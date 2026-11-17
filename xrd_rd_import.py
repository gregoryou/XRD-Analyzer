"""Self-validating reader for Philips / PANalytical ``.rd`` files (``V3RD``).

STATUS: checked against ONE real file (PW3064 / PW3050/60 Theta-Theta, polycrystalline Si standard):
810-byte header, 16-bit little-endian counts, float64 values for K-alpha1 / K-alpha2 / ratio
(offsets 94 / 102 / 110) and for step / start / end (offsets 214 / 222 / 230). Nothing else is
assumed, and nothing is taken on trust: the scan axis is accepted only when three numbers stored
next to each other in the header (start, step, end) agree with the number of data points in the
file, and the counts must look like detector counts.

In that real file the header describes 2992 points but only 2711 are stored, i.e. exactly 281
points are missing. The silicon reflections show that the stored points end at the header's "end"
angle, so when the discrepancy is a whole number of steps the axis is anchored at "end" and a
warning says so. This rests on a single file: verify the axis on a known reflection.
When any check fails a ``ValueError`` explains what was found and nothing is returned. ``V5RD``
files are recognised but not read (the layout is unknown without a sample).

Run ``python xrd_rd_import.py file.rd`` to print the diagnostics for a file.
"""
from __future__ import annotations

import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

V3_MAGIC = b"V3RD"
V5_MAGIC = b"V5RD"
V3_DATA_OFFSET = 810            # size of the V3 header; 16-bit little-endian counts follow
MAX_TRAILER_POINTS = 16         # a few extra values after the scan are tolerated
BLOCK_WINDOW_BYTES = 24         # start, step and end are stored next to each other
MIN_POINTS = 20

# anode: (K-alpha1, K-alpha2) in angstrom
KNOWN_LINES: Dict[str, Tuple[float, float]] = {
    "Cu": (1.540598, 1.544426), "Co": (1.788965, 1.792850), "Cr": (2.289700, 2.293606),
    "Fe": (1.936042, 1.939980), "Mo": (0.709300, 0.713590), "Ni": (1.657910, 1.661750),
    "Ag": (0.559421, 0.563798),
}
LINE_TOLERANCE = 3e-4


@dataclass
class RdScan:
    x: np.ndarray
    y: np.ndarray
    version: str
    wavelength_A: Optional[float]
    wavelength_kalpha2_A: Optional[float]
    anode: str
    header_text: str
    layout: Dict[str, object] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)


# ----------------------------------------------------------------------------- header scanning
def _scan_numbers(header: bytes, dtype: str) -> Tuple[np.ndarray, np.ndarray]:
    """Return (offsets, values) of every finite number of the given dtype at every byte offset."""
    size = np.dtype(dtype).itemsize
    offsets: List[np.ndarray] = []
    values: List[np.ndarray] = []
    for shift in range(size):
        count = (len(header) - shift) // size
        if count <= 0:
            continue
        array = np.frombuffer(header, dtype=dtype, count=count, offset=shift)
        index = np.nonzero(np.isfinite(array))[0]
        offsets.append(shift + index * size)
        values.append(array[index].astype(float))
    if not offsets:
        return np.zeros(0, dtype=int), np.zeros(0)
    offs, vals = np.concatenate(offsets), np.concatenate(values)
    order = np.argsort(offs, kind="stable")
    return offs[order], vals[order]


def find_axis_blocks(header: bytes, point_counts) -> List[Tuple[float, float, int, int, str]]:
    """Find (start, step, n_points, offset, dtype) whose start/step/end are consistent.

    start, step and end must lie within BLOCK_WINDOW_BYTES of each other and satisfy
    end = start + (n - 1) * step (or start + n * step) for one of the candidate point counts.
    Identical blocks (a float32 and a float64 copy, or a repeated copy) are merged; blocks that
    disagree are all returned so that the caller can refuse to guess. When the end value could
    mean either "last point" or "end of the last step", the longer reading (last point) is used.
    """
    point_counts = [int(n) for n in point_counts if int(n) >= MIN_POINTS]
    found: List[Tuple[float, float, int, int, str]] = []
    for dtype in ("<f4", "<f8"):
        offs, vals = _scan_numbers(header, dtype)
        if len(offs) == 0:
            continue
        starts = np.nonzero((vals > 0.1) & (vals < 170.0))[0]
        for i in starts:
            near = np.nonzero(np.abs(offs - offs[i]) <= BLOCK_WINDOW_BYTES)[0]
            steps = [j for j in near if j != i and 5e-4 <= vals[j] <= 5.0]
            ends = [k for k in near if k != i and 0.1 < vals[k] <= 180.0 and vals[k] > vals[i]]
            for j in steps:
                for k in ends:
                    if k == j:
                        continue
                    s, d, e = vals[i], vals[j], vals[k]
                    tol = max(2e-4, 0.002 * d)
                    for n in point_counts:
                        if abs(s + (n - 1) * d - e) <= tol or abs(s + n * d - e) <= tol:
                            found.append((float(s), float(d), n, int(offs[i]), dtype))
    # The same start/step/end matches N points (end = last point) and N-1 points (end = start + N*step
    # for N-1 points). Blocks with the same start and step are one scan: keep the longest reading,
    # i.e. "end = last point". Blocks with different start or step are genuinely different.
    merged: List[Tuple[float, float, int, int, str]] = []
    for block in sorted(found, key=lambda b: -b[2]):
        if not any(abs(block[0] - m[0]) < 1e-3 and abs(block[1] - m[1]) < 1e-5 for m in merged):
            merged.append(block)
    return merged


def find_trimmed_blocks(header: bytes, n_points: int) -> List[Tuple[float, float, float, int, int, str]]:
    """Find (start, step, end) describing MORE points than the file holds, by a whole number of steps.

    Returns (start, step, end, n_header, offset, dtype). Blocks with the same start and step are merged.
    """
    found: List[Tuple[float, float, float, int, int, str]] = []
    for dtype in ("<f4", "<f8"):
        offs, vals = _scan_numbers(header, dtype)
        if len(offs) == 0:
            continue
        for i in np.nonzero((vals > 0.1) & (vals < 170.0))[0]:
            near = np.nonzero(np.abs(offs - offs[i]) <= BLOCK_WINDOW_BYTES)[0]
            for j in near:
                if j == i or not (5e-4 <= vals[j] <= 5.0):
                    continue
                for k in near:
                    if k in (i, j) or not (vals[i] < vals[k] <= 180.0):
                        continue
                    n_header = (vals[k] - vals[i]) / vals[j] + 1.0
                    if abs(n_header - round(n_header)) <= 0.02 and n_points < round(n_header) <= 2 * n_points:
                        found.append((float(vals[i]), float(vals[j]), float(vals[k]), int(round(n_header)),
                                      int(offs[i]), dtype))
    merged: List[Tuple[float, float, float, int, int, str]] = []
    for block in found:
        if not any(abs(block[0] - m[0]) < 1e-3 and abs(block[1] - m[1]) < 1e-5 for m in merged):
            merged.append(block)
    return merged


def find_wavelengths(header: bytes) -> Optional[Tuple[str, float, Optional[float]]]:
    """Look for K-alpha wavelengths of a known anode stored as float32/float64 in the header."""
    hits: Dict[str, Dict[str, float]] = {}
    for dtype in ("<f4", "<f8"):
        offs, vals = _scan_numbers(header, dtype)
        for anode, (ka1, ka2) in KNOWN_LINES.items():
            for name, nominal in (("ka1", ka1), ("ka2", ka2)):
                match = np.nonzero(np.abs(vals - nominal) <= LINE_TOLERANCE)[0]
                if len(match) and name not in hits.setdefault(anode, {}):
                    hits[anode][name] = float(vals[match[0]])
    hits = {anode: found for anode, found in hits.items() if "ka1" in found}
    if not hits:
        return None
    anode = max(hits, key=lambda a: (len(hits[a]), a == "Cu"))
    return anode, hits[anode]["ka1"], hits[anode].get("ka2")


def header_strings(header: bytes, minimum: int = 4, limit: int = 8) -> List[str]:
    runs = re.findall(rb"[ -~]{%d,}" % minimum, header[4:])
    cleaned = [r.decode("ascii").strip() for r in runs]
    return [c for c in cleaned if c][:limit]


def count_plausibility(y) -> Tuple[bool, float, str]:
    """Do the values look like detector counts (non-negative, Poisson-like point-to-point scatter)?"""
    y = np.asarray(y, dtype=float)
    if len(y) < MIN_POINTS or not np.all(np.isfinite(y)):
        return False, float("nan"), "zbyt mało punktów lub wartości nieskończone"
    if y.min() < 0:
        return False, float("nan"), "ujemne zliczenia"
    if float(np.ptp(y)) == 0.0:
        return False, float("nan"), "wszystkie wartości są równe"
    level = float(np.median(y))
    ratio = float(np.median(np.abs(np.diff(y))) / math.sqrt(max(level, 1.0)))
    if level >= 20.0 and not (0.05 <= ratio <= 30.0):
        return False, ratio, (f"rozrzut sąsiednich punktów ({ratio:.3g}·√poziom) nie przypomina zliczeń "
                              "(oczekiwano 0,05–30)")
    return True, ratio, "ok"


# ----------------------------------------------------------------------------- reading
def _diagnostics_hint(path: Path) -> str:
    return (f"Uruchom `python xrd_rd_import.py \"{path.name}\"` i prześlij wynik (lub sam plik) autorowi programu, "
            "aby dopisać obsługę tego wariantu.")


def parse_philips_rd(path) -> RdScan:
    path = Path(path)
    data = path.read_bytes()
    if len(data) < 64:
        raise ValueError("Plik .rd jest zbyt krótki.")
    magic = data[:4]
    if magic == V5_MAGIC:
        raise ValueError("Rozpoznano plik Philips RD w wersji V5RD, ale jego układ nie jest jeszcze obsługiwany "
                         "(bez pliku wzorcowego nie zgaduję osi kątowej). " + _diagnostics_hint(path))
    if magic != V3_MAGIC:
        raise ValueError("To nie jest plik Philips RD (oczekiwano znacznika V3RD na początku pliku).")
    if len(data) < V3_DATA_OFFSET + 2 * MIN_POINTS:
        raise ValueError("Plik V3RD nie zawiera danych za nagłówkiem (810 bajtów).")

    header = data[:V3_DATA_OFFSET]
    n_available = (len(data) - V3_DATA_OFFSET) // 2
    candidates = range(n_available, max(n_available - MAX_TRAILER_POINTS, MIN_POINTS) - 1, -1)
    blocks = find_axis_blocks(header, candidates)
    anchor, missing = "start", 0
    axis_end = None
    if not blocks:
        trimmed = find_trimmed_blocks(header, n_available)
        if len(trimmed) == 1:
            t_start, t_step, axis_end, n_header, offset, dtype = trimmed[0]
            blocks = [(t_start, t_step, n_available, offset, dtype)]
            anchor, missing = "end", n_header - n_available
        elif len(trimmed) > 1:
            listing = "; ".join(f"start {b[0]:.5g}°, krok {b[1]:.6g}°, koniec {b[2]:.6g}°, N={b[3]}" for b in trimmed[:4])
            raise ValueError("Nagłówek V3RD opisuje więcej punktów niż zawiera plik, a zakresów pasujących do tego jest "
                             f"kilka, więc nie wybieram żadnego: {listing}. " + _diagnostics_hint(path))
    if not blocks:
        raise ValueError(
            "Nie znaleziono w nagłówku V3RD trzech sąsiadujących liczb (start, krok, koniec) zgodnych z liczbą "
            f"punktów w pliku ({n_available}). Oś kątowa nie została odgadnięta. " + _diagnostics_hint(path))
    if len(blocks) > 1:
        listing = "; ".join(f"start {b[0]:.5g}°, krok {b[1]:.6g}°, N={b[2]} (przesunięcie {b[3]}, {b[4]})" for b in blocks[:4])
        raise ValueError("Nagłówek V3RD zawiera kilka różnych, wewnętrznie spójnych zakresów osi, więc nie wybieram "
                         f"żadnego: {listing}. " + _diagnostics_hint(path))
    start, step, n, offset, dtype = blocks[0]
    if anchor == "end":
        start = axis_end - (n - 1) * step

    counts = np.frombuffer(data, dtype="<u2", count=n, offset=V3_DATA_OFFSET).astype(float)
    plausible, ratio, reason = count_plausibility(counts)
    if not plausible:
        raise ValueError(f"Oś kątowa odczytana z nagłówka jest spójna, ale dane nie wyglądają na zliczenia: {reason}. "
                         + _diagnostics_hint(path))

    warnings: List[str] = ["Plik .rd odczytano heurystycznie (układ potwierdzono spójnością zakresu, kroku i liczby "
                           "punktów); sprawdź oś 2θ i długość fali."]
    if n_available - n > 0:
        warnings.append(f"Pominięto {n_available - n} dodatkowych wartości za skanem.")
    if anchor == "end":
        warnings.append(f"Nagłówek opisuje o {missing} punktów więcej, niż zawiera plik. Przyjęto, że brakuje ich na "
                        "początku skanu (oś zakotwiczona na kącie końcowym) — potwierdzono na jednym pliku wzorca Si; "
                        "zweryfikuj oś na znanym refleksie.")
    found = find_wavelengths(header)
    anode, wavelength, wavelength2 = ("", None, None) if found is None else found
    if found is None:
        warnings.append("Nie znaleziono długości fali w nagłówku; przyjęto domyślną (Cu Kα1) — zweryfikuj ją.")
    elif wavelength2 is None:
        warnings.append("W nagłówku znaleziono tylko Kα1; sprawdź długość fali.")
    text = "; ".join(header_strings(header))[:300]
    layout = {"version": "3", "data_offset": V3_DATA_OFFSET, "dtype": "uint16 LE", "n_points": n,
              "axis_block_offset": offset, "axis_block_dtype": dtype, "start_deg": start, "step_deg": step,
              "anchor": anchor, "missing_points": missing,
              "poisson_ratio": ratio}
    return RdScan(start + step * np.arange(n), counts, "3", wavelength, wavelength2, anode, text, layout, warnings)


def describe_philips_rd(path) -> str:
    """Human-readable diagnostics, meant to be attached to a bug report."""
    path = Path(path)
    data = path.read_bytes()
    lines = [f"Plik: {path.name}  rozmiar: {len(data)} B  znacznik: {data[:4]!r}"]
    header = data[:V3_DATA_OFFSET] if data[:4] == V3_MAGIC else data[:4096]
    lines.append(f"Teksty w nagłówku: {header_strings(header, limit=20)}")
    found = find_wavelengths(header)
    lines.append(f"Długość fali (znane linie anod): {found}")
    if data[:4] == V3_MAGIC and len(data) > V3_DATA_OFFSET:
        n_available = (len(data) - V3_DATA_OFFSET) // 2
        lines.append(f"Punktów 16-bitowych za 810-bajtowym nagłówkiem: {n_available}")
        blocks = find_axis_blocks(header, range(n_available, max(n_available - MAX_TRAILER_POINTS, MIN_POINTS) - 1, -1))
        lines.append(f"Spójne bloki (start, krok, N, przesunięcie, typ): {blocks}")
        lines.append("Bloki opisujące więcej punktów niż plik (start, krok, koniec, N, przesunięcie, typ): "
                     f"{find_trimmed_blocks(header, n_available)}")
        counts = np.frombuffer(data, dtype="<u2", count=n_available, offset=V3_DATA_OFFSET).astype(float)
        lines.append(f"Zliczenia: min {counts.min():.0f} max {counts.max():.0f} mediana {np.median(counts):.0f} "
                     f"wiarygodność: {count_plausibility(counts)}")
    else:
        lines.append("Układ danych nie jest znany; dołącz ten wydruk do zgłoszenia.")
        lines.append("Początek pliku (hex): " + data[:96].hex(" "))
    try:
        scan = parse_philips_rd(path)
        lines.append(f"ODCZYT POWIÓDŁ SIĘ: {len(scan.x)} punktów, 2θ {scan.x[0]:.4f}–{scan.x[-1]:.4f}°, "
                     f"λ={scan.wavelength_A}, anoda={scan.anode or '?'}")
    except ValueError as exc:
        lines.append(f"ODCZYT NIE POWIÓDŁ SIĘ: {exc}")
    return "\n".join(lines)


if __name__ == "__main__":
    for argument in sys.argv[1:]:
        print(describe_philips_rd(argument))
        print()
