"""Calibration of a measurement series with an external powder standard (silicon).

A polycrystalline Si standard scanned together with the samples gives two things:

* the 2theta ERROR of the instrument as a function of angle. Observed minus theoretical positions of
  the Si reflections are fitted with  d(2theta) = p0 + p1*cos(theta) + p2*sin(2theta)
  (zero shift, sample displacement, transparency); the terms are chosen by AICc. A sample's peak
  positions are corrected as  2theta_true = 2theta_obs - d(2theta_obs).
* the instrumental line width (Gaussian and Lorentzian FWHM at each Si reflection), which is the IPF
  needed for Scherrer / Williamson-Hall.

Caveats that matter physically (they are also written into the stored notes):
  - The cos(theta) term is the sample displacement of the STANDARD. It transfers to a sample only if
    the standard and the sample were mounted at the same height (same holder, same spinner). The zero
    shift transfers whenever the instrument alignment did not change between the scans.
  - An uncertified standard has an uncertain lattice parameter: the literature value is used unless
    the certificate value is entered, and any error in a shifts the zero term.
  - The correction is only valid inside the angular range covered by the Si reflections.
"""
from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field
from datetime import date
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from xrd_profile_analysis import InstrumentBroadeningModel, fit_peak_profile

SI_LATTICE_A = 5.4310205          # angstrom, literature value near 22.5 degC (used for uncertified Si)
SYSTEMATIC_FLOOR_DEG = 0.002      # added in quadrature to the fit uncertainty of each reflection
MIN_REFLECTIONS = 3
MODELS: Dict[str, Tuple[str, ...]] = {
    "zero": ("zero",),
    "zero+cos": ("zero", "cos_theta"),
    "zero+sin2": ("zero", "sin_2theta"),
    "zero+cos+sin2": ("zero", "cos_theta", "sin_2theta"),
}
MODEL_LABELS = {
    "auto": "Automatyczny (AICc)",
    "zero": "Tylko zero",
    "zero+cos": "Zero + cosθ (przesunięcie próbki)",
    "zero+sin2": "Zero + sin2θ (transparencja)",
    "zero+cos+sin2": "Zero + cosθ + sin2θ",
}
TERM_LABELS = {"zero": "zero", "cos_theta": "cosθ", "sin_2theta": "sin2θ"}


# ----------------------------------------------------------------------------- silicon geometry
def si_reflection_allowed(h: int, k: int, l: int) -> bool:
    """Diamond structure: all odd, or all even with h+k+l divisible by 4."""
    parity = {h % 2, k % 2, l % 2}
    if parity == {1}:
        return True
    return parity == {0} and (h + k + l) % 4 == 0 and (h, k, l) != (0, 0, 0)


def si_reflections(wavelength_A: float, a_A: float = SI_LATTICE_A, two_theta_min: float = 5.0,
                   two_theta_max: float = 170.0) -> List[Tuple[str, float]]:
    """(hkl label, theoretical Kalpha1 2theta) of every distinct Si reflection in the range."""
    by_n: Dict[int, Tuple[int, int, int]] = {}
    for h in range(0, 13):
        for k in range(0, h + 1):
            for l in range(0, k + 1):
                if si_reflection_allowed(h, k, l):
                    n = h * h + k * k + l * l
                    if n not in by_n or (h, k, l) > by_n[n]:
                        by_n[n] = (h, k, l)
    rows = []
    for n, (h, k, l) in sorted(by_n.items()):
        s = wavelength_A * math.sqrt(n) / (2.0 * a_A)
        if s >= 1.0:
            continue
        tt = 2.0 * math.degrees(math.asin(s))
        if two_theta_min <= tt <= two_theta_max:
            rows.append((f"{h}{k}{l}", tt))
    return rows


# ----------------------------------------------------------------------------- reflection fits
@dataclass
class StandardReflection:
    hkl: str
    two_theta_theory: float
    two_theta_obs: float
    u_two_theta: float
    fwhm_obs: float
    gaussian_fwhm: float
    lorentzian_fwhm: float
    u_gaussian_fwhm: Optional[float]
    u_lorentzian_fwhm: Optional[float]
    r_squared: Optional[float]
    used: bool
    note: str = ""

    @property
    def delta(self) -> float:
        return self.two_theta_obs - self.two_theta_theory

    def to_dict(self) -> dict:
        """JSON-safe: a failed fit has NaN values, which are stored as null."""
        data = {key: (None if isinstance(value, float) and not math.isfinite(value) else value)
                for key, value in self.__dict__.items()}
        delta = self.delta
        data["delta"] = delta if math.isfinite(delta) else None
        return data

    @classmethod
    def from_dict(cls, data: dict) -> "StandardReflection":
        numeric = ("two_theta_theory", "two_theta_obs", "u_two_theta", "fwhm_obs", "gaussian_fwhm", "lorentzian_fwhm")
        values = {key: data.get(key) for key in cls.__dataclass_fields__ if key in data}
        for key in numeric:
            if values.get(key) is None:
                values[key] = math.nan
        return cls(**values)


def _finite(value) -> Optional[float]:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def fit_standard_reflections(x, y, wavelength_A: float, wavelength_kalpha2_A: float, kalpha2_ratio: float = 0.5,
                             a_A: float = SI_LATTICE_A, min_r_squared: float = 0.95,
                             max_offset_deg: float = 0.35) -> List[StandardReflection]:
    """Fit every Si reflection that lies fully inside the scan; flag the ones that are usable."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    rows: List[StandardReflection] = []
    for hkl, theory in si_reflections(wavelength_A, a_A, float(x.min()), float(x.max())):
        theta = math.radians(theory / 2.0)
        doublet = math.degrees(2.0 * math.tan(theta) * (wavelength_kalpha2_A - wavelength_A) / wavelength_A)
        half = max(0.6, 0.4 + 3.5 * doublet)
        if theory - half < x[0] or theory + half > x[-1]:
            continue
        window = (x >= theory - half) & (x <= theory + half)
        if int(window.sum()) < 12:
            continue
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fit = fit_peak_profile(x[window], y[window], center_hint=theory, wavelength_A=wavelength_A,
                                       wavelength_kalpha2_A=wavelength_kalpha2_A, kalpha2_ratio=kalpha2_ratio, hkl=hkl)
        except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
            rows.append(StandardReflection(hkl, theory, math.nan, math.nan, math.nan, math.nan, math.nan, None, None,
                                           None, False, f"dopasowanie nieudane: {exc}"[:120]))
            continue
        data = fit.to_dict()
        cov = data.get("covariance")
        u_g = math.sqrt(max(0.0, cov[2][2])) if cov else None
        u_l = math.sqrt(max(0.0, cov[3][3])) if cov else None
        u_tt = _finite(data.get("u_two_theta_deg"))
        centre = float(data["two_theta_kalpha1_deg"])
        r2 = _finite(data.get("r_squared"))
        problems = []
        if not data.get("fit_success"):
            problems.append("dopasowanie niezbieżne")
        if r2 is None or r2 < min_r_squared:
            problems.append(f"R² < {min_r_squared:g}")
        if abs(centre - theory) > max_offset_deg:
            problems.append(f"położenie odbiega o {centre - theory:+.2f}°")
        if not 0.02 <= float(data["fwhm_obs_deg"]) <= 1.2:
            problems.append("nieprawdopodobna szerokość")
        rows.append(StandardReflection(
            hkl, theory, centre, u_tt if u_tt is not None else 0.005, float(data["fwhm_obs_deg"]),
            float(data["gaussian_fwhm_deg"]), float(data["lorentzian_fwhm_deg"]), u_g, u_l, r2,
            not problems, "; ".join(problems)))
    return rows


# ----------------------------------------------------------------------------- axis correction
def _basis(terms: Sequence[str], two_theta) -> np.ndarray:
    tt = np.atleast_1d(np.asarray(two_theta, dtype=float))
    columns = []
    for term in terms:
        if term == "zero":
            columns.append(np.ones_like(tt))
        elif term == "cos_theta":
            columns.append(np.cos(np.radians(tt / 2.0)))
        elif term == "sin_2theta":
            columns.append(np.sin(np.radians(tt)))
        else:
            raise ValueError(f"Nieznany człon korekty: {term}")
    return np.column_stack(columns)


@dataclass
class AxisCorrection:
    terms: List[str]
    params: List[float]
    u_params: List[Optional[float]]
    covariance: Optional[List[List[float]]]
    rms_deg: float
    chi2_reduced: Optional[float]
    n_points: int
    two_theta_min: float
    two_theta_max: float
    model: str = ""
    aicc: Dict[str, float] = field(default_factory=dict)

    def delta(self, two_theta):
        """Instrument error d(2theta) = observed - true, same shape as the input."""
        tt = np.asarray(two_theta, dtype=float)
        value = _basis(self.terms, tt.reshape(-1)) @ np.asarray(self.params)
        return value.reshape(tt.shape) if tt.ndim else float(value[0])

    def u_delta(self, two_theta):
        tt = np.asarray(two_theta, dtype=float)
        if not self.covariance:
            return np.full(tt.shape, math.nan) if tt.ndim else math.nan
        design = _basis(self.terms, tt.reshape(-1))
        variance = np.einsum("ij,jk,ik->i", design, np.asarray(self.covariance), design)
        value = np.sqrt(np.maximum(variance, 0.0))
        return value.reshape(tt.shape) if tt.ndim else float(value[0])

    def corrected(self, two_theta):
        tt = np.asarray(two_theta, dtype=float)
        return tt - self.delta(tt)

    def raw_from_corrected(self, corrected):
        """Inverse of corrected(): the instrument error is smooth and tiny, so fixed-point iteration converges."""
        value = np.asarray(corrected, dtype=float)
        raw = value.copy()
        for _ in range(4):
            raw = value + self.delta(raw)
        return raw

    def covers(self, two_theta, margin_deg: float = 3.0):
        tt = np.asarray(two_theta, dtype=float)
        return (tt >= self.two_theta_min - margin_deg) & (tt <= self.two_theta_max + margin_deg)

    def describe(self) -> str:
        parts = []
        for term, value, u in zip(self.terms, self.params, self.u_params):
            parts.append(f"{value:+.4f}" + (f"±{u:.4f}" if u is not None else "") + ("" if term == "zero" else f"·{TERM_LABELS[term]}"))
        return "Δ(2θ) = " + " ".join(parts) + f" [°];  RMS reszt {self.rms_deg:.4f}°, n={self.n_points}"

    def to_dict(self) -> dict:
        return {"terms": list(self.terms), "params": [float(p) for p in self.params],
                "u_params": [None if u is None else float(u) for u in self.u_params],
                "covariance": None if self.covariance is None else [[float(v) for v in row] for row in self.covariance],
                "rms_deg": float(self.rms_deg), "chi2_reduced": None if self.chi2_reduced is None else float(self.chi2_reduced),
                "n_points": int(self.n_points), "two_theta_min": float(self.two_theta_min),
                "two_theta_max": float(self.two_theta_max), "model": self.model,
                "aicc": {k: float(v) for k, v in self.aicc.items()}}

    @classmethod
    def from_dict(cls, data: dict) -> "AxisCorrection":
        return cls(terms=list(data["terms"]), params=[float(p) for p in data["params"]],
                   u_params=[None if u is None else float(u) for u in data.get("u_params", [None] * len(data["params"]))],
                   covariance=data.get("covariance"), rms_deg=float(data.get("rms_deg", 0.0)),
                   chi2_reduced=data.get("chi2_reduced"), n_points=int(data.get("n_points", 0)),
                   two_theta_min=float(data["two_theta_min"]), two_theta_max=float(data["two_theta_max"]),
                   model=str(data.get("model", "")), aicc=dict(data.get("aicc", {})))


def fit_axis_correction(two_theta_obs, two_theta_theory, u_two_theta=None, model: str = "auto") -> AxisCorrection:
    """Weighted least squares of (observed - theoretical) against the chosen terms.

    model="auto" picks the term set with the lowest AICc among those the data can support
    (n >= terms + 2). A systematic floor is added to every uncertainty so that one very sharp
    reflection cannot dominate the fit.
    """
    obs = np.asarray(two_theta_obs, dtype=float)
    theory = np.asarray(two_theta_theory, dtype=float)
    n = len(obs)
    if n < MIN_REFLECTIONS or len(theory) != n:
        raise ValueError(f"Korekta osi wymaga co najmniej {MIN_REFLECTIONS} użytecznych refleksów wzorca (jest {n}).")
    delta = obs - theory
    u = np.full(n, 0.005) if u_two_theta is None else np.asarray(u_two_theta, dtype=float)
    u = np.sqrt(np.where(np.isfinite(u) & (u > 0), u, 0.005) ** 2 + SYSTEMATIC_FLOOR_DEG ** 2)
    weights = 1.0 / u ** 2

    def solve(terms):
        design = _basis(terms, obs)
        a = design * np.sqrt(weights)[:, None]
        b = delta * np.sqrt(weights)
        params, *_ = np.linalg.lstsq(a, b, rcond=None)
        residual = delta - design @ params
        chi2 = float(np.sum(weights * residual ** 2))
        return params, residual, chi2, a

    eligible = {name: terms for name, terms in MODELS.items() if n >= len(terms) + 2}
    if model != "auto":
        if model not in MODELS:
            raise ValueError(f"Nieznany model korekty: {model}")
        if model not in eligible:
            raise ValueError(f"Model „{MODEL_LABELS[model]}” wymaga co najmniej {len(MODELS[model]) + 2} refleksów (jest {n}).")
        chosen = model
        aicc = {}
    else:
        aicc = {}
        for name, terms in eligible.items():
            _p, _r, chi2, _a = solve(terms)
            k = len(terms)
            aicc[name] = chi2 + 2 * k + 2 * k * (k + 1) / (n - k - 1)
        best = min(aicc.values())
        # prefer the simplest model within 2 AICc units of the best one
        chosen = next(name for name in MODELS if name in aicc and aicc[name] <= best + 2.0)
    terms = MODELS[chosen]
    params, residual, chi2, a = solve(terms)
    dof = n - len(terms)
    chi2_red = chi2 / dof if dof > 0 else None
    cov = np.linalg.inv(a.T @ a) * (max(1.0, chi2_red) if chi2_red is not None else 1.0)
    u_params = [float(math.sqrt(max(0.0, cov[i, i]))) for i in range(len(terms))]
    return AxisCorrection(list(terms), [float(p) for p in params], u_params, cov.tolist(),
                          float(math.sqrt(np.mean(residual ** 2))), chi2_red, n, float(obs.min()), float(obs.max()),
                          chosen, aicc)


# ----------------------------------------------------------------------------- the whole calibration
@dataclass
class StandardCalibration:
    axis: AxisCorrection
    ipf: Optional[dict]
    reflections: List[StandardReflection]
    standard: str
    a_nominal_A: float
    certified: bool
    wavelength_A: float
    wavelength_kalpha2_A: float
    kalpha2_ratio: float
    source_name: str = ""
    source_path: str = ""
    created: str = ""
    notes: List[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {"axis": self.axis.to_dict(), "ipf": self.ipf, "reflections": [r.to_dict() for r in self.reflections],
                "standard": self.standard, "a_nominal_A": float(self.a_nominal_A), "certified": bool(self.certified),
                "wavelength_A": float(self.wavelength_A), "wavelength_kalpha2_A": float(self.wavelength_kalpha2_A),
                "kalpha2_ratio": float(self.kalpha2_ratio), "source_name": self.source_name,
                "source_path": self.source_path, "created": self.created, "notes": list(self.notes)}

    @classmethod
    def from_dict(cls, data: dict) -> "StandardCalibration":
        return cls(axis=AxisCorrection.from_dict(data["axis"]), ipf=data.get("ipf"),
                   reflections=[StandardReflection.from_dict(r) for r in data.get("reflections", [])],
                   standard=str(data.get("standard", "Si")), a_nominal_A=float(data.get("a_nominal_A", SI_LATTICE_A)),
                   certified=bool(data.get("certified", False)), wavelength_A=float(data["wavelength_A"]),
                   wavelength_kalpha2_A=float(data.get("wavelength_kalpha2_A", 0.0)),
                   kalpha2_ratio=float(data.get("kalpha2_ratio", 0.5)), source_name=str(data.get("source_name", "")),
                   source_path=str(data.get("source_path", "")), created=str(data.get("created", "")),
                   notes=list(data.get("notes", [])))

    def compact(self) -> dict:
        """What is stored in every sample: model, IPF, and the reflection table without fit arrays."""
        return self.to_dict()


def build_instrument_profile(reflections: Sequence[StandardReflection], name: str, wavelength_A: float,
                             source_file: str = "") -> Optional[InstrumentBroadeningModel]:
    """Instrument broadening model (Gaussian / Lorentzian FWHM vs 2theta) from the usable reflections."""
    used = sorted((r for r in reflections if r.used), key=lambda r: r.two_theta_obs)
    unique = []
    for r in used:
        if not unique or r.two_theta_obs - unique[-1].two_theta_obs > 1e-6:
            unique.append(r)
    if len(unique) < 2:
        return None
    values = dict(two_theta_deg=[r.two_theta_obs for r in unique], gaussian_fwhm_deg=[r.gaussian_fwhm for r in unique],
                  lorentzian_fwhm_deg=[r.lorentzian_fwhm for r in unique], name=name, standard=name, optics_id=name,
                  wavelength_A=wavelength_A, source_file=source_file)
    if all(r.u_gaussian_fwhm is not None and r.u_lorentzian_fwhm is not None for r in unique):
        values.update(u_gaussian_fwhm_deg=[r.u_gaussian_fwhm for r in unique],
                      u_lorentzian_fwhm_deg=[r.u_lorentzian_fwhm for r in unique])
    return InstrumentBroadeningModel(**values)


def calibrate_from_reflections(reflections: Sequence[StandardReflection], wavelength_A: float, wavelength_kalpha2_A: float,
                               kalpha2_ratio: float = 0.5, a_A: float = SI_LATTICE_A, certified: bool = False,
                               model: str = "auto", source_name: str = "", source_path: str = "") -> StandardCalibration:
    used = [r for r in reflections if r.used]
    axis = fit_axis_correction([r.two_theta_obs for r in used], [r.two_theta_theory for r in used],
                               [r.u_two_theta for r in used], model)
    label = f"Si {source_name}".strip()
    ipf = build_instrument_profile(used, label, wavelength_A, source_path)
    notes = []
    if not certified:
        notes.append(f"Wzorzec niecertyfikowany: przyjęto parametr sieci Si a = {a_A:.5f} Å (wartość literaturowa). "
                     "Błąd a przesuwa człon zerowy.")
    if "cos_theta" in axis.terms:
        notes.append("Człon cosθ to przesunięcie próbki WZORCA. Dla próbki ma sens tylko wtedy, gdy wzorzec i próbka "
                     "były zamocowane na tej samej wysokości.")
    if ipf is None:
        notes.append("Za mało refleksów do profilu aparatury (potrzeba co najmniej 2).")
    else:
        notes.append(f"Profil aparatury obejmuje 2θ = {ipf.two_theta_deg[0]:.1f}–{ipf.two_theta_deg[-1]:.1f}°; "
                     "poza tym zakresem nie jest ekstrapolowany.")
    if len(used) < len(reflections):
        skipped = ", ".join(f"{r.hkl} ({r.note or 'wyłączony'})" for r in reflections if not r.used)
        notes.append(f"Pominięte refleksy: {skipped}.")
    return StandardCalibration(axis, None if ipf is None else ipf.to_dict(), list(reflections), "Si", float(a_A),
                               bool(certified), float(wavelength_A), float(wavelength_kalpha2_A), float(kalpha2_ratio),
                               source_name, source_path, date.today().isoformat(), notes)


def calibrate_from_standard(x, y, wavelength_A: float, wavelength_kalpha2_A: float, kalpha2_ratio: float = 0.5,
                            a_A: float = SI_LATTICE_A, certified: bool = False, model: str = "auto",
                            source_name: str = "", source_path: str = "") -> StandardCalibration:
    reflections = fit_standard_reflections(x, y, wavelength_A, wavelength_kalpha2_A, kalpha2_ratio, a_A)
    usable = sum(r.used for r in reflections)
    if usable < MIN_REFLECTIONS:
        found = ", ".join(f"{r.hkl}: {r.note or 'ok'}" for r in reflections) or "żadnego refleksu Si w zakresie skanu"
        raise ValueError(f"Znaleziono {usable} użytecznych refleksów Si, a potrzeba co najmniej {MIN_REFLECTIONS}. "
                         f"Szczegóły: {found}.")
    return calibrate_from_reflections(reflections, wavelength_A, wavelength_kalpha2_A, kalpha2_ratio, a_A, certified,
                                      model, source_name, source_path)


# ----------------------------------------------------------------------------- using the calibration
ROLE_KEY = "calibration_role"          # "standard" marks the Si scan itself
CAL_KEY = "calibration"                # compact StandardCalibration stored in every sample
APPLIED_KEY = "calibration_applied"    # what has actually been applied to the stored peak positions
WAVELENGTH_TOLERANCE_A = 1e-4


def calibration_status(dataset) -> Tuple[Optional[StandardCalibration], str]:
    """(calibration, "") when it applies to this dataset, otherwise (calibration or None, reason)."""
    md = dataset.metadata
    raw = md.get(CAL_KEY)
    if not raw:
        return None, "brak kalibracji"
    try:
        calibration = StandardCalibration.from_dict(raw)
    except (KeyError, TypeError, ValueError) as exc:
        return None, f"uszkodzone dane kalibracji: {exc}"
    if md.get(ROLE_KEY) == "standard":
        return calibration, "to jest sam wzorzec"
    if abs(float(dataset.wavelength_A) - calibration.wavelength_A) > WAVELENGTH_TOLERANCE_A:
        return calibration, (f"inna długość fali (próbka {dataset.wavelength_A:.5f} Å, wzorzec "
                             f"{calibration.wavelength_A:.5f} Å)")
    return calibration, ""


def calibrated_axis(dataset) -> Tuple[np.ndarray, Optional[StandardCalibration]]:
    """The 2theta axis with the standard's correction applied (or the raw axis and None)."""
    calibration, reason = calibration_status(dataset)
    x = np.asarray(dataset.x, dtype=float)
    if calibration is None or reason:
        return x, None
    return calibration.axis.corrected(x), calibration


def _shift_peaks(peaks, wavelength_A: float, axis: AxisCorrection, undo: bool) -> int:
    count = 0
    for p in peaks:
        if undo:
            p.two_theta = float(axis.raw_from_corrected(p.two_theta))
        else:
            p.two_theta = float(axis.corrected(p.two_theta))
        s = math.sin(math.radians(p.two_theta / 2.0))
        if s > 0:
            p.d_A = wavelength_A / (2.0 * s)
        count += 1
    return count


def apply_calibration_to_fresh_peaks(dataset) -> Optional[dict]:
    """Peaks were just detected on the raw axis: correct their positions and d. Returns what was applied."""
    dataset.metadata.pop(APPLIED_KEY, None)
    dataset.metadata.pop("calibration_skipped", None)
    calibration, reason = calibration_status(dataset)
    if calibration is None:
        if reason not in ("brak kalibracji", "to jest sam wzorzec") and dataset.metadata.get(CAL_KEY):
            dataset.metadata["calibration_skipped"] = reason
        return None
    if reason:
        if reason != "to jest sam wzorzec":
            dataset.metadata["calibration_skipped"] = reason
        return None
    return _apply(dataset, calibration)


def _apply(dataset, calibration: StandardCalibration) -> dict:
    outside = [p.two_theta for p in dataset.peaks if not bool(calibration.axis.covers(p.two_theta))]
    n = _shift_peaks(dataset.peaks, dataset.wavelength_A, calibration.axis, undo=False)
    record = {"axis": calibration.axis.to_dict(), "wavelength_A": float(dataset.wavelength_A), "n_peaks": n,
              "standard": calibration.source_name or calibration.standard, "outside_range": [float(v) for v in outside]}
    dataset.metadata[APPLIED_KEY] = record
    return record


def reapply_calibration(dataset) -> Optional[dict]:
    """Peaks already carry a correction (possibly an older one): undo it, then apply the current one.

    Idempotent: calling it twice leaves the peak positions unchanged.
    """
    previous = dataset.metadata.pop(APPLIED_KEY, None)
    if previous:
        # d is always recomputed with the CURRENT wavelength (it may have changed since the correction was applied)
        _shift_peaks(dataset.peaks, dataset.wavelength_A, AxisCorrection.from_dict(previous["axis"]), undo=True)
    dataset.metadata.pop("calibration_skipped", None)
    calibration, reason = calibration_status(dataset)
    if calibration is None or reason:
        if reason and reason not in ("brak kalibracji", "to jest sam wzorzec"):
            dataset.metadata["calibration_skipped"] = reason
        return None
    return _apply(dataset, calibration)


def install_calibration(datasets: Sequence, standard, calibration: StandardCalibration) -> Dict[str, List[str]]:
    """Store the calibration in the standard (role) and in every other dataset. Returns {"applied": [...], "skipped": [...]}."""
    standard.metadata[ROLE_KEY] = "standard"
    standard.metadata_sources[ROLE_KEY] = "user"
    standard.metadata[CAL_KEY] = calibration.compact()
    standard.metadata_sources[CAL_KEY] = "calibration from standard"
    applied, skipped = [], []
    for ds in datasets:
        if ds is standard:
            continue
        ds.metadata.pop(ROLE_KEY, None)
        ds.metadata[CAL_KEY] = calibration.compact()
        ds.metadata_sources[CAL_KEY] = "calibration from standard"
        _calibration, reason = calibration_status(ds)
        ds.metadata.pop("calibration_skipped", None)
        if reason:
            ds.metadata["calibration_skipped"] = reason
        (skipped if reason else applied).append(ds.name if not reason else f"{ds.name}: {reason}")
    return {"applied": applied, "skipped": skipped}


def remove_calibration(datasets: Sequence) -> int:
    count = 0
    for ds in datasets:
        previous = ds.metadata.pop(APPLIED_KEY, None)
        if previous and ds.peaks:
            _shift_peaks(ds.peaks, ds.wavelength_A, AxisCorrection.from_dict(previous["axis"]), undo=True)
        for key in (CAL_KEY, ROLE_KEY, "calibration_skipped"):
            if ds.metadata.pop(key, None) is not None:
                count += 1
            ds.metadata_sources.pop(key, None)
    return count


def calibration_summary(dataset) -> str:
    """One line for the sample summary panel."""
    md = dataset.metadata
    if md.get(ROLE_KEY) == "standard":
        calibration, _ = calibration_status(dataset)
        return f"Wzorzec kalibracyjny Si: {calibration.axis.describe()}" if calibration else "Wzorzec kalibracyjny Si"
    applied = md.get(APPLIED_KEY)
    if applied:
        axis = AxisCorrection.from_dict(applied["axis"])
        text = f"2θ skorygowane wzorcem „{applied.get('standard', 'Si')}”: {axis.describe()}"
        if applied.get("outside_range"):
            text += f"; UWAGA: {len(applied['outside_range'])} pików poza zakresem wzorca"
        return text
    if md.get("calibration_skipped"):
        return f"Kalibracja wzorcem Si NIE została zastosowana: {md['calibration_skipped']}"
    if md.get(CAL_KEY):
        return "Kalibracja wzorcem Si zapisana; wykonaj analizę, aby skorygować położenia pików."
    return ""
