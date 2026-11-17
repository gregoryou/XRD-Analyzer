from __future__ import annotations

import sys
import os
import csv
import json
import math
import re
import sqlite3
import traceback
import importlib
import io
import zipfile
import uuid
import tempfile
from copy import deepcopy
from datetime import datetime
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Tuple, Dict, Optional

import numpy as np
from scipy.signal import find_peaks, peak_widths
from scipy import sparse
from scipy.sparse.linalg import spsolve

from PyQt6.QtCore import Qt, QSettings, QSignalBlocker, QThread, pyqtSignal
from PyQt6.QtGui import QColor, QAction, QKeySequence
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QFileDialog, QMessageBox, QVBoxLayout,
    QHBoxLayout, QFormLayout, QLineEdit, QPushButton, QLabel, QComboBox,
    QDoubleSpinBox, QSpinBox, QCheckBox, QTableWidget, QTableWidgetItem,
    QSplitter, QGroupBox, QHeaderView, QDialog, QDialogButtonBox, QTextEdit,
    QInputDialog, QAbstractItemView, QTabWidget, QProgressDialog, QColorDialog, QToolButton,
    QGridLayout, QFrame, QScrollArea, QSizePolicy, QAbstractSpinBox
)

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.transforms import blended_transform_factory
from xrd_rd_import import parse_philips_rd
from xrd_kalpha import dataset_kalpha2_candidates, traceability_rows, verified_wavelength_pair, cod_assignment_for_peak
from xrd_cod_offline import (
    default_cod_path, parse_elements, parse_ids, sync_api, import_zip_archive,
    download_full_cod_zip, CODSyncError,
)
from xrd_detection_support import find_candidate_peaks, noise_sigma_from_differences
from xrd_calibration import apply_calibration_to_fresh_peaks, calibration_summary, reapply_calibration
from xrd_plot_support import (
    display_smoothing_settings, smooth_for_display, sanitize_plot_state,
    SMOOTH_METHODS, SMOOTH_METHOD_LABELS,
    valid_color, optional_bool, safe_window_rect,
)

try:
    from openpyxl import Workbook
    OPENPYXL_OK = True
except Exception:
    OPENPYXL_OK = False

gemmi = None
GEMMI_OK = False
GEMMI_IMPORT_ERROR = ""

def ensure_gemmi() -> bool:
    """Lazy-load gemmi so an installation made after the application was first built is detected.

    The diagnostic deliberately reports sys.executable because on Windows it is common to install
    gemmi into a different Python installation than the one used to launch this GUI.
    """
    global gemmi, GEMMI_OK, GEMMI_IMPORT_ERROR
    if GEMMI_OK and gemmi is not None:
        return True
    try:
        gemmi = importlib.import_module("gemmi")
        GEMMI_OK = True
        GEMMI_IMPORT_ERROR = ""
        return True
    except Exception as exc:
        gemmi = None
        GEMMI_OK = False
        GEMMI_IMPORT_ERROR = f"{type(exc).__name__}: {exc}"
        return False

ensure_gemmi()

APP_NAME = "XRD Analyzer"
APP_VERSION = "20.3.5.1"
DEFAULT_WAVELENGTH_A = 1.5406  # Cu Kα, Å

DEFAULT_SEPARATOR_COLOR = "#808080"

def normalizuj_separator(item: dict, liczba_probek: int) -> dict:
    """Zwróć separator w aktualnym schemacie niezależnie od wersji projektu.

    Migracja jest wykonywana wyłącznie w pamięci. Stary pojedynczy klucz ``label``
    staje się etykietą NAD separatorem, a etykieta pod separatorem pozostaje pusta.
    Oryginalny plik projektu nie jest modyfikowany do czasu kolejnego zapisu.
    """
    item = dict(item or {})
    try:
        position = int(item.get("position", 0))
    except (TypeError, ValueError):
        position = 0
    try:
        offset_pct = float(item.get("offset_pct", 0.0))
    except (TypeError, ValueError):
        offset_pct = 0.0

    top_label = item.get("top_label", None)
    if top_label is None:
        top_label = item.get("label", "")

    return {
        "position": max(0, min(int(max(0, liczba_probek)), position)),
        "top_label": str(top_label or ""),
        "bottom_label": str(item.get("bottom_label", "") or ""),
        "offset_pct": max(-48.0, min(48.0, offset_pct)),
        "color": str(item.get("color", DEFAULT_SEPARATOR_COLOR) or DEFAULT_SEPARATOR_COLOR),
        "labels_visible": bool(item.get("labels_visible", True)),
    }


@dataclass
class PeakResult:
    index: int
    two_theta: float
    d_A: float
    height: float
    rel_intensity: float
    fwhm_deg: float
    phase: str = "unassigned"
    hkl: str = ""
    delta_2theta: float = math.nan
    custom_label: str = ""
    prominence: float = math.nan
    snr: float = math.nan
    phase_source: str = ""  # "internal", "cod", or empty/unassigned

    @property
    def width_half_prominence_deg(self) -> float:
        """Detection width; legacy fwhm_deg is kept only for file/API compatibility."""
        return self.fwhm_deg


def angular_distance_selection(indices, x, strength, min_distance_deg):
    """Return positions in candidate arrays after suppression in actual degrees.

    Highest prominence wins; ties prefer the lower angle deterministically.
    The comparison tolerance only covers floating point round-off at the boundary.
    """
    indices = np.asarray(indices, dtype=int)
    strength = np.asarray(strength, dtype=float)
    if not math.isfinite(min_distance_deg) or min_distance_deg < 0:
        raise ValueError("Minimalny odstęp kątowy musi być skończony i nieujemny.")
    order = np.lexsort((indices, -strength))
    selected = []
    for pos in order:
        angle = float(x[indices[pos]])
        if all(abs(angle - float(x[indices[j]])) + 1e-12 >= min_distance_deg for j in selected):
            selected.append(int(pos))
    return np.asarray(sorted(selected, key=lambda pos: indices[pos]), dtype=int)


@dataclass
class Dataset:
    path: Path
    name: str
    material_type: str
    x: np.ndarray
    y: np.ndarray
    wavelength_A: float = DEFAULT_WAVELENGTH_A
    corrected_y: Optional[np.ndarray] = None
    baseline: Optional[np.ndarray] = None
    peaks: List[PeakResult] = field(default_factory=list)
    # Instrument/measurement metadata. Values may come from the source file,
    # an instrument profile, or manual user input. metadata_sources stores
    # that provenance explicitly (e.g. "file", "profile: X'Pert", "user").
    metadata: Dict[str, object] = field(default_factory=dict)
    metadata_sources: Dict[str, str] = field(default_factory=dict)
    # Project identity survives renaming, relocation and missing source files.
    dataset_id: str = field(default_factory=lambda: uuid.uuid4().hex)

    @property
    def max_total_intensity(self) -> float:
        return float(np.nanmax(self.y)) if len(self.y) else float("nan")

    @property
    def scan_start_deg(self) -> float:
        return float(np.nanmin(self.x)) if len(self.x) else float("nan")

    @property
    def scan_end_deg(self) -> float:
        return float(np.nanmax(self.x)) if len(self.x) else float("nan")

    @property
    def scan_step_deg(self) -> float:
        if len(self.x) < 2:
            return float("nan")
        return float(np.nanmedian(np.diff(self.x)))


DEFAULT_INSTRUMENT_PROFILE = {
    "name": "PANalytical PW3040/60 X'Pert Pro – Cu Kα (35 kV / 40 mA)",
    "instrument": "PANalytical PW3040/60 X'Pert Pro",
    "manufacturer": "Malvern Panalytical, Malvern, UK",
    "radiation": "Cu Kα",
    "anode": "Cu",
    "wavelength_A": 1.5406,
    "tube_voltage_kV": 35.0,
    "tube_current_mA": 40.0,
    "scan_rate_deg_s": 0.031,
    "geometry": "θ–2θ",
    "notes": "User-defined laboratory profile. Scan range and step are read from each data file whenever available.",
}


class InstrumentProfileStore:
    """Persistent, editable instrument profiles stored as JSON in XRDAnalyzerData."""
    def __init__(self):
        data_dir = Path.home() / "XRDAnalyzerData"
        try:
            data_dir.mkdir(parents=True, exist_ok=True)
        except OSError:
            # A read-only profile directory must not prevent the program from starting.
            data_dir = Path(tempfile.gettempdir()) / "XRDAnalyzerData"
            data_dir.mkdir(parents=True, exist_ok=True)
        self.path = data_dir / "instrument_profiles.json"
        self.profiles: Dict[str, dict] = {}
        self._load()

    def _load(self):
        try:
            if self.path.exists():
                raw = json.loads(self.path.read_text(encoding="utf-8"))
                if isinstance(raw, list):
                    for rec in raw:
                        if isinstance(rec, dict) and str(rec.get("name", "")).strip():
                            self.profiles[str(rec["name"]).strip()] = dict(rec)
                elif isinstance(raw, dict):
                    for name, rec in raw.items():
                        if isinstance(rec, dict):
                            r = dict(rec); r["name"] = str(r.get("name", name))
                            self.profiles[r["name"]] = r
        except Exception:
            self.profiles = {}
        if not self.profiles:
            self.profiles[DEFAULT_INSTRUMENT_PROFILE["name"]] = dict(DEFAULT_INSTRUMENT_PROFILE)
            try:
                self.save()
            except OSError:
                pass

    def save(self):
        records = [self.profiles[k] for k in sorted(self.profiles, key=str.lower)]
        self.path.write_text(json.dumps(records, ensure_ascii=False, indent=2), encoding="utf-8")

    def names(self) -> List[str]:
        return sorted(self.profiles.keys(), key=str.lower)

    def get(self, name: str) -> Optional[dict]:
        rec = self.profiles.get(str(name))
        return dict(rec) if rec else None

    def upsert(self, profile: dict):
        name = str(profile.get("name", "")).strip()
        if not name:
            raise ValueError("Nazwa profilu aparatury nie może być pusta.")
        rec = dict(profile); rec["name"] = name
        self.profiles[name] = rec
        self.save()

    def delete(self, name: str):
        if name in self.profiles:
            del self.profiles[name]
            if not self.profiles:
                self.profiles[DEFAULT_INSTRUMENT_PROFILE["name"]] = dict(DEFAULT_INSTRUMENT_PROFILE)
            self.save()


def _meta_missing(ds: Dataset, key: str) -> bool:
    if key == "wavelength_A":
        src = str(ds.metadata_sources.get("wavelength_A", "default")).lower()
        return src in ("", "default", "unknown", "unverified")
    value = ds.metadata.get(key, None)
    return value is None or (isinstance(value, str) and not value.strip())


def set_dataset_wavelength(ds: Dataset, wavelength_A: float, source: str) -> bool:
    """Recalculate measured d-spacings and invalidate results that used the old lambda.

    Return whether an existing identification was invalidated. Measured peak
    positions, intensities, widths, and user peak text remain unchanged.
    """
    wavelength_A = float(wavelength_A)
    if not math.isfinite(wavelength_A) or wavelength_A <= 0:
        raise ValueError("Długość fali musi być skończoną liczbą większą od zera.")
    if wavelength_A == ds.wavelength_A:
        ds.metadata["wavelength_A"] = wavelength_A
        ds.metadata_sources["wavelength_A"] = source
        return False

    # Compute before modifying the dataset so invalid input cannot leave mixed results.
    spacings = []
    for peak in ds.peaks:
        if not math.isfinite(peak.two_theta) or not 0 < peak.two_theta <= 180:
            raise ValueError("Nie można przeliczyć d: pozycja piku musi należeć do zakresu 0 < 2θ ≤ 180°.")
        spacings.append(bragg_d(peak.two_theta, wavelength_A))
    dependent_keys = (
        "phase_fits", "structural_family_fits", "manual_cod_peak_assignments",
        "ai_phase_ranking", "ai_user_confirmed_candidate", "publication_phase_aliases",
        "scientific_analysis",
        "cod_auto_candidate_count", "ai_phase_model", "ai_phase_warning", "cod_search_warning",
    )
    invalidated = any(
        p.phase not in ("", "unassigned") or p.hkl or p.phase_source
        for p in ds.peaks
    ) or any(bool(ds.metadata.get(key)) for key in dependent_keys[:6])
    ds.wavelength_A = wavelength_A
    ds.metadata["wavelength_A"] = wavelength_A
    ds.metadata_sources["wavelength_A"] = source
    for peak, spacing in zip(ds.peaks, spacings):
        peak.d_A = spacing
        peak.phase = "unassigned"
        peak.hkl = ""
        peak.delta_2theta = math.nan
        peak.phase_source = ""
    for key in dependent_keys:
        ds.metadata.pop(key, None)
        ds.metadata_sources.pop(key, None)
    if invalidated:
        ds.metadata["phase_identification_source"] = "wymaga ponownej identyfikacji po zmianie λ"
        ds.metadata_sources["phase_identification_source"] = "zmiana długości fali"
    if ds.metadata.get("calibration"):
        reapply_calibration(ds)      # a calibration made for another wavelength is no longer valid for this scan
    return bool(invalidated)


def apply_instrument_profile(ds: Dataset, profile: dict, fill_missing_only: bool = True) -> bool:
    """Apply an instrument profile while preserving file metadata by default."""
    if not profile:
        return False
    pname = str(profile.get("name", "Profil aparatury")).strip()
    wave = profile.get("wavelength_A", None)
    invalidated = False
    if wave not in (None, "") and ((not fill_missing_only) or _meta_missing(ds, "wavelength_A")):
        invalidated = set_dataset_wavelength(ds, wave, f"profile: {pname}")
    keys = (
        "instrument", "manufacturer", "radiation", "anode",
        "tube_voltage_kV", "tube_current_mA", "scan_rate_deg_s",
        "geometry", "notes",
    )
    for key in keys:
        value = profile.get(key, None)
        if value is None or (isinstance(value, str) and not value.strip()):
            continue
        if (not fill_missing_only) or _meta_missing(ds, key):
            ds.metadata[key] = value
            ds.metadata_sources[key] = f"profile: {pname}"
    ds.metadata["instrument_profile"] = pname
    ds.metadata_sources["instrument_profile"] = "profile"
    return invalidated


def metadata_source(ds: Dataset, key: str) -> str:
    return str(ds.metadata_sources.get(key, "unknown"))


# --------------------------- phase database ---------------------------

def cubic_d(a: float, h: int, k: int, l: int) -> float:
    return a / math.sqrt(h*h + k*k + l*l)


def hex_d(a: float, c: float, h: int, k: int, l: int) -> float:
    inv_d2 = (4.0/3.0) * (h*h + h*k + k*k) / (a*a) + (l*l)/(c*c)
    return 1.0 / math.sqrt(inv_d2)


def tetragonal_d(a: float, c: float, h: int, k: int, l: int) -> float:
    inv_d2 = (h*h + k*k)/(a*a) + (l*l)/(c*c)
    return 1.0 / math.sqrt(inv_d2)


def formula_elements(formula: str) -> set[str]:
    """Extract element symbols from formulas such as Fe-Cr-Ni, Ni3Al or Ti-6Al-4V."""
    if not formula:
        return set()
    return set(re.findall(r"[A-Z][a-z]?", str(formula)))


# Metallic elements used only to build concise, human-readable COD labels.
# This does not change phase matching or chemistry filtering.
COD_METALLIC_ELEMENTS = set(
    "Li Be Na Mg Al K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn "
    "Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Th Pa U".split()
)


def _cod_clean_formula(formula: str) -> str:
    s = str(formula or "").strip()
    s = re.sub(r"^[-\s]+|[-\s]+$", "", s)
    return re.sub(r"\s+", " ", s).strip()


def cod_structure_short(space_group: str = "", sg_no=None) -> str:
    """Return a familiar structure label when COD space-group data justify it."""
    try:
        n = int(float(sg_no))
    except Exception:
        n = 0
    sg = re.sub(r"\s+", "", str(space_group or "")).lower()
    if n == 225 or sg in {"fm-3m", "fm3m"}:
        return "FCC"
    if n == 229 or sg in {"im-3m", "im3m"}:
        return "BCC"
    if n == 194 or sg in {"p63/mmc", "p6_3/mmc"}:
        return "HCP"
    if n:
        if 195 <= n <= 230:
            return "cubic"
        if 168 <= n <= 194:
            return "hexagonal"
        if 143 <= n <= 167:
            return "trigonal"
        if 75 <= n <= 142:
            return "tetragonal"
        if 16 <= n <= 74:
            return "orthorhombic"
        if 3 <= n <= 15:
            return "monoclinic"
        if 1 <= n <= 2:
            return "triclinic"
    return str(space_group or "").strip()


def cod_friendly_phase_names(formula: str, chemname: str = "", mineral: str = "", commonname: str = "",
                             space_group: str = "", sg_no=None):
    """Create publication-friendly and compact labels from COD metadata.

    COD often has no alloy trade/phase label such as 'Cu-Ni FCC'.  For simple metallic
    entries the formula and space group are sufficient to generate a neutral descriptive
    label.  The COD ID remains available in tables/details but is intentionally omitted
    from plot labels.
    """
    formula_c = _cod_clean_formula(formula)
    elems_ordered = []
    for e in re.findall(r"[A-Z][a-z]?", formula_c):
        if e not in elems_ordered:
            elems_ordered.append(e)
    struct = cod_structure_short(space_group, sg_no)
    simple_metal = bool(elems_ordered) and len(elems_ordered) <= 5 and all(e in COD_METALLIC_ELEMENTS for e in elems_ordered)

    # For FCC/BCC/HCP metallic entries, this is much clearer than a raw COD number.
    if simple_metal and struct in {"FCC", "BCC", "HCP"}:
        if len(elems_ordered) == 1:
            base = str(commonname or mineral or chemname or elems_ordered[0]).strip()
            base = base if base and base not in {"\\N", "?", "."} else elems_ordered[0]
            return f"{base} ({struct})", f"{base} {struct}"
        base = "–".join(elems_ordered)
        return f"{base} solid solution ({struct})", f"{base} {struct}"

    cod_name = str(commonname or mineral or chemname or "").strip()
    if cod_name in {"", "\\N", "?", "."}:
        cod_name = formula_c or "COD structure"
    long_name = f"{cod_name} ({struct})" if struct and struct.lower() not in cod_name.lower() else cod_name
    # Plot labels should stay short even for long systematic chemical names.
    compact_formula = re.sub(r"\s+", "", formula_c) or cod_name
    short_name = f"{compact_formula} {struct}".strip() if struct else compact_formula
    if len(short_name) > 30:
        short_name = short_name[:29] + "…"
    return long_name, short_name


PEAK_LABEL_STYLE_SHORT = "Publikacyjny – krótki"
PEAK_LABEL_STYLE_FULL = "Publikacyjny – pełny"
PEAK_LABEL_STYLE_DEBUG = "Baza danych / diagnostyka"
PEAK_LABEL_STYLES = [PEAK_LABEL_STYLE_SHORT, PEAK_LABEL_STYLE_FULL, PEAK_LABEL_STYLE_DEBUG]

HKL_DISPLAY_AUTO = "Pokaż hkl – z przypisanego źródła (zalecane)"
HKL_DISPLAY_OFF = "Ukryj hkl"
HKL_DISPLAY_INTERNAL = "Pokaż hkl – tylko baza wewnętrzna"
HKL_DISPLAY_COD = "Pokaż hkl – tylko COD"
HKL_DISPLAY_MODES = [HKL_DISPLAY_AUTO, HKL_DISPLAY_OFF, HKL_DISPLAY_INTERNAL, HKL_DISPLAY_COD]


def _migrate_peak_label_style(text: str) -> str:
    return {
        "Publication – short": PEAK_LABEL_STYLE_SHORT,
        "Publication – full": PEAK_LABEL_STYLE_FULL,
        "Database / debug": PEAK_LABEL_STYLE_DEBUG,
    }.get(str(text), str(text))

def _migrate_hkl_display(text: str) -> str:
    return {
        "Show hkl – assigned source (recommended)": HKL_DISPLAY_AUTO,
        "Hide hkl": HKL_DISPLAY_OFF,
        "Show hkl – Internal database only": HKL_DISPLAY_INTERNAL,
        "Show hkl – COD only": HKL_DISPLAY_COD,
    }.get(str(text), str(text))

def _migrate_phase_source(text: str) -> str:
    return {
        "Internal database only": "Tylko baza wewnętrzna",
        "External COD only": "Tylko zewnętrzna baza COD",
        "Internal + external COD": "Baza wewnętrzna + zewnętrzna COD",
        "No automatic phase identification": "Bez automatycznej identyfikacji faz",
    }.get(str(text), str(text))

def _migrate_label_mode(text: str) -> str:
    return {
        "Matched phases only": "Tylko dopasowane fazy",
        "Top N strongest peaks": "Top N najsilniejszych pików",
        "All detected peaks": "Wszystkie wykryte piki",
    }.get(str(text), str(text))


def _clean_hkl_label(hkl: str) -> str:
    h = str(hkl or "").split(" / ")[0].strip()
    if not h:
        return ""
    h = h.strip("()[]{} ")
    return f"({h})" if h else ""


def _publication_phase_names(phase: str, material_type: str = "") -> tuple[str, str]:
    """Return (short, full) publication labels for a phase name.

    This is deliberately descriptive rather than a chemical-identification claim.
    It normalizes common crystallographic names returned by the internal library or COD.
    """
    raw = str(phase or "").split(" / ")[0].strip()
    low = raw.lower().replace("–", "-")
    mat = str(material_type or "").lower()
    if not raw or low == "unassigned":
        return "", ""
    if "austen" in low or "gamma" in low or "γ" in raw:
        return "γ", "γ-austenite FCC"
    if "delta" in low or "δ" in raw:
        return "δ", "δ-ferrite BCC"
    if "ferrite" in low or "ferryt" in low:
        # Stainless-steel weld/base-metal ferrite is commonly δ; generic steels use α.
        if any(k in mat for k in ("316", "304", "stainless", "1.4404", "1.4301", "x2cr", "x5cr")):
            return "δ", "δ-ferrite BCC"
        return "α-Fe", "α-ferrite BCC"
    if ("cu-ni" in low or "cu–ni" in raw.lower() or ("cu" in low and "ni" in low)):
        return "Cu–Ni", "Cu–Ni FCC"
    if ("ti" in low or "titan" in low) and ("hcp" in low or "alpha" in low or "α" in raw):
        return "α-Ti", "α-Ti HCP"
    if ("ti" in low or "titan" in low) and ("bcc" in low or "beta" in low or "β" in raw):
        return "β-Ti", "β-Ti BCC"
    if "aluminium" in low or "aluminum" in low or raw.strip() == "Al" or low.startswith("al "):
        return "Al", "Al FCC"
    if raw.strip() == "Cu" or low.startswith("copper"):
        return "Cu", "Cu FCC"
    if raw.strip() == "Ni" or low.startswith("nickel"):
        return "Ni", "Ni FCC"
    if raw.strip() == "Fe" or low.startswith("iron"):
        return "Fe", "Fe"
    # Remove COD wording and verbose solid-solution wording while retaining structure type.
    compact = re.sub(r"\bsolid solution\b", "", raw, flags=re.I)
    compact = re.sub(r"\s+", " ", compact).strip(" -–")
    short = compact
    for token in (" (FCC)", " (BCC)", " (HCP)"):
        short = short.replace(token, "")
    short = short.replace(" FCC", "").replace(" BCC", "").replace(" HCP", "").strip() or compact
    if len(short) > 22:
        short = short[:21] + "…"
    full = compact
    return short, full


def _peak_cod_id(dataset: Dataset, peak: PeakResult):
    fits = dataset.metadata.get("phase_fits", []) if isinstance(dataset.metadata, dict) else []
    if not isinstance(fits, list):
        return None
    phase = str(peak.phase or "").strip().lower()
    for rec in fits:
        if not isinstance(rec, dict):
            continue
        names = [rec.get("plot_phase"), rec.get("phase"), rec.get("name"), rec.get("formula")]
        names = [str(x).strip().lower() for x in names if x not in (None, "")]
        if phase and phase in names and rec.get("cod_id") not in (None, ""):
            try:
                return int(rec.get("cod_id"))
            except Exception:
                return rec.get("cod_id")
    # If there is only one accepted COD phase, it is a useful diagnostic fallback.
    cod_fits = [r for r in fits if isinstance(r, dict) and r.get("cod_id") not in (None, "")]
    if len(cod_fits) == 1:
        try:
            return int(cod_fits[0].get("cod_id"))
        except Exception:
            return cod_fits[0].get("cod_id")
    return None


def _hkl_for_display(peak: PeakResult, hkl_mode: str) -> str:
    """Return a formatted Miller-index label according to the selected provenance filter."""
    mode = str(hkl_mode or HKL_DISPLAY_AUTO)
    source = str(getattr(peak, "phase_source", "") or "").strip().lower()
    if mode == HKL_DISPLAY_OFF:
        return ""
    if mode == HKL_DISPLAY_INTERNAL and not source.startswith("internal"):
        return ""
    if mode == HKL_DISPLAY_COD and source != "cod":
        return ""
    return _clean_hkl_label(peak.hkl)


def format_peak_label(dataset: Dataset, peak: PeakResult, style: str, hkl_mode: str = HKL_DISPLAY_AUTO) -> str:
    """Format one peak label for publication or database/debug display.

    Miller indices are optional. In combined COD + Internal mode, provenance stored on
    each assigned peak decides whether the hkl belongs to Internal or COD.
    User-confirmed publication aliases from the AI-assisted review are presentation-only:
    they never change the underlying COD/Internal assignment stored on the peak.
    """
    if str(peak.custom_label or "").strip():
        return str(peak.custom_label).strip()
    hkl = _hkl_for_display(peak, hkl_mode)
    if hkl and str(getattr(peak,'phase_source','')).lower() == 'cod':
        phase_name = str(peak.phase or '').lower()
        cod_rec = cod_assignment_for_peak(dataset,peak)
        space_group = str(cod_rec.get('space_group','')).replace(' ','').lower()
        is_cubic = ('cubic' in phase_name or '-3' in space_group or 'm3' in space_group)
        if is_cubic:
            hkl = _clean_hkl_label(canonical_hkl_for_structure(peak.hkl,'cubic'))
    phase_for_label = str(peak.phase or "")
    aliases = dataset.metadata.get("publication_phase_aliases", {}) if isinstance(dataset.metadata, dict) else {}
    if isinstance(aliases, dict) and phase_for_label:
        alias = aliases.get(phase_for_label)
        if not alias:
            low = phase_for_label.strip().lower()
            for k, v in aliases.items():
                if str(k).strip().lower() == low:
                    alias = v; break
        if alias:
            phase_for_label = str(alias)
    short, full = _publication_phase_names(phase_for_label, dataset.material_type)
    if str(getattr(peak, "phase_source", "") or "") == TENTATIVE_SOURCE and short:
        short, full = short + "?", f"{full} ({TENTATIVE_SUFFIX})"
    style = str(style or PEAK_LABEL_STYLE_FULL)
    if style == PEAK_LABEL_STYLE_SHORT:
        if short and hkl:
            return f"{short}{hkl}"
        return short or hkl or ""
    if style == PEAK_LABEL_STYLE_DEBUG:
        base = full or str(peak.phase or "unassigned")
        cid = _peak_cod_id(dataset, peak)
        parts = [base + (f" {hkl}" if hkl else "")]
        if cid is not None:
            parts.append(f"COD {cid}")
        source = str(getattr(peak, "phase_source", "") or "").strip()
        if source:
            parts.append(f"source={source}")
        parts.append(f"2θ={peak.two_theta:.2f}°")
        return "\n".join(parts)
    if full and hkl:
        return f"{full} {hkl}"
    return full or hkl or ""


def default_elements_for_material(material_text: str, sample_name: str = "") -> str:
    """Suggest a conservative chemistry filter for automatic phase/COD search.

    The preset intentionally contains the principal alloying elements only. Trace elements
    (O, N, C, Pb, P, S, etc.) can be added manually when oxide/nitride/carbide or inclusion
    searches are required.  A narrower COD chemistry filter greatly reduces false matches.
    """
    t = f"{material_text} {sample_name}".lower()
    is_x5_1810 = ("x5crni18-10" in t or "1.4301" in t or "aisi 304" in t or "aisi304" in t)
    is_316 = ("316" in t or "x2crnimo17-12-2" in t or "1.4404" in t)
    is_stainless = ("stainless" in t or "nierdzew" in t)
    is_cuni = ("cuni" in t or "monel" in t or "ercuni" in t or "cunifer" in t)
    is_weld = ("weld" in t or "spoina" in t or "joint" in t)

    # Dissimilar 316L/CuNi weld metal: principal elements contributed by both sheets/fillers.
    if is_weld and (is_316 or "316l" in t) and is_cuni:
        return "Fe Cr Ni Mo Mn Si Cu"
    if is_x5_1810:
        return "Fe Cr Ni Mn Si"
    if is_316:
        return "Fe Cr Ni Mo Mn Si"
    if is_stainless:
        return "Fe Cr Ni Mn Si Mo"
    if is_cuni or "copper alloy" in t:
        return "Cu Ni Fe Mn"
    if "brass" in t:
        return "Cu Zn"
    if "bronze" in t:
        return "Cu Sn"
    if "ti6al4v" in t or "ti-6al-4v" in t or "titanium alloy" in t or re.search(r"(^|[^a-z])ti([^a-z]|$)", t):
        return "Ti Al V Fe"
    if "aluminium" in t or "aluminum" in t or re.search(r"(^|[^a-z])al([^a-z]|$)", t):
        return "Al Mg Si Cu Zn Mn Fe Cr Ti Sc Zr"
    if "zircaloy" in t or "zirconium" in t or re.search(r"(^|[^a-z])zr([^a-z]|$)", t):
        return "Zr Sn Fe Cr Nb"
    if "magnesium" in t or "az31" in t or "az91" in t or re.search(r"(^|[^a-z])mg([^a-z]|$)", t):
        return "Mg Al Zn Mn Si Cu Fe"
    if "steel" in t or "stal" in t:
        return "Fe C Mn Si Cr Ni Mo V W Nb Ti Al N B"
    return ""


def _cubic_family_hkls(lattice: str = "P", max_index: int = 10, d_min_A: float = 0.60,
                       a_A: float = 4.0) -> List[Tuple[int,int,int]]:
    """Generate unique cubic reflection families h>=k>=l with basic extinction rules."""
    lattice = lattice.upper()
    out = []
    for h in range(0, max_index + 1):
        for k in range(0, h + 1):
            for l in range(0, k + 1):
                if h == k == l == 0:
                    continue
                allowed = True
                if lattice == "FCC":
                    allowed = ((h % 2 == k % 2 == l % 2))
                elif lattice == "BCC":
                    allowed = ((h + k + l) % 2 == 0)
                elif lattice == "DIAMOND":
                    # Diamond: FCC condition plus extinction of all-even reflections
                    # for which h+k+l is not divisible by 4.
                    fcc = (h % 2 == k % 2 == l % 2)
                    if not fcc:
                        allowed = False
                    elif h % 2 == 0:  # all even
                        allowed = ((h + k + l) % 4 == 0)
                elif lattice == "P":
                    allowed = True
                if not allowed:
                    continue
                d = cubic_d(a_A, h, k, l)
                if d >= d_min_A:
                    out.append((h,k,l))
    out.sort(key=lambda q: (q[0]*q[0]+q[1]*q[1]+q[2]*q[2], q), reverse=False)
    return out


def _tetragonal_family_hkls(a_A: float, c_A: float, max_index: int = 9,
                            d_min_A: float = 0.60, centering: str = "P") -> List[Tuple[int,int,int]]:
    out = []
    centering = centering.upper()
    for h in range(0, max_index + 1):
        for k in range(0, h + 1):
            for l in range(0, max_index + 1):
                if h == k == l == 0:
                    continue
                if centering == "I" and ((h+k+l) % 2 != 0):
                    continue
                d = tetragonal_d(a_A, c_A, h, k, l)
                if d >= d_min_A:
                    out.append((h,k,l))
    # Collapse near-identical d values while keeping one representative hkl.
    unique = []
    seen = []
    for hkl in sorted(out, key=lambda q: tetragonal_d(a_A,c_A,*q), reverse=True):
        d = tetragonal_d(a_A,c_A,*hkl)
        if any(abs(d-x) < 1e-5 for x in seen):
            continue
        seen.append(d); unique.append(hkl)
    return unique


def seed_phase_definitions():
    """
    Engineering reference library for common metallic materials.

    For simple substitutional alloys the alloy name is a material preset/reference,
    while XRD identifies the underlying crystallographic phase (e.g. alpha-Fe BCC,
    gamma-austenite FCC, alpha-Al FCC, alpha-Ti HCP). Lattice parameters are
    representative and composition/temperature dependent. COD search remains the
    preferred route for uncommon intermetallics and precise crystallographic work.
    """
    hcp_hkls = [
        (1,0,0),(0,0,2),(1,0,1),(1,0,2),(1,1,0),(1,0,3),(2,0,0),(1,1,2),
        (2,0,1),(0,0,4),(2,0,2),(1,0,4),(2,0,3),(2,1,0),(2,1,1),(1,1,4),
        (2,1,2),(1,0,5),(2,0,4),(3,0,0),(2,1,3),(0,0,6),(3,0,2),(2,2,0),
        (3,1,0),(2,1,4),(3,0,3),(2,2,2),(3,1,1),(3,0,4)
    ]

    phases = []

    def add_cubic(name, formula, structure, space_group, a, category, notes="", hkls=None, lattice=None, fit_lattice=None):
        if lattice is None:
            su = structure.upper()
            if "DIAMOND" in su:
                lattice = "DIAMOND"
            elif "FCC" in su:
                lattice = "FCC"
            elif "BCC" in su:
                lattice = "BCC"
            elif "B2" in su or "L1" in su or "PRIMITIVE" in su:
                lattice = "P"
            else:
                lattice = "P"
        hkls = hkls or _cubic_family_hkls(lattice, max_index=10, d_min_A=0.60, a_A=a)
        peaks = [(cubic_d(a,*hkl), f"{hkl[0]}{hkl[1]}{hkl[2]}") for hkl in hkls]
        if fit_lattice is None:
            fit_lattice = (lattice in {"FCC", "BCC"} and ("FCC" in structure.upper() or "BCC" in structure.upper()))
        model = f"CUBIC_{lattice}" if fit_lattice and lattice in {"FCC", "BCC"} else "FIXED"
        note2 = (notes + " " if notes else "") + (
            "Reference lattice parameter is an initial search seed only; XRD Analyzer fits the cubic lattice parameter to the measured reflections."
            if fit_lattice else ""
        )
        phases.append((name, formula, structure, space_group, category, note2.strip(), peaks, model, float(a), None, int(bool(fit_lattice))))

    def add_hex(name, formula, structure, space_group, a, c, category, notes="", hkls=None, fit_lattice=None):
        hs = hkls or hcp_hkls
        peaks = []
        seen = []
        for hkl in hs:
            d = hex_d(a,c,*hkl)
            if d < 0.60 or any(abs(d-x) < 1e-5 for x in seen):
                continue
            seen.append(d)
            peaks.append((d, f"{hkl[0]}{hkl[1]}{hkl[2]}"))
        if fit_lattice is None:
            fit_lattice = (structure.upper().strip() == "HCP")
        model = "HCP" if fit_lattice else "FIXED"
        note2 = (notes + " " if notes else "") + (
            "Reference a and c are initial search seeds only; XRD Analyzer fits the HCP lattice parameters to the measured reflections."
            if fit_lattice else ""
        )
        phases.append((name, formula, structure, space_group, category, note2.strip(), peaks, model, float(a), float(c), int(bool(fit_lattice))))

    def add_tetragonal(name, formula, structure, space_group, a, c, category, notes="", centering="P"):
        hs = _tetragonal_family_hkls(a, c, max_index=9, d_min_A=0.60, centering=centering)
        peaks = [(tetragonal_d(a,c,*hkl), f"{hkl[0]}{hkl[1]}{hkl[2]}") for hkl in hs]
        phases.append((name, formula, structure, space_group, category, notes, peaks, "FIXED", float(a), float(c), 0))

    # ------------------------------------------------------------------
    # IRON, UNALLOYED / LOW-ALLOY / ALLOYED STEELS
    # ------------------------------------------------------------------
    add_cubic("α-Fe ferrite – unalloyed/carbon steel", "Fe", "BCC", "Im-3m", 2.8665, "steel – unalloyed",
              "Main ferritic matrix reference for unalloyed and carbon steels; alloying and stress shift peak positions.")
    add_cubic("α-Fe ferrite – low-alloy steel", "Fe-C-Mn-Si-Cr-Mo-Ni", "BCC", "Im-3m", 2.868, "steel – low alloy",
              "Representative ferritic matrix of low-alloy steels; use chemistry filtering and COD for exact precipitates.")
    add_cubic("γ-austenite – alloy/stainless steel", "Fe-Cr-Ni", "FCC", "Fm-3m", 3.595, "steel – alloyed",
              "Representative austenitic phase; lattice parameter varies with C, N, Cr, Ni, Mo and residual stress.")
    add_cubic("X5CrNi18-10 / EN 1.4301 – γ austenite", "Fe-Cr-Ni", "FCC", "Fm-3m", 3.586,
              "steel – X5CrNi18-10",
              "Dedicated X5CrNi18-10 austenitic stainless-steel reference. The 3.586 Å seed is consistent with the supplied 0/750/950/1100 °C air-cooled XRD series; the actual a is fitted from each measured pattern.")
    add_cubic("X5CrNi18-10 / EN 1.4301 – δ ferrite", "Fe-Cr-Ni", "BCC", "Im-3m", 2.875,
              "steel – X5CrNi18-10",
              "Possible minor BCC ferritic constituent. Require multiple independent BCC reflections before reporting it.")
    add_cubic("Ferritic stainless steel α (Fe-Cr)", "Fe-Cr", "BCC", "Im-3m", 2.875, "steel – stainless ferritic")
    add_cubic("Duplex stainless steel γ-austenite", "Fe-Cr-Ni-Mo-N", "FCC", "Fm-3m", 3.600, "steel – duplex")
    add_cubic("Duplex stainless steel δ-ferrite", "Fe-Cr-Ni-Mo-N", "BCC", "Im-3m", 2.880, "steel – duplex")

    # Exact grade requested by the user: EN 1.4404 / AISI 316L.
    add_cubic("X2CrNiMo17-12-2 / 1.4404 / AISI 316L – γ austenite", "Fe-Cr-Ni-Mo", "FCC", "Fm-3m", 3.594,
              "steel – X2CrNiMo17-12-2",
              "Dedicated practical reference for X2CrNiMo17-12-2 (EN 1.4404 / AISI 316L). Representative a varies with chemistry and stress.")
    add_cubic("X2CrNiMo17-12-2 / 1.4404 / AISI 316L – δ ferrite", "Fe-Cr-Ni-Mo", "BCC", "Im-3m", 2.876,
              "steel – X2CrNiMo17-12-2",
              "Possible minor ferritic phase in 316L/weld metal; use multiple BCC reflections for confirmation.")
    # Keep backward-compatible names used by earlier projects.
    add_cubic("Austenite γ (Fe-Cr-Ni)", "Fe-Cr-Ni", "FCC", "Fm-3m", 3.595, "steel")
    add_cubic("Ferrite α / δ (Fe-rich)", "Fe", "BCC", "Im-3m", 2.8665, "steel")
    add_cubic("316L austenite (reference)", "Fe-Cr-Ni-Mo", "FCC", "Fm-3m", 3.594, "316L")

    # Martensite: representative BCT; c/a varies strongly with carbon.
    add_tetragonal("Martensite α' – carbon/low-alloy steel", "Fe-C", "BCT", "I4/mmm", 2.866, 2.90,
                   "steel – martensitic", "Representative BCT martensite; exact tetragonality is carbon dependent.", centering="I")
    add_tetragonal("Martensitic stainless steel α'", "Fe-Cr-C", "BCT", "I4/mmm", 2.87, 2.90,
                   "steel – stainless martensitic", centering="I")
    add_tetragonal("Tool/high-speed steel martensite", "Fe-C-Cr-W-Mo-V", "BCT", "I4/mmm", 2.87, 2.91,
                   "steel – tool/HSS", centering="I")

    # ------------------------------------------------------------------
    # ALUMINIUM AND COMMON ALUMINIUM-ALLOY PHASES
    # ------------------------------------------------------------------
    add_cubic("Aluminium / α-Al alloy matrix", "Al", "FCC", "Fm-3m", 4.0495, "aluminium and alloys",
              "Primary FCC aluminium matrix for 1xxx and the matrix of most wrought/cast Al alloys.")
    add_cubic("2xxx Al-Cu alloy matrix α-Al", "Al-Cu", "FCC", "Fm-3m", 4.045, "aluminium alloy – 2xxx")
    add_cubic("5xxx Al-Mg alloy matrix α-Al", "Al-Mg", "FCC", "Fm-3m", 4.055, "aluminium alloy – 5xxx")
    add_cubic("6xxx Al-Mg-Si alloy matrix α-Al", "Al-Mg-Si", "FCC", "Fm-3m", 4.050, "aluminium alloy – 6xxx")
    add_cubic("7xxx Al-Zn-Mg-Cu alloy matrix α-Al", "Al-Zn-Mg-Cu", "FCC", "Fm-3m", 4.045, "aluminium alloy – 7xxx")
    add_tetragonal("θ-Al2Cu", "Al2Cu", "tetragonal", "I4/mcm", 6.067, 4.877, "aluminium alloy – precipitate",
                   "Common equilibrium Al-Cu phase; exact intensities require a structural model.", centering="I")
    add_cubic("β-Mg2Si", "Mg2Si", "cubic antifluorite", "Fm-3m", 6.35, "aluminium alloy – precipitate", lattice="FCC")
    add_cubic("Silicon – Al-Si alloys", "Si", "diamond cubic", "Fd-3m", 5.431, "aluminium alloy – Al-Si", lattice="DIAMOND")
    add_hex("η-MgZn2 (C14 reference)", "MgZn2", "hexagonal Laves C14", "P63/mmc", 5.15, 8.48,
            "aluminium alloy – 7xxx precipitate", "Candidate positional reference; use COD/Rietveld for rigorous identification.")
    add_cubic("Al3Sc / coherent L12 reference", "Al3Sc", "L12 ordered cubic", "Pm-3m", 4.103,
              "aluminium alloy – precipitate", lattice="P")
    add_cubic("Al3Zr / coherent L12 reference", "Al3Zr", "L12 ordered cubic", "Pm-3m", 4.08,
              "aluminium alloy – precipitate", lattice="P")

    # ------------------------------------------------------------------
    # COPPER AND COPPER ALLOYS
    # ------------------------------------------------------------------
    add_cubic("Copper", "Cu", "FCC", "Fm-3m", 3.6149, "copper and alloys")
    add_cubic("α-brass Cu-Zn (FCC reference)", "Cu-Zn", "FCC", "Fm-3m", 3.68, "copper alloy – brass",
              "Representative alpha-brass matrix; lattice parameter increases with Zn content.")
    add_cubic("β-brass CuZn (B2)", "Cu-Zn", "B2 ordered cubic", "Pm-3m", 2.945, "copper alloy – brass", lattice="P")
    add_cubic("α-bronze Cu-Sn (FCC reference)", "Cu-Sn", "FCC", "Fm-3m", 3.66, "copper alloy – bronze")
    add_cubic("Cu-Al bronze α (FCC reference)", "Cu-Al", "FCC", "Fm-3m", 3.64, "copper alloy – aluminium bronze")
    add_cubic("Cu-Al bronze β/B2 reference", "Cu-Al", "B2/BCC", "Pm-3m", 2.93, "copper alloy – aluminium bronze", lattice="P")
    add_cubic("Cu-Ni solid solution", "Cu-Ni", "FCC", "Fm-3m", 3.59, "copper alloy – Cu-Ni",
              "Substitutional FCC solid solution; the actual lattice parameter depends strongly on Ni/Fe/Mn content.")
    add_cubic("CuNi10Fe1Mn reference", "Cu-Ni-Fe-Mn", "FCC", "Fm-3m", 3.60, "copper alloy – Cu-Ni",
              "Composition-specific peak positions are not fixed; the program fits a to the measured pattern.")
    add_cubic("CuNi30Mn1Fe reference", "Cu-Ni-Fe-Mn", "FCC", "Fm-3m", 3.59, "copper alloy – Cu-Ni",
              "Reference value is only a starting guess; use fitted a from the measured pattern.")
    add_cubic("67ERCuNi / ERCuNi reference", "Cu-Ni-Fe-Mn", "FCC", "Fm-3m", 3.58, "copper alloy – filler wire",
              "Reference value is only a starting guess; ERCuNi chemistry and residual stress shift the FCC reflections.")

    # ------------------------------------------------------------------
    # TITANIUM AND TITANIUM ALLOYS
    # ------------------------------------------------------------------
    add_hex("Titanium α / commercially pure Ti", "Ti", "HCP", "P63/mmc", 2.951, 4.683, "titanium and alloys")
    add_cubic("Titanium β", "Ti", "BCC", "Im-3m", 3.306, "titanium and alloys")
    add_hex("Ti-6Al-4V α (HCP reference)", "Ti-Al-V", "HCP", "P63/mmc", 2.924, 4.665, "titanium alloy – Ti-6Al-4V",
            "Seed updated from the supplied 0/750/950/1100 °C air-cooled Ti-6Al-4V XRD series. These values are only starting parameters; a and c are fitted independently for every measured sample.")
    add_cubic("Ti-6Al-4V β (BCC reference)", "Ti-Al-V", "BCC", "Im-3m", 3.24, "titanium alloy – Ti-6Al-4V")
    add_hex("α2-Ti3Al (D019 positional reference)", "Ti3Al", "hexagonal ordered D019", "P63/mmc", 5.78, 4.64,
            "titanium alloy – intermetallic", "Use COD for exact ordered-phase intensities.")
    add_tetragonal("γ-TiAl (L10 reference)", "TiAl", "tetragonal L10", "P4/mmm", 3.998, 4.07,
                   "titanium alloy – intermetallic", centering="P")

    # ------------------------------------------------------------------
    # ZIRCONIUM AND ZIRCONIUM ALLOYS
    # ------------------------------------------------------------------
    add_hex("Zirconium α", "Zr", "HCP", "P63/mmc", 3.231, 5.147, "zirconium and alloys")
    add_cubic("Zirconium β", "Zr", "BCC", "Im-3m", 3.61, "zirconium and alloys")
    add_hex("Zircaloy / Zr-alloy α matrix", "Zr-Sn-Fe-Cr", "HCP", "P63/mmc", 3.23, 5.15,
            "zirconium alloy", "Representative alpha-Zr matrix for zirconium alloys; exact lattice parameters depend on chemistry.")

    # ------------------------------------------------------------------
    # MAGNESIUM AND MAGNESIUM ALLOYS
    # ------------------------------------------------------------------
    add_hex("Magnesium / α-Mg matrix", "Mg", "HCP", "P63/mmc", 3.2094, 5.210, "magnesium and alloys")
    add_hex("AZ31/AZ91 α-Mg matrix reference", "Mg-Al-Zn-Mn", "HCP", "P63/mmc", 3.209, 5.211,
            "magnesium alloy – AZ series")
    add_cubic("β-Mg17Al12 positional reference", "Mg17Al12", "complex cubic", "I-43m", 10.56,
              "magnesium alloy – precipitate", "Positional screening reference; complex basis means calculated intensities are not represented.", lattice="BCC")
    # Mg2Si and MgZn2 are already present above and are also relevant to Mg alloys.

    # Additional common elemental metallic references retained for general use.
    add_cubic("Nickel", "Ni", "FCC", "Fm-3m", 3.5238, "metal")
    add_cubic("Chromium", "Cr", "BCC", "Im-3m", 2.8840, "metal")
    add_cubic("Molybdenum", "Mo", "BCC", "Im-3m", 3.1470, "metal")
    add_cubic("Tungsten", "W", "BCC", "Im-3m", 3.1652, "metal")
    add_cubic("Vanadium", "V", "BCC", "Im-3m", 3.0300, "metal")
    add_cubic("Niobium", "Nb", "BCC", "Im-3m", 3.3004, "metal")
    add_cubic("Tantalum", "Ta", "BCC", "Im-3m", 3.3058, "metal")
    add_hex("Zinc", "Zn", "HCP", "P63/mmc", 2.6649, 4.9468, "metal")
    add_hex("Cobalt (HCP)", "Co", "HCP", "P63/mmc", 2.5071, 4.0695, "metal")
    add_cubic("Cobalt (FCC)", "Co", "FCC", "Fm-3m", 3.5447, "metal")
    add_cubic("Silver", "Ag", "FCC", "Fm-3m", 4.0862, "metal")
    add_cubic("Gold", "Au", "FCC", "Fm-3m", 4.0782, "metal")
    add_cubic("Lead", "Pb", "FCC", "Fm-3m", 4.9502, "metal")
    add_cubic("Platinum", "Pt", "FCC", "Fm-3m", 3.9231, "metal")
    add_cubic("Palladium", "Pd", "FCC", "Fm-3m", 3.8898, "metal")

    # Selected ordered/intermetallic phases useful in alloy work.
    add_cubic("Ni3Al γ'", "Ni3Al", "L12 ordered cubic", "Pm-3m", 3.572, "intermetallic", lattice="P")
    add_cubic("Fe3Al DO3", "Fe3Al", "ordered cubic", "Fm-3m", 5.79, "intermetallic",
              "Superlattice reflections depend on chemical order.", lattice="P")

    return phases


def seed_experimental_reference_patterns():
    """Measured reference profiles bundled with this build.

    These are empirical benchmark patterns from the user-supplied XRDML series, not
    universal phase standards and not temperature classifiers.  They are stored in a
    separate internal table and are intentionally excluded from automatic phase
    assignment.  Automatic identification continues to use crystallographic phase
    models (with fitted lattice parameters), while these records provide an auditable
    comparison/reference history.
    """
    return [
        {
            "name": "X5CrNi18-10 reference – 0 °C / unheated",
            "material": "X5CrNi18-10 / EN 1.4301", "treatment": "0 °C reference",
            "phase_label": "γ-austenite FCC", "structure": "FCC", "wavelength_A": 1.540598,
            "a_fit_A": 3.5865, "c_fit_A": None,
            "notes": "Measured XRDML benchmark. Dominant FCC austenite; no independent strong BCC ferrite series detected.",
            "peaks": [(43.7196,"111",100.0),(50.9884,"200",21.8),(74.7469,"220",9.6),(90.8238,"311",13.8),(95.9846,"222",7.2),(118.3155,"400",1.6)],
        },
        {
            "name": "X5CrNi18-10 reference – 750 °C air-cooled",
            "material": "X5CrNi18-10 / EN 1.4301", "treatment": "750 °C, free air cooling",
            "phase_label": "γ-austenite FCC", "structure": "FCC", "wavelength_A": 1.540598,
            "a_fit_A": 3.5854, "c_fit_A": None,
            "notes": "Measured XRDML benchmark. FCC sequence (111),(200),(220),(311),(222); temperature history is metadata, not inferred from XRD peak position alone.",
            "peaks": [(43.7189,"111",100.0),(50.9657,"200",21.6),(74.8768,"220",9.8),(90.8191,"311",13.5),(95.9830,"222",6.4)],
        },
        {
            "name": "X5CrNi18-10 reference – 950 °C air-cooled",
            "material": "X5CrNi18-10 / EN 1.4301", "treatment": "950 °C, free air cooling",
            "phase_label": "γ-austenite FCC", "structure": "FCC", "wavelength_A": 1.540598,
            "a_fit_A": 3.5862, "c_fit_A": None,
            "notes": "Measured XRDML benchmark. Dominant FCC austenite.",
            "peaks": [(43.7191,"111",100.0),(50.9700,"200",21.9),(74.7342,"220",9.4),(90.9061,"311",13.6),(96.0237,"222",5.5),(118.3235,"400",1.8)],
        },
        {
            "name": "X5CrNi18-10 reference – 1100 °C air-cooled",
            "material": "X5CrNi18-10 / EN 1.4301", "treatment": "1100 °C, free air cooling",
            "phase_label": "γ-austenite FCC", "structure": "FCC", "wavelength_A": 1.540598,
            "a_fit_A": 3.5876, "c_fit_A": None,
            "notes": "Measured XRDML benchmark. Dominant FCC austenite after free cooling.",
            "peaks": [(43.7019,"111",100.0),(50.9371,"200",23.3),(74.7440,"220",9.8),(90.7923,"311",16.2),(96.0413,"222",8.3),(118.1785,"400",2.2)],
        },
        {
            "name": "Ti-6Al-4V reference – 0 °C / unheated",
            "material": "Ti-6Al-4V", "treatment": "0 °C reference",
            "phase_label": "α-Ti HCP (dominant)", "structure": "HCP", "wavelength_A": 1.540598,
            "a_fit_A": 2.92235, "c_fit_A": 4.66088,
            "notes": "Measured XRDML benchmark. Dominant α-Ti HCP. A weak/overlapped β-Ti contribution cannot be quantified from these positional data alone.",
            "peaks": [(35.4876,"100",18.2),(38.6408,"002",29.7),(40.5549,"101",100.0),(53.4493,"102",16.6),(63.6733,"110",25.0),(71.1958,"103",19.5),(76.9927,"112",20.2),(78.2657,"201",4.2),(82.7305,"004",2.6),(87.6550,"202",1.8),(93.4852,"104",2.3),(103.4423,"203",3.9),(110.5465,"211",8.0),(115.3679,"114",8.1)],
        },
        {
            "name": "Ti-6Al-4V reference – 750 °C air-cooled",
            "material": "Ti-6Al-4V", "treatment": "750 °C, free air cooling",
            "phase_label": "α-Ti HCP (dominant)", "structure": "HCP", "wavelength_A": 1.540598,
            "a_fit_A": 2.92320, "c_fit_A": 4.66381,
            "notes": "Measured XRDML benchmark. Dominant α-Ti HCP; peak-intensity changes can reflect texture/microstructure and are not a unique thermometer.",
            "peaks": [(35.4727,"100",18.5),(38.6340,"002",31.1),(40.5327,"101",100.0),(53.4360,"102",15.4),(63.6282,"110",21.4),(71.1660,"103",24.5),(76.9368,"112",16.9),(78.2631,"201",4.2),(82.7809,"004",3.0),(87.5068,"202",1.9),(93.2610,"104",3.0),(103.2846,"203",3.3),(110.5894,"211",7.2),(115.3685,"114",7.9)],
        },
        {
            "name": "Ti-6Al-4V reference – 950 °C air-cooled",
            "material": "Ti-6Al-4V", "treatment": "950 °C, free air cooling",
            "phase_label": "α-Ti HCP (dominant)", "structure": "HCP", "wavelength_A": 1.540598,
            "a_fit_A": 2.92503, "c_fit_A": 4.66788,
            "notes": "Measured XRDML benchmark. Dominant α-Ti HCP after free cooling. Do not report β-Ti from a single overlapped reflection.",
            "peaks": [(35.4732,"100",16.0),(38.6228,"002",29.5),(40.5398,"101",100.0),(53.3781,"102",16.0),(63.5755,"110",26.1),(71.0876,"103",20.8),(76.8768,"112",19.8),(78.0755,"201",3.9),(82.5862,"004",3.0),(87.5647,"202",1.9),(93.2281,"104",2.7),(103.3043,"203",3.9),(110.4218,"211",9.2),(115.1777,"114",8.8)],
        },
        {
            "name": "Ti-6Al-4V reference – 1100 °C air-cooled",
            "material": "Ti-6Al-4V", "treatment": "1100 °C, free air cooling",
            "phase_label": "α-Ti HCP (dominant after cooling)", "structure": "HCP", "wavelength_A": 1.540598,
            "a_fit_A": 2.92623, "c_fit_A": 4.66694,
            "notes": "Measured XRDML benchmark after heating above the usual Ti-6Al-4V beta-transus region and free cooling. Final room-temperature pattern is still dominated by HCP α/α′-type reflections; the prior high-temperature β state must not be inferred as retained β solely from heat-treatment temperature.",
            "peaks": [(35.4417,"100",16.9),(38.5801,"002",21.8),(40.5073,"101",100.0),(53.3713,"102",7.7),(63.5514,"110",4.7),(71.0378,"103",9.7),(76.8798,"112",8.6),(78.1545,"201",5.4),(82.6384,"004",2.5),(87.4960,"202",2.5),(93.1770,"104",0.9),(103.2094,"203",3.6),(110.3172,"211",3.1),(115.2838,"114",1.8)],
        },
    ]


class PhaseDatabase:
    def __init__(self):
        data_dir = Path.home() / "XRDAnalyzerData"
        data_dir.mkdir(parents=True, exist_ok=True)
        self.path = data_dir / "xrd_phase_database.sqlite"
        self.conn = sqlite3.connect(self.path)
        self.conn.execute("PRAGMA foreign_keys=ON")
        self._create_schema()
        self._seed_if_empty()

    def _create_schema(self):
        self.conn.execute("""
            CREATE TABLE IF NOT EXISTS phases (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                name TEXT NOT NULL UNIQUE,
                formula TEXT,
                structure TEXT,
                space_group TEXT,
                category TEXT,
                notes TEXT,
                lattice_model TEXT DEFAULT 'FIXED',
                a_ref REAL,
                c_ref REAL,
                fit_lattice INTEGER DEFAULT 0
            )
        """)
        # Migrate databases created by older XRD Analyzer versions.
        cols = {r[1] for r in self.conn.execute("PRAGMA table_info(phases)").fetchall()}
        for col, ddl in [
            ("lattice_model", "TEXT DEFAULT 'FIXED'"),
            ("a_ref", "REAL"),
            ("c_ref", "REAL"),
            ("fit_lattice", "INTEGER DEFAULT 0"),
        ]:
            if col not in cols:
                self.conn.execute(f"ALTER TABLE phases ADD COLUMN {col} {ddl}")
        self.conn.execute("""
            CREATE TABLE IF NOT EXISTS phase_peaks (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                phase_id INTEGER NOT NULL,
                d_A REAL NOT NULL,
                hkl TEXT,
                rel_intensity REAL DEFAULT 100,
                FOREIGN KEY(phase_id) REFERENCES phases(id) ON DELETE CASCADE
            )
        """)
        self.conn.execute("""
            CREATE TABLE IF NOT EXISTS experimental_references (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                name TEXT NOT NULL UNIQUE,
                material TEXT,
                treatment TEXT,
                phase_label TEXT,
                structure TEXT,
                wavelength_A REAL,
                a_fit_A REAL,
                c_fit_A REAL,
                notes TEXT
            )
        """)
        self.conn.execute("""
            CREATE TABLE IF NOT EXISTS experimental_reference_peaks (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                reference_id INTEGER NOT NULL,
                two_theta REAL NOT NULL,
                d_A REAL,
                hkl TEXT,
                rel_intensity REAL DEFAULT 100,
                FOREIGN KEY(reference_id) REFERENCES experimental_references(id) ON DELETE CASCADE
            )
        """)
        self.conn.commit()

    def _seed_if_empty(self):
        # Insert new built-in references and migrate older installations to the
        # expanded reflection lists. User-added/COD entries have different names
        # and are therefore left untouched.
        for name, formula, structure, sg, category, notes, peaks, lattice_model, a_ref, c_ref, fit_lattice in seed_phase_definitions():
            row = self.conn.execute("SELECT id FROM phases WHERE name=?", (name,)).fetchone()
            if row:
                pid = int(row[0])
                old_count = int(self.conn.execute(
                    "SELECT COUNT(*) FROM phase_peaks WHERE phase_id=?", (pid,)
                ).fetchone()[0])
                # Update metadata and replace the old short peak list only when
                # the new built-in definition is more complete.
                self.conn.execute(
                    "UPDATE phases SET formula=?,structure=?,space_group=?,category=?,notes=?,lattice_model=?,a_ref=?,c_ref=?,fit_lattice=? WHERE id=?",
                    (formula, structure, sg, category, notes, lattice_model, a_ref, c_ref, fit_lattice, pid)
                )
                # Built-in references are regenerated on every upgrade. Their stored
                # d-spacings are SEARCH SEEDS only for fit-enabled alloy matrices, not
                # fixed experimental peak positions. User/COD entries have other names.
                self.conn.execute("DELETE FROM phase_peaks WHERE phase_id=?", (pid,))
                for dA, hkl in peaks:
                    self.conn.execute(
                        "INSERT INTO phase_peaks(phase_id,d_A,hkl,rel_intensity) VALUES (?,?,?,?)",
                        (pid, float(dA), hkl, 100.0)
                    )
                continue
            cur = self.conn.execute(
                "INSERT INTO phases(name, formula, structure, space_group, category, notes, lattice_model, a_ref, c_ref, fit_lattice) VALUES (?,?,?,?,?,?,?,?,?,?)",
                (name, formula, structure, sg, category, notes, lattice_model, a_ref, c_ref, fit_lattice)
            )
            pid = cur.lastrowid
            for dA, hkl in peaks:
                self.conn.execute(
                    "INSERT INTO phase_peaks(phase_id,d_A,hkl,rel_intensity) VALUES (?,?,?,?)",
                    (pid, float(dA), hkl, 100.0)
                )

        # Experimental benchmark patterns are stored separately from phase models.
        # They document measured results but are NOT used as automatic phase labels,
        # because a heat-treatment temperature cannot be identified uniquely from XRD.
        for rec in seed_experimental_reference_patterns():
            row = self.conn.execute("SELECT id FROM experimental_references WHERE name=?", (rec["name"],)).fetchone()
            if row:
                rid = int(row[0])
                self.conn.execute(
                    "UPDATE experimental_references SET material=?,treatment=?,phase_label=?,structure=?,wavelength_A=?,a_fit_A=?,c_fit_A=?,notes=? WHERE id=?",
                    (rec.get("material",""), rec.get("treatment",""), rec.get("phase_label",""),
                     rec.get("structure",""), rec.get("wavelength_A"), rec.get("a_fit_A"),
                     rec.get("c_fit_A"), rec.get("notes",""), rid)
                )
                self.conn.execute("DELETE FROM experimental_reference_peaks WHERE reference_id=?", (rid,))
            else:
                cur = self.conn.execute(
                    "INSERT INTO experimental_references(name,material,treatment,phase_label,structure,wavelength_A,a_fit_A,c_fit_A,notes) VALUES (?,?,?,?,?,?,?,?,?)",
                    (rec["name"], rec.get("material",""), rec.get("treatment",""), rec.get("phase_label",""),
                     rec.get("structure",""), rec.get("wavelength_A"), rec.get("a_fit_A"), rec.get("c_fit_A"), rec.get("notes",""))
                )
                rid = int(cur.lastrowid)
            wave = float(rec.get("wavelength_A") or DEFAULT_WAVELENGTH_A)
            for tt, hkl, rel in rec.get("peaks", []):
                try:
                    dA = bragg_d(float(tt), wave)
                except Exception:
                    dA = None
                self.conn.execute(
                    "INSERT INTO experimental_reference_peaks(reference_id,two_theta,d_A,hkl,rel_intensity) VALUES (?,?,?,?,?)",
                    (rid, float(tt), dA, str(hkl), float(rel))
                )
        self.conn.commit()

    def all_phases(self):
        return self.conn.execute(
            "SELECT id,name,formula,structure,space_group,category,notes FROM phases ORDER BY category,name"
        ).fetchall()

    def phase_peaks(self, phase_id: int):
        return self.conn.execute(
            "SELECT d_A,hkl,rel_intensity FROM phase_peaks WHERE phase_id=? ORDER BY d_A DESC", (phase_id,)
        ).fetchall()

    def phase_lattice_model(self, phase_id: int):
        row = self.conn.execute(
            "SELECT lattice_model,a_ref,c_ref,fit_lattice FROM phases WHERE id=?", (phase_id,)
        ).fetchone()
        if not row:
            return "FIXED", None, None, False
        return (str(row[0] or "FIXED"), row[1], row[2], bool(row[3]))

    def all_experimental_references(self):
        return self.conn.execute(
            "SELECT id,name,material,treatment,phase_label,structure,wavelength_A,a_fit_A,c_fit_A,notes FROM experimental_references ORDER BY material,treatment,name"
        ).fetchall()

    def experimental_reference_peaks(self, reference_id: int):
        return self.conn.execute(
            "SELECT two_theta,d_A,hkl,rel_intensity FROM experimental_reference_peaks WHERE reference_id=? ORDER BY two_theta",
            (int(reference_id),)
        ).fetchall()

    def add_phase(self, name, formula, structure, sg, category, notes, peaks: List[Tuple[float,str,float]]):
        cur = self.conn.execute(
            "INSERT INTO phases(name,formula,structure,space_group,category,notes) VALUES (?,?,?,?,?,?)",
            (name,formula,structure,sg,category,notes)
        )
        pid = cur.lastrowid
        for dA,hkl,rel in peaks:
            self.conn.execute(
                "INSERT INTO phase_peaks(phase_id,d_A,hkl,rel_intensity) VALUES (?,?,?,?)",
                (pid,dA,hkl,rel)
            )
        self.conn.commit()

    def delete_phase(self, phase_id: int):
        self.conn.execute("DELETE FROM phases WHERE id=?", (phase_id,))
        self.conn.commit()


class CODDatabase:
    """Read-only connector for the Crystallography Open Database SQLite .db3 export."""
    REQUIRED_COLUMNS = {
        "file", "a", "b", "c", "alpha", "beta", "gamma", "sg", "sgNumber",
        "formula", "chemname", "mineral", "commonname", "doi", "journal", "year", "method"
    }

    def __init__(self, path: Optional[Path] = None):
        self.path: Optional[Path] = None
        self.conn: Optional[sqlite3.Connection] = None
        self.record_count: int = 0
        if path:
            self.open(path)

    @property
    def is_open(self) -> bool:
        return self.conn is not None

    def close(self):
        if self.conn is not None:
            try:
                self.conn.close()
            except Exception:
                pass
        self.conn = None
        self.path = None
        self.record_count = 0

    def open(self, path: Path):
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(path)
        if path.suffix.lower() not in (".db3", ".sqlite", ".sqlite3", ".db"):
            raise ValueError("Wskaż plik bazy COD SQLite (.db3/.sqlite/.db).")
        self.close()
        uri = f"file:{path.resolve().as_posix()}?mode=ro"
        conn = sqlite3.connect(uri, uri=True, timeout=20)
        conn.execute("PRAGMA query_only=ON")
        tables = {r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        if "data" not in tables:
            conn.close()
            raise ValueError("This SQLite file does not contain the COD 'data' table.")
        cols = {r[1] for r in conn.execute("PRAGMA table_info(data)")}
        missing = sorted(self.REQUIRED_COLUMNS - cols)
        if missing:
            conn.close()
            raise ValueError("COD database schema is missing columns: " + ", ".join(missing))
        self.conn = conn
        self.path = path.resolve()
        self.record_count = int(conn.execute("SELECT COUNT(*) FROM data").fetchone()[0])

    def search(self, text: str, limit: int = 200):
        if self.conn is None:
            raise RuntimeError("COD database is not loaded.")
        q = (text or "").strip()
        cols = "file,formula,chemname,mineral,commonname,sg,sgNumber,a,b,c,alpha,beta,gamma,doi,journal,year,method"
        if not q:
            sql = f"SELECT {cols} FROM data WHERE a IS NOT NULL AND b IS NOT NULL AND c IS NOT NULL LIMIT ?"
            return self.conn.execute(sql, (int(limit),)).fetchall()
        if q.isdigit():
            row = self.conn.execute(f"SELECT {cols} FROM data WHERE file=?", (int(q),)).fetchall()
            if row:
                return row
        like = f"%{q}%"
        sql = f"""
            SELECT {cols} FROM data
            WHERE formula LIKE ? OR chemname LIKE ? OR mineral LIKE ? OR commonname LIKE ?
               OR sg LIKE ? OR doi LIKE ?
            LIMIT ?
        """
        return self.conn.execute(sql, (like, like, like, like, like, like, int(limit))).fetchall()

    def get_entry(self, cod_id: int):
        if self.conn is None:
            raise RuntimeError("COD database is not loaded.")
        cols = "file,formula,chemname,mineral,commonname,sg,sgNumber,a,b,c,alpha,beta,gamma,doi,journal,year,method"
        return self.conn.execute(f"SELECT {cols} FROM data WHERE file=?", (int(cod_id),)).fetchone()

    def get_entry_details(self, cod_id: int) -> dict:
        """Return the useful descriptive/bibliographic metadata stored in the COD SQLite export."""
        if self.conn is None:
            raise RuntimeError("COD database is not loaded.")
        cols = [
            "file", "formula", "calcformula", "commonname", "chemname", "mineral",
            "sg", "sgHall", "sgNumber", "a", "b", "c", "alpha", "beta", "gamma", "vol",
            "celltemp", "diffrtemp", "cellpressure", "diffrpressure", "compoundsource", "nel",
            "authors", "title", "journal", "year", "volume", "issue", "firstpage", "lastpage", "doi",
            "method", "radiation", "wavelength", "radType", "radSymbol", "Rall", "Robs", "wRall", "wRobs",
            "gofall", "gofobs", "status", "flags", "date"
        ]
        row = self.conn.execute(f"SELECT {','.join(cols)} FROM data WHERE file=?", (int(cod_id),)).fetchone()
        if row is None:
            return {}
        out = {}
        for k, v in zip(cols, row):
            if isinstance(v, str):
                v = self._clean_cod_text(v)
            out[k] = v
        return out

    @staticmethod
    def _clean_cod_text(v) -> str:
        if v is None:
            return ""
        s = str(v).strip()
        return "" if s in ("\\N", "?", ".") else s

    def candidate_structures(self, allowed_elements: set[str], max_candidates: int = 5000,
                             max_cell_A: float = 25.0,
                             priority_weights: Optional[Dict[str, float]] = None):
        """
        Return COD structures whose chemical formula contains only the allowed elements.
        This first-stage chemistry/cell filter keeps the measured-pattern search fast enough
        for the ~500k-record COD SQLite export.
        """
        if self.conn is None:
            raise RuntimeError("COD database is not loaded.")
        allowed = {e.strip() for e in allowed_elements if e.strip()}
        if not allowed:
            raise ValueError(
                "Enter at least one allowed element (for example: Fe Cr Ni Mo Cu Mn). "
                "A completely unrestricted scan of the full COD is intentionally disabled."
            )
        patterns = [f"%{e}%" for e in sorted(allowed)]
        formula_where = " OR ".join(["formula LIKE ?"] * len(patterns))
        max_nel = max(1, len(allowed))
        sql = f"""
            SELECT file,formula,chemname,mineral,commonname,sg,sgNumber,
                   a,b,c,alpha,beta,gamma,year,method,nel
            FROM data
            WHERE a IS NOT NULL AND b IS NOT NULL AND c IS NOT NULL
              AND alpha IS NOT NULL AND beta IS NOT NULL AND gamma IS NOT NULL
              AND a>0 AND b>0 AND c>0
              AND a<? AND b<? AND c<?
              AND CAST(nel AS INTEGER)<=?
              AND ({formula_where})
        """
        rows = self.conn.execute(sql, (float(max_cell_A), float(max_cell_A), float(max_cell_A), max_nel, *patterns)).fetchall()
        kept = []
        for row in rows:
            els = formula_elements(self._clean_cod_text(row[1]))
            if not els or not els.issubset(allowed):
                continue
            # Prefer powder-diffraction records, chemically richer entries and newer records
            # when the chemistry filter still leaves more structures than requested.
            method = self._clean_cod_text(row[14]).lower()
            powder_bonus = 1 if "powder" in method else 0
            try:
                year = int(float(row[13])) if row[13] not in (None, "", "\\N", "?", ".") else 0
            except Exception:
                year = 0
            chem_capture = float(sum((priority_weights or {}).get(e, 0.0) for e in els))
            dominant_bonus = 0
            if priority_weights:
                dom = max(priority_weights, key=priority_weights.get)
                dominant_bonus = 1 if dom in els else 0
            kept.append((powder_bonus, dominant_bonus, chem_capture, len(els), year, row, els))
        # Preserve experimentally useful powder records, but keep chemistry relevance ahead of
        # mere formula richness/year so the 5k-candidate cap does not discard the likely matrix.
        kept.sort(key=lambda x: (x[0], x[1], x[2], x[3], x[4]), reverse=True)
        return [(x[5], x[6]) for x in kept[:int(max_candidates)]]

    @staticmethod
    def _row_to_reflections(row, wavelength_A: float, two_theta_min: float, two_theta_max: float):
        """Generate reflection positions from COD unit-cell/space-group data (no intensities)."""
        if not ensure_gemmi():
            raise RuntimeError(gemmi_diagnostic_message("Measured-pattern COD search"))
        (file_id, formula, chemname, mineral, commonname, sg_text, sg_no,
         a,b,c,alpha,beta,gamma,year,method,nel) = row
        sg_text = CODDatabase._clean_cod_text(sg_text)
        try:
            sg = gemmi.SpaceGroup(sg_text) if sg_text else gemmi.SpaceGroup(int(sg_no))
            cell = gemmi.UnitCell(float(a), float(b), float(c), float(alpha), float(beta), float(gamma))
        except Exception as e:
            raise ValueError(f"Invalid COD unit cell / space group: {e}")

        theta_min = math.radians(max(0.001, float(two_theta_min)) / 2.0)
        theta_max = math.radians(min(179.0, float(two_theta_max)) / 2.0)
        d_max = wavelength_A / (2.0 * math.sin(theta_min))
        d_min = wavelength_A / (2.0 * math.sin(theta_max))
        hkls = gemmi.make_miller_array(cell, sg, float(d_min), float(d_max), True)
        refs = []
        for hkl in hkls:
            try:
                d = float(cell.calculate_d(hkl))
                tt = d_to_2theta(d, wavelength_A)
            except Exception:
                continue
            if tt is None or tt < two_theta_min or tt > two_theta_max:
                continue
            h,k,l = (int(hkl[0]), int(hkl[1]), int(hkl[2]))
            refs.append((float(tt), float(d), f"{h}{k}{l}"))
        refs.sort(key=lambda x: x[0])
        return refs

    def search_measured_pattern(self, dataset: Dataset, allowed_elements: set[str],
                                tolerance_deg: float = 0.35, max_candidates: int = 5000,
                                top_results: int = 25, max_observed_peaks: int = 8,
                                min_relative_intensity: float = 5.0,
                                progress_callback=None, cancel_callback=None,
                                exclude_resolved_kalpha2: bool = True):
        """
        Szukaj w COD unit-cell/space-group records against measured peak positions.

        The score uses weighted observed-peak coverage, RMS 2theta error and a penalty
        for overly dense predicted reflection sets. Because COD .db3 does not provide a
        ready powder-intensity pattern for every record, this is a candidate-screening
        tool, not a Rietveld/refinement result.
        """
        if self.conn is None:
            raise RuntimeError("COD database is not loaded.")
        if not dataset.peaks:
            raise ValueError("Brak wykrytych pików. Najpierw przeanalizuj bieżącą próbkę.")

        excluded_indices = set(dataset_kalpha2_candidates(dataset)) if exclude_resolved_kalpha2 else set()
        obs = [p for idx,p in enumerate(dataset.peaks)
               if idx not in excluded_indices and p.rel_intensity >= float(min_relative_intensity)]
        obs.sort(key=lambda p: p.height, reverse=True)
        obs = obs[:max(2, int(max_observed_peaks))]
        if len(obs) < 2:
            raise ValueError("Wyszukiwanie wzorca COD wymaga co najmniej dwóch zmierzonych pików.")
        obs_tt = np.asarray([p.two_theta for p in obs], dtype=float)
        weights = np.sqrt(np.asarray([max(p.rel_intensity, 0.1) for p in obs], dtype=float) / 100.0)
        weights = weights / weights.sum()

        priority = dataset_element_priority_weights(dataset)
        candidates = self.candidate_structures(allowed_elements, max_candidates=max_candidates, priority_weights=priority)
        results = []
        n = len(candidates)
        tol = max(1e-6, float(tolerance_deg))
        for i, (row, row_elements) in enumerate(candidates):
            if cancel_callback and cancel_callback():
                break
            if progress_callback and (i % 10 == 0 or i == n-1):
                progress_callback(i+1, n)
            try:
                refs = self._row_to_reflections(row, dataset.wavelength_A, float(dataset.x.min()), float(dataset.x.max()))
                if not refs:
                    continue
                ref_tt = np.asarray([r[0] for r in refs], dtype=float)
                obs_diff = np.asarray([float(np.min(np.abs(ref_tt - o))) for o in obs_tt])
                matched = obs_diff <= tol
                match_count = int(np.sum(matched))
                if match_count < 2:
                    continue

                weighted_coverage = float(np.sum(weights[matched]))
                rms = float(np.sqrt(np.average(obs_diff[matched]**2, weights=weights[matched])))
                # Dense low-symmetry patterns can match almost anything by chance. Penalize them
                # gently because generated reflection positions do not include structure-factor intensities.
                expected_sparse = max(8.0, 2.2 * len(obs_tt))
                density_penalty = 1.0 / (1.0 + max(0.0, len(refs)-expected_sparse) / (4.0*expected_sparse))
                rms_factor = math.exp(-0.5 * (rms / tol)**2)
                chem_factor, chem_capture, chem_reason = phase_chemistry_factor(dataset, row_elements)
                sg_text_for_factor = self._clean_cod_text(row[5])
                sg_no_for_factor = row[6]
                structure_factor, structure_reason = phase_structure_factor(
                    dataset, space_group=sg_text_for_factor, sg_no=sg_no_for_factor
                )
                method = self._clean_cod_text(row[14]).lower()
                method_factor = 1.03 if "powder" in method else 1.0
                score = min(100.0, 100.0 * weighted_coverage * rms_factor * density_penalty
                            * chem_factor * structure_factor * method_factor)

                matches = []
                for p, diff, ok in zip(obs, obs_diff, matched):
                    nearest_idx = int(np.argmin(np.abs(ref_tt - p.two_theta)))
                    rtt, rd, rhkl = refs[nearest_idx]
                    matches.append({
                        "obs_2theta": p.two_theta, "obs_d": p.d_A, "obs_rel": p.rel_intensity,
                        "ref_2theta": rtt, "ref_d": rd, "hkl": canonical_hkl_for_structure(rhkl, cod_structure_short(sg_text_for_factor, sg_no_for_factor)),
                        "delta": float(diff), "matched": bool(ok),
                    })
                (file_id, formula, chemname, mineral, commonname, sg_text, sg_no,
                 a,b,c,alpha,beta,gamma,year,method_raw,nel) = row
                display_name = self._clean_cod_text(chemname) or self._clean_cod_text(mineral) or self._clean_cod_text(commonname) or self._clean_cod_text(formula) or "COD structure"
                friendly_name, plot_name = cod_friendly_phase_names(
                    self._clean_cod_text(formula), self._clean_cod_text(chemname),
                    self._clean_cod_text(mineral), self._clean_cod_text(commonname),
                    self._clean_cod_text(sg_text), sg_no
                )
                contextual_name = contextual_cod_phase_name(
                    dataset, self._clean_cod_text(formula), plot_name,
                    self._clean_cod_text(sg_text), sg_no
                )
                if contextual_name != plot_name:
                    friendly_name = contextual_name
                    plot_name = contextual_name
                results.append({
                    "score": score, "matched": match_count, "used": len(obs), "rms": rms,
                    "chem_factor": float(chem_factor), "chem_capture": float(chem_capture),
                    "chem_reason": chem_reason, "structure_factor": float(structure_factor),
                    "structure_reason": structure_reason,
                    "ref_count": len(refs), "cod_id": int(file_id),
                    "formula": self._clean_cod_text(formula), "name": friendly_name, "plot_name": plot_name,
                    "cod_name": display_name, "chemname": self._clean_cod_text(chemname),
                    "mineral": self._clean_cod_text(mineral), "commonname": self._clean_cod_text(commonname),
                    "space_group": self._clean_cod_text(sg_text), "sg_no": sg_no,
                    "a": a, "b": b, "c": c, "alpha": alpha, "beta": beta, "gamma": gamma,
                    "year": year, "method": self._clean_cod_text(method_raw),
                    "matches": matches,
                })
            except Exception:
                continue

        results.sort(key=lambda r: (r["score"], r["matched"], -r["rms"]), reverse=True)

        # Candidate clustering: established search-match programs suppress repeated entries
        # representing practically the same calculated pattern. Keep the highest-ranked member
        # of each geometry cluster so chemistry-relevant alternatives remain visible instead of
        # dozens of near-duplicate COD cards.
        unique = []
        fingerprints = set()
        for r in results:
            struct = cod_structure_short(r.get("space_group", ""), r.get("sg_no"))
            geom = tuple(round(float(m.get("ref_2theta",0.0))/0.18) for m in r.get("matches",[]) if m.get("matched"))
            fp = (struct, geom)
            if fp in fingerprints:
                continue
            fingerprints.add(fp)
            unique.append(r)
            if len(unique) >= int(top_results):
                break
        return unique, obs, len(candidates)

    def entry_to_reference(self, cod_id: int, d_min_A: float = 0.70, d_max_A: float = 15.0):
        """
        Convert a selected COD unit cell + space group into a positional powder reference.
        Intensities are not calculated here; matching uses peak positions (d / 2theta).
        """
        if not ensure_gemmi():
            raise RuntimeError(gemmi_diagnostic_message("COD reflection generation"))
        row = self.get_entry(cod_id)
        if row is None:
            raise ValueError(f"COD entry {cod_id} was not found.")
        (file_id, formula, chemname, mineral, commonname, sg_text, sg_no,
         a,b,c,alpha,beta,gamma,doi,journal,year,method) = row
        vals = (a,b,c,alpha,beta,gamma)
        if any(v is None for v in vals):
            raise ValueError(f"COD {file_id} has incomplete unit-cell parameters.")
        a,b,c,alpha,beta,gamma = map(float, vals)
        sg_text = self._clean_cod_text(sg_text)
        try:
            sg = gemmi.SpaceGroup(sg_text) if sg_text else gemmi.SpaceGroup(int(sg_no))
        except Exception as e:
            raise ValueError(f"Cannot interpret COD space group '{sg_text}' ({sg_no}): {e}")
        cell = gemmi.UnitCell(a,b,c,alpha,beta,gamma)
        ops = sg.operations()

        # A conservative index limit works well for metals while still handling moderate cells.
        # The cap avoids pathological runtimes for very large molecular cells.
        max_axis = max(a,b,c)
        hmax = max(8, min(28, int(math.ceil(max_axis / max(d_min_A, 0.35))) + 2))
        refl_by_d = {}
        for h in range(-hmax, hmax+1):
            for k in range(-hmax, hmax+1):
                for l in range(-hmax, hmax+1):
                    if h == 0 and k == 0 and l == 0:
                        continue
                    hkl = [h,k,l]
                    try:
                        if ops.is_systematically_absent(hkl):
                            continue
                        d = float(cell.calculate_d(hkl))
                    except Exception:
                        continue
                    if not math.isfinite(d) or d < d_min_A or d > d_max_A:
                        continue
                    key = round(d, 5)
                    # Prefer a compact all-positive/low-index label when equivalent d values occur.
                    score = (sum(abs(v) for v in hkl), sum(1 for v in hkl if v < 0))
                    old = refl_by_d.get(key)
                    if old is None or score < old[0]:
                        def idx_text(v):
                            return f"-{abs(v)}" if v < 0 else str(v)
                        refl_by_d[key] = (score, f"{idx_text(h)}{idx_text(k)}{idx_text(l)}")

        if not refl_by_d:
            raise ValueError(f"No usable reflections could be generated for COD {file_id}.")
        peaks = [(float(d), hkl, 100.0) for d,(_,hkl) in sorted(refl_by_d.items(), reverse=True)]

        formula_c = self._clean_cod_text(formula)
        chem_c = self._clean_cod_text(chemname)
        mineral_c = self._clean_cod_text(mineral)
        common_c = self._clean_cod_text(commonname)
        display = chem_c or mineral_c or common_c or formula_c or "COD structure"
        name = f"COD {file_id} – {display}"
        structure = f"COD crystal structure; SG {sg.hm}"
        notes_parts = [
            f"Imported from COD record {file_id}.",
            f"Cell: a={a:g}, b={b:g}, c={c:g} Å; α={alpha:g}, β={beta:g}, γ={gamma:g}°.",
            "Reference reflections were generated from the unit cell and systematic absences; intensities are not COD powder intensities."
        ]
        if doi and self._clean_cod_text(doi): notes_parts.append(f"DOI: {self._clean_cod_text(doi)}")
        if journal and self._clean_cod_text(journal): notes_parts.append(f"Journal: {self._clean_cod_text(journal)} {year or ''}".strip())
        return {
            "name": name,
            "formula": formula_c,
            "structure": structure,
            "space_group": sg.hm,
            "category": "COD",
            "notes": " ".join(notes_parts),
            "peaks": peaks,
        }


# --------------------------- parsers ---------------------------

def parse_cif(path: Path) -> Dataset:
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    meta = {}
    for line in lines:
        st = line.strip()
        if st.startswith("_pd_meas_2theta_range_min"):
            try: meta["min"] = float(st.split()[-1])
            except Exception: pass
        elif st.startswith("_pd_meas_2theta_range_inc"):
            try: meta["inc"] = float(st.split()[-1])
            except Exception: pass
        elif st.startswith("_diffrn_radiation_wavelength"):
            try: meta["wavelength"] = float(st.split()[-1])
            except Exception: pass
        elif st.startswith("_diffrn_source_target"):
            meta["anode"] = st.split(maxsplit=1)[1].strip(" '\"") if len(st.split(maxsplit=1)) > 1 else ""
        elif st.startswith("_diffrn_radiation_type"):
            meta["radiation"] = st.split(maxsplit=1)[1].strip(" '\"") if len(st.split(maxsplit=1)) > 1 else ""
        elif st.startswith("_pd_meas_datetime_initiated"):
            meta["measurement_datetime"] = st.split(maxsplit=1)[1].strip(" '\"") if len(st.split(maxsplit=1)) > 1 else ""

    y = None
    x = None
    i = 0
    while i < len(lines):
        if lines[i].strip() == "loop_":
            tags = []
            j = i + 1
            while j < len(lines) and lines[j].lstrip().startswith("_"):
                tags.append(lines[j].strip())
                j += 1
            if "_pd_meas_intensity_total" in tags:
                cols = len(tags)
                idx = tags.index("_pd_meas_intensity_total")
                vals = []
                k = j
                while k < len(lines):
                    st = lines[k].strip()
                    if not st or st == "loop_" or st.startswith("_"):
                        break
                    parts = st.split()
                    if len(parts) < cols:
                        break
                    try:
                        vals.append(float(parts[idx]))
                    except Exception:
                        break
                    k += 1
                y = np.asarray(vals, dtype=float)
                if "min" in meta and "inc" in meta:
                    x = meta["min"] + np.arange(len(y)) * meta["inc"]
                break
        i += 1

    if x is None or y is None or len(y) == 0:
        raise ValueError("W pliku CIF nie znaleziono profilu intensywności dyfrakcji proszkowej.")
    wave = float(meta.get("wavelength", DEFAULT_WAVELENGTH_A))
    ds = Dataset(path, path.stem, "", x, y, wave)
    ds.metadata["source_format"] = "CIF/pdCIF"
    ds.metadata_sources["source_format"] = "file"
    if "wavelength" in meta:
        ds.metadata["wavelength_A"] = wave
        ds.metadata_sources["wavelength_A"] = "file"
    else:
        ds.metadata["wavelength_A"] = wave
        ds.metadata_sources["wavelength_A"] = "default/unverified"
    for key in ("anode", "radiation", "measurement_datetime"):
        if meta.get(key):
            ds.metadata[key] = meta[key]
            ds.metadata_sources[key] = "file"
    if meta.get("anode") and not meta.get("radiation"):
        ds.metadata["radiation"] = f"{meta['anode']} Kα"
        ds.metadata_sources["radiation"] = "derived from file"
    return ds


def parse_udf(path: Path) -> Dataset:
    text = path.read_text(encoding="utf-8", errors="ignore")
    lines = text.splitlines()
    start = end = step = None
    wavelength = DEFAULT_WAVELENGTH_A
    wavelength_found = False
    raw_start = None
    meta: Dict[str, object] = {"source_format": "PANalytical UDF"}

    for i, line in enumerate(lines):
        st = line.strip()
        if st.startswith("LabdaAlpha1"):
            try:
                wavelength = float(st.split(",")[1]); wavelength_found = True
            except Exception:
                pass
        elif st.startswith("DataAngleRange"):
            try:
                parts = [q.strip(" /,") for q in st.split(",")]
                start, end = float(parts[1]), float(parts[2])
            except Exception:
                pass
        elif st.startswith("ScanStepSize"):
            try:
                parts = [q.strip(" /,") for q in st.split(",")]
                step = float(parts[1])
            except Exception:
                pass
        elif st.startswith("GeneratorVoltage"):
            try: meta["tube_voltage_kV"] = float(st.split(",")[1].strip(" /"))
            except Exception: pass
        elif st.startswith("TubeCurrent"):
            try: meta["tube_current_mA"] = float(st.split(",")[1].strip(" /"))
            except Exception: pass
        elif st.startswith("Anode"):
            try: meta["anode"] = st.split(",", 1)[1].strip(" /,")
            except Exception: pass
        elif st.startswith("DiffrType"):
            try: meta["instrument_family"] = st.split(",", 1)[1].strip(" /,")
            except Exception: pass
        elif st.startswith("FileDateTime"):
            try: meta["measurement_datetime"] = st.split(",", 1)[1].strip(" /,")
            except Exception: pass
        elif st.startswith("ScanType"):
            try: meta["scan_mode"] = st.split(",", 1)[1].strip(" /,")
            except Exception: pass
        elif st.startswith("Title2") and "Goniometer=" in st:
            try:
                desc = st.split(",", 1)[1]
                m = re.search(r"Goniometer=([^;]+(?:\([^)]*\))?)", desc, flags=re.I)
                if m:
                    meta["goniometer"] = m.group(1).strip()
                    if "theta/theta" in m.group(1).lower():
                        meta["geometry"] = "Theta/Theta"
            except Exception:
                pass
        elif st.startswith("RawScan"):
            raw_start = i + 1
            break

    if raw_start is None:
        raise ValueError("Nie znaleziono sekcji RawScan w pliku UDF.")
    vals = []
    for line in lines[raw_start:]:
        for token in line.replace("/", "").replace(",", " ").split():
            try: vals.append(float(token))
            except Exception: pass
    y = np.asarray(vals, dtype=float)
    if start is None or step is None:
        raise ValueError("W pliku UDF brakuje zakresu kątów lub kroku skanu.")
    x = start + np.arange(len(y)) * step
    ds = Dataset(path, path.stem, "", x, y, wavelength)
    for key, value in meta.items():
        if value not in (None, ""):
            ds.metadata[key] = value
            ds.metadata_sources[key] = "file"
    ds.metadata["wavelength_A"] = wavelength
    ds.metadata_sources["wavelength_A"] = "file" if wavelength_found else "default/unverified"
    if meta.get("anode"):
        ds.metadata["radiation"] = f"{meta['anode']} Kα"
        ds.metadata_sources["radiation"] = "derived from file"
    return ds


def supplement_legacy_rd_provenance(dataset: Dataset) -> None:
    """Enrich old projects only when an unchanged source .rd is still available.

    Never replace the project scan axis or its measured intensities. If the source
    is unavailable/changed or inconsistent with the saved axis, do nothing.
    """
    md = dataset.metadata if isinstance(dataset.metadata, dict) else {}
    if md.get('rd_source_checked_2035_1'):
        return
    md['rd_source_checked_2035_1'] = True
    path = Path(dataset.path)
    if path.suffix.lower() != '.rd' or not path.is_file():
        return
    try:
        scan = parse_philips_rd(path)
        if len(scan.x) != len(dataset.x) or not np.allclose(scan.x,dataset.x,rtol=0,atol=1e-4):
            return
        if scan.wavelength_A and scan.wavelength_kalpha2_A:
            md.setdefault('wavelength_kalpha1_A',float(scan.wavelength_A))
            md.setdefault('wavelength_kalpha2_A',float(scan.wavelength_kalpha2_A))
            dataset.metadata_sources.setdefault('wavelength_kalpha2_A','file rechecked (heuristic header scan)')
        md.setdefault('rd_axis_missing_points',int(scan.layout.get('missing_points',0)))
        md.setdefault('rd_axis_anchor',str(scan.layout.get('anchor','start')))
        dataset.metadata_sources.setdefault('rd_axis_missing_points','file rechecked / header comparison')
    except (OSError,ValueError):
        return


def parse_rd(path: Path) -> Dataset:
    """Philips/PANalytical .rd (V3RD). The reader validates the layout it finds and refuses to guess."""
    scan = parse_philips_rd(path)
    wavelength = scan.wavelength_A if scan.wavelength_A else DEFAULT_WAVELENGTH_A
    ds = Dataset(path, path.stem, "", scan.x, scan.y, wavelength)
    ds.metadata["source_format"] = f"Philips RD v{scan.version} (odczyt heurystyczny)"
    ds.metadata_sources["source_format"] = "file"
    ds.metadata["wavelength_A"] = wavelength
    ds.metadata_sources["wavelength_A"] = "file (heurystyka)" if scan.wavelength_A else "default/unverified"
    if scan.anode:
        ds.metadata["anode"] = scan.anode
        ds.metadata_sources["anode"] = "derived from file"
        ds.metadata["radiation"] = f"{scan.anode} Kα"
        ds.metadata_sources["radiation"] = "derived from file"
    if scan.wavelength_A and scan.wavelength_kalpha2_A:
        ds.metadata["wavelength_kalpha1_A"] = float(scan.wavelength_A)
        ds.metadata["wavelength_kalpha2_A"] = float(scan.wavelength_kalpha2_A)
        ds.metadata_sources["wavelength_kalpha1_A"] = "file (heuristic header scan)"
        ds.metadata_sources["wavelength_kalpha2_A"] = "file (heuristic header scan)"
    ds.metadata["rd_axis_missing_points"] = int(scan.layout.get("missing_points", 0))
    ds.metadata["rd_axis_anchor"] = str(scan.layout.get("anchor", "start"))
    ds.metadata_sources["rd_axis_missing_points"] = "file/header comparison"
    notes = list(scan.warnings)
    if scan.wavelength_kalpha2_A:
        notes.append(f"Kα2 z nagłówka: {scan.wavelength_kalpha2_A:.6g} Å.")
    if scan.header_text:
        notes.append(f"Teksty z nagłówka: {scan.header_text}")
    ds.metadata["notes"] = " ".join(notes)
    ds.metadata_sources["notes"] = "derived"
    return ds


def parse_xrdml(path: Path) -> Dataset:
    """Read PANalytical/Malvern Panalytical XRDML files across schema versions.

    XRDML 1.3 and 2.x are supported by matching XML elements by local name
    instead of hard-coding the namespace URI.
    """
    root = ET.parse(path).getroot()

    def local_name(tag):
        return tag.split("}", 1)[-1] if "}" in tag else tag

    def first_descendant(name, parent=None):
        base = parent if parent is not None else root
        for el in base.iter():
            if local_name(el.tag) == name:
                return el
        return None

    def descendants(name):
        return [el for el in root.iter() if local_name(el.tag) == name]

    def child_by_name(parent, name):
        if parent is None:
            return None
        for el in list(parent):
            if local_name(el.tag) == name:
                return el
        return None

    schema_uri = ""
    if root.tag.startswith("{"):
        schema_uri = root.tag[1:].split("}", 1)[0]
    schema_version = "XRDML"
    mver = re.search(r"XRDMeasurement/([0-9.]+)", schema_uri)
    if mver:
        schema_version = f"XRDML {mver.group(1)}"

    meta: Dict[str, object] = {"source_format": schema_version}

    ka1 = first_descendant("kAlpha1")
    wavelength_found = ka1 is not None and bool((ka1.text or "").strip())
    try:
        wavelength = float(ka1.text.strip()) if wavelength_found else DEFAULT_WAVELENGTH_A
    except Exception:
        wavelength = DEFAULT_WAVELENGTH_A; wavelength_found = False

    node = first_descendant("intensities")
    if node is None:
        node = first_descendant("counts")
    if node is None or not (node.text or "").strip():
        raise ValueError(
            "Nie znaleziono danych intensities/counts w pliku XRDML. "
            f"Wykryta przestrzeń nazw XML: {schema_uri or 'brak'}."
        )
    try:
        y = np.asarray([float(v) for v in node.text.split()], dtype=float)
    except Exception as exc:
        raise ValueError(f"Nie udało się odczytać intensywności z pliku XRDML: {exc}")
    if y.size < 2:
        raise ValueError("Plik XRDML zawiera mniej niż dwa punkty intensywności.")

    x = None
    for pos in descendants("positions"):
        axis = (pos.attrib.get("axis", "") or "").lower().replace(" ", "")
        if "2theta" not in axis:
            continue
        start_node = child_by_name(pos, "startPosition")
        end_node = child_by_name(pos, "endPosition")
        list_pos = child_by_name(pos, "listPositions")
        if list_pos is not None and (list_pos.text or "").strip():
            vals = np.asarray([float(v) for v in list_pos.text.split()], dtype=float)
            if vals.size == y.size:
                x = vals
            elif vals.size > 1:
                x = np.linspace(float(vals[0]), float(vals[-1]), y.size)
        elif start_node is not None and end_node is not None and start_node.text and end_node.text:
            x = np.linspace(float(start_node.text), float(end_node.text), y.size)
        if x is not None:
            break
    if x is None:
        raise ValueError("Nie znaleziono pozycji 2θ w pliku XRDML.")
    if x.size != y.size:
        raise ValueError(f"Niezgodna liczba punktów w XRDML: pozycje {x.size}, intensywności {y.size}.")

    # Instrument and acquisition metadata available in many XRDML exports.
    tension = first_descendant("tension")
    current = first_descendant("current")
    anode = first_descendant("anodeMaterial")
    timestamp = first_descendant("startTimeStamp")
    ctime = first_descendant("commonCountingTime")
    app = first_descendant("applicationSoftware")
    scan = first_descendant("scan")
    used_wl = first_descendant("usedWavelength")

    def float_text(el):
        try: return float((el.text or "").strip()) if el is not None else None
        except Exception: return None

    tv = float_text(tension); tc = float_text(current); ct = float_text(ctime)
    if tv is not None: meta["tube_voltage_kV"] = tv
    if tc is not None: meta["tube_current_mA"] = tc
    if anode is not None and (anode.text or "").strip(): meta["anode"] = anode.text.strip()
    if timestamp is not None and (timestamp.text or "").strip(): meta["measurement_datetime"] = timestamp.text.strip()
    if ct is not None: meta["counting_time_s_as_stored"] = ct
    if app is not None and (app.text or "").strip():
        ver = app.attrib.get("version", "")
        meta["application_software"] = (app.text.strip() + (f" {ver}" if ver else "")).strip()
    if scan is not None:
        if scan.attrib.get("mode"): meta["scan_mode"] = scan.attrib.get("mode")
        if scan.attrib.get("scanAxis"): meta["scan_axis"] = scan.attrib.get("scanAxis")
    if used_wl is not None and used_wl.attrib.get("intended"):
        meta["radiation"] = used_wl.attrib.get("intended")
    elif meta.get("anode"):
        meta["radiation"] = f"{meta['anode']} Kα"

    # Comments often contain the goniometer/configuration string.
    entries = []
    for el in descendants("entry"):
        if (el.text or "").strip(): entries.append(el.text.strip())
    for entry in entries:
        if entry.lower().startswith("configuration="):
            meta["configuration"] = entry.split("=", 1)[1].strip()
        if "goniometer=" in entry.lower():
            m = re.search(r"Goniometer=([^;]+(?:\([^)]*\))?)", entry, flags=re.I)
            if m:
                meta["goniometer"] = m.group(1).strip()
                if "theta/theta" in m.group(1).lower():
                    meta["geometry"] = "Theta/Theta"

    ds = Dataset(path, path.stem, "", x, y, wavelength)
    for key, value in meta.items():
        if value not in (None, ""):
            ds.metadata[key] = value
            ds.metadata_sources[key] = "file"
    ds.metadata["wavelength_A"] = wavelength
    ds.metadata_sources["wavelength_A"] = "file" if wavelength_found else "default/unverified"
    return ds


def _leading_numbers(tokens, decimal_comma: bool) -> List[float]:
    """Parse consecutive numeric fields from the start of a line; stop at the first text field."""
    numbers: List[float] = []
    for token in tokens:
        if decimal_comma:
            token = token.replace(",", ".")
        try:
            value = float(token)
        except ValueError:
            break
        if not math.isfinite(value):
            break
        numbers.append(value)
    return numbers


def parse_xy_line(line: str) -> List[float]:
    """Return the leading numeric fields of one ASCII data line (2theta, intensity, ...).

    Supported separators: space/tab, semicolon (decimal point or decimal comma) and comma
    (decimal point). A comma is a field separator when it touches a field boundary
    ("35.00,1000" or "35.00, 1000") and a decimal comma when fields are separated by
    whitespace ("35,02 1000"). Header and comment lines yield an empty list because a data
    line must *start* with numbers.
    """
    text = line.strip().lstrip("\ufeff")
    if not text:
        return []
    if ";" in text or "\t" in text:
        return _leading_numbers([t for t in re.split(r"[;\t\s]+", text) if t], decimal_comma=True)
    if "," not in text:
        return _leading_numbers(text.split(), decimal_comma=False)
    fields = text.split()
    boundary_comma = any(t.startswith(",") or t.endswith(",") for t in fields)
    if len(fields) >= 2 and not boundary_comma:
        numbers = _leading_numbers(fields, decimal_comma=True)  # whitespace separated, decimal comma
        if len(numbers) >= 2:
            return numbers
    return _leading_numbers([t for t in re.split(r"[,\s]+", text) if t], decimal_comma=False)


def parse_csv_generic(path: Path) -> Dataset:
    """Read generic two-column ASCII/CSV/ASC XRD data.

    The first two numeric values on each usable line are interpreted as
    2theta [deg] and intensity. Header/text lines are ignored. Comma, semicolon,
    tab and space separated files are supported (see parse_xy_line). A file whose
    2theta axis strictly decreases is reversed; other orderings are not changed.
    """
    raw = path.read_text(encoding="utf-8-sig", errors="ignore")
    xs, ys = [], []
    for line in raw.splitlines():
        nums = parse_xy_line(line)
        if len(nums) >= 2:
            xs.append(nums[0]); ys.append(nums[1])

    if len(xs) < 5:
        raise ValueError("Plik ASCII/CSV musi zawierać co najmniej 5 wierszy z dwiema kolumnami liczbowymi: "
                         "2θ oraz intensywność (separator: spacja, tabulator, średnik lub przecinek).")
    x = np.asarray(xs, dtype=float); y = np.asarray(ys, dtype=float)
    if len(x) > 1:
        keep = np.ones(len(x), dtype=bool)
        keep[1:] = ~((x[1:] == x[:-1]) & (y[1:] == y[:-1]))
        x, y = x[keep], y[keep]
    reversed_axis = bool(len(x) > 1 and np.all(np.diff(x) < 0))
    if reversed_axis:
        x, y = x[::-1].copy(), y[::-1].copy()

    ds = Dataset(path, path.stem, "", x, y, DEFAULT_WAVELENGTH_A)
    ds.metadata["source_format"] = path.suffix.upper().lstrip(".") or "ASCII"
    ds.metadata_sources["source_format"] = "file"
    ds.metadata["wavelength_A"] = DEFAULT_WAVELENGTH_A
    ds.metadata_sources["wavelength_A"] = "default/unverified"
    note = "Plik tekstowy zawiera wyłącznie dane liczbowe dyfraktogramu; parametry aparatury nie są w nim zapisane."
    if reversed_axis:
        note += " Oś 2θ w pliku malała, więc kolejność punktów odwrócono (wartości bez zmian)."
    ds.metadata["notes"] = note
    ds.metadata_sources["notes"] = "derived"
    return ds


def infer_material_type(filename: str) -> str:
    t = filename.lower()
    if "spoina" in t or "weld" in t:
        if "monel" in t or "3a" in t or "3b" in t:
            return "Spoina różnoimienna 316L/CuNi30Mn1Fe – 67ERCuNi filler"
        if "316l" in t or "4a" in t or "4b" in t:
            return "Spoina różnoimienna 316L/CuNi30Mn1Fe – 316L filler"
        return "Spoina różnoimienna 316L/CuNi30Mn1Fe"
    if "67ercuni" in t or "ercuni" in t:
        return "67ERCuNi filler wire"
    if "cuni30" in t or "cunifer" in t:
        return "CuNi30Mn1Fe base material"
    if "x2crnimo17-12-2" in t or "1.4404" in t or "316l" in t:
        return "X2CrNiMo17-12-2 / 1.4404 / AISI 316L"
    if "x5crni18-10" in t or "1.4301" in t or "aisi304" in t or "aisi_304" in t:
        return "X5CrNi18-10 / EN 1.4301 austenitic stainless steel"
    if re.search(r"(^|[^a-z])ti([^a-z]|$)", t) or "titan" in t:
        return "Titanium alloy / Ti-6Al-4V"
    if "zircaloy" in t or "zircon" in t or re.search(r"(^|[^a-z])zr([^a-z]|$)", t):
        return "Zirconium / zirconium alloy"
    if "magnesium" in t or "az31" in t or "az91" in t or re.search(r"(^|[^a-z])mg([^a-z]|$)", t):
        return "Magnesium alloy"
    if "aluminium" in t or "aluminum" in t:
        return "Aluminium alloy"
    if "copper" in t or "brass" in t or "bronze" in t:
        return "Copper alloy"
    return ""


def load_dataset(path: Path) -> Dataset:
    suffix = path.suffix.lower()
    if suffix in (".cif",):
        ds = parse_cif(path)
    elif suffix in (".udf",):
        ds = parse_udf(path)
    elif suffix in (".rd",):
        ds = parse_rd(path)
    elif suffix in (".xrdml", ".xml"):
        ds = parse_xrdml(path)
    elif suffix in (".csv", ".txt", ".dat", ".asc"):
        ds = parse_csv_generic(path)
    else:
        raise ValueError(f"Nieobsługiwany typ pliku: {suffix}")
    # Safety: never treat a filename as verified material chemistry.  Keep a possible
    # filename-based guess only as an explicitly UNVERIFIED suggestion for the wizard/UI.
    # It is not copied into material_type and therefore cannot bias phase ranking unless
    # the user consciously accepts/enters it.
    suggestion = infer_material_type(path.stem)
    ds.material_type = ""
    if suggestion:
        ds.metadata["filename_material_suggestion"] = suggestion
        ds.metadata_sources["filename_material_suggestion"] = "unverified filename suggestion (not used for analysis)"
    return ds


# --------------------------- analysis ---------------------------

BASELINE_REFERENCE_STEP_DEG = 0.05  # scan step of the 316L / CuNi / Monel reference scans for which lam=2e6 is known to work


def baseline_asls(y: np.ndarray, lam=2e6, p=0.01, niter=10, x=None) -> np.ndarray:
    """Asymmetric least-squares baseline (Eilers & Boelens).

    The penalty acts on second differences between *points*, so a fixed lam makes the baseline
    stiffer or more flexible in degrees depending on the scan step (stiffness ~ step^-4): with a
    0.005 deg step the baseline followed broad peaks and removed them. When the 2theta axis x is
    given, lam is rescaled by (0.05 / median step)^4 so that the smoothness in degrees is the same
    for every step (scans with step 0.05 deg, the reference, are unchanged); without x the legacy
    per-point behaviour is kept.
    """
    y = np.asarray(y, dtype=float)
    L = len(y)
    if x is not None:
        axis = np.asarray(x, dtype=float)
        if axis.ndim == 1 and len(axis) == L and L > 1:
            step = float(np.median(np.diff(axis)))
            if math.isfinite(step) and step > 0:
                lam = float(np.clip(lam * (BASELINE_REFERENCE_STEP_DEG / step) ** 4, 1e2, 1e10))
    D = sparse.diags([1.0, -2.0, 1.0], [0, -1, -2], shape=(L, L-2), format="csc")
    w = np.ones(L)
    for _ in range(niter):
        W = sparse.spdiags(w, 0, L, L).tocsc()
        Z = W + lam * (D @ D.T)
        z = spsolve(Z, w * y)
        w = p * (y > z) + (1-p) * (y < z)
    return np.asarray(z)


def bragg_d(two_theta_deg: float, wavelength_A: float) -> float:
    theta = math.radians(two_theta_deg / 2.0)
    return wavelength_A / (2.0 * math.sin(theta))


def d_to_2theta(d_A: float, wavelength_A: float) -> Optional[float]:
    v = wavelength_A / (2.0 * d_A)
    if v <= 0 or v >= 1:
        return None
    return math.degrees(2.0 * math.asin(v))


def trusted_material_context(dataset: Dataset) -> str:
    """Return only user-confirmed material context used for ranking/naming.

    The sample filename/name is deliberately excluded.  Material/grade text is trusted only
    when the user has saved Sample setup / chemistry or explicitly confirmed a built-in preset.
    This prevents a misleading filename from steering phase identification.
    """
    md = dataset.metadata if isinstance(dataset.metadata, dict) else {}
    srcs = dataset.metadata_sources if isinstance(dataset.metadata_sources, dict) else {}
    user_confirmed = any(
        "user" in str(srcs.get(k, "") or "").lower()
        for k in ("material_type", "chemical_composition", "allowed_elements", "expected_matrix_structure")
    ) or bool(str(md.get("material_preset_confirmed", "") or "").strip())
    if not user_confirmed:
        return ""
    return " ".join([
        str(dataset.material_type or ""),
        str(md.get("base_material_A", "")),
        str(md.get("base_material_B", "")),
        str(md.get("filler_material", "")),
    ]).strip()


def material_bias(phase_name: str, category: str, material_text: str) -> float:
    t = material_text.lower()
    n = phase_name.lower()
    c = (category or "").lower()
    bonus = 0.0
    if "x5crni18-10" in t or "1.4301" in t or "aisi 304" in t or "aisi304" in t:
        if "x5crni18-10" in n or "1.4301" in n: bonus += 0.55
        elif "austen" in n or "ferrite" in n or "steel" in c: bonus += 0.25
    elif "316" in t or "x2crnimo17-12-2" in t or "1.4404" in t or "stainless" in t or "stal" in t:
        if "x2crnimo17-12-2" in n or "1.4404" in n or "316l" in n: bonus += 0.55
        elif "316" in n or "austen" in n or "ferrite" in n or "steel" in c: bonus += 0.25
    if "cuni" in t or "cu-ni" in t or "monel" in t or "ercuni" in t:
        if ("67ercuni" in t or "ercuni" in t) and "ercuni" in n:
            bonus += 0.55
        elif ("cuni30" in t or "cuni30mn1fe" in t) and "cuni30" in n:
            bonus += 0.55
        elif "cu-ni" in n or "cuni" in n or "ercuni" in n or "cu-ni" in c:
            bonus += 0.25
    if "ti-6al-4v" in t or "ti6al4v" in t or "titanium" in t:
        if "ti-6al-4v" in n: bonus += 0.45
        elif "titanium" in n or "titanium" in c: bonus += 0.20
    if "aluminium" in t or "aluminum" in t:
        if "aluminium" in n or "α-al" in n or "aluminium" in c: bonus += 0.30
    if "copper" in t or "brass" in t or "bronze" in t or "cuni" in t:
        if "copper" in c or "brass" in c or "bronze" in c or "cu-ni" in c: bonus += 0.30
    if "zirconium" in t or "zircaloy" in t:
        if "zirconium" in n or "zirconium" in c: bonus += 0.35
    if "magnesium" in t or "az31" in t or "az91" in t:
        if "magnesium" in n or "magnesium" in c: bonus += 0.35
    return bonus


def dataset_allowed_elements(dataset: Dataset) -> set[str]:
    raw = dataset.metadata.get("allowed_elements", "") if isinstance(dataset.metadata, dict) else ""
    if isinstance(raw, (list, tuple, set)):
        return {str(x).strip().capitalize() for x in raw if str(x).strip()}
    return {t[0].upper() + t[1:].lower() for t in re.split(r"[,;\s]+", str(raw).strip()) if re.fullmatch(r"[A-Za-z]{1,2}", t)}


EXPECTED_MATRIX_ITEMS = [
    ("Automatycznie na podstawie zatwierdzonej chemii / materiału", "auto"),
    ("Nieznana / bez preferencji struktury", "unknown"),
    ("γ-austenit – osnowa Fe, FCC", "gamma_fcc"),
    ("Ferryt / martenzyt – osnowa Fe, BCC", "ferrite_bcc"),
    ("Roztwór stały Cu–Ni – FCC", "cuni_fcc"),
    ("α-Ti – HCP", "alpha_ti_hcp"),
    ("β-Ti – BCC", "beta_ti_bcc"),
    ("Osnowa aluminium – FCC", "al_fcc"),
    ("FCC – ogólnie", "fcc"),
    ("BCC – ogólnie", "bcc"),
    ("HCP – ogólnie", "hcp"),
]


MATERIAL_SETUP_PRESETS = {
    "X5CrNi18-10 / EN 1.4301 – austenitic stainless steel": {
        "material": "X5CrNi18-10 / EN 1.4301 austenitic stainless steel",
        "elements": ["Fe","Cr","Ni","Mn","Si"], "matrix": "gamma_fcc",
    },
    "X2CrNiMo17-12-2 / EN 1.4404 / AISI 316L": {
        "material": "X2CrNiMo17-12-2 / 1.4404 / AISI 316L",
        "elements": ["Fe","Cr","Ni","Mo","Mn","Si"], "matrix": "gamma_fcc",
    },
    "Ti-6Al-4V": {
        "material": "Ti-6Al-4V", "elements": ["Ti","Al","V","Fe"], "matrix": "auto",
    },
    "CuNi30Mn1Fe": {
        "material": "CuNi30Mn1Fe", "elements": ["Cu","Ni","Mn","Fe"], "matrix": "cuni_fcc",
    },
    "67ERCuNi / ERCuNi filler": {
        "material": "67ERCuNi / ERCuNi filler wire", "elements": ["Cu","Ni","Fe","Mn"], "matrix": "cuni_fcc",
    },
}


def _material_priority_preset(dataset: Dataset) -> Dict[str, float]:
    """Return a broad material prior ONLY after explicit user confirmation of a preset.

    A filename or free-text material name is never enough to activate this prior.  The flag
    ``material_preset_confirmed`` is written only when the user presses *Zastosuj preset* in the
    Kreator konfiguracji próbki and then saves the setup.  Exact/range/order chemistry remains the
    preferred source whenever available.
    """
    md = dataset.metadata if isinstance(dataset.metadata, dict) else {}
    if not str(md.get("material_preset_confirmed", "") or "").strip():
        return {}
    t = trusted_material_context(dataset).lower()
    allowed = dataset_allowed_elements(dataset)

    def filt(vals):
        vals = {e: float(v) for e, v in vals.items() if (not allowed or e in allowed)}
        total = sum(vals.values())
        return {e: v/total for e, v in vals.items()} if total > 0 else {}

    is_x5_1810 = ("x5crni18-10" in t or "1.4301" in t or "aisi 304" in t or "aisi304" in t)
    is_316 = ("316" in t or "x2crnimo17-12-2" in t or "1.4404" in t)
    is_cuni = ("cuni" in t or "ercuni" in t or "monel" in t or "cunifer" in t)
    is_weld = ("weld" in t or "spoina" in t or "joint" in t)
    if is_weld and is_316 and is_cuni:
        # Spoina różnoimienna metal can contain substantial contributions from both sides/filler.
        # Do not pretend that one matrix element is known a priori.
        return filt({"Fe": 34.0, "Cu": 34.0, "Ni": 18.0, "Cr": 9.0, "Mo": 2.0, "Mn": 2.0, "Si": 1.0})
    if is_x5_1810:
        return filt({"Fe": 70.5, "Cr": 18.0, "Ni": 9.5, "Mn": 1.5, "Si": 0.5})
    if is_316:
        return filt({"Fe": 66.0, "Cr": 17.0, "Ni": 11.0, "Mo": 2.2, "Mn": 1.5, "Si": 0.7})
    if "cuni30" in t or "cuni30mn1fe" in t:
        return filt({"Cu": 67.5, "Ni": 30.5, "Mn": 1.0, "Fe": 1.0})
    if "67ercuni" in t or "ercuni" in t or "monel" in t:
        return filt({"Cu": 67.0, "Ni": 30.0, "Fe": 2.0, "Mn": 1.0})
    if "ti6al4v" in t or "ti-6al-4v" in t:
        return filt({"Ti": 89.5, "Al": 6.0, "V": 4.0, "Fe": 0.5})
    return {}


def dataset_element_priority_weights(dataset: Dataset) -> Dict[str, float]:
    """Return chemistry priorities for phase ranking.

    Exact/range/order chemistry entered in the Kreator konfiguracji próbki has the highest priority.
    A broad alloy prior is used only when the user explicitly confirmed a built-in material
    preset in Sample setup.  A plain element list is deliberately treated as weak, equal
    information; filenames never create a chemical prior.
    """
    rows = dataset.metadata.get("chemical_composition", []) if isinstance(dataset.metadata, dict) else []
    mode = str(dataset.metadata.get("composition_mode", "elements")) if isinstance(dataset.metadata, dict) else "elements"
    vals: Dict[str, float] = {}

    if isinstance(rows, list) and rows:
        if mode == "order":
            n = len(rows)
            # A steeper decay than the old linear weighting makes the first user-entered
            # element genuinely represent the matrix/dominant constituent.
            for i, rec in enumerate(rows):
                if isinstance(rec, dict) and rec.get("element") in ELEMENT_Z:
                    vals[str(rec["element"])] = float((n - i) ** 1.7)
        elif mode in ("exact", "ranges"):
            numeric_sum = 0.0; balance_els = []
            for rec in rows:
                if not isinstance(rec, dict):
                    continue
                el = str(rec.get("element", "")); raw = str(rec.get("value", "")).strip().lower()
                if el not in ELEMENT_Z:
                    continue
                if raw in ("balance", "bal", "rest", "remainder"):
                    balance_els.append(el); continue
                v = None
                try:
                    v = float(str(rec.get("value", "")).replace(",", "."))
                except Exception:
                    try:
                        lo = float(str(rec.get("min", "")).replace(",", "."))
                        hi = float(str(rec.get("max", "")).replace(",", "."))
                        v = (lo + hi) / 2.0
                    except Exception:
                        v = None
                if v is not None and v >= 0:
                    vals[el] = max(0.001, float(v)); numeric_sum += float(v)
            if balance_els:
                remain = max(1.0, 100.0 - numeric_sum)
                for el in balance_els:
                    vals[el] = remain / len(balance_els)
        elif mode == "elements":
            preset = _material_priority_preset(dataset)
            if preset:
                return preset
            for rec in rows:
                if isinstance(rec, dict) and rec.get("element") in ELEMENT_Z:
                    vals[str(rec["element"])] = 1.0

    if not vals:
        preset = _material_priority_preset(dataset)
        if preset:
            return preset
        els = sorted(dataset_allowed_elements(dataset), key=lambda e: ELEMENT_Z.get(e, 999))
        if not els:
            return {}
        vals = {e: 1.0 for e in els}

    total = sum(vals.values())
    return {e: v / total for e, v in vals.items()} if total > 0 else {}


def dataset_dominant_element(dataset: Dataset) -> tuple[Optional[str], float, bool]:
    """Return (element, normalized weight, strong_evidence).

    strong_evidence means that dominance came from exact/range/order chemistry or an explicitly
    confirmed material preset, rather than from an unordered element list.
    """
    weights = dataset_element_priority_weights(dataset)
    if not weights:
        return None, 0.0, False
    dom = max(weights, key=weights.get)
    md = dataset.metadata if isinstance(dataset.metadata, dict) else {}
    mode = str(md.get("composition_mode", "elements"))
    rows = md.get("chemical_composition", [])
    strong = mode in {"exact", "ranges", "order"} and isinstance(rows, list) and bool(rows)
    if not strong and _material_priority_preset(dataset):
        text = trusted_material_context(dataset).lower()
        dissimilar = (("weld" in text or "spoina" in text or "joint" in text)
                      and ("316" in text or "1.4404" in text)
                      and ("cuni" in text or "ercuni" in text or "monel" in text))
        strong = not dissimilar
    return dom, float(weights[dom]), strong


def dataset_material_signature(dataset: Dataset) -> str:
    """Infer a broad alloy family from user chemistry/material metadata, never from XRD alone.

    The result is used for contextual phase naming and ranking.  It does not force a phase:
    the diffraction pattern still has to support FCC/BCC/HCP reflections.
    """
    md = dataset.metadata if isinstance(dataset.metadata, dict) else {}
    t = trusted_material_context(dataset).lower()
    if "x5crni18-10" in t or "1.4301" in t or "aisi 304" in t or "aisi304" in t:
        return "fe_cr_ni_austenitic"
    if "316" in t or "x2crnimo17-12-2" in t or "1.4404" in t:
        return "fe_cr_ni_austenitic"
    if "ti6al4v" in t or "ti-6al-4v" in t:
        return "ti64"
    weights = dataset_element_priority_weights(dataset)
    if weights:
        fe, cr, ni = weights.get("Fe", 0.0), weights.get("Cr", 0.0), weights.get("Ni", 0.0)
        ti, al, v = weights.get("Ti", 0.0), weights.get("Al", 0.0), weights.get("V", 0.0)
        if fe >= 0.35 and cr >= 0.08 and ni >= 0.04:
            return "fe_cr_ni_austenitic"
        if ti >= 0.45 and al >= 0.05 and v >= 0.02:
            return "ti64"
    # An unordered user-confirmed element list is weak quantitative information, but its
    # membership is still useful for broad *context* once the diffraction itself has
    # established FCC/HCP/BCC.  This never comes from the filename.
    allowed = dataset_allowed_elements(dataset)
    if {"Ti","Al","V"}.issubset(allowed):
        return "ti64_elements"
    if {"Fe","Cr","Ni","Cu"}.issubset(allowed):
        return "mixed_fe_cu_weld"
    if {"Fe","Cr","Ni"}.issubset(allowed):
        return "fe_cr_ni_elements"
    if {"Cu","Ni"}.issubset(allowed) and "Cr" not in allowed:
        return "cuni_elements"
    return ""


def dataset_expected_matrix_code(dataset: Dataset) -> tuple[str, bool]:
    """Return (expected matrix code, explicit_user_choice)."""
    md = dataset.metadata if isinstance(dataset.metadata, dict) else {}
    code = str(md.get("expected_matrix_structure", "auto") or "auto")
    if code != "auto":
        return code, True

    t = trusted_material_context(dataset).lower()
    is_x5_1810 = ("x5crni18-10" in t or "1.4301" in t or "aisi 304" in t or "aisi304" in t)
    is_316 = ("316" in t or "x2crnimo17-12-2" in t or "1.4404" in t)
    is_cuni = ("cuni" in t or "ercuni" in t or "monel" in t or "cunifer" in t)
    is_weld = ("weld" in t or "spoina" in t or "joint" in t)
    if is_weld and is_316 and is_cuni:
        return "fcc", False
    if is_x5_1810 or is_316:
        return "gamma_fcc", False
    if is_cuni:
        return "cuni_fcc", False
    if "ti6al4v" in t or "ti-6al-4v" in t:
        return "ti_dual", False
    if "aluminium" in t or "aluminum" in t:
        return "al_fcc", False
    sig = dataset_material_signature(dataset)
    # Chemistry alone may describe duplex/ferritic as well as austenitic Fe-Cr-Ni steel,
    # so it is used for contextual naming but does not force an FCC structure preference.
    if sig == "ti64":
        return "ti_dual", False
    return "unknown", False


def phase_structure_factor(dataset: Dataset, structure: str = "", space_group: str = "", sg_no=None) -> tuple[float, str]:
    """Soft compatibility factor for an expected matrix crystal structure."""
    code, explicit = dataset_expected_matrix_code(dataset)
    if code in {"", "unknown", "auto"}:
        return 1.0, "no preference"
    struct = str(structure or "").upper().strip()
    if not struct:
        struct = cod_structure_short(space_group, sg_no).upper()
    targets = {
        "gamma_fcc": {"FCC"}, "cuni_fcc": {"FCC"}, "al_fcc": {"FCC"}, "fcc": {"FCC"},
        "ferrite_bcc": {"BCC"}, "beta_ti_bcc": {"BCC"}, "bcc": {"BCC"},
        "alpha_ti_hcp": {"HCP"}, "hcp": {"HCP"}, "ti_dual": {"HCP", "BCC"},
    }.get(code, set())
    if not targets:
        return 1.0, "no preference"
    if struct in targets:
        return (1.12 if explicit else 1.08), f"{struct} matches expected matrix"
    # Explicit user choices are intentionally stronger than an automatic material-name hint.
    return (0.45 if explicit else 0.72), f"{struct or 'unknown'} differs from expected {'/'.join(sorted(targets))}"


def phase_chemistry_factor(dataset: Dataset, phase_elements: set[str]) -> tuple[float, float, str]:
    """Return (multiplicative factor, captured bulk weight, explanation).

    This is the key safeguard against false matrix assignments such as 'Cr FCC' for 316L.
    A minor-element-only COD entry may remain visible as a secondary candidate, but it is
    strongly down-ranked when it attempts to explain the dominant whole-pattern peaks.
    """
    elems = {str(e) for e in (phase_elements or set()) if str(e)}
    weights = dataset_element_priority_weights(dataset)
    if not weights or not elems:
        return 1.0, 0.0, "chemistry unavailable"

    captured = float(sum(weights.get(e, 0.0) for e in elems))
    factor = 0.25 + 0.75 * math.sqrt(max(0.0, min(1.0, captured)))
    dom, dom_w, strong_dom = dataset_dominant_element(dataset)
    reasons = [f"captures {captured*100:.1f}% of bulk-priority weight"]

    if dom and strong_dom:
        if dom in elems:
            factor *= 1.05
            reasons.append(f"contains dominant {dom}")
        else:
            # Missing a known matrix element is a major warning for a whole-pattern matrix fit.
            miss = 0.22 if dom_w >= 0.45 else (0.35 if dom_w >= 0.30 else 0.55)
            factor *= miss
            reasons.append(f"missing dominant {dom}")

        # A single minor constituent should not become the bulk matrix solely through a
        # coincident lattice parameter (the exact failure seen for Cr FCC in 316L).
        if len(elems) == 1:
            only = next(iter(elems))
            if only != dom and weights.get(only, 0.0) < 0.35:
                factor *= 0.25
                reasons.append(f"single minor constituent {only}")

    # If a multi-element COD record claims to represent the alloy matrix, missing a major
    # alloying element is meaningful.  This specifically separates Fe-Cr-Ni stainless steel
    # from geometrically similar Fe-Mn-Ni records.  A pure dominant-element record (e.g. Ti
    # used as an HCP structural proxy for Ti-6Al-4V) is penalized only mildly, because COD
    # often stores an end-member structure rather than the exact substitutional alloy.
    major = [(e, w) for e, w in weights.items() if w >= 0.07]
    missing_major = [(e, w) for e, w in major if e not in elems]
    if missing_major:
        if len(elems) == 1 and dom in elems:
            factor *= 0.88
            reasons.append("dominant-element structural proxy; missing alloying elements " + ",".join(e for e,_ in missing_major))
        else:
            for e, w in missing_major:
                factor *= max(0.45, 1.0 - 2.0*w)
            reasons.append("missing major alloying element(s) " + ",".join(e for e,_ in missing_major))

    # For an unordered but explicitly user-confirmed element list, use alloy-system membership
    # as a weak qualitative restraint. This is deliberately softer than exact percentages, but
    # prevents a single allowed minor element from becoming the whole matrix merely because its
    # lattice happens to coincide.
    md = dataset.metadata if isinstance(dataset.metadata, dict) else {}
    if str(md.get("composition_mode", "elements")) == "elements":
        allowed = dataset_allowed_elements(dataset)
        if {"Fe","Cr","Ni"}.issubset(allowed) and "Cu" not in allowed:
            nsys = len(elems & {"Fe","Cr","Ni"})
            if nsys >= 2:
                factor *= 1.08; reasons.append("matches Fe-Cr-Ni alloy system")
            elif len(elems) == 1 and next(iter(elems)) in {"Cr","Ni","Mn","Mo","Si"}:
                factor *= 0.38; reasons.append("single alloying-element proxy for Fe-Cr-Ni matrix")
            elif "Fe" in elems:
                factor *= 0.82; reasons.append("Fe end-member proxy")
        elif {"Cu","Ni"}.issubset(allowed) and "Cr" not in allowed:
            if {"Cu","Ni"}.issubset(elems):
                factor *= 1.10; reasons.append("matches Cu-Ni alloy system")
            elif len(elems)==1 and ("Cu" in elems or "Ni" in elems):
                factor *= 0.82; reasons.append("Cu/Ni end-member proxy")

    # Reward chemically richer solid-solution records when they capture the principal alloying
    # elements, without requiring every trace element to appear in the COD formula.
    top = sorted(weights.items(), key=lambda kv: kv[1], reverse=True)[:3]
    top_cov = sum(w for e, w in top if e in elems) / max(1e-12, sum(w for _, w in top))
    factor *= 0.92 + 0.16 * top_cov
    factor = max(0.02, min(1.15, float(factor)))
    return factor, captured, "; ".join(reasons)


def contextual_cod_phase_name(dataset: Dataset, formula: str, current_plot_name: str,
                              space_group: str = "", sg_no=None) -> str:
    """Translate a chemically compatible COD structural reference into a sample-context phase label.

    The COD formula/ID is still preserved in the result table and detail dialog.  This only changes
    the human-facing phase label after chemistry-aware ranking has selected the candidate.
    """
    code, _explicit = dataset_expected_matrix_code(dataset)
    struct = cod_structure_short(space_group, sg_no)
    elems = formula_elements(formula)
    allowed = dataset_allowed_elements(dataset)
    signature = dataset_material_signature(dataset)
    if (code == "gamma_fcc" or signature in {"fe_cr_ni_austenitic","fe_cr_ni_elements"}) and struct == "FCC" and ("Fe" in elems or {"Fe","Cr","Ni"}.issubset(allowed)):
        return "γ-austenite FCC"
    if (code == "cuni_fcc" or signature == "cuni_elements") and struct == "FCC" and ({"Cu", "Ni"} & elems) and {"Cu", "Ni"}.issubset(allowed | elems):
        return "Cu–Ni FCC"
    if signature == "mixed_fe_cu_weld" and struct == "FCC":
        return "FCC solid solution (Fe–Cr–Ni–Cu)"
    if (code in {"alpha_ti_hcp", "ti_dual"} or signature == "ti64") and struct == "HCP" and "Ti" in elems:
        return "α-Ti HCP"
    if (code in {"beta_ti_bcc", "ti_dual"} or signature == "ti64") and struct == "BCC" and "Ti" in elems:
        return "β-Ti BCC"
    if code == "al_fcc" and struct == "FCC" and "Al" in elems:
        return "α-Al FCC"
    return current_plot_name


def canonical_hkl_for_structure(hkl: str, structure: str) -> str:
    """Canonicalize symmetry-equivalent cubic Miller labels for publication.

    COD/gemmi may return e.g. 002, 022 or 113.  For cubic materials these are the
    same reflection families conventionally written 200, 220 and 311.
    """
    raw = str(hkl or "").strip().strip("()[]{} ")
    st = str(structure or "").upper()
    if st not in {"FCC", "BCC", "CUBIC"}:
        return raw
    vals = None
    if re.fullmatch(r"-?\d+[, ]+-?\d+[, ]+-?\d+", raw):
        try:
            vals = [abs(int(v)) for v in re.split(r"[, ]+", raw)]
        except Exception:
            vals = None
    elif re.fullmatch(r"\d{3}", raw):
        vals = [int(raw[0]), int(raw[1]), int(raw[2])]
    if vals is None or len(vals) != 3:
        return raw
    vals.sort(reverse=True)
    return "".join(str(v) for v in vals)


def _greedy_unique_peak_matches(obs_peaks, refs, tolerance_deg: float):
    """One-to-one matching between measured peaks and reference reflections."""
    pairs = []
    for oi, p in enumerate(obs_peaks):
        for ri, ref in enumerate(refs):
            delta = abs(float(p.two_theta) - float(ref[0]))
            if delta <= tolerance_deg:
                pairs.append((delta, oi, ri))
    pairs.sort(key=lambda q: q[0])
    used_o, used_r, matches = set(), set(), []
    for delta, oi, ri in pairs:
        if oi in used_o or ri in used_r:
            continue
        used_o.add(oi); used_r.add(ri)
        matches.append((oi, ri, delta))
    return matches


def _parse_simple_hkl(hkl_text: str):
    """Parse the simple non-negative hkl labels used by the built-in metal library."""
    t = str(hkl_text or "").strip().replace("(", "").replace(")", "")
    if any(sep in t for sep in [",", " ", ";"]):
        parts = [q for q in re.split(r"[,;\s]+", t) if q]
        if len(parts) == 3:
            try:
                return tuple(int(q) for q in parts)
            except Exception:
                return None
    if re.fullmatch(r"\d{3}", t):
        return int(t[0]), int(t[1]), int(t[2])
    return None


def _refs_for_cubic_hkls(hkl_rows, a_A: float, wavelength_A: float, xmin: float, xmax: float):
    refs = []
    for hkl, rel in hkl_rows:
        h, k, l = hkl
        d = cubic_d(float(a_A), h, k, l)
        tt = d_to_2theta(d, wavelength_A)
        if tt is not None and xmin - 0.5 <= tt <= xmax + 0.5:
            refs.append((float(tt), f"{h}{k}{l}", float(d), float(rel)))
    refs.sort(key=lambda r: r[0])
    return refs


def _fit_cubic_lattice(pattern_peaks, hkl_rows, wavelength_A: float, a_ref: float,
                       xmin: float, xmax: float, tolerance_deg: float):
    """Fit cubic a from the measured peak positions while preserving allowed hkl families.

    The database a_ref is only an initial search seed. Candidate lattice parameters are
    generated directly from measured d-spacings and allowed hkl values, then the whole
    pattern is matched one-to-one and a is refined from all matched reflections.
    """
    if not hkl_rows or not a_ref or not math.isfinite(float(a_ref)):
        return None
    a_ref = float(a_ref)
    # ±7 % safely covers normal alloying/temperature/stress shifts while preventing
    # a low-order measured peak from being mapped onto a completely different hkl.
    lo, hi = 0.93 * a_ref, 1.07 * a_ref
    low_order = sorted(hkl_rows, key=lambda q: sum(v*v for v in q[0]))[:12]
    candidates = [a_ref]
    for p in pattern_peaks[:15]:
        for hkl, _rel in low_order:
            n = sum(v*v for v in hkl)
            if n <= 0:
                continue
            a = float(p.d_A) * math.sqrt(n)
            if lo <= a <= hi:
                candidates.append(a)
    # de-duplicate candidates without losing physically meaningful shifts
    candidates = sorted({round(a, 6) for a in candidates})
    best = None
    for a in candidates:
        refs = _refs_for_cubic_hkls(hkl_rows, a, wavelength_A, xmin, xmax)
        matches = _greedy_unique_peak_matches(pattern_peaks, refs, tolerance_deg)
        if not matches:
            continue
        rms = float(np.sqrt(np.mean([m[2]**2 for m in matches])))
        icov = sum(pattern_peaks[m[0]].rel_intensity for m in matches)
        key = (len(matches), icov, -rms)
        if best is None or key > best[0]:
            best = (key, a, refs, matches, rms)
    if best is None:
        return None

    a = float(best[1])
    # Refine twice from matched d-spacings. Weight strong measured peaks more heavily.
    for _ in range(2):
        refs = _refs_for_cubic_hkls(hkl_rows, a, wavelength_A, xmin, xmax)
        matches = _greedy_unique_peak_matches(pattern_peaks, refs, tolerance_deg)
        if not matches:
            break
        vals, weights = [], []
        for oi, ri, _delta in matches:
            hkl = _parse_simple_hkl(refs[ri][1])
            if hkl is None:
                continue
            n = sum(v*v for v in hkl)
            vals.append(float(pattern_peaks[oi].d_A) * math.sqrt(n))
            weights.append(max(1.0, float(pattern_peaks[oi].rel_intensity)))
        if vals:
            a = float(np.average(np.asarray(vals), weights=np.asarray(weights)))
    refs = _refs_for_cubic_hkls(hkl_rows, a, wavelength_A, xmin, xmax)
    matches = _greedy_unique_peak_matches(pattern_peaks, refs, tolerance_deg)
    rms = float(np.sqrt(np.mean([m[2]**2 for m in matches]))) if matches else 9.0
    return refs, matches, rms, {"model": "cubic", "a_A": a}


def _refs_for_hcp_hkls(hkl_rows, a_A: float, c_A: float, wavelength_A: float, xmin: float, xmax: float):
    refs = []
    for hkl, rel in hkl_rows:
        h, k, l = hkl
        try:
            d = hex_d(float(a_A), float(c_A), h, k, l)
        except Exception:
            continue
        tt = d_to_2theta(d, wavelength_A)
        if tt is not None and xmin - 0.5 <= tt <= xmax + 0.5:
            refs.append((float(tt), f"{h}{k}{l}", float(d), float(rel)))
    refs.sort(key=lambda r: r[0])
    return refs


def _fit_hcp_lattice(pattern_peaks, hkl_rows, wavelength_A: float, a_ref: float, c_ref: float,
                     xmin: float, xmax: float, tolerance_deg: float):
    """Fit HCP a and c with a coarse seed grid followed by weighted linear refinement."""
    if not hkl_rows or not a_ref or not c_ref:
        return None
    a_ref, c_ref = float(a_ref), float(c_ref)
    best = None
    # Coarse ±4% grid is used only to establish the reflection correspondence.
    for sa in np.linspace(0.96, 1.04, 9):
        for sc in np.linspace(0.96, 1.04, 9):
            a, c = a_ref * float(sa), c_ref * float(sc)
            refs = _refs_for_hcp_hkls(hkl_rows, a, c, wavelength_A, xmin, xmax)
            matches = _greedy_unique_peak_matches(pattern_peaks, refs, tolerance_deg)
            if not matches:
                continue
            rms = float(np.sqrt(np.mean([m[2]**2 for m in matches])))
            icov = sum(pattern_peaks[m[0]].rel_intensity for m in matches)
            key = (len(matches), icov, -rms)
            if best is None or key > best[0]:
                best = (key, a, c, refs, matches, rms)
    if best is None:
        return None
    a, c = float(best[1]), float(best[2])
    for _ in range(2):
        refs = _refs_for_hcp_hkls(hkl_rows, a, c, wavelength_A, xmin, xmax)
        matches = _greedy_unique_peak_matches(pattern_peaks, refs, tolerance_deg)
        X, Y, W = [], [], []
        for oi, ri, _delta in matches:
            hkl = _parse_simple_hkl(refs[ri][1])
            if hkl is None:
                continue
            h, k, l = hkl
            A = (4.0/3.0) * (h*h + h*k + k*k)
            B = l*l
            X.append([A, B])
            Y.append(1.0 / (float(pattern_peaks[oi].d_A)**2))
            W.append(max(1.0, float(pattern_peaks[oi].rel_intensity)))
        if len(X) >= 2:
            X = np.asarray(X, float); Y = np.asarray(Y, float); W = np.asarray(W, float)
            sw = np.sqrt(W)[:, None]
            try:
                coef, *_ = np.linalg.lstsq(X * sw, Y * sw[:,0], rcond=None)
                if coef[0] > 0 and coef[1] > 0:
                    aa, cc = 1.0/math.sqrt(float(coef[0])), 1.0/math.sqrt(float(coef[1]))
                    if 0.90*a_ref <= aa <= 1.10*a_ref and 0.90*c_ref <= cc <= 1.10*c_ref:
                        a, c = aa, cc
            except Exception:
                pass
    refs = _refs_for_hcp_hkls(hkl_rows, a, c, wavelength_A, xmin, xmax)
    matches = _greedy_unique_peak_matches(pattern_peaks, refs, tolerance_deg)
    rms = float(np.sqrt(np.mean([m[2]**2 for m in matches]))) if matches else 9.0
    return refs, matches, rms, {"model": "hcp", "a_A": a, "c_A": c}


def _phase_reference_or_fit(dataset: Dataset, db: PhaseDatabase, pid: int, pattern_peaks, tolerance_deg: float):
    model, a_ref, c_ref, fit_enabled = db.phase_lattice_model(pid)
    raw_rows = db.phase_peaks(pid)
    parsed = []
    for dA, hkl, rel in raw_rows:
        hp = _parse_simple_hkl(hkl)
        if hp is not None:
            parsed.append((hp, float(rel) if rel is not None else 100.0))
    xmin, xmax = float(dataset.x.min()), float(dataset.x.max())
    if fit_enabled and model in {"CUBIC_FCC", "CUBIC_BCC"} and parsed:
        fit = _fit_cubic_lattice(pattern_peaks, parsed, dataset.wavelength_A, a_ref, xmin, xmax, tolerance_deg)
        if fit is not None:
            return fit
    if fit_enabled and model == "HCP" and parsed:
        fit = _fit_hcp_lattice(pattern_peaks, parsed, dataset.wavelength_A, a_ref, c_ref, xmin, xmax, tolerance_deg)
        if fit is not None:
            return fit
    refs = []
    for dA, hkl, rel in raw_rows:
        tt = d_to_2theta(dA, dataset.wavelength_A)
        if tt is not None and xmin - 0.5 <= tt <= xmax + 0.5:
            refs.append((float(tt), str(hkl), float(dA), float(rel) if rel is not None else 0.0))
    matches = _greedy_unique_peak_matches(pattern_peaks, refs, tolerance_deg)
    rms = float(np.sqrt(np.mean([m[2]**2 for m in matches]))) if matches else 9.0
    return refs, matches, rms, {"model": "fixed"}


def assign_phases(dataset: Dataset, db: PhaseDatabase, tolerance_deg: float = 0.35,
                  allowed_elements: Optional[set[str]] = None,
                  min_phase_matches: int = 2,
                  min_pattern_rel_intensity: float = 5.0,
                  max_pattern_peaks: int = 20):
    """Whole-pattern phase identification with lattice fitting for simple alloy matrices.

    Built-in FCC/BCC/HCP alloy references store structure + allowed hkl families. Their
    database lattice parameters are only initial search seeds. The actual lattice
    parameter(s) are fitted from the measured pattern before phase scoring. This avoids
    forcing 316L, CuNi30Mn1Fe and ERCuNi peaks onto artificially identical 2theta values.
    """
    if not dataset.peaks:
        dataset.metadata["phase_fits"] = []
        return []

    if allowed_elements is None:
        allowed_elements = dataset_allowed_elements(dataset)
    else:
        allowed_elements = {str(e).strip().capitalize() for e in allowed_elements if str(e).strip()}

    pattern_peaks = [p for p in dataset.peaks if p.rel_intensity >= float(min_pattern_rel_intensity)]
    pattern_peaks = sorted(pattern_peaks, key=lambda p: p.rel_intensity, reverse=True)[:max(2, int(max_pattern_peaks))]
    if not pattern_peaks:
        pattern_peaks = sorted(dataset.peaks, key=lambda p: p.rel_intensity, reverse=True)[:max(2, int(max_pattern_peaks))]

    required_matches = max(2, int(min_phase_matches))
    required_matches = min(required_matches, len(pattern_peaks))

    ranked = []
    phase_cache = {}
    material_text = trusted_material_context(dataset)

    for row in db.all_phases():
        pid, name, formula, structure, sg, category, notes = row
        phase_els = formula_elements(formula)
        if allowed_elements and phase_els and not phase_els.issubset(allowed_elements):
            continue

        refs, matches, rms, fit_info = _phase_reference_or_fit(dataset, db, pid, pattern_peaks, tolerance_deg)
        if not refs:
            continue
        nmatch = len(matches)
        if nmatch < required_matches:
            continue
        expected_code, explicit_matrix = dataset_expected_matrix_code(dataset)
        signature = dataset_material_signature(dataset)
        struct_short = str(structure or "").upper().strip()
        # Do not promote a minor/overlapped BCC pattern to a separate beta/ferrite phase
        # on only two coincident reflections.  Three independent matches are required
        # unless the user explicitly selected that BCC phase as the expected matrix.
        if struct_short == "BCC" and nmatch < 3:
            if signature == "ti64" and expected_code != "beta_ti_bcc":
                continue
            if expected_code == "gamma_fcc" and not explicit_matrix:
                continue

        obs_coverage = nmatch / max(1, len(pattern_peaks))
        intensity_coverage = sum(pattern_peaks[m[0]].rel_intensity for m in matches) / max(
            1e-12, sum(p.rel_intensity for p in pattern_peaks)
        )
        base_score = (
            0.52 * obs_coverage
            + 0.38 * intensity_coverage
            + 0.10 * max(0.0, 1.0 - rms / max(tolerance_deg, 1e-9))
            + material_bias(name, category, material_text)
        )
        chem_factor, chem_capture, chem_reason = phase_chemistry_factor(dataset, phase_els)
        structure_factor, structure_reason = phase_structure_factor(dataset, structure=structure, space_group=sg)
        score = base_score * chem_factor * structure_factor
        if nmatch == 2:
            score -= 0.06
        elif nmatch >= 4:
            score += 0.05
        fit_info = dict(fit_info or {})
        fit_info.update({
            "chem_factor": float(chem_factor), "chem_capture": float(chem_capture),
            "chem_reason": chem_reason, "structure_factor": float(structure_factor),
            "structure_reason": structure_reason,
        })
        ranked.append((score, pid, name, category, nmatch, formula, rms, matches, fit_info))
        phase_cache[pid] = refs

    ranked.sort(reverse=True, key=lambda q: q[0])
    accepted = ranked[:8]

    for p in dataset.peaks:
        candidates = []
        for score, pid, name, category, nmatch, formula, rms, matches, fit_info in accepted:
            refs = phase_cache[pid]
            tt, hkl, dA, rel = min(refs, key=lambda r: abs(r[0] - p.two_theta))
            delta = abs(tt - p.two_theta)
            if delta <= tolerance_deg:
                candidates.append((-score, delta, name, hkl, tt, nmatch, rms))
        if not candidates:
            p.phase = "unassigned"
            p.hkl = ""
            p.delta_2theta = math.nan
            p.phase_source = ""
            continue
        candidates.sort()
        best = candidates[0]
        p.phase = best[2]
        p.hkl = best[3]
        p.delta_2theta = float(best[1])
        p.phase_source = "internal"

    fit_summary = []
    for score, pid, name, category, nmatch, formula, rms, matches, fit_info in accepted:
        item = {
            "phase": name,
            "formula": formula,
            "score": float(score),
            "matched_reflections": int(nmatch),
            "rms_2theta_deg": float(rms),
            "model": fit_info.get("model", "fixed"),
            "chemistry_fit_percent": 100.0 * float(fit_info.get("chem_factor", 1.0)),
            "chemistry_capture_percent": 100.0 * float(fit_info.get("chem_capture", 0.0)),
            "chemistry_note": fit_info.get("chem_reason", ""),
            "structure_fit_percent": 100.0 * float(fit_info.get("structure_factor", 1.0)),
            "structure_note": fit_info.get("structure_reason", ""),
        }
        if "a_A" in fit_info:
            item["a_A"] = float(fit_info["a_A"])
        if "c_A" in fit_info:
            item["c_A"] = float(fit_info["c_A"])
        fit_summary.append(item)
    dataset.metadata["phase_fits"] = fit_summary
    dataset.metadata_sources["phase_fits"] = "calculated from measured peak positions"
    dataset.metadata["phase_identification_source"] = "internal database"
    dataset.metadata_sources["phase_identification_source"] = "analysis setting"
    return accepted


TENTATIVE_SOURCE = "internal-tentative"
TENTATIVE_SUFFIX = "niepotwierdzony, 1 refleks"


def assign_tentative_single_reflections(dataset: Dataset, db: PhaseDatabase, tolerance_deg: float = 0.35,
                                        allowed_elements: Optional[set] = None, min_snr: float = 10.0,
                                        min_rel_intensity: float = 2.0) -> list:
    """Mark a strong, otherwise unexplained peak as a TENTATIVE ferrite / martensite (110).

    The regular assignment demands several reflections of a BCC phase, because a lone peak cannot
    prove a phase. A weak 110 of a ferrite next to a strong austenite pattern is nevertheless a
    common, real finding in stainless steel, so such a peak gets a clearly marked provisional
    label (phase_source "internal-tentative", label suffix "?") instead of staying anonymous.
    Only peaks left unassigned are considered, only steel BCC/BCT phases whose elements are all
    allowed, only peaks with SNR >= min_snr, and only within the angular tolerance of the 110
    reflection. Whether the weaker 200 reflection is visible is recorded as a note, not a condition.
    """
    dataset.metadata.pop("tentative_phases", None)
    dataset.metadata_sources.pop("tentative_phases", None)
    allowed = {str(e).strip().capitalize() for e in (allowed_elements or set()) if str(e).strip()}
    if not allowed:
        return []
    peaks = [p for p in dataset.peaks
             if p.phase in ("", "unassigned") and not p.phase_source
             and math.isfinite(p.snr) and p.snr >= float(min_snr) and p.rel_intensity >= float(min_rel_intensity)]
    if not peaks:
        return []
    material_text = trusted_material_context(dataset)
    has_austenite = any("austen" in str(p.phase).lower() or "γ" in str(p.phase) for p in dataset.peaks)
    options = []
    for pid, name, formula, structure, sg, category, notes in db.all_phases():
        struct = str(structure or "").upper()
        if not ("BCC" in struct or "BCT" in struct) or not str(category or "").lower().startswith("steel"):
            continue
        phase_els = formula_elements(formula)
        if not phase_els or not phase_els.issubset(allowed):
            continue
        d110 = next((float(d) for d, hkl, rel in db.phase_peaks(pid) if str(hkl).strip() == "110"), None)
        tt110 = d_to_2theta(d110, dataset.wavelength_A) if d110 else None
        if tt110 is not None:
            options.append((name, category, float(tt110), d110))
    if not options:
        return []
    sig = None
    try:
        base = getattr(dataset, "baseline", None)
        if base is not None and len(base) == len(dataset.y):
            sig = np.asarray(dataset.y, dtype=float) - np.asarray(base, dtype=float)
    except Exception:
        sig = None
    records = []
    for p in sorted(peaks, key=lambda q: -q.rel_intensity):
        # Positions closer than 0.05 deg count as a tie; then chemistry/material context decides, and in a
        # sample that already contains austenite the delta-ferrite of austenitic steels is preferred.
        ranked = sorted(((int(abs(tt - p.two_theta) / 0.05),
                          -(material_bias(name, category, material_text) + (1.0 if has_austenite and ("δ" in name or "delta" in name.lower()) else 0.0)),
                          abs(tt - p.two_theta), name, tt, d110)
                         for name, category, tt, d110 in options if abs(tt - p.two_theta) <= tolerance_deg))
        if not ranked:
            continue
        _tie, _bias, delta, name, tt, d110 = ranked[0]
        p.phase, p.hkl, p.delta_2theta, p.phase_source = name, "110", float(delta), TENTATIVE_SOURCE
        note = "jeden refleks 110; potwierdź refleksem 200 (ok. {:.1f}°) lub 211 (ok. {:.1f}°)".format(
            2 * math.degrees(math.asin(min(1.0, dataset.wavelength_A / (2 * d110 / math.sqrt(2))))),
            2 * math.degrees(math.asin(min(1.0, dataset.wavelength_A / (2 * d110 / math.sqrt(3))))))
        if sig is not None:
            a_bcc = d110 * math.sqrt(2)
            tt200 = 2 * math.degrees(math.asin(min(1.0, dataset.wavelength_A / (2 * a_bcc / 2.0))))
            window = (dataset.x > tt200 - 0.5) & (dataset.x < tt200 + 0.5)
            if window.any() and len(sig) >= 25:
                # 5-point mean: noise falls by sqrt(5); "visible" needs 4 sigma of that AND half the expected height
                smooth = np.convolve(sig, np.ones(5) / 5.0, mode="same")
                excess = float(np.max(smooth[window]))
                expected = 0.15 * float(p.height)
                noise5 = noise_sigma_from_differences(sig) / math.sqrt(5.0)
                visible = excess >= max(0.5 * expected, 4.0 * noise5)
                note += (f"; w 200 zmierzono {excess:.0f} zliczeń (średnia z 5 punktów) wobec ok. {expected:.0f} "
                         f"oczekiwanych ({'zgodne' if visible else 'brak potwierdzenia'})")
        records.append({"two_theta": float(p.two_theta), "phase": name, "hkl": "110",
                        "delta_2theta": float(delta), "note": note})
    if records:
        dataset.metadata["tentative_phases"] = records
        dataset.metadata_sources["tentative_phases"] = "single-reflection ferrite check"
    return records


def clear_phase_assignments(dataset: Dataset):
    """Wyczyść automatic phase assignments while preserving measured peaks and custom labels."""
    dataset.metadata.pop("tentative_phases", None)
    for p in dataset.peaks:
        p.phase = "unassigned"
        p.hkl = ""
        p.delta_2theta = math.nan
        p.phase_source = ""
    dataset.metadata["phase_fits"] = []


def assign_phases_from_cod(dataset: Dataset, cod_db: CODDatabase, tolerance_deg: float = 0.35,
                           allowed_elements: Optional[set[str]] = None,
                           min_phase_matches: int = 2,
                           min_pattern_rel_intensity: float = 5.0,
                           max_pattern_peaks: int = 20,
                           max_candidates: int = 5000,
                           top_results: int = 8):
    """Assign phases directly from the external COD .db3 database.

    This path deliberately does NOT consult the internal SQLite phase library. COD unit-cell
    and space-group records are screened against the measured whole pattern; the top
    candidates are then used to label measured peaks from generated allowed reflections.
    COD powder intensities are generally unavailable, so this remains a positional
    candidate-screening/identification workflow rather than a Rietveld refinement.
    """
    clear_phase_assignments(dataset)
    if not dataset.peaks:
        dataset.metadata["phase_identification_source"] = "external COD only"
        dataset.metadata_sources["phase_identification_source"] = "analysis setting"
        return []
    if not cod_db.is_open:
        raise RuntimeError("External COD-only identification is selected, but no COD .db3 database is loaded.")
    if not ensure_gemmi():
        raise RuntimeError(gemmi_diagnostic_message("External COD phase identification"))

    if allowed_elements is None:
        allowed_elements = dataset_allowed_elements(dataset)
    else:
        allowed_elements = {str(e).strip().capitalize() for e in allowed_elements if str(e).strip()}
    if not allowed_elements:
        raise ValueError(
            "Identyfikacja wyłącznie w zewnętrznej bazie COD wymaga filtra dozwolonych pierwiastków. "
            "Najpierw wybierz skład próbki przyciskiem „Układ okresowy…”."
        )

    # Whole-pattern search first.  Limiting the observed set to the strongest reflections
    # suppresses noise and makes the large COD scan tractable.
    results, used_obs, candidate_count = cod_db.search_measured_pattern(
        dataset, allowed_elements, tolerance_deg=tolerance_deg,
        max_candidates=max(100, int(max_candidates)),
        top_results=max(3, int(top_results)),
        max_observed_peaks=max(2, min(int(max_pattern_peaks), 12)),
        min_relative_intensity=float(min_pattern_rel_intensity),
    )

    required = max(2, min(int(min_phase_matches), len(used_obs) if used_obs else int(min_phase_matches)))
    accepted = [r for r in results if int(r.get("matched", 0)) >= required]
    expected_code, explicit_matrix = dataset_expected_matrix_code(dataset)
    signature = dataset_material_signature(dataset)
    filtered_secondary = []
    for r in accepted:
        st = cod_structure_short(r.get("space_group", ""), r.get("sg_no"))
        nm = int(r.get("matched", 0))
        if st == "BCC" and nm < 3:
            if signature == "ti64" and expected_code != "beta_ti_bcc":
                continue
            if expected_code == "gamma_fcc" and not explicit_matrix:
                continue
        filtered_secondary.append(r)
    accepted = filtered_secondary
    # When the wizard gives strong evidence for a dominant matrix element, do not let an
    # extremely chemistry-incompatible whole-pattern candidate automatically label all strong
    # peaks. Such records remain visible in the COD search table for expert inspection.
    _dom, _dom_w, _strong_dom = dataset_dominant_element(dataset)
    if _strong_dom:
        chemically_plausible = [r for r in accepted if float(r.get("chem_factor", 1.0)) >= 0.18]
        if chemically_plausible:
            accepted = chemically_plausible
    accepted = accepted[:max(1, int(top_results))]

    ref_cache = {}
    xmin, xmax = float(dataset.x.min()), float(dataset.x.max())
    for r in accepted:
        try:
            ref = cod_db.entry_to_reference(int(r["cod_id"]))
        except Exception:
            continue
        refs = []
        for dA, hkl, rel in ref.get("peaks", []):
            tt = d_to_2theta(float(dA), dataset.wavelength_A)
            if tt is not None and xmin - 0.5 <= tt <= xmax + 0.5:
                refs.append((float(tt), str(hkl), float(dA), float(rel)))
        if refs:
            ref_cache[int(r["cod_id"])] = refs

    # Use whole-pattern score as the primary discriminator and peak distance second.
    for p in dataset.peaks:
        candidates = []
        for r in accepted:
            cid = int(r["cod_id"])
            refs = ref_cache.get(cid)
            if not refs:
                continue
            tt, hkl, dA, rel = min(refs, key=lambda q: abs(q[0] - p.two_theta))
            delta = abs(float(tt) - float(p.two_theta))
            if delta <= float(tolerance_deg):
                candidates.append((-float(r["score"]), delta, r, hkl, tt))
        if not candidates:
            continue
        candidates.sort(key=lambda q: (q[0], q[1]))
        _negscore, delta, r, hkl, tt = candidates[0]
        display = str(r.get("plot_name") or r.get("name") or r.get("formula") or "COD structure")
        p.phase = display
        p.hkl = canonical_hkl_for_structure(str(hkl), cod_structure_short(r.get("space_group", ""), r.get("sg_no")))
        p.delta_2theta = float(delta)
        p.phase_source = "cod"

    fit_summary = []
    for r in accepted:
        fit_summary.append({
            "phase": r.get("name") or r.get("formula") or "COD structure",
            "plot_phase": r.get("plot_name") or r.get("name") or r.get("formula") or "COD structure",
            "formula": r.get("formula", ""),
            "score": float(r.get("score", 0.0)) / 100.0,
            "score_percent": float(r.get("score", 0.0)),
            "matched_reflections": int(r.get("matched", 0)),
            "rms_2theta_deg": float(r.get("rms", math.nan)),
            "model": "COD external unit cell / space group",
            "chemistry_fit_percent": 100.0 * float(r.get("chem_factor", 1.0)),
            "chemistry_capture_percent": 100.0 * float(r.get("chem_capture", 0.0)),
            "chemistry_note": r.get("chem_reason", ""),
            "structure_fit_percent": 100.0 * float(r.get("structure_factor", 1.0)),
            "structure_note": r.get("structure_reason", ""),
            "cod_id": int(r["cod_id"]),
            "space_group": r.get("space_group", ""),
            "a_A": r.get("a", ""),
            "b_A": r.get("b", ""),
            "c_A": r.get("c", ""),
        })
    dataset.metadata["phase_fits"] = fit_summary
    dataset.metadata["cod_auto_candidate_count"] = int(candidate_count)
    dataset.metadata["phase_identification_source"] = "external COD only"
    dataset.metadata_sources["phase_fits"] = "external COD .db3 whole-pattern search"
    dataset.metadata_sources["phase_identification_source"] = "analysis setting"
    return accepted


def apply_manual_cod_peak_overrides(dataset: Dataset) -> int:
    """Reapply user-approved COD peak assignments after peak detection/automatic matching.

    Assignments are stored by measured 2θ and intentionally override automatic labels.
    They remain valid across project save/reopen and ordinary re-analysis as long as a
    corresponding detected peak is still present within the saved tolerance.
    """
    rows = dataset.metadata.get("manual_cod_peak_assignments", []) if isinstance(dataset.metadata, dict) else []
    if not isinstance(rows, list) or not rows or not dataset.peaks:
        return 0
    applied = 0
    for rec in rows:
        if not isinstance(rec, dict):
            continue
        try:
            tt = float(rec.get("peak_two_theta"))
        except Exception:
            continue
        p = min(dataset.peaks, key=lambda q: abs(q.two_theta - tt))
        tol = max(0.08, float(rec.get("tolerance_deg", 0.35) or 0.35))
        if abs(p.two_theta - tt) > tol:
            continue
        p.phase = str(rec.get("phase", "COD structure") or "COD structure")
        p.hkl = str(rec.get("hkl", "") or "")
        try:
            p.delta_2theta = float(rec.get("delta_2theta", math.nan))
        except Exception:
            p.delta_2theta = math.nan
        p.phase_source = "cod"
        applied += 1
    return applied


def combine_internal_and_cod_assignments(dataset: Dataset, db: PhaseDatabase, cod_db: CODDatabase,
                                         tolerance_deg: float, allowed_elements: set[str],
                                         min_phase_matches: int, min_pattern_rel_intensity: float,
                                         max_pattern_peaks: int):
    """Use both sources; whole-pattern scores decide between conflicting assignments.

    Internal simple-alloy references are lattice-fitted to the measured pattern. COD is then
    searched independently.  A COD label replaces an internal one only when its normalized
    whole-pattern score is clearly stronger; otherwise the internal fitted-alloy label stays.
    """
    assign_phases(dataset, db, tolerance_deg, allowed_elements, min_phase_matches,
                  min_pattern_rel_intensity, max_pattern_peaks)
    internal_state = [(p.phase, p.hkl, p.delta_2theta, p.phase_source) for p in dataset.peaks]
    internal_fits = list(dataset.metadata.get("phase_fits", []) or [])
    internal_scores = {str(f.get("phase", "")): float(f.get("score", 0.0))*100.0 for f in internal_fits if isinstance(f, dict)}

    try:
        assign_phases_from_cod(dataset, cod_db, tolerance_deg, allowed_elements, min_phase_matches,
                               min_pattern_rel_intensity, max_pattern_peaks)
        cod_state = [(p.phase, p.hkl, p.delta_2theta, p.phase_source) for p in dataset.peaks]
        cod_fits = list(dataset.metadata.get("phase_fits", []) or [])
        cod_scores = {str(f.get("phase", "")): float(f.get("score_percent", 0.0)) for f in cod_fits if isinstance(f, dict)}
    except Exception:
        # If COD is unavailable, preserve the internal analysis in combined mode.
        cod_state = [("unassigned", "", math.nan, "") for _ in dataset.peaks]
        cod_fits = []
        cod_scores = {}

    for p, i_state, c_state in zip(dataset.peaks, internal_state, cod_state):
        iph, ihkl, idelta, isource = i_state
        cph, chkl, cdelta, csource = c_state
        if cph != "unassigned":
            cscore = cod_scores.get(cph, 0.0)
            iscore = internal_scores.get(iph, 0.0)
            # External COD must be materially stronger to override a fitted internal alloy model.
            if iph == "unassigned" or cscore >= iscore + 7.5:
                p.phase, p.hkl, p.delta_2theta, p.phase_source = cph, chkl, cdelta, (csource or "cod")
            else:
                p.phase, p.hkl, p.delta_2theta, p.phase_source = iph, ihkl, idelta, (isource or "internal")
        else:
            p.phase, p.hkl, p.delta_2theta, p.phase_source = iph, ihkl, idelta, (isource or ("internal" if iph != "unassigned" else ""))

    dataset.metadata["phase_fits"] = internal_fits + cod_fits
    dataset.metadata["phase_identification_source"] = "internal + external COD"
    dataset.metadata_sources["phase_fits"] = "combined internal lattice-fit + external COD search"
    dataset.metadata_sources["phase_identification_source"] = "analysis setting"
    return dataset.metadata["phase_fits"]


# --------------------------- robust structure-family fallback ---------------------------

GENERIC_FCC_HKLS = [
    ((1,1,1),100.0), ((2,0,0),85.0), ((2,2,0),65.0), ((3,1,1),55.0),
    ((2,2,2),35.0), ((4,0,0),28.0), ((3,3,1),24.0), ((4,2,0),22.0),
]
GENERIC_BCC_HKLS = [
    ((1,1,0),100.0), ((2,0,0),55.0), ((2,1,1),70.0), ((2,2,0),35.0),
    ((3,1,0),28.0), ((2,2,2),24.0), ((3,2,1),20.0), ((4,0,0),18.0),
]


def _generic_cubic_fit(dataset: Dataset, structure: str, tolerance_deg: float = 0.35,
                       min_rel_intensity: float = 2.0, max_pattern_peaks: int = 16):
    """Fit a generic FCC/BCC reflection family without assuming a material grade.

    This is intentionally database-independent.  It is used as a safe fallback when COD
    cannot be searched (e.g. chemistry is unknown) or when an exact COD compound is less
    informative than the measured lattice family.  It identifies *structure*, not chemistry.
    """
    if not dataset.peaks:
        return None
    st = str(structure).upper()
    hkl_rows = GENERIC_FCC_HKLS if st == 'FCC' else GENERIC_BCC_HKLS
    obs = [q for q in dataset.peaks if q.rel_intensity >= float(min_rel_intensity)]
    obs = sorted(obs, key=lambda q: q.rel_intensity, reverse=True)[:max(3, int(max_pattern_peaks))]
    if len(obs) < 2:
        return None
    # Broad but physically realistic range for common metallic cubic matrices.
    lo_a, hi_a = ((2.8, 5.6) if st == 'FCC' else (2.35, 4.5))
    xmin, xmax = float(np.nanmin(dataset.x)), float(np.nanmax(dataset.x))
    candidates = []
    for q in obs[:12]:
        for hkl, _rel in hkl_rows[:6]:
            n = sum(v*v for v in hkl)
            a = float(q.d_A) * math.sqrt(n)
            if lo_a <= a <= hi_a:
                candidates.append(round(a, 6))
    if not candidates:
        return None
    best = None
    total_int = max(1e-12, sum(q.rel_intensity for q in obs))
    for a in sorted(set(candidates)):
        refs = _refs_for_cubic_hkls(hkl_rows, a, dataset.wavelength_A, xmin, xmax)
        matches = _greedy_unique_peak_matches(obs, refs, float(tolerance_deg))
        if len(matches) < 2:
            continue
        icov = sum(obs[m[0]].rel_intensity for m in matches) / total_int
        rms = float(np.sqrt(np.mean([m[2]**2 for m in matches])))
        coverage = len(matches) / max(1, min(len(obs), len(refs)))
        # Strong patterns with 3+ reflections are preferred; two-line coincidences remain weak.
        score = 0.48*icov + 0.34*coverage + 0.18*max(0.0, 1.0-rms/max(float(tolerance_deg),1e-9))
        if len(matches) == 2:
            score *= 0.72
        elif len(matches) >= 4:
            score *= 1.06
        key = (score, len(matches), icov, -rms)
        if best is None or key > best[0]:
            best = (key, a, refs, matches, rms, icov)
    if best is None:
        return None
    a = float(best[1])
    # Refine a from all matched d-spacings, weighted by measured intensity.
    for _ in range(2):
        refs = _refs_for_cubic_hkls(hkl_rows, a, dataset.wavelength_A, xmin, xmax)
        matches = _greedy_unique_peak_matches(obs, refs, float(tolerance_deg))
        vals=[]; w=[]
        for oi, ri, _dlt in matches:
            hp = _parse_simple_hkl(refs[ri][1])
            if hp is None:
                continue
            vals.append(float(obs[oi].d_A)*math.sqrt(sum(v*v for v in hp)))
            w.append(max(1.0, float(obs[oi].rel_intensity)))
        if vals:
            a = float(np.average(vals, weights=w))
    refs = _refs_for_cubic_hkls(hkl_rows, a, dataset.wavelength_A, xmin, xmax)
    matches = _greedy_unique_peak_matches(obs, refs, float(tolerance_deg))
    if len(matches) < 2:
        return None
    rms = float(np.sqrt(np.mean([m[2]**2 for m in matches])))
    icov = sum(obs[m[0]].rel_intensity for m in matches) / total_int
    coverage = len(matches) / max(1, min(len(obs), len(refs)))
    score = 0.48*icov + 0.34*coverage + 0.18*max(0.0, 1.0-rms/max(float(tolerance_deg),1e-9))
    if len(matches)==2: score*=0.72
    elif len(matches)>=4: score*=1.06
    return {
        'structure': st, 'a_A': a, 'score': float(score), 'rms': rms,
        'matched': len(matches), 'intensity_coverage': float(icov),
        'refs': refs, 'matches': matches, 'observed': obs,
    }


def _generic_hcp_ti_fit(dataset: Dataset, tolerance_deg: float = 0.35,
                        min_rel_intensity: float = 2.0, max_pattern_peaks: int = 18):
    """Conservative HCP fallback for user-confirmed Ti-containing samples."""
    allowed = dataset_allowed_elements(dataset)
    if 'Ti' not in allowed or not dataset.peaks:
        return None
    obs = [q for q in dataset.peaks if q.rel_intensity >= float(min_rel_intensity)]
    obs = sorted(obs, key=lambda q:q.rel_intensity, reverse=True)[:max(4,int(max_pattern_peaks))]
    if len(obs)<3:
        return None
    hcp = [((1,0,0),80),((0,0,2),70),((1,0,1),100),((1,0,2),45),((1,1,0),40),((1,0,3),35),((1,1,2),30),((2,0,1),25),((0,0,4),20)]
    # Ti / Ti-alloy seeds; fitting routine refines both a and c.
    fit = _fit_hcp_lattice(obs, hcp, dataset.wavelength_A, 2.95, 4.68,
                           float(np.nanmin(dataset.x)), float(np.nanmax(dataset.x)), float(tolerance_deg))
    if fit is None:
        return None
    refs,matches,rms,info=fit
    if len(matches)<3:
        return None
    total=max(1e-12,sum(q.rel_intensity for q in obs))
    icov=sum(obs[m[0]].rel_intensity for m in matches)/total
    coverage=len(matches)/max(1,min(len(obs),len(refs)))
    score=0.50*icov+0.32*coverage+0.18*max(0.0,1.0-rms/max(float(tolerance_deg),1e-9))
    return {'structure':'HCP','a_A':info.get('a_A'),'c_A':info.get('c_A'),'score':float(score),
            'rms':float(rms),'matched':len(matches),'intensity_coverage':float(icov),
            'refs':refs,'matches':matches,'observed':obs}


def _structure_context_label(dataset: Dataset, structure: str) -> str:
    """Human-facing phase-family name based only on user-confirmed chemistry/structure hints."""
    st=str(structure).upper()
    code,_explicit=dataset_expected_matrix_code(dataset)
    els=dataset_allowed_elements(dataset)
    has=lambda *e: set(e).issubset(els)
    if st=='FCC':
        if code=='gamma_fcc': return 'γ-austenite FCC'
        if code=='cuni_fcc': return 'Cu–Ni FCC'
        if code=='al_fcc': return 'α-Al FCC'
        if has('Fe','Cr','Ni') and 'Cu' not in els: return 'γ-austenite FCC'
        if has('Cu','Ni') and 'Cr' not in els and 'Mo' not in els: return 'Cu–Ni FCC'
        if has('Fe','Cr','Ni','Cu'): return 'FCC solid solution (Fe–Cr–Ni–Cu)'
        return 'FCC structure (candidate)'
    if st=='BCC':
        if code=='beta_ti_bcc': return 'β-Ti BCC'
        if code=='ferrite_bcc': return 'ferrite / δ-ferrite BCC'
        if 'Ti' in els: return 'β-Ti BCC (candidate)'
        if 'Fe' in els: return 'ferrite / BCC (candidate)'
        return 'BCC structure (candidate)'
    if st=='HCP':
        if code=='alpha_ti_hcp' or 'Ti' in els: return 'α-Ti HCP'
        return 'HCP structure (candidate)'
    return f'{st} structure (candidate)'


def apply_structural_family_fallback(dataset: Dataset, tolerance_deg: float = 0.35,
                                     min_rel_intensity: float = 2.0,
                                     max_pattern_peaks: int = 18,
                                     overwrite_weak_cod: bool = True):
    """Fill unassigned peaks using robust lattice-family recognition.

    The routine does not use filenames and does not claim an exact compound.  It fits FCC/BCC
    (and Ti-HCP when Ti was explicitly allowed), assigns conventional hkl families and stores
    an auditable structure-only result.  This mirrors the conservative workflow in search-match
    software: use composition restraints when known, but keep structural evidence distinct from
    database chemistry.
    """
    fits=[]
    for st in ('FCC','BCC'):
        f=_generic_cubic_fit(dataset,st,tolerance_deg,min_rel_intensity,max_pattern_peaks)
        if f is not None: fits.append(f)
    h=_generic_hcp_ti_fit(dataset,tolerance_deg,min_rel_intensity,max_pattern_peaks)
    if h is not None: fits.append(h)
    if not fits:
        dataset.metadata['structural_family_fits']=[]
        return []
    fits.sort(key=lambda z:(z['score'],z['matched'],z['intensity_coverage'],-z['rms']), reverse=True)
    # Avoid a 2-line alternative overriding a well-supported family.
    accepted=[]
    top=fits[0]
    if top['matched']>=3 and top['score']>=0.48:
        accepted=[top]
    elif top['matched']>=2 and top['score']>=0.58:
        accepted=[top]
    # A secondary structure must explain at least two measured peaks that are NOT already
    # explained by the dominant family. This suppresses the common false BCC solution where
    # FCC (200)/(220)/(222) are merely re-indexed as BCC (110)/(200)/(211).
    if accepted:
        top_tts={round(float(top['observed'][m[0]].two_theta),3) for m in top['matches']}
        for cand in fits[1:]:
            unique=[]
            for m in cand['matches']:
                tt=round(float(cand['observed'][m[0]].two_theta),3)
                if all(abs(tt-x)>max(0.12,0.45*float(tolerance_deg)) for x in top_tts):
                    unique.append(tt)
            if cand['matched']>=3 and cand['score']>=0.50 and len(unique)>=2:
                cand['unique_secondary_peaks']=unique
                accepted.append(cand)
                break
    cod_fit_by_name={}
    for r in (dataset.metadata.get('phase_fits',[]) or []):
        if isinstance(r,dict):
            for k in ('plot_phase','phase','name'):
                if r.get(k): cod_fit_by_name[str(r.get(k))]=r
    for fam in accepted:
        label=_structure_context_label(dataset,fam['structure'])
        fam['label']=label
        for oi,ri,delta in fam['matches']:
            p=fam['observed'][oi]
            ref=fam['refs'][ri]
            # Determine whether an existing COD label is weak enough to replace.
            replace=(str(p.phase or '').strip().lower() in ('','unassigned'))
            if not replace and overwrite_weak_cod and p.phase_source=='cod':
                rec=cod_fit_by_name.get(str(p.phase),{})
                chem=float(rec.get('chemistry_fit_percent',100.0) or 100.0) if isinstance(rec,dict) else 100.0
                struct=float(rec.get('structure_fit_percent',100.0) or 100.0) if isinstance(rec,dict) else 100.0
                replace=(chem<45.0 or struct<65.0)
            if not replace:
                continue
            p.phase=label
            p.hkl=canonical_hkl_for_structure(str(ref[1]),fam['structure'])
            p.delta_2theta=float(delta)
            p.phase_source='structure'
    summary=[]
    for fam in accepted:
        item={'phase':fam.get('label',fam['structure']), 'plot_phase':fam.get('label',fam['structure']),
              'score':fam['score'], 'score_percent':100*fam['score'], 'matched_reflections':fam['matched'],
              'rms_2theta_deg':fam['rms'], 'model':f"structure-only {fam['structure']} lattice fit",
              'a_A':fam.get('a_A',''), 'c_A':fam.get('c_A',''), 'source':'structure fallback'}
        summary.append(item)
    dataset.metadata['structural_family_fits']=summary
    dataset.metadata_sources['structural_family_fits']='measured peak lattice-family fit'
    return summary



# --------------------------- AI-assisted evidence ranking ---------------------------
AI_EVIDENCE_MODEL_NAME = "Hybrid XRD evidence ranker v1"


def _clamp01(value, default: float = 0.0) -> float:
    try:
        v = float(value)
        if not math.isfinite(v):
            return float(default)
        return max(0.0, min(1.0, v))
    except Exception:
        return float(default)


def _phase_candidate_source(rec: dict) -> str:
    """Return a compact provenance label for one phase-fit record."""
    if rec.get("cod_id") not in (None, ""):
        return "COD"
    src = str(rec.get("source", "") or "").lower()
    model = str(rec.get("model", "") or "").lower()
    if "structure" in src or "structure-only" in model:
        return "Structure fit"
    return "Internal"


def _candidate_publication_label(dataset: Dataset, phase_name: str) -> str:
    """Convert a database/technical phase name into a cautious publication label."""
    _short, full = _publication_phase_names(str(phase_name or ""), dataset.material_type)
    return full or str(phase_name or "candidate")


def build_ai_assisted_ranking(dataset: Dataset, tolerance_deg: float = 0.35) -> List[dict]:
    """Fuse XRD, chemistry and database evidence into an advisory candidate ranking.

    This is intentionally *not* an autonomous phase classifier and it never changes peak
    assignments.  It combines the evidence already produced by the physics/database
    workflow (COD, internal references and lattice-family fitting) into one auditable score.
    The score is an evidence index, not a statistical probability.  User confirmation is
    required before a suggested publication label is adopted.
    """
    md = dataset.metadata if isinstance(dataset.metadata, dict) else {}
    raw_records: List[dict] = []
    for key in ("phase_fits", "structural_family_fits"):
        rows = md.get(key, [])
        if isinstance(rows, list):
            for rec in rows:
                if isinstance(rec, dict):
                    raw_records.append(dict(rec))

    if not raw_records:
        md["ai_phase_ranking"] = []
        md["ai_phase_model"] = AI_EVIDENCE_MODEL_NAME
        dataset.metadata_sources["ai_phase_ranking"] = "derived from current XRD/database evidence"
        return []

    strong = [p for p in dataset.peaks if float(getattr(p, "rel_intensity", 0.0)) >= 5.0]
    n_strong = max(1, len(strong))
    chemistry_known = bool(dataset_allowed_elements(dataset))
    tol = max(0.05, float(tolerance_deg or 0.35))
    ranked: List[dict] = []
    seen = set()

    for rec in raw_records:
        source = _phase_candidate_source(rec)
        phase = str(rec.get("plot_phase") or rec.get("phase") or rec.get("name") or rec.get("formula") or "candidate")
        formula = str(rec.get("formula", "") or "")
        dedupe = (source, phase.strip().lower(), formula.strip().lower())
        if dedupe in seen:
            continue
        seen.add(dedupe)

        raw_score = rec.get("score_percent", None)
        if raw_score in (None, ""):
            try:
                raw_score = 100.0 * float(rec.get("score", 0.0))
            except Exception:
                raw_score = 0.0
        raw_norm = _clamp01(float(raw_score) / 100.0)

        try:
            matched = int(rec.get("matched_reflections", rec.get("matched", 0)) or 0)
        except Exception:
            matched = 0
        match_target = max(3, min(5, n_strong))
        match_norm = _clamp01(matched / max(1.0, float(match_target)))

        try:
            rms = float(rec.get("rms_2theta_deg", rec.get("rms", math.nan)))
        except Exception:
            rms = math.nan
        rms_norm = math.exp(-0.5 * (rms / tol) ** 2) if math.isfinite(rms) else 0.45

        if chemistry_known:
            chem_pct = rec.get("chemistry_fit_percent", None)
            chem_norm = _clamp01(float(chem_pct) / 100.0, 0.55) if chem_pct not in (None, "") else 0.55
        else:
            chem_pct = None
            # Unknown chemistry must remain neutral, not accidentally become a positive signal.
            chem_norm = 0.68

        struct_pct = rec.get("structure_fit_percent", None)
        struct_norm = _clamp01(float(struct_pct) / 100.0, 0.72) if struct_pct not in (None, "") else 0.72
        source_prior = {"COD": 0.90, "Internal": 0.94, "Structure fit": 0.82}.get(source, 0.85)

        # Transparent evidence fusion.  Whole-pattern/matching evidence carries most weight;
        # chemistry and matrix expectations refine the ordering but cannot create a phase on
        # their own.  Two-line coincidences are deliberately down-weighted.
        evidence = (
            0.27 * raw_norm
            + 0.24 * match_norm
            + 0.18 * rms_norm
            + 0.16 * chem_norm
            + 0.09 * struct_norm
            + 0.06 * source_prior
        )
        penalty_reasons = []
        if matched < 2:
            evidence *= 0.45; penalty_reasons.append("fewer than 2 matched reflections")
        elif matched == 2:
            evidence *= 0.78; penalty_reasons.append("only 2 matched reflections")
        if chemistry_known and chem_norm < 0.35:
            evidence *= 0.66; penalty_reasons.append("weak chemistry compatibility")
        if math.isfinite(rms) and rms > tol:
            evidence *= 0.75; penalty_reasons.append("RMS mismatch exceeds tolerance")
        evidence = _clamp01(evidence)

        if evidence >= 0.82 and matched >= 4:
            band = "high"
        elif evidence >= 0.68 and matched >= 3:
            band = "moderate-high"
        elif evidence >= 0.52:
            band = "moderate"
        else:
            band = "low"

        pub = _candidate_publication_label(dataset, phase)
        reasons = [
            f"{matched} matched reflection(s)",
            (f"RMS Δ2θ={rms:.3f}°" if math.isfinite(rms) else "RMS Δ2θ unavailable"),
            (f"chemistry={100*chem_norm:.0f}%" if chemistry_known else "chemistry unknown/neutral"),
            f"structure={100*struct_norm:.0f}%",
            f"source={source}",
        ]
        reasons.extend(penalty_reasons)
        ranked.append({
            "phase": phase,
            "plot_phase": str(rec.get("plot_phase") or phase),
            "publication_label": pub,
            "formula": formula,
            "source": source,
            "cod_id": rec.get("cod_id", ""),
            "matched_reflections": matched,
            "rms_2theta_deg": rms,
            "base_score_percent": 100.0 * raw_norm,
            "chemistry_fit_percent": (100.0 * chem_norm if chemistry_known else None),
            "structure_fit_percent": 100.0 * struct_norm,
            "ai_evidence_score_percent": 100.0 * evidence,
            "confidence_band": band,
            "reason": "; ".join(reasons),
            "model": AI_EVIDENCE_MODEL_NAME,
        })

    ranked.sort(key=lambda r: (
        float(r.get("ai_evidence_score_percent", 0.0)),
        int(r.get("matched_reflections", 0)),
        -float(r.get("rms_2theta_deg", 99.0)) if math.isfinite(float(r.get("rms_2theta_deg", math.nan))) else -99.0,
    ), reverse=True)

    # Relative support helps the user see whether the first two candidates are essentially tied.
    if ranked:
        ev = np.asarray([float(r["ai_evidence_score_percent"]) / 100.0 for r in ranked], dtype=float)
        z = np.exp(5.0 * (ev - float(np.max(ev))))
        z = z / max(float(np.sum(z)), 1e-12)
        for i, (r, support) in enumerate(zip(ranked, z), start=1):
            r["rank"] = i
            r["relative_support_percent"] = 100.0 * float(support)
        if len(ranked) > 1:
            margin = float(ranked[0]["ai_evidence_score_percent"]) - float(ranked[1]["ai_evidence_score_percent"])
        else:
            margin = float(ranked[0]["ai_evidence_score_percent"])
        ranked[0]["top_margin_percent"] = margin

    md["ai_phase_ranking"] = ranked
    md["ai_phase_model"] = AI_EVIDENCE_MODEL_NAME
    dataset.metadata_sources["ai_phase_ranking"] = "advisory fusion of XRD/COD/internal/chemistry evidence"
    return ranked


def publication_description_from_ranking(dataset: Dataset, ranking: Optional[List[dict]] = None) -> str:
    """Generate cautious publication-style English wording from the current evidence."""
    if ranking is None:
        ranking = dataset.metadata.get("ai_phase_ranking", []) if isinstance(dataset.metadata, dict) else []
    if not isinstance(ranking, list) or not ranking:
        return ("No phase assignment is proposed because the current analysis does not contain "
                "sufficient database/structural evidence. Peak positions can still be reported independently.")

    confirmed = dataset.metadata.get("ai_user_confirmed_candidate", {}) if isinstance(dataset.metadata, dict) else {}
    top = None
    if isinstance(confirmed, dict) and confirmed.get("phase"):
        cphase = str(confirmed.get("phase")).strip().lower()
        for r in ranking:
            if str(r.get("phase", "")).strip().lower() == cphase:
                top = r; break
    if top is None:
        top = ranking[0]

    label = str(top.get("publication_label") or top.get("phase") or "the leading phase candidate")
    matched = int(top.get("matched_reflections", 0) or 0)
    rms = top.get("rms_2theta_deg", math.nan)
    rms_txt = f" with an RMS 2θ mismatch of {float(rms):.3f}°" if isinstance(rms, (int, float)) and math.isfinite(float(rms)) else ""
    confirmed_text = bool(isinstance(confirmed, dict) and confirmed.get("phase"))
    if confirmed_text:
        lead = f"The diffraction pattern was assigned to {label}"
        tail = "The assignment was user-confirmed after database/structure review."
    else:
        lead = f"The diffraction pattern is most consistent with {label}"
        tail = ("This is an advisory candidate ranking and should be confirmed against the known chemistry "
                "and, where necessary, whole-pattern/Rietveld refinement before being reported as definitive.")
    return f"{lead}, based on {matched} matched reflection(s){rms_txt}. {tail}"


class AIAssistedReviewDialog(QDialog):
    """Advisory comparison of COD, internal and structure-only phase candidates."""
    def __init__(self, dataset: Dataset, tolerance_deg: float = 0.35, parent=None):
        super().__init__(parent)
        self.dataset = dataset
        self.tolerance_deg = float(tolerance_deg)
        self.setWindowTitle("Wspomagana AI ocena faz – ranking dowodów")
        self.resize(1320, 780)
        lay = QVBoxLayout(self)

        info = QLabel(
            "Lokalny model wspomagający: dowody z pików i całego profilu XRD + wyniki COD/bazy wewnętrznej + chemia zatwierdzona przez użytkownika. "
            "Wynik jest indeksem siły dowodów, a nie prawdopodobieństwem, i nigdy nie zatwierdza fazy automatycznie. "
            "Wybierz kandydata i zatwierdź go dopiero po ocenie naukowej."
        )
        info.setWordWrap(True); lay.addWidget(info)

        self.table = QTableWidget(0, 12)
        self.table.setHorizontalHeaderLabels([
            "Pozycja", "Indeks dowodów [0–100]", "Poziom", "Wsparcie względne [%]", "Źródło", "Kandydat",
            "Etykieta publikacyjna", "Wzór", "Dopasowane", "RMS Δ2θ [°]", "Chemia [%]", "Struktura [%]"
        ])
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(5, QHeaderView.ResizeMode.Stretch)
        self.table.horizontalHeader().setSectionResizeMode(6, QHeaderView.ResizeMode.Stretch)
        self.table.itemSelectionChanged.connect(self._selection_changed)
        lay.addWidget(self.table, 1)

        self.detail = QTextEdit(); self.detail.setReadOnly(True); self.detail.setMinimumHeight(130)
        lay.addWidget(self.detail)
        self.description = QTextEdit(); self.description.setReadOnly(False); self.description.setMinimumHeight(90)
        self.description.setToolTip("Sformułowanie w stylu publikacyjnym wygenerowane na podstawie aktualnie ocenianych dowodów. Można je edytować i kopiować; program nie wstawia go automatycznie do manuskryptu.")
        lay.addWidget(QLabel("Opis w stylu publikacyjnym (po angielsku):")); lay.addWidget(self.description)

        bl = QHBoxLayout()
        self.rebuild_btn = QPushButton("Przelicz ranking")
        self.confirm_btn = QPushButton("Zatwierdź wybranego kandydata + etykietę publikacyjną")
        self.clear_btn = QPushButton("Usuń zatwierdzenie użytkownika / aliasy")
        self.copy_btn = QPushButton("Kopiuj opis")
        close = QPushButton("Zamknij")
        for b in (self.rebuild_btn, self.confirm_btn, self.clear_btn, self.copy_btn): bl.addWidget(b)
        bl.addStretch(1); bl.addWidget(close); lay.addLayout(bl)
        self.rebuild_btn.clicked.connect(self.rebuild)
        self.confirm_btn.clicked.connect(self.confirm_selected)
        self.clear_btn.clicked.connect(self.clear_confirmation)
        self.copy_btn.clicked.connect(self.copy_description)
        close.clicked.connect(self.accept)
        self.ranking: List[dict] = []
        self.rebuild()

    def rebuild(self):
        self.ranking = build_ai_assisted_ranking(self.dataset, self.tolerance_deg)
        self.table.setRowCount(len(self.ranking))
        confirmed = self.dataset.metadata.get("ai_user_confirmed_candidate", {}) if isinstance(self.dataset.metadata, dict) else {}
        cphase = str(confirmed.get("phase", "") if isinstance(confirmed, dict) else "").strip().lower()
        for row, r in enumerate(self.ranking):
            vals = [
                r.get("rank", row+1), f"{float(r.get('ai_evidence_score_percent',0)):.1f}", r.get("confidence_band", ""),
                f"{float(r.get('relative_support_percent',0)):.1f}", r.get("source", ""), r.get("phase", ""),
                r.get("publication_label", ""), r.get("formula", ""), r.get("matched_reflections", 0),
                ("" if not math.isfinite(float(r.get("rms_2theta_deg", math.nan))) else f"{float(r.get('rms_2theta_deg')):.3f}"),
                ("nieznana" if r.get("chemistry_fit_percent") is None else f"{float(r.get('chemistry_fit_percent')):.0f}"),
                f"{float(r.get('structure_fit_percent',0)):.0f}",
            ]
            for col, val in enumerate(vals):
                item = QTableWidgetItem(str(val)); item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(row, col, item)
            if cphase and str(r.get("phase", "")).strip().lower() == cphase:
                for col in range(self.table.columnCount()):
                    it = self.table.item(row, col)
                    if it: it.setToolTip("Kandydat zatwierdzony przez użytkownika")
        if self.ranking:
            self.table.selectRow(0)
        else:
            self.detail.setPlainText("Brak dostępnych kandydatów. Najpierw przeanalizuj próbkę. Jeśli chemia jest znana, wprowadź ją w Kreatorze konfiguracji próbki / Chemii serii i ponów analizę.")
        self.description.setPlainText(publication_description_from_ranking(self.dataset, self.ranking))

    def selected_candidate(self) -> Optional[dict]:
        row = self.table.currentRow()
        return self.ranking[row] if 0 <= row < len(self.ranking) else None

    def _selection_changed(self):
        r = self.selected_candidate()
        if not r:
            return
        lines = [
            f"Kandydat: {r.get('phase','')}",
            f"Etykieta publikacyjna: {r.get('publication_label','')}",
            f"Źródło: {r.get('source','')}",
            f"Indeks dowodów: {float(r.get('ai_evidence_score_percent',0)):.1f}/100 ({r.get('confidence_band','')}); nie jest prawdopodobieństwem",
            f"Wsparcie względne wśród pokazanych kandydatów: {float(r.get('relative_support_percent',0)):.1f}%",
            f"Dopasowane refleksy: {r.get('matched_reflections',0)}",
            f"RMS Δ2θ: {r.get('rms_2theta_deg','')}",
            f"Wzór / kontekst COD: {r.get('formula','')}",
            f"Uzasadnienie: {r.get('reason','')}",
            "\nWażne: indeks dowodów AI nie jest prawdopodobieństwem fazy i nie zastępuje potwierdzenia eksperckiego ani analizy Rietvelda."
        ]
        self.detail.setPlainText("\n".join(lines))

    def confirm_selected(self):
        r = self.selected_candidate()
        if not r:
            QMessageBox.information(self, "Ocena faz wspomagana AI", "Najpierw wybierz jednego kandydata."); return
        phase = str(r.get("phase", "") or "").strip()
        pub = str(r.get("publication_label", "") or phase).strip()
        if not phase:
            return
        self.dataset.metadata["ai_user_confirmed_candidate"] = {
            "phase": phase, "publication_label": pub, "source": r.get("source", ""),
            "ai_evidence_score_percent": r.get("ai_evidence_score_percent", 0.0),
            "model": AI_EVIDENCE_MODEL_NAME,
        }
        self.dataset.metadata_sources["ai_user_confirmed_candidate"] = "user confirmation in AI-assisted review"
        aliases = self.dataset.metadata.get("publication_phase_aliases", {})
        if not isinstance(aliases, dict): aliases = {}
        aliases[str(r.get("phase", phase))] = pub
        if r.get("plot_phase"): aliases[str(r.get("plot_phase"))] = pub
        self.dataset.metadata["publication_phase_aliases"] = aliases
        self.dataset.metadata_sources["publication_phase_aliases"] = "user confirmation in AI-assisted review"
        self.description.setPlainText(publication_description_from_ranking(self.dataset, self.ranking))
        QMessageBox.information(self, "Ocena faz wspomagana AI", f"Zatwierdzono do wyświetlania publikacyjnego:\n{pub}\n\nBazowe przypisania pików COD/bazy wewnętrznej nie zostały nadpisane.")
        self.rebuild()
        parent = self.parent()
        if parent is not None:
            try:
                parent.redraw(); parent.update_table(); parent.update_summary(); parent.mark_project_dirty()
            except Exception:
                pass

    def clear_confirmation(self):
        self.dataset.metadata.pop("ai_user_confirmed_candidate", None)
        self.dataset.metadata.pop("publication_phase_aliases", None)
        self.dataset.metadata_sources.pop("ai_user_confirmed_candidate", None)
        self.dataset.metadata_sources.pop("publication_phase_aliases", None)
        self.rebuild()
        parent = self.parent()
        if parent is not None:
            try:
                parent.redraw(); parent.update_summary(); parent.mark_project_dirty()
            except Exception:
                pass

    def copy_description(self):
        QApplication.clipboard().setText(self.description.toPlainText())


# ---------------------- optional COD Offline sync worker --------------------
class CODOfflineWorker(QThread):
    progress_changed = pyqtSignal(int, str)
    task_completed = pyqtSignal(object)
    task_failed = pyqtSignal(str)

    def __init__(self, task: str, db_path: Path, args: dict, parent=None):
        super().__init__(parent)
        self.task = task
        self.db_path = Path(db_path)
        self.args = args
        self._cancelled = False

    def request_cancel(self):
        self._cancelled = True

    def _cancel_requested(self):
        return self._cancelled

    def run(self):
        try:
            kwargs = {"cancelled": self._cancel_requested, "progress": self.progress_changed.emit}
            if self.task == "api_elements":
                result = sync_api(self.db_path, elements=self.args["elements"],
                                  strict=self.args.get("strict", False), **kwargs)
            elif self.task == "api_ids":
                result = sync_api(self.db_path, ids=self.args["ids"], **kwargs)
            elif self.task == "archive":
                result = import_zip_archive(self.db_path, self.args["archive"], **kwargs)
            elif self.task == "bulk_download":
                archive = Path(self.args["archive"])
                if not archive.exists():
                    download_full_cod_zip(archive, **kwargs)
                if self._cancel_requested():
                    raise CODSyncError("Anulowano po pobraniu archiwum COD.")
                result = import_zip_archive(self.db_path, archive, **kwargs)
            else:
                raise CODSyncError("Nieznany rodzaj zadania COD Offline.")
            self.task_completed.emit(result)
        except Exception as exc:
            self.task_failed.emit(f"{type(exc).__name__}: {exc}")


# --------------------------- database dialog ---------------------------
class PhaseDatabaseDialog(QDialog):
    def __init__(self, db: PhaseDatabase, cod_db: CODDatabase, settings: QSettings, parent=None):
        super().__init__(parent)
        self.db = db
        self.cod_db = cod_db
        self.settings = settings
        self.setWindowTitle("Bazy faz – wzorce lokalne i COD")
        self.resize(1250, 760)
        root = QVBoxLayout(self)

        self.tabs = QTabWidget()
        root.addWidget(self.tabs, 1)

        # ---------------- local reference DB ----------------
        local_tab = QWidget(); layout = QVBoxLayout(local_tab)
        info = QLabel(
            "Lokalna, edytowalna baza wzorców faz. Jest używana tylko wtedy, gdy w oknie głównym jako źródło identyfikacji wybrano bazę wewnętrzną lub bazę wewnętrzną + COD. "
            "W trybie „Tylko zewnętrzna baza COD” cała biblioteka wewnętrzna jest pomijana przy automatycznym przypisywaniu faz. Rekordy COD można nadal importować tutaj do późniejszego użycia."
        )
        info.setWordWrap(True); layout.addWidget(info)
        self.search = QLineEdit(); self.search.setPlaceholderText("Szukaj lokalnej fazy, wzoru, struktury, kategorii…")
        self.search.textChanged.connect(self.refresh); layout.addWidget(self.search)
        self.table = QTableWidget(0, 7)
        self.table.setHorizontalHeaderLabels(["ID","Faza","Wzór","Struktura","Grupa przestrzenna","Kategoria","Uwagi"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(6, QHeaderView.ResizeMode.Stretch)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        layout.addWidget(self.table, 1)
        btns = QHBoxLayout()
        add_btn = QPushButton("Dodaj fazę")
        del_btn = QPushButton("Usuń zaznaczone")
        peaks_btn = QPushButton("Pokaż piki referencyjne")
        add_btn.clicked.connect(self.add_phase); del_btn.clicked.connect(self.delete_selected); peaks_btn.clicked.connect(self.show_peaks)
        for b in (add_btn, del_btn, peaks_btn): btns.addWidget(b)
        btns.addStretch(1); layout.addLayout(btns)
        self.tabs.addTab(local_tab, "Lokalne wzorce faz")

        # ---------------- bundled measured reference patterns ----------------
        exp_tab = QWidget(); el = QVBoxLayout(exp_tab)
        exp_info = QLabel(
            "Zmierzono wzorce referencyjne zapisane w bazie wewnętrznej. Rekordy dokumentują dostarczone serie XRDML stali X5CrNi18-10 oraz Ti-6Al-4V. "
            "Nie są używane jako automatyczne etykiety faz i nie należy interpretować ich jako jednoznacznego klasyfikatora temperatury obróbki; automatyczna identyfikacja korzysta z modeli krystalograficznych faz."
        )
        exp_info.setWordWrap(True); el.addWidget(exp_info)
        self.exp_table = QTableWidget(0, 10)
        self.exp_table.setHorizontalHeaderLabels([
            "ID","Wzorzec","Materiał","Obróbka","Zinterpretowana faza","Struktura","λ [Å]","dopas. a [Å]","dopas. c [Å]","Uwagi"
        ])
        self.exp_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.exp_table.horizontalHeader().setSectionResizeMode(9, QHeaderView.ResizeMode.Stretch)
        self.exp_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        el.addWidget(self.exp_table, 1)
        exp_buttons = QHBoxLayout()
        exp_show = QPushButton("Pokaż zmierzone piki / hkl")
        exp_show.clicked.connect(self.show_experimental_reference)
        exp_buttons.addWidget(exp_show); exp_buttons.addStretch(1); el.addLayout(exp_buttons)
        self.tabs.addTab(exp_tab, "Wzorce pomiarowe")

        # ---------------- COD external DB ----------------
        cod_tab = QWidget(); cl = QVBoxLayout(cod_tab)
        top = QHBoxLayout()
        self.cod_status = QLabel("Baza COD: niewczytana")
        self.cod_status.setWordWrap(True)
        load_cod_btn = QPushButton("Wczytaj bazę COD .db3")
        disconnect_btn = QPushButton("Odłącz")
        load_cod_btn.clicked.connect(self.load_cod_database)
        disconnect_btn.clicked.connect(self.disconnect_cod)
        top.addWidget(load_cod_btn); top.addWidget(disconnect_btn); top.addWidget(self.cod_status, 1)
        cl.addLayout(top)

        # COD Offline synchronization; read-only CODDatabase connector remains unchanged.
        self._cod_worker = None
        self._cod_progress = None
        self.cod_filter_mode = QComboBox()
        self.cod_filter_mode.addItems(["Pierwiastki (wymagane)", "Numery COD"])
        self.cod_filter = QLineEdit()
        self.cod_filter.setPlaceholderText("np. Fe Cr Ni albo 1000000, 2002926")
        self.cod_filter.setText(str(self.settings.value("cod_sync_filter", "")))
        self.cod_strict = QCheckBox("Tylko podane pierwiastki")
        self.cod_strict.setChecked(True)
        self.cod_online_btn = QPushButton("Pobierz / aktualizuj COD Online → SQLite")
        self.cod_online_btn.clicked.connect(self.start_cod_online)
        sync_row = QHBoxLayout()
        sync_row.addWidget(self.cod_filter_mode)
        sync_row.addWidget(self.cod_filter, 1)
        sync_row.addWidget(self.cod_strict)
        sync_row.addWidget(self.cod_online_btn)
        cl.addLayout(sync_row)

        archive_row = QHBoxLayout()
        self.cod_archive_btn = QPushButton("Importuj lokalne archiwum CIF ZIP (pełna baza)")
        self.cod_archive_btn.clicked.connect(self.start_cod_archive)
        self.cod_bulk_btn = QPushButton("Pobierz pełne archiwum COD (duże, opcjonalne)")
        self.cod_bulk_btn.clicked.connect(self.start_cod_bulk_download)
        archive_row.addWidget(self.cod_archive_btn)
        archive_row.addWidget(self.cod_bulk_btn)
        archive_row.addStretch(1)
        cl.addLayout(archive_row)
        note = QLabel(
            "COD Offline 20.3.5: pierwsze pobranie wymaga internetu; ponowna praca jest lokalna. "
            "Filtr online ogranicza wynik do 30 000 struktur / 48 MB. "
            "Aktualizacja ponownie pobiera wybrane metadane (nie pełny indeks COD). "
            "Pełne archiwum jest bardzo duże; import może zająć dużo czasu i dysku. "
            "Tylko dane struktury — bez eksperymentalnych intensywności proszkowych."
        )
        note.setWordWrap(True)
        cl.addWidget(note)
        cod_info_btn = QPushButton("Skąd pobrać COD? (.db3 / COD Offline)")
        cod_info_btn.clicked.connect(self.show_cod_help)
        cl.addWidget(cod_info_btn)

        help_label = QLabel(
            "Search by COD ID, formula, chemical/mineral/common name, space group or DOI. "
            "The COD file is opened read-only. In Tylko zewnętrzna baza COD mode the program searches this .db3 directly and does not consult the internal phase database. "
            "Importing a COD record is optional and only makes it available for later internal-database workflows."
        )
        help_label.setWordWrap(True); cl.addWidget(help_label)

        search_row = QHBoxLayout()
        self.cod_search = QLineEdit(); self.cod_search.setPlaceholderText("np. 1524232, Cu Ni, iron, F m -3 m…")
        self.cod_limit = QSpinBox(); self.cod_limit.setRange(10, 1000); self.cod_limit.setValue(200)
        cod_search_btn = QPushButton("Szukaj w COD")
        cod_reset_btn = QPushButton("Przywróć domyślne")
        cod_search_btn.clicked.connect(self.search_cod)
        cod_reset_btn.clicked.connect(self.reset_cod_search_defaults)
        self.cod_search.returnPressed.connect(self.search_cod)
        search_row.addWidget(self.cod_search, 1); search_row.addWidget(QLabel("Maks. liczba wyników:")); search_row.addWidget(self.cod_limit); search_row.addWidget(cod_search_btn); search_row.addWidget(cod_reset_btn)
        cl.addLayout(search_row)

        self.cod_table = QTableWidget(0, 14)
        self.cod_table.setHorizontalHeaderLabels([
            "COD ID","Wzór","Nazwa chemiczna","Nazwa mineralna / zwyczajowa","Grupa przestrzenna","Nr GP",
            "a [Å]","b [Å]","c [Å]","α [°]","β [°]","γ [°]","Rok","DOI"
        ])
        self.cod_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.cod_table.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeMode.Stretch)
        self.cod_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.cod_table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.cod_table.setSortingEnabled(True)
        cl.addWidget(self.cod_table, 1)

        cod_buttons = QHBoxLayout()
        import_btn = QPushButton("Importuj wybrane struktury COD do wzorców lokalnych")
        details_btn = QPushButton("Pokaż szczegóły COD")
        import_btn.clicked.connect(self.import_selected_cod)
        details_btn.clicked.connect(self.show_cod_details)
        cod_buttons.addWidget(import_btn); cod_buttons.addWidget(details_btn); cod_buttons.addStretch(1)
        cl.addLayout(cod_buttons)
        self.tabs.addTab(cod_tab, "COD .db3")

        close_row = QHBoxLayout(); close_row.addStretch(1)
        close_btn = QPushButton("Zamknij"); close_btn.clicked.connect(self.accept); close_row.addWidget(close_btn)
        root.addLayout(close_row)

        self.refresh()
        self.refresh_experimental_references()
        self.update_cod_status()

    def show_cod_help(self):
        from xrd_help_content import COD_HELP_HTML, show_help_dialog
        show_help_dialog(self, "COD — pobieranie i baza offline", COD_HELP_HTML)

    def _cod_local_target(self) -> Path:
        # Never modify user-supplied external SQLite databases. Always store
        # downloaded data under the application-specific COD Offline location.
        return default_cod_path()

    def start_cod_online(self):
        if self._cod_worker is not None and self._cod_worker.isRunning():
            return
        raw = self.cod_filter.text().strip()
        try:
            if self.cod_filter_mode.currentIndex() == 0:
                symbols = parse_elements(raw)
                task, args = "api_elements", {"elements": symbols, "strict": self.cod_strict.isChecked()}
            else:
                ids = parse_ids(raw)
                task, args = "api_ids", {"ids": ids}
            self.settings.setValue("cod_sync_filter", raw)
        except ValueError as exc:
            QMessageBox.warning(self, "Filtr COD", str(exc))
            return
        self._start_cod_worker(task, args)

    def start_cod_archive(self):
        archive, _ = QFileDialog.getOpenFileName(
            self, "Wybierz archiwum COD z plikami CIF", str(self._cod_local_target().parent),
            "Archiwa ZIP z CIF (*.zip);;Wszystkie pliki (*.*)")
        if archive:
            self._start_cod_worker("archive", {"archive": str(Path(archive))})

    def start_cod_bulk_download(self):
        if QMessageBox.question(
            self, "Pełne archiwum COD",
            "Oficjalne archiwum COD może liczyć wiele GB. Pobranie i indeksowanie "
            "wymagają dużo miejsca na dysku i mogą długo trwać. "
            "Dane zostaną zapisane lokalnie, bez nadpisania zdjęć ani projektów. "
            "Dla zwykłej pracy zalecamy podzbiór według składu. Kontynuować?"
        ) != QMessageBox.StandardButton.Yes:
            return
        name = f"cod-cifs-mysql_{datetime.now():%Y%m%d}.zip"
        self._start_cod_worker("bulk_download", {"archive": str(self._cod_local_target().parent / name)})

    def _start_cod_worker(self, task, args):
        if self._cod_worker is not None and self._cod_worker.isRunning():
            return
        destination = self._cod_local_target()
        self._cod_prior_path = Path(self.cod_db.path) if self.cod_db.is_open else None
        # The old COD connector is read-only; close it while writing the same file.
        if self.cod_db.is_open and Path(self.cod_db.path).resolve() == destination.resolve():
            self.cod_db.close()
        self._cod_worker = CODOfflineWorker(task, destination, args, self)
        self._cod_progress = QProgressDialog("Przygotowanie COD Offline...", "Anuluj", 0, 100, self)
        self._cod_progress.setWindowTitle("Aktualizacja bazy COD")
        self._cod_progress.setMinimumDuration(0)
        self._cod_progress.setWindowModality(Qt.WindowModality.WindowModal)
        self._cod_progress.setAutoClose(False)
        self._cod_progress.canceled.connect(self._cod_worker.request_cancel)
        self._cod_worker.progress_changed.connect(self._cod_progress_event)
        self._cod_worker.task_completed.connect(self._cod_done)
        self._cod_worker.task_failed.connect(self._cod_failure)
        self._cod_worker.finished.connect(self._cod_finished)
        self._cod_progress.show()
        self._cod_worker.start()

    def _cod_progress_event(self, value, message):
        if self._cod_progress is not None:
            self._cod_progress.setLabelText(message)
            if value >= 0:
                self._cod_progress.setValue(value)

    def _cod_done(self, result):
        try:
            self.cod_db.open(Path(result["path"]))
            self.settings.setValue("cod_db_path", str(Path(result["path"]).resolve()))
            self.settings.setValue("cod_db_dir", str(Path(result["path"]).resolve().parent))
            self.update_cod_status()
            QMessageBox.information(
                self, "COD Offline – aktualizacja zakończona",
                f"Źródło: {result['source']}\n"
                f"Zaimportowano / zaktualizowano: {result['imported']:,}\n"
                f"Łączna liczba struktur w SQLite: {result['total']:,}\n"
                f"Baza: {result['path']}\n\n"
                "Dane można teraz przeszukiwać bez połączenia z internetem.")
        except Exception as exc:
            QMessageBox.critical(self, "COD Offline", f"Zapis udał się, ale nie można otworzyć SQLite: {exc}")

    def _cod_failure(self, message):
        # Restore previously loaded reference database after a rejected update.
        if not self.cod_db.is_open and self._cod_prior_path is not None:
            try:
                self.cod_db.open(self._cod_prior_path)
            except Exception:
                pass
        self.update_cod_status()
        QMessageBox.warning(self, "COD Offline – problem", message)

    def _cod_finished(self):
        if self._cod_progress is not None:
            self._cod_progress.close()
            self._cod_progress.deleteLater()
            self._cod_progress = None
        if self._cod_worker is not None:
            self._cod_worker.deleteLater()
            self._cod_worker = None
        self.update_cod_status()

    def accept(self):
        if self._cod_worker is not None and self._cod_worker.isRunning():
            QMessageBox.information(self, "COD Offline", "Trwa pobieranie lub import. Anuluj zadanie i poczekaj na jego zakończenie przed zamknięciem okna.")
            return
        super().accept()

    def reject(self):
        if self._cod_worker is not None and self._cod_worker.isRunning():
            self._cod_worker.request_cancel()
            return
        super().reject()

    def closeEvent(self, event):
        if self._cod_worker is not None and self._cod_worker.isRunning():
            self._cod_worker.request_cancel()
            event.ignore()
            return
        super().closeEvent(event)

    def reset_cod_search_defaults(self):
        self.cod_search.clear()
        self.cod_limit.setValue(200)
        self.cod_table.setRowCount(0)

    def refresh(self):
        q = self.search.text().strip().lower()
        rows = []
        for r in self.db.all_phases():
            if not q or q in " ".join(str(v or "") for v in r).lower(): rows.append(r)
        self.table.setRowCount(len(rows))
        for i,row in enumerate(rows):
            for j,v in enumerate(row): self.table.setItem(i,j,QTableWidgetItem(str(v or "")))

    def refresh_experimental_references(self):
        rows = self.db.all_experimental_references()
        self.exp_table.setRowCount(len(rows))
        for i, row in enumerate(rows):
            for j, v in enumerate(row):
                if j in (6,7,8) and v not in (None, ""):
                    try:
                        txt = f"{float(v):.6g}"
                    except Exception:
                        txt = str(v)
                else:
                    txt = str(v or "")
                self.exp_table.setItem(i, j, QTableWidgetItem(txt))

    def selected_experimental_reference_id(self):
        row = self.exp_table.currentRow()
        if row < 0:
            return None
        try:
            return int(self.exp_table.item(row, 0).text())
        except Exception:
            return None

    def show_experimental_reference(self):
        rid = self.selected_experimental_reference_id()
        if rid is None:
            return
        meta = None
        for row in self.db.all_experimental_references():
            if int(row[0]) == int(rid):
                meta = row; break
        peaks = self.db.experimental_reference_peaks(rid)
        if meta is None:
            return
        _, name, material, treatment, phase_label, structure, wave, a_fit, c_fit, notes = meta
        lines = [
            f"Wzorzec: {name}", f"Materiał: {material}", f"Obróbka: {treatment}",
            f"Interpretacja: {phase_label}", f"Struktura: {structure}", f"Długość fali: {wave} Å",
            f"Dopasowana sieć: a={a_fit if a_fit is not None else '-'} Å, c={c_fit if c_fit is not None else '-'} Å",
            "", str(notes or ""), "", "2θ [°]\td [Å]\thkl\tIrel [%]"
        ]
        for tt, dA, hkl, rel in peaks:
            lines.append(f"{float(tt):.4f}\t{float(dA):.5f}\t{hkl}\t{float(rel):.1f}")
        dlg = QDialog(self); dlg.setWindowTitle("Wewnętrzny wzorzec pomiarowy")
        lay = QVBoxLayout(dlg); te = QTextEdit(); te.setReadOnly(True); te.setPlainText("\n".join(lines)); lay.addWidget(te)
        bb = QDialogButtonBox(QDialogButtonBox.StandardButton.Close); bb.button(QDialogButtonBox.StandardButton.Close).setText("Zamknij"); bb.rejected.connect(dlg.reject); bb.accepted.connect(dlg.accept); lay.addWidget(bb)
        dlg.resize(760, 620); dlg.exec()

    def selected_phase_id(self):
        row = self.table.currentRow()
        if row < 0: return None
        try: return int(self.table.item(row,0).text())
        except Exception: return None

    def show_peaks(self):
        pid = self.selected_phase_id()
        if pid is None: return
        rows = self.db.phase_peaks(pid)
        model, a_ref, c_ref, fit_enabled = self.db.phase_lattice_model(pid)
        if fit_enabled:
            seed = f"Początkowy parametr sieci jest tylko wartością startową: model={model}, a_ref={a_ref if a_ref is not None else '-'} Å"
            if c_ref is not None:
                seed += f", c_ref={c_ref} Å"
            seed += "\nRzeczywiste parametry sieci są dopasowywane do zmierzonych pozycji pików XRD przed oceną fazy.\n\n"
        else:
            seed = f"Stały wzorzec pozycyjny: model={model}.\n\n"
        text = seed + "d [Å]\thkl\tIrel\n" + "\n".join(f"{d:.5f}\t{hkl}\t{rel:g}" for d,hkl,rel in rows)
        dlg = QDialog(self); dlg.setWindowTitle("Piki referencyjne")
        lay = QVBoxLayout(dlg); te = QTextEdit(); te.setReadOnly(True); te.setPlainText(text); lay.addWidget(te)
        bb = QDialogButtonBox(QDialogButtonBox.StandardButton.Close); bb.button(QDialogButtonBox.StandardButton.Close).setText("Zamknij"); bb.rejected.connect(dlg.reject); bb.accepted.connect(dlg.accept); lay.addWidget(bb)
        dlg.resize(550,500); dlg.exec()

    def add_phase(self):
        name, ok = QInputDialog.getText(self, "Dodaj fazę", "Nazwa fazy:")
        if not ok or not name.strip(): return
        formula, _ = QInputDialog.getText(self, "Dodaj fazę", "Wzór / chemia:")
        structure, _ = QInputDialog.getText(self, "Dodaj fazę", "Struktura krystaliczna (np. FCC, BCC, HCP):")
        sg, _ = QInputDialog.getText(self, "Dodaj fazę", "Grupa przestrzenna:")
        category, _ = QInputDialog.getText(self, "Dodaj fazę", "Kategoria:")
        notes, _ = QInputDialog.getMultiLineText(self, "Dodaj fazę", "Uwagi:")
        peak_text, ok = QInputDialog.getMultiLineText(self, "Piki referencyjne", "Wprowadź jeden pik w każdym wierszu jako: d_A,hkl,względna_intensywność\nPrzykład: 2.071,111,100")
        if not ok: return
        peaks=[]
        try:
            for line in peak_text.splitlines():
                if not line.strip(): continue
                p=[x.strip() for x in line.replace(";",",").split(",")]
                dA=float(p[0]); hkl=p[1] if len(p)>1 else ""; rel=float(p[2]) if len(p)>2 else 100.0
                peaks.append((dA,hkl,rel))
            if not peaks: raise ValueError("Wymagany jest co najmniej jeden pik referencyjny.")
            self.db.add_phase(name.strip(), formula, structure, sg, category, notes, peaks); self.refresh()
        except Exception as e: QMessageBox.critical(self,"Baza danych",str(e))

    def delete_selected(self):
        pid = self.selected_phase_id()
        if pid is None: return
        if QMessageBox.question(self,"Usuń fazę","Usunąć zaznaczoną fazę z lokalnej bazy danych?") == QMessageBox.StandardButton.Yes:
            self.db.delete_phase(pid); self.refresh()

    def update_cod_status(self):
        if self.cod_db.is_open:
            self.cod_status.setText(f"Baza COD: {self.cod_db.path}  |  rekordy: {self.cod_db.record_count:,}")
        else:
            self.cod_status.setText("Baza COD: niewczytana")

    def load_cod_database(self):
        start_dir = self.settings.value("cod_db_dir", "")
        path, _ = QFileDialog.getOpenFileName(self, "Wczytaj bazę COD SQLite", start_dir, "COD SQLite (*.db3 *.sqlite *.sqlite3 *.db);;Wszystkie pliki (*.*)")
        if not path: return
        try:
            self.cod_db.open(Path(path))
            self.settings.setValue("cod_db_path", str(Path(path).resolve()))
            self.settings.setValue("cod_db_dir", str(Path(path).resolve().parent))
            self.update_cod_status()
            QMessageBox.information(self, "Baza COD", f"Wczytano w trybie tylko do odczytu:\n{path}\n\nRekordy: {self.cod_db.record_count:,}")
        except Exception as e:
            QMessageBox.critical(self, "Baza COD", str(e))

    def disconnect_cod(self):
        self.cod_db.close(); self.settings.remove("cod_db_path"); self.cod_table.setRowCount(0); self.update_cod_status()

    @staticmethod
    def _clean(v):
        if v is None: return ""
        s=str(v).strip(); return "" if s in ("\\N","?",".") else s

    def search_cod(self):
        if not self.cod_db.is_open:
            QMessageBox.information(self, "Baza COD", "Najpierw wczytaj bazę COD .db3."); return
        try:
            rows = self.cod_db.search(self.cod_search.text(), self.cod_limit.value())
            self.cod_table.setSortingEnabled(False); self.cod_table.setRowCount(len(rows))
            for r,row in enumerate(rows):
                (file_id, formula, chemname, mineral, commonname, sg, sgno, a,b,c,alpha,beta,gamma,doi,journal,year,method)=row
                min_common = self._clean(mineral) or self._clean(commonname)
                vals=[file_id,self._clean(formula),self._clean(chemname),min_common,self._clean(sg),sgno,a,b,c,alpha,beta,gamma,year,self._clean(doi)]
                for cidx,v in enumerate(vals): self.cod_table.setItem(r,cidx,QTableWidgetItem("" if v is None else str(v)))
            self.cod_table.setSortingEnabled(True)
            if not rows: QMessageBox.information(self, "Wyszukiwanie COD", "Nie znaleziono pasujących rekordów COD.")
        except Exception as e:
            QMessageBox.critical(self, "Wyszukiwanie COD", f"{e}\n\n{traceback.format_exc()}")

    def selected_cod_ids(self):
        ids=set()
        for item in self.cod_table.selectedItems():
            row=item.row(); cell=self.cod_table.item(row,0)
            if cell:
                try: ids.add(int(cell.text()))
                except Exception: pass
        return sorted(ids)

    def import_selected_cod(self):
        ids=self.selected_cod_ids()
        if not ids:
            QMessageBox.information(self, "Import COD", "Najpierw zaznacz jeden lub więcej wierszy COD."); return
        imported=[]; skipped=[]
        for cod_id in ids:
            try:
                ref=self.cod_db.entry_to_reference(cod_id)
                try:
                    self.db.add_phase(ref["name"],ref["formula"],ref["structure"],ref["space_group"],ref["category"],ref["notes"],ref["peaks"])
                    imported.append(f"{cod_id} ({len(ref['peaks'])} reflections)")
                except sqlite3.IntegrityError:
                    skipped.append(f"{cod_id}: already imported")
            except Exception as e:
                skipped.append(f"{cod_id}: {e}")
        self.refresh()
        msg=""
        if imported: msg += "Imported:\n" + "\n".join(imported)
        if skipped: msg += ("\n\n" if msg else "") + "Not imported:\n" + "\n".join(skipped)
        QMessageBox.information(self, "Import COD", msg or "Nic nie zaimportowano.")

    def show_cod_details(self):
        ids=self.selected_cod_ids()
        if not ids: return
        row=self.cod_db.get_entry(ids[0])
        if not row: return
        labels=["COD ID","Formula","Chemical name","Mineral","Common name","Space group","SG number","a [Å]","b [Å]","c [Å]","α [°]","β [°]","γ [°]","DOI","Journal","Year","Method"]
        text="\n".join(f"{lab}: {self._clean(val)}" for lab,val in zip(labels,row))
        dlg=QDialog(self); dlg.setWindowTitle(f"COD {ids[0]} – szczegóły"); lay=QVBoxLayout(dlg)
        te=QTextEdit(); te.setReadOnly(True); te.setPlainText(text); lay.addWidget(te)
        bb=QDialogButtonBox(QDialogButtonBox.StandardButton.Close); bb.rejected.connect(dlg.reject); bb.accepted.connect(dlg.accept); lay.addWidget(bb)
        dlg.resize(650,520); dlg.exec()


# --------------------------- periodic table selector ---------------------------
ELEMENT_SYMBOLS = [
    "H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar",
    "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
    "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
    "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
    "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac",
    "Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh",
    "Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"
]
ELEMENT_Z = {sym: i + 1 for i, sym in enumerate(ELEMENT_SYMBOLS)}

# Main-table coordinates (period, group). Lantanowce/actinides are shown in separate rows.
PERIODIC_POSITIONS = {
    "H":(1,1), "He":(1,18),
    "Li":(2,1), "Be":(2,2), "B":(2,13), "C":(2,14), "N":(2,15), "O":(2,16), "F":(2,17), "Ne":(2,18),
    "Na":(3,1), "Mg":(3,2), "Al":(3,13), "Si":(3,14), "P":(3,15), "S":(3,16), "Cl":(3,17), "Ar":(3,18),
    "K":(4,1), "Ca":(4,2), "Sc":(4,3), "Ti":(4,4), "V":(4,5), "Cr":(4,6), "Mn":(4,7), "Fe":(4,8), "Co":(4,9), "Ni":(4,10), "Cu":(4,11), "Zn":(4,12), "Ga":(4,13), "Ge":(4,14), "As":(4,15), "Se":(4,16), "Br":(4,17), "Kr":(4,18),
    "Rb":(5,1), "Sr":(5,2), "Y":(5,3), "Zr":(5,4), "Nb":(5,5), "Mo":(5,6), "Tc":(5,7), "Ru":(5,8), "Rh":(5,9), "Pd":(5,10), "Ag":(5,11), "Cd":(5,12), "In":(5,13), "Sn":(5,14), "Sb":(5,15), "Te":(5,16), "I":(5,17), "Xe":(5,18),
    "Cs":(6,1), "Ba":(6,2), "Hf":(6,4), "Ta":(6,5), "W":(6,6), "Re":(6,7), "Os":(6,8), "Ir":(6,9), "Pt":(6,10), "Au":(6,11), "Hg":(6,12), "Tl":(6,13), "Pb":(6,14), "Bi":(6,15), "Po":(6,16), "At":(6,17), "Rn":(6,18),
    "Fr":(7,1), "Ra":(7,2), "Rf":(7,4), "Db":(7,5), "Sg":(7,6), "Bh":(7,7), "Hs":(7,8), "Mt":(7,9), "Ds":(7,10), "Rg":(7,11), "Cn":(7,12), "Nh":(7,13), "Fl":(7,14), "Mc":(7,15), "Lv":(7,16), "Ts":(7,17), "Og":(7,18),
}
LANTHANIDES = ["La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu"]
ACTINIDES = ["Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr"]

METAL_ELEMENTS = set(
    "Li Be Na Mg Al K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd Cs Ba "
    "La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn "
    "Ga In Sn Tl Pb Bi Po Nh Fl Mc Lv".split()
)
COMMON_WELD_ELEMENTS = set("Fe Cr Ni Mo Cu Mn Si C Al Ti Nb V W Co Zn Mg N O S P".split())


def gemmi_diagnostic_message(context: str = "Wyszukiwanie COD") -> str:
    ensure_gemmi()
    if GEMMI_OK:
        version = getattr(gemmi, "__version__", "unknown")
        return f"{context}: biblioteka gemmi {version} jest dostępna."
    return (
        f"{context} wymaga biblioteki „gemmi”, ale uruchomiony Python nie może jej zaimportować.\n\n"
        f"Python używany przez XRD Analyzer:\n{sys.executable}\n\n"
        f"Błąd importu:\n{GEMMI_IMPORT_ERROR or 'nieznany błąd'}\n\n"
        "Zainstaluj gemmi w TYM interpreterze i uruchom ponownie XRD Analyzer:\n"
        f'"{sys.executable}" -m pip install --upgrade gemmi\n\n'
        "Jeżeli gemmi jest już zainstalowane, prawdopodobnie zostało zainstalowane w innym środowisku Pythona."
    )


class PeriodicTableDialog(QDialog):
    """Interactive periodic table used to build COD chemistry filters."""
    def __init__(self, selected=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Wybierz dozwolone pierwiastki – układ okresowy")
        self.resize(1150, 610)
        self.buttons = {}
        selected = set(selected or [])

        root = QVBoxLayout(self)
        help_label = QLabel(
            "Zaznacz wszystkie pierwiastki, które mogą występować w rozważanych fazach COD. "
            "Pierwiastek zaznaczony jest dozwolony; brak zaznaczenia wyklucza rekordy COD zawierające ten pierwiastek."
        )
        help_label.setWordWrap(True)
        root.addWidget(help_label)

        quick = QHBoxLayout()
        clear_btn = QPushButton("Wyczyść")
        metals_btn = QPushButton("Wybierz metale")
        weld_btn = QPushButton("Typowe pierwiastki spawalnicze")
        select_current_btn = QPushButton("Zachowaj bieżący wybór")
        clear_btn.clicked.connect(lambda: self.set_selection(set()))
        metals_btn.clicked.connect(lambda: self.set_selection(METAL_ELEMENTS))
        weld_btn.clicked.connect(lambda: self.set_selection(COMMON_WELD_ELEMENTS))
        select_current_btn.clicked.connect(lambda: None)
        quick.addWidget(clear_btn)
        quick.addWidget(metals_btn)
        quick.addWidget(weld_btn)
        quick.addStretch(1)
        self.count_label = QLabel()
        quick.addWidget(self.count_label)
        root.addLayout(quick)

        frame = QFrame()
        grid = QGridLayout(frame)
        grid.setHorizontalSpacing(4)
        grid.setVerticalSpacing(4)

        # group numbers
        for g in range(1, 19):
            lab = QLabel(str(g)); lab.setAlignment(Qt.AlignmentFlag.AlignCenter)
            grid.addWidget(lab, 0, g)
        # period numbers
        for period in range(1, 8):
            lab = QLabel(str(period)); lab.setAlignment(Qt.AlignmentFlag.AlignCenter)
            grid.addWidget(lab, period, 0)

        def add_element(sym, row, col):
            btn = QToolButton()
            btn.setText(f"{ELEMENT_Z[sym]}\n{sym}")
            btn.setCheckable(True)
            btn.setChecked(sym in selected)
            btn.setMinimumSize(52, 48)
            btn.setToolTip(f"Liczba atomowa {ELEMENT_Z[sym]} — {sym}")
            btn.toggled.connect(self.update_count)
            grid.addWidget(btn, row, col)
            self.buttons[sym] = btn

        for sym, (row, col) in PERIODIC_POSITIONS.items():
            add_element(sym, row, col)

        # placeholders for f-block in the main table
        for row, text in ((6, "57–71"), (7, "89–103")):
            lab = QLabel(text); lab.setAlignment(Qt.AlignmentFlag.AlignCenter)
            grid.addWidget(lab, row, 3)

        # f-block rows
        grid.addWidget(QLabel("Lantanowce"), 8, 0, 1, 3)
        grid.addWidget(QLabel("Aktynowce"), 9, 0, 1, 3)
        for i, sym in enumerate(LANTHANIDES):
            add_element(sym, 8, 3 + i)
        for i, sym in enumerate(ACTINIDES):
            add_element(sym, 9, 3 + i)

        root.addWidget(frame, 1)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.button(QDialogButtonBox.StandardButton.Ok).setText("Zastosuj")
        buttons.button(QDialogButtonBox.StandardButton.Cancel).setText("Anuluj")
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        root.addWidget(buttons)
        self.update_count()

    def selected_elements(self):
        return {sym for sym, btn in self.buttons.items() if btn.isChecked()}

    def set_selection(self, symbols):
        symbols = set(symbols)
        for sym, btn in self.buttons.items():
            btn.blockSignals(True)
            btn.setChecked(sym in symbols)
            btn.blockSignals(False)
        self.update_count()

    def update_count(self):
        chosen = self.selected_elements()
        self.count_label.setText(f"Wybrano: {len(chosen)}")



# --------------------------- sample setup / chemistry wizard ---------------------------
class SampleSetupWizard(QDialog):
    """User-guided sample description and chemistry input.

    The wizard deliberately does not infer elemental composition from the diffraction
    pattern itself. Material-name suggestions are optional presets only; user-entered
    chemistry has priority and is saved with the dataset/project.
    """
    MODE_ITEMS = [
        ("Dokładny skład", "exact"),
        ("Zakresy składu", "ranges"),
        ("Znana kolejność: od największej do najmniejszej zawartości", "order"),
        ("Znana tylko lista pierwiastków", "elements"),
        ("Skład nieznany – tylko detekcja pików", "unknown"),
    ]

    def __init__(self, dataset: Dataset, parent=None):
        super().__init__(parent)
        self.dataset = dataset
        self._preset_was_applied = str((dataset.metadata or {}).get("material_preset_confirmed", "") or "")
        self.setWindowTitle("Kreator konfiguracji próbki")
        self.resize(900, 680)
        root = QVBoxLayout(self)

        intro = QLabel(
            "Podaj to, co wiadomo o próbce. Dokładne udziały procentowe nie są obowiązkowe. "
            "Jeżeli znasz tylko kolejność, umieść pierwiastek dominujący jako pierwszy, a najmniejszy znany składnik jako ostatni. "
            "Dane te ograniczają i porządkują kandydatów faz; XRD nie jest traktowane jako pomiar składu pierwiastkowego."
        )
        intro.setWordWrap(True); root.addWidget(intro)

        form_box = QGroupBox("1. Próbka")
        form = QFormLayout(form_box)
        self.name_edit = QLineEdit(dataset.name)
        self.material_edit = QLineEdit(dataset.material_type)
        self.preset_combo = QComboBox()
        self.preset_combo.addItem("Bez wbudowanego presetu materiału", "")
        for preset_name in MATERIAL_SETUP_PRESETS:
            self.preset_combo.addItem(preset_name, preset_name)
        self.apply_preset_btn = QPushButton("Zastosuj preset")
        preset_row = QWidget(); preset_lay = QHBoxLayout(preset_row); preset_lay.setContentsMargins(0,0,0,0)
        preset_lay.addWidget(self.preset_combo, 1); preset_lay.addWidget(self.apply_preset_btn)
        self.kind_combo = QComboBox(); self.kind_combo.addItems([
            "Nieznana / próbka ogólna", "Materiał podstawowy", "Drut spawalniczy / spoiwo", "Stopiwo",
            "Spoina różnoimienna", "Proszek / powłoka", "Inne"
        ])
        self.base_a_edit = QLineEdit(str(dataset.metadata.get("base_material_A", "")))
        self.base_b_edit = QLineEdit(str(dataset.metadata.get("base_material_B", "")))
        self.filler_edit = QLineEdit(str(dataset.metadata.get("filler_material", "")))
        self.matrix_combo = QComboBox()
        for label, code in EXPECTED_MATRIX_ITEMS:
            self.matrix_combo.addItem(label, code)
        saved_matrix = str(dataset.metadata.get("expected_matrix_structure", "auto") or "auto")
        for i in range(self.matrix_combo.count()):
            if self.matrix_combo.itemData(i) == saved_matrix:
                self.matrix_combo.setCurrentIndex(i); break
        self.matrix_combo.setToolTip(
            "Opcjonalna wskazówka krystalograficzna do rankingu. Jest to miękka preferencja, a nie wymuszone przypisanie fazy. "
            "Dla znanego 316L wybierz γ-austenit – osnowa Fe, FCC."
        )
        form.addRow("Nazwa próbki (opcjonalnie):", self.name_edit)
        form.addRow("Materiał / gatunek (opcjonalnie):", self.material_edit)
        form.addRow("Wbudowany preset materiału:", preset_row)
        form.addRow("Typ próbki:", self.kind_combo)
        form.addRow("Materiał podstawowy A (opcjonalnie):", self.base_a_edit)
        form.addRow("Materiał podstawowy B (opcjonalnie):", self.base_b_edit)
        form.addRow("Spoiwo / drut (opcjonalnie):", self.filler_edit)
        form.addRow("Oczekiwana osnowa / struktura:", self.matrix_combo)
        root.addWidget(form_box)

        chem_box = QGroupBox("2. Skład chemiczny")
        chem = QVBoxLayout(chem_box)
        mode_row = QHBoxLayout()
        self.mode_combo = QComboBox()
        for label, code in self.MODE_ITEMS:
            self.mode_combo.addItem(label, code)
        saved_mode = str(dataset.metadata.get("composition_mode", "elements" if dataset.metadata.get("allowed_elements") else "unknown"))
        for i in range(self.mode_combo.count()):
            if self.mode_combo.itemData(i) == saved_mode:
                self.mode_combo.setCurrentIndex(i); break
        self.unit_combo = QComboBox(); self.unit_combo.addItems(["wt.%", "at.%"])
        self.unit_combo.setCurrentText(str(dataset.metadata.get("composition_unit", "wt.%")))
        mode_row.addWidget(QLabel("Poziom wiedzy:")); mode_row.addWidget(self.mode_combo, 1)
        mode_row.addWidget(QLabel("Jednostka:")); mode_row.addWidget(self.unit_combo)
        chem.addLayout(mode_row)

        self.comp_table = QTableWidget(0, 5)
        self.comp_table.setHorizontalHeaderLabels(["Pierwiastek", "Wartość", "Min.", "Maks.", "Uwagi"])
        self.comp_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        self.comp_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        self.comp_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
        self.comp_table.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeMode.ResizeToContents)
        self.comp_table.horizontalHeader().setSectionResizeMode(4, QHeaderView.ResizeMode.Stretch)
        self.comp_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        chem.addWidget(self.comp_table, 1)

        br = QHBoxLayout()
        self.periodic_btn = QPushButton("Dodaj z układu okresowego…")
        self.add_btn = QPushButton("Dodaj wiersz")
        self.del_btn = QPushButton("Usuń wiersz")
        self.up_btn = QPushButton("W górę")
        self.down_btn = QPushButton("W dół")
        self.suggest_btn = QPushButton("Użyj niezweryfikowanej sugestii z nazwy…")
        self.suggest_btn.setToolTip("Sugestia opiera się wyłącznie na tekście/nazwie. Nic nie zostanie zastosowane, dopóki nie użyjesz tego przycisku i nie zapiszesz kreatora.")
        for b in (self.periodic_btn, self.add_btn, self.del_btn, self.up_btn, self.down_btn, self.suggest_btn): br.addWidget(b)
        chem.addLayout(br)
        order_help = QLabel("W trybie „Znana kolejność” kolejność wierszy oznacza ranking chemiczny: pierwszy wiersz = największa znana zawartość.")
        order_help.setWordWrap(True); chem.addWidget(order_help)
        root.addWidget(chem_box, 1)

        self.summary = QLabel(); self.summary.setWordWrap(True); root.addWidget(self.summary)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.button(QDialogButtonBox.StandardButton.Ok).setText("Zakończ / zapisz konfigurację próbki")
        buttons.button(QDialogButtonBox.StandardButton.Cancel).setText("Anuluj")
        buttons.accepted.connect(self._accept); buttons.rejected.connect(self.reject)
        root.addWidget(buttons)

        self.periodic_btn.clicked.connect(self._periodic)
        self.add_btn.clicked.connect(lambda: self._append_row(""))
        self.del_btn.clicked.connect(self._remove_selected)
        self.up_btn.clicked.connect(lambda: self._move(-1))
        self.down_btn.clicked.connect(lambda: self._move(+1))
        self.suggest_btn.clicked.connect(self._suggest)
        self.apply_preset_btn.clicked.connect(self._apply_material_preset)
        self.mode_combo.currentIndexChanged.connect(self._mode_changed)
        self.comp_table.itemChanged.connect(lambda *_: self._update_summary())
        self.name_edit.textChanged.connect(lambda *_: self._update_summary())
        self.material_edit.textChanged.connect(lambda *_: self._update_summary())
        self.base_a_edit.textChanged.connect(lambda *_: self._update_summary())
        self.base_b_edit.textChanged.connect(lambda *_: self._update_summary())
        self.filler_edit.textChanged.connect(lambda *_: self._update_summary())
        self.matrix_combo.currentIndexChanged.connect(lambda *_: self._update_summary())

        self._load_existing()
        self._mode_changed()
        self._update_summary()

    def _append_row(self, element: str, value: str = "", vmin: str = "", vmax: str = "", notes: str = ""):
        r = self.comp_table.rowCount(); self.comp_table.insertRow(r)
        for c, v in enumerate((element, value, vmin, vmax, notes)):
            self.comp_table.setItem(r, c, QTableWidgetItem(str(v)))

    def _load_existing(self):
        rows = self.dataset.metadata.get("chemical_composition", [])
        if isinstance(rows, list) and rows:
            for rec in rows:
                if not isinstance(rec, dict): continue
                self._append_row(rec.get("element", ""), rec.get("value", ""), rec.get("min", ""), rec.get("max", ""), rec.get("notes", ""))
        else:
            raw = str(self.dataset.metadata.get("allowed_elements", "") or "")
            for el in re.split(r"[,;\s]+", raw.strip()):
                if el:
                    self._append_row(el)
        saved_kind = str(self.dataset.metadata.get("sample_kind", ""))
        if saved_kind and self.kind_combo.findText(saved_kind) >= 0:
            self.kind_combo.setCurrentText(saved_kind)

    def _rows(self):
        out=[]
        for r in range(self.comp_table.rowCount()):
            vals=[]
            for c in range(5):
                item=self.comp_table.item(r,c); vals.append(item.text().strip() if item else "")
            el=vals[0]
            if not el: continue
            el=el[0].upper()+el[1:].lower() if len(el)>1 else el.upper()
            out.append({"element":el,"value":vals[1],"min":vals[2],"max":vals[3],"notes":vals[4],"rank":len(out)+1})
        return out

    def _periodic(self):
        current={r["element"] for r in self._rows() if r.get("element")}
        dlg=PeriodicTableDialog(current,self)
        if dlg.exec()==QDialog.DialogCode.Accepted:
            chosen=dlg.selected_elements(); existing=[r["element"] for r in self._rows()]
            # Preserve existing row order; append newly selected elements by atomic number.
            for el in sorted(chosen-set(existing), key=lambda e:ELEMENT_Z.get(e,999)):
                self._append_row(el)
            # Remove elements that were explicitly unselected.
            for r in range(self.comp_table.rowCount()-1,-1,-1):
                item=self.comp_table.item(r,0); el=item.text().strip() if item else ""
                el=el[0].upper()+el[1:].lower() if len(el)>1 else el.upper()
                if el and el not in chosen: self.comp_table.removeRow(r)
            self._update_summary()

    def _remove_selected(self):
        rows=sorted({i.row() for i in self.comp_table.selectionModel().selectedRows()}, reverse=True)
        if not rows and self.comp_table.currentRow()>=0: rows=[self.comp_table.currentRow()]
        for r in rows: self.comp_table.removeRow(r)
        self._update_summary()

    def _move(self, direction:int):
        r=self.comp_table.currentRow(); nr=r+direction
        if r<0 or nr<0 or nr>=self.comp_table.rowCount(): return
        vals=[]
        for c in range(5):
            a=self.comp_table.takeItem(r,c); b=self.comp_table.takeItem(nr,c); vals.append((a,b))
        for c,(a,b) in enumerate(vals):
            self.comp_table.setItem(r,c,b or QTableWidgetItem("")); self.comp_table.setItem(nr,c,a or QTableWidgetItem(""))
        self.comp_table.selectRow(nr); self._update_summary()

    def _apply_material_preset(self):
        key = str(self.preset_combo.currentData() or "")
        if not key or key not in MATERIAL_SETUP_PRESETS:
            return
        rec = MATERIAL_SETUP_PRESETS[key]
        self._preset_was_applied = key
        self.material_edit.setText(str(rec.get("material", key)))
        # Presets intentionally encode order, not invented exact percentages.
        for i in range(self.mode_combo.count()):
            if self.mode_combo.itemData(i) == "order":
                self.mode_combo.setCurrentIndex(i); break
        self.comp_table.setRowCount(0)
        for el in rec.get("elements", []):
            self._append_row(str(el))
        matrix = str(rec.get("matrix", "auto"))
        for i in range(self.matrix_combo.count()):
            if self.matrix_combo.itemData(i) == matrix:
                self.matrix_combo.setCurrentIndex(i); break
        self._update_summary()

    def _suggest(self):
        filename_hint = str((self.dataset.metadata or {}).get("filename_material_suggestion", "") or "")
        text=" ".join([self.material_edit.text(), self.name_edit.text(), self.base_a_edit.text(), self.base_b_edit.text(), self.filler_edit.text(), filename_hint])
        suggestion=default_elements_for_material(text,text)
        els=[e for e in re.split(r"[,;\s]+", suggestion) if e]
        if not els:
            QMessageBox.information(self,"Konfiguracja próbki","Brak niezweryfikowanej sugestii na podstawie tekstu/nazwy. Wprowadź znaną chemię ręcznie albo samodzielnie wybierz wbudowany preset.")
            return
        existing={r["element"] for r in self._rows()}
        for el in els:
            if el not in existing:
                self._append_row(el); existing.add(el)
        if self.mode_combo.currentData()=="unknown":
            for i in range(self.mode_combo.count()):
                if self.mode_combo.itemData(i)=="elements": self.mode_combo.setCurrentIndex(i); break
        self._update_summary()

    def _mode_changed(self):
        code=str(self.mode_combo.currentData())
        unknown=code=="unknown"
        self.comp_table.setEnabled(not unknown)
        for b in (self.periodic_btn,self.add_btn,self.del_btn,self.up_btn,self.down_btn,self.suggest_btn): b.setEnabled(not unknown)
        numeric = code in ("exact","ranges")
        self.unit_combo.setEnabled(numeric)
        self._update_summary()

    def _update_summary(self):
        code=str(self.mode_combo.currentData())
        rows=self._rows()
        if code=="unknown":
            txt="Filtr chemiczny nie zostanie zastosowany. Detekcja pików pozostaje dostępna; identyfikacja faz wyłącznie z zewnętrznej bazy COD będzie pomijana do czasu podania składu chemicznego."
        elif code=="order":
            txt="Znana kolejność: " + (" > ".join(r["element"] for r in rows) if rows else "(brak)")
        else:
            txt="Dozwolone pierwiastki: " + (", ".join(r["element"] for r in rows) if rows else "(brak)")
        matrix_label = self.matrix_combo.currentText()
        txt += f"\nWskazówka dotycząca oczekiwanej osnowy: {matrix_label}. Zgodność chemii i pików nadal pozostaje częścią rankingu."
        self.summary.setText(txt)

    def _accept(self):
        code=str(self.mode_combo.currentData())
        rows=self._rows()
        invalid=[r["element"] for r in rows if r["element"] not in ELEMENT_Z]
        if invalid:
            QMessageBox.warning(self,"Konfiguracja próbki","Nieznany symbol / symbole pierwiastków: "+", ".join(invalid)); return
        if code!="unknown" and not rows:
            QMessageBox.warning(self,"Konfiguracja próbki","Dodaj co najmniej jeden znany pierwiastek albo wybierz „Skład nieznany”. "); return
        ds=self.dataset
        ds.name=self.name_edit.text().strip() or ds.name
        ds.material_type=self.material_edit.text().strip()
        ds.metadata_sources["material_type"]="user / sample setup wizard"
        ds.metadata["sample_kind"]=self.kind_combo.currentText()
        ds.metadata["base_material_A"]=self.base_a_edit.text().strip()
        ds.metadata["base_material_B"]=self.base_b_edit.text().strip()
        ds.metadata["filler_material"]=self.filler_edit.text().strip()
        ds.metadata_sources["base_material_A"]="user / sample setup wizard"
        ds.metadata_sources["base_material_B"]="user / sample setup wizard"
        ds.metadata_sources["filler_material"]="user / sample setup wizard"
        ds.metadata["expected_matrix_structure"]=str(self.matrix_combo.currentData() or "auto")
        ds.metadata_sources["expected_matrix_structure"]="user / sample setup wizard"
        if self._preset_was_applied:
            ds.metadata["material_preset_confirmed"] = self._preset_was_applied
            ds.metadata_sources["material_preset_confirmed"] = "user / sample setup wizard"
        else:
            ds.metadata.pop("material_preset_confirmed", None)
            ds.metadata_sources.pop("material_preset_confirmed", None)
        ds.metadata["composition_mode"]=code
        ds.metadata["composition_unit"]=self.unit_combo.currentText()
        ds.metadata["chemical_composition"]=rows if code!="unknown" else []
        ds.metadata_sources["chemical_composition"]="user / sample setup wizard"
        if code=="unknown":
            ds.metadata.pop("allowed_elements",None); ds.metadata_sources.pop("allowed_elements",None)
        else:
            # Preserve the user's order; it has meaning in order mode.
            ds.metadata["allowed_elements"]=" ".join(r["element"] for r in rows)
            ds.metadata_sources["allowed_elements"]="user / sample setup wizard"
        self.accept()

# --------------------------- measured-pattern COD search ---------------------------

class BatchChemistryDialog(QDialog):
    """Ręczna konfiguracja ograniczeń chemicznych dla serii próbek.

    Wersja v20.1 celowo NIE wiąże na żywo pola „pierwiastki wspólne” z każdym
    rekordem. Pole wspólne jest jedynie buforem pomocniczym: użytkownik może
    jawnie dodać je do zaznaczonych próbek lub zastąpić nim ich skład. Każda
    próbka ma własny, edytowalny zestaw pierwiastków i może zostać niezależnie
    wyczyszczona. Dzięki temu zmiana jednego rekordu nie zmienia pozostałych.
    """
    def __init__(self, datasets: List[Dataset], parent=None):
        super().__init__(parent)
        self.datasets = datasets
        self.setWindowTitle('Chemia serii / ograniczenia składu')
        self.resize(1280, 760)
        root = QVBoxLayout(self)

        info = QLabel(
            'Każda próbka ma własny zestaw pierwiastków. Pole „pierwiastki wspólne” nie zmienia próbek automatycznie — '
            'najpierw wpisz/wybierz pierwiastki wspólne, a następnie użyj przycisku „Dodaj wspólne do zaznaczonych” albo '
            '„Zastąp zaznaczone wspólnymi”. Możesz też edytować lub wyczyścić każdy rekord osobno. '
            'Nazwy plików nigdy nie są traktowane jako źródło składu chemicznego.'
        )
        info.setWordWrap(True)
        root.addWidget(info)

        common_box = QGroupBox('Pierwiastki wspólne — bufor pomocniczy')
        ch = QGridLayout(common_box)
        self.common_edit = QLineEdit()
        self.common_edit.setPlaceholderText('np. Fe Cr Ni')
        self.common_pt = QPushButton('Układ okresowy…')
        self.intersection_btn = QPushButton('Wczytaj część wspólną istniejącej chemii')
        self.add_common_btn = QPushButton('Dodaj wspólne do zaznaczonych ✓')
        self.replace_common_btn = QPushButton('Zastąp zaznaczone wspólnymi ✓')
        ch.addWidget(self.common_edit, 0, 0, 1, 3)
        ch.addWidget(self.common_pt, 0, 3)
        ch.addWidget(self.intersection_btn, 1, 0)
        ch.addWidget(self.add_common_btn, 1, 1)
        ch.addWidget(self.replace_common_btn, 1, 2)
        root.addWidget(common_box)

        self.table = QTableWidget(len(datasets), 6)
        self.table.setHorizontalHeaderLabels([
            'Użyj', 'Próbka', 'Pierwiastki próbki', 'Dodatkowe pierwiastki',
            'Efektywne pierwiastki', 'Oczekiwana osnowa / struktura'
        ])
        hh = self.table.horizontalHeader()
        hh.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        hh.setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        hh.setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        hh.setSectionResizeMode(3, QHeaderView.ResizeMode.Stretch)
        hh.setSectionResizeMode(4, QHeaderView.ResizeMode.Stretch)
        hh.setSectionResizeMode(5, QHeaderView.ResizeMode.Stretch)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.verticalHeader().setDefaultSectionSize(34)

        for r, ds in enumerate(datasets):
            use = QTableWidgetItem()
            use.setFlags(use.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            use.setCheckState(Qt.CheckState.Checked)
            self.table.setItem(r, 0, use)

            nm = QTableWidgetItem(ds.name)
            nm.setFlags(nm.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.table.setItem(r, 1, nm)

            allowed = set(dataset_allowed_elements(ds))
            old_extra = self._parse(ds.metadata.get('batch_extra_elements', ''))
            base = allowed - old_extra if old_extra else allowed
            self.table.setItem(r, 2, QTableWidgetItem(self._fmt(base)))
            self.table.setItem(r, 3, QTableWidgetItem(self._fmt(old_extra)))

            eff = QTableWidgetItem('')
            eff.setFlags(eff.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.table.setItem(r, 4, eff)

            combo = QComboBox()
            combo.addItems([label for label, _code in EXPECTED_MATRIX_ITEMS])
            code = str(ds.metadata.get('expected_matrix_structure', 'unknown') or 'unknown')
            label = next((lab for lab, c in EXPECTED_MATRIX_ITEMS if c == code),
                         'Nieznana / bez preferencji struktury')
            combo.setCurrentText(label)
            self.table.setCellWidget(r, 5, combo)

        root.addWidget(self.table, 1)

        # Narzędzia dla pojedynczego rekordu.
        row_one = QHBoxLayout()
        self.sample_pt = QPushButton('Układ okresowy — pierwiastki wybranej próbki…')
        self.extra_pt = QPushButton('Układ okresowy — dodatkowe dla wybranej próbki…')
        self.clear_current_base = QPushButton('Wyczyść pierwiastki wybranej próbki')
        self.clear_current_extra = QPushButton('Wyczyść dodatkowe wybranej próbki')
        for b in (self.sample_pt, self.extra_pt, self.clear_current_base, self.clear_current_extra):
            row_one.addWidget(b)
        root.addLayout(row_one)

        # Narzędzia dla grupy / wszystkich rekordów.
        row_many = QHBoxLayout()
        self.select_all = QPushButton('Zaznacz wszystko ✓')
        self.clear_sel = QPushButton('Odznacz wszystko')
        self.clear_checked = QPushButton('Wyczyść chemię zaznaczonych ✓')
        self.clear_all = QPushButton('Wyczyść chemię wszystkich próbek')
        self.preview_btn = QPushButton('Odśwież efektywne pierwiastki')
        for b in (self.select_all, self.clear_sel, self.clear_checked, self.clear_all, self.preview_btn):
            row_many.addWidget(b)
        root.addLayout(row_many)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.button(QDialogButtonBox.StandardButton.Ok).setText('Zastosuj do zaznaczonych próbek ✓')
        buttons.button(QDialogButtonBox.StandardButton.Cancel).setText('Anuluj')
        root.addWidget(buttons)

        buttons.accepted.connect(self.apply_and_accept)
        buttons.rejected.connect(self.reject)
        self.select_all.clicked.connect(lambda: self._set_checks(True))
        self.clear_sel.clicked.connect(lambda: self._set_checks(False))
        self.common_pt.clicked.connect(self._common_periodic)
        self.sample_pt.clicked.connect(lambda: self._row_periodic(2))
        self.extra_pt.clicked.connect(lambda: self._row_periodic(3))
        self.clear_current_base.clicked.connect(lambda: self._clear_current_cell(2))
        self.clear_current_extra.clicked.connect(lambda: self._clear_current_cell(3))
        self.clear_checked.clicked.connect(self._clear_checked_rows)
        self.clear_all.clicked.connect(self._clear_all_rows)
        self.preview_btn.clicked.connect(self.refresh_effective)
        self.intersection_btn.clicked.connect(self._load_intersection)
        self.add_common_btn.clicked.connect(lambda: self._apply_common_to_checked(replace=False))
        self.replace_common_btn.clicked.connect(lambda: self._apply_common_to_checked(replace=True))
        self.table.itemChanged.connect(self._table_item_changed)

        # Uzupełnij pole wspólne tylko informacyjnie. Nie modyfikuje ono żadnej próbki.
        self._load_intersection(initial_only=True)
        self.refresh_effective()

    @staticmethod
    def _parse(text):
        return CODPatternSearchDialog._parse_elements(str(text or ''))

    @staticmethod
    def _fmt(elements):
        return ' '.join(sorted(set(elements or []), key=lambda e: ELEMENT_Z.get(e, 999)))

    def _checked(self, r):
        it = self.table.item(r, 0)
        return bool(it and it.checkState() == Qt.CheckState.Checked)

    def _checked_rows(self):
        return [r for r in range(self.table.rowCount()) if self._checked(r)]

    def _mark_row_for_apply(self, r):
        if 0 <= r < self.table.rowCount():
            it = self.table.item(r, 0)
            if it is not None and it.checkState() != Qt.CheckState.Checked:
                self.table.blockSignals(True)
                try:
                    it.setCheckState(Qt.CheckState.Checked)
                finally:
                    self.table.blockSignals(False)

    def _set_checks(self, val):
        self.table.blockSignals(True)
        try:
            for r in range(self.table.rowCount()):
                self.table.item(r, 0).setCheckState(Qt.CheckState.Checked if val else Qt.CheckState.Unchecked)
        finally:
            self.table.blockSignals(False)
        self.refresh_effective()

    def _table_item_changed(self, item):
        if item is None:
            return
        if item.column() in (2, 3):
            self._mark_row_for_apply(item.row())
        if item.column() in (0, 2, 3):
            self.refresh_effective()

    def _load_intersection(self, initial_only=False):
        sets = [set(dataset_allowed_elements(ds)) for ds in self.datasets if dataset_allowed_elements(ds)]
        inter = set.intersection(*sets) if sets else set()
        if inter or not initial_only:
            self.common_edit.setText(self._fmt(inter))

    def _common_periodic(self):
        dlg = PeriodicTableDialog(self._parse(self.common_edit.text()), self)
        if dlg.exec() == QDialog.DialogCode.Accepted:
            self.common_edit.setText(self._fmt(dlg.selected_elements()))

    def _current_row_or_warn(self):
        r = self.table.currentRow()
        if r < 0:
            QMessageBox.information(self, 'Chemia serii', 'Najpierw zaznacz wiersz próbki.')
            return -1
        return r

    def _row_periodic(self, column):
        r = self._current_row_or_warn()
        if r < 0:
            return
        cur = self._parse(self.table.item(r, column).text() if self.table.item(r, column) else '')
        dlg = PeriodicTableDialog(cur, self)
        if dlg.exec() == QDialog.DialogCode.Accepted:
            self._mark_row_for_apply(r)
            self.table.item(r, column).setText(self._fmt(dlg.selected_elements()))

    def _clear_current_cell(self, column):
        r = self._current_row_or_warn()
        if r >= 0 and self.table.item(r, column):
            self._mark_row_for_apply(r)
            self.table.item(r, column).setText('')

    def _clear_checked_rows(self):
        rows = self._checked_rows()
        if not rows:
            QMessageBox.information(self, 'Chemia serii', 'Nie zaznaczono żadnych próbek symbolem ✓.')
            return
        self.table.blockSignals(True)
        try:
            for r in rows:
                self.table.item(r, 2).setText('')
                self.table.item(r, 3).setText('')
        finally:
            self.table.blockSignals(False)
        self.refresh_effective()

    def _clear_all_rows(self):
        ans = QMessageBox.question(
            self, 'Chemia serii',
            'Wyczyścić pierwiastki we wszystkich widocznych próbkach?\n\nZmiana zostanie zapisana dopiero po naciśnięciu „Zastosuj”.',
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        if ans != QMessageBox.StandardButton.Yes:
            return
        self.table.blockSignals(True)
        try:
            for r in range(self.table.rowCount()):
                self.table.item(r, 0).setCheckState(Qt.CheckState.Checked)
                self.table.item(r, 2).setText('')
                self.table.item(r, 3).setText('')
        finally:
            self.table.blockSignals(False)
        self.refresh_effective()

    def _apply_common_to_checked(self, replace=False):
        common = self._parse(self.common_edit.text())
        rows = self._checked_rows()
        if not rows:
            QMessageBox.information(self, 'Chemia serii', 'Nie zaznaczono żadnych próbek symbolem ✓.')
            return
        if not common and not replace:
            QMessageBox.information(self, 'Chemia serii', 'Pole pierwiastków wspólnych jest puste.')
            return
        self.table.blockSignals(True)
        try:
            for r in rows:
                current = self._parse(self.table.item(r, 2).text() if self.table.item(r, 2) else '')
                new = set(common) if replace else (current | common)
                self.table.item(r, 2).setText(self._fmt(new))
        finally:
            self.table.blockSignals(False)
        self.refresh_effective()

    def refresh_effective(self):
        # Każdy rekord jest niezależny: NIE dołączamy tutaj pola common_edit.
        self.table.blockSignals(True)
        try:
            for r in range(self.table.rowCount()):
                base = self._parse(self.table.item(r, 2).text() if self.table.item(r, 2) else '')
                extra = self._parse(self.table.item(r, 3).text() if self.table.item(r, 3) else '')
                eff = base | extra
                text = self._fmt(eff)
                if self.table.item(r, 4) is None:
                    it = QTableWidgetItem(text)
                    it.setFlags(it.flags() & ~Qt.ItemFlag.ItemIsEditable)
                    self.table.setItem(r, 4, it)
                else:
                    self.table.item(r, 4).setText(text)
        finally:
            self.table.blockSignals(False)

    def apply_and_accept(self):
        rows = self._checked_rows()
        if not rows:
            QMessageBox.warning(self, 'Chemia serii', 'Zaznacz co najmniej jedną próbkę symbolem ✓.')
            return

        common_buffer = self._parse(self.common_edit.text())
        for r in rows:
            ds = self.datasets[r]
            base = self._parse(self.table.item(r, 2).text() if self.table.item(r, 2) else '')
            extra = self._parse(self.table.item(r, 3).text() if self.table.item(r, 3) else '')
            eff = base | extra

            if eff:
                txt = self._fmt(eff)
                ds.metadata['allowed_elements'] = txt
                ds.metadata_sources['allowed_elements'] = 'użytkownik / chemia serii'
                ds.metadata['composition_mode'] = 'elements'
                ds.metadata['chemical_composition'] = [
                    {'element': e, 'value': ''} for e in sorted(eff, key=lambda e: ELEMENT_Z.get(e, 999))
                ]
                ds.metadata_sources['chemical_composition'] = 'użytkownik / chemia serii (lista pierwiastków)'
            else:
                ds.metadata.pop('allowed_elements', None)
                ds.metadata_sources.pop('allowed_elements', None)
                ds.metadata['composition_mode'] = 'unknown'
                ds.metadata['chemical_composition'] = []
                ds.metadata_sources['chemical_composition'] = 'użytkownik / chemia serii — wyczyszczono'

            # Pola pomocnicze są zapisywane wyłącznie do odtworzenia układu dialogu.
            ds.metadata['batch_common_elements'] = self._fmt(common_buffer)
            ds.metadata['batch_sample_elements'] = self._fmt(base)
            ds.metadata['batch_extra_elements'] = self._fmt(extra)

            combo = self.table.cellWidget(r, 5)
            label = combo.currentText() if isinstance(combo, QComboBox) else 'Nieznana / bez preferencji struktury'
            code = next((c for lab, c in EXPECTED_MATRIX_ITEMS if lab == label), 'unknown')
            ds.metadata['expected_matrix_structure'] = code
            ds.metadata_sources['expected_matrix_structure'] = 'użytkownik / chemia serii'

        self.accept()

class CODPatternSearchDialog(QDialog):
    def __init__(self, dataset: Dataset, cod_db: CODDatabase, phase_db: PhaseDatabase,
                 tolerance_deg: float = 0.35, parent=None):
        super().__init__(parent)
        self.dataset = dataset
        self.cod_db = cod_db
        self.phase_db = phase_db
        self.parent_window = parent
        self.default_tolerance_deg = float(tolerance_deg)
        self.results = []
        self.used_obs = []
        self.setWindowTitle("Wyszukaj COD dla zmierzonego dyfraktogramu")
        self.resize(1320, 780)
        root = QVBoxLayout(self)

        info = QLabel(
            f"Zmierzona próbka: {dataset.name} | wykryte piki: {len(dataset.peaks)} | "
            f"λ = {dataset.wavelength_A:.6f} Å\n"
            "Wyszukiwanie porównuje zmierzone położenia pików z refleksami wygenerowanymi z komórek elementarnych i grup przestrzennych COD. "
            "Jest to narzędzie do przesiewowej oceny kandydatów: większość rekordów COD .db3 nie zawiera gotowego eksperymentalnego profilu intensywności proszkowej. "
            "Kα₂ jest oznaczane heurystycznie tylko przy dwóch źródłowych długościach fal; nie usuwa to zmierzonych pików."
        )
        info.setWordWrap(True)
        root.addWidget(info)

        # Show the exact Python environment used by the GUI. This makes Windows
        # multi-Python installations much easier to diagnose.
        if ensure_gemmi():
            gver = getattr(gemmi, "__version__", "unknown")
            self.gemmi_status = QLabel(f"gemmi: {gver}  |  Python: {sys.executable}")
        else:
            self.gemmi_status = QLabel(f"gemmi: NIE WYKRYTO  |  Python: {sys.executable}")
            self.gemmi_status.setToolTip(GEMMI_IMPORT_ERROR or "nie udało się zaimportować gemmi")
        self.gemmi_status.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        root.addWidget(self.gemmi_status)

        settings_box = QGroupBox("Ustawienia wyszukiwania")
        sf = QFormLayout(settings_box)
        saved_elements = str(dataset.metadata.get("allowed_elements", "") or "").strip()
        self.elements_edit = QLineEdit(saved_elements)
        filename_hint = str(dataset.metadata.get("filename_material_suggestion", "") or "").strip()
        manual_hint = default_elements_for_material(dataset.material_type, dataset.name) if dataset.material_type else ""
        hint = manual_hint or (default_elements_for_material(filename_hint, "") if filename_hint else "")
        self.elements_edit.setPlaceholderText(("Niezweryfikowana sugestia (niezastosowana): " + hint) if hint else "Dozwolone pierwiastki, np. Fe Cr Ni Mo Cu Mn Si C")
        self.elements_edit.setToolTip(
            "Uwzględniane są tylko wzory COD złożone z tych pierwiastków. Dodaj O, N, C itd., jeśli spodziewasz się takich faz."
        )
        self.periodic_btn = QPushButton("Układ okresowy…")
        self.periodic_btn.setToolTip("Wybierz dozwolone pierwiastki na interaktywnym układzie okresowym")
        self.periodic_btn.clicked.connect(self.open_periodic_table)
        elements_widget = QWidget()
        elements_layout = QHBoxLayout(elements_widget)
        elements_layout.setContentsMargins(0, 0, 0, 0)
        elements_layout.setSpacing(6)
        elements_layout.addWidget(self.elements_edit, 1)
        elements_layout.addWidget(self.periodic_btn)
        self.tolerance = QDoubleSpinBox(); self.tolerance.setRange(0.02, 2.0); self.tolerance.setDecimals(3); self.tolerance.setValue(float(tolerance_deg)); self.tolerance.setSuffix(" °")
        self.max_peaks = QSpinBox(); self.max_peaks.setRange(2, 30); self.max_peaks.setValue(min(8, max(2, len(dataset.peaks))))
        self.min_rel = QDoubleSpinBox(); self.min_rel.setRange(0.0, 100.0); self.min_rel.setValue(5.0); self.min_rel.setSuffix(" %")
        self.max_candidates = QSpinBox(); self.max_candidates.setRange(100, 50000); self.max_candidates.setSingleStep(500); self.max_candidates.setValue(5000)
        self.top_results = QSpinBox(); self.top_results.setRange(5, 200); self.top_results.setValue(30)
        sf.addRow("Dozwolone pierwiastki:", elements_widget)
        sf.addRow("Tolerancja dopasowania 2θ:", self.tolerance)
        sf.addRow("Liczba użytych najsilniejszych pików:", self.max_peaks)
        sf.addRow("Minimalna względna intensywność:", self.min_rel)
        self.kalpha2_filter_cb = QCheckBox("Nie traktuj rozdzielonych kandydatów Kα₂ jako niezależnych refleksów (tylko ranking COD)")
        self.kalpha2_filter_cb.setToolTip("Działa tylko przy dwóch długościach fal odczytanych z pliku. Nie usuwa pików; umożliwia wyłączenie filtra.")
        self.kalpha2_filter_cb.setChecked(verified_wavelength_pair(dataset) is not None)
        sf.addRow("Dublet Kα:", self.kalpha2_filter_cb)
        sf.addRow("Maks. liczba kandydatów COD po filtrze chemicznym:", self.max_candidates)
        sf.addRow("Liczba pokazywanych najlepszych wyników:", self.top_results)
        self.reset_search_defaults_btn = QPushButton("Przywróć domyślne")
        self.reset_search_defaults_btn.setToolTip("Przywróć zalecane ustawienia wyszukiwania wzorca COD dla tej próbki")
        self.reset_search_defaults_btn.clicked.connect(self.reset_search_defaults)
        sf.addRow(self.reset_search_defaults_btn)
        root.addWidget(settings_box)

        search_row = QHBoxLayout()
        self.search_btn = QPushButton("Wyszukaj COD dla zmierzonego dyfraktogramu")
        self.search_btn.clicked.connect(self.run_search)
        self.status = QLabel("Gotowy.")
        self.status.setWordWrap(True)
        search_row.addWidget(self.search_btn)
        search_row.addWidget(self.status, 1)
        root.addLayout(search_row)

        self.table = QTableWidget(0, 15)
        self.table.setHorizontalHeaderLabels([
            "Rank", "Score [%]", "COD ID", "Formula", "Czytelna nazwa fazy", "Space group",
            "Dopasowane piki", "RMS Δ2θ [°]", "Generated refs", "a [Å]", "b [Å]", "c [Å]", "Year",
            "Mnożnik chemii [×]", "Mnożnik struktury [×]"
        ])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(4, QHeaderView.ResizeMode.Stretch)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setSortingEnabled(False)
        self.table.doubleClicked.connect(self.show_match_details)
        self.table.itemSelectionChanged.connect(self.populate_peak_choices)
        root.addWidget(self.table, 1)

        peak_box = QGroupBox("Proponowane przypisania pików dla wybranego kandydata COD")
        peak_lay = QVBoxLayout(peak_box)
        self.peak_choice_info = QLabel(
            "Wybierz powyżej kandydata COD, a następnie zaznacz tylko te zmierzone piki, które chcesz zaakceptować. "
            "Jest to ręczne przypisanie; rekord COD nie jest przez to importowany do bazy wewnętrznej."
        )
        self.peak_choice_info.setWordWrap(True)
        peak_lay.addWidget(self.peak_choice_info)
        self.peak_choice_table = QTableWidget(0, 9)
        self.peak_choice_table.setHorizontalHeaderLabels([
            "Użyj", "Zmierzony 2θ [°]", "Irel [%]", "COD 2θ [°]", "hkl", "Δ2θ [°]",
            "W tolerancji", "Bieżąca faza", "Bieżące hkl"
        ])
        self.peak_choice_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.peak_choice_table.horizontalHeader().setSectionResizeMode(7, QHeaderView.ResizeMode.Stretch)
        self.peak_choice_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.peak_choice_table.setMinimumHeight(170)
        self.peak_choice_table.setMaximumHeight(260)
        peak_lay.addWidget(self.peak_choice_table)
        peak_buttons = QHBoxLayout()
        self.select_matched_btn = QPushButton("Zaznacz dopasowane piki")
        self.clear_peak_selection_btn = QPushButton("Wyczyść wybór pików")
        self.apply_selected_peaks_btn = QPushButton("Zastosuj wybrane piki COD do próbki")
        self.clear_manual_cod_btn = QPushButton("Usuń ręczne przypisania pików COD")
        self.select_matched_btn.clicked.connect(self.select_all_matched_peak_choices)
        self.clear_peak_selection_btn.clicked.connect(self.clear_peak_choices)
        self.apply_selected_peaks_btn.clicked.connect(self.apply_selected_peak_assignments)
        self.clear_manual_cod_btn.clicked.connect(self.clear_manual_cod_assignments)
        for b in (self.select_matched_btn, self.clear_peak_selection_btn, self.apply_selected_peaks_btn, self.clear_manual_cod_btn):
            peak_buttons.addWidget(b)
        peak_lay.addLayout(peak_buttons)
        root.addWidget(peak_box)

        note = QLabel(
            "Interpretacja: wysoki wynik oznacza, że wybrane zmierzone piki leżą blisko refleksów dozwolonych przez strukturę z COD. "
            "Ranking łączy zgodność pików, priorytet wynikający ze składu chemicznego oraz opcjonalnie oczekiwaną osnowę/strukturę z Kreatora konfiguracji próbki. "
            "Rekordy opisujące wyłącznie pierwiastek minorowy (np. sam Cr dla próbki 316L z dominującym Fe) są silnie obniżane w rankingu. "
            "Dla stopów substytucyjnych wybrany rekord COD może być wzorcem strukturalnym, a nie dokładnym składem stopu; krytyczne identyfikacje należy potwierdzić składem (np. EDS) i w razie potrzeby analizą całego profilu/Rietvelda."
        )
        note.setWordWrap(True)
        root.addWidget(note)

        buttons = QHBoxLayout()
        details = QPushButton("Pokaż dopasowane piki")
        cod_details_btn = QPushButton("Pokaż szczegóły rekordu COD")
        import_btn = QPushButton("Importuj wybranego kandydata COD")
        export_btn = QPushButton("Eksportuj wyniki wyszukiwania")
        close_btn = QPushButton("Zamknij")
        details.clicked.connect(self.show_match_details)
        cod_details_btn.clicked.connect(self.show_cod_record_details)
        import_btn.clicked.connect(self.import_selected)
        export_btn.clicked.connect(self.export_results)
        close_btn.clicked.connect(self.accept)
        for b in (details, cod_details_btn, import_btn, export_btn): buttons.addWidget(b)
        buttons.addStretch(1); buttons.addWidget(close_btn)
        root.addLayout(buttons)

    def reset_search_defaults(self):
        """Restore recommended COD-search defaults without touching measured data."""
        saved_elements = str(self.dataset.metadata.get("allowed_elements", "") or "").strip()
        self.elements_edit.setText(saved_elements)
        self.tolerance.setValue(self.default_tolerance_deg)
        self.max_peaks.setValue(min(8, max(2, len(self.dataset.peaks))))
        self.min_rel.setValue(5.0)
        self.max_candidates.setValue(5000)
        self.top_results.setValue(30)
        self.kalpha2_filter_cb.setChecked(verified_wavelength_pair(self.dataset) is not None)
        self.status.setText("Przywrócono zalecane ustawienia wyszukiwania COD.")

    @staticmethod
    def _parse_elements(text: str) -> set[str]:
        tokens = re.split(r"[,;\s]+", text.strip())
        out = set()
        for t in tokens:
            if not t:
                continue
            t = t[0].upper() + t[1:].lower()
            if re.fullmatch(r"[A-Z][a-z]?", t):
                out.add(t)
        return out

    def open_periodic_table(self):
        current = self._parse_elements(self.elements_edit.text())
        dlg = PeriodicTableDialog(current, self)
        if dlg.exec() == QDialog.DialogCode.Accepted:
            selected = dlg.selected_elements()
            ordered = sorted(selected, key=lambda s: ELEMENT_Z.get(s, 999))
            value = " ".join(ordered)
            self.elements_edit.setText(value)
            self.dataset.metadata["allowed_elements"] = value
            self.dataset.metadata_sources["allowed_elements"] = "user / periodic table"

    def run_search(self):
        if not self.cod_db.is_open:
            QMessageBox.information(self, "Wyszukiwanie wzorca COD", "Najpierw wczytaj bazę COD .db3.")
            return
        if not ensure_gemmi():
            QMessageBox.critical(self, "Wyszukiwanie wzorca COD", gemmi_diagnostic_message("Wyszukiwanie wzorca COD"))
            return
        elements = self._parse_elements(self.elements_edit.text())
        if elements:
            ordered = sorted(elements, key=lambda e: ELEMENT_Z.get(e, 999))
            self.dataset.metadata["allowed_elements"] = " ".join(ordered)
            self.dataset.metadata_sources["allowed_elements"] = "user / COD search"
        if not elements:
            QMessageBox.information(
                self, "Wyszukiwanie wzorca COD",
                "Wprowadź dozwolone pierwiastki chemiczne. Dla spoin różnoimiennych można np. podać: Fe Cr Ni Mo Cu Mn Si C"
            )
            return

        progress = QProgressDialog("Filtrowanie COD i porównywanie wygenerowanych refleksów…", "Anuluj", 0, self.max_candidates.value(), self)
        progress.setWindowTitle("Wyszukaj COD dla zmierzonego dyfraktogramu")
        progress.setWindowModality(Qt.WindowModality.WindowModal)
        progress.setMinimumDuration(0)
        progress.show()

        def progress_cb(done, total):
            progress.setMaximum(max(1, total))
            progress.setValue(min(done, total))
            progress.setLabelText(f"Porównywanie kandydatów COD: {done:,} / {total:,}")
            QApplication.processEvents()

        try:
            self.search_btn.setEnabled(False)
            results, obs, candidate_count = self.cod_db.search_measured_pattern(
                self.dataset,
                elements,
                tolerance_deg=self.tolerance.value(),
                max_candidates=self.max_candidates.value(),
                top_results=self.top_results.value(),
                max_observed_peaks=self.max_peaks.value(),
                min_relative_intensity=self.min_rel.value(),
                exclude_resolved_kalpha2=self.kalpha2_filter_cb.isChecked(),
                progress_callback=progress_cb,
                cancel_callback=progress.wasCanceled,
            )
            self.results = results
            self.used_obs = obs
            self.populate_results()
            used_txt = ", ".join(f"{p.two_theta:.3f}°" for p in obs)
            suppressed = len(dataset_kalpha2_candidates(self.dataset)) if self.kalpha2_filter_cb.isChecked() else 0
            self.status.setText(
                f"Kandydaci COD po filtrze chemicznym: {candidate_count:,}. "
                f"Użyte zmierzone piki ({len(obs)}): {used_txt}. Wyniki: {len(results)}. "
                f"Pominięte kandydaty Kα₂: {suppressed} (piki zachowane w projekcie)."
            )
        except Exception as e:
            QMessageBox.critical(self, "Wyszukiwanie wzorca COD", f"{e}\n\n{traceback.format_exc()}")
        finally:
            progress.close()
            self.search_btn.setEnabled(True)

    def populate_results(self):
        self.table.setRowCount(len(self.results))
        for i, r in enumerate(self.results):
            vals = [
                i+1, f"{r['score']:.1f}", r["cod_id"], r["formula"], r["name"], r["space_group"],
                f"{r['matched']}/{r['used']}", f"{r['rms']:.4f}", r["ref_count"],
                f"{float(r['a']):.5g}", f"{float(r['b']):.5g}", f"{float(r['c']):.5g}", r["year"] or "",
                f"{float(r.get('chem_factor',1.0)):.3f}",
                f"{float(r.get('structure_factor',1.0)):.3f}"
            ]
            for c, v in enumerate(vals):
                item = QTableWidgetItem(str(v))
                if c == 2:
                    item.setData(Qt.ItemDataRole.UserRole, int(r["cod_id"]))
                self.table.setItem(i, c, item)
        if self.results:
            self.table.selectRow(0)
            self.populate_peak_choices()
        else:
            self.peak_choice_table.setRowCount(0)

    def selected_result(self):
        row = self.table.currentRow()
        if row < 0:
            return None
        cell = self.table.item(row, 2)
        if not cell:
            return None
        try:
            cod_id = int(cell.text())
        except Exception:
            return None
        return next((r for r in self.results if int(r["cod_id"]) == cod_id), None)

    def populate_peak_choices(self):
        """Show selectable measured↔COD reflection suggestions for the chosen candidate."""
        r = self.selected_result()
        self.peak_choice_table.setRowCount(0)
        if not r:
            self.peak_choice_info.setText("Wybierz powyżej kandydata COD, aby przejrzeć i wybrać proponowane przypisania pików.")
            return
        matches = list(r.get("matches", []) or [])
        self.peak_choice_table.setRowCount(len(matches))
        self.peak_choice_info.setText(
            f"COD {r.get('cod_id')} – {r.get('name') or r.get('formula') or 'struktura'}. "
            "Dopasowane piki są wstępnie zaznaczone; odznacz każde przypisanie, którego nie chcesz zaakceptować."
        )
        for row, m in enumerate(matches):
            use_item = QTableWidgetItem("")
            use_item.setFlags((use_item.flags() | Qt.ItemFlag.ItemIsUserCheckable | Qt.ItemFlag.ItemIsEnabled) & ~Qt.ItemFlag.ItemIsEditable)
            use_item.setCheckState(Qt.CheckState.Checked if bool(m.get("matched")) else Qt.CheckState.Unchecked)
            use_item.setData(Qt.ItemDataRole.UserRole, float(m.get("obs_2theta", math.nan)))
            self.peak_choice_table.setItem(row, 0, use_item)
            obs_tt = float(m.get("obs_2theta", math.nan))
            current = min(self.dataset.peaks, key=lambda p: abs(p.two_theta - obs_tt)) if self.dataset.peaks else None
            vals = [
                f"{obs_tt:.4f}", f"{float(m.get('obs_rel', 0.0)):.1f}",
                f"{float(m.get('ref_2theta', math.nan)):.4f}", str(m.get("hkl", "")),
                f"{float(m.get('delta', math.nan)):.4f}", "yes" if bool(m.get("matched")) else "no",
                str(current.phase if current is not None else ""), str(current.hkl if current is not None else ""),
            ]
            for c, v in enumerate(vals, start=1):
                item = QTableWidgetItem(v)
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.peak_choice_table.setItem(row, c, item)

    def select_all_matched_peak_choices(self):
        r = self.selected_result()
        if not r:
            return
        matches = list(r.get("matches", []) or [])
        for row, m in enumerate(matches):
            item = self.peak_choice_table.item(row, 0)
            if item is not None:
                item.setCheckState(Qt.CheckState.Checked if bool(m.get("matched")) else Qt.CheckState.Unchecked)

    def clear_peak_choices(self):
        for row in range(self.peak_choice_table.rowCount()):
            item = self.peak_choice_table.item(row, 0)
            if item is not None:
                item.setCheckState(Qt.CheckState.Unchecked)

    def apply_selected_peak_assignments(self):
        r = self.selected_result()
        if not r:
            QMessageBox.information(self, "Wybór pików COD", "Najpierw wybierz kandydata COD.")
            return
        selected = []
        matches = list(r.get("matches", []) or [])
        for row, m in enumerate(matches):
            item = self.peak_choice_table.item(row, 0)
            if item is not None and item.checkState() == Qt.CheckState.Checked:
                selected.append(m)
        if not selected:
            QMessageBox.information(self, "Wybór pików COD", "Najpierw zaznacz co najmniej jedno proponowane przypisanie piku.")
            return

        phase_name = str(r.get("plot_name") or r.get("name") or r.get("formula") or "COD structure")
        cod_id = int(r.get("cod_id"))
        existing = list(self.dataset.metadata.get("manual_cod_peak_assignments", []) or [])
        # Replace previous manual assignments for the same observed peaks, but preserve other
        # manually accepted peaks (possibly from another COD candidate).
        selected_tts = [float(m.get("obs_2theta")) for m in selected]
        kept = []
        for rec in existing:
            try:
                old_tt = float(rec.get("peak_two_theta"))
            except Exception:
                continue
            if not any(abs(old_tt - tt) <= max(0.08, self.tolerance.value()) for tt in selected_tts):
                kept.append(rec)

        applied = 0
        for m in selected:
            obs_tt = float(m.get("obs_2theta"))
            if not self.dataset.peaks:
                continue
            p = min(self.dataset.peaks, key=lambda q: abs(q.two_theta - obs_tt))
            if abs(p.two_theta - obs_tt) > max(0.08, self.tolerance.value()):
                continue
            p.phase = phase_name
            p.hkl = canonical_hkl_for_structure(str(m.get("hkl", "")), cod_structure_short(r.get("space_group", ""), r.get("sg_no")))
            p.delta_2theta = float(m.get("delta", math.nan))
            p.phase_source = "cod"
            kept.append({
                "peak_two_theta": float(p.two_theta), "cod_id": cod_id, "phase": phase_name,
                "formula": str(r.get("formula", "")), "space_group": str(r.get("space_group", "")),
                "hkl": p.hkl, "ref_two_theta": float(m.get("ref_2theta", math.nan)),
                "delta_2theta": float(m.get("delta", math.nan)), "tolerance_deg": float(self.tolerance.value()),
            })
            applied += 1
        self.dataset.metadata["manual_cod_peak_assignments"] = kept
        self.dataset.metadata_sources["manual_cod_peak_assignments"] = "user / COD peak selection"

        fits = list(self.dataset.metadata.get("phase_fits", []) or [])
        if not any(isinstance(f, dict) and int(f.get("cod_id", -1)) == cod_id for f in fits):
            fits.append({
                "phase": str(r.get("name") or phase_name), "plot_phase": phase_name,
                "formula": str(r.get("formula", "")), "score_percent": float(r.get("score", 0.0)),
                "matched_reflections": int(r.get("matched", 0)), "rms_2theta_deg": float(r.get("rms", math.nan)),
                "model": "COD manual peak selection", "cod_id": cod_id,
                "space_group": str(r.get("space_group", "")), "manual_selection": True,
            })
            self.dataset.metadata["phase_fits"] = fits

        if self.parent_window is not None:
            try:
                self.parent_window.update_table(); self.parent_window.redraw(); self.parent_window.update_summary(); self.parent_window.mark_project_dirty()
            except Exception:
                pass
        self.populate_peak_choices()
        QMessageBox.information(self, "Wybór pików COD", f"Zastosowano {applied} ręcznie wybranych przypisań pików COD do próbki.")

    def clear_manual_cod_assignments(self):
        if not self.dataset.metadata.get("manual_cod_peak_assignments"):
            QMessageBox.information(self, "Wybór pików COD", "Ta próbka nie ma ręcznych przypisań pików COD.")
            return
        self.dataset.metadata.pop("manual_cod_peak_assignments", None)
        self.dataset.metadata_sources.pop("manual_cod_peak_assignments", None)
        if self.parent_window is not None and hasattr(self.parent_window, "_assign_phases_for_dataset"):
            try:
                self.parent_window._assign_phases_for_dataset(self.dataset)
                self.parent_window.update_table(); self.parent_window.redraw(); self.parent_window.update_summary(); self.parent_window.mark_project_dirty()
            except Exception:
                clear_phase_assignments(self.dataset)
        self.populate_peak_choices()

    def show_match_details(self):
        r = self.selected_result()
        if not r:
            return
        lines = [
            f"COD ID: {r['cod_id']}",
            f"Wzór: {r['formula']}",
            f"Czytelna nazwa fazy: {r['name']}",
            f"Nazwa COD: {r.get('cod_name','')}",
            f"Grupa przestrzenna: {r['space_group']} (nr {r.get('sg_no','')})",
            f"Wynik: {r['score']:.1f}% | dopasowane: {r['matched']}/{r['used']} | RMS Δ2θ: {r['rms']:.4f}°",
            f"Mnożnik rankingu chemii: ×{float(r.get('chem_factor',1.0)):.3f} (może przekroczyć 1) | "
            f"pokrycie priorytetu składu: {100.0*float(r.get('chem_capture',0.0)):.1f}%",
            f"Uwagi do chemii: {r.get('chem_reason','')}",
            f"Mnożnik rankingu struktury osnowy: ×{float(r.get('structure_factor',1.0)):.3f} | "
            f"{r.get('structure_reason','')}",
            "",
            "Zmierzony 2θ\tZmierzone d [Å]\tIrel [%]\tCOD 2θ\tCOD d [Å]\thkl\tΔ2θ\tstatus"
        ]
        for m in r["matches"]:
            status = "DOPASOWANY" if m["matched"] else "poza tolerancją"
            lines.append(
                f"{m['obs_2theta']:.4f}\t{m['obs_d']:.5f}\t{m['obs_rel']:.1f}\t"
                f"{m['ref_2theta']:.4f}\t{m['ref_d']:.5f}\t{m['hkl']}\t{m['delta']:.4f}\t{status}"
            )
        dlg = QDialog(self); dlg.setWindowTitle(f"COD {r['cod_id']} – dopasowanie do zmierzonego dyfraktogramu")
        lay = QVBoxLayout(dlg); te = QTextEdit(); te.setReadOnly(True); te.setPlainText("\n".join(lines)); lay.addWidget(te)
        bb = QDialogButtonBox(QDialogButtonBox.StandardButton.Close); bb.button(QDialogButtonBox.StandardButton.Close).setText("Zamknij"); bb.rejected.connect(dlg.reject); bb.accepted.connect(dlg.accept); lay.addWidget(bb)
        dlg.resize(900, 560); dlg.exec()

    def show_cod_record_details(self):
        r = self.selected_result()
        if not r:
            QMessageBox.information(self, "Rekord COD", "Najpierw wybierz kandydata COD.")
            return
        try:
            d = self.cod_db.get_entry_details(int(r["cod_id"]))
        except Exception as e:
            QMessageBox.critical(self, "Rekord COD", str(e))
            return
        if not d:
            QMessageBox.information(self, "Rekord COD", "Nie znaleziono metadanych COD dla tego wpisu.")
            return

        friendly = r.get("name") or r.get("formula") or "struktura COD"
        fields = [
            ("Czytelna nazwa fazy", friendly),
            ("COD ID", d.get("file")),
            ("Wzór", d.get("formula")),
            ("Wzór obliczony", d.get("calcformula")),
            ("Nazwa zwyczajowa", d.get("commonname")),
            ("Nazwa chemiczna", d.get("chemname")),
            ("Nazwa minerału", d.get("mineral")),
            ("Grupa przestrzenna (H-M)", d.get("sg")),
            ("Grupa przestrzenna (Hall)", d.get("sgHall")),
            ("Numer grupy przestrzennej", d.get("sgNumber")),
            ("a [Å]", d.get("a")), ("b [Å]", d.get("b")), ("c [Å]", d.get("c")),
            ("alpha [°]", d.get("alpha")), ("beta [°]", d.get("beta")), ("gamma [°]", d.get("gamma")),
            ("Objętość komórki [Å³]", d.get("vol")),
            ("Temperatura komórki", d.get("celltemp")), ("Temperatura dyfrakcji", d.get("diffrtemp")),
            ("Ciśnienie komórki", d.get("cellpressure")), ("Ciśnienie podczas dyfrakcji", d.get("diffrpressure")),
            ("Źródło związku", d.get("compoundsource")), ("Liczba różnych pierwiastków", d.get("nel")),
            ("Metoda", d.get("method")), ("Promieniowanie", d.get("radiation")),
            ("Długość fali [Å]", d.get("wavelength")), ("Typ promieniowania", d.get("radType")),
            ("Autorzy", d.get("authors")), ("Tytuł publikacji", d.get("title")),
            ("Czasopismo", d.get("journal")), ("Rok", d.get("year")), ("Tom", d.get("volume")),
            ("Numer", d.get("issue")), ("Strony", f"{d.get('firstpage') or ''}–{d.get('lastpage') or ''}".strip("–")),
            ("DOI", d.get("doi")), ("R(all)", d.get("Rall")), ("R(obs)", d.get("Robs")),
            ("wR(all)", d.get("wRall")), ("wR(obs)", d.get("wRobs")),
            ("GOF(all)", d.get("gofall")), ("GOF(obs)", d.get("gofobs")),
            ("Status", d.get("status")), ("Flagi", d.get("flags")), ("Data COD", d.get("date")),
        ]
        lines = [f"{k}: {v}" for k, v in fields if v not in (None, "", "\\N", "?", ".")]
        dlg = QDialog(self)
        dlg.setWindowTitle(f"COD {r['cod_id']} – szczegóły rekordu")
        lay = QVBoxLayout(dlg)
        te = QTextEdit(); te.setReadOnly(True); te.setPlainText("\n".join(lines)); lay.addWidget(te)
        bb = QDialogButtonBox(QDialogButtonBox.StandardButton.Close); bb.button(QDialogButtonBox.StandardButton.Close).setText("Zamknij"); bb.rejected.connect(dlg.reject); bb.accepted.connect(dlg.accept); lay.addWidget(bb)
        dlg.resize(920, 680); dlg.exec()

    def import_selected(self):
        r = self.selected_result()
        if not r:
            QMessageBox.information(self, "Import COD", "Najpierw wybierz kandydata COD.")
            return
        try:
            ref = self.cod_db.entry_to_reference(int(r["cod_id"]))
            try:
                self.phase_db.add_phase(ref["name"], ref["formula"], ref["structure"], ref["space_group"], ref["category"], ref["notes"], ref["peaks"])
                imported_msg = f"Imported COD {r['cod_id']} to local phase references."
            except sqlite3.IntegrityError:
                imported_msg = f"COD {r['cod_id']} is already present in local phase references."
            if self.parent_window is not None and hasattr(self.parent_window, "_assign_phases_for_dataset"):
                try:
                    self.parent_window._assign_phases_for_dataset(self.dataset)
                except Exception:
                    clear_phase_assignments(self.dataset)
            else:
                assign_phases(self.dataset, self.phase_db, self.tolerance.value())
            if self.parent_window is not None:
                try:
                    self.parent_window.update_table(); self.parent_window.redraw(); self.parent_window.update_summary()
                except Exception:
                    pass
            QMessageBox.information(self, "Import COD", imported_msg + "\nPeak assignments were refreshed using the currently selected identification source.")
        except Exception as e:
            QMessageBox.critical(self, "Import COD", f"{e}\n\n{traceback.format_exc()}")

    def export_results(self):
        if not self.results:
            QMessageBox.information(self, "Eksport wyszukiwania COD", "Najpierw uruchom wyszukiwanie dla zmierzonego dyfraktogramu.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Eksportuj wyniki wyszukiwania COD", f"{self.dataset.name}_COD_search.xlsx", "Excel (*.xlsx);;CSV (*.csv)")
        if not path:
            return
        headers = ["Rank","Score [%]","COD ID","Formula","Czytelna nazwa fazy","Original COD name","Space group","Dopasowane piki","Used peaks","RMS delta 2theta [deg]","Generated reflections","a [A]","b [A]","c [A]","alpha [deg]","beta [deg]","gamma [deg]","Year","Method","Mnożnik chemii [×]","Pokrycie składu [%]"]
        rows = []
        for i,r in enumerate(self.results):
            rows.append([i+1,r["score"],r["cod_id"],r["formula"],r["name"],r.get("cod_name",""),r["space_group"],r["matched"],r["used"],r["rms"],r["ref_count"],r["a"],r["b"],r["c"],r["alpha"],r["beta"],r["gamma"],r["year"],r["method"],r.get('chem_factor',1.0),100*float(r.get('chem_capture',0.0))])
        ext = Path(path).suffix.lower()
        if ext == ".xlsx":
            if not OPENPYXL_OK:
                QMessageBox.warning(self, "Eksport", "Biblioteka openpyxl nie jest zainstalowana. Zapisz zamiast tego jako CSV.")
                return
            wb = Workbook(); ws = wb.active; ws.title = "Wyszukiwanie wzorca COD"; ws.append(headers)
            for row in rows: ws.append(row)
            wb.save(path)
        else:
            if not ext: path += ".csv"
            with open(path, "w", newline="", encoding="utf-8-sig") as f:
                w = csv.writer(f, delimiter=";"); w.writerow(headers); w.writerows(rows)
        QMessageBox.information(self, "Eksport wyszukiwania COD", f"Zapisano:\n{path}")



# --------------------------- multi-pattern plot creator ---------------------------


class PlotOrderConfigDialog(QDialog):
    """Konfiguracja kolejności próbek oraz dowolnej liczby linii podziału.

    Separator jest pełnoprawnym wierszem listy. Ma dwa niezależne, opcjonalne
    napisy: nad i pod linią. Checkbox separatora steruje wyłącznie widocznością
    napisów; sama linia pozostaje widoczna.
    """
    def __init__(self, creator, parent=None):
        super().__init__(parent or creator)
        self.creator = creator
        self.setWindowTitle("Konfiguracja kolejności materiałów i linii podziału")
        self.resize(1240, 760)
        self.entries = self._build_entries()

        root = QVBoxLayout(self)
        info = QLabel(
            "Próbki i linie podziału są w jednej kolejności. Dla wiersza LINIA PODZIAŁU pole po lewej stronie koloru "
            "oznacza napis NAD linią, a pole po prawej — napis POD linią. Oba napisy są opcjonalne. "
            "Checkbox w kolumnie „Użyj / etykiety” ukrywa lub pokazuje oba napisy, ale nie usuwa samej linii."
        )
        info.setWordWrap(True)
        root.addWidget(info)

        self.table = QTableWidget(0, 9)
        self.table.setHorizontalHeaderLabels([
            "Typ", "Użyj / etykiety", "Nazwa próbki / Napis nad linią", "Kolor",
            "Plik źródłowy / Napis pod linią", "Przesunięcie [%]",
            "Wygładź wykres", "Siła wygładzenia [%]", "Metoda wygładzania"
        ])
        hh = self.table.horizontalHeader()
        hh.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        hh.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        hh.setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        hh.setSectionResizeMode(3, QHeaderView.ResizeMode.ResizeToContents)
        hh.setSectionResizeMode(4, QHeaderView.ResizeMode.Stretch)
        hh.setSectionResizeMode(5, QHeaderView.ResizeMode.ResizeToContents)
        hh.setSectionResizeMode(6, QHeaderView.ResizeMode.ResizeToContents)
        hh.setSectionResizeMode(7, QHeaderView.ResizeMode.ResizeToContents)
        hh.setSectionResizeMode(8, QHeaderView.ResizeMode.ResizeToContents)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setEditTriggers(
            QAbstractItemView.EditTrigger.DoubleClicked
            | QAbstractItemView.EditTrigger.EditKeyPressed
            | QAbstractItemView.EditTrigger.AnyKeyPressed
        )
        self.table.verticalHeader().setDefaultSectionSize(36)
        root.addWidget(self.table, 1)

        row = QHBoxLayout()
        self.add_sep_btn = QPushButton("Dodaj linię podziału")
        self.del_sep_btn = QPushButton("Usuń wybraną linię")
        self.up_btn = QPushButton("W górę")
        self.down_btn = QPushButton("W dół")
        self.all_btn = QPushButton("Zaznacz wszystkie próbki")
        self.none_btn = QPushButton("Odznacz wszystkie próbki")
        self.reset_btn = QPushButton("Przywróć nazwy i kolory")
        self.up_btn.setToolTip("Przesuń zaznaczoną próbkę lub linię o jeden wiersz wyżej.")
        self.down_btn.setToolTip("Przesuń zaznaczoną próbkę lub linię o jeden wiersz niżej.")
        for b in (self.add_sep_btn, self.del_sep_btn, self.up_btn, self.down_btn,
                  self.all_btn, self.none_btn, self.reset_btn):
            row.addWidget(b)
        root.addLayout(row)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.button(QDialogButtonBox.StandardButton.Ok).setText("Zastosuj")
        buttons.button(QDialogButtonBox.StandardButton.Cancel).setText("Anuluj")
        root.addWidget(buttons)

        self.add_sep_btn.clicked.connect(self._add_separator)
        self.del_sep_btn.clicked.connect(self._delete_separator)
        self.up_btn.clicked.connect(lambda: self._move(-1))
        self.down_btn.clicked.connect(lambda: self._move(1))
        self.all_btn.clicked.connect(lambda: self._set_samples(True))
        self.none_btn.clicked.connect(lambda: self._set_samples(False))
        self.reset_btn.clicked.connect(self._reset_sample_styles)
        buttons.accepted.connect(self._accept_changes)
        buttons.rejected.connect(self.reject)
        self._populate()

    def _build_entries(self):
        by_pos = {}
        for sep in getattr(self.creator, 'separators', []):
            norm = normalizuj_separator(sep, len(self.creator.rows))
            by_pos.setdefault(norm['position'], []).append(norm)
        out = []
        for i, row in enumerate(self.creator.rows):
            for sep in by_pos.get(i, []):
                out.append({'kind': 'separator', 'separator': dict(sep)})
            out.append({'kind': 'sample', 'row': dict(row)})
        for sep in by_pos.get(len(self.creator.rows), []):
            out.append({'kind': 'separator', 'separator': dict(sep)})
        return out

    def _populate(self):
        self.table.blockSignals(True)
        self.table.clearContents()
        self.table.setRowCount(len(self.entries))
        for r, e in enumerate(self.entries):
            if e['kind'] == 'sample':
                info = e['row']
                typ = QTableWidgetItem("Próbka")
                typ.setFlags(typ.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(r, 0, typ)

                cb = QCheckBox()
                cb.setChecked(bool(info.get('enabled', True)))
                cb.setToolTip("Włącz lub wyłącz próbkę na wykresie.")
                wrap = QWidget(); lay = QHBoxLayout(wrap); lay.setContentsMargins(6,0,6,0)
                lay.addWidget(cb); lay.addStretch(1); self.table.setCellWidget(r, 1, wrap)

                name = QTableWidgetItem(str(info.get('label', info['dataset'].name)))
                name.setToolTip("Nazwa próbki wyświetlana na wykresie. Dwuklik pozwala ją edytować.")
                self.table.setItem(r, 2, name)

                color = str(info.get('color', '#1f77b4'))
                cbtn = QPushButton(color); cbtn.setStyleSheet(f"background-color: {color};")
                cbtn.setToolTip("Wybierz kolor krzywej tej próbki.")
                cbtn.clicked.connect(lambda checked=False, rr=r: self._choose_color(rr))
                self.table.setCellWidget(r, 3, cbtn)

                src = QTableWidgetItem(info['dataset'].path.name)
                src.setFlags(src.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(r, 4, src)
                off = QTableWidgetItem("—"); off.setFlags(off.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(r, 5, off)
                cb.stateChanged.connect(lambda state, rr=r: self._sample_checked(rr, state))
                smoothing = display_smoothing_settings(info)
                info.update(smoothing)
                smooth_cb = QCheckBox()
                smooth_cb.setChecked(smoothing['smooth_enabled'])
                smooth_cb.setToolTip(
                    "Wygładzenie tylko wizualne. Nie zmienia danych pomiarowych, wykrytych pików, "
                    "FWHM ani wyników analizy. Procent oznacza udział krzywej wygładzonej: "
                    "y = (1 − p/100)·y + (p/100)·y_wygładzone. Metodę wybierasz w ostatniej kolumnie."
                )
                wrap = QWidget(); lay = QHBoxLayout(wrap); lay.setContentsMargins(6,0,6,0)
                lay.addWidget(smooth_cb); lay.addStretch(1); self.table.setCellWidget(r, 6, wrap)
                smooth_spin = QDoubleSpinBox()
                smooth_spin.setRange(0.0, 100.0); smooth_spin.setDecimals(1)
                smooth_spin.setSuffix(" %"); smooth_spin.setValue(smoothing['smooth_strength_percent'])
                smooth_spin.setEnabled(smoothing['smooth_enabled'])
                smooth_spin.setToolTip(
                    "0%: oryginalna krzywa; 100%: pełne wygładzenie. Nominalne okno "
                    f"{smoothing['smooth_window_deg']:g}° 2θ, wielomian stopnia {smoothing['smooth_polyorder']}. "
                    "Okno jest nieparzyste i ograniczone długością danych. "
                    "Dla nierównej osi używana jest pomocnicza interpolacja. Mniej niż 5 punktów: bez zmian."
                )
                self.table.setCellWidget(r, 7, smooth_spin)
                method_combo = QComboBox()
                for method_key in SMOOTH_METHODS:
                    method_combo.addItem(SMOOTH_METHOD_LABELS[method_key], method_key)
                method_combo.setCurrentIndex(max(0, SMOOTH_METHODS.index(smoothing['smooth_method'])))
                method_combo.setEnabled(smoothing['smooth_enabled'])
                method_combo.setToolTip(
                    "Z zachowaniem pików: to samo wygładzenie tła, ale wąskie piki zachowują wysokość "
                    "i szerokość (klasyczny filtr spłaszcza piki węższe niż ok. 0,25°). "
                    "Klasyczny Savitzky–Golay: jednakowe okno dla całego wykresu, jak w wersjach do 20.3.1."
                )
                self.table.setCellWidget(r, 8, method_combo)
                smooth_cb.toggled.connect(smooth_spin.setEnabled)
                smooth_cb.toggled.connect(method_combo.setEnabled)
                smooth_cb.toggled.connect(lambda checked, rr=r: self._smoothing_changed(rr, 'smooth_enabled', checked))
                smooth_spin.valueChanged.connect(lambda value, rr=r: self._smoothing_changed(rr, 'smooth_strength_percent', value))
                method_combo.currentIndexChanged.connect(
                    lambda index, rr=r, combo=method_combo: self._smoothing_changed(rr, 'smooth_method', combo.itemData(index)))
            else:
                sep = normalizuj_separator(e['separator'], len(self.creator.rows))
                e['separator'] = sep
                typ = QTableWidgetItem("LINIA PODZIAŁU")
                typ.setFlags(typ.flags() & ~Qt.ItemFlag.ItemIsEditable)
                typ.setToolTip("Separator pomiędzy grupami materiałów.")
                self.table.setItem(r, 0, typ)

                cb = QCheckBox()
                cb.setChecked(bool(sep.get('labels_visible', True)))
                cb.setToolTip("Pokaż lub ukryj napis nad i pod linią. Sama linia podziału pozostaje widoczna.")
                wrap = QWidget(); lay = QHBoxLayout(wrap); lay.setContentsMargins(6,0,6,0)
                lay.addWidget(cb); lay.addStretch(1); self.table.setCellWidget(r, 1, wrap)
                cb.stateChanged.connect(lambda state, rr=r: self._separator_labels_checked(rr, state))

                top_item = QTableWidgetItem(str(sep.get('top_label', '')))
                top_item.setToolTip("Napis wyświetlany nad linią podziału, np. „Materiały odniesienia”. Pole może pozostać puste.")
                self.table.setItem(r, 2, top_item)

                color = str(sep.get('color', DEFAULT_SEPARATOR_COLOR))
                cbtn = QPushButton(color); cbtn.setStyleSheet(f"background-color: {color};")
                cbtn.setToolTip("Wybierz kolor linii podziału i jej napisów.")
                cbtn.clicked.connect(lambda checked=False, rr=r: self._choose_color(rr))
                self.table.setCellWidget(r, 3, cbtn)

                bottom_item = QTableWidgetItem(str(sep.get('bottom_label', '')))
                bottom_item.setToolTip("Napis wyświetlany pod linią podziału, np. „Materiały badane”. Pole może pozostać puste.")
                self.table.setItem(r, 4, bottom_item)

                spin = QSpinBox(); spin.setRange(-48, 48)
                spin.setValue(int(round(float(sep.get('offset_pct', 0.0)))))
                spin.setSuffix(" %")
                spin.setToolTip(
                    "Przesuwa separator w pionie w obrębie luki pomiędzy sąsiednimi dyfraktogramami. "
                    "0% oznacza środek luki; zakres ±48% pozwala zbliżyć linię do sąsiedniego przebiegu."
                )
                spin.valueChanged.connect(lambda val, rr=r: self._separator_offset(rr, val))
                self.table.setCellWidget(r, 5, spin)
                for column in (6, 7, 8):
                    self.table.removeCellWidget(r, column)
                    empty = QTableWidgetItem("—")
                    empty.setFlags(empty.flags() & ~Qt.ItemFlag.ItemIsEditable)
                    self.table.setItem(r, column, empty)
        self.table.blockSignals(False)

    def _sync_texts(self):
        for r, e in enumerate(self.entries):
            if e['kind'] == 'sample':
                it = self.table.item(r, 2)
                if it:
                    e['row']['label'] = it.text().strip() or e['row']['dataset'].name
            else:
                a = self.table.item(r, 2); b = self.table.item(r, 4)
                e['separator']['top_label'] = a.text().strip() if a else ''
                e['separator']['bottom_label'] = b.text().strip() if b else ''

    def _sample_checked(self, row, state):
        if 0 <= row < len(self.entries) and self.entries[row]['kind'] == 'sample':
            self.entries[row]['row']['enabled'] = (state == Qt.CheckState.Checked.value)

    def _smoothing_changed(self, row, key, value):
        if 0 <= row < len(self.entries) and self.entries[row]['kind'] == 'sample':
            self.entries[row]['row'][key] = value

    def _separator_labels_checked(self, row, state):
        if 0 <= row < len(self.entries) and self.entries[row]['kind'] == 'separator':
            self.entries[row]['separator']['labels_visible'] = (state == Qt.CheckState.Checked.value)

    def _separator_offset(self, row, val):
        if 0 <= row < len(self.entries) and self.entries[row]['kind'] == 'separator':
            self.entries[row]['separator']['offset_pct'] = float(val)

    def _choose_color(self, row):
        if not (0 <= row < len(self.entries)):
            return
        entry = self.entries[row]
        if entry['kind'] == 'sample':
            target = entry['row']; default = '#1f77b4'; title = "Wybierz kolor krzywej"
        else:
            target = entry['separator']; default = DEFAULT_SEPARATOR_COLOR; title = "Wybierz kolor linii podziału"
        current = str(target.get('color', default) or default)
        col = QColorDialog.getColor(QColor(current), self, title)
        if col.isValid():
            self._sync_texts()
            target['color'] = col.name()
            self._populate()
            self.table.selectRow(row)

    def _add_separator(self):
        self._sync_texts()
        r = self.table.currentRow(); ins = (r + 1 if r >= 0 else len(self.entries))
        self.entries.insert(ins, {
            'kind': 'separator',
            'separator': normalizuj_separator({
                'top_label': '', 'bottom_label': '', 'offset_pct': 0.0,
                'color': DEFAULT_SEPARATOR_COLOR, 'labels_visible': True,
            }, len(self.creator.rows))
        })
        self._populate(); self.table.selectRow(ins)

    def _delete_separator(self):
        self._sync_texts(); r = self.table.currentRow()
        if 0 <= r < len(self.entries) and self.entries[r]['kind'] == 'separator':
            self.entries.pop(r); self._populate()
            if self.entries:
                self.table.selectRow(min(r, len(self.entries)-1))
        else:
            QMessageBox.information(self, "Linia podziału", "Zaznacz wiersz linii podziału, którą chcesz usunąć.")

    def _move(self, direction):
        self._sync_texts(); r = self.table.currentRow(); nr = r + direction
        if r < 0 or nr < 0 or nr >= len(self.entries):
            return
        self.entries[r], self.entries[nr] = self.entries[nr], self.entries[r]
        self._populate(); self.table.selectRow(nr)

    def _set_samples(self, enabled):
        self._sync_texts()
        for e in self.entries:
            if e['kind'] == 'sample':
                e['row']['enabled'] = bool(enabled)
        self._populate()

    def _reset_sample_styles(self):
        self._sync_texts(); i = 0
        for e in self.entries:
            if e['kind'] == 'sample':
                e['row']['label'] = e['row']['dataset'].name
                e['row']['color'] = self.creator.DEFAULT_COLORS[i % len(self.creator.DEFAULT_COLORS)]
                i += 1
        self._populate()

    def _accept_changes(self):
        self._sync_texts(); rows = []; seps = []; count = 0
        for e in self.entries:
            if e['kind'] == 'sample':
                rows.append(e['row']); count += 1
            else:
                sep = dict(e['separator']); sep['position'] = count
                seps.append(normalizuj_separator(sep, len(self.creator.rows)))
        self.creator.rows = rows
        self.creator.separators = seps
        self.creator._update_order_summary()
        self.creator.redraw()
        self.accept()


# Kreator wykresów – konfiguracja publikacyjna
class PlotCreatorDialog(QDialog):
    """Publication-oriented creator for plotting any subset of loaded XRD datasets.

    The dialog keeps its state serializable so names, colours, ordering and figure
    settings can be stored in an XRD Analyzer project.
    """

    DEFAULT_COLORS = [
        "#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e",
        "#17becf", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22",
    ]

    def __init__(self, datasets: List[Dataset], parent=None, state: Optional[dict] = None):
        super().__init__(parent)
        self.setWindowTitle("Kreator wykresów – wiele dyfraktogramów XRD")
        self.setWindowFlags(
            self.windowFlags()
            | Qt.WindowType.WindowMaximizeButtonHint
            | Qt.WindowType.WindowMinimizeButtonHint
        )
        self.resize(1500, 920)
        self.datasets = datasets
        self.parent_window = parent
        state = sanitize_plot_state(state)
        self._restoring_state = True
        self._last_synced_state = None
        self.rows = self._build_rows_from_state(state or {})
        self.separators = self._build_separators_from_state(state or {})
        self._last_axes = []
        # Manual label positions in Kreator wykresów. Keys are stable per layout/sample/peak
        # and values store the text-anchor position. They are serialized with the project.
        self._label_positions: Dict[str, dict] = {}
        self._annotation_entries: List[dict] = []
        self._drag_info: Optional[dict] = None
        self._legend_artists: List[object] = []
        self._manual_legend_loc = None
        self._manual_legend_positions = {}

        outer = QVBoxLayout(self)
        split = QSplitter(Qt.Orientation.Horizontal)
        outer.addWidget(split, 1)

        # ---- Left: sample selection + styling ----
        left = QWidget(); ll = QVBoxLayout(left)
        split.addWidget(left)

        order_box = QGroupBox("Kolejność materiałów i linie podziału")
        order_lay = QVBoxLayout(order_box)
        self.order_summary_label = QLabel("")
        self.order_summary_label.setWordWrap(True)
        self.order_config_btn = QPushButton("Konfiguracja kolejności materiałów…")
        self.order_config_btn.setToolTip("Otwiera dużą tabelę próbek. W tej samej tabeli można dodawać wiele linii podziału i przesuwać je przyciskami W górę / W dół.")
        order_lay.addWidget(self.order_summary_label)
        order_lay.addWidget(self.order_config_btn)
        ll.addWidget(order_box)

        cg = QGroupBox("Konfiguracja wykresu")
        cf = QFormLayout(cg)
        self.mode_combo = QComboBox()
        self.mode_combo.addItems(["Nałożone – wszystkie na jednym wykresie", "Ułożone pionowo", "Oddzielne panele"])
        self.data_combo = QComboBox()
        self.data_combo.addItems([
            "Surowa intensywność pomiarowa",
            "Intensywność po korekcji tła",
            "Surowa znormalizowana (0–1)",
            "Po korekcji tła + normalizacja (0–1)",
        ])
        self.offset_spin = QDoubleSpinBox()
        self.offset_spin.setRange(0.00, 20.0)
        self.offset_spin.setDecimals(2)
        self.offset_spin.setSingleStep(0.10)
        self.offset_spin.setValue(1.15)
        self.offset_spin.setToolTip(
            "Używane tylko w układzie pionowym. Zmienia odległość między poziomami bazowymi kolejnych dyfraktogramów. "
            "Przy włączonym zachowaniu wysokości pików figura rośnie w pionie, dzięki czemu piki nie stają się wizualnie mniejsze."
        )
        self.keep_peak_height_cb = QCheckBox("Zachowaj wysokość pików przy zmianie odstępu")
        self.keep_peak_height_cb.setChecked(True)
        self.keep_peak_height_cb.setToolTip(
            "Zalecane. W układzie pionowym wysokość podglądu i eksportu rośnie automatycznie wraz z odstępem pionowym. "
            "Zwiększa to fizyczną odległość między dyfraktogramami bez ściskania znormalizowanych pików."
        )
        self.effective_height_label = QLabel("Efektywna wysokość figury: —")
        self.effective_height_label.setToolTip(
            "Rzeczywista wysokość używana w podglądzie i eksporcie po kompensacji odstępów. "
            "Pole wysokości figury pozostaje wysokością odniesienia dla domyślnego odstępu."
        )
        self.linewidth_spin = QDoubleSpinBox(); self.linewidth_spin.setRange(0.2, 5.0); self.linewidth_spin.setDecimals(2); self.linewidth_spin.setValue(1.2)
        self.legend_cb = QCheckBox("Pokaż legendę"); self.legend_cb.setChecked(True)
        self.legend_pos_combo = QComboBox()
        self.legend_pos_combo.addItems([
            "Najlepsze automatyczne", "Prawy górny", "Lewy górny", "Lewy dolny", "Prawy dolny",
            "Prawy środek", "Lewy środek", "Dolny środek", "Górny środek", "Środek",
            "Na zewnątrz po prawej", "Ręcznie – przeciągnij w podglądzie",
        ])
        self.legend_pos_combo.setCurrentText("Najlepsze automatyczne")
        self.legend_pos_combo.setToolTip("Wybierz stałe położenie legendy, umieść ją poza wykresem albo wybierz tryb ręczny i przeciągnij legendę w podglądzie.")
        self.legend_font_spin = QSpinBox(); self.legend_font_spin.setRange(5, 24); self.legend_font_spin.setValue(8); self.legend_font_spin.setSuffix(" pt")
        self.sample_names_cb = QCheckBox("Pokaż nazwy próbek przy dyfraktogramach"); self.sample_names_cb.setChecked(True)
        self.sample_names_cb.setToolTip("Przydatne dla wykresów ułożonych pionowo i oddzielnych paneli. Gdy nazwy próbek są pokazane przy krzywych, legendę można zwykle wyłączyć.")
        self.grid_cb = QCheckBox("Pokaż siatkę"); self.grid_cb.setChecked(True)
        self.peaks_cb = QCheckBox("Pokaż znaczniki wykrytych pików")
        self.peak_labels_cb = QCheckBox("Pokaż etykiety pików")
        self.peak_label_style_combo = QComboBox(); self.peak_label_style_combo.addItems(PEAK_LABEL_STYLES); self.peak_label_style_combo.setCurrentText(PEAK_LABEL_STYLE_FULL)
        self.peak_label_style_combo.setToolTip("Publikacyjny – krótki: γ(111), Cu–Ni(200), α-Ti(101). Publikacyjny – pełny: γ-austenite FCC (111). Baza / diagnostyka: pełna faza + ID COD + hkl + 2θ.")
        self.hkl_display_combo = QComboBox(); self.hkl_display_combo.addItems(HKL_DISPLAY_MODES); self.hkl_display_combo.setCurrentText(HKL_DISPLAY_AUTO)
        self.hkl_display_combo.setToolTip("Wybierz, czy indeksy Millera mają być ukryte, pokazane z przypisanego źródła fazy, czy ograniczone do bazy wewnętrznej/COD.")
        self.peak_label_position_combo = QComboBox()
        self.peak_label_position_combo.addItems(["Nad pikiem", "Pod pikiem", "Po prawej stronie piku", "Po lewej stronie piku"])
        self.peak_label_position_combo.setCurrentText("Nad pikiem")
        self.peak_label_position_combo.setToolTip(
            "Automatyczne położenie początkowe etykiet pików. Ręcznie przeciągnięte etykiety zachowują pozycje do użycia funkcji „Przywróć położenie etykiet”."
        )
        self.peak_label_orientation_combo = QComboBox()
        self.peak_label_orientation_combo.addItems(["Pionowo", "Poziomo"])
        self.peak_label_orientation_combo.setCurrentText("Pionowo")
        self.peak_label_orientation_combo.setToolTip("Orientacja tekstu etykiet pików. Linia prowadząca pozostaje połączona z pikiem podczas przesuwania etykiety.")
        self.peak_font_spin = QSpinBox(); self.peak_font_spin.setRange(5, 24); self.peak_font_spin.setValue(8); self.peak_font_spin.setSuffix(" pt")
        cf.addRow("Układ:", self.mode_combo)
        cf.addRow("Dane:", self.data_combo)
        cf.addRow("Odstęp pionowy:", self.offset_spin)
        cf.addRow(self.keep_peak_height_cb)
        cf.addRow(self.effective_height_label)
        cf.addRow("Grubość linii:", self.linewidth_spin)
        cf.addRow(self.legend_cb)
        cf.addRow("Położenie legendy:", self.legend_pos_combo)
        cf.addRow("Czcionka legendy:", self.legend_font_spin)
        cf.addRow(self.sample_names_cb)
        cf.addRow(self.grid_cb)
        cf.addRow(self.peaks_cb)
        cf.addRow(self.peak_labels_cb)
        cf.addRow("Styl etykiet pików:", self.peak_label_style_combo)
        cf.addRow("Indeksy Millera (hkl):", self.hkl_display_combo)
        cf.addRow("Położenie etykiet pików:", self.peak_label_position_combo)
        cf.addRow("Orientacja etykiet pików:", self.peak_label_orientation_combo)
        cf.addRow("Czcionka etykiet pików:", self.peak_font_spin)
        self.reset_label_positions_btn = QPushButton("Przywróć położenie etykiet")
        self.reset_label_positions_btn.setToolTip("Usuń wszystkie ręczne przesunięcia etykiet pików i przywróć rozmieszczenie automatyczne")
        cf.addRow(self.reset_label_positions_btn)
        self.reset_figure_defaults_btn = QPushButton("Przywróć domyślne")
        self.reset_figure_defaults_btn.setToolTip("Przywróć domyślny układ, tryb danych oraz ustawienia linii i pików")
        cf.addRow(self.reset_figure_defaults_btn)
        self.publication_template_btn = QPushButton("Zastosuj szablon porównawczy do publikacji")
        self.publication_template_btn.setToolTip("Utwórz porównanie ułożone pionowo podobne do XRD_all_8_publication_EN: nazwy próbek przy krzywych, opcjonalne linie podziału i pasek znaczników faz, bez zbędnej legendy.")
        cf.addRow(self.publication_template_btn)
        self.phase_marker_strip_cb = QCheckBox("Pokaż pasek znaczników faz na dole")
        cf.addRow(self.phase_marker_strip_cb)
        lg = QGroupBox("Etykiety i rozmiar")
        lf = QFormLayout(lg)
        self.title_edit = QLineEdit("Dyfraktogramy XRD")
        self.subtitle_edit = QLineEdit("")
        self.subtitle_edit.setPlaceholderText("Opcjonalny drugi wiersz tytułu, np. Cu Kα, λ = 0,15406 nm")
        self.xlabel_edit = QLineEdit("2θ [°]")
        self.ylabel_edit = QLineEdit("Intensywność [j.u.]")
        self.font_family_combo = QComboBox(); self.font_family_combo.setEditable(True)
        self.font_family_combo.addItems(["DejaVu Sans", "Arial", "Times New Roman", "Calibri"])
        self.font_family_combo.setCurrentText("DejaVu Sans")
        self.font_family_combo.setToolTip("Krój czcionki używany dla tytułu, osi, wartości osi, nazw próbek, legendy i etykiet pików. Jeśli czcionka jest niedostępna, Matplotlib użyje czcionki zastępczej.")
        self.title_font_spin = QSpinBox(); self.title_font_spin.setRange(6, 36); self.title_font_spin.setValue(12); self.title_font_spin.setSuffix(" pt")
        self.axis_font_spin = QSpinBox(); self.axis_font_spin.setRange(6, 28); self.axis_font_spin.setValue(10); self.axis_font_spin.setSuffix(" pt")
        self.tick_font_spin = QSpinBox(); self.tick_font_spin.setRange(5, 24); self.tick_font_spin.setValue(9); self.tick_font_spin.setSuffix(" pt")
        self.sample_font_spin = QSpinBox(); self.sample_font_spin.setRange(5, 24); self.sample_font_spin.setValue(9); self.sample_font_spin.setSuffix(" pt")
        self.width_spin = QDoubleSpinBox(); self.width_spin.setRange(4.0, 30.0); self.width_spin.setValue(11.0); self.width_spin.setSuffix(" in")
        self.height_spin = QDoubleSpinBox(); self.height_spin.setRange(3.0, 30.0); self.height_spin.setValue(7.5); self.height_spin.setSuffix(" in")
        self.dpi_combo = QComboBox(); self.dpi_combo.addItems(["300", "600", "900", "1200"]); self.dpi_combo.setCurrentText("600")
        lf.addRow("Tytuł:", self.title_edit)
        lf.addRow("Podtytuł:", self.subtitle_edit)
        lf.addRow("Oś X:", self.xlabel_edit)
        lf.addRow("Oś Y:", self.ylabel_edit)
        lf.addRow("Krój czcionki:", self.font_family_combo)
        lf.addRow("Czcionka tytułu:", self.title_font_spin)
        lf.addRow("Czcionka opisów osi:", self.axis_font_spin)
        lf.addRow("Czcionka wartości osi:", self.tick_font_spin)
        lf.addRow("Czcionka nazw próbek:", self.sample_font_spin)
        lf.addRow("Szerokość figury:", self.width_spin)
        lf.addRow("Wysokość figury:", self.height_spin)
        lf.addRow("DPI rastra:", self.dpi_combo)
        self.reset_labels_defaults_btn = QPushButton("Przywróć domyślne")
        self.reset_labels_defaults_btn.setToolTip("Przywróć domyślny tytuł, opisy osi, rozmiar figury i DPI rastra")
        lf.addRow(self.reset_labels_defaults_btn)
        # Keep the control panel readable: Figure configuration on the left and
        # Labels and size on the right. The sample table remains above both columns.
        settings_columns = QHBoxLayout()
        settings_columns.setSpacing(10)
        settings_columns.addWidget(cg, 1)
        settings_columns.addWidget(lg, 1)
        ll.addLayout(settings_columns)

        hint = QLabel("Podgląd Kreatora wykresów ma stałą skalę: kółko myszy nie powiększa wykresu. Gdy figura jest większa od podglądu, użyj pionowego lub poziomego paska przewijania. Etykiety pików można przeciągać lewym przyciskiem myszy; linia prowadząca pozostaje połączona z pikiem. Dwuklik przy piku pozwala edytować tekst.")
        hint.setWordWrap(True)
        ll.addWidget(hint)

        export_row = QHBoxLayout()
        self.refresh_btn = QPushButton("Odśwież podgląd")
        self.save_btn = QPushButton("Zapisz figurę…")
        self.close_btn = QPushButton("Zamknij")
        export_row.addWidget(self.refresh_btn); export_row.addWidget(self.save_btn); export_row.addWidget(self.close_btn)
        ll.addLayout(export_row)

        # ---- Right: live Matplotlib preview ----
        # The canvas lives inside a non-resizing QScrollArea. This is important for
        # publication figures with a large physical height: the preview keeps the
        # requested aspect/size and a vertical scrollbar lets the user reach lower
        # panels and drag their annotations instead of Qt shrinking the whole figure.
        right = QWidget(); rl = QVBoxLayout(right)
        split.addWidget(right)
        self.fig = Figure(figsize=(11, 7.5))
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self.preview_scroll = QScrollArea()
        self.preview_scroll.setWidgetResizable(False)
        self.preview_scroll.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop)
        self.preview_scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.preview_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.preview_scroll.setFrameShape(QFrame.Shape.StyledPanel)
        self.preview_scroll.setWidget(self.canvas)
        self.preview_scroll.setToolTip(
            "Użyj pasków przewijania, aby poruszać się po dużym podglądzie. Przewijanie nie zmienia powiększenia ani osi wykresu."
        )
        rl.addWidget(self.preview_scroll, 1)
        split.setSizes([780, 1000]); split.setStretchFactor(0, 0); split.setStretchFactor(1, 1)

        self.refresh_btn.clicked.connect(self.redraw)
        self.save_btn.clicked.connect(self.save_figure)
        self.close_btn.clicked.connect(self.accept)
        self.reset_figure_defaults_btn.clicked.connect(self.reset_figure_defaults)
        self.reset_labels_defaults_btn.clicked.connect(self.reset_labels_defaults)
        self.reset_label_positions_btn.clicked.connect(self.reset_label_positions)
        self.publication_template_btn.clicked.connect(self.apply_publication_template)
        self.phase_marker_strip_cb.stateChanged.connect(lambda *_: self.redraw())
        self.mode_combo.currentIndexChanged.connect(self._mode_changed)
        self.data_combo.currentIndexChanged.connect(self.redraw)
        self.offset_spin.valueChanged.connect(self.redraw)
        self.keep_peak_height_cb.stateChanged.connect(self.redraw)
        self.linewidth_spin.valueChanged.connect(self.redraw)
        self.legend_cb.stateChanged.connect(self._legend_setting_changed)
        self.legend_pos_combo.currentTextChanged.connect(self._legend_setting_changed)
        self.legend_font_spin.valueChanged.connect(self.redraw)
        self.sample_names_cb.stateChanged.connect(self.redraw)
        self.grid_cb.stateChanged.connect(self.redraw)
        self.peaks_cb.stateChanged.connect(self.redraw)
        self.peak_labels_cb.stateChanged.connect(self.redraw)
        self.peak_label_style_combo.currentTextChanged.connect(self.redraw)
        self.hkl_display_combo.currentTextChanged.connect(self.redraw)
        self.peak_label_position_combo.currentTextChanged.connect(self.redraw)
        self.peak_label_orientation_combo.currentTextChanged.connect(self.redraw)
        self.peak_font_spin.valueChanged.connect(self.redraw)
        self.title_edit.editingFinished.connect(self.redraw)
        self.subtitle_edit.editingFinished.connect(self.redraw)
        self.xlabel_edit.editingFinished.connect(self.redraw)
        self.ylabel_edit.editingFinished.connect(self.redraw)
        self.font_family_combo.currentTextChanged.connect(self.redraw)
        self.title_font_spin.valueChanged.connect(self.redraw)
        self.axis_font_spin.valueChanged.connect(self.redraw)
        self.tick_font_spin.valueChanged.connect(self.redraw)
        self.sample_font_spin.valueChanged.connect(self.redraw)
        self.width_spin.valueChanged.connect(self.redraw)
        self.height_spin.valueChanged.connect(self.redraw)
        # Deliberately no scroll_event connection here: Kreator wykresów must keep fixed positioning.
        self.canvas.mpl_connect("button_press_event", self._on_mouse_press)
        self.canvas.mpl_connect("motion_notify_event", self._on_mouse_motion)
        self.canvas.mpl_connect("button_release_event", self._on_mouse_release)

        self.order_config_btn.clicked.connect(self._open_order_config)
        self._apply_state(state or {})
        self._update_order_summary()
        self.legend_pos_combo.setEnabled(self.legend_cb.isChecked())
        self.legend_font_spin.setEnabled(self.legend_cb.isChecked())
        self._restoring_state = False
        self._mode_changed()
        self._last_synced_state = deepcopy(self.get_state())
        if self.parent_window is not None:
            self.parent_window._plot_creator_dialog = self
        # Text and DPI edits must also count, even without Enter or a redraw.
        for widget in self.findChildren(QLineEdit):
            widget.textChanged.connect(self._sync_to_project)
        for widget in self.findChildren(QComboBox):
            widget.currentTextChanged.connect(self._sync_to_project)
        for widget in self.findChildren(QAbstractSpinBox):
            widget.editingFinished.connect(self._sync_to_project)

    def _sync_to_project(self, *_args):
        if self._restoring_state or self._last_synced_state is None:
            return
        state = self.get_state()
        parent = self.parent_window
        if parent is not None:
            parent.plot_creator_state = deepcopy(state)
            if state != self._last_synced_state:
                parent.project_dirty = True
        self._last_synced_state = deepcopy(state)

    def done(self, result):
        self._sync_to_project()
        parent = self.parent_window
        if parent is not None:
            if getattr(parent, "_plot_creator_dialog", None) is self:
                parent._plot_creator_dialog = None
            if getattr(parent, "project_dirty", False):
                parent.mark_project_dirty()
        super().done(result)

    def accept(self):
        self.done(QDialog.DialogCode.Accepted)

    def reject(self):
        # X / Escape close this editor; they do not discard the project edits.
        self.done(QDialog.DialogCode.Rejected)

    def reset_sample_styles(self):
        for i, row in enumerate(self.rows):
            row["label"] = row["dataset"].name
            row["color"] = self.DEFAULT_COLORS[i % len(self.DEFAULT_COLORS)]
        self._update_order_summary()
        self.redraw()

    def reset_label_positions(self):
        """Restore automatic Kreator wykresów peak-label placement."""
        self._label_positions.clear()
        self._drag_info = None
        self.redraw()

    def reset_figure_defaults(self):
        self.mode_combo.setCurrentIndex(0)
        self.data_combo.setCurrentIndex(0)
        self.offset_spin.setValue(1.15)
        self.keep_peak_height_cb.setChecked(True)
        self.linewidth_spin.setValue(1.20)
        self.legend_cb.setChecked(True)
        self.legend_pos_combo.setCurrentText("Najlepsze automatyczne")
        self.legend_font_spin.setValue(8)
        self._manual_legend_loc = None
        self.sample_names_cb.setChecked(True)
        self.grid_cb.setChecked(True)
        self.peaks_cb.setChecked(False)
        self.peak_labels_cb.setChecked(False)
        self.peak_label_style_combo.setCurrentText(PEAK_LABEL_STYLE_FULL)
        self.hkl_display_combo.setCurrentText(HKL_DISPLAY_AUTO)
        self.peak_label_position_combo.setCurrentText("Nad pikiem")
        self.peak_label_orientation_combo.setCurrentText("Pionowo")
        self.peak_font_spin.setValue(8)
        self._mode_changed()
        self.redraw()

    def reset_labels_defaults(self):
        self.title_edit.setText("Dyfraktogramy XRD")
        self.xlabel_edit.setText("2θ [°]")
        self.ylabel_edit.setText("Intensywność [j.u.]")
        self.font_family_combo.setCurrentText("DejaVu Sans")
        self.title_font_spin.setValue(12)
        self.axis_font_spin.setValue(10)
        self.tick_font_spin.setValue(9)
        self.sample_font_spin.setValue(9)
        self.width_spin.setValue(11.0)
        self.height_spin.setValue(7.5)
        self.dpi_combo.setCurrentText("600")
        self.redraw()

    def _build_rows_from_state(self, state: dict):
        saved = state.get("rows", []) if isinstance(state, dict) else []
        saved = [r for r in saved if isinstance(r, dict)] if isinstance(saved, list) else []
        rows = []
        used = set()
        ds_by_id = {ds.dataset_id: ds for ds in self.datasets}
        ds_by_path = {}
        for ds in self.datasets:
            paths = [str(ds.path.resolve())] + list(getattr(ds, "_project_source_paths", []))
            for path in paths:
                ds_by_path.setdefault(str(Path(path).resolve()), []).append(ds)
        for r in saved:
            dataset_id = r.get("dataset_id")
            ds = ds_by_id.get(dataset_id) if isinstance(dataset_id, str) else None
            # Older projects only have a path. Do not guess ambiguous matches.
            if ds is None and not r.get("dataset_id") and isinstance(r.get("path"), str):
                candidates = {d.dataset_id: d for d in ds_by_path.get(str(Path(r["path"]).resolve()), [])}
                if len(candidates) == 1:
                    ds = next(iter(candidates.values()))
            if ds is None or ds.dataset_id in used:
                continue
            used.add(ds.dataset_id)
            rows.append({
                "dataset": ds,
                "enabled": optional_bool(r.get("enabled", True), True),
                "label": str(r.get("label", ds.name)),
                "color": valid_color(r.get("color"), self.DEFAULT_COLORS[len(rows) % len(self.DEFAULT_COLORS)]),
                **display_smoothing_settings(r),
            })
        for ds in self.datasets:
            if ds.dataset_id in used:
                continue
            rows.append({
                "dataset": ds,
                "enabled": True,
                "label": ds.name,
                "color": self.DEFAULT_COLORS[len(rows) % len(self.DEFAULT_COLORS)],
                **display_smoothing_settings({}),
            })
        return rows

    def _build_separators_from_state(self, state: dict):
        raw = state.get("separators", []) if isinstance(state, dict) else []
        out = []
        if isinstance(raw, list):
            for item in raw:
                if isinstance(item, dict):
                    out.append(normalizuj_separator(item, len(self.rows)))

        # Najstarszy format: pojedynczy separator grupy poza listą separators.
        if not out and bool(state.get("group_separator", False)):
            try:
                pos = int(state.get("group_break_after", min(4, max(1, len(self.rows)))))
            except Exception:
                pos = min(4, max(1, len(self.rows)))
            out.append(normalizuj_separator({
                "position": pos,
                "top_label": state.get("top_group_label", "Reference materials"),
                "bottom_label": state.get("bottom_group_label", "Dissimilar welded joints"),
                "offset_pct": 0.0,
                "color": DEFAULT_SEPARATOR_COLOR,
                "labels_visible": True,
            }, len(self.rows)))
        return out

    def _open_order_config(self):
        dlg=PlotOrderConfigDialog(self,self)
        dlg.exec()

    def _update_order_summary(self):
        enabled=sum(1 for r in self.rows if r.get("enabled",True))
        smoothed = sum(1 for r in self._selected_rows() if r.get("smooth_enabled") and r.get("smooth_strength_percent", 0) > 0)
        self.order_summary_label.setText(
            f"Próbki: {len(self.rows)} (używane: {enabled})   •   Linie podziału: {len(self.separators)}"
            + (f"\nWygładzenie wizualne: {smoothed} próbek. Dane analityczne pozostają niezmienione." if smoothed else "")
        )

    def _selected_rows(self):
        """Zwróć aktywne próbki Kreatora wykresów w bieżącej kolejności.

        Od v20 tabela została przeniesiona do osobnego okna konfiguracji.
        Starsza implementacja pobierała zaznaczenie bezpośrednio z tabeli,
        dlatego po refaktoryzacji metoda przypadkowo zniknęła i Kreator wykresów
        kończył się błędem AttributeError. Źródłem prawdy jest teraz self.rows.
        """
        return [
            row for row in self.rows
            if bool(row.get("enabled", True)) and row.get("dataset") is not None
        ]

    def _apply_state(self, state: dict):
        restoring = self._restoring_state
        self._restoring_state = True
        blockers = [QSignalBlocker(widget) for widget in self.findChildren(QWidget)]
        try:
            # Old artists belong to the previous state and must not overwrite it.
            self._legend_artists = []
            self._apply_state_values(sanitize_plot_state(state))
        finally:
            del blockers
            self._restoring_state = restoring

    def _apply_state_values(self, state: dict):
        if not state:
            return
        raw_positions = state.get("label_positions", {})
        if isinstance(raw_positions, dict):
            self._label_positions = {str(k): dict(v) for k, v in raw_positions.items() if isinstance(v, dict)}
        self.mode_combo.setCurrentIndex(int(state.get("mode", 0)))
        self.data_combo.setCurrentIndex(int(state.get("data_mode", 0)))
        self.offset_spin.setValue(float(state.get("vertical_spacing", 1.15)))
        self.keep_peak_height_cb.setChecked(bool(state.get("keep_peak_height_with_spacing", True)))
        self.linewidth_spin.setValue(float(state.get("line_width", 1.2)))
        self.legend_cb.setChecked(bool(state.get("legend", True)))
        legend_pos = str(state.get("legend_position", "Najlepsze automatyczne"))
        legend_pos = {"Best":"Najlepsze automatyczne","Upper right":"Prawy górny","Upper left":"Lewy górny","Lower left":"Lewy dolny","Lower right":"Prawy dolny","Center right":"Prawy środek","Center left":"Lewy środek","Lower center":"Dolny środek","Upper center":"Górny środek","Center":"Środek","Outside right":"Na zewnątrz po prawej","Manual – drag in preview":"Ręcznie – przeciągnij w podglądzie"}.get(legend_pos,legend_pos)
        if self.legend_pos_combo.findText(legend_pos) >= 0:
            self.legend_pos_combo.setCurrentText(legend_pos)
        self.legend_font_spin.setValue(int(state.get("legend_font_size", 8)))
        self.sample_names_cb.setChecked(bool(state.get("sample_names", True)))
        manual_loc = state.get("manual_legend_loc", None)
        if isinstance(manual_loc, (list, tuple)) and len(manual_loc) == 2:
            try:
                self._manual_legend_loc = (float(manual_loc[0]), float(manual_loc[1]))
            except Exception:
                self._manual_legend_loc = None
        positions = state.get("manual_legend_positions", {})
        if isinstance(positions, dict):
            for key, value in positions.items():
                try:
                    if len(value) == 2 and all(math.isfinite(float(v)) for v in value):
                        self._manual_legend_positions[str(key)] = [float(v) for v in value]
                except (TypeError, ValueError, OverflowError):
                    continue
        self.grid_cb.setChecked(bool(state.get("grid", True)))
        self.peaks_cb.setChecked(bool(state.get("peak_markers", False)))
        self.peak_labels_cb.setChecked(bool(state.get("peak_labels", False)))
        lbl_style = _migrate_peak_label_style(str(state.get("peak_label_style", PEAK_LABEL_STYLE_FULL)))
        if self.peak_label_style_combo.findText(lbl_style) >= 0:
            self.peak_label_style_combo.setCurrentText(lbl_style)
        hkl_mode = _migrate_hkl_display(str(state.get("hkl_display", HKL_DISPLAY_AUTO)))
        if self.hkl_display_combo.findText(hkl_mode) >= 0:
            self.hkl_display_combo.setCurrentText(hkl_mode)
        peak_pos = str(state.get("peak_label_position", "Nad pikiem"))
        peak_pos = {"Above peak":"Nad pikiem","Below peak":"Pod pikiem","Right of peak":"Po prawej stronie piku","Left of peak":"Po lewej stronie piku"}.get(peak_pos,peak_pos)
        if self.peak_label_position_combo.findText(peak_pos) >= 0:
            self.peak_label_position_combo.setCurrentText(peak_pos)
        peak_orient = str(state.get("peak_label_orientation", "Pionowo"))
        peak_orient = {"Vertical":"Pionowo","Horizontal":"Poziomo"}.get(peak_orient,peak_orient)
        if self.peak_label_orientation_combo.findText(peak_orient) >= 0:
            self.peak_label_orientation_combo.setCurrentText(peak_orient)
        self.peak_font_spin.setValue(int(state.get("peak_font_size", 8)))
        self.separators = self._build_separators_from_state(state)
        self.phase_marker_strip_cb.setChecked(bool(state.get("phase_marker_strip", False)))
        self.title_edit.setText(str(state.get("title", "Dyfraktogramy XRD")))
        self.subtitle_edit.setText(str(state.get("subtitle", "")))
        self.xlabel_edit.setText(str(state.get("xlabel", "2θ [°]")))
        self.ylabel_edit.setText(str(state.get("ylabel", "Intensywność [j.u.]")))
        self.font_family_combo.setCurrentText(str(state.get("font_family", "DejaVu Sans")))
        self.title_font_spin.setValue(int(state.get("title_font_size", 12)))
        self.axis_font_spin.setValue(int(state.get("axis_font_size", 10)))
        self.tick_font_spin.setValue(int(state.get("tick_font_size", 9)))
        self.sample_font_spin.setValue(int(state.get("sample_font_size", 9)))
        self.width_spin.setValue(float(state.get("width_in", 11.0)))
        self.height_spin.setValue(float(state.get("height_in", 7.5)))
        self.dpi_combo.setCurrentText(str(state.get("dpi", "600")))

    def get_state(self) -> dict:
        self._capture_manual_legend_position()
        return {
            "rows": [
                {
                    "path": str(r["dataset"].path.resolve()),
                    "dataset_id": r["dataset"].dataset_id,
                    "enabled": bool(r["enabled"]),
                    "label": r["label"],
                    "color": r["color"],
                    **display_smoothing_settings(r),
                }
                for r in self.rows
            ],
            "mode": self.mode_combo.currentIndex(),
            "data_mode": self.data_combo.currentIndex(),
            "vertical_spacing": self.offset_spin.value(),
            "keep_peak_height_with_spacing": self.keep_peak_height_cb.isChecked(),
            "line_width": self.linewidth_spin.value(),
            "legend": self.legend_cb.isChecked(),
            "legend_position": self.legend_pos_combo.currentText(),
            "legend_font_size": self.legend_font_spin.value(),
            "manual_legend_loc": list(self._manual_legend_loc) if isinstance(self._manual_legend_loc, tuple) else None,
            "manual_legend_positions": deepcopy(self._manual_legend_positions),
            "sample_names": self.sample_names_cb.isChecked(),
            "grid": self.grid_cb.isChecked(),
            "peak_markers": self.peaks_cb.isChecked(),
            "peak_labels": self.peak_labels_cb.isChecked(),
            "peak_label_style": self.peak_label_style_combo.currentText(),
            "hkl_display": self.hkl_display_combo.currentText(),
            "peak_label_position": self.peak_label_position_combo.currentText(),
            "peak_label_orientation": self.peak_label_orientation_combo.currentText(),
            "peak_font_size": self.peak_font_spin.value(),
            "separators": [normalizuj_separator(x, len(self.rows)) for x in self.separators],
            "phase_marker_strip": self.phase_marker_strip_cb.isChecked(),
            "title": self.title_edit.text(),
            "subtitle": self.subtitle_edit.text(),
            "xlabel": self.xlabel_edit.text(),
            "ylabel": self.ylabel_edit.text(),
            "font_family": self.font_family_combo.currentText(),
            "title_font_size": self.title_font_spin.value(),
            "axis_font_size": self.axis_font_spin.value(),
            "tick_font_size": self.tick_font_spin.value(),
            "sample_font_size": self.sample_font_spin.value(),
            "width_in": self.width_spin.value(),
            "height_in": self.height_spin.value(),
            "dpi": self.dpi_combo.currentText(),
            "label_positions": deepcopy(self._label_positions),
        }

    def _figure_title_text(self):
        """Return the Kreator wykresów title with an optional subtitle on a second line.

        Kept as a PlotCreatorDialog method because redraw() can be called during dialog
        initialization before the user interacts with any control.
        """
        title = self.title_edit.text().strip()
        subtitle = self.subtitle_edit.text().strip()
        if title and subtitle:
            return title + "\n" + subtitle
        return title or subtitle

    def apply_publication_template(self):
        """Configure a clean stacked comparison inspired by the user's reference figure."""
        self.mode_combo.setCurrentText('Ułożone pionowo')
        self.data_combo.setCurrentText('Po korekcji tła + normalizacja (0–1)')
        self.offset_spin.setValue(1.15)
        self.keep_peak_height_cb.setChecked(True)
        self.legend_cb.setChecked(False)
        self.sample_names_cb.setChecked(True)
        self.grid_cb.setChecked(True)
        self.peaks_cb.setChecked(False)
        self.peak_labels_cb.setChecked(False)
        if len(self._selected_rows()) >= 4 and not self.separators:
            self.separators = [normalizuj_separator({
                "position": min(4, len(self.rows)-1),
                "top_label": "Reference materials",
                "bottom_label": "Dissimilar welded joints",
                "offset_pct": 0.0,
                "color": DEFAULT_SEPARATOR_COLOR,
                "labels_visible": True,
            }, len(self.rows))]
            self._update_order_summary()
        self.phase_marker_strip_cb.setChecked(True)
        self.ylabel_edit.setText('Normalized intensity + vertical offset')
        sel=self._selected_rows()
        waves={round(float(i['dataset'].wavelength_A),6) for i in sel if getattr(i['dataset'],'wavelength_A',None)}
        if len(waves)==1:
            lam_A=next(iter(waves))
            self.title_edit.setText("XRD patterns")
            self.subtitle_edit.setText(f"Cu Kα, λ = {lam_A/10.0:.5f} nm")
        self.redraw()

    def _draw_phase_marker_strip(self, ax, selected, spacing, bottom_rows=3):
        if not self.phase_marker_strip_cb.isChecked():
            return 0.0
        # Collect median measured positions by publication phase/hkl. This produces the same
        # visual language as reference tick strips without pretending that a fixed database
        # angle applies to every alloy composition.
        phase_map={}
        phase_colors={}
        for info in selected:
            ds=info['dataset']
            for p in ds.peaks:
                short,full=_publication_phase_names(p.phase,ds.material_type)
                h=_clean_hkl_label(p.hkl)
                if not full or not h: continue
                key=full
                phase_map.setdefault(key,{}).setdefault(h,[]).append(float(p.two_theta))
                phase_colors.setdefault(key,info['color'])
        if not phase_map:
            return 0.0
        # Prefer familiar matrix phases first and cap the strip to four rows.
        def pri(k):
            kl=k.lower()
            return (0 if 'austenite' in kl else 1 if 'cu–ni' in kl or 'cu-ni' in kl else 2 if 'ferrite' in kl else 3, k)
        keys=sorted(phase_map,key=pri)[:4]
        x0,x1=ax.get_xlim(); xr=max(1e-9,x1-x0)
        row_step=0.18
        base=-0.20
        for j,key in enumerate(keys):
            yy=base-j*row_step
            col=phase_colors.get(key,'0.25')
            ax.text(x0+0.007*xr,yy,key,ha='left',va='center',fontsize=max(6,self.tick_font_spin.value()-1),fontfamily=self.font_family_combo.currentText())
            for h,vals in sorted(phase_map[key].items(),key=lambda kv:np.median(kv[1])):
                xx=float(np.median(vals))
                ax.vlines(xx,yy-0.055,yy+0.055,color=col,lw=max(1.0,self.linewidth_spin.value()))
                ax.text(xx,yy-0.07,h.strip('()'),ha='center',va='top',rotation=90,fontsize=max(5,self.tick_font_spin.value()-2),fontfamily=self.font_family_combo.currentText())
        return abs(base-(len(keys)-1)*row_step)+0.18

    def _mode_changed(self):
        # Vertical spacing has a meaning only for stacked plots. In overlay mode
        # all curves deliberately share the same zero level.
        mode = self.mode_combo.currentIndex()
        stacked = mode == 1
        self.offset_spin.setEnabled(stacked)
        self.keep_peak_height_cb.setEnabled(stacked)
        self.effective_height_label.setVisible(stacked)
        self.sample_names_cb.setEnabled(mode in (1, 2))
        self.redraw()

    def _legend_setting_changed(self, *_args):
        if self._restoring_state:
            return
        enabled = self.legend_cb.isChecked()
        self.legend_pos_combo.setEnabled(enabled)
        self.legend_font_spin.setEnabled(enabled)
        if self.legend_pos_combo.currentText() != "Ręcznie – przeciągnij w podglądzie":
            self._manual_legend_loc = None
            self._manual_legend_positions.clear()
        self.redraw()

    def _capture_manual_legend_position(self):
        if self.legend_pos_combo.currentText() != "Ręcznie – przeciągnij w podglądzie":
            return
        for index, leg in enumerate(list(getattr(self, "_legend_artists", []))):
            try:
                loc = getattr(leg, "_loc", None)
                if isinstance(loc, (tuple, list)) and len(loc) == 2:
                    x, y = float(loc[0]), float(loc[1])
                    if math.isfinite(x) and math.isfinite(y):
                        if index == 0:
                            self._manual_legend_loc = (x, y)
                        key = getattr(leg, "_xrd_legend_key", None)
                        if key is not None:
                            self._manual_legend_positions[key] = [x, y]
            except Exception:
                pass

    def _legend_kwargs(self):
        text = self.legend_pos_combo.currentText()
        mapping = {
            "Najlepsze automatyczne": "best", "Prawy górny": "upper right", "Lewy górny": "upper left",
            "Lewy dolny": "lower left", "Prawy dolny": "lower right",
            "Prawy środek": "center right", "Lewy środek": "center left",
            "Dolny środek": "lower center", "Górny środek": "upper center", "Środek": "center",
        }
        if text == "Na zewnątrz po prawej":
            return {"loc": "upper left", "bbox_to_anchor": (1.01, 1.0), "borderaxespad": 0.0}
        if text == "Ręcznie – przeciągnij w podglądzie":
            return {"loc": self._manual_legend_loc if self._manual_legend_loc is not None else "upper right"}
        return {"loc": mapping.get(text, "best")}

    def _add_legend(self, ax, dataset=None):
        if not self.legend_cb.isChecked():
            return None
        kwargs = self._legend_kwargs()
        key = f"m{self.mode_combo.currentIndex()}"
        if dataset is not None:
            key += f"|id:{dataset.dataset_id}"
        if (self.legend_pos_combo.currentText() == "Ręcznie – przeciągnij w podglądzie"
                and key in self._manual_legend_positions):
            kwargs["loc"] = tuple(self._manual_legend_positions[key])
        prop = {"size": self.legend_font_spin.value(), "family": self.font_family_combo.currentText()}
        leg = ax.legend(prop=prop, **kwargs)
        if leg is not None:
            leg._xrd_legend_key = key
            self._legend_artists.append(leg)
            if self.legend_pos_combo.currentText() == "Ręcznie – przeciągnij w podglądzie":
                try:
                    leg.set_draggable(True, use_blit=False, update="loc")
                except TypeError:
                    try:
                        leg.set_draggable(True)
                    except Exception:
                        pass
                except Exception:
                    pass
        return leg

    def _apply_axis_typography(self, ax):
        family = self.font_family_combo.currentText()
        try:
            ax.tick_params(axis="both", labelsize=self.tick_font_spin.value())
            for lab in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
                lab.set_fontfamily(family)
        except Exception:
            pass

    def _curve(self, ds: Dataset, row=None):
        mode = self.data_combo.currentIndex()
        raw = np.asarray(ds.y, float)
        if mode in (1, 3):
            if ds.corrected_y is not None:
                y = np.asarray(ds.corrected_y, float)
            else:
                try:
                    base = baseline_asls(raw, x=ds.x)
                    y = np.clip(raw - base, 0, None)
                except Exception:
                    y = raw - np.nanmin(raw)
        else:
            y = raw.copy()
        # Work only on display values. Dataset arrays and scientific results are
        # never assigned here; resampling is interpolated back to original x.
        if row is None:
            row = next((item for item in self.rows if item['dataset'] is ds), {})
        smoothing = display_smoothing_settings(row)
        if smoothing['smooth_enabled']:
            y = smooth_for_display(ds.x, y, smoothing['smooth_strength_percent'],
                                   smoothing['smooth_window_deg'], smoothing['smooth_polyorder'],
                                   method=smoothing['smooth_method'])
        if mode in (2, 3):
            lo = float(np.nanmin(y)); hi = float(np.nanmax(y))
            y = (y - lo) / (hi - lo) if hi > lo else np.zeros_like(y)
        return y

    def _peak_text(self, ds: Dataset, p: PeakResult) -> str:
        return format_peak_label(ds, p, self.peak_label_style_combo.currentText(), self.hkl_display_combo.currentText())

    def _peak_label_key(self, ds: Dataset, p: PeakResult) -> str:
        # Keep the key independent of the displayed phase text so an edited label or
        # database-name change does not discard the user's manual position.
        key = f"m{self.mode_combo.currentIndex()}|id:{ds.dataset_id}|{p.two_theta:.12g}"
        legacy = f"m{self.mode_combo.currentIndex()}|{ds.name}|{p.two_theta:.4f}"
        if key not in self._label_positions and legacy in self._label_positions:
            matches = [d for d in self.datasets if d.name == ds.name]
            peaks = [peak for peak in ds.peaks if f"{peak.two_theta:.4f}" == f"{p.two_theta:.4f}"]
            if len(matches) == 1 and len(peaks) == 1:
                self._label_positions[key] = self._label_positions.pop(legacy)
        return key

    def _register_annotation(self, ann, ax, ds: Dataset, p: PeakResult, coord_mode: str,
                             band_low: Optional[float] = None, band_high: Optional[float] = None):
        entry = {
            "annotation": ann, "ax": ax, "dataset": ds, "peak": p,
            "key": self._peak_label_key(ds, p), "coord_mode": coord_mode,
            "band_low": band_low, "band_high": band_high,
        }
        self._annotation_entries.append(entry)
        try:
            ann.set_picker(True)
        except Exception:
            pass
        return entry

    def _expand_for_peak_labels(self, ax):
        if not self.peak_labels_cb.isChecked():
            return
        try:
            y0, y1 = ax.get_ylim()
            span = y1 - y0
            if span > 0:
                # A small amount of headroom is enough because tall-peak labels are
                # positioned in axes coordinates and forced to grow downwards.
                ax.set_ylim(y0, y1 + 0.10 * span)
        except Exception:
            pass

    def _draw_peak_labels(self, ax, ds: Dataset, y: np.ndarray, offset: float = 0.0,
                          color: Optional[str] = None, band_low: Optional[float] = None,
                          band_high: Optional[float] = None):
        """Draw draggable peak labels with user-selectable placement and orientation.

        ``Peak-label position`` determines the automatic starting point (above, below,
        right or left of the measured peak). ``Peak-label orientation`` controls only
        text rotation. The Matplotlib annotation leader line always connects the label
        to the peak and therefore follows the label during manual dragging.

        Manually saved positions take precedence over automatic placement. In stacked
        mode the y coordinate is clamped to the pattern's own band so labels cannot be
        dragged into another diffractogram.
        """
        if not self.peak_labels_cb.isChecked() or not ds.peaks:
            return
        valid = [p for p in ds.peaks if 0 <= p.index < len(y)]
        if not valid:
            return

        y0, y1 = ax.get_ylim()
        span = max(1e-12, y1 - y0)
        mixed = blended_transform_factory(ax.transData, ax.transAxes)
        ordered = sorted(valid, key=lambda q: q.two_theta)
        last_x = None
        lane_flip = False
        stacked_band = band_low is not None and band_high is not None
        placement = self.peak_label_position_combo.currentText()
        orientation = self.peak_label_orientation_combo.currentText()
        rotation = 90 if orientation == "Pionowo" else 0
        x_min = float(np.nanmin(ds.x)); x_max = float(np.nanmax(ds.x))
        x_span = max(1e-9, x_max - x_min)
        dx = 0.018 * x_span

        def alignment_for(where: str):
            if where == "Po prawej stronie piku":
                return "left", "center"
            if where == "Po lewej stronie piku":
                return "right", "center"
            if where == "Pod pikiem":
                return "center", "top"
            return "center", "bottom"  # Above peak

        ha, va = alignment_for(placement)

        for p in ordered:
            py = float(y[p.index]) + offset
            label = self._peak_text(ds, p)
            if not label:
                continue
            if last_x is not None and abs(p.two_theta - last_x) < 3.0:
                lane_flip = not lane_flip
            else:
                lane_flip = False
            last_x = p.two_theta
            key = self._peak_label_key(ds, p)
            saved = self._label_positions.get(key, {})

            if stacked_band:
                lo = float(band_low); hi = float(band_high)
                if hi <= lo:
                    hi = lo + max(0.2, abs(span) * 0.02)
                bspan = max(1e-9, hi - lo)
                margin = max(0.025, 0.035 * bspan)
                dy = max(0.055, 0.085 * bspan)

                auto_x = float(p.two_theta)
                auto_y = py
                if placement == "Po prawej stronie piku":
                    auto_x += dx
                    auto_y += (0.05 * bspan if lane_flip else 0.0)
                elif placement == "Po lewej stronie piku":
                    auto_x -= dx
                    auto_y += (0.05 * bspan if lane_flip else 0.0)
                elif placement == "Pod pikiem":
                    auto_y -= dy * (1.45 if lane_flip else 1.0)
                else:  # Above peak
                    auto_y += dy * (1.45 if lane_flip else 1.0)

                tx = float(saved.get("x", auto_x)) if saved.get("coord") == "data" else auto_x
                ty = float(saved.get("y", auto_y)) if saved.get("coord") == "data" else auto_y
                tx = min(max(tx, x_min), x_max)
                ty = min(max(ty, lo + margin), hi - margin)

                ann = ax.annotate(
                    label, xy=(p.two_theta, py), xycoords="data",
                    xytext=(tx, ty), textcoords="data",
                    ha=ha, va=va, fontsize=self.peak_font_spin.value(),
                    fontfamily=self.font_family_combo.currentText(),
                    rotation=rotation, rotation_mode="anchor",
                    arrowprops=dict(arrowstyle="-", lw=0.55, color=color or "0.35"),
                    annotation_clip=True,
                )
                self._register_annotation(ann, ax, ds, p, "data", lo, hi)
            else:
                frac = (py - y0) / span
                frac = min(max(frac, 0.04), 0.96)
                dy_frac = 0.085 * (1.35 if lane_flip else 1.0)
                auto_x = float(p.two_theta)
                auto_frac = frac
                if placement == "Po prawej stronie piku":
                    auto_x += dx
                    auto_frac += (0.035 if lane_flip else 0.0)
                elif placement == "Po lewej stronie piku":
                    auto_x -= dx
                    auto_frac += (0.035 if lane_flip else 0.0)
                elif placement == "Pod pikiem":
                    auto_frac -= dy_frac
                else:  # Above peak
                    auto_frac += dy_frac

                auto_x = min(max(auto_x, x_min), x_max)
                auto_frac = min(max(auto_frac, 0.05), 0.97)
                tx = float(saved.get("x", auto_x)) if saved.get("coord") == "axes_y" else auto_x
                ty = float(saved.get("y", auto_frac)) if saved.get("coord") == "axes_y" else auto_frac
                tx = min(max(tx, x_min), x_max)
                ty = min(max(ty, 0.05), 0.97)

                ann = ax.annotate(
                    label, xy=(p.two_theta, py), xycoords="data",
                    xytext=(tx, ty), textcoords=mixed,
                    ha=ha, va=va, fontsize=self.peak_font_spin.value(),
                    fontfamily=self.font_family_combo.currentText(),
                    rotation=rotation, rotation_mode="anchor",
                    arrowprops=dict(arrowstyle="-", lw=0.55, color=color or "0.35"),
                    annotation_clip=True,
                )
                self._register_annotation(ann, ax, ds, p, "axes_y")
            try:
                ann.set_clip_on(True)
                ann.set_clip_path(ax.patch)
            except Exception:
                pass

    def _effective_stacked_height(self, selected_count: int, mode: int) -> float:
        """Return the physical figure height used for preview/export.

        A stacked plot stores every normalized trace in a common data axis. Merely
        increasing the numerical offset expands the Y data range; on a fixed-height
        canvas that necessarily makes every unit-high peak look smaller. When the
        preservation option is enabled, scale the physical canvas by the same change
        in stacked Y range. Thus a normalized peak keeps approximately the same
        physical height while the baseline-to-baseline distance really increases.
        """
        base_h = float(self.height_spin.value())
        if mode != 1 or not self.keep_peak_height_cb.isChecked() or selected_count <= 1:
            return base_h
        spacing = max(0.01, float(self.offset_spin.value()))
        ref_spacing = 1.15
        # 1.30 represents one normalized trace plus the small fixed top/bottom
        # headroom used by the stacked axis. It prevents spacing from scaling margins.
        current_span = (selected_count - 1) * spacing + 1.30
        reference_span = (selected_count - 1) * ref_spacing + 1.30
        if reference_span <= 0:
            return base_h
        return max(3.0, base_h * current_span / reference_span)

    def _update_effective_height_label(self, effective_h: float, mode: int):
        if not hasattr(self, "effective_height_label"):
            return
        if mode == 1:
            if self.keep_peak_height_cb.isChecked():
                self.effective_height_label.setText(f"Effective figure height: {effective_h:.2f} in")
            else:
                self.effective_height_label.setText(
                    f"Effective figure height: {float(self.height_spin.value()):.2f} in (fixed)"
                )
        else:
            self.effective_height_label.setText("Effective figure height: —")

    def _sync_preview_canvas_size(self):
        """Keep the preview canvas at the requested physical figure size.

        FigureCanvas normally follows the available Qt layout space. In Kreator wykresów
        that made a tall figure shrink to the viewport, so lower panels could not be
        reached for annotation editing. A fixed-size canvas inside QScrollArea gives
        the preview a real scrollable working surface while leaving Matplotlib zoom
        disabled.
        """
        try:
            dpi = float(self.fig.dpi or 100.0)
            # Use the figure's ACTUAL size. In stacked mode the height may be larger
            # than the reference Figure height in order to preserve peak size.
            fig_w, fig_h = self.fig.get_size_inches()
            w_px = max(320, int(round(float(fig_w) * dpi)))
            h_px = max(240, int(round(float(fig_h) * dpi)))
            self.canvas.setFixedSize(w_px, h_px)
        except Exception:
            pass

    @staticmethod
    def _bbox_koliduje_z_linia(bbox, line, margines_px=2.0):
        """Sprawdź lokalną kolizję prostokąta tekstu z krzywą Line2D w pikselach."""
        try:
            path = line.get_path().transformed(line.get_transform())
            xy = np.asarray(path.vertices, float)
            if xy.size == 0:
                return False
            mask = ((xy[:, 0] >= bbox.x0 - margines_px) & (xy[:, 0] <= bbox.x1 + margines_px))
            if not np.any(mask):
                return False
            local_y = xy[mask, 1]
            return bool(np.any((local_y >= bbox.y0 - margines_px) & (local_y <= bbox.y1 + margines_px)))
        except Exception:
            return False

    def _draw_separator(self, ax, sep: dict, sep_y: float):
        """Narysuj linię oraz dwie niezależne etykiety z odstępem w punktach typograficznych."""
        sep = normalizuj_separator(sep, len(self.rows))
        color = str(sep.get("color", DEFAULT_SEPARATOR_COLOR))
        line = ax.axhline(sep_y, color=color, lw=0.8, alpha=0.75, zorder=0.4)
        result = {"line": line, "texts": [], "sep_y": float(sep_y), "sep": sep}
        if not bool(sep.get("labels_visible", True)):
            return result
        top = str(sep.get("top_label", "")).strip()
        bottom = str(sep.get("bottom_label", "")).strip()
        if not top and not bottom:
            return result

        mix = blended_transform_factory(ax.transAxes, ax.transData)
        fontsize = max(7, self.tick_font_spin.value())
        pad_pt = max(3.0, float(fontsize) * 0.35)
        if top:
            ann_top = ax.annotate(
                top, xy=(0.985, sep_y), xycoords=mix,
                xytext=(0, pad_pt), textcoords="offset points",
                ha="right", va="bottom", fontsize=fontsize,
                fontfamily=self.font_family_combo.currentText(), color=color,
                zorder=7, annotation_clip=True,
            )
            result["texts"].append(ann_top)
        if bottom:
            ann_bottom = ax.annotate(
                bottom, xy=(0.985, sep_y), xycoords=mix,
                xytext=(0, -pad_pt), textcoords="offset points",
                ha="right", va="top", fontsize=fontsize,
                fontfamily=self.font_family_combo.currentText(), color=color,
                zorder=7, annotation_clip=True,
            )
            result["texts"].append(ann_bottom)
        return result

    def _resolve_separator_collisions(self, separator_entries, trace_lines):
        """Przenieś etykiety separatorów w poziomie, gdy kolidują z pikami lub krzywymi."""
        if not separator_entries:
            return
        try:
            renderer = self.canvas.get_renderer()
        except Exception:
            return

        peak_boxes = []
        for entry in self._annotation_entries:
            ann = entry.get("annotation")
            if ann is None:
                continue
            try:
                peak_boxes.append(ann.get_window_extent(renderer).expanded(1.04, 1.08))
            except Exception:
                pass

        occupied = []
        x_candidates = [0.985, 0.82, 0.68, 0.52]
        for record in separator_entries:
            texts = [t for t in record.get("texts", []) if t is not None]
            if not texts:
                continue
            chosen = None
            for xfrac in x_candidates:
                for ann in texts:
                    ann.xy = (xfrac, float(record.get("sep_y", 0.0)))
                try:
                    boxes = [ann.get_window_extent(renderer).expanded(1.02, 1.06) for ann in texts]
                except Exception:
                    boxes = []
                if not boxes:
                    chosen = xfrac
                    break
                collision = False
                all_other = peak_boxes + occupied
                for box in boxes:
                    if any(box.overlaps(other) for other in all_other):
                        collision = True; break
                    if any(self._bbox_koliduje_z_linia(box, line, 2.0) for line in trace_lines):
                        collision = True; break
                if not collision:
                    chosen = xfrac
                    occupied.extend(boxes)
                    break
            if chosen is None:
                # Ostatni fallback: pozostaw położenie, ale dodaj subtelne tło poprawiające czytelność.
                chosen = x_candidates[-1]
                for ann in texts:
                    ann.xy = (chosen, float(record.get("sep_y", 0.0)))
                    try:
                        ann.set_bbox(dict(facecolor=ann.axes.get_facecolor(), edgecolor="none", alpha=0.88, pad=1.0))
                    except Exception:
                        pass
                try:
                    occupied.extend([ann.get_window_extent(renderer).expanded(1.02, 1.06) for ann in texts])
                except Exception:
                    pass

    @staticmethod
    def _data_dy_from_points(ax, points: float) -> float:
        """Przelicz pionowe przesunięcie w punktach typograficznych na jednostki danych osi Y."""
        try:
            px = float(points) * float(ax.figure.dpi) / 72.0
            inv = ax.transData.inverted()
            y0 = inv.transform((0.0, 0.0))[1]
            y1 = inv.transform((0.0, px))[1]
            return float(y1 - y0)
        except Exception:
            return 0.0

    def redraw(self):
        if self._restoring_state:
            return
        selected = self._selected_rows()
        # Preserve the user's viewport when a setting causes a redraw. This is
        # especially useful while editing labels near the bottom of a tall figure.
        vpos = self.preview_scroll.verticalScrollBar().value() if hasattr(self, "preview_scroll") else 0
        hpos = self.preview_scroll.horizontalScrollBar().value() if hasattr(self, "preview_scroll") else 0
        self._capture_manual_legend_position()
        self._sync_to_project()
        self.fig.clear()
        self._legend_artists = []
        self._annotation_entries = []
        self._drag_info = None
        mode = self.mode_combo.currentIndex()
        effective_h = self._effective_stacked_height(len(selected), mode)
        self._update_effective_height_label(effective_h, mode)
        self.fig.set_size_inches(float(self.width_spin.value()), effective_h, forward=True)
        self._sync_preview_canvas_size()
        self._last_axes = []
        if not selected:
            ax = self.fig.add_subplot(111)
            ax.text(0.5, 0.5, "Wybierz co najmniej jeden dyfraktogram", ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off(); self.canvas.draw_idle(); return

        lw = float(self.linewidth_spin.value())
        grid = self.grid_cb.isChecked()
        trace_lines = []
        separator_entries = []

        if mode == 2:  # separate panels
            axes = self.fig.subplots(len(selected), 1, sharex=True, squeeze=False).ravel()
            self._last_axes = list(axes)
            pending = []
            for ax, info in zip(axes, selected):
                ds = info["dataset"]; y = self._curve(ds, info)
                line, = ax.plot(ds.x, y, lw=lw, color=info["color"], label=info["label"])
                trace_lines.append(line)
                if self.peaks_cb.isChecked() and ds.peaks:
                    valid = [p for p in ds.peaks if 0 <= p.index < len(y)]
                    ax.scatter([p.two_theta for p in valid], [y[p.index] for p in valid], s=12, color=info["color"], zorder=4)
                ax.set_ylabel(info["label"] if self.sample_names_cb.isChecked() else "",
                              fontsize=self.sample_font_spin.value(), fontfamily=self.font_family_combo.currentText())
                self._apply_axis_typography(ax)
                ax.grid(alpha=0.18 if grid else 0)
                self._add_legend(ax, ds)
                pending.append((ax, ds, y, 0.0, info["color"]))
            for ax, *_ in pending:
                self._expand_for_peak_labels(ax)
            for args in pending:
                self._draw_peak_labels(*args)
            axes[-1].set_xlabel(self.xlabel_edit.text(), fontsize=self.axis_font_spin.value(), fontfamily=self.font_family_combo.currentText())
            self.fig.suptitle(self._figure_title_text(), fontsize=self.title_font_spin.value(), fontfamily=self.font_family_combo.currentText())
            self.fig.text(0.01, 0.5, self.ylabel_edit.text(), rotation=90, va="center",
                          fontsize=self.axis_font_spin.value(), fontfamily=self.font_family_combo.currentText())
        else:
            ax = self.fig.add_subplot(111)
            self._last_axes = [ax]
            pending = []
            if mode == 0:  # overlay
                for info in selected:
                    ds = info["dataset"]; y = self._curve(ds, info)
                    line, = ax.plot(ds.x, y, lw=lw, color=info["color"], label=info["label"])
                    trace_lines.append(line)
                    if self.peaks_cb.isChecked() and ds.peaks:
                        valid = [p for p in ds.peaks if 0 <= p.index < len(y)]
                        ax.scatter([p.two_theta for p in valid], [y[p.index] for p in valid], s=12, color=info["color"], zorder=4)
                    pending.append((ax, ds, y, 0.0, info["color"]))
                ax.set_ylabel(self.ylabel_edit.text(), fontsize=self.axis_font_spin.value(), fontfamily=self.font_family_combo.currentText())
            else:  # stacked
                spacing = float(self.offset_spin.value())
                for i, info in enumerate(selected):
                    ds = info["dataset"]; y = self._curve(ds, info)
                    lo = float(np.nanmin(y)); hi = float(np.nanmax(y))
                    yn = (y - lo) / (hi - lo) if hi > lo else np.zeros_like(y)
                    offset = (len(selected) - 1 - i) * spacing
                    line, = ax.plot(ds.x, yn + offset, lw=lw, color=info["color"], label=info["label"])
                    trace_lines.append(line)
                    if self.peaks_cb.isChecked() and ds.peaks:
                        valid = [p for p in ds.peaks if 0 <= p.index < len(yn)]
                        ax.scatter([p.two_theta for p in valid], [yn[p.index] + offset for p in valid], s=12, color=info["color"], zorder=4)
                    if self.sample_names_cb.isChecked():
                        ax.text(float(ds.x.min()) + 0.01 * (float(ds.x.max()) - float(ds.x.min())), offset + 0.43,
                                info["label"], fontsize=self.sample_font_spin.value(), va="center",
                                fontfamily=self.font_family_combo.currentText())
                    band_low = offset - 0.02
                    # Every pattern gets its own annotation band. For the TOP pattern
                    # there is no next baseline, so do NOT let the band grow with
                    # Vertical spacing; otherwise a large spacing creates a large empty
                    # area above the whole figure. Lower patterns may use the actual gap.
                    if i == 0:
                        # Reserve a finite editing band for the top pattern.  Vertical text needs
                        # more room than horizontal text, but the margin stays independent of spacing.
                        extra_top = 0.58 if (self.peak_labels_cb.isChecked() and self.peak_label_orientation_combo.currentText()=="Pionowo") else 0.38
                        band_high = offset + 1.0 + extra_top
                    else:
                        band_high = offset + max(1.04, spacing - 0.04)
                    pending.append((ax, ds, yn, offset, info["color"], band_low, band_high))
                top_offset = (len(selected) - 1) * spacing
                if self.peak_labels_cb.isChecked():
                    if self.peak_label_position_combo.currentText()=="Nad pikiem":
                        top_pad = 0.62 if self.peak_label_orientation_combo.currentText()=="Pionowo" else 0.42
                    else:
                        top_pad = 0.24
                else:
                    top_pad = 0.12
                # Publication marker strip lives below the lowest trace and therefore gets its
                # own fixed data-space margin; it never scales with Vertical spacing.
                ax.set_xlim(float(min(info['dataset'].x.min() for info in selected)), float(max(info['dataset'].x.max() for info in selected)))
                marker_depth = self._draw_phase_marker_strip(ax, selected, spacing)
                bottom_pad = max(0.08, marker_depth)
                ax.set_ylim(-bottom_pad, top_offset + 1.0 + top_pad)
                # Wiele linii podziału. Napisy są odsuwane w punktach typograficznych,
                # a nie w jednostkach danych, dzięki czemu zachowują stały odstęp przy zmianie DPI.
                if self.separators and len(selected) > 1:
                    prepared = []
                    for sep0 in self.separators:
                        sep = normalizuj_separator(sep0, len(self.rows))
                        pos = sep["position"]
                        br = sum(1 for r in self.rows[:pos] if r.get("enabled", True))
                        if br <= 0 or br >= len(selected):
                            continue
                        sep_y = (len(selected) - br - 0.5) * spacing
                        sep_y += (float(sep.get("offset_pct", 0.0)) / 100.0) * spacing
                        prepared.append([sep, float(sep_y)])

                    # Jeżeli dwa separatory trafią w dokładnie tę samą granicę, zachowaj oba
                    # i rozsuń je tylko wizualnie o kilka punktów — bez zmiany danych projektu.
                    groups = {}
                    for idx, (_sep, yy) in enumerate(prepared):
                        groups.setdefault(round(yy, 9), []).append(idx)
                    for indices in groups.values():
                        if len(indices) <= 1:
                            continue
                        center = (len(indices) - 1) / 2.0
                        for j, idx in enumerate(indices):
                            prepared[idx][1] += self._data_dy_from_points(ax, (j - center) * 8.0)

                    for sep, sep_y in prepared:
                        separator_entries.append(self._draw_separator(ax, sep, sep_y))
                ax.set_ylabel(self.ylabel_edit.text(), fontsize=self.axis_font_spin.value(), fontfamily=self.font_family_combo.currentText())
                ax.set_yticks([])
            if mode != 1:
                self._expand_for_peak_labels(ax)
            for args in pending:
                self._draw_peak_labels(*args)
            ax.set_title(self._figure_title_text(), fontsize=self.title_font_spin.value(), fontfamily=self.font_family_combo.currentText())
            ax.set_xlabel(self.xlabel_edit.text(), fontsize=self.axis_font_spin.value(), fontfamily=self.font_family_combo.currentText())
            self._apply_axis_typography(ax)
            ax.grid(alpha=0.18 if grid else 0)
            self._add_legend(ax)

        try:
            self.fig.tight_layout(rect=[0.025, 0.025, 0.995, 0.965] if mode == 2 else None)
        except Exception:
            pass
        # Pierwsze pełne rysowanie ustala rzeczywiste rozmiary tekstów po tight_layout.
        # Dopiero wtedy można wiarygodnie wykrywać kolizje separatorów z etykietami pików.
        try:
            self.canvas.draw()
            self._resolve_separator_collisions(separator_entries, trace_lines)
        except Exception:
            pass
        self.canvas.draw_idle()
        if hasattr(self, "preview_scroll"):
            # Clamp automatically if the new canvas became smaller.
            vs = self.preview_scroll.verticalScrollBar()
            hs = self.preview_scroll.horizontalScrollBar()
            vs.setValue(min(vpos, vs.maximum()))
            hs.setValue(min(hpos, hs.maximum()))

    def _on_scroll(self, event):
        ax = event.inaxes
        if ax is None:
            return
        base_scale = 1.25
        scale = 1.0 / base_scale if event.button == "up" else base_scale
        x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
        xc = event.xdata if event.xdata is not None else (x0 + x1) / 2.0
        yc = event.ydata if event.ydata is not None else (y0 + y1) / 2.0
        new_w = (x1 - x0) * scale; new_h = (y1 - y0) * scale
        rx = (xc - x0) / (x1 - x0) if x1 != x0 else 0.5
        ry = (yc - y0) / (y1 - y0) if y1 != y0 else 0.5
        ax.set_xlim(xc - new_w * rx, xc + new_w * (1 - rx))
        ax.set_ylim(yc - new_h * ry, yc + new_h * (1 - ry))
        self.canvas.draw_idle()

    @staticmethod
    def _left_button(event) -> bool:
        b = getattr(event, "button", None)
        return getattr(b, "value", b) == 1

    def _on_mouse_press(self, event):
        if getattr(event, "dblclick", False):
            self._on_double_click(event)
            return
        if not self._left_button(event) or event.inaxes is None:
            return
        # Prefer the visually topmost annotation when labels overlap.
        for entry in reversed(self._annotation_entries):
            ann = entry.get("annotation")
            if entry.get("ax") is not event.inaxes or ann is None:
                continue
            try:
                hit, _ = ann.contains(event)
            except Exception:
                hit = False
            if hit:
                self._drag_info = entry
                try:
                    ann.set_alpha(0.65)
                except Exception:
                    pass
                self.canvas.draw_idle()
                return

    def _on_mouse_motion(self, event):
        info = self._drag_info
        if not info or event.inaxes is not info.get("ax") or event.xdata is None:
            return
        ann = info.get("annotation")
        ax = info.get("ax")
        if ann is None or ax is None:
            return
        x = float(event.xdata)
        mode = info.get("coord_mode")
        if mode == "data":
            if event.ydata is None:
                return
            y = float(event.ydata)
            lo = info.get("band_low"); hi = info.get("band_high")
            if lo is not None and hi is not None:
                margin = max(0.01, 0.025 * max(1e-6, float(hi)-float(lo)))
                y = min(max(y, float(lo)+margin), float(hi)-margin)
            ann.xyann = (x, y)
        else:
            try:
                yfrac = float(ax.transAxes.inverted().transform((event.x, event.y))[1])
            except Exception:
                return
            yfrac = min(max(yfrac, 0.04), 0.97)
            ann.xyann = (x, yfrac)
        self.canvas.draw_idle()

    def _on_mouse_release(self, event):
        info = self._drag_info
        if not info:
            return
        ann = info.get("annotation")
        try:
            if ann is not None:
                pos = ann.xyann
                self._label_positions[str(info.get("key"))] = {
                    "x": float(pos[0]), "y": float(pos[1]),
                    "coord": str(info.get("coord_mode")),
                }
                ann.set_alpha(1.0)
        except Exception:
            pass
        self._drag_info = None
        self.canvas.draw_idle()
        if self.parent_window is not None:
            try:
                self.parent_window.plot_creator_state = self.get_state()
                self.parent_window.mark_project_dirty()
            except Exception:
                pass

    def _on_double_click(self, event):
        if not getattr(event, "dblclick", False) or event.inaxes is None or event.xdata is None:
            return
        candidates = []
        xspan = abs(event.inaxes.get_xlim()[1] - event.inaxes.get_xlim()[0])
        max_dx = max(0.10, xspan * 0.008)
        for info in self._selected_rows():
            ds = info["dataset"]
            for p in ds.peaks:
                dx = abs(p.two_theta - float(event.xdata))
                if dx <= max_dx:
                    candidates.append((dx, ds, p))
        if not candidates:
            return
        _, ds, p = min(candidates, key=lambda t: t[0])
        current = p.custom_label or self._peak_text(ds, p)
        text, ok = QInputDialog.getText(self, "Etykieta piku", f"Etykieta dla {ds.name}, 2θ = {p.two_theta:.3f}°:", text=current)
        if ok:
            p.custom_label = text.strip()
            self.redraw()
            if self.parent_window is not None:
                try:
                    self.parent_window.update_table()
                    self.parent_window.redraw()
                    self.parent_window.mark_project_dirty()
                except Exception:
                    pass

    def save_figure(self):
        if not self._selected_rows():
            QMessageBox.information(self, "Zapis figury", "Najpierw wybierz co najmniej jeden dyfraktogram.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Zapisz figurę jakości publikacyjnej", "XRD_multi_pattern.png",
            "PNG (*.png);;TIFF (*.tif *.tiff);;SVG vector (*.svg);;PDF vector (*.pdf);;EPS vector (*.eps)"
        )
        if not path:
            return
        ext = Path(path).suffix.lower()
        if not ext:
            path += ".png"; ext = ".png"
        dpi = int(self.dpi_combo.currentText())
        try:
            self.redraw()
            kwargs = {"bbox_inches": "tight"}
            if ext in (".png", ".tif", ".tiff"):
                kwargs["dpi"] = dpi
            self.fig.savefig(path, **kwargs)
            QMessageBox.information(self, "Zapis figury", f"Zapisano:\n{path}")
        except Exception as e:
            QMessageBox.critical(self, "Zapis figury", f"Nie udało się zapisać figury:\n{e}")


class InstrumentMeasurementDialog(QDialog):
    """Edit measurement metadata and reusable instrument profiles.

    File-derived scan start/end/step are shown read-only. User-supplied values are
    stored together with provenance so project files preserve what came from the
    source file and what was entered manually or via a saved profile.
    """
    FIELD_SPECS = [
        ("instrument", "Aparat", "text"),
        ("manufacturer", "Producent", "text"),
        ("radiation", "Promieniowanie", "text"),
        ("anode", "Anoda", "text"),
        ("wavelength_A", "Długość fali λ [Å]", "float"),
        ("tube_voltage_kV", "Napięcie lampy [kV]", "float"),
        ("tube_current_mA", "Prąd lampy [mA]", "float"),
        ("scan_rate_deg_s", "Szybkość skanowania [°/s]", "float"),
        ("geometry", "Geometria", "text"),
        ("notes", "Uwagi", "text"),
    ]

    def __init__(self, dataset: Dataset, all_datasets: List[Dataset], profiles: InstrumentProfileStore, parent=None):
        super().__init__(parent)
        self.dataset = dataset
        self.all_datasets = all_datasets
        self.profiles = profiles
        self.setWindowTitle("Aparatura / parametry pomiaru")
        self.setWindowFlags(self.windowFlags() | Qt.WindowType.WindowMaximizeButtonHint | Qt.WindowType.WindowMinimizeButtonHint)
        self.resize(980, 760)

        root = QVBoxLayout(self)

        # Profile controls
        pg = QGroupBox("Profil aparatury")
        pl = QGridLayout(pg)
        self.profile_combo = QComboBox()
        self.profile_combo.addItems(self.profiles.names())
        selected = str(dataset.metadata.get("instrument_profile", ""))
        if selected and self.profile_combo.findText(selected) >= 0:
            self.profile_combo.setCurrentText(selected)
        self.profile_load_btn = QPushButton("Wczytaj profil do pól")
        self.profile_save_btn = QPushButton("Zapisz pola jako profil…")
        self.profile_delete_btn = QPushButton("Usuń profil")
        self.profile_reset_btn = QPushButton("Przywróć wbudowany profil domyślny")
        self.apply_current_btn = QPushButton("Zastosuj profil do bieżącej próbki")
        self.apply_all_btn = QPushButton("Zastosuj profil do wszystkich próbek z brakującymi metadanymi")
        self.fill_missing_cb = QCheckBox("Uzupełniaj tylko brakujące wartości (zachowaj metadane odczytane z plików)")
        self.fill_missing_cb.setChecked(True)
        pl.addWidget(QLabel("Profil:"), 0, 0)
        pl.addWidget(self.profile_combo, 0, 1, 1, 3)
        pl.addWidget(self.profile_load_btn, 1, 0)
        pl.addWidget(self.profile_save_btn, 1, 1)
        pl.addWidget(self.profile_delete_btn, 1, 2)
        pl.addWidget(self.profile_reset_btn, 1, 3)
        pl.addWidget(self.fill_missing_cb, 2, 0, 1, 4)
        pl.addWidget(self.apply_current_btn, 3, 0, 1, 2)
        pl.addWidget(self.apply_all_btn, 3, 2, 1, 2)
        root.addWidget(pg)

        # Editable metadata
        mg = QGroupBox(f"Metadane pomiaru – {dataset.name}")
        grid = QGridLayout(mg)
        grid.addWidget(QLabel("Parametr"), 0, 0)
        grid.addWidget(QLabel("Wartość"), 0, 1)
        grid.addWidget(QLabel("Źródło / pochodzenie"), 0, 2)
        self.fields: Dict[str, QWidget] = {}
        self.source_labels: Dict[str, QLabel] = {}
        row = 1
        for key, label, typ in self.FIELD_SPECS:
            grid.addWidget(QLabel(label), row, 0)
            if typ == "float":
                w = QDoubleSpinBox()
                w.setDecimals(6 if key == "wavelength_A" else 4)
                w.setRange(0.0, 100000.0)
                w.setSpecialValueText("nie ustawiono")
                w.setKeyboardTracking(False)
            else:
                w = QLineEdit()
            self.fields[key] = w
            sl = QLabel("")
            sl.setWordWrap(True)
            self.source_labels[key] = sl
            grid.addWidget(w, row, 1)
            grid.addWidget(sl, row, 2)
            row += 1

        # Read-only values derived directly from diffraction axis.
        for key, label, value in [
            ("scan_start_deg", "Początek skanu 2θ [°]", dataset.scan_start_deg),
            ("scan_end_deg", "Koniec skanu 2θ [°]", dataset.scan_end_deg),
            ("scan_step_deg", "Krok skanu [°]", dataset.scan_step_deg),
            ("points", "Liczba punktów", len(dataset.x)),
            ("max_intensity", "Maksymalna zmierzona intensywność", dataset.max_total_intensity),
        ]:
            grid.addWidget(QLabel(label), row, 0)
            val = QLineEdit(f"{value:.6g}" if isinstance(value, float) else str(value))
            val.setReadOnly(True)
            src = QLabel("wyznaczone z danych dyfrakcyjnych")
            grid.addWidget(val, row, 1)
            grid.addWidget(src, row, 2)
            row += 1

        self.reload_sample_btn = QPushButton("Przywróć pola z metadanych bieżącej próbki")
        self.reload_sample_btn.setToolTip("Odrzuć niezapisane zmiany w tym oknie i ponownie wczytaj wartości zapisane dla bieżącej próbki")
        grid.addWidget(self.reload_sample_btn, row, 0, 1, 3)
        root.addWidget(mg, 1)

        note = QLabel(
            "Pliki ASCII/ASC zwykle zawierają jedynie 2θ i intensywność. "
            "Dla długości fali, napięcia/prądu lampy i szybkości skanowania użyj zweryfikowanego profilu aparatury. "
            "Wartości odczytane z XRDML/UDF/CIF są domyślnie zachowywane i oznaczane jako metadane pochodzące z pliku."
        )
        note.setWordWrap(True)
        root.addWidget(note)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Save | QDialogButtonBox.StandardButton.Close)
        buttons.button(QDialogButtonBox.StandardButton.Save).setText("Zapisz metadane bieżącej próbki")
        buttons.button(QDialogButtonBox.StandardButton.Close).setText("Zamknij")
        buttons.accepted.connect(self.save_current_metadata)
        buttons.rejected.connect(self.reject)
        root.addWidget(buttons)

        self.profile_load_btn.clicked.connect(self.load_profile_fields)
        self.profile_save_btn.clicked.connect(self.save_profile_from_fields)
        self.profile_delete_btn.clicked.connect(self.delete_profile)
        self.profile_reset_btn.clicked.connect(self.reset_to_builtin_profile)
        self.apply_current_btn.clicked.connect(self.apply_profile_current)
        self.apply_all_btn.clicked.connect(self.apply_profile_all)
        self.reload_sample_btn.clicked.connect(self._load_dataset_fields)
        self._load_dataset_fields()

    def reset_to_builtin_profile(self):
        name = DEFAULT_INSTRUMENT_PROFILE["name"]
        # Ensure the built-in profile exists even if the user profile file was modified.
        if self.profiles.get(name) is None:
            self.profiles.upsert(dict(DEFAULT_INSTRUMENT_PROFILE))
            self.profile_combo.clear()
            self.profile_combo.addItems(self.profiles.names())
        if self.profile_combo.findText(name) >= 0:
            self.profile_combo.setCurrentText(name)
        profile = self.profiles.get(name) or dict(DEFAULT_INSTRUMENT_PROFILE)
        for key, _, _ in self.FIELD_SPECS:
            if key in profile:
                self._set_widget_value(key, profile.get(key))
                self.source_labels[key].setText("built-in default profile")
        self.fill_missing_cb.setChecked(True)

    def _set_widget_value(self, key: str, value):
        w = self.fields[key]
        if isinstance(w, QDoubleSpinBox):
            try:
                w.setValue(float(value))
            except Exception:
                w.setValue(0.0)
        else:
            w.setText("" if value is None else str(value))

    def _widget_value(self, key: str):
        w = self.fields[key]
        if isinstance(w, QDoubleSpinBox):
            v = float(w.value())
            return None if v <= 0 else v
        return w.text().strip()

    def _load_dataset_fields(self):
        ds = self.dataset
        for key, _, _ in self.FIELD_SPECS:
            if key == "wavelength_A":
                value = ds.wavelength_A
            else:
                value = ds.metadata.get(key, "")
            self._set_widget_value(key, value)
            self.source_labels[key].setText(metadata_source(ds, key))

    def _selected_profile(self) -> Optional[dict]:
        return self.profiles.get(self.profile_combo.currentText())

    def load_profile_fields(self):
        profile = self._selected_profile()
        if not profile:
            return
        for key, _, _ in self.FIELD_SPECS:
            self._set_widget_value(key, profile.get(key, ""))

    def _fields_profile_dict(self, name: str) -> dict:
        rec = {"name": name}
        for key, _, _ in self.FIELD_SPECS:
            rec[key] = self._widget_value(key)
        return rec

    def save_profile_from_fields(self):
        current = self.profile_combo.currentText().strip() or "Profil aparatury"
        name, ok = QInputDialog.getText(self, "Zapisz profil aparatury", "Profile name:", text=current)
        if not ok or not name.strip():
            return
        try:
            self.profiles.upsert(self._fields_profile_dict(name.strip()))
            self.profile_combo.clear(); self.profile_combo.addItems(self.profiles.names())
            self.profile_combo.setCurrentText(name.strip())
        except Exception as e:
            QMessageBox.critical(self, "Profil aparatury", f"Nie udało się zapisać profilu:\n{e}")

    def delete_profile(self):
        name = self.profile_combo.currentText().strip()
        if not name:
            return
        ans = QMessageBox.question(self, "Usuń profil aparatury", f"Usuń profil '{name}'?",
                                   QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
        if ans != QMessageBox.StandardButton.Yes:
            return
        self.profiles.delete(name)
        self.profile_combo.clear(); self.profile_combo.addItems(self.profiles.names())

    def apply_profile_current(self):
        profile = self._selected_profile()
        if not profile:
            return
        try:
            invalidated = apply_instrument_profile(self.dataset, profile, self.fill_missing_cb.isChecked())
        except ValueError as exc:
            QMessageBox.warning(self, "Długość fali", str(exc))
            return
        self._load_dataset_fields()
        if self.parent() is not None:
            try:
                self.parent().refresh_current_metadata_view()
                self.parent().mark_project_dirty()
            except Exception:
                pass

        if invalidated:
            QMessageBox.information(self, "Zmiana długości fali",
                "Przeliczono odległości d dla nowej długości fali. Poprzednie przypisania faz "
                "i zależne wyniki zostały unieważnione. Uruchom ponowną identyfikację faz.")

    def apply_profile_all(self):
        profile = self._selected_profile()
        if not profile:
            return
        count = 0
        invalidated_count = 0
        for ds in self.all_datasets:
            before = dict(ds.metadata_sources)
            try:
                invalidated_count += int(apply_instrument_profile(ds, profile, True))
            except ValueError as exc:
                QMessageBox.warning(self, "Długość fali", f"{ds.name}: {exc}")
                continue
            if ds.metadata_sources != before:
                count += 1
        if self.dataset in self.all_datasets:
            self._load_dataset_fields()
        if self.parent() is not None:
            try:
                self.parent().refresh_current_metadata_view()
                self.parent().mark_project_dirty()
            except Exception:
                pass
        message = f"Profil zastosowano do {count} wczytanych próbek z brakującymi metadanymi."
        if invalidated_count:
            message += (f"\nPrzeliczono d i unieważniono poprzednie wyniki fazowe dla {invalidated_count} próbek. "
                        "Uruchom ponowną identyfikację faz.")
        QMessageBox.information(self, "Profil aparatury", message)

    def save_current_metadata(self):
        ds = self.dataset
        wavelength = self._widget_value("wavelength_A")
        invalidated = False
        if wavelength is not None:
            try:
                invalidated = set_dataset_wavelength(ds, wavelength, "user")
            except ValueError as exc:
                QMessageBox.warning(self, "Długość fali", str(exc))
                return
        for key, _, _ in self.FIELD_SPECS:
            value = self._widget_value(key)
            if key == "wavelength_A":
                continue
            if value in (None, ""):
                ds.metadata.pop(key, None)
                ds.metadata_sources.pop(key, None)
            else:
                ds.metadata[key] = value
                ds.metadata_sources[key] = "user"
        self._load_dataset_fields()
        if self.parent() is not None:
            try:
                self.parent().refresh_current_metadata_view()
                self.parent().mark_project_dirty()
            except Exception:
                pass
        message = "Metadane bieżącej próbki zapisano w stanie projektu."
        if invalidated:
            message += ("\nPrzeliczono d dla nowej długości fali i unieważniono poprzednie wyniki fazowe. "
                        "Uruchom ponowną identyfikację faz.")
        QMessageBox.information(self, "Metadane pomiaru", message)

# --------------------------- main GUI ---------------------------
# --------------------------- small UI helpers ---------------------------

NUMERIC_PEAK_COLUMNS = (2, 3, 4, 5, 6, 7, 8, 12)  # numeric columns of the peak table


class NumericSortItem(QTableWidgetItem):
    """Table cell that sorts by its numeric value instead of alphabetically ("9.5" < "10.2")."""

    @staticmethod
    def _value(item):
        try:
            return float(item.text().replace(",", "."))
        except (ValueError, AttributeError):
            return None

    def __lt__(self, other):
        a, b = self._value(self), self._value(other)
        if a is None or b is None:
            if a is None and b is None:
                return super().__lt__(other)
            return b is None  # numbers sort before empty cells
        return a < b


def _number_or_text(value):
    """Return a float for numeric text (also with decimal comma), otherwise the text unchanged."""
    try:
        number = float(str(value).replace(",", "."))
    except ValueError:
        return value
    return number if math.isfinite(number) else value


def _safe_float(value, default: float = 0.0) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    return number if math.isfinite(number) else default


def _safe_filename(name: str, fallback: str = "XRD") -> str:
    """Make a sample name usable as a file name (no path separators or reserved characters)."""
    cleaned = re.sub(r'[\\/:*?"<>|\r\n\t]+', "_", str(name or "")).strip(" ._")
    return cleaned or fallback


def _display_text(value, missing: str = "nie ustawiono") -> str:
    text = "" if value is None else str(value).strip()
    return missing if text in ("", "not set") else text


def _display_number(value, unit: str, missing: str = "nie ustawiono") -> str:
    if value in (None, "", "not set"):
        return missing
    try:
        return f"{float(value):g} {unit}"
    except (TypeError, ValueError):
        return str(value)


_IN_EXCEPTION_HOOK = False


def install_exception_hook():
    """Report unexpected exceptions instead of letting PyQt6 abort the whole program.

    Since PyQt 5.5 an unhandled Python exception inside a Qt slot terminates the process
    unless sys.excepthook is replaced. The hook appends the traceback to
    ~/XRDAnalyzerData/XRD_Analyzer_error.log and shows a dialog; the application keeps running.
    """
    def hook(exc_type, exc, tb):
        global _IN_EXCEPTION_HOOK
        if issubclass(exc_type, KeyboardInterrupt) or _IN_EXCEPTION_HOOK:
            sys.__excepthook__(exc_type, exc, tb)
            return
        _IN_EXCEPTION_HOOK = True
        try:
            text = "".join(traceback.format_exception(exc_type, exc, tb))
            try:
                sys.__stderr__.write(text)
            except Exception:
                pass
            log_note = ""
            try:
                log_dir = Path.home() / "XRDAnalyzerData"
                log_dir.mkdir(parents=True, exist_ok=True)
                log_path = log_dir / "XRD_Analyzer_error.log"
                with log_path.open("a", encoding="utf-8") as stream:
                    stream.write(f"\n=== {datetime.now():%Y-%m-%d %H:%M:%S} | XRD Analyzer {APP_VERSION} ===\n{text}")
                log_note = f"\n\nSzczegóły zapisano w pliku:\n{log_path}"
            except Exception:
                pass
            try:
                QMessageBox.critical(
                    QApplication.activeWindow(), "Nieoczekiwany błąd",
                    f"{exc_type.__name__}: {exc}\n\n"
                    "Program nadal działa, ale ostatnia operacja nie została dokończona. "
                    "Zapisz projekt i, jeśli błąd się powtarza, prześlij plik dziennika autorowi." + log_note)
            except Exception:
                pass
        finally:
            _IN_EXCEPTION_HOOK = False

    sys.excepthook = hook


class MainWindow(QMainWindow):
    PROJECT_VERSION = 13

    def __init__(self):
        super().__init__()
        self.setWindowTitle("XRD Analyzer – analiza wsadowa, wyszukiwanie COD i kreator wykresów")
        self.resize(1500, 900)
        self.db = PhaseDatabase()
        self.settings = QSettings(QSettings.defaultFormat(), QSettings.Scope.UserScope,
                                  "OpenAI-XRDTools", "XRD Analyzer")
        self.instrument_profiles = InstrumentProfileStore()
        self.cod_db = CODDatabase()
        saved_cod = self.settings.value("cod_db_path", "")
        if saved_cod and Path(str(saved_cod)).exists():
            try:
                self.cod_db.open(Path(str(saved_cod)))
            except Exception:
                pass
        self.datasets: List[Dataset] = []
        self.project_path: Optional[Path] = None
        self.project_dirty = False
        self._loading_project = False
        self.plot_creator_state: dict = {}
        self._plot_creator_dialog = None
        self._last_save_error = None

        self._build_menus()
        self.statusBar().showMessage("Gotowy. Wczytaj pliki XRD (Ctrl+Shift+O) lub otwórz projekt (Ctrl+O).")

        central = QWidget(); self.setCentralWidget(central)
        outer = QVBoxLayout(central)

        top = QGridLayout()
        self.open_project_btn = QPushButton("Otwórz projekt")
        self.save_project_btn = QPushButton("Zapisz projekt")
        self.save_project_as_btn = QPushButton("Zapisz projekt jako…")
        self.load_btn = QPushButton("Otwórz jeden lub wiele plików XRD")
        self.analyze_btn = QPushButton("Analizuj bieżący")
        self.analyze_all_btn = QPushButton("Analizuj wszystkie wczytane")
        self.batch_chemistry_btn = QPushButton("Chemia serii…")
        self.batch_chemistry_btn.setToolTip("Ustaw pierwiastki wspólne dla wybranej serii oraz dodatkowe pierwiastki i wskazówki osnowy dla poszczególnych próbek. Nazwy plików nie są używane jako informacja chemiczna.")
        self.plot_creator_btn = QPushButton("Kreator wykresów")
        self.save_plot_btn = QPushButton("Zapisz bieżący wykres")
        self.export_btn = QPushButton("Eksportuj tabelę")
        self.db_btn = QPushButton("Bazy faz / COD")
        self.cod_pattern_btn = QPushButton("Wyszukaj COD dla zmierzonego dyfraktogramu")
        self.ai_review_btn = QPushButton("Ocena faz wspomagana AI…")
        self.ai_review_btn.setToolTip("Doradcze porównanie kandydatów z COD, bazy wewnętrznej i dopasowania samej struktury. Program nie zatwierdza faz automatycznie.")
        self.instrument_btn = QPushButton("Aparatura / pomiar")
        self.close_program_btn = QPushButton("Zamknij program")
        for tip_button, tip_text in (
            (self.open_project_btn, "Otwórz zapisany projekt (.xrdproj).  [Ctrl+O]"),
            (self.save_project_btn, "Zapisz bieżący projekt.  [Ctrl+S]"),
            (self.save_project_as_btn, "Zapisz projekt pod nową nazwą.  [Ctrl+Shift+S]"),
            (self.load_btn, "Wczytaj jeden lub wiele dyfraktogramów (CIF, UDF, RD, XRDML, CSV, TXT, DAT, ASC).  [Ctrl+Shift+O]"),
            (self.instrument_btn, "Metadane aparatury i pomiaru: długość fali, lampa, geometria, profile aparatury."),
            (self.analyze_btn, "Wykryj piki i zidentyfikuj fazy dla bieżącej próbki.  [F5]"),
            (self.analyze_all_btn, "Przeanalizuj wszystkie wczytane próbki z bieżącymi ustawieniami.  [Shift+F5]"),
            (self.db_btn, "Zarządzaj wewnętrzną bazą faz i zewnętrzną bazą COD."),
            (self.cod_pattern_btn, "Wyszukaj w bazie COD struktury pasujące do zmierzonych pików."),
            (self.plot_creator_btn, "Wykres porównawczy wielu próbek z kolejnością, wygładzaniem i eksportem.  [Ctrl+K]"),
            (self.save_plot_btn, "Zapisz bieżący wykres jako PNG, TIFF, PDF, SVG lub EPS.  [Ctrl+E]"),
            (self.export_btn, "Zapisz tabelę pików do pliku Excel lub CSV.  [Ctrl+T]"),
            (self.close_program_btn, "Zamknij program (zapyta o niezapisane zmiany).  [Ctrl+Q]"),
        ):
            if not tip_button.toolTip():
                tip_button.setToolTip(tip_text)
        for button_index, b in enumerate((
            # row 1: project and data
            self.open_project_btn, self.save_project_btn, self.save_project_as_btn,
            self.load_btn, self.instrument_btn,
            # row 2: analysis and databases
            self.analyze_btn, self.analyze_all_btn, self.batch_chemistry_btn,
            self.db_btn, self.cod_pattern_btn,
            # row 3: results and output
            self.ai_review_btn, self.plot_creator_btn, self.save_plot_btn, self.export_btn,
            self.close_program_btn,
        )):
            top.addWidget(b, button_index // 5, button_index % 5)
        # A single row of translated toolbar buttons otherwise imposes a
        # >3900 px minimum width and defeats every geometry restoration check.
        toolbar_widget = QWidget()
        toolbar_widget.setLayout(top)
        toolbar_scroll = QScrollArea()
        toolbar_scroll.setWidgetResizable(True)
        toolbar_scroll.setWidget(toolbar_widget)
        toolbar_scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        toolbar_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        toolbar_scroll.setFrameShape(QFrame.Shape.NoFrame)
        toolbar_scroll.setFixedHeight(toolbar_widget.sizeHint().height() + 22)
        outer.addWidget(toolbar_scroll)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        outer.addWidget(splitter, 1)

        # Left controls
        left = QWidget(); llay = QVBoxLayout(left)
        left_scroll = QScrollArea()
        left_scroll.setWidgetResizable(True)
        left_scroll.setWidget(left)
        left_scroll.setFrameShape(QFrame.Shape.NoFrame)
        left_scroll.setMinimumWidth(330)
        splitter.addWidget(left_scroll)

        file_group = QGroupBox("Wczytane próbki")
        fl = QVBoxLayout(file_group)
        nav = QHBoxLayout()
        self.prev_btn = QPushButton("◀ Poprzedni XRD")
        self.sample_combo = QComboBox()
        self.next_btn = QPushButton("Następny XRD ▶")
        nav.addWidget(self.prev_btn)
        nav.addWidget(self.sample_combo, 1)
        nav.addWidget(self.next_btn)
        fl.addLayout(nav)
        form = QFormLayout()
        form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)  # narrow panel: label above field
        self.sample_name = QLineEdit(); self.material_type = QLineEdit()
        self.material_type.setPlaceholderText("np. 316L, CuNi30Mn1Fe, Ti-6Al-4V, stopiwo")
        self.allowed_elements_edit = QLineEdit()
        self.allowed_elements_edit.setPlaceholderText("np. Ti Al V Fe O N")
        self.allowed_elements_btn = QPushButton("Układ okresowy…")
        elem_widget = QWidget(); elem_lay = QHBoxLayout(elem_widget); elem_lay.setContentsMargins(0,0,0,0); elem_lay.setSpacing(4)
        elem_lay.addWidget(self.allowed_elements_edit, 1); elem_lay.addWidget(self.allowed_elements_btn)
        self.sample_setup_btn = QPushButton("Kreator konfiguracji próbki…")
        self.sample_setup_btn.setToolTip("Wprowadź chemię i opcjonalnie oczekiwaną osnowę/strukturę (np. γ-austenit FCC dla 316L). Dane te są używane w rankingu faz.")
        self.wavelength = QDoubleSpinBox(); self.wavelength.setDecimals(6); self.wavelength.setRange(0.1, 5.0); self.wavelength.setValue(DEFAULT_WAVELENGTH_A)
        # Emit valueChanged only when editing is finished. Otherwise every intermediate value typed
        # (1, 1.5, 1.54 ...) would recompute d and invalidate saved phase / scientific results.
        self.wavelength.setKeyboardTracking(False); self.wavelength.setSingleStep(0.001)
        self.wavelength.setToolTip("Zmiana λ przelicza odstępy d i unieważnia wyniki identyfikacji faz oraz analizy naukowej "
                                   "tej próbki. Zatwierdź wartość klawiszem Enter lub przejściem do innego pola.")
        self.phase_source_combo = QComboBox()
        self.phase_source_combo.addItems([
            "Tylko baza wewnętrzna",
            "Tylko zewnętrzna baza COD",
            "Baza wewnętrzna + zewnętrzna COD",
            "Bez automatycznej identyfikacji faz",
        ])
        saved_phase_source = _migrate_phase_source(str(self.settings.value("phase_identification_source", "Tylko baza wewnętrzna")))
        if self.phase_source_combo.findText(saved_phase_source) >= 0:
            self.phase_source_combo.setCurrentText(saved_phase_source)
        self.phase_source_combo.setToolTip(
            "Wybierz bazę dopuszczoną do identyfikacji faz. Tryb „Tylko zewnętrzna baza COD” całkowicie wyłącza wewnętrzną bibliotekę faz z automatycznej identyfikacji."
        )
        form.addRow("Nazwa próbki:", self.sample_name)
        form.addRow("Typ materiału:", self.material_type)
        form.addRow("Konfiguracja próbki:", self.sample_setup_btn)
        form.addRow("Dozwolone pierwiastki:", elem_widget)
        form.addRow("Źródło identyfikacji faz:", self.phase_source_combo)
        form.addRow("Długość fali λ [Å]:", self.wavelength)
        self.reset_identification_btn = QPushButton("Przywróć domyślne identyfikacji")
        self.reset_identification_btn.setToolTip("Wyczyść zatwierdzoną chemię, przywróć domyślne źródło bazy oraz długość fali z pliku/domyślną")
        form.addRow(self.reset_identification_btn)
        fl.addLayout(form)
        llay.addWidget(file_group)

        analysis_group = QGroupBox("Detekcja pików")
        af = QFormLayout(analysis_group)
        af.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
        self.baseline_cb = QCheckBox("Korekcja tła (AsLS)"); self.baseline_cb.setChecked(True)
        self.prominence = QDoubleSpinBox(); self.prominence.setRange(0.05, 50); self.prominence.setValue(3.0); self.prominence.setSuffix(" %")
        self.min_distance = QDoubleSpinBox(); self.min_distance.setRange(0.01, 10); self.min_distance.setValue(0.25); self.min_distance.setSuffix(" °")
        self.min_rel_height = QDoubleSpinBox(); self.min_rel_height.setRange(0.0, 100.0); self.min_rel_height.setValue(5.0); self.min_rel_height.setSuffix(" %")
        self.min_snr = QDoubleSpinBox(); self.min_snr.setRange(0.0, 100.0); self.min_snr.setValue(8.0); self.min_snr.setDecimals(1)
        self.min_fwhm = QDoubleSpinBox(); self.min_fwhm.setRange(0.0, 5.0); self.min_fwhm.setValue(0.10); self.min_fwhm.setDecimals(3); self.min_fwhm.setSuffix(" °")
        self.max_peaks = QSpinBox(); self.max_peaks.setRange(5, 500); self.max_peaks.setValue(50)
        self.phase_tol = QDoubleSpinBox(); self.phase_tol.setRange(0.02, 2.0); self.phase_tol.setValue(0.35); self.phase_tol.setSuffix(" °")
        self.min_phase_matches = QSpinBox(); self.min_phase_matches.setRange(2, 8); self.min_phase_matches.setValue(2)
        for spin in (self.prominence, self.min_distance, self.min_rel_height, self.min_snr, self.min_fwhm,
                     self.max_peaks, self.phase_tol, self.min_phase_matches):
            spin.setKeyboardTracking(False)  # no re-save of the project on every typed digit
        af.addRow(self.baseline_cb)
        af.addRow("Min. prominencja:", self.prominence)
        af.addRow("Min. względna wysokość piku:", self.min_rel_height)
        af.addRow("Min. sygnał/szum (SNR):", self.min_snr)
        self.min_fwhm.setToolTip("Szerokość detekcyjna przy połowie prominencji. FWHM modelowe jest dostępne w menu Analiza naukowa.")
        af.addRow("Min. szerokość detekcyjna:", self.min_fwhm)
        af.addRow("Maks. liczba wykrytych pików:", self.max_peaks)
        af.addRow("Min. odległość między pikami:", self.min_distance)
        af.addRow("Tolerancja dopasowania fazy:", self.phase_tol)
        af.addRow("Min. liczba dopasowanych refleksów / fazę:", self.min_phase_matches)
        self.reset_peak_defaults_btn = QPushButton("Przywróć domyślne")
        self.reset_peak_defaults_btn.setToolTip("Przywróć zalecane ustawienia detekcji pików i dopasowania faz")
        af.addRow(self.reset_peak_defaults_btn)
        llay.addWidget(analysis_group)

        plot_group = QGroupBox("Etykiety wykresu i eksport")
        pf = QFormLayout(plot_group)
        pf.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
        self.title_edit = QLineEdit("Dyfraktogram XRD – profil pomiarowy i prawdopodobne fazy")
        self.xlabel_edit = QLineEdit("2θ [°]")
        self.ylabel_edit = QLineEdit("Intensywność [zliczenia]")
        self.annotate_cb = QCheckBox("Opisuj wykryte piki"); self.annotate_cb.setChecked(True)
        self.label_mode_combo = QComboBox(); self.label_mode_combo.addItems([
            "Tylko dopasowane fazy", "Top N najsilniejszych pików", "Wszystkie wykryte piki"
        ])
        self.peak_label_style_combo = QComboBox(); self.peak_label_style_combo.addItems(PEAK_LABEL_STYLES); self.peak_label_style_combo.setCurrentText(PEAK_LABEL_STYLE_FULL)
        self.peak_label_style_combo.setToolTip("Publikacyjny – krótki: γ(111), Cu–Ni(200), α-Ti(101). Publikacyjny – pełny: γ-austenite FCC (111). Baza / diagnostyka: pełna faza + ID COD + hkl + 2θ.")
        self.hkl_display_combo = QComboBox(); self.hkl_display_combo.addItems(HKL_DISPLAY_MODES); self.hkl_display_combo.setCurrentText(HKL_DISPLAY_AUTO)
        self.hkl_display_combo.setToolTip("Wybierz, czy indeksy Millera mają być ukryte, pokazane z przypisanego źródła fazy, czy ograniczone do bazy wewnętrznej/COD.")
        self.label_top_n = QSpinBox(); self.label_top_n.setRange(1, 200); self.label_top_n.setValue(12)
        self.peak_font_spin = QSpinBox(); self.peak_font_spin.setRange(5, 24); self.peak_font_spin.setValue(7); self.peak_font_spin.setSuffix(" pt")
        for spin in (self.label_top_n, self.peak_font_spin):
            spin.setKeyboardTracking(False)  # redraw once per edit, not per typed digit
        self.show_baseline_cb = QCheckBox("Pokaż obliczone tło")
        self.dpi_combo = QComboBox(); self.dpi_combo.addItems(["300", "600", "1200"]); self.dpi_combo.setCurrentText("600")
        pf.addRow("Tytuł:", self.title_edit)
        pf.addRow("Oś X:", self.xlabel_edit)
        pf.addRow("Oś Y:", self.ylabel_edit)
        pf.addRow(self.annotate_cb)
        pf.addRow("Które piki etykietować:", self.label_mode_combo)
        pf.addRow("Styl etykiet pików:", self.peak_label_style_combo)
        pf.addRow("Indeksy Millera (hkl):", self.hkl_display_combo)
        pf.addRow("Liczba etykiet Top N:", self.label_top_n)
        pf.addRow("Czcionka etykiet pików:", self.peak_font_spin)
        pf.addRow(self.show_baseline_cb)
        pf.addRow("DPI PNG/TIFF:", self.dpi_combo)
        self.reset_plot_defaults_btn = QPushButton("Przywróć domyślne")
        self.reset_plot_defaults_btn.setToolTip("Przywróć domyślny tytuł, opisy osi, adnotacje pików i DPI eksportu")
        self.reset_view_btn = QPushButton("Przywróć widok wykresu")
        self.reset_view_btn.setToolTip("Przywróć pełny automatycznie wyskalowany widok głównego wykresu po powiększaniu kółkiem myszy")
        pf.addRow(self.reset_plot_defaults_btn)
        pf.addRow(self.reset_view_btn)
        llay.addWidget(plot_group)

        for combo in (self.sample_combo, self.phase_source_combo, self.label_mode_combo,
                      self.peak_label_style_combo, self.hkl_display_combo):
            # long entries must not force a horizontal scroll bar in the narrow left panel
            combo.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
            combo.setMinimumContentsLength(12)
        help_label = QLabel("Główny wykres: kółko myszy powiększa wokół kursora; użyj „Przywróć widok wykresu”, aby wrócić do pełnego zakresu. Dwuklik przy wykrytym piku pozwala edytować lub dodać jego nazwę. Podgląd Kreatora wykresów nie zmienia powiększenia kółkiem myszy.")
        help_label.setWordWrap(True)
        llay.addWidget(help_label)

        self.summary_label = QLabel("Nie wczytano próbki."); self.summary_label.setWordWrap(True)
        llay.addWidget(self.summary_label)
        llay.addStretch(1)

        # Right area with plot + table tabs
        self.tabs = QTabWidget(); splitter.addWidget(self.tabs)
        plot_tab = QWidget(); plot_lay = QVBoxLayout(plot_tab)
        self.fig = Figure(figsize=(9, 6)); self.canvas = FigureCanvas(self.fig); self.ax = self.fig.add_subplot(111)
        plot_lay.addWidget(self.canvas)
        self.tabs.addTab(plot_tab, "Wykres")

        table_tab = QWidget(); table_lay = QVBoxLayout(table_tab)
        self.table = QTableWidget(0, 13)
        self.table.setHorizontalHeaderLabels([
            "Próbka", "Typ materiału", "2θ [°]", "d [Å]", "Wysokość piku",
            "Intensywność względna [%]", "Szerokość ½ prominencji [°2θ]", "Prominencja", "SNR",
            "Prawdopodobna faza krystalograficzna", "hkl", "Etykieta piku", "Maks. całkowita intensywność"
        ])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(9, QHeaderView.ResizeMode.Stretch)
        self.table.horizontalHeader().setSectionResizeMode(11, QHeaderView.ResizeMode.Stretch)
        self.table.setSortingEnabled(True)
        table_lay.addWidget(self.table)
        self.tabs.addTab(table_tab, "Tabela pików")

        splitter.setStretchFactor(0, 0); splitter.setStretchFactor(1, 1); splitter.setSizes([450, 1050])

        # Signals
        self.open_project_btn.clicked.connect(self.open_project)
        self.save_project_btn.clicked.connect(self.save_project)
        self.save_project_as_btn.clicked.connect(self.save_project_as)
        self.load_btn.clicked.connect(self.open_files)
        self.analyze_btn.clicked.connect(self.analyze_current)
        self.analyze_all_btn.clicked.connect(self.analyze_all)
        self.batch_chemistry_btn.clicked.connect(self.open_batch_chemistry)
        self.plot_creator_btn.clicked.connect(self.open_plot_creator)
        self.save_plot_btn.clicked.connect(self.save_plot)
        self.export_btn.clicked.connect(self.export_table)
        self.db_btn.clicked.connect(self.open_database)
        self.cod_pattern_btn.clicked.connect(self.search_cod_measured_pattern)
        self.ai_review_btn.clicked.connect(self.open_ai_phase_review)
        self.instrument_btn.clicked.connect(self.open_instrument_metadata)
        self.close_program_btn.clicked.connect(self.close)
        self.allowed_elements_btn.clicked.connect(self.open_current_periodic_table)
        self.sample_setup_btn.clicked.connect(self.open_sample_setup_wizard)
        self.allowed_elements_edit.editingFinished.connect(self.update_allowed_elements)
        self.reset_identification_btn.clicked.connect(self.reset_identification_defaults)
        self.reset_peak_defaults_btn.clicked.connect(self.reset_peak_detection_defaults)
        self.reset_plot_defaults_btn.clicked.connect(self.reset_main_plot_defaults)
        self.reset_view_btn.clicked.connect(self.reset_plot_view)
        self.phase_source_combo.currentTextChanged.connect(self.phase_source_changed)
        self.sample_combo.currentIndexChanged.connect(self.switch_dataset)
        self.prev_btn.clicked.connect(self.previous_dataset)
        self.next_btn.clicked.connect(self.next_dataset)
        self.sample_name.editingFinished.connect(self.update_metadata)
        self.material_type.editingFinished.connect(self.update_metadata)
        self.wavelength.valueChanged.connect(self.update_metadata)
        self.annotate_cb.stateChanged.connect(self._plot_setting_changed)
        self.label_mode_combo.currentIndexChanged.connect(self._plot_setting_changed)
        self.peak_label_style_combo.currentTextChanged.connect(self._plot_setting_changed)
        self.hkl_display_combo.currentTextChanged.connect(self._plot_setting_changed)
        self.label_top_n.valueChanged.connect(self._plot_setting_changed)
        self.peak_font_spin.valueChanged.connect(self._plot_setting_changed)
        self.show_baseline_cb.stateChanged.connect(self._plot_setting_changed)
        self.title_edit.editingFinished.connect(self._plot_setting_changed)
        self.xlabel_edit.editingFinished.connect(self._plot_setting_changed)
        self.ylabel_edit.editingFinished.connect(self._plot_setting_changed)
        self.dpi_combo.currentTextChanged.connect(lambda *_: self.mark_project_dirty())
        self.baseline_cb.stateChanged.connect(lambda *_: self.mark_project_dirty())
        self.prominence.valueChanged.connect(lambda *_: self.mark_project_dirty())
        self.min_distance.valueChanged.connect(lambda *_: self.mark_project_dirty())
        self.phase_tol.valueChanged.connect(lambda *_: self.mark_project_dirty())
        self.min_rel_height.valueChanged.connect(lambda *_: self.mark_project_dirty())
        self.min_snr.valueChanged.connect(lambda *_: self.mark_project_dirty())
        self.min_fwhm.valueChanged.connect(lambda *_: self.mark_project_dirty())
        self.max_peaks.valueChanged.connect(lambda *_: self.mark_project_dirty())
        self.min_phase_matches.valueChanged.connect(lambda *_: self.mark_project_dirty())
        self.table.itemChanged.connect(self._table_item_changed)
        self.canvas.mpl_connect("scroll_event", self._on_scroll)
        self.canvas.mpl_connect("button_press_event", self._on_plot_double_click)
        self._update_nav_buttons()
        # Small screens still expose every control by scrolling, without forcing
        # the top-level window outside the available desktop area.
        central = self.takeCentralWidget()
        central_scroll = QScrollArea()
        central_scroll.setWidgetResizable(True)
        central_scroll.setWidget(central)
        central_scroll.setFrameShape(QFrame.Shape.NoFrame)
        self.setCentralWidget(central_scroll)
        self._restore_window_geometry()

    def _build_menus(self):
        """Menu bar mirroring the toolbar, with keyboard shortcuts.

        Slots are wrapped in a lambda so that the 'checked' flag of QAction.triggered is never
        passed into optional arguments (e.g. open_sample_setup_wizard(dataset=None)).
        """
        def add(menu, text, slot, shortcut=None, tip=""):
            action = QAction(text, self)
            if shortcut:
                action.setShortcut(QKeySequence(shortcut))
            if tip:
                action.setStatusTip(tip)
            action.triggered.connect(lambda _checked=False, target=slot: target())
            menu.addAction(action)
            return action

        file_menu = self.menuBar().addMenu("&Plik")
        add(file_menu, "Otwórz projekt…", self.open_project, "Ctrl+O", "Otwórz zapisany projekt .xrdproj")
        add(file_menu, "Zapisz projekt", self.save_project, "Ctrl+S", "Zapisz bieżący projekt")
        add(file_menu, "Zapisz projekt jako…", self.save_project_as, "Ctrl+Shift+S")
        file_menu.addSeparator()
        add(file_menu, "Otwórz pliki XRD…", self.open_files, "Ctrl+Shift+O", "Wczytaj jeden lub wiele dyfraktogramów")
        file_menu.addSeparator()
        add(file_menu, "Zamknij program", self.close, "Ctrl+Q")

        analysis_menu = self.menuBar().addMenu("&Analiza")
        add(analysis_menu, "Analizuj bieżący", self.analyze_current, "F5", "Wykryj piki i zidentyfikuj fazy bieżącej próbki")
        add(analysis_menu, "Analizuj wszystkie wczytane", self.analyze_all, "Shift+F5")
        analysis_menu.addSeparator()
        add(analysis_menu, "Konfiguracja próbki…", self.open_sample_setup_wizard)
        add(analysis_menu, "Chemia serii…", self.open_batch_chemistry)
        add(analysis_menu, "Dozwolone pierwiastki (układ okresowy)…", self.open_current_periodic_table)
        add(analysis_menu, "Aparatura / pomiar…", self.open_instrument_metadata)

        science_menu = self.menuBar().addMenu("Analiza naukowa")
        add(science_menu, "Profile / Scherrer / WH / sin²ψ / Wilkens…", self.open_scientific_analysis, "Ctrl+Shift+A")
        add(science_menu, "Kalibracja z wzorca Si…", self.open_calibration, "Ctrl+Shift+K",
            "Popraw położenia pików i profil aparatury skanem wzorca krzemowego zmierzonym razem z próbkami")

        plot_menu = self.menuBar().addMenu("&Wykres")
        add(plot_menu, "Kreator wykresów…", self.open_plot_creator, "Ctrl+K")
        add(plot_menu, "Zapisz bieżący wykres…", self.save_plot, "Ctrl+E")
        add(plot_menu, "Eksportuj tabelę pików…", self.export_table, "Ctrl+T")
        plot_menu.addSeparator()
        add(plot_menu, "Przywróć widok wykresu", self.reset_plot_view, "Ctrl+0")

        db_menu = self.menuBar().addMenu("&Bazy faz")
        add(db_menu, "Bazy faz / COD…", self.open_database)
        add(db_menu, "Wyszukaj COD dla zmierzonego dyfraktogramu…", self.search_cod_measured_pattern)
        add(db_menu, "Ocena faz wspomagana AI…", self.open_ai_phase_review)

        help_menu = self.menuBar().addMenu("Pomo&c")
        add(help_menu, "Skróty klawiszowe i wskazówki…", self.show_shortcuts_help, "F1")
        add(help_menu, "Badania XRD — analiza, wykresy i formaty plików…", self.show_xrd_overview)
        add(help_menu, "Bazy COD — skąd pobrać i jak wczytać?…", self.show_cod_download_help)
        add(help_menu, "O programie…", self.show_about)

    def _status(self, text: str, timeout_ms: int = 8000):
        """Non-blocking feedback in the status bar (never raises)."""
        try:
            self.statusBar().showMessage(str(text), timeout_ms)
        except Exception:
            pass

    def _last_dir(self) -> str:
        value = str(self.settings.value("ui/last_dir", "") or "")
        return value if value and Path(value).is_dir() else ""

    def _remember_dir(self, path):
        try:
            self.settings.setValue("ui/last_dir", str(Path(str(path)).parent))
        except Exception:
            pass

    def _default_save_path(self, filename: str) -> str:
        folder = self._last_dir()
        return str(Path(folder) / filename) if folder else filename

    def show_shortcuts_help(self):
        QMessageBox.information(
            self, "Skróty klawiszowe i wskazówki",
            "Skróty klawiszowe\n"
            "  Ctrl+O – otwórz projekt        Ctrl+S – zapisz projekt        Ctrl+Shift+S – zapisz jako…\n"
            "  Ctrl+Shift+O – wczytaj pliki XRD\n"
            "  F5 – analizuj bieżący           Shift+F5 – analizuj wszystkie\n"
            "  Ctrl+K – kreator wykresów     Ctrl+E – zapisz wykres        Ctrl+T – eksportuj tabelę pików\n"
            "  Ctrl+Shift+A – analiza naukowa    Ctrl+Shift+K – kalibracja z wzorca Si\n  Ctrl+0 – przywróć widok wykresu    Ctrl+Q – zamknij\n\n"
            "Typowy przebieg pracy\n"
            "  1. Wczytaj pliki XRD i sprawdź długość fali λ oraz metadane aparatury.\n"
            "  2. Podaj skład próbki (Konfiguracja próbki / Chemia serii) – bez tego fazy nie są dopasowywane chemicznie.\n"
            "  3. Ustaw parametry detekcji pików i naciśnij Analizuj (F5).\n"
            "  4. Zapisz projekt – po nadaniu nazwy kolejne zmiany zapisują się automatycznie.\n\n"
            "Wskazówki\n"
            "  • Kółko myszy powiększa wykres wokół kursora; „Przywróć widok wykresu” wraca do pełnego zakresu.\n"
            "  • Dwuklik przy piku pozwala nadać własną etykietę.\n"
            "  • Kliknięcie nagłówka kolumny sortuje tabelę pików (liczby sortowane są numerycznie).\n"
            "  • Wartości w polach liczbowych zatwierdzaj klawiszem Enter lub przejściem do innego pola.")

    def show_xrd_overview(self):
        from xrd_help_content import XRD_OVERVIEW_HTML, show_help_dialog
        show_help_dialog(self, "Badania XRD — zastosowania, wykresy i formaty", XRD_OVERVIEW_HTML)

    def show_cod_download_help(self):
        from xrd_help_content import COD_HELP_HTML, show_help_dialog
        show_help_dialog(self, "COD — pobieranie i baza offline", COD_HELP_HTML)

    def show_about(self):
        QMessageBox.information(
            self, "O programie",
            f"{APP_NAME} {APP_VERSION}\nFormat projektu: {self.PROJECT_VERSION}\n\n"
            "Analiza dyfraktogramów XRD: detekcja pików, identyfikacja faz (baza wewnętrzna i COD), "
            "profile linii, Scherrer / Williamson–Hall, naprężenia sin²ψ oraz kreator wykresów.\n\n"
            "Import: .udf, .rd (V3RD), .xrdml/.xml, .cif (pdCIF), .csv/.txt/.dat/.asc.\n"
            "Eksport: .png/.tif/.tiff/.pdf/.svg/.eps, .csv/.xlsx. Projekty: .xrdproj/.json.\n"
            "Więcej: Pomoc → Badania XRD — analiza, wykresy i formaty plików.\n\n"
            "Wyniki materiałowe mają status warunkowy i wymagają walidacji na wzorcu tej samej aparatury.")

    def _available_window_screens(self):
        # Keep space for a normal window's title bar and frame. Qt coordinates
        # are logical pixels, so the same code works across DPI scaling changes.
        screens = []
        for screen in QApplication.screens():
            rect = screen.availableGeometry()
            margin_x = min(8, max(0, rect.width() // 20))
            margin_y = min(32, max(0, rect.height() // 20))
            screens.append((rect.x() + margin_x, rect.y() + margin_y,
                            rect.width() - 2 * margin_x, rect.height() - margin_y - margin_x))
        return screens

    def _restore_window_geometry(self):
        """Restore machine-local UI state; projects never own window geometry."""
        restored = False
        geometry = self.settings.value("ui/main_geometry")
        if geometry is not None:
            try:
                restored = bool(self.restoreGeometry(geometry))
            except (TypeError, ValueError):
                pass
        saved_rect = self.settings.value("ui/main_normal_geometry")
        if saved_rect is None and restored:
            rect = self.normalGeometry()
            saved_rect = [rect.x(), rect.y(), rect.width(), rect.height()]
        self.setWindowState(Qt.WindowState.WindowNoState)
        self.setGeometry(*safe_window_rect(saved_rect, self._available_window_screens()))
        saved_state = self.settings.value("ui/main_state")
        if saved_state is not None:
            try:
                self.restoreState(saved_state, 1)
            except (TypeError, ValueError):
                pass
        # Preserve the historical maximized first launch, but respect subsequent
        # normal-window launches. Never restore a minimized or fullscreen state.
        maximized = optional_bool(self.settings.value("ui/main_maximized", True), True)
        if maximized:
            self.setWindowState(Qt.WindowState.WindowMaximized)
        return restored

    def _save_window_geometry(self):
        rect = self.normalGeometry() if self.isMaximized() or self.isMinimized() else self.geometry()
        safe = safe_window_rect([rect.x(), rect.y(), rect.width(), rect.height()],
                                self._available_window_screens())
        self.settings.setValue("ui/main_geometry", self.saveGeometry())
        self.settings.setValue("ui/main_normal_geometry", list(safe))
        self.settings.setValue("ui/main_state", self.saveState(1))
        self.settings.setValue("ui/main_maximized", bool(self.isMaximized()))
        self.settings.sync()

    def open_calibration(self):
        if not self.datasets:
            QMessageBox.information(self, "Kalibracja", "Najpierw wczytaj skan wzorca Si i skany próbek.")
            return
        from xrd_calibration_dialog import CalibrationDialog
        try:
            dialog = CalibrationDialog(self.datasets, self.sample_combo.currentIndex(), self)
            dialog.exec()
        except Exception as exc:
            QMessageBox.critical(self, "Kalibracja", f"{exc}\n\n{traceback.format_exc()}")

    def open_scientific_analysis(self):
        ds = self.current_dataset()
        if ds is None:
            QMessageBox.information(self, "Analiza naukowa", "Najpierw wczytaj dyfraktogram.")
            return
        from xrd_analysis_dialog import ScientificAnalysisDialog
        try:
            dialog = ScientificAnalysisDialog(ds, self)
            dialog.exec()
        except Exception as exc:
            QMessageBox.critical(self, "Analiza naukowa", f"{exc}\n\n{traceback.format_exc()}")

    # --------------------------- reset helpers ---------------------------
    def reset_identification_defaults(self):
        ds = self.current_dataset()
        if ds is None:
            return
        self.allowed_elements_edit.clear()
        filename_hint = str(ds.metadata.get("filename_material_suggestion", "") or "")
        text_hint = self.material_type.text().strip() or ds.material_type or filename_hint or self.sample_name.text().strip() or ds.name
        suggested = default_elements_for_material(text_hint, "")
        self.allowed_elements_edit.setPlaceholderText(("Niezweryfikowana sugestia (niezastosowana): " + suggested) if suggested else "Ustaw chemię w Kreatorze konfiguracji próbki lub w układzie okresowym")
        self.phase_source_combo.setCurrentText("Tylko baza wewnętrzna")
        try:
            wl = float(ds.metadata.get("wavelength_A", ds.wavelength_A))
        except Exception:
            wl = DEFAULT_WAVELENGTH_A
        if not math.isfinite(wl) or wl <= 0:
            wl = DEFAULT_WAVELENGTH_A
        self.wavelength.setValue(wl)
        ds.metadata.pop("allowed_elements", None)
        ds.metadata_sources.pop("allowed_elements", None)
        ds.metadata.pop("material_preset_confirmed", None)
        ds.metadata_sources.pop("material_preset_confirmed", None)
        ds.metadata["composition_mode"] = "unknown"
        ds.metadata["chemical_composition"] = []
        ds.metadata["phase_identification_source"] = "internal database"
        ds.metadata_sources["phase_identification_source"] = "ustawienie domyślne"
        self.mark_project_dirty()
        self.refresh_current_metadata_view()

    def reset_peak_detection_defaults(self):
        self.baseline_cb.setChecked(True)
        self.prominence.setValue(3.0)
        self.min_rel_height.setValue(5.0)
        self.min_snr.setValue(8.0)
        self.min_fwhm.setValue(0.10)
        self.max_peaks.setValue(50)
        self.min_distance.setValue(0.25)
        self.phase_tol.setValue(0.35)
        self.min_phase_matches.setValue(2)
        self.mark_project_dirty()

    def reset_main_plot_defaults(self):
        self.title_edit.setText("Dyfraktogram XRD – profil pomiarowy i prawdopodobne fazy")
        self.xlabel_edit.setText("2θ [°]")
        self.ylabel_edit.setText("Intensywność [zliczenia]")
        self.annotate_cb.setChecked(True)
        self.label_mode_combo.setCurrentText("Tylko dopasowane fazy")
        self.peak_label_style_combo.setCurrentText(PEAK_LABEL_STYLE_FULL)
        self.hkl_display_combo.setCurrentText(HKL_DISPLAY_AUTO)
        self.label_top_n.setValue(12)
        self.peak_font_spin.setValue(7)
        self.show_baseline_cb.setChecked(False)
        self.dpi_combo.setCurrentText("600")
        self.mark_project_dirty()
        self.redraw()

    def reset_plot_view(self):
        """Restore the main plot to the full autoscaled range after interactive zoom."""
        self.redraw()

    # --------------------------- project persistence ---------------------------
    def mark_project_dirty(self):
        if self._loading_project:
            return
        self.project_dirty = True
        if self.project_path is not None:
            # Once a project has a path, edits are immediately persisted.
            try:
                self._write_project(self.project_path, silent=True)
                self.project_dirty = False
            except Exception as exc:
                self.project_dirty = True
                message = str(exc)
                if message != self._last_save_error:
                    QMessageBox.warning(self, "Automatyczny zapis projektu",
                                        "Nie udało się zapisać zmian. Pozostają w pamięci; użyj „Zapisz projekt” "
                                        "lub „Zapisz projekt jako…”.\n\n" + message)
                self._last_save_error = message

    def _serialize_peak(self, p: PeakResult) -> dict:
        return {
            "two_theta": p.two_theta,
            "d_A": p.d_A,
            "height": p.height,
            "rel_intensity": p.rel_intensity,
            "fwhm_deg": p.fwhm_deg,
            "width_half_prominence_deg": p.width_half_prominence_deg,
            "width_definition": "scipy_peak_width_half_prominence",
            "phase": p.phase,
            "hkl": p.hkl,
            "delta_2theta": None if not math.isfinite(p.delta_2theta) else p.delta_2theta,
            "custom_label": p.custom_label,
            "prominence": None if not math.isfinite(p.prominence) else p.prominence,
            "snr": None if not math.isfinite(p.snr) else p.snr,
            "phase_source": str(getattr(p, "phase_source", "") or ""),
        }

    def _project_dict(self, project_path: Path) -> dict:
        if self._plot_creator_dialog is not None:
            self._plot_creator_dialog._sync_to_project()
        base = project_path.parent.resolve()
        data = []
        for ds in self.datasets:
            try:
                rel = os.path.relpath(ds.path.resolve(), base)
            except Exception:
                rel = str(ds.path.resolve())
            data.append({
                "dataset_id": ds.dataset_id,
                "path": rel,
                "absolute_path": str(ds.path.resolve()),
                "embedded_data": f"datasets/{len(data):03d}.npz",
                "name": ds.name,
                "material_type": ds.material_type,
                "wavelength_A": ds.wavelength_A,
                "metadata": ds.metadata,
                "metadata_sources": ds.metadata_sources,
                "peaks": [self._serialize_peak(p) for p in ds.peaks],
            })
        return {
            "format": "XRDAnalyzerProject",
            "version": self.PROJECT_VERSION,
            "datasets": data,
            "current_index": self.sample_combo.currentIndex(),
            "analysis_settings": {
                "baseline": self.baseline_cb.isChecked(),
                "prominence_percent": self.prominence.value(),
                "min_relative_height_percent": self.min_rel_height.value(),
                "min_snr": self.min_snr.value(),
                "min_fwhm_deg": self.min_fwhm.value(),
                "max_detected_peaks": self.max_peaks.value(),
                "min_peak_distance_deg": self.min_distance.value(),
                "phase_tolerance_deg": self.phase_tol.value(),
                "min_phase_matches": self.min_phase_matches.value(),
                "phase_identification_source": self.phase_source_combo.currentText(),
            },
            "main_plot": {
                "title": self.title_edit.text(),
                "xlabel": self.xlabel_edit.text(),
                "ylabel": self.ylabel_edit.text(),
                "annotate_peaks": self.annotate_cb.isChecked(),
                "label_mode": self.label_mode_combo.currentText(),
                "peak_label_style": self.peak_label_style_combo.currentText(),
                "hkl_display": self.hkl_display_combo.currentText(),
                "label_top_n": self.label_top_n.value(),
                "peak_font_size": self.peak_font_spin.value(),
                "show_baseline": self.show_baseline_cb.isChecked(),
                "dpi": self.dpi_combo.currentText(),
            },
            "plot_creator": deepcopy(self.plot_creator_state),
            "cod_db_path": str(self.cod_db.path) if getattr(self.cod_db, "path", None) else "",
        }

    def _write_project(self, path: Path, silent: bool = False):
        """Write a self-contained compressed .xrdproj archive.

        The archive contains project.json plus one compressed NumPy data block per
        sample (2theta, raw intensity, baseline and corrected intensity). Therefore
        reopening the project does not require the original XRD/CIF/ASC files.
        """
        path = Path(path)
        if path.suffix.lower() != ".xrdproj":
            path = path.with_suffix(".xrdproj")
        payload = self._project_dict(path)
        tmp = path.with_suffix(path.suffix + ".tmp")
        try:
            with zipfile.ZipFile(tmp, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
                for i, (ds, rec) in enumerate(zip(self.datasets, payload.get("datasets", []))):
                    entry = str(rec.get("embedded_data") or f"datasets/{i:03d}.npz")
                    bio = io.BytesIO()
                    baseline = np.asarray(ds.baseline, dtype=float) if ds.baseline is not None else np.asarray([], dtype=float)
                    corrected = np.asarray(ds.corrected_y, dtype=float) if ds.corrected_y is not None else np.asarray([], dtype=float)
                    np.savez_compressed(
                        bio, x=np.asarray(ds.x, dtype=float), y=np.asarray(ds.y, dtype=float),
                        baseline=baseline, corrected_y=corrected,
                    )
                    zf.writestr(entry, bio.getvalue())
                zf.writestr("project.json", json.dumps(payload, ensure_ascii=False, indent=2))
            os.replace(tmp, path)
        finally:
            try:
                if tmp.exists():
                    tmp.unlink()
            except Exception:
                pass
        self.project_path = path
        self.project_dirty = False
        self._last_save_error = None
        self.setWindowTitle(f"XRD Analyzer – {path.name}")
        self.settings.setValue("last_project", str(path))
        if not silent:
            QMessageBox.information(self, "Zapisz projekt",
                                    f"Samowystarczalny projekt został zapisany:\n{path}\n\nOryginalne pliki wejściowe XRD nie są już potrzebne do ponownego otwarcia projektu.")

    def save_project(self):
        if self.project_path is None:
            return self.save_project_as()
        try:
            self._write_project(self.project_path, silent=False)
            return True
        except Exception as e:
            QMessageBox.critical(self, "Zapisz projekt", f"Nie udało się zapisać projektu:\n{e}\n\n{traceback.format_exc()}")
            return False

    def save_project_as(self):
        default_name = str(self.project_path) if self.project_path else self._default_save_path("XRD_analysis.xrdproj")
        path, _ = QFileDialog.getSaveFileName(self, "Zapisz projekt XRD Analyzer", default_name, "Projekt XRD Analyzer (*.xrdproj);;JSON (*.json)")
        if not path:
            return False
        self._remember_dir(path)
        try:
            self._write_project(Path(path), silent=False)
            return True
        except Exception as e:
            QMessageBox.critical(self, "Zapisz projekt", f"Nie udało się zapisać projektu:\n{e}\n\n{traceback.format_exc()}")
            return False

    def _can_replace_project(self):
        if self._plot_creator_dialog is not None:
            self._plot_creator_dialog._sync_to_project()
        if not self.datasets or not self.project_dirty:
            return True
        answer = QMessageBox.question(self, "Niezapisane zmiany",
            "Zapisać bieżący projekt przed otwarciem innego?",
            QMessageBox.StandardButton.Save | QMessageBox.StandardButton.Discard | QMessageBox.StandardButton.Cancel,
            QMessageBox.StandardButton.Save)
        if answer == QMessageBox.StandardButton.Discard:
            return True
        if answer == QMessageBox.StandardButton.Save:
            return bool(self.save_project())
        return False

    def _capture_project_memory(self):
        names = ("baseline_cb", "prominence", "min_rel_height", "min_snr", "min_fwhm", "max_peaks",
                 "min_distance", "phase_tol", "min_phase_matches", "phase_source_combo", "title_edit", "xlabel_edit",
                 "ylabel_edit", "annotate_cb", "label_mode_combo", "peak_label_style_combo", "hkl_display_combo",
                 "label_top_n", "peak_font_spin", "show_baseline_cb", "dpi_combo")
        widgets = []
        for name in names:
            w = getattr(self, name)
            if isinstance(w, QCheckBox): setter, value = "setChecked", w.isChecked()
            elif isinstance(w, QComboBox): setter, value = "setCurrentText", w.currentText()
            elif isinstance(w, QLineEdit): setter, value = "setText", w.text()
            else: setter, value = "setValue", w.value()
            widgets.append((w, setter, value))
        return dict(datasets=list(self.datasets), index=self.sample_combo.currentIndex(),
                    path=self.project_path, dirty=self.project_dirty, title=self.windowTitle(),
                    creator=deepcopy(self.plot_creator_state), widgets=widgets)

    def _restore_project_memory(self, snapshot):
        self._loading_project = True
        blockers = [QSignalBlocker(w) for w, _, _ in snapshot["widgets"]] + [QSignalBlocker(self.sample_combo)]
        for widget, setter, value in snapshot["widgets"]:
            getattr(widget, setter)(value)
        self.datasets = snapshot["datasets"]
        self.sample_combo.clear(); self.sample_combo.addItems([d.name for d in self.datasets])
        self.sample_combo.setCurrentIndex(snapshot["index"])
        self.project_path = snapshot["path"]; self.project_dirty = snapshot["dirty"]
        self.plot_creator_state = snapshot["creator"]; self.setWindowTitle(snapshot["title"])
        del blockers
        self.switch_dataset(snapshot["index"])
        self._loading_project = False

    def open_project(self):
        path, _ = QFileDialog.getOpenFileName(self, "Otwórz projekt XRD Analyzer", self._last_dir(), "Projekt XRD Analyzer (*.xrdproj *.json);;Wszystkie pliki (*.*)")
        if not path:
            return
        if not self._can_replace_project():
            return
        project_zip = None
        snapshot = None
        try:
            project_file = Path(path)
            if zipfile.is_zipfile(project_file):
                project_zip = zipfile.ZipFile(project_file, "r")
                payload = json.loads(project_zip.read("project.json").decode("utf-8"))
            else:
                # Backward compatibility with legacy JSON .xrdproj files.
                payload = json.loads(project_file.read_text(encoding="utf-8"))
            if payload.get("format") != "XRDAnalyzerProject":
                raise ValueError("To nie jest plik projektu XRD Analyzer.")
            if self._plot_creator_dialog is not None:
                self._plot_creator_dialog.reject()
            snapshot = self._capture_project_memory()
            self._loading_project = True
            self.datasets.clear(); self.sample_combo.clear()
            project_dir = Path(path).parent.resolve()
            missing = []

            aset = payload.get("analysis_settings", {})
            self.baseline_cb.setChecked(bool(aset.get("baseline", True)))
            self.prominence.setValue(float(aset.get("prominence_percent", 3.0)))
            self.min_rel_height.setValue(float(aset.get("min_relative_height_percent", 5.0)))
            self.min_snr.setValue(float(aset.get("min_snr", 8.0)))
            self.min_fwhm.setValue(float(aset.get("min_fwhm_deg", 0.10)))
            self.max_peaks.setValue(int(aset.get("max_detected_peaks", 50)))
            self.min_distance.setValue(float(aset.get("min_peak_distance_deg", 0.25)))
            self.phase_tol.setValue(float(aset.get("phase_tolerance_deg", 0.35)))
            self.min_phase_matches.setValue(int(aset.get("min_phase_matches", 2)))
            phase_source = _migrate_phase_source(str(aset.get("phase_identification_source", self.settings.value("phase_identification_source", "Tylko baza wewnętrzna"))))
            if self.phase_source_combo.findText(phase_source) >= 0:
                self.phase_source_combo.setCurrentText(phase_source)

            pset = payload.get("main_plot", {})
            self.title_edit.setText(str(pset.get("title", "Dyfraktogram XRD – profil pomiarowy i prawdopodobne fazy")))
            self.xlabel_edit.setText(str(pset.get("xlabel", "2θ [°]")))
            self.ylabel_edit.setText(str(pset.get("ylabel", "Intensywność [zliczenia]")))
            self.annotate_cb.setChecked(bool(pset.get("annotate_peaks", True)))
            mode = _migrate_label_mode(str(pset.get("label_mode", "Tylko dopasowane fazy")))
            if self.label_mode_combo.findText(mode) >= 0:
                self.label_mode_combo.setCurrentText(mode)
            lbl_style = _migrate_peak_label_style(str(pset.get("peak_label_style", PEAK_LABEL_STYLE_FULL)))
            if self.peak_label_style_combo.findText(lbl_style) >= 0:
                self.peak_label_style_combo.setCurrentText(lbl_style)
            hkl_mode = _migrate_hkl_display(str(pset.get("hkl_display", HKL_DISPLAY_AUTO)))
            if self.hkl_display_combo.findText(hkl_mode) >= 0:
                self.hkl_display_combo.setCurrentText(hkl_mode)
            self.label_top_n.setValue(int(pset.get("label_top_n", 12)))
            self.peak_font_spin.setValue(int(pset.get("peak_font_size", 7)))
            self.show_baseline_cb.setChecked(bool(pset.get("show_baseline", False)))
            self.dpi_combo.setCurrentText(str(pset.get("dpi", "600")))

            loaded_ids = set()
            for rec in payload.get("datasets", []):
                rel = Path(str(rec.get("path", "")))
                candidate = rel if rel.is_absolute() else project_dir / rel
                if not candidate.exists() and rec.get("absolute_path"):
                    candidate = Path(str(rec.get("absolute_path")))

                ds = None
                embedded = str(rec.get("embedded_data", "") or "")
                if project_zip is not None and embedded and embedded in project_zip.namelist():
                    try:
                        with np.load(io.BytesIO(project_zip.read(embedded)), allow_pickle=False) as arr:
                            x = np.asarray(arr["x"], dtype=float)
                            y = np.asarray(arr["y"], dtype=float)
                            base = np.asarray(arr["baseline"], dtype=float) if "baseline" in arr.files else np.asarray([])
                            corr = np.asarray(arr["corrected_y"], dtype=float) if "corrected_y" in arr.files else np.asarray([])
                        display_path = candidate if str(candidate) else Path(str(rec.get("absolute_path", rec.get("name", "embedded.xrd"))))
                        ds = Dataset(
                            path=display_path, name=str(rec.get("name", display_path.stem)),
                            material_type=str(rec.get("material_type", "")), x=x, y=y,
                            wavelength_A=float(rec.get("wavelength_A", DEFAULT_WAVELENGTH_A)),
                        )
                        ds.baseline = base if len(base) == len(y) and len(base) else None
                        ds.corrected_y = corr if len(corr) == len(y) and len(corr) else None
                    except Exception as e:
                        missing.append(f"embedded data for {rec.get('name', '?')}: {e}")
                        ds = None

                if ds is None:
                    if not candidate.exists():
                        missing.append(str(candidate))
                        continue
                    ds = load_dataset(candidate)

                ds.name = str(rec.get("name", ds.name))
                saved_id = rec.get("dataset_id")
                if isinstance(saved_id, str) and saved_id and saved_id not in loaded_ids:
                    ds.dataset_id = saved_id
                loaded_ids.add(ds.dataset_id)
                ds._project_source_paths = [
                    str(value) for value in (rec.get("absolute_path"), candidate)
                    if value
                ]
                ds.material_type = str(rec.get("material_type", ds.material_type))
                ds.wavelength_A = float(rec.get("wavelength_A", ds.wavelength_A))
                if isinstance(rec.get("metadata"), dict):
                    ds.metadata.update(rec.get("metadata", {}))
                if isinstance(rec.get("metadata_sources"), dict):
                    ds.metadata_sources.update({str(k): str(v) for k, v in rec.get("metadata_sources", {}).items()})
                # v18 safety migration: old projects may contain chemistry that v16 injected
                # automatically from filenames/project context.  Remove only those automatic
                # records; user-entered chemistry remains untouched.
                self._discard_legacy_automatic_chemistry(ds)
                ds.metadata["wavelength_A"] = ds.wavelength_A
                ds.metadata_sources.setdefault("wavelength_A", "project")
                # Legacy projects do not contain processed arrays. New self-contained
                # projects restore them exactly and only recompute if they are absent.
                if self.baseline_cb.isChecked() and (ds.baseline is None or ds.corrected_y is None):
                    try:
                        ds.baseline = baseline_asls(ds.y, x=ds.x)
                        ds.corrected_y = np.clip(ds.y - ds.baseline, 0, None)
                    except Exception:
                        ds.baseline = np.zeros_like(ds.y)
                        ds.corrected_y = ds.y.copy()
                saved_peaks = []
                for pr in rec.get("peaks", []):
                    tt = float(pr.get("two_theta", 0.0))
                    idx = int(np.argmin(np.abs(ds.x - tt))) if len(ds.x) else 0
                    delta = pr.get("delta_2theta", None)
                    saved_peaks.append(PeakResult(
                        index=idx,
                        two_theta=tt,
                        d_A=float(pr.get("d_A", bragg_d(tt, ds.wavelength_A))),
                        height=float(pr.get("height", 0.0)),
                        rel_intensity=float(pr.get("rel_intensity", 0.0)),
                        fwhm_deg=float(pr.get("width_half_prominence_deg", pr.get("fwhm_deg", 0.0))),
                        phase=str(pr.get("phase", "unassigned")),
                        hkl=str(pr.get("hkl", "")),
                        delta_2theta=float(delta) if delta is not None else math.nan,
                        custom_label=str(pr.get("custom_label", "")),
                        prominence=float(pr.get("prominence")) if pr.get("prominence") is not None else math.nan,
                        snr=float(pr.get("snr")) if pr.get("snr") is not None else math.nan,
                        phase_source=str(pr.get("phase_source", "") or ""),
                    ))
                ds.peaks = saved_peaks
                self.datasets.append(ds)
                self.sample_combo.addItem(ds.name)

            if payload.get("datasets") and not self.datasets:
                raise ValueError("Nie odtworzono żadnej próbki. Bieżący projekt pozostaje zachowany.")
            raw_creator = payload.get("plot_creator", {})
            self.plot_creator_state = deepcopy(raw_creator) if isinstance(raw_creator, dict) else {}
            cod_path = str(payload.get("cod_db_path", ""))
            if cod_path and Path(cod_path).exists():
                try:
                    self.cod_db.open(Path(cod_path))
                    self.settings.setValue("cod_db_path", cod_path)
                except Exception:
                    pass

            self.project_path = Path(path)
            idx = int(payload.get("current_index", 0))
            if self.datasets:
                self.sample_combo.setCurrentIndex(min(max(idx, 0), len(self.datasets)-1))
            self.project_dirty = False
            self.setWindowTitle(f"XRD Analyzer – {self.project_path.name}")
            self.switch_dataset(self.sample_combo.currentIndex())
            self._loading_project = False
            if project_zip is not None:
                project_zip.close(); project_zip = None
            if missing:
                QMessageBox.warning(self, "Otwórz projekt", "Projekt otwarto, ale części danych starszego formatu lub danych osadzonych nie udało się odtworzyć:\n" + "\n".join(missing))
        except Exception as e:
            try:
                if project_zip is not None:
                    project_zip.close()
            except Exception:
                pass
            if snapshot is not None:
                self._restore_project_memory(snapshot)
            self._loading_project = False
            QMessageBox.critical(self, "Otwórz projekt", f"Nie udało się otworzyć projektu:\n{e}\n\n{traceback.format_exc()}")

    # --------------------------- datasets / navigation ---------------------------
    def current_dataset(self) -> Optional[Dataset]:
        i = self.sample_combo.currentIndex()
        return self.datasets[i] if 0 <= i < len(self.datasets) else None

    def _update_nav_buttons(self):
        i = self.sample_combo.currentIndex(); n = len(self.datasets)
        self.prev_btn.setEnabled(n > 1 and i > 0)
        self.next_btn.setEnabled(n > 1 and 0 <= i < n - 1)

    def previous_dataset(self):
        i = self.sample_combo.currentIndex()
        if i > 0:
            self.sample_combo.setCurrentIndex(i - 1)

    def next_dataset(self):
        i = self.sample_combo.currentIndex()
        if 0 <= i < len(self.datasets) - 1:
            self.sample_combo.setCurrentIndex(i + 1)

    def open_files(self):
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Otwórz pliki XRD", self._last_dir(),
            "Pliki XRD (*.cif *.CIF *.udf *.UDF *.rd *.RD *.xrdml *.XRDML *.xml *.csv *.CSV *.txt *.dat *.asc *.ASC);;Wszystkie pliki (*.*)"
        )
        if not paths:
            return
        self._remember_dir(paths[0])
        errors = []
        newly_loaded = []
        for p in paths:
            try:
                ds = load_dataset(Path(p))
                self.datasets.append(ds); self.sample_combo.addItem(ds.name); newly_loaded.append(ds)
            except Exception as e:
                errors.append(f"{Path(p).name}: {e}")
        if self.datasets:
            self.sample_combo.setCurrentIndex(len(self.datasets) - 1)
        self._update_nav_buttons()
        self.mark_project_dirty()
        if errors:
            QMessageBox.warning(self, "Niektórych plików nie udało się otworzyć", "\n".join(errors))
        rd_cautions = [(ds.name,int(ds.metadata.get("rd_axis_missing_points",0) or 0))
                       for ds in newly_loaded if ds.metadata.get("rd_axis_missing_points")]
        if rd_cautions:
            lines = [f"{name}: niezgodność {missing} punktów" for name,missing in rd_cautions]
            QMessageBox.warning(self, "Weryfikacja osi 2θ dla plików .rd",
                "Importer .rd zastosował heurystyczne zakotwiczenie osi na końcu skanu:\n"
                + "\n".join(lines) + "\n\nSprawdź pozycje wzorcowych refleksów, zanim użyjesz wyników naukowo. "
                "Oryginalny plik nie został zmieniony.")
        if newly_loaded:
            self._status(f"Wczytano plików XRD: {len(newly_loaded)}. Sprawdź λ i skład próbki, następnie naciśnij Analizuj (F5).")
            ans = QMessageBox.question(
                self, "Konfiguracja próbki",
                f"Czy skonfigurować teraz znany skład chemiczny dla {len(newly_loaded)} nowo wczytanych próbek?\n\n"
                "Możesz podać dokładny skład, zakresy, kolejność pierwiastków, prostą listę albo oznaczyć skład jako nieznany.",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
            )
            if ans == QMessageBox.StandardButton.Yes:
                for ds in newly_loaded:
                    self.open_sample_setup_wizard(ds)
        # Text/ASC files and some converted formats do not contain verified instrument metadata.
        needs_meta = [ds for ds in newly_loaded if _meta_missing(ds, "wavelength_A") or _meta_missing(ds, "instrument")]
        if needs_meta:
            ans = QMessageBox.question(
                self, "Metadane aparatury",
                f"{len(needs_meta)} nowo wczytanych plików nie zawiera kompletnych, zweryfikowanych metadanych aparatury.\n\n"
                "Czy otworzyć okno „Aparatura / pomiar” i zastosować zapisany profil?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
            )
            if ans == QMessageBox.StandardButton.Yes:
                self.open_instrument_metadata()

    def switch_dataset(self, idx):
        ds = self.current_dataset()
        self._update_nav_buttons()
        if not ds:
            return
        supplement_legacy_rd_provenance(ds)
        self.sample_name.blockSignals(True); self.material_type.blockSignals(True); self.wavelength.blockSignals(True); self.allowed_elements_edit.blockSignals(True)
        try:
            self.sample_name.setText(ds.name)
            self.material_type.setText(ds.material_type)
            self.wavelength.setValue(ds.wavelength_A)
            saved = str(ds.metadata.get("allowed_elements", "") or "").strip()
            self.allowed_elements_edit.setText(saved)
            filename_hint = str(ds.metadata.get("filename_material_suggestion", "") or "")
            hint = default_elements_for_material(ds.material_type, "") if ds.material_type else (default_elements_for_material(filename_hint, "") if filename_hint else "")
            self.allowed_elements_edit.setPlaceholderText(("Niezweryfikowana sugestia (niezastosowana): " + hint) if hint else "Ustaw chemię w Kreatorze konfiguracji próbki lub w układzie okresowym")
        finally:
            self.sample_name.blockSignals(False); self.material_type.blockSignals(False); self.wavelength.blockSignals(False); self.allowed_elements_edit.blockSignals(False)
        self.update_table(); self.redraw(); self.update_summary()

    def refresh_current_metadata_view(self):
        ds = self.current_dataset()
        if not ds:
            return
        self.wavelength.blockSignals(True)
        try:
            self.wavelength.setValue(ds.wavelength_A)
        finally:
            self.wavelength.blockSignals(False)
        self.update_table(); self.redraw(); self.update_summary()

    def phase_source_changed(self, text: str):
        self.settings.setValue("phase_identification_source", str(text))
        self.mark_project_dirty()
        ds = self.current_dataset()
        if ds and ds.peaks:
            try:
                self._assign_phases_for_dataset(ds)
                self.update_table(); self.redraw(); self.update_summary()
            except Exception as e:
                # Changing the selector should not destroy the already measured peak table.
                clear_phase_assignments(ds)
                self.update_table(); self.redraw(); self.update_summary()
                if "External COD" in str(text):
                    QMessageBox.warning(self, "Źródło identyfikacji faz", str(e))

    def _discard_legacy_automatic_chemistry(self, ds: Dataset) -> bool:
        """Remove chemistry injected by older versions from filenames/project context.

        User-entered chemistry is never touched.  This migration prevents a v16 project from
        silently continuing to use an automatic name-based preset after upgrading to v18.
        """
        src_allowed = str(ds.metadata_sources.get("allowed_elements", "") or "").lower()
        src_comp = str(ds.metadata_sources.get("chemical_composition", "") or "").lower()
        is_auto = ("automatic" in src_allowed or "auto preset" in src_allowed or
                   "automatic" in src_comp or "auto preset" in src_comp)
        if not is_auto:
            return False
        ds.metadata.pop("allowed_elements", None)
        ds.metadata_sources.pop("allowed_elements", None)
        if "automatic" in src_comp or "auto preset" in src_comp:
            ds.metadata.pop("chemical_composition", None)
            ds.metadata_sources.pop("chemical_composition", None)
            ds.metadata["composition_mode"] = "unknown"
        ds.metadata["legacy_auto_chemistry_removed"] = True
        ds.metadata_sources["legacy_auto_chemistry_removed"] = "v18 safety migration"
        return True

    def _ensure_allowed_elements_for_dataset(self, ds: Dataset) -> set[str]:
        """Return only chemistry that has been explicitly supplied/confirmed by the user.

        v18 deliberately does not infer chemistry from filenames, sample names or other loaded
        project files.  Built-in presets remain available in Sample setup, but require an
        explicit user action.
        """
        self._discard_legacy_automatic_chemistry(ds)
        return dataset_allowed_elements(ds)

    def _assign_phases_for_dataset(self, ds: Dataset):
        mode = self.phase_source_combo.currentText()
        elements = self._ensure_allowed_elements_for_dataset(ds)
        min_rel = max(2.0, self.min_rel_height.value())
        maxp = min(30, self.max_peaks.value())
        result = []
        if mode == "Bez automatycznej identyfikacji faz":
            clear_phase_assignments(ds)
            ds.metadata["phase_identification_source"] = "disabled"
            ds.metadata_sources["phase_identification_source"] = "analysis setting"
        elif mode == "Tylko zewnętrzna baza COD":
            if not elements:
                # Safe structural result instead of a misleading chemical guess or a page of
                # 'unassigned' labels. COD itself is not searched without user chemistry.
                clear_phase_assignments(ds)
                result = apply_structural_family_fallback(
                    ds, self.phase_tol.value(), min_rel_intensity=min(3.0,min_rel), max_pattern_peaks=maxp
                )
                ds.metadata["phase_identification_source"] = "structure-only fallback – chemistry unknown"
                ds.metadata_sources["phase_identification_source"] = "measured peak lattice-family fit"
            else:
                try:
                    result = assign_phases_from_cod(
                        ds, self.cod_db, self.phase_tol.value(), elements,
                        self.min_phase_matches.value(), min_rel, maxp,
                        max_candidates=8000, top_results=12,
                    )
                    ds.metadata.pop("cod_search_warning", None)
                except Exception as exc:
                    clear_phase_assignments(ds)
                    ds.metadata["cod_search_warning"] = str(exc)
                    ds.metadata_sources["cod_search_warning"] = "automatic analysis"
                    result = []
                apply_structural_family_fallback(
                    ds, self.phase_tol.value(), min_rel_intensity=min(3.0,min_rel), max_pattern_peaks=maxp,
                    overwrite_weak_cod=True
                )
        elif mode == "Baza wewnętrzna + zewnętrzna COD":
            if not elements:
                clear_phase_assignments(ds)
                result = apply_structural_family_fallback(
                    ds, self.phase_tol.value(), min_rel_intensity=min(3.0,min_rel), max_pattern_peaks=maxp
                )
                ds.metadata["phase_identification_source"] = "structure-only fallback – chemistry unknown"
                ds.metadata_sources["phase_identification_source"] = "measured peak lattice-family fit"
            else:
                result = combine_internal_and_cod_assignments(
                    ds, self.db, self.cod_db, self.phase_tol.value(), elements,
                    self.min_phase_matches.value(), min_rel, maxp,
                )
                apply_structural_family_fallback(
                    ds, self.phase_tol.value(), min_rel_intensity=min(3.0,min_rel), max_pattern_peaks=maxp,
                    overwrite_weak_cod=True
                )
                assign_tentative_single_reflections(ds, self.db, self.phase_tol.value(), elements)
        else:
            if not elements:
                clear_phase_assignments(ds)
                result = apply_structural_family_fallback(
                    ds, self.phase_tol.value(), min_rel_intensity=min(3.0,min_rel), max_pattern_peaks=maxp
                )
                ds.metadata["phase_identification_source"] = "structure-only fallback – chemistry unknown"
                ds.metadata_sources["phase_identification_source"] = "measured peak lattice-family fit"
            else:
                result = assign_phases(
                    ds, self.db, self.phase_tol.value(), elements,
                    min_phase_matches=self.min_phase_matches.value(),
                    min_pattern_rel_intensity=min_rel,
                    max_pattern_peaks=maxp,
                )
                apply_structural_family_fallback(
                    ds, self.phase_tol.value(), min_rel_intensity=min(3.0,min_rel), max_pattern_peaks=maxp,
                    overwrite_weak_cod=False
                )
                assign_tentative_single_reflections(ds, self.db, self.phase_tol.value(), elements)
        # Explicit user-approved COD peak choices always have final priority.
        apply_manual_cod_peak_overrides(ds)
        # Build an advisory, source-comparative ranking. This never changes assignments.
        try:
            build_ai_assisted_ranking(ds, self.phase_tol.value())
        except Exception as exc:
            ds.metadata["ai_phase_ranking"] = []
            ds.metadata["ai_phase_warning"] = f"{type(exc).__name__}: {exc}"
            ds.metadata_sources["ai_phase_warning"] = "AI-assisted review"
        return result

    def open_batch_chemistry(self):
        if not self.datasets:
            QMessageBox.information(self, 'Chemia serii', 'Najpierw otwórz jeden lub więcej plików XRD.')
            return
        # Store current edits first, then let the user explicitly configure the series.
        self.update_metadata()
        dlg=BatchChemistryDialog(self.datasets,self)
        if dlg.exec()==QDialog.DialogCode.Accepted:
            cur=self.current_dataset()
            if cur is not None:
                self.allowed_elements_edit.blockSignals(True)
                try: self.allowed_elements_edit.setText(str(cur.metadata.get('allowed_elements','') or ''))
                finally: self.allowed_elements_edit.blockSignals(False)
            self.mark_project_dirty(); self.update_summary(); self.redraw()

    def open_sample_setup_wizard(self, dataset: Optional[Dataset] = None):
        ds = dataset or self.current_dataset()
        if ds is None:
            QMessageBox.information(self, "Konfiguracja próbki", "Najpierw otwórz plik XRD.")
            return False
        dlg = SampleSetupWizard(ds, self)
        if dlg.exec() == QDialog.DialogCode.Accepted:
            # Refresh combo text for a renamed sample.
            try:
                idx = self.datasets.index(ds)
                if 0 <= idx < self.sample_combo.count():
                    self.sample_combo.setItemText(idx, ds.name)
                    if idx == self.sample_combo.currentIndex():
                        self.switch_dataset(idx)
            except Exception:
                pass
            self.mark_project_dirty()
            self.update_summary(); self.redraw()
            return True
        return False

    def open_current_periodic_table(self):
        ds = self.current_dataset()
        if not ds:
            QMessageBox.information(self, "Dozwolone pierwiastki", "Najpierw otwórz próbkę XRD.")
            return
        current = CODPatternSearchDialog._parse_elements(self.allowed_elements_edit.text())
        dlg = PeriodicTableDialog(current, self)
        if dlg.exec() == QDialog.DialogCode.Accepted:
            selected = dlg.selected_elements()
            ordered = sorted(selected, key=lambda e: ELEMENT_Z.get(e, 999))
            self.allowed_elements_edit.setText(" ".join(ordered))
            self.update_allowed_elements()

    def update_allowed_elements(self):
        ds = self.current_dataset()
        if not ds:
            return
        elements = CODPatternSearchDialog._parse_elements(self.allowed_elements_edit.text())
        ordered = sorted(elements, key=lambda e: ELEMENT_Z.get(e, 999))
        value = " ".join(ordered)
        self.allowed_elements_edit.setText(value)
        ds.metadata["allowed_elements"] = value
        ds.metadata_sources["allowed_elements"] = "user / periodic table"
        if ds.peaks:
            try:
                self._assign_phases_for_dataset(ds)
            except Exception as e:
                clear_phase_assignments(ds)
                if self.phase_source_combo.currentText() == "Tylko zewnętrzna baza COD":
                    QMessageBox.warning(self, "Identyfikacja zewnętrzna COD", str(e))
        self.update_table(); self.redraw(); self.update_summary(); self.mark_project_dirty()

    def open_instrument_metadata(self):
        ds = self.current_dataset()
        if not ds:
            QMessageBox.information(self, "Aparatura / pomiar", "Najpierw otwórz próbkę XRD.")
            return
        dlg = InstrumentMeasurementDialog(ds, self.datasets, self.instrument_profiles, self)
        dlg.exec()
        self.refresh_current_metadata_view()
        self.mark_project_dirty()

    def update_metadata(self):
        ds = self.current_dataset()
        if not ds:
            return
        ds.name = self.sample_name.text().strip() or ds.path.stem
        ds.material_type = self.material_type.text().strip()
        ds.metadata_sources["material_type"] = "user / main window"
        elements = CODPatternSearchDialog._parse_elements(self.allowed_elements_edit.text())
        if elements:
            ordered = sorted(elements, key=lambda e: ELEMENT_Z.get(e, 999))
            ds.metadata["allowed_elements"] = " ".join(ordered)
            ds.metadata_sources["allowed_elements"] = "user"
        else:
            ds.metadata.pop("allowed_elements", None)
            ds.metadata_sources.pop("allowed_elements", None)
        new_wave = float(self.wavelength.value())
        invalidated = False
        if abs(new_wave - ds.wavelength_A) > 1e-12:
            try:
                invalidated = set_dataset_wavelength(ds, new_wave, "user")
            except ValueError as exc:
                QMessageBox.warning(self, "Długość fali", str(exc))
                return
        i = self.sample_combo.currentIndex()
        if i >= 0:
            self.sample_combo.setItemText(i, ds.name)
        self.update_table(); self.redraw(); self.update_summary(); self.mark_project_dirty()
        if invalidated:
            QMessageBox.information(self, "Zmiana długości fali",
                "Przeliczono d dla nowej długości fali i unieważniono poprzednie wyniki fazowe. "
                "Uruchom ponowną identyfikację faz.")

    # --------------------------- analysis ---------------------------
    def _analyze_dataset_peaks_only(self, ds: Dataset):
        """Detect/quantify peaks without attempting phase identification."""
        x = np.asarray(ds.x, dtype=float)
        if (x.ndim != 1 or len(x) < 3 or len(x) != len(ds.y)
                or not np.all(np.isfinite(x)) or np.any(np.diff(x) <= 0)):
            raise ValueError("Oś 2θ musi zawierać co najmniej 3 skończone, ściśle rosnące wartości "
                             "i mieć tyle samo punktów co intensywność. Dane nie są automatycznie sortowane.")
        if self.baseline_cb.isChecked():
            ds.baseline = baseline_asls(ds.y, x=ds.x)
            sig = ds.y - ds.baseline
        else:
            ds.baseline = np.zeros_like(ds.y)
            sig = ds.y.copy()
        sig = np.asarray(sig, dtype=float)
        sig[~np.isfinite(sig)] = 0.0
        positive = sig[sig > 0]
        scale = float(np.nanmax(sig)) if len(sig) else 1.0
        if not math.isfinite(scale) or scale <= 0:
            scale = 1.0
        # Robust noise estimate from first differences suppresses slowly varying background.
        noise_sigma = noise_sigma_from_differences(sig)
        # Broad reflections (also on the flank of a strong one) need less than min_snr x point noise.
        idx, prominences = find_candidate_peaks(sig, noise_sigma, scale, self.prominence.value(),
                                                self.min_snr.value(), self.min_rel_height.value())
        if len(idx):
            _widths_pts, _width_heights, left_ips, right_ips = peak_widths(sig, idx, rel_height=0.5)
            sample_indices = np.arange(len(x), dtype=float)
            widths = np.interp(right_ips, sample_indices, x) - np.interp(left_ips, sample_indices, x)
            keep = widths >= self.min_fwhm.value()
            idx = idx[keep]; widths = widths[keep]; prominences = prominences[keep]
            selected = angular_distance_selection(idx, x, prominences, self.min_distance.value())
            idx = idx[selected]; widths = widths[selected]; prominences = prominences[selected]
            if len(idx) > self.max_peaks.value():
                order_strength = np.argsort(prominences)[::-1][:self.max_peaks.value()]
                idx = idx[order_strength]; widths = widths[order_strength]; prominences = prominences[order_strength]
            order = np.argsort(idx)
            idx = idx[order]; widths = widths[order]; prominences = prominences[order]
        maxh = float(np.nanmax(sig[idx])) if len(idx) else 1.0
        old_labels = [(p.two_theta, p.custom_label) for p in ds.peaks if p.custom_label]
        peaks = []
        for n, pi in enumerate(idx):
            tt = float(ds.x[pi]); h = float(sig[pi]); dA = bragg_d(tt, ds.wavelength_A)
            custom = ""
            if old_labels:
                nearest = min(old_labels, key=lambda q: abs(q[0] - tt))
                if abs(nearest[0] - tt) <= max(0.15, self.phase_tol.value()):
                    custom = nearest[1]
            snr = h / max(noise_sigma, 1e-12)
            peaks.append(PeakResult(
                int(pi), tt, dA, h, 100*h/maxh, float(widths[n]),
                custom_label=custom, prominence=float(prominences[n]), snr=float(snr)
            ))
        ds.peaks = peaks
        apply_calibration_to_fresh_peaks(ds)       # silicon-standard correction of 2theta and d (if the series has one)

    def _analyze_dataset(self, ds: Dataset):
        self._analyze_dataset_peaks_only(ds)
        self._assign_phases_for_dataset(ds)

    def analyze_current(self):
        ds = self.current_dataset()
        if not ds:
            QMessageBox.information(self, "Analiza", "Najpierw otwórz plik XRD.")
            return
        self.update_metadata()
        try:
            self._analyze_dataset(ds)
            # Synchronize any user-confirmed chemistry stored by Sample setup/COD search.
            self.allowed_elements_edit.blockSignals(True)
            try:
                self.allowed_elements_edit.setText(str(ds.metadata.get("allowed_elements", "") or ""))
            finally:
                self.allowed_elements_edit.blockSignals(False)
            self.update_table(); self.redraw(); self.update_summary(); self.mark_project_dirty()
            self._status(f"Przeanalizowano „{ds.name}”: wykryto pików {len(ds.peaks)}.")
        except Exception as e:
            QMessageBox.critical(self, "Błąd analizy", f"{e}\n\n{traceback.format_exc()}")

    def analyze_all(self):
        if not self.datasets:
            QMessageBox.information(self, "Analiza wsadowa", "Najpierw otwórz jeden lub więcej plików XRD.")
            return
        self.update_metadata()
        progress = QProgressDialog("Analizowanie wczytanych dyfraktogramów XRD…", "Anuluj", 0, len(self.datasets), self)
        progress.setWindowTitle("Analiza wsadowa XRD")
        progress.setWindowModality(Qt.WindowModality.WindowModal)
        errors = []
        skipped_chemistry = []
        canceled = False
        processed = len(self.datasets)
        for i, ds in enumerate(self.datasets):
            progress.setValue(i)
            progress.setLabelText(f"Analizowanie {ds.name} ({i+1}/{len(self.datasets)})")
            QApplication.processEvents()
            if progress.wasCanceled():
                canceled = True
                processed = i
                break
            try:
                elements = self._ensure_allowed_elements_for_dataset(ds)
                self._analyze_dataset(ds)
                if not elements and self.phase_source_combo.currentText() in ("Tylko zewnętrzna baza COD", "Baza wewnętrzna + zewnętrzna COD", "Tylko baza wewnętrzna"):
                    skipped_chemistry.append(ds.name)
            except Exception as e:
                errors.append(f"{ds.name}: {e}")
        progress.setValue(len(self.datasets))
        cur = self.current_dataset()
        if cur is not None:
            self.allowed_elements_edit.blockSignals(True)
            try:
                self.allowed_elements_edit.setText(str(cur.metadata.get("allowed_elements", "") or ""))
            finally:
                self.allowed_elements_edit.blockSignals(False)
        self.update_table(); self.redraw(); self.update_summary(); self.mark_project_dirty()
        if errors:
            msg = "Analiza wsadowa zakończona z błędami:\n" + "\n".join(errors)
            if canceled:
                msg += f"\n\nAnalizę przerwano na prośbę użytkownika (przeanalizowano {processed} z {len(self.datasets)} próbek)."
            if skipped_chemistry:
                msg += ("\n\nRozpoznanie rodziny strukturalnej wykonano, ale nie uruchomiono chemicznego dopasowania COD/bazy wewnętrznej dla "
                        f"{len(skipped_chemistry)} próbek, ponieważ skład chemiczny nie został zatwierdzony:\n" +
                        "\n".join(skipped_chemistry))
            QMessageBox.warning(self, "Analiza wsadowa", msg)
        elif canceled:
            QMessageBox.information(
                self, "Analiza wsadowa",
                f"Analizę przerwano na prośbę użytkownika. Przeanalizowano {processed} z {len(self.datasets)} próbek; "
                "pozostałe zachowały poprzednie wyniki.")
        elif skipped_chemistry:
            QMessageBox.information(
                self, "Analiza wsadowa",
                "Detekcja pików i rozpoznanie rodziny strukturalnej zostały wykonane dla wszystkich wczytanych próbek.\n\n"
                f"Nie wykonano chemicznego dopasowania faz dla {len(skipped_chemistry)} próbek. "
                "Użyj „Chemia serii…” lub „Kreator konfiguracji próbki…”, aby dodać pierwiastki zatwierdzone przez użytkownika, a następnie ponów analizę w celu uzyskania chemicznych nazw faz.\n\n" +
                "\n".join(skipped_chemistry)
            )
        else:
            QMessageBox.information(self, "Analiza wsadowa", "Analiza wsadowa została zakończona dla wszystkich wczytanych próbek.")

    # --------------------------- Kreator wykresów ---------------------------
    def open_plot_creator(self):
        if not self.datasets:
            QMessageBox.information(self, "Kreator wykresów", "Najpierw otwórz jeden lub więcej plików XRD.")
            return
        try:
            dlg = PlotCreatorDialog(self.datasets, self, self.plot_creator_state)
            dlg.setWindowState(dlg.windowState() | Qt.WindowState.WindowMaximized)
            dlg.exec()
            # done() has already committed both accepted and rejected closures.
            dlg.deleteLater()
        except Exception as exc:
            # A GUI callback exception must not terminate the entire analyzer.  Save a
            # traceback beside the program and show the user a concise diagnostic.
            import traceback
            try:
                log_path = Path(__file__).with_name("XRD_Analyzer_error.log")
                log_path.write_text(traceback.format_exc(), encoding="utf-8")
                log_note = f"\n\nDziennik diagnostyczny: {log_path}"
            except Exception:
                log_note = ""
            QMessageBox.critical(
                self,
                "Błąd Kreatora wykresów",
                "Nie udało się otworzyć Kreatora wykresów, ale program główny nadal działa.\n\n"
                f"{type(exc).__name__}: {exc}" + log_note,
            )

    # --------------------------- current plot / labels ---------------------------
    def _peak_text(self, ds: Dataset, p: PeakResult) -> str:
        return format_peak_label(ds, p, self.peak_label_style_combo.currentText(), self.hkl_display_combo.currentText())

    def update_summary(self):
        ds = self.current_dataset()
        if not ds:
            self.summary_label.setText("Nie wczytano próbki."); return
        inst = _display_text(ds.metadata.get("instrument"))
        rad = _display_text(ds.metadata.get("radiation"))
        tvtxt = _display_number(ds.metadata.get("tube_voltage_kV"), "kV")
        tctxt = _display_number(ds.metadata.get("tube_current_mA"), "mA")
        ratetxt = _display_number(ds.metadata.get("scan_rate_deg_s"), "°/s")
        fits = ds.metadata.get("phase_fits", []) if isinstance(ds.metadata, dict) else []
        fit_bits = []
        if isinstance(fits, list):
            for f in fits[:3]:
                if not isinstance(f, dict) or f.get("model") == "fixed":
                    continue
                nm = str(f.get("phase", "phase"))
                if len(nm) > 32:
                    nm = nm[:30] + "…"
                vals = []
                if f.get("a_A") is not None:
                    vals.append(f"a={float(f['a_A']):.5f} Å")
                if f.get("c_A") is not None:
                    vals.append(f"c={float(f['c_A']):.5f} Å")
                fit_bits.append(f"{nm}: " + ", ".join(vals))
        fit_line = ("\nDopasowana sieć z pików pomiarowych: " + " | ".join(fit_bits)) if fit_bits else ""
        ai_rows = ds.metadata.get("ai_phase_ranking", []) if isinstance(ds.metadata, dict) else []
        ai_line = ""
        if isinstance(ai_rows, list) and ai_rows:
            r0 = ai_rows[0] if isinstance(ai_rows[0], dict) else {}
            if r0:
                ai_line = (f"\nOcena wspomagana AI: {r0.get('publication_label') or r0.get('phase','kandydat')} "
                           f"| indeks dowodów {_safe_float(r0.get('ai_evidence_score_percent')):.1f}/100 "
                           f"| {r0.get('confidence_band','')} | źródło {r0.get('source','')}")
                confirmed = ds.metadata.get("ai_user_confirmed_candidate", {})
                if isinstance(confirmed, dict) and confirmed.get("phase"):
                    ai_line += " | zatwierdzone przez użytkownika"
                else:
                    ai_line += " | niezatwierdzone"
        cal_text = calibration_summary(ds)
        cal_line = f"\n{cal_text}" if cal_text else ""
        rd_missing = int(ds.metadata.get("rd_axis_missing_points", 0) or 0)
        rd_warning = (f"\nUWAGA .rd: nagłówek zawiera o {rd_missing} punktów więcej niż plik; "
                      "oś wyznaczono heurystycznie od końca skanu. Zweryfikuj wzorzec!") if rd_missing else ""
        alpha_note = (f"\nKandydaci rozdzielonych dubletów Kα₂: {len(dataset_kalpha2_candidates(ds))} "
                      "(wymagają weryfikacji)") if verified_wavelength_pair(ds) else ""
        self.summary_label.setText(
            f"Plik: {ds.path.name} | Format: {ds.metadata.get('source_format', ds.path.suffix)}\n"
            f"Punkty: {len(ds.x)} | Zakres: {ds.scan_start_deg:.3f}–{ds.scan_end_deg:.3f}° 2θ | Krok: {ds.scan_step_deg:.5f}°\n"
            f"Maks. zmierzona intensywność całkowita: {ds.max_total_intensity:.1f} | Wykryte piki: {len(ds.peaks)}\n"
            f"λ = {ds.wavelength_A:.6f} Å [{metadata_source(ds, 'wavelength_A')}] | Promieniowanie: {rad}\n"
            f"Aparatura: {inst} | Lampa: {tvtxt}, {tctxt} | Szybkość skanowania: {ratetxt}\n"
            f"Dozwolone pierwiastki do identyfikacji faz: {ds.metadata.get('allowed_elements', 'bez ograniczeń')} "
            f"[{ds.metadata_sources.get('allowed_elements', 'brak')}]\n"
            f"Źródło identyfikacji faz: {ds.metadata.get('phase_identification_source', self.phase_source_combo.currentText())}"
            f"{fit_line}{ai_line}{cal_line}{rd_warning}{alpha_note}\n"
            f"Filtr pików: rel≥{self.min_rel_height.value():g}%, SNR≥{self.min_snr.value():g}, "
            f"Szerokość ½ prominencji≥{self.min_fwhm.value():g}° | faza wymaga ≥{self.min_phase_matches.value()} dopasowanych refleksów"
        )

    def _plot_setting_changed(self, *args):
        self.redraw(); self.mark_project_dirty()

    def redraw(self):
        ds = self.current_dataset(); self.ax.clear()
        if not ds:
            self.canvas.draw_idle(); return
        yplot = ds.corrected_y if ds.corrected_y is not None else ds.y
        self.ax.plot(ds.x, yplot, lw=1.1, label="Pomiar")
        if self.show_baseline_cb.isChecked() and ds.baseline is not None and ds.corrected_y is None:
            self.ax.plot(ds.x, ds.baseline, lw=0.9, label="Linia tła")
        valid = []
        py = np.array([])
        if ds.peaks:
            valid = [p for p in ds.peaks if 0 <= p.index < len(yplot)]
            px = np.array([p.two_theta for p in valid]); py = np.array([yplot[p.index] for p in valid])
            self.ax.scatter(px, py, s=22, zorder=5, label="Wykryte piki")
        self.ax.set_title(self.title_edit.text())
        self.ax.set_xlabel(self.xlabel_edit.text()); self.ax.set_ylabel(self.ylabel_edit.text())
        self.ax.grid(alpha=0.18)
        self.ax.legend(loc="best", fontsize=8)
        self.ax.relim(); self.ax.autoscale_view()

        if valid and self.annotate_cb.isChecked():
            label_pairs = list(zip(valid, py))
            mode = self.label_mode_combo.currentText()
            if mode == "Tylko dopasowane fazy":
                label_pairs = [(p, yv) for p, yv in label_pairs if p.custom_label.strip() or (p.phase and p.phase.lower() != "unassigned")]
            elif mode == "Top N najsilniejszych pików":
                top_ids = {id(p) for p in sorted(valid, key=lambda q: q.rel_intensity, reverse=True)[:self.label_top_n.value()]}
                label_pairs = [(p, yv) for p, yv in label_pairs if id(p) in top_ids or p.custom_label.strip()]
            # In all modes, a custom label is always respected.
            if label_pairs:
                y0, y1 = self.ax.get_ylim(); span = max(1e-12, y1 - y0)
                self.ax.set_ylim(y0, y1 + 0.10 * span)
                y0, y1 = self.ax.get_ylim(); span = max(1e-12, y1 - y0)
                mixed = blended_transform_factory(self.ax.transData, self.ax.transAxes)
            kalpha2_indices = dataset_kalpha2_candidates(ds)
            pair_angle_ids = {id(ds.peaks[idx]) for idx in kalpha2_indices}
            close_lane = {}
            positions = sorted([p for p, _ in label_pairs], key=lambda p: p.two_theta)
            for pos, p0 in enumerate(positions):
                if pos and float(p0.two_theta)-float(positions[pos-1].two_theta) < 1.15:
                    close_lane[id(p0)] = 0.70 if pos % 2 else 0.91
                    close_lane[id(positions[pos-1])] = 0.91 if pos % 2 else 0.70
            for p, yv in label_pairs:
                yv = float(yv)
                frac = (yv - y0) / span
                label = ("Kα₂?" if id(p) in pair_angle_ids and not p.custom_label
                         else self._peak_text(ds, p))
                if not label:
                    continue
                if frac >= 0.68:
                    target_frac = 0.965
                else:
                    target_frac = min(0.82, max(frac + 0.08, 0.18))
                target_frac = close_lane.get(id(p),target_frac)
                ann = self.ax.annotate(
                    label,
                    xy=(p.two_theta, yv), xycoords="data",
                    xytext=(p.two_theta, target_frac), textcoords=mixed,
                    ha="right", va="center",
                    fontsize=self.peak_font_spin.value(), rotation=90,
                    rotation_mode="anchor",
                    arrowprops=dict(arrowstyle="-", lw=0.45, color="0.35"),
                    annotation_clip=True,
                )
                try:
                    ann.set_clip_on(True)
                    ann.set_clip_path(self.ax.patch)
                except Exception:
                    pass
        try:
            self.fig.tight_layout()
        except Exception:
            pass
        self.canvas.draw_idle()

    def _on_scroll(self, event):
        if event.inaxes is None:
            return
        ax = event.inaxes
        base_scale = 1.25
        scale = 1.0 / base_scale if event.button == "up" else base_scale
        x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
        xc = event.xdata if event.xdata is not None else (x0 + x1) / 2.0
        yc = event.ydata if event.ydata is not None else (y0 + y1) / 2.0
        nw = (x1 - x0) * scale; nh = (y1 - y0) * scale
        rx = (xc - x0) / (x1 - x0) if x1 != x0 else 0.5
        ry = (yc - y0) / (y1 - y0) if y1 != y0 else 0.5
        ax.set_xlim(xc - nw * rx, xc + nw * (1 - rx))
        ax.set_ylim(yc - nh * ry, yc + nh * (1 - ry))
        self.canvas.draw_idle()

    def _on_plot_double_click(self, event):
        if not getattr(event, "dblclick", False) or event.inaxes != self.ax or event.xdata is None:
            return
        ds = self.current_dataset()
        if not ds or not ds.peaks:
            return
        xspan = abs(self.ax.get_xlim()[1] - self.ax.get_xlim()[0])
        p = min(ds.peaks, key=lambda q: abs(q.two_theta - float(event.xdata)))
        if abs(p.two_theta - float(event.xdata)) > max(0.10, xspan * 0.008):
            return
        current = p.custom_label or self._peak_text(ds, p)
        text, ok = QInputDialog.getText(self, "Etykieta piku", f"Etykieta dla piku przy 2θ = {p.two_theta:.3f}°:", text=current)
        if ok:
            p.custom_label = text.strip()
            self.update_table(); self.redraw(); self.mark_project_dirty()

    def update_table(self):
        ds = self.current_dataset()
        self.table.blockSignals(True)
        self.table.setSortingEnabled(False)
        if not ds:
            self.table.setRowCount(0); self.table.setSortingEnabled(True); self.table.blockSignals(False); return
        self.table.setRowCount(len(ds.peaks))
        for r, p in enumerate(ds.peaks):
            vals = [
                ds.name, ds.material_type, f"{p.two_theta:.4f}", f"{p.d_A:.5f}", f"{p.height:.2f}",
                f"{p.rel_intensity:.2f}", f"{p.fwhm_deg:.4f}",
                "" if not math.isfinite(p.prominence) else f"{p.prominence:.2f}",
                "" if not math.isfinite(p.snr) else f"{p.snr:.2f}",
                (f"{p.phase} ({TENTATIVE_SUFFIX})" if p.phase_source == TENTATIVE_SOURCE else p.phase),
                p.hkl, p.custom_label, f"{ds.max_total_intensity:.2f}"
            ]
            for c, v in enumerate(vals):
                item = NumericSortItem(str(v)) if c in NUMERIC_PEAK_COLUMNS else QTableWidgetItem(str(v))
                if c == 2:
                    item.setData(Qt.ItemDataRole.UserRole, r)
                if c == 10:
                    src = str(getattr(p, "phase_source", "") or "").strip()
                    src_name = "zewnętrzna baza COD" if src == "cod" else "baza wewnętrzna" if src == "internal" else "nieprzypisane"
                    item.setToolTip("Źródło indeksów Millera: " + src_name)
                if c != 11:
                    item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(r, c, item)
        self.table.setSortingEnabled(True)
        self.table.blockSignals(False)

    def _table_item_changed(self, item: QTableWidgetItem):
        if self._loading_project or item.column() != 11:
            return
        ds = self.current_dataset()
        if not ds:
            return
        # Because sorting may rearrange rows, recover the peak from the 2θ item.
        tt_item = self.table.item(item.row(), 2)
        if tt_item is None:
            return
        try:
            tt = float(tt_item.text().replace(",", "."))
            p = min(ds.peaks, key=lambda q: abs(q.two_theta - tt))
            p.custom_label = item.text().strip()
            self.redraw(); self.mark_project_dirty()
        except Exception:
            return

    # --------------------------- export / database ---------------------------
    def save_plot(self):
        try:
            self._save_plot()
        except Exception as exc:
            QMessageBox.critical(
                self, "Zapis wykresu",
                f"Nie udało się zapisać wykresu:\n{exc}\n\n"
                "Jeśli plik jest otwarty w innym programie, zamknij go i spróbuj ponownie.")

    def _save_plot(self):
        ds = self.current_dataset()
        if not ds:
            QMessageBox.information(self, "Zapis wykresu", "Najpierw otwórz próbkę."); return
        path, filt = QFileDialog.getSaveFileName(
            self, "Zapisz wykres wysokiej rozdzielczości", self._default_save_path(f"{_safe_filename(ds.name)}_XRD.png"),
            "PNG (*.png);;TIFF (*.tif *.tiff);;PDF (wektorowy) (*.pdf);;SVG (wektorowy) (*.svg);;EPS (wektorowy) (*.eps)"
        )
        if not path:
            return
        self._remember_dir(path)
        dpi = int(self.dpi_combo.currentText())
        ext = Path(path).suffix.lower()
        if not ext:
            path += ".png"; ext = ".png"
        kwargs = {"bbox_inches": "tight"}
        if ext in (".png", ".tif", ".tiff"):
            kwargs["dpi"] = dpi
        self.fig.savefig(path, **kwargs)
        self._status(f"Zapisano wykres: {path}", 12000)

    def export_table(self):
        try:
            self._export_table()
        except Exception as exc:
            QMessageBox.critical(
                self, "Eksport",
                f"Nie udało się zapisać tabeli:\n{exc}\n\n"
                "Jeśli plik jest otwarty w Excelu lub innym programie, zamknij go i spróbuj ponownie.")

    def _export_table(self):
        ds = self.current_dataset()
        if not ds:
            QMessageBox.information(self, "Eksport", "Najpierw otwórz i przeanalizuj próbkę."); return
        path, selected_filter = QFileDialog.getSaveFileName(
            self, "Eksportuj tabelę pików", self._default_save_path(f"{_safe_filename(ds.name)}_peaks.xlsx"),
            "Excel (*.xlsx);;CSV (*.csv)")
        if not path:
            return
        self._remember_dir(path)
        if not Path(path).suffix:
            path += ".csv" if "csv" in str(selected_filter).lower() else ".xlsx"
        headers = [self.table.horizontalHeaderItem(c).text() for c in range(self.table.columnCount())]
        rows = [[self.table.item(r, c).text() if self.table.item(r, c) else "" for c in range(self.table.columnCount())] for r in range(self.table.rowCount())]
        ext = Path(path).suffix.lower()
        if ext == ".xlsx":
            if not OPENPYXL_OK:
                QMessageBox.warning(self, "Eksport", "Biblioteka openpyxl nie jest zainstalowana; zapisz zamiast tego jako CSV."); return
            # numbers must be real numbers in Excel (right-aligned, summable), not text
            rows = [[_number_or_text(v) if c in NUMERIC_PEAK_COLUMNS else v for c, v in enumerate(row)] for row in rows]
            wb = Workbook(); ws = wb.active; ws.title = "XRD peaks"; ws.append(headers)
            for row in rows:
                ws.append(row)
            meta_ws = wb.create_sheet("Metadane pomiaru")
            meta_ws.append(["Parametr", "Wartość", "Źródło / pochodzenie"])
            meta_rows = [
                ("Nazwa próbki", ds.name, "projekt/użytkownik"),
                ("Typ materiału", ds.material_type, "projekt/użytkownik"),
                ("Plik źródłowy", str(ds.path), "plik"),
                ("Format źródłowy", ds.metadata.get("source_format", ds.path.suffix), metadata_source(ds, "source_format")),
                ("Aparatura", ds.metadata.get("instrument", ""), metadata_source(ds, "instrument")),
                ("Producent", ds.metadata.get("manufacturer", ""), metadata_source(ds, "manufacturer")),
                ("Promieniowanie", ds.metadata.get("radiation", ""), metadata_source(ds, "radiation")),
                ("Anoda", ds.metadata.get("anode", ""), metadata_source(ds, "anode")),
                ("Długość fali [Å]", ds.wavelength_A, metadata_source(ds, "wavelength_A")),
                ("Długość fali Kα₂ [Å]", ds.metadata.get("wavelength_kalpha2_A", ""), metadata_source(ds,"wavelength_kalpha2_A")),
                ("Rekonstrukcja osi RD (anchor)", ds.metadata.get("rd_axis_anchor",""), metadata_source(ds,"rd_axis_missing_points")),
                ("Brakujące punkty nagłówka RD", ds.metadata.get("rd_axis_missing_points",""), metadata_source(ds,"rd_axis_missing_points")),
                ("Napięcie lampy [kV]", ds.metadata.get("tube_voltage_kV", ""), metadata_source(ds, "tube_voltage_kV")),
                ("Prąd lampy [mA]", ds.metadata.get("tube_current_mA", ""), metadata_source(ds, "tube_current_mA")),
                ("Początek skanu 2θ [°]", ds.scan_start_deg, "wyznaczone z danych dyfrakcyjnych"),
                ("Koniec skanu 2θ [°]", ds.scan_end_deg, "wyznaczone z danych dyfrakcyjnych"),
                ("Krok skanu [°]", ds.scan_step_deg, "wyznaczone z danych dyfrakcyjnych"),
                ("Szybkość skanowania [°/s]", ds.metadata.get("scan_rate_deg_s", ""), metadata_source(ds, "scan_rate_deg_s")),
                ("Geometria", ds.metadata.get("geometry", ""), metadata_source(ds, "geometry")),
                ("Data/czas pomiaru", ds.metadata.get("measurement_datetime", ""), metadata_source(ds, "measurement_datetime")),
                ("Uwagi", ds.metadata.get("notes", ""), metadata_source(ds, "notes")),
            ]
            for r in meta_rows:
                meta_ws.append(list(r))
            fits = ds.metadata.get("phase_fits", []) if isinstance(ds.metadata, dict) else []
            if isinstance(fits, list) and fits:
                fit_ws = wb.create_sheet("Dopasowanie sieci")
                fit_ws.append(["Prawdopodobna faza", "Model", "a [Å]", "c [Å]", "Dopasowane refleksy", "RMS Δ2θ [°]", "Wynik"])
                for f in fits:
                    if not isinstance(f, dict):
                        continue
                    fit_ws.append([
                        f.get("phase", ""), f.get("model", ""), f.get("a_A", ""), f.get("c_A", ""),
                        f.get("matched_reflections", ""), f.get("rms_2theta_deg", ""), f.get("score", "")
                    ])
            trace_ws = wb.create_sheet("Pochodzenie pikow")
            trace_headers, trace_rows = traceability_rows(ds)
            trace_ws.append(trace_headers)
            for trace_row in trace_rows:
                trace_ws.append(trace_row)
            cod_phase_ws = wb.create_sheet("Fazy COD")
            cod_phase_ws.append(["COD ID", "Faza", "Grupa przestrzenna", "Dopasowane refleksy",
                                 "RMS Δ2θ [°]", "Wynik rankingu [%]", "Tryb przypisania",
                                 "Mnożnik chemii [×]", "Pokrycie składu [%]"])
            seen_cod = set()
            for fit in (fits if isinstance(fits,list) else []):
                if not isinstance(fit,dict) or fit.get('cod_id') in (None,''):
                    continue
                if str(fit.get('cod_id')) in seen_cod:
                    continue
                seen_cod.add(str(fit.get('cod_id')))
                cod_phase_ws.append([fit.get('cod_id'),fit.get('phase',''),fit.get('space_group',''),
                  fit.get('matched_reflections',''), fit.get('rms_2theta_deg',''),
                  fit.get('score_percent',fit.get('score','')),
                  'ręcznie zatwierdzone' if fit.get('manual_selection') else 'automatyczna sugestia',
                  fit.get('chem_factor',''),fit.get('chemistry_capture_percent','')])
            for rec in (ds.metadata.get('manual_cod_peak_assignments',[]) or []):
                if not isinstance(rec,dict) or rec.get('cod_id') in (None,''):
                    continue
                if str(rec['cod_id']) not in seen_cod:
                    seen_cod.add(str(rec['cod_id']))
                    cod_phase_ws.append([rec['cod_id'],rec.get('phase',''),rec.get('space_group',''),
                                         '', '', '', 'ręcznie zatwierdzone', '', ''])
            wb.save(path)
        else:
            if not ext:
                path += ".csv"
            with open(path, "w", newline="", encoding="utf-8-sig") as f:
                w = csv.writer(f, delimiter=";"); w.writerow(headers); w.writerows(rows)
        self._status(f"Zapisano tabelę pików: {path}", 12000)

    def open_database(self):
        PhaseDatabaseDialog(self.db, self.cod_db, self.settings, self).exec()
        self.mark_project_dirty()

    def search_cod_measured_pattern(self):
        ds = self.current_dataset()
        if not ds:
            QMessageBox.information(self, "Wyszukiwanie wzorca COD", "Najpierw otwórz próbkę XRD.")
            return
        if not self.cod_db.is_open:
            ans = QMessageBox.question(
                self, "Wyszukiwanie wzorca COD",
                "Nie wczytano bazy COD .db3. Czy otworzyć teraz menedżer baz?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
            )
            if ans == QMessageBox.StandardButton.Yes:
                self.open_database()
            if not self.cod_db.is_open:
                return
        if not ds.peaks:
            self.analyze_current()
            if not ds.peaks:
                QMessageBox.information(self, "Wyszukiwanie wzorca COD", "Przy bieżących ustawieniach detekcji nie wykryto żadnych pików.")
                return
        dlg = CODPatternSearchDialog(ds, self.cod_db, self.db, self.phase_tol.value(), self)
        dlg.setWindowFlags(dlg.windowFlags() | Qt.WindowType.WindowMaximizeButtonHint | Qt.WindowType.WindowMinimizeButtonHint)
        dlg.exec()
        saved = str(ds.metadata.get("allowed_elements", "") or "").strip()
        if saved:
            self.allowed_elements_edit.setText(saved)
            if ds.peaks:
                try:
                    self._assign_phases_for_dataset(ds)
                except Exception as e:
                    clear_phase_assignments(ds)
                    if self.phase_source_combo.currentText() == "Tylko zewnętrzna baza COD":
                        QMessageBox.warning(self, "Identyfikacja zewnętrzna COD", str(e))
                self.update_table(); self.redraw(); self.update_summary()
        self.mark_project_dirty()

    def open_ai_phase_review(self):
        ds = self.current_dataset()
        if not ds:
            QMessageBox.information(self, "Ocena faz wspomagana AI", "Najpierw otwórz próbkę XRD.")
            return
        if not ds.peaks:
            self.analyze_current()
            if not ds.peaks:
                return
        try:
            build_ai_assisted_ranking(ds, self.phase_tol.value())
            dlg = AIAssistedReviewDialog(ds, self.phase_tol.value(), self)
            dlg.setWindowFlags(dlg.windowFlags() | Qt.WindowType.WindowMaximizeButtonHint | Qt.WindowType.WindowMinimizeButtonHint)
            dlg.exec()
            self.update_table(); self.redraw(); self.update_summary(); self.mark_project_dirty()
        except Exception as exc:
            QMessageBox.critical(self, "Ocena faz wspomagana AI", f"{type(exc).__name__}: {exc}\n\n{traceback.format_exc()}")

    def closeEvent(self, event):
        if self._plot_creator_dialog is not None:
            self._plot_creator_dialog._sync_to_project()
        if self.project_path is not None and self.project_dirty:
            try:
                self._write_project(self.project_path, silent=True)
            except Exception as exc:
                ans = QMessageBox.question(
                    self, "Nie zapisano projektu",
                    f"Zapis nie powiódł się:\n{exc}\n\nZamknąć program bez zapisania zmian?",
                    QMessageBox.StandardButton.Discard | QMessageBox.StandardButton.Cancel,
                    QMessageBox.StandardButton.Cancel,
                )
                if ans != QMessageBox.StandardButton.Discard:
                    event.ignore(); return
        elif self.project_path is None and self.datasets:
            ans = QMessageBox.question(
                self, "Zapisz projekt",
                "Czy przed zamknięciem zapisać bieżącą pracę XRD jako projekt?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No | QMessageBox.StandardButton.Cancel
            )
            if ans == QMessageBox.StandardButton.Cancel:
                event.ignore(); return
            if ans == QMessageBox.StandardButton.Yes:
                if not self.save_project_as():
                    event.ignore(); return
        self._save_window_geometry()
        event.accept()


def main():
    app = QApplication(sys.argv)
    install_exception_hook()
    w = MainWindow()
    w.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
