"""Scientific workspace: all fits consume original Dataset.x / Dataset.y."""
from __future__ import annotations

import csv
import json
import math
import re
from copy import deepcopy
from pathlib import Path

import numpy as np
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QFormLayout, QLabel,
    QPushButton, QTabWidget, QWidget, QDoubleSpinBox, QLineEdit, QCheckBox,
    QTableWidget, QTableWidgetItem, QHeaderView, QTextEdit, QFileDialog,
    QMessageBox, QInputDialog, QApplication)
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

from xrd_profile_analysis import (PeakFitResult, InstrumentBroadeningModel,
    fit_peak_profile, correct_instrument_broadening, scherrer_size,
    fit_williamson_hall, fit_cubic_lattice, _json_safe)
from xrd_advanced_analysis import fit_sin2psi, fit_wilkens_profiles
from xrd_calibration import calibrated_axis, calibration_status


PROFILE_HEADER_TIPS = {
    0: "Odznacz, aby pominąć refleks w dopasowaniu, WH i innych obliczeniach.",
    1: "Nazwa fazy. Do WH i sieci wszystkie użyte refleksy muszą należeć do jednej fazy.",
    2: "Indeksy Millera refleksu, np. 111 lub 2 0 0.",
    3: "Przybliżona pozycja refleksu 2θ (Kα1) – punkt startowy dopasowania.",
    4: "Początek zakresu dopasowania [°2θ]. Zakres obejmuje refleks (i dublet Kα) oraz tło po obu stronach.",
    5: "Koniec zakresu dopasowania [°2θ].",
    6: "Współczynnik kontrastu dyslokacji C_hkl (tylko dla modelu Wilkensa); w pozostałych obliczeniach może być pusty.",
    7: "Zmierzona szerokość połówkowa (FWHM) modelu pseudo-Voigta Kα1, przed korektą instrumentalną.",
    8: "Rozmiar domen koherentnego rozpraszania wg Scherrera (po korekcie IPF). To nie jest rozmiar ziarna.",
    9: "Status dopasowania: ACCEPTED/CONDITIONAL – można używać z zastrzeżeniami, REJECTED – wynik odrzucony.",
}

BUTTON_TIPS = {
    "Dodaj zakres": "Dodaj nowy zakres 2θ wokół refleksu, który chcesz dopasować.",
    "Usuń zakres": "Usuń zaznaczony wiersz tabeli.",
    "Dopasuj zaznaczone profile": "Dopasuj pseudo-Voigta Kα1/Kα2 z liniowym tłem do surowych danych w każdym zaznaczonym zakresie.",
    "Oblicz WH": "Williamson–Hall: wymaga ≥3 różnych hkl jednej fazy oraz wczytanego wzorca IPF.",
    "Parametr sieci kubicznej": "Dopasuj parametr sieci a i przesunięcie zera dla struktury kubicznej (wymaga hkl).",
    "IPF z wzorca Si": "Użyj profilu aparatury zmierzonego na wzorcu krzemowym tej serii (Analiza naukowa → Kalibracja z wzorca Si).",
    "Wczytaj wzorzec IPF": "Wczytaj profil instrumentalny (JSON lub CSV) zmierzony na wzorcu, np. LaB6.",
    "Eksportuj fity jako wzorzec": "Zapisz bieżące dopasowania wzorca jako plik IPF (min. dwa refleksy).",
    "Usuń wzorzec": "Usuń wczytany profil instrumentalny; korekta szerokości zostanie wyłączona.",
    "Dodaj pomiar": "Dodaj pusty wiersz do serii pomiarów sin²ψ.",
    "Usuń pomiar": "Usuń zaznaczony wiersz serii.",
    "Wczytaj CSV serii": "Wczytaj serię z CSV (kolumny: psi_deg, phi_deg, two_theta_deg, u_two_theta_deg, hkl, phase).",
    "Dodaj bieżący fit": "Dodaj pozycję refleksu zaznaczonego w zakładce profili.",
    "Zbierz hkl z próbek": "Zbierz dopasowania wybranego hkl ze wszystkich wczytanych próbek.",
}

STRESS_HEADERS = ("ψ [°]", "φ [°]", "2θ [°]", "u(2θ) [°]")
STRESS_FIELD_LABELS = {
    "xec_half_s2_per_MPa": "S2/2 [1/MPa]", "u_half_s2_per_MPa": "u(S2/2) [1/MPa]",
    "d0_A": "d0 [Å]", "u_d0_A": "u(d0) [Å]",
}
_DECIMAL_COMMA = re.compile(r"^\s*[-+]?\d+,\d+(?:[eE][-+]?\d+)?\s*$")


def number_spin(value, low, high, decimals=6, step=None):
    widget = QDoubleSpinBox()
    widget.setDecimals(decimals)
    widget.setRange(low, high)
    widget.setValue(value)
    # valueChanged only after editing is finished: every intermediate typed value
    # would otherwise invalidate all fits and results.
    widget.setKeyboardTracking(False)
    if step is not None:
        widget.setSingleStep(step)
    return widget


def parse_number(text, label, allow_empty=False):
    """Parse a user-typed number (decimal point or comma) with a readable Polish error."""
    cleaned = str(text).strip().replace(",", ".")
    if not cleaned:
        if allow_empty:
            return None
        raise ValueError(f"Uzupełnij pole: {label}.")
    try:
        value = float(cleaned)
    except ValueError:
        raise ValueError(f"{label}: „{str(text).strip()}” nie jest liczbą.") from None
    if not math.isfinite(value):
        raise ValueError(f"{label}: wartość musi być skończoną liczbą.")
    return value


def read_csv_records(path):
    """Read a CSV whose separator (; , or tab) is taken from the header line.

    Decimal commas are accepted in files that use ';' or tab as the separator.
    """
    text = Path(path).read_text(encoding="utf-8-sig")
    lines = text.splitlines()
    header = next((line for line in lines if line.strip()), "")
    counts = {delimiter: header.count(delimiter) for delimiter in (";", "\t", ",")}
    delimiter = max(counts, key=lambda key: counts[key])
    if counts[delimiter] == 0:
        raise ValueError("Nie rozpoznano separatora kolumn w pierwszym wierszu pliku "
                         "(oczekiwano średnika, przecinka lub tabulatora).")
    records = list(csv.DictReader(lines, delimiter=delimiter, skipinitialspace=True))
    if delimiter != ",":
        for record in records:
            for key, value in record.items():
                if isinstance(value, str) and _DECIMAL_COMMA.match(value):
                    record[key] = value.replace(",", ".")
    return records


def format_value(value, digits=7):
    if value is None:
        return "—"
    if isinstance(value, (float, np.floating)):
        return f"{value:.{digits}g}" if math.isfinite(value) else "—"
    if isinstance(value, (list, tuple)):
        return "[" + "; ".join(format_value(item, digits) for item in value) + "]"
    return str(value)


def format_with_uncertainty(value, uncertainty):
    """'value ± u' with the uncertainty rounded to two significant digits (no exponent notation)."""
    if value is None or not math.isfinite(float(value)):
        return "—"
    if uncertainty is None or not math.isfinite(float(uncertainty)) or float(uncertainty) <= 0:
        return format_value(float(value), 6)
    decimals = min(10, max(0, 1 - math.floor(math.log10(float(uncertainty)))))
    return f"{float(value):.{decimals}f} ± {float(uncertainty):.{decimals}f}"


class ScientificAnalysisDialog(QDialog):
    """Persist logical input/results inside dataset.metadata, with provenance."""

    def __init__(self, dataset, parent=None):
        super().__init__(parent)
        self.dataset = dataset
        self.setWindowTitle(f"Analiza naukowa — {dataset.name}")
        self.resize(1180, 860)
        self.setWindowFlags(self.windowFlags() | Qt.WindowType.WindowMaximizeButtonHint)
        old = dataset.metadata.get("scientific_analysis", {})
        self.state = deepcopy(old) if isinstance(old, dict) else {}
        self.state.setdefault("schema", 1)
        self.state.setdefault("fits", {})
        self.state.setdefault("results", {})
        self._building = True
        root = QVBoxLayout(self)
        intro = QLabel("Dopasowanie korzysta z surowych intensywności i własnego liniowego tła. "
            "Wybierz izolowane refleksy, zakresy i tę samą fazę do WH. "
            "Szerokości detekcyjne nie są używane w Scherrerze ani WH.")
        intro.setWordWrap(True)
        root.addWidget(intro)
        self.tabs = QTabWidget()
        root.addWidget(self.tabs, 1)
        self._build_profiles()
        self._build_stress()
        self._build_wilkens()
        diagnostics = QWidget()
        dl = QVBoxLayout(diagnostics)
        self.details = QTextEdit(); self.details.setReadOnly(True)
        dl.addWidget(self.details)
        self.tabs.addTab(diagnostics, "Diagnostyka i niepewności")
        row = QHBoxLayout()
        self.export_btn = QPushButton("Eksportuj wyniki JSON")
        self.csv_btn = QPushButton("Eksportuj tabelę profili CSV")
        self.close_btn = QPushButton("Zapisz w projekcie i zamknij")
        row.addWidget(self.export_btn); row.addWidget(self.csv_btn); row.addStretch(); row.addWidget(self.close_btn)
        root.addLayout(row)
        self.export_btn.clicked.connect(self.export_json)
        self.csv_btn.clicked.connect(self.export_csv)
        self.close_btn.clicked.connect(self.accept)
        self._building = False
        self.stress_table.itemChanged.connect(lambda *_: self._invalidate_advanced("sin2psi"))
        for widget in self.stress_fields.values():
            widget.textChanged.connect(lambda *_: self._invalidate_advanced("sin2psi"))
        self.stress_assumptions.toggled.connect(lambda *_: self._invalidate_advanced("sin2psi"))
        for widget in self.wilkens_fields.values():
            widget.textChanged.connect(lambda *_: self._invalidate_advanced("wilkens"))
        for button in self.findChildren(QPushButton):
            # Enter confirms a typed value; it must not also trigger e.g. "Dodaj zakres"
            button.setAutoDefault(False)
        try:
            self.use_standard_ipf(auto=True)         # series calibrated with a silicon standard: adopt its IPF
        except Exception:
            pass
        self._refresh_results()

    def _build_profiles(self):
        tab = QWidget(); layout = QVBoxLayout(tab)
        form = QHBoxLayout()
        settings = self.state.get("profile_settings", {})
        cu = abs(self.dataset.wavelength_A-1.5406) < .003
        self.wave1 = number_spin(settings.get("wavelength_A", 1.540598 if cu else self.dataset.wavelength_A), .01, 10, step=.0001)
        self.wave2 = number_spin(settings.get("wavelength_kalpha2_A", 1.544426 if cu else self.dataset.wavelength_A), .01, 10, step=.0001)
        self.ratio = number_spin(settings.get("kalpha2_ratio", .5 if cu else 0), 0, 1, 3, step=.01)
        self.shape_k = number_spin(settings.get("K", .9), .01, 3, 3, step=.01)
        for label, widget in (("λ1 [Å]", self.wave1), ("λ2 [Å]", self.wave2),
                              ("Kα2/Kα1", self.ratio), ("K Scherrera", self.shape_k)):
            form.addWidget(QLabel(label)); form.addWidget(widget)
        self.ratio.setToolTip("Stosunek pól Kα2/Kα1. Ustaw 0 dla linii monochromatycznej; sprawdź źródło promieniowania.")
        layout.addLayout(form)
        self.profile_table = QTableWidget(0, 10)
        self.profile_table.setHorizontalHeaderLabels(["Użyj", "Faza", "hkl", "2θ start [°]", "Od [°]", "Do [°]",
                                                     "C_hkl", "FWHM obs. [°]", "D Scherrer [nm]", "Status"])
        for column, tip in PROFILE_HEADER_TIPS.items():
            header_item = self.profile_table.horizontalHeaderItem(column)
            if header_item is not None:
                header_item.setToolTip(tip)
        self.profile_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.profile_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        records = self.state.get("regions")
        if not isinstance(records, list):
            records = []
            for p in self.dataset.peaks:
                half = max(.8, 4*p.width_half_prominence_deg)
                records.append(dict(use=True, phase=p.phase if p.phase != "unassigned" else "", hkl=p.hkl,
                    center=p.two_theta, low=max(float(self.dataset.x[0]), p.two_theta-half),
                    high=min(float(self.dataset.x[-1]), p.two_theta+half), contrast=""))
        for record in records:
            self._append_region(record)
        layout.addWidget(self.profile_table, 2)
        row = QHBoxLayout()
        for label, callback in (("Dodaj zakres", self.add_region), ("Usuń zakres", self.remove_region),
                                ("Dopasuj zaznaczone profile", self.fit_profiles),
                                ("Oblicz WH", self.calculate_wh), ("Parametr sieci kubicznej", self.calculate_lattice)):
            b = QPushButton(label); b.setToolTip(BUTTON_TIPS.get(label, "")); b.clicked.connect(callback); row.addWidget(b)
        layout.addLayout(row)
        instrument_row = QHBoxLayout()
        for label, callback in (("IPF z wzorca Si", self.use_standard_ipf),
                                ("Wczytaj wzorzec IPF", self.load_instrument),
                                ("Eksportuj fity jako wzorzec", self.export_instrument),
                                ("Usuń wzorzec", self.clear_instrument)):
            b = QPushButton(label); b.setToolTip(BUTTON_TIPS.get(label, "")); b.clicked.connect(callback); instrument_row.addWidget(b)
        self.instrument_label = QLabel(); self.instrument_label.setWordWrap(True)
        instrument_row.addWidget(self.instrument_label, 1)
        layout.addLayout(instrument_row)
        self.summary = QLabel(); self.summary.setWordWrap(True); layout.addWidget(self.summary)
        self.figure = Figure(figsize=(9, 3.5), layout="constrained")
        self.canvas = FigureCanvasQTAgg(self.figure); layout.addWidget(self.canvas, 3)
        self.tabs.addTab(tab, "Profile / Scherrer / WH")
        self.profile_table.itemChanged.connect(self._profile_edited)
        self.profile_table.currentCellChanged.connect(lambda *_: self.draw_profile())
        for widget in (self.wave1, self.wave2, self.ratio):
            widget.valueChanged.connect(self._radiation_edited)
        self.shape_k.valueChanged.connect(self._k_edited)

    def _append_region(self, rec):
        table = self.profile_table; row = table.rowCount(); table.insertRow(row)
        checkbox = QCheckBox(); checkbox.setChecked(bool(rec.get("use", True)))
        checkbox.toggled.connect(self._selection_changed)
        table.setCellWidget(row, 0, checkbox)
        for col, key in enumerate(("phase", "hkl", "center", "low", "high", "contrast"), 1):
            value = rec.get(key, "")
            if key in ("center", "low", "high") and isinstance(value, (int, float)):
                value = f"{value:.5f}"
            table.setItem(row, col, QTableWidgetItem(str(value)))
        for col in (7, 8, 9):
            item = QTableWidgetItem("—"); item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            table.setItem(row, col, item)

    def regions(self):
        result = []
        for row in range(self.profile_table.rowCount()):
            result.append(dict(use=self.profile_table.cellWidget(row, 0).isChecked(), **{
                key: self.profile_table.item(row, col).text().strip()
                for col, key in enumerate(("phase", "hkl", "center", "low", "high", "contrast"), 1)}))
        return result

    def add_region(self):
        center, ok = QInputDialog.getDouble(self, "Nowy zakres", "Przybliżone 2θ [°]:",
            float(np.mean(self.dataset.x)), float(self.dataset.x[0]), float(self.dataset.x[-1]), 4)
        if ok:
            self._building = True
            self._append_region(dict(use=True, center=center, low=max(float(self.dataset.x[0]), center-1),
                                     high=min(float(self.dataset.x[-1]), center+1)))
            self._building = False
            self._save_state()

    def remove_region(self):
        row = self.profile_table.currentRow()
        if row < 0: return
        self.profile_table.removeRow(row)
        self.state["fits"] = {str(int(key)-(int(key)>row)): val for key,val in self.state["fits"].items() if int(key)!=row}
        self.state["results"].clear(); self._save_state(); self._refresh_results()

    def _profile_edited(self, item):
        if self._building or item.column() >= 7: return
        if item.column() != 6:
            self.state["fits"].pop(str(item.row()), None)
        self.state["results"].clear(); self._save_state(); self._refresh_results()

    def _selection_changed(self):
        if self._building: return
        for key in ("williamson_hall", "cubic_lattice", "wilkens"):
            self.state["results"].pop(key, None)
        self._save_state(); self._refresh_results()

    def _invalidate_advanced(self, key):
        if self._building: return
        self.state["results"].pop(key, None)
        if key == "sin2psi":
            self.stress_summary.setText("Dane zmienione — oblicz ponownie.")
            self.stress_figure.clear(); self.stress_canvas.draw_idle()
        else:
            self.wilkens_summary.setText("Dane zmienione — oblicz ponownie.")
            self.wilkens_figure.clear(); self.wilkens_canvas.draw_idle()
        self._save_state()

    def _radiation_edited(self):
        if self._building: return
        self.state["fits"].clear(); self.state["results"].clear()
        self._save_state(); self._refresh_results()

    def _k_edited(self):
        if self._building: return
        self.state["results"].clear(); self._recalculate_corrections(); self._save_state(); self._refresh_results()

    def _profile_settings(self):
        return dict(wavelength_A=self.wave1.value(), wavelength_kalpha2_A=self.wave2.value(),
                    kalpha2_ratio=self.ratio.value(), K=self.shape_k.value())

    def fit_profiles(self):
        errors = []
        settings = self._profile_settings(); settings.pop("K")
        # With a silicon-standard calibration the fits use the corrected 2theta axis, so every position
        # (and the lattice, WH and stress results built on them) is corrected consistently.
        axis, calibration = calibrated_axis(self.dataset)
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            for row, rec in enumerate(self.regions()):
                if not rec["use"]: continue
                self.state["fits"].pop(str(row), None)
                try:
                    low = parse_number(rec["low"], "Od [°]"); high = parse_number(rec["high"], "Do [°]")
                    center = parse_number(rec["center"], "2θ start [°]")
                    if not low < center < high: raise ValueError("Zakres musi obejmować centrum refleksu.")
                    mask = (axis >= low) & (axis <= high)
                    fit = fit_peak_profile(axis[mask], self.dataset.y[mask], center_hint=center,
                                           hkl=rec["hkl"], **settings)
                    record = dict(fit=fit.to_dict(), phase=rec["phase"],
                                  source="raw Dataset.x/y" if calibration is None else "Dataset.y, oś 2θ skorygowana wzorcem Si",
                                  low=low, high=high)
                    if calibration is not None:
                        record["axis_correction"] = calibration.axis.describe()
                    self.state["fits"][str(row)] = record
                except Exception as exc:
                    errors.append(f"Wiersz {row+1}: {exc}")
            self.state["results"].clear()
            self._recalculate_corrections()
            self._save_state(); self._refresh_results()
            self.profile_table.selectRow(next((int(k) for k in self.state["fits"]), 0))
        finally:
            QApplication.restoreOverrideCursor()
        if errors: QMessageBox.warning(self, "Dopasowanie profili", "\n".join(errors))

    def _recalculate_corrections(self):
        model = self.state.get("instrument")
        instrument = InstrumentBroadeningModel(**model) if model else None
        for record in self.state["fits"].values():
            record.pop("corrected", None); record.pop("scherrer", None)
            if instrument is None: continue
            fit = PeakFitResult(**record["fit"])
            correction = correct_instrument_broadening(fit, instrument)
            record["corrected"] = correction.to_dict()
            if correction.fwhm_corrected_deg is not None:
                record["scherrer"] = scherrer_size(fit.two_theta_kalpha1_deg, correction.fwhm_corrected_deg,
                    fit.wavelength_A, self.shape_k.value(), correction.u_fwhm_corrected_deg,
                    fit.u_two_theta_deg).to_dict()

    def load_instrument(self):
        path, _ = QFileDialog.getOpenFileName(self, "Profil instrumentalny", "", "Wzorzec (*.json *.csv)")
        if not path: return
        try:
            if Path(path).suffix.lower() == ".json":
                values = json.loads(Path(path).read_text(encoding="utf-8-sig"))
            else:
                rows = read_csv_records(path)
                values = {key: [float(r[key]) for r in rows] for key in
                          ("two_theta_deg", "gaussian_fwhm_deg", "lorentzian_fwhm_deg")}
                for key in ("u_gaussian_fwhm_deg", "u_lorentzian_fwhm_deg"):
                    if all(r.get(key) for r in rows): values[key] = [float(r[key]) for r in rows]
                values.update(wavelength_A=self.wave1.value(), name=Path(path).stem, source_file=str(path))
            instrument = InstrumentBroadeningModel(**values)
            self.state["instrument"] = instrument.to_dict()
            self.state["results"].clear(); self._recalculate_corrections(); self._save_state(); self._refresh_results()
        except Exception as exc:
            QMessageBox.warning(self, "Wzorzec", str(exc))

    def use_standard_ipf(self, auto: bool = False):
        """Take the instrument profile measured on the silicon standard of this series."""
        calibration, reason = calibration_status(self.dataset)
        if calibration is None or not calibration.ipf:
            if not auto:
                QMessageBox.information(self, "Wzorzec Si", "Ta próbka nie ma profilu aparatury z wzorca. "
                                        "Wykonaj kalibrację: Analiza naukowa → Kalibracja z wzorca Si.")
            return
        if reason:
            if not auto:
                QMessageBox.information(self, "Wzorzec Si", f"Profilu wzorca nie można użyć: {reason}.")
            return
        if auto and self.state.get("instrument"):
            return
        self.state["instrument"] = deepcopy(calibration.ipf)
        self.state["results"].clear(); self._recalculate_corrections(); self._save_state(); self._refresh_results()

    def clear_instrument(self):
        self.state.pop("instrument", None); self.state["results"].clear()
        self._recalculate_corrections(); self._save_state(); self._refresh_results()

    def _selected_fits(self, corrected=False):
        selected = []
        for row, rec in enumerate(self.regions()):
            if not rec["use"]: continue
            record = self.state["fits"].get(str(row))
            if not record or record["fit"]["status"] == "REJECTED":
                raise ValueError(f"Wiersz {row+1}: wymagany zaakceptowany fit; odznacz nieużywany refleks.")
            if corrected and record.get("corrected", {}).get("fwhm_corrected_deg") is None:
                raise ValueError(f"Wiersz {row+1}: brak skorygowanej szerokości; wczytaj zgodny wzorzec IPF.")
            selected.append((rec, record))
        if not selected: raise ValueError("Zaznacz co najmniej jeden dopasowany refleks.")
        return selected

    def export_instrument(self):
        try:
            selected = self._selected_fits()
            if len(selected) < 2: raise ValueError("Wzorzec wymaga co najmniej dwóch refleksów obejmujących zakres próbki.")
            name, ok = QInputDialog.getText(self, "Wzorzec", "Nazwa wzorca / konfiguracja optyki:")
            if not ok or not name.strip(): return
            fits = sorted([r["fit"] for _,r in selected], key=lambda f: f["two_theta_kalpha1_deg"])
            values = dict(two_theta_deg=[f["two_theta_kalpha1_deg"] for f in fits],
                gaussian_fwhm_deg=[f["gaussian_fwhm_deg"] for f in fits],
                lorentzian_fwhm_deg=[f["lorentzian_fwhm_deg"] for f in fits],
                name=name, standard=name, optics_id=name, wavelength_A=self.wave1.value(), source_file=str(self.dataset.path))
            if all(f["covariance"] is not None for f in fits):
                values.update(u_gaussian_fwhm_deg=[math.sqrt(max(0,f["covariance"][2][2])) for f in fits],
                              u_lorentzian_fwhm_deg=[math.sqrt(max(0,f["covariance"][3][3])) for f in fits])
            model = InstrumentBroadeningModel(**values)
            path, _ = QFileDialog.getSaveFileName(self, "Zapisz wzorzec IPF", "instrument_profile.json", "JSON (*.json)")
            if path: Path(path).write_text(json.dumps(model.to_dict(), ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
        except Exception as exc: QMessageBox.warning(self, "Wzorzec", str(exc))

    def _same_phase(self, selected):
        phases = {rec["phase"].strip() for rec,_ in selected}
        if len(phases) != 1 or not next(iter(phases)):
            raise ValueError("Podaj i wybierz refleksy jednej fazy. Nie łącz różnych faz w jednym dopasowaniu.")
        if any(not rec["hkl"].strip() for rec,_ in selected): raise ValueError("Uzupełnij hkl wszystkich użytych refleksów.")
        return next(iter(phases))

    def calculate_wh(self):
        try:
            selected = self._selected_fits(corrected=True); self._same_phase(selected)
            widths = [r["corrected"]["fwhm_corrected_deg"] for _,r in selected]
            sigmas = [r["corrected"]["u_fwhm_corrected_deg"] for _,r in selected]
            result = fit_williamson_hall([r["fit"]["two_theta_kalpha1_deg"] for _,r in selected], widths,
                self.wave1.value(), self.shape_k.value(), sigmas if all(s is not None and s>0 for s in sigmas) else None,
                [rec["hkl"] for rec,_ in selected]).to_dict()
            self.state["results"]["williamson_hall"] = result
            self._save_state(); self._show_diagnostic(result)
            self.figure.clear(); ax, res = self.figure.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios":[3,1]})
            ax.errorbar(result["x"], result["y"], yerr=result["sigma_y"], fmt="o")
            x=np.asarray(result["x"]); ax.plot(x, result["intercept"]+result["epsilon"]*x)
            ax.set_ylabel("β cosθ [rad]"); res.plot(x,result["residuals"],"o"); res.axhline(0,color="gray")
            res.set(xlabel="4 sinθ",ylabel="Residua"); self.canvas.draw()
            self.summary.setText(f"WH: ε_app={format_value(result['epsilon'], 5)}; D_WH={format_value(result['D_nm'], 4)} nm. " + " ".join(result["warnings"]))
        except Exception as exc: QMessageBox.warning(self, "Williamson–Hall", str(exc))

    def calculate_lattice(self):
        try:
            selected = self._selected_fits(); self._same_phase(selected)
            sigmas = [r["fit"]["u_two_theta_deg"] for _,r in selected]
            result = fit_cubic_lattice([r["fit"]["two_theta_kalpha1_deg"] for _,r in selected],
                [rec["hkl"] for rec,_ in selected], self.wave1.value(), fit_zero_shift=True,
                sigma_2theta_deg=sigmas if all(s is not None and s>0 for s in sigmas) else None).to_dict()
            self.state["results"]["cubic_lattice"] = result; self._save_state(); self._show_diagnostic(result)
            self.summary.setText(f"Sieć kubiczna: a={format_with_uncertainty(result['a_A'], result['u_a_A'])} Å; "
                f"zero shift={format_value(result['zero_shift_deg'], 5)}°. " + " ".join(result["warnings"]))
        except Exception as exc: QMessageBox.warning(self, "Parametr sieci", str(exc))

    def _refresh_results(self):
        self._building = True
        for row in range(self.profile_table.rowCount()):
            record = self.state["fits"].get(str(row), {})
            fit, sch = record.get("fit", {}), record.get("scherrer", {})
            for col, value in ((7,fit.get("fwhm_obs_deg")),(8,sch.get("D_nm")),(9,fit.get("status"))):
                self.profile_table.item(row,col).setText(format_value(value, 5 if col == 7 else 4))
        self._building = False
        inst = self.state.get("instrument")
        self.instrument_label.setText("IPF: " + (inst.get("name","") if inst else "brak — Scherrer/WH wymagają wzorca"))
        self.summary.setText("Wybierz wiersz, aby obejrzeć profil i residua. Wyniki materiałowe są warunkowe.")
        self.draw_profile()

    def draw_profile(self):
        if self._building: return
        row = self.profile_table.currentRow()
        record = self.state["fits"].get(str(row))
        self.figure.clear()
        if record:
            fit=record["fit"]; x=np.asarray(fit["x_deg"]); y=np.asarray(fit["fitted_y"]); residual=np.asarray(fit["residuals"])
            ax, res = self.figure.subplots(2,1,sharex=True,gridspec_kw={"height_ratios":[3,1]})
            ax.plot(x,y+residual,".",ms=3,label="Pomiar surowy"); ax.plot(x,y,label="TCH / Kα")
            ax.legend(); ax.set_ylabel("Intensywność"); res.plot(x,residual); res.axhline(0,color="gray")
            res.set(xlabel="2θ [°]",ylabel="Residua")
            self.summary.setText(f"FWHM={format_with_uncertainty(fit['fwhm_obs_deg'], fit['u_fwhm_obs_deg'])}° (1σ); "
                f"R²={format_value(fit['r_squared'])}; RMSE={format_value(fit['rmse'], 4)}. " + " ".join(fit["warnings"]))
            self._show_diagnostic(record)
        else:
            ax=self.figure.add_subplot(111); ax.plot(self.dataset.x,self.dataset.y,lw=.6)
            ax.set(xlabel="2θ [°]",ylabel="Intensywność surowa")
        self.canvas.draw_idle()

    def _build_stress(self):
        tab=QWidget(); layout=QVBoxLayout(tab)
        label=QLabel("Seria tego samego hkl i fazy przy jednym φ. Potrzebne ≥3 niezależne sin²ψ, "
                     "rentgenowskie S2/2 [1/MPa] oraz niezależne d0 [Å]. ±ψ służą kontroli rozszczepienia.")
        label.setWordWrap(True); layout.addWidget(label)
        self.stress_table=QTableWidget(0,7)
        self.stress_table.setHorizontalHeaderLabels(["ψ [°]","φ [°]","2θ [°]","u(2θ) [°]","hkl","Faza","Próbka"])
        self.stress_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.stress_keys=("psi_deg","phi_deg","two_theta_deg","u_two_theta_deg","hkl","phase","sample")
        for rec in self.state.get("stress_measurements",[]): self._append_stress(rec)
        layout.addWidget(self.stress_table,1)
        row=QHBoxLayout()
        for title, callback in (("Dodaj pomiar",lambda:self._append_stress({})),("Usuń pomiar",self.remove_stress),
                               ("Wczytaj CSV serii",self.load_stress), ("Dodaj bieżący fit",self.stress_from_fit),
                               ("Zbierz hkl z próbek", self.stress_from_datasets)):
            b=QPushButton(title); b.setToolTip(BUTTON_TIPS.get(title, "")); b.clicked.connect(callback); row.addWidget(b)
        layout.addLayout(row)
        settings=self.state.get("stress_settings",{})
        form=QFormLayout(); self.stress_fields={}
        for key,label,default in (("xec_half_s2_per_MPa","S2/2 [1/MPa]",""),("u_half_s2_per_MPa","u(S2/2) [1/MPa]","0"),
                                 ("d0_A","d0 [Å]",""),("u_d0_A","u(d0) [Å]","0"),
                                 ("xec_hkl","hkl XEC",""),("xec_phase","Faza XEC",""),("xec_source","Źródło XEC","")):
            widget=QLineEdit(str(settings.get(key,default))); self.stress_fields[key]=widget; form.addRow(label,widget)
        self.stress_assumptions=QCheckBox("Potwierdzam założenie płaskiego stanu naprężenia i zgodność XEC z fazą/hkl")
        self.stress_assumptions.setChecked(bool(settings.get("assumptions_confirmed",False)))
        form.addRow(self.stress_assumptions); layout.addLayout(form)
        b=QPushButton("Oblicz naprężenie sin²ψ"); b.clicked.connect(self.calculate_stress); layout.addWidget(b)
        self.stress_summary=QLabel("Brak wyniku."); self.stress_summary.setWordWrap(True); layout.addWidget(self.stress_summary)
        self.stress_figure=Figure(figsize=(8,3),layout="constrained"); self.stress_canvas=FigureCanvasQTAgg(self.stress_figure)
        layout.addWidget(self.stress_canvas,1); self.tabs.addTab(tab,"Naprężenia sin²ψ")

    def _append_stress(self,rec):
        blocked=self.stress_table.blockSignals(True)
        row=self.stress_table.rowCount(); self.stress_table.insertRow(row)
        for col,key in enumerate(self.stress_keys):
            value=rec.get(key)
            self.stress_table.setItem(row,col,QTableWidgetItem("" if value is None else str(value)))
        self.stress_table.blockSignals(blocked)
        if not self._building:self._invalidate_advanced("sin2psi")

    def remove_stress(self):
        if self.stress_table.currentRow()>=0:
            self.stress_table.removeRow(self.stress_table.currentRow()); self._invalidate_advanced("sin2psi")

    def load_stress(self):
        path,_=QFileDialog.getOpenFileName(self,"Seria ψ","","CSV (*.csv)")
        if not path:return
        try:
            records=read_csv_records(path)
            if not all(all(key in rec for key in self.stress_keys[:6]) for rec in records):
                raise ValueError("Kolumny CSV: "+", ".join(self.stress_keys[:6]))
            self.stress_table.setRowCount(0)
            for rec in records:self._append_stress(rec)
        except Exception as exc:QMessageBox.warning(self,"Seria ψ",str(exc))

    def stress_from_fit(self):
        record=self.state["fits"].get(str(self.profile_table.currentRow()))
        if not record:
            QMessageBox.information(self,"Seria ψ","Zaznacz w zakładce „Profile” wiersz z dopasowanym refleksem.")
            return
        f=record["fit"]
        self._append_stress(dict(two_theta_deg=f["two_theta_kalpha1_deg"],u_two_theta_deg=f["u_two_theta_deg"],
                                hkl=f["hkl"],phase=record["phase"],sample=self.dataset.name))

    def stress_from_datasets(self):
        from xrd_profile_analysis import _parse_hkl
        text, ok=QInputDialog.getText(self,"Seria ψ","hkl do zebrania, np. 311:",text=self.stress_fields["xec_hkl"].text())
        if not ok:return
        try:
            target=_parse_hkl(text); count=0
            self._save_state()
            for dataset in getattr(self.parent(),"datasets",[self.dataset]):
                for record in dataset.metadata.get("scientific_analysis",{}).get("fits",{}).values():
                    fit=record["fit"]
                    if fit.get("status")=="REJECTED" or not fit.get("hkl"):continue
                    if _parse_hkl(fit["hkl"])!=target:continue
                    if not math.isclose(fit["wavelength_A"],self.wave1.value(),rel_tol=1e-6):continue
                    self._append_stress(dict(psi_deg=dataset.metadata.get("psi_deg",""),phi_deg=dataset.metadata.get("phi_deg",""),
                        two_theta_deg=fit["two_theta_kalpha1_deg"],u_two_theta_deg=fit["u_two_theta_deg"],
                        hkl=fit["hkl"],phase=record["phase"],sample=dataset.name)); count+=1
            self.stress_summary.setText(f"Dodano {count} dopasowanych pozycji. Uzupełnij ψ i φ, sprawdź fazę oraz usuń powtórzenia.")
        except Exception as exc:QMessageBox.warning(self,"Seria ψ",str(exc))

    def _stress_input(self,parse=False):
        measurements=[]
        for row in range(self.stress_table.rowCount()):
            rec={key:self.stress_table.item(row,col).text().strip() for col,key in enumerate(self.stress_keys)}
            if parse:
                for col,key in enumerate(self.stress_keys[:4]):
                    rec[key]=parse_number(rec[key],f"Pomiar {row+1}, kolumna „{STRESS_HEADERS[col]}”")
                rec["wavelength_A"]=self.wave1.value()
            measurements.append(rec)
        settings={key:widget.text().strip() for key,widget in self.stress_fields.items()}
        settings["assumptions_confirmed"]=self.stress_assumptions.isChecked()
        if parse:
            for key in ("xec_half_s2_per_MPa","u_half_s2_per_MPa","d0_A","u_d0_A"):
                value=parse_number(settings[key],STRESS_FIELD_LABELS[key],allow_empty=True)
                settings[key]=value if value is not None else (None if key in ("xec_half_s2_per_MPa","d0_A") else 0)
        return measurements,settings

    def calculate_stress(self):
        try:
            measurements,settings=self._stress_input(True)
            result=fit_sin2psi(measurements,**settings)
            self.state["results"]["sin2psi"]=result; self._save_state(); self._show_diagnostic(result)
            self.stress_summary.setText(f"{result['status']}: σ={format_value(result['stress_MPa'])} MPa; "
                f"95% CI={format_value(result.get('stress_ci95_MPa'))}. "+" ".join(result["warnings"]))
            self.stress_figure.clear()
            if "sin2psi" in result:
                ax,res=self.stress_figure.subplots(2,1,sharex=True,gridspec_kw={"height_ratios":[3,1]})
                x=np.asarray(result["sin2psi"]); order=np.argsort(x)
                ax.errorbar(x,result["d_A"],yerr=result["u_d_A"],fmt="o")
                ax.plot(x[order],np.asarray(result["fitted_d_A"])[order]); ax.set_ylabel("d [Å]")
                res.plot(x,result["residuals_A"],"o"); res.axhline(0,color="gray");res.set(xlabel="sin²ψ",ylabel="Residua [Å]")
            self.stress_canvas.draw()
        except Exception as exc:QMessageBox.warning(self,"sin²ψ",str(exc))

    def _build_wilkens(self):
        tab=QWidget(); layout=QVBoxLayout(tab)
        info=QLabel("Eksperymentalny model całych, rozdzielonych profili Wilkensa: ρ, Re, lognormalny rozkład domen. "
            "To ograniczony model, wymagający porównania z CMWP przed interpretacją bezwzględnej ρ. "
            "W zakładce profili zaznacz ≥3 różne hkl jednej fazy, wpisz C_hkl i wczytaj IPF.")
        info.setWordWrap(True); layout.addWidget(info)
        settings=self.state.get("wilkens_settings",{})
        form=QFormLayout();self.wilkens_fields={}
        for key,label,default in (("burgers_vector_nm","Wektor Burgersa b [nm]",""),("structure","Struktura (np. FCC)",""),
                                 ("instrument_source","Identyfikacja wzorca / optyki",""),
                                 ("contrast_source","Źródło C_hkl",""),("slip_system","Układ poślizgu","")):
            widget=QLineEdit(str(settings.get(key,default)));self.wilkens_fields[key]=widget;form.addRow(label,widget)
        layout.addLayout(form)
        row=QHBoxLayout()
        b=QPushButton("Dopasuj model Wilkensa");b.clicked.connect(self.calculate_wilkens);row.addWidget(b)
        b=QPushButton("Eksportuj dane whole-profile JSON");b.clicked.connect(self.export_whole_profile);row.addWidget(b)
        layout.addLayout(row)
        self.wilkens_summary=QLabel("Brak wyniku. Program nie wyznacza ρ z 1/D².")
        self.wilkens_summary.setWordWrap(True);layout.addWidget(self.wilkens_summary)
        self.wilkens_figure=Figure(figsize=(9,5),layout="constrained");self.wilkens_canvas=FigureCanvasQTAgg(self.wilkens_figure)
        layout.addWidget(self.wilkens_canvas,1);self.tabs.addTab(tab,"Wilkens — whole-profile")

    def _whole_profile_input(self):
        selected=self._selected_fits();phase=self._same_phase(selected)
        if not self.state.get("instrument"):raise ValueError("Wczytaj profil instrumentalny IPF.")
        instrument=InstrumentBroadeningModel(**self.state["instrument"])
        if not math.isclose(instrument.wavelength_A,self.wave1.value(),rel_tol=1e-6):raise ValueError("Niezgodna λ wzorca.")
        profiles=[]
        for rec,record in selected:
            fit=record["fit"];mask=(self.dataset.x>=record["low"])&(self.dataset.x<=record["high"])
            g,l,_,_=instrument.at(fit["two_theta_kalpha1_deg"])
            profiles.append(dict(x_deg=self.dataset.x[mask].tolist(),y=self.dataset.y[mask].tolist(),
                hkl=rec["hkl"],phase=phase,center_deg=fit["two_theta_kalpha1_deg"],contrast_factor=parse_number(rec["contrast"],f"C_hkl dla refleksu {rec['hkl']}"),
                instrument_gaussian_fwhm_deg=g,instrument_lorentzian_fwhm_deg=l))
        settings={key:widget.text().strip() for key,widget in self.wilkens_fields.items()}
        settings["burgers_vector_nm"]=parse_number(settings["burgers_vector_nm"],"Wektor Burgersa b [nm]",allow_empty=True)
        settings.update(phase=phase,wavelength_A=self.wave1.value(),kalpha2_wavelength_A=self.wave2.value(),kalpha2_ratio=self.ratio.value())
        return profiles,settings

    def calculate_wilkens(self):
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            profiles,settings=self._whole_profile_input()
            result=fit_wilkens_profiles(profiles,**settings)
            self.state["results"]["wilkens"]=result;self._save_state();self._show_diagnostic(result)
            self.wilkens_summary.setText(f"{result['status']}: ρ={format_value(result.get('rho_m2'))} m⁻²; "
                f"95% CI={format_value(result.get('rho_ci95_m2'))}. "+" ".join(result["warnings"]))
            self.wilkens_figure.clear()
            records=result.get("profiles",[])
            if records:
                axes=self.wilkens_figure.subplots(len(records),2,squeeze=False)
                for rec,aa in zip(records,axes):
                    aa[0].plot(rec["x_deg"],rec["observed_y"],".",ms=2);aa[0].plot(rec["x_deg"],rec["fitted_y"])
                    aa[0].set_title(str(rec["hkl"]));aa[1].plot(rec["x_deg"],rec["residuals_y"]);aa[1].axhline(0,color="gray")
                    aa[0].set_xlabel("2θ [°]");aa[1].set_title("Residua")
            self.wilkens_canvas.draw()
        except Exception as exc:QMessageBox.warning(self,"Wilkens",str(exc))
        finally:QApplication.restoreOverrideCursor()

    def export_whole_profile(self):
        try:
            profiles,settings=self._whole_profile_input()
            path,_=QFileDialog.getSaveFileName(self,"Dane whole-profile","whole_profile.json","JSON (*.json)")
            if path:Path(path).write_text(json.dumps(dict(schema="xrd-analyzer-whole-profile/1",profiles=profiles,settings=settings),
                ensure_ascii=False,indent=2,allow_nan=False),encoding="utf-8")
        except Exception as exc:QMessageBox.warning(self,"Whole-profile",str(exc))

    def _show_diagnostic(self,result):
        if hasattr(self,"details"):self.details.setPlainText(json.dumps(_json_safe(result),ensure_ascii=False,indent=2,allow_nan=False))

    def _save_state(self):
        if self._building:return
        self.state["regions"]=self.regions();self.state["profile_settings"]=self._profile_settings()
        self.state["stress_measurements"],self.state["stress_settings"]=self._stress_input()
        self.state["wilkens_settings"]={key:widget.text().strip() for key,widget in self.wilkens_fields.items()}
        self.state["dataset_id"]=self.dataset.dataset_id
        value=_json_safe(self.state)
        if self.dataset.metadata.get("scientific_analysis")!=value:
            self.dataset.metadata["scientific_analysis"]=deepcopy(value)
            if self.parent() is not None:self.parent().mark_project_dirty()

    def export_json(self):
        self._save_state()
        path,_=QFileDialog.getSaveFileName(self,"Eksport wyników","scientific_analysis.json","JSON (*.json)")
        if path:
            try:Path(path).write_text(json.dumps(_json_safe(self.state),ensure_ascii=False,indent=2,allow_nan=False),encoding="utf-8")
            except Exception as exc:QMessageBox.warning(self,"Eksport",str(exc))

    def export_csv(self):
        path,_=QFileDialog.getSaveFileName(self,"Tabela profili","profile_results.csv","CSV (*.csv)")
        if not path:return
        try:
            with open(path,"w",encoding="utf-8-sig",newline="") as stream:
                writer=csv.writer(stream,delimiter=";")
                writer.writerow(["phase","hkl","two_theta_deg","u_two_theta_deg","d_A","u_d_A","fwhm_obs_deg",
                    "u_fwhm_obs_deg","fwhm_instr_deg","fwhm_corrected_deg","u_fwhm_corrected_deg","D_nm","u_D_nm","status","warnings"])
                for rec in self.state["fits"].values():
                    f=rec["fit"];c=rec.get("corrected",{});s=rec.get("scherrer",{})
                    writer.writerow([rec["phase"],f["hkl"],f["two_theta_kalpha1_deg"],f["u_two_theta_deg"],f["d_A"],f["u_d_A"],
                        f["fwhm_obs_deg"],f["u_fwhm_obs_deg"],c.get("fwhm_instr_deg"),c.get("fwhm_corrected_deg"),
                        c.get("u_fwhm_corrected_deg"),s.get("D_nm"),s.get("u_D_nm"),f["status"],
                        " | ".join(f["warnings"]+c.get("warnings",[])+s.get("warnings",[]))])
        except Exception as exc:QMessageBox.warning(self,"Eksport",str(exc))

    def done(self,result):
        self._save_state()
        super().done(result)
