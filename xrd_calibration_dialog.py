"""Dialog: calibrate the measurement series with the silicon standard scanned together with the samples."""
from __future__ import annotations

import math

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (QApplication, QCheckBox, QComboBox, QDialog, QDoubleSpinBox, QFormLayout, QHBoxLayout,
                             QHeaderView, QLabel, QMessageBox, QPushButton, QTableWidget, QTableWidgetItem,
                             QVBoxLayout)

from xrd_calibration import (MODEL_LABELS, ROLE_KEY, SI_LATTICE_A, StandardReflection, calibrate_from_reflections,
                             fit_standard_reflections, install_calibration, reapply_calibration, remove_calibration)
from xrd_rd_import import KNOWN_LINES

INTRO = ("Wzorzec krzemowy zeskanowany razem z próbkami służy do poprawienia ich wyników. Z położeń refleksów Si "
         "program wyznacza błąd osi 2θ aparatu (zero, opcjonalnie przesunięcie próbki) oraz szerokość linii aparatury "
         "(IPF). Położenia pików i odstępy d we wszystkich próbkach o tej samej długości fali są poprawiane, a profil "
         "aparatury można od razu użyć w analizie naukowej (Scherrer, Williamson–Hall). Człon cosθ (przesunięcie próbki) "
         "przenosi się na próbki tylko wtedy, gdy wzorzec i próbki były zamocowane na tej samej wysokości — domyślnie "
         "stosowane jest samo zero.")


class CalibrationDialog(QDialog):
    HEADERS = ["Użyj", "hkl", "2θ teoret. [°]", "2θ zmierz. [°]", "Δ [°]", "u(2θ) [°]", "FWHM [°]", "G / L [°]", "R²", "Uwagi"]

    def __init__(self, datasets, current_index: int = 0, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Kalibracja z wzorca krzemowego (Si)")
        self.resize(1120, 860)
        self.datasets = datasets
        self.reflections: list = []
        self.calibration = None
        self._building = False
        layout = QVBoxLayout(self)
        intro = QLabel(INTRO)
        intro.setWordWrap(True)
        layout.addWidget(intro)

        form = QFormLayout()
        form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
        self.standard_combo = QComboBox()
        start = current_index if 0 <= current_index < len(datasets) else 0
        for index, ds in enumerate(datasets):
            marked = ds.metadata.get(ROLE_KEY) == "standard"
            self.standard_combo.addItem(f"{ds.name}{'  ✓ wzorzec kalibracyjny' if marked else ''}")
            if marked:
                start = index
        self.standard_combo.setCurrentIndex(start)
        self.standard_combo.setToolTip("Skan krzemu zmierzony razem z próbkami.")
        self.wave1_label = QLabel()
        self.wave2 = self._spin(1.544426, 0.01, 10.0, 6, 1e-4)
        self.ratio = self._spin(0.5, 0.0, 1.0, 3, 0.01)
        self.a_spin = self._spin(SI_LATTICE_A, 5.0, 5.9, 6, 1e-5)
        self.certified = QCheckBox("Parametr sieci pochodzi z certyfikatu wzorca")
        self.model_combo = QComboBox()
        for key, label in MODEL_LABELS.items():
            self.model_combo.addItem(label, key)
        self.model_combo.setCurrentIndex(self.model_combo.findData("zero"))
        self.model_combo.setToolTip("Zero: poprawka stała. Zero + cosθ: dodatkowo przesunięcie próbki wzorca (patrz opis u góry). "
                                    "Automatyczny: wybór kryterium AICc.")
        form.addRow("Wzorzec (skan Si):", self.standard_combo)
        form.addRow("Długość fali Kα1:", self.wave1_label)
        form.addRow("Długość fali Kα2 [Å]:", self.wave2)
        form.addRow("Stosunek Kα2/Kα1:", self.ratio)
        form.addRow("Parametr sieci Si a [Å]:", self.a_spin)
        form.addRow("", self.certified)
        form.addRow("Model korekty 2θ:", self.model_combo)
        layout.addLayout(form)
        self.fit_btn = QPushButton("Dopasuj refleksy wzorca")
        self.fit_btn.setToolTip("Dopasowuje profile wszystkich refleksów Si leżących w zakresie skanu.")
        layout.addWidget(self.fit_btn)

        self.table = QTableWidget(0, len(self.HEADERS))
        self.table.setHorizontalHeaderLabels(self.HEADERS)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(len(self.HEADERS) - 1, QHeaderView.ResizeMode.Stretch)
        self.table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        layout.addWidget(self.table, 2)
        self.figure = Figure(figsize=(9, 4.2), layout="constrained")
        self.canvas = FigureCanvasQTAgg(self.figure)
        layout.addWidget(self.canvas, 3)
        self.summary = QLabel("Wybierz skan wzorca i naciśnij „Dopasuj refleksy wzorca”.")
        self.summary.setWordWrap(True)
        layout.addWidget(self.summary)

        buttons = QHBoxLayout()
        self.apply_btn = QPushButton("Zastosuj do próbek")
        self.apply_btn.setToolTip("Zapisuje kalibrację w projekcie i poprawia położenia pików oraz d we wszystkich próbkach "
                                  "o tej samej długości fali.")
        self.remove_btn = QPushButton("Usuń kalibrację ze wszystkich")
        self.close_btn = QPushButton("Zamknij")
        for button in (self.apply_btn, self.remove_btn):
            buttons.addWidget(button)
        buttons.addStretch(1)
        buttons.addWidget(self.close_btn)
        layout.addLayout(buttons)
        self.apply_btn.setEnabled(False)

        self.standard_combo.currentIndexChanged.connect(self._standard_changed)
        self.fit_btn.clicked.connect(self.fit)
        self.model_combo.currentIndexChanged.connect(self._recompute)
        self.apply_btn.clicked.connect(self.apply)
        self.remove_btn.clicked.connect(self.remove)
        self.close_btn.clicked.connect(self.accept)
        for button in self.findChildren(QPushButton):
            button.setAutoDefault(False)
        self._standard_changed()

    @staticmethod
    def _spin(value, low, high, decimals, step):
        widget = QDoubleSpinBox()
        widget.setDecimals(decimals)
        widget.setRange(low, high)
        widget.setSingleStep(step)
        widget.setValue(value)
        widget.setKeyboardTracking(False)
        return widget

    # ------------------------------------------------------------------ helpers
    def _standard(self):
        return self.datasets[self.standard_combo.currentIndex()]

    def _standard_changed(self, *_):
        if not self.datasets:
            return
        ds = self._standard()
        self.wave1_label.setText(f"{ds.wavelength_A:.6f} Å (z danych wzorca)")
        for _anode, (ka1, ka2) in KNOWN_LINES.items():
            if abs(ka1 - ds.wavelength_A) <= 3e-4:
                self.wave2.setValue(ka2)
                break
        else:
            self.wave2.setValue(ds.wavelength_A * 1.0025)
        self.reflections, self.calibration = [], None
        self.table.setRowCount(0)
        self.apply_btn.setEnabled(False)

    def _refresh_parent(self):
        parent = self.parent()
        for name in ("update_table", "redraw", "update_summary", "mark_project_dirty"):
            method = getattr(parent, name, None)
            if callable(method):
                method()

    # ------------------------------------------------------------------ fitting
    def fit(self):
        ds = self._standard()
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            self.reflections = fit_standard_reflections(ds.x, ds.y, ds.wavelength_A, self.wave2.value(), self.ratio.value(),
                                                        self.a_spin.value())
        except Exception as exc:
            self.reflections = []
            QMessageBox.warning(self, "Kalibracja", f"Nie udało się dopasować refleksów wzorca:\n{exc}")
        finally:
            QApplication.restoreOverrideCursor()
        self._fill_table()
        self._recompute()

    def _fill_table(self):
        self._building = True
        self.table.setRowCount(len(self.reflections))
        for row, r in enumerate(self.reflections):
            check = QTableWidgetItem()
            check.setFlags(Qt.ItemFlag.ItemIsUserCheckable | Qt.ItemFlag.ItemIsEnabled)
            check.setCheckState(Qt.CheckState.Checked if r.used else Qt.CheckState.Unchecked)
            self.table.setItem(row, 0, check)
            values = [r.hkl, f"{r.two_theta_theory:.4f}", f"{r.two_theta_obs:.4f}", f"{r.delta:+.4f}", f"{r.u_two_theta:.4f}",
                      f"{r.fwhm_obs:.4f}", f"{r.gaussian_fwhm:.3f} / {r.lorentzian_fwhm:.3f}",
                      "—" if r.r_squared is None else f"{r.r_squared:.4f}", r.note]
            for column, text in enumerate(values, start=1):
                self.table.setItem(row, column, QTableWidgetItem(text))
        self._building = False
        try:
            self.table.itemChanged.disconnect(self._item_changed)
        except TypeError:
            pass
        self.table.itemChanged.connect(self._item_changed)

    def _item_changed(self, item):
        if not self._building and item.column() == 0:
            self._recompute()

    def _used_reflections(self):
        rows = []
        for row, r in enumerate(self.reflections):
            item = self.table.item(row, 0)
            checked = item is not None and item.checkState() == Qt.CheckState.Checked
            usable = math.isfinite(r.two_theta_obs) and math.isfinite(r.fwhm_obs)
            rows.append(StandardReflection(**{**r.__dict__, "used": bool(checked and usable)}))
        return rows

    def _recompute(self, *_):
        self.calibration = None
        self.apply_btn.setEnabled(False)
        if not self.reflections:
            self.figure.clear()
            self.canvas.draw_idle()
            return
        ds = self._standard()
        rows = self._used_reflections()
        try:
            self.calibration = calibrate_from_reflections(
                rows, ds.wavelength_A, self.wave2.value(), self.ratio.value(), self.a_spin.value(),
                self.certified.isChecked(), self.model_combo.currentData() or "auto", ds.name, str(ds.path))
        except ValueError as exc:
            self.summary.setText(f"Kalibracja niemożliwa: {exc}")
            self._draw(rows)
            return
        axis = self.calibration.axis
        self.summary.setText(
            f"Model: {MODEL_LABELS.get(axis.model, axis.model)}.  {axis.describe()}\n"
            + "\n".join("• " + note for note in self.calibration.notes))
        self.apply_btn.setEnabled(True)
        self._draw(rows)

    def _draw(self, rows):
        self.figure.clear()
        ax1, ax2 = self.figure.subplots(1, 2)
        used = [r for r in rows if r.used]
        other = [r for r in rows if not r.used and math.isfinite(r.two_theta_obs)]
        if used:
            ax1.errorbar([r.two_theta_obs for r in used], [r.delta for r in used], yerr=[r.u_two_theta for r in used],
                         fmt="o", color="tab:blue", label="użyte refleksy Si")
        if other:
            ax1.plot([r.two_theta_obs for r in other], [r.delta for r in other], "x", color="tab:red", label="pominięte")
        if self.calibration is not None:
            axis = self.calibration.axis
            tt = np.linspace(axis.two_theta_min - 3, axis.two_theta_max + 3, 200)
            ax1.plot(tt, axis.delta(tt), color="tab:green", label="model Δ(2θ)")
            band = axis.u_delta(tt)
            if np.all(np.isfinite(band)):
                ax1.fill_between(tt, axis.delta(tt) - band, axis.delta(tt) + band, color="tab:green", alpha=0.2)
        ax1.axhline(0.0, color="0.6", lw=0.8)
        ax1.set_xlabel("2θ [°]")
        ax1.set_ylabel("Δ = zmierzone − teoretyczne [°]")
        ax1.set_title("Błąd osi 2θ aparatu")
        ax1.legend(fontsize=8)
        if used:
            ax2.plot([r.two_theta_obs for r in used], [r.fwhm_obs for r in used], "o-", color="tab:purple")
        ax2.set_xlabel("2θ [°]")
        ax2.set_ylabel("FWHM linii wzorca [°]")
        ax2.set_title("Poszerzenie aparaturowe")
        self.canvas.draw_idle()

    # ------------------------------------------------------------------ apply / remove
    def apply(self):
        if self.calibration is None:
            return
        standard = self._standard()
        result = install_calibration(self.datasets, standard, self.calibration)
        parent = self.parent()
        reassigned = []
        for ds in self.datasets:
            if ds is standard:
                continue
            reapply_calibration(ds)
            had_phases = any(p.phase not in ("", "unassigned") for p in ds.peaks)
            if had_phases and ds.metadata.get("allowed_elements") and hasattr(parent, "_assign_phases_for_dataset"):
                try:
                    parent._assign_phases_for_dataset(ds)
                    reassigned.append(ds.name)
                except Exception:
                    pass
        self._refresh_parent()
        text = [f"Wzorzec: {standard.name}.", f"Zastosowano do próbek: {len(result['applied'])}."]
        if result["skipped"]:
            text.append("Pominięto (inna długość fali):\n  " + "\n  ".join(result["skipped"]))
        text.append("Próbki, które nie miały jeszcze wykrytych pików, zostaną skorygowane podczas analizy (F5 / Shift+F5).")
        if reassigned:
            text.append("Ponownie zidentyfikowano fazy: " + ", ".join(reassigned) + ".")
        QMessageBox.information(self, "Kalibracja", "\n".join(text))

    def remove(self):
        answer = QMessageBox.question(self, "Kalibracja", "Usunąć kalibrację wzorcem Si ze wszystkich próbek?\n"
                                      "Położenia pików wrócą do wartości zmierzonych.")
        if answer != QMessageBox.StandardButton.Yes:
            return
        count = remove_calibration(self.datasets)
        self._refresh_parent()
        self.summary.setText(f"Usunięto kalibrację ({count} wpisów).")
