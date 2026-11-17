"""Behavioral project/Qt persistence tests; no real preferences or user data touched.

Run: python -m unittest discover -s tests -p test_plot_creator_persistence.py -v
XRD_ANALYZER_TEST_SOURCE optionally selects an unchanged baseline module.
"""
import copy
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import zipfile

os.environ["QT_QPA_PLATFORM"] = "offscreen"
import numpy as np
from PyQt6.QtCore import QSettings, QSignalBlocker, Qt
from PyQt6.QtTest import QTest
from PyQt6.QtWidgets import QApplication, QCheckBox, QComboBox, QDoubleSpinBox, QLineEdit, QSpinBox

SOURCE = Path(os.environ.get("XRD_ANALYZER_TEST_SOURCE", str(Path(__file__).resolve().parents[1] / "XRD_Analyzer_GUI.py")))
sys.path.insert(0, str(SOURCE.parent))
spec = importlib.util.spec_from_file_location("xrd_persistence_under_test", SOURCE)
xrd = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = xrd
spec.loader.exec_module(xrd)
APP = QApplication.instance() or QApplication([])


def isolate_settings(directory):
    QSettings.setDefaultFormat(QSettings.Format.IniFormat)
    for scope in (QSettings.Scope.UserScope, QSettings.Scope.SystemScope):
        QSettings.setPath(QSettings.Format.IniFormat, scope, str(directory))


def make_datasets(directory):
    result = []
    x = np.linspace(20.0, 80.0, 121)
    for i in range(3):
        y = 5.0 + (100.0 + i * 10) * np.exp(-((x - 40.0) / 1.3) ** 2)
        source = directory / f"sample-{i}.dat"
        np.savetxt(source, np.column_stack((x, y)))
        ds = xrd.Dataset(path=source, name=f"Próbka {i}", material_type="", x=x.copy(), y=y)
        ds.baseline = np.full_like(x, 5.0)
        ds.corrected_y = y - ds.baseline
        ds.peaks = [xrd.PeakResult(index=40, two_theta=40.0, d_A=2.25, height=float(y[40]),
                                  rel_intensity=100.0, fwhm_deg=0.3, custom_label=f"Pik {i}")]
        result.append(ds)
    return result


def configure(dialog, mode=1):
    """Set every exposed configuration family, including actual annotation drag."""
    widgets = [
        w for w in dialog.findChildren(xrd.QWidget)
        if isinstance(w, (QCheckBox, QComboBox, QDoubleSpinBox, QLineEdit, QSpinBox))]
    blockers = [QSignalBlocker(w) for w in widgets]
    dialog.mode_combo.setCurrentIndex(mode)
    dialog.data_combo.setCurrentIndex(3)
    dialog.offset_spin.setValue(1.72)
    dialog.keep_peak_height_cb.setChecked(False)
    dialog.linewidth_spin.setValue(2.15)
    dialog.legend_cb.setChecked(True)
    dialog.legend_pos_combo.setCurrentText("Ręcznie – przeciągnij w podglądzie")
    dialog.legend_font_spin.setValue(11)
    dialog.sample_names_cb.setChecked(False)
    dialog.grid_cb.setChecked(False)
    dialog.peaks_cb.setChecked(True)
    dialog.peak_labels_cb.setChecked(True)
    dialog.peak_label_style_combo.setCurrentIndex(0)
    dialog.hkl_display_combo.setCurrentIndex(1)
    dialog.peak_label_position_combo.setCurrentText("Po lewej stronie piku")
    dialog.peak_label_orientation_combo.setCurrentText("Poziomo")
    dialog.peak_font_spin.setValue(10)
    dialog.phase_marker_strip_cb.setChecked(True)
    dialog.title_edit.setText("Test trwałości — α/β")
    dialog.subtitle_edit.setText("Cu Kα, λ = 0,15406 nm")
    dialog.xlabel_edit.setText("Kąt 2θ [°]")
    dialog.ylabel_edit.setText("Natężenie względne")
    dialog.font_family_combo.setCurrentText("DejaVu Serif")
    for name, value in (("title_font_spin", 16), ("axis_font_spin", 13),
                        ("tick_font_spin", 12), ("sample_font_spin", 11)):
        getattr(dialog, name).setValue(value)
    dialog.width_spin.setValue(8.3)
    dialog.height_spin.setValue(5.2)
    dialog.dpi_combo.setCurrentText("900")
    del blockers
    dialog.rows[:] = [dialog.rows[2], dialog.rows[0], dialog.rows[1]]
    for i, row in enumerate(dialog.rows):
        row.update(label=f"Opis {i} — α", color=["#123456", "#cc5500", "#009955"][i], enabled=i != 2)
    dialog.separators = [xrd.normalizuj_separator({
        "position": 1, "top_label": "Referencja", "bottom_label": "Badane", "offset_pct": 14.5,
        "color": "#664422", "labels_visible": True}, len(dialog.rows))]
    dialog.redraw()
    for legend in dialog._legend_artists:
        legend.set_loc((0.23, 0.62))
    dialog.redraw()
    if dialog._annotation_entries:
        entry = dialog._annotation_entries[0]
        ann = entry["annotation"]
        old_x, old_y = ann.xyann
        ann.xyann = (float(old_x) + 1.1, float(old_y))
        dialog._drag_info = entry
        dialog._on_mouse_release(SimpleNamespace())
    return copy.deepcopy(dialog.get_state())


def project_payload(path):
    with zipfile.ZipFile(path) as archive:
        return json.loads(archive.read("project.json"))


def read_in_fresh_process(project, output, home):
    """Worker entry point: entirely new QApplication/MainWindow and no cached state."""
    home.mkdir(parents=True, exist_ok=True)
    isolate_settings(home)
    with patch.object(Path, "home", return_value=home), \
         patch.object(xrd.QMessageBox, "information"), \
         patch.object(xrd.QMessageBox, "warning") as warning, \
         patch.object(xrd.QMessageBox, "critical") as critical:
        window = xrd.MainWindow()
        with patch.object(xrd.QFileDialog, "getOpenFileName", return_value=(str(project), "")):
            window.open_project()
        if critical.called or warning.called:
            raise AssertionError(f"Unexpected open error: {critical.call_args}, {warning.call_args}")
        dialog = xrd.PlotCreatorDialog(window.datasets, window, window.plot_creator_state)
        result = {"state": dialog.get_state(), "datasets": [
            {"name": d.name, "x": d.x.tolist(), "y": d.y.tolist(), "baseline": d.baseline.tolist(),
             "corrected_y": d.corrected_y.tolist(), "label": d.peaks[0].custom_label}
            for d in window.datasets]}
        output.write_text(json.dumps(result, ensure_ascii=False), encoding="utf-8")
        dialog.deleteLater()
        window.project_path = None
        window.db.conn.close()
        window.deleteLater()


class PlotCreatorPersistenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="xrd-plot-persistence-")
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        isolate_settings(self.base / "settings")
        self.home_patch = patch.object(Path, "home", return_value=self.base)
        self.home_patch.start()
        self.addCleanup(self.home_patch.stop)
        self.messages = []
        for name in ("information", "warning", "critical"):
            p = patch.object(xrd.QMessageBox, name)
            self.messages.append(p.start())
            self.addCleanup(p.stop)
        self.window = xrd.MainWindow()
        self.window.datasets = make_datasets(self.base)
        self.window.sample_combo.addItems([d.name for d in self.window.datasets])
        self.window.project_dirty = False
        self.dialogs = []

    def tearDown(self):
        self.window.project_path = None
        self.window.project_dirty = False
        for dialog in self.dialogs:
            dialog.deleteLater()
        self.window.db.conn.close()
        self.window.deleteLater()
        APP.processEvents()

    def dialog(self, state=None):
        dialog = xrd.PlotCreatorDialog(self.window.datasets, self.window,
                                       self.window.plot_creator_state if state is None else state)
        self.dialogs.append(dialog)
        return dialog

    def modal(self, callback):
        def execute(dialog):
            self.dialogs.append(dialog)
            callback(dialog)
            return dialog.result()
        with patch.object(xrd.PlotCreatorDialog, "exec", execute):
            self.window.open_plot_creator()
        self.assertFalse(self.messages[2].called, str(self.messages[2].call_args))

    def restore_fresh(self, project):
        output = self.base / "fresh-result.json"
        result = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--restore-project",
                                 str(project), str(output), str(self.base / "fresh-home")],
                                text=True, capture_output=True, encoding="utf-8", errors="replace", timeout=90)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return json.loads(output.read_text(encoding="utf-8"))

    def test_all_nondefaults_survive_restore_in_each_layout(self):
        for mode in (0, 1, 2):
            with self.subTest(mode=mode):
                expected = configure(self.dialog({}), mode)
                restored = self.dialog(expected)
                self.assertEqual(restored.get_state(), expected)
                self.assertTrue(expected["label_positions"])

    def test_legend_drag_is_captured_without_redraw(self):
        dialog = self.dialog({})
        dialog.legend_pos_combo.setCurrentText("Ręcznie – przeciągnij w podglądzie")
        self.assertTrue(dialog._legend_artists)
        dialog._legend_artists[0].set_loc((0.17, 0.43))
        self.assertEqual(dialog.get_state()["manual_legend_loc"], [0.17, 0.43])

    def test_individual_panel_legends_survive_reopen(self):
        dialog = self.dialog({})
        configure(dialog, mode=2)
        positions = [(0.11, 0.23), (0.61, 0.47)]
        self.assertEqual(len(dialog._legend_artists), len(positions))
        for legend, position in zip(dialog._legend_artists, positions):
            legend.set_loc(position)
        expected = dialog.get_state()
        dialog.reject()
        restored = self.dialog()
        self.assertEqual(restored.get_state(), expected)
        self.assertEqual([tuple(leg._loc) for leg in restored._legend_artists], positions)

    def test_duplicate_names_and_rename_keep_distinct_label_positions(self):
        for ds in self.window.datasets:
            ds.name = "Ta sama nazwa"
        dialog = self.dialog({})
        configure(dialog)
        entries = dialog._annotation_entries
        self.assertEqual(len(entries), 2)
        self.assertNotEqual(entries[0]["key"], entries[1]["key"])
        for index, entry in enumerate(entries):
            x, y = entry["annotation"].xyann
            entry["annotation"].xyann = (float(x) + index + 1, float(y))
            dialog._drag_info = entry
            dialog._on_mouse_release(SimpleNamespace())
        expected = dialog.get_state()["label_positions"]
        dialog.accept()
        for ds in self.window.datasets:
            ds.name = "Nowa nazwa"
        restored = self.dialog()
        self.assertEqual(restored.get_state()["label_positions"], expected)
        for entry in restored._annotation_entries:
            self.assertEqual(tuple(entry["annotation"].xyann),
                             (expected[entry["key"]]["x"], expected[entry["key"]]["y"]))

    def test_open_close_without_changes_does_not_autosave(self):
        self.window._write_project(self.base / "nochange.xrdproj", silent=True)
        with patch.object(self.window, "_write_project") as write:
            self.modal(lambda dialog: dialog.reject())
        write.assert_not_called()
        self.assertFalse(self.window.project_dirty)

    def test_cancel_main_close_after_failed_save_retains_state(self):
        from PyQt6.QtGui import QCloseEvent
        self.window._write_project(self.base / "close-error.xrdproj", silent=True)
        dialog = self.dialog({})
        dialog.title_edit.setText("Niezapisana praca")
        event = QCloseEvent()
        with patch.object(self.window, "_write_project", side_effect=OSError("disk full")), \
             patch.object(xrd.QMessageBox, "question", return_value=xrd.QMessageBox.StandardButton.Cancel):
            self.window.closeEvent(event)
        self.assertFalse(event.isAccepted())
        self.assertTrue(self.window.project_dirty)
        self.assertEqual(self.window.plot_creator_state["title"], "Niezapisana praca")

    def test_close_reopen_accept_reject_window_close_and_escape(self):
        for variant in ("accept", "reject", "close", "escape"):
            with self.subTest(variant=variant):
                expected = []
                self.window.plot_creator_state = {}
                def change_and_close(dialog):
                    expected.append(configure(dialog))
                    if variant == "escape":
                        dialog.show()
                        QTest.keyClick(dialog, Qt.Key.Key_Escape)
                    else:
                        getattr(dialog, variant)()
                self.modal(change_and_close)
                restored = []
                self.modal(lambda d: (restored.append(d.get_state()), d.accept()))
                self.assertEqual(restored[0], expected[0])

    def test_save_while_dialog_active_includes_uncommitted_text_and_dpi(self):
        project = self.base / "active.xrdproj"
        expected = []
        def save_active(dialog):
            configure(dialog)
            dialog.title_edit.setText("Ostatni wpis bez utraty fokusu")
            dialog.dpi_combo.setCurrentText("1200")
            expected.append(copy.deepcopy(dialog.get_state()))
            with patch.object(xrd.QFileDialog, "getSaveFileName", return_value=(str(project), "")):
                self.assertTrue(self.window.save_project_as())
            self.assertEqual(project_payload(project)["plot_creator"], expected[0])
            dialog.reject()
        self.modal(save_active)
        self.assertEqual(self.restore_fresh(project)["state"], expected[0])

    def test_save_close_restart_open_restores_state_and_embedded_arrays(self):
        project = self.base / "restart.xrdproj"
        expected = []
        self.modal(lambda d: (expected.append(configure(d)), d.accept()))
        with patch.object(xrd.QFileDialog, "getSaveFileName", return_value=(str(project), "")):
            self.assertTrue(self.window.save_project_as())
        for dataset in self.window.datasets:
            dataset.path.unlink()
        result = self.restore_fresh(project)
        self.assertEqual(result["state"], expected[0])
        for actual, dataset in zip(result["datasets"], self.window.datasets):
            for field in ("x", "y", "baseline", "corrected_y"):
                np.testing.assert_array_equal(actual[field], getattr(dataset, field))
            self.assertEqual(actual["label"], dataset.peaks[0].custom_label)

    def test_save_as_and_relocation_preserve_row_identity(self):
        original = self.base / "original.xrdproj"
        expected = []
        self.modal(lambda d: (expected.append(configure(d)), d.accept()))
        self.window._write_project(original, silent=True)
        other = self.base / "saved-as"
        other.mkdir()
        saved_as = other / "renamed.xrdproj"
        with patch.object(xrd.QFileDialog, "getSaveFileName", return_value=(str(saved_as), "")):
            self.assertTrue(self.window.save_project_as())
        # Relocate the whole source/project folder, then remove absolute source paths
        # from the copied archive to exercise relative resolution and embedded arrays.
        moved = self.base / "moved"
        moved.mkdir()
        moved_project = moved / "relocated.xrdproj"
        with zipfile.ZipFile(saved_as) as src, zipfile.ZipFile(moved_project, "w") as dst:
            payload = json.loads(src.read("project.json"))
            for index, record in enumerate(payload["datasets"]):
                record["path"] = f"relocated-source-{index}.xy"
                record["absolute_path"] = ""
            for name in src.namelist():
                dst.writestr(name, json.dumps(payload, ensure_ascii=False) if name == "project.json" else src.read(name))
        for dataset in self.window.datasets:
            dataset.path.unlink()
        restored = self.restore_fresh(moved_project)["state"]
        expected_semantic = copy.deepcopy(expected[0])
        restored_semantic = copy.deepcopy(restored)
        for state in (expected_semantic, restored_semantic):
            for row in state["rows"]:
                row.pop("path", None)
        self.assertEqual(restored_semantic, expected_semantic)

    def test_project_without_plot_section_has_defaults(self):
        project = self.base / "legacy-no-plot.xrdproj"
        self.window._write_project(project, silent=True)
        payload = project_payload(project)
        payload.pop("plot_creator", None)
        legacy = self.base / "legacy-no-plot.json"
        legacy.write_text(json.dumps(payload), encoding="utf-8")
        expected = self.dialog({}).get_state()
        result = self.restore_fresh(legacy)
        actual = result["state"]
        # New optional per-dataset identifiers may be initialized on first use.
        for state in (actual, expected):
            for row in state["rows"]:
                row.pop("dataset_id", None)
        self.assertEqual(actual, expected)

    def test_legacy_english_settings_and_separator_migrate(self):
        legacy = {"mode": 1, "legend_position": "Upper left", "peak_label_position": "Below peak",
                  "peak_label_orientation": "Horizontal", "group_separator": True,
                  "group_break_after": 1, "top_group_label": "Reference", "bottom_group_label": "Test"}
        state = self.dialog(legacy).get_state()
        self.assertEqual(state["legend_position"], "Lewy górny")
        self.assertEqual(state["peak_label_position"], "Pod pikiem")
        self.assertEqual(state["peak_label_orientation"], "Poziomo")
        self.assertEqual(state["separators"][0]["top_label"], "Reference")
        self.assertEqual(self.dialog(state).get_state(), state)

    def test_failed_save_keeps_dirty_state_and_previous_archive(self):
        project = self.base / "write-failure.xrdproj"
        self.window._write_project(project, silent=True)
        previous_bytes = project.read_bytes()
        dialog = self.dialog({})
        with patch.object(xrd.os, "replace", side_effect=OSError("simulated full disk")):
            dialog.title_edit.setText("Niezapisana zmiana")
            dialog.title_edit.editingFinished.emit()
            dialog.accept()
            self.window.save_project()
        self.assertTrue(self.window.project_dirty)
        self.assertEqual(project.read_bytes(), previous_bytes)
        self.assertFalse(project.with_suffix(".xrdproj.tmp").exists())
        self.assertTrue(self.messages[2].called)
        self.assertEqual(self.window.plot_creator_state["title"], "Niezapisana zmiana")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--restore-project":
        read_in_fresh_process(*(Path(arg) for arg in sys.argv[2:5]))
    else:
        unittest.main()
