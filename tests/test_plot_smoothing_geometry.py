"""Presentation smoothing, project round trips and display geometry regressions."""
import copy
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

os.environ["QT_QPA_PLATFORM"] = "offscreen"
import numpy as np
from scipy.signal import savgol_filter
from PyQt6.QtCore import QSettings, Qt
from PyQt6.QtWidgets import QApplication, QCheckBox

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import XRD_Analyzer_GUI as xrd
from xrd_plot_support import display_smoothing_window, safe_window_rect, smooth_for_display

APP = QApplication.instance() or QApplication([])


class SmoothingNumericalTests(unittest.TestCase):
    def test_percentage_is_exact_mixture_and_inputs_unchanged(self):
        x = np.linspace(20, 22, 201)
        y = np.sin(x) + 0.03 * np.cos(100 * x)
        original_x, original_y = x.copy(), y.copy()
        x.setflags(write=False); y.setflags(write=False)
        full = smooth_for_display(x, y, 100)
        expected = savgol_filter(y, 25, 2, mode="interp")
        np.testing.assert_allclose(full, expected)
        np.testing.assert_allclose(smooth_for_display(x, y, 15), .85 * y + .15 * full)
        untouched = smooth_for_display(x, y, 0)
        np.testing.assert_array_equal(untouched, y)
        self.assertFalse(np.shares_memory(untouched, y))
        np.testing.assert_array_equal(x, original_x)
        np.testing.assert_array_equal(y, original_y)

    def test_irregular_axis_keeps_original_alignment_and_mixture(self):
        x = np.r_[np.linspace(20, 20.3, 101), np.linspace(20.31, 23, 80)]
        y = np.sin(x * 3) + np.cos(x * 120) / 10
        full = smooth_for_display(x, y, 100)
        self.assertEqual(full.shape, y.shape)
        self.assertTrue(np.all(np.isfinite(full)))
        np.testing.assert_allclose(smooth_for_display(x, y, 70), .3 * y + .7 * full)

    def test_short_arrays_and_invalid_inputs_are_stable(self):
        for n in range(0, 5):
            y = np.arange(n, dtype=float)
            np.testing.assert_array_equal(smooth_for_display(y, y, 100), y)
        y = np.arange(8, dtype=float)
        for x in (y[::-1], np.zeros(8), np.r_[np.nan, y[1:]]):
            np.testing.assert_array_equal(smooth_for_display(x, y, 100), y)
        y[2] = np.nan
        np.testing.assert_array_equal(smooth_for_display(np.arange(8), y, 100), y)

    def test_window_is_odd_bounded_and_explicit(self):
        self.assertEqual(display_smoothing_window(201, .01), 25)
        self.assertEqual(display_smoothing_window(8, .01), 7)
        for n in range(5, 60):
            for order in (1, 2, 3, 4):
                window = display_smoothing_window(n, .001, .25, order)
                self.assertLessEqual(window, n)
                self.assertGreater(window, order)
                self.assertEqual(window % 2, 1)


class GeometryLogicTests(unittest.TestCase):
    def assert_inside(self, rect, screen):
        x, y, w, h = rect; sx, sy, sw, sh = screen
        self.assertGreaterEqual(x, sx); self.assertGreaterEqual(y, sy)
        self.assertLessEqual(x + w, sx + sw); self.assertLessEqual(y + h, sy + sh)

    def test_removed_monitor_and_absurd_sizes_fit_available_area(self):
        screen = (0, 0, 1366, 728)
        for saved in ((2500, 800, 1500, 900), (0, 0, 20, 10), (0, 0, 20000, 9000),
                      None, "corrupt", [0, float('nan'), 1200, 800], [0, 0, -3, 500]):
            self.assert_inside(safe_window_rect(saved, [screen]), screen)

    def test_negative_monitor_origin_and_overlap_preserve_location(self):
        screens = [(-1920, 0, 1920, 1040), (0, 0, 1366, 728)]
        saved = (-1700, 100, 1200, 800)
        self.assertEqual(safe_window_rect(saved, screens), saved)

    def test_small_screen_overrides_normal_minimum(self):
        screen = (0, 20, 640, 460)
        self.assert_inside(safe_window_rect([0, 0, 200, 100], [screen]), screen)


class PlotSmoothingGuiTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='xrd-display-')
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        QSettings.setDefaultFormat(QSettings.Format.IniFormat)
        for scope in (QSettings.Scope.UserScope, QSettings.Scope.SystemScope):
            QSettings.setPath(QSettings.Format.IniFormat, scope, str(self.base / 'settings'))
        self.home_patch = patch.object(Path, 'home', return_value=self.base)
        self.home_patch.start(); self.addCleanup(self.home_patch.stop)
        self.window = xrd.MainWindow()
        x = np.linspace(20, 22, 201)
        for index in range(2):
            y = 5 + 10 * np.exp(-((x - 21) / .13) ** 2) + np.cos(x * 200)
            ds = xrd.Dataset(path=self.base / f'sample{index}.dat', name=f'Sample {index}',
                             material_type='', x=x.copy(), y=y.copy())
            ds.baseline = np.full_like(y, 5)
            ds.corrected_y = y - ds.baseline
            self.window.datasets.append(ds)
        self.dialogs = []

    def tearDown(self):
        self.window.project_path = None; self.window.project_dirty = False
        for dialog in self.dialogs:
            dialog.deleteLater()
        self.window.db.conn.close()
        self.window.deleteLater()
        APP.processEvents()

    def dialog(self, state=None):
        dialog = xrd.PlotCreatorDialog(self.window.datasets, self.window, state or {})
        self.dialogs.append(dialog)
        return dialog

    def test_checkbox_enables_only_own_strength_and_survives_reorder(self):
        creator = self.dialog()
        order = xrd.PlotOrderConfigDialog(creator)
        self.dialogs.append(order)
        self.assertFalse(order.table.cellWidget(0, 7).isEnabled())
        self.assertFalse(order.table.cellWidget(1, 7).isEnabled())
        cb = order.table.cellWidget(0, 6).findChild(QCheckBox)
        cb.setChecked(True)
        order.table.cellWidget(0, 7).setValue(70)
        self.assertTrue(order.table.cellWidget(0, 7).isEnabled())
        self.assertFalse(order.table.cellWidget(1, 7).isEnabled())
        order.table.selectRow(0); order._move(1); order._accept_changes()
        self.assertFalse(creator.rows[0]['smooth_enabled'])
        self.assertEqual(creator.rows[1]['smooth_strength_percent'], 70)
        self.assertIs(creator.rows[1]['dataset'], self.window.datasets[0])

    def test_method_selector_is_per_row_and_defaults_to_peak_preserving(self):
        creator = self.dialog()
        order = xrd.PlotOrderConfigDialog(creator)
        self.dialogs.append(order)
        combo = order.table.cellWidget(0, 8)
        self.assertEqual(combo.currentData(), 'peak_preserving')
        self.assertFalse(combo.isEnabled())
        order.table.cellWidget(0, 6).findChild(QCheckBox).setChecked(True)
        self.assertTrue(combo.isEnabled())
        self.assertFalse(order.table.cellWidget(1, 8).isEnabled())
        combo.setCurrentIndex(combo.findData('savgol'))
        order._accept_changes()
        self.assertEqual(creator.rows[0]['smooth_method'], 'savgol')
        self.assertEqual(creator.rows[1]['smooth_method'], 'peak_preserving')

    def test_stored_classic_method_reproduces_the_plain_filter(self):
        creator = self.dialog()
        creator.rows[0].update(smooth_enabled=True, smooth_strength_percent=100, smooth_method='savgol')
        creator.data_combo.setCurrentIndex(0)
        ds = self.window.datasets[0]
        np.testing.assert_allclose(creator._curve(ds), savgol_filter(ds.y, 25, 2, mode='interp'))

    def test_rendering_all_layouts_changes_only_display_copy(self):
        creator = self.dialog()
        snapshots = [[getattr(ds, k).copy() for k in ('x', 'y', 'baseline', 'corrected_y')]
                     for ds in self.window.datasets]
        original_peaks = copy.deepcopy([ds.peaks for ds in self.window.datasets])
        creator.rows[0].update(smooth_enabled=True, smooth_strength_percent=70)
        for mode in (0, 1, 2):
            creator.mode_combo.setCurrentIndex(mode)
            for data_mode in range(4):
                creator.data_combo.setCurrentIndex(data_mode)
                creator.redraw()
        for ds, arrays in zip(self.window.datasets, snapshots):
            for key, expected in zip(('x', 'y', 'baseline', 'corrected_y'), arrays):
                np.testing.assert_array_equal(getattr(ds, key), expected)
        self.assertEqual([ds.peaks for ds in self.window.datasets], original_peaks)
        creator.data_combo.setCurrentIndex(0)
        self.assertFalse(np.array_equal(creator._curve(self.window.datasets[0]), self.window.datasets[0].y))
        np.testing.assert_array_equal(creator._curve(self.window.datasets[1]), self.window.datasets[1].y)

    def test_save_restart_project_restores_smoothing_by_dataset_id(self):
        creator = self.dialog()
        creator.rows[0].update(smooth_enabled=True, smooth_strength_percent=70)
        creator.rows.reverse(); creator.redraw()
        expected = creator.get_state()
        project = self.base / 'smooth.xrdproj'
        self.window._write_project(project, silent=True)
        output = self.base / 'restored.json'
        helper = Path(__file__).with_name('test_plot_creator_persistence.py')
        # Existing restart helper validates the complete plot state and arrays.
        # Its result reader expects one peak in each sample.
        for ds in self.window.datasets:
            ds.peaks = [xrd.PeakResult(index=100, two_theta=21., d_A=2., height=10.,
                                      rel_intensity=100., fwhm_deg=.2)]
        self.window._write_project(project, silent=True)
        result = subprocess.run([sys.executable, str(helper), '--restore-project', str(project),
                                 str(output), str(self.base / 'new-home')], capture_output=True,
                                text=True, encoding='utf-8', errors='replace', timeout=90)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        restored = json.loads(output.read_text(encoding='utf-8'))
        self.assertEqual(restored['state'], expected)
        for actual, ds in zip(restored['datasets'], self.window.datasets):
            np.testing.assert_array_equal(actual['y'], ds.y)

    def test_bad_optional_fields_do_not_discard_valid_siblings(self):
        ds = self.window.datasets[0]
        state = {'title': 'keep this', 'mode': 'bad', 'width_in': None,
                 'legend_font_size': float('inf'), 'line_width': 2.4,
                 'rows': [{'dataset_id': ds.dataset_id, 'label': 'keep row', 'color': 'bad color',
                           'smooth_enabled': True, 'smooth_strength_percent': 70,
                           'smooth_polyorder': None}],
                 'separators': [{'position': float('inf'), 'offset_pct': 'bad', 'color': 'bad'}],
                 'label_positions': {'bad': {'coord': 'data', 'x': 'bad', 'y': 1}},
                 'manual_legend_loc': [float('nan'), 0]}
        restored = self.dialog(state).get_state()
        self.assertEqual(restored['title'], 'keep this')
        self.assertEqual(restored['line_width'], 2.4)
        self.assertEqual(restored['rows'][0]['label'], 'keep row')
        self.assertEqual(restored['rows'][0]['smooth_strength_percent'], 70)
        self.assertEqual(restored['rows'][0]['smooth_polyorder'], 2)
        self.assertEqual(restored['label_positions'], {})

    def test_old_projects_default_to_no_smoothing(self):
        creator = self.dialog({'rows': [{'dataset_id': self.window.datasets[0].dataset_id,
                                        'label': 'legacy'}]})
        for row in creator.rows:
            self.assertFalse(row['smooth_enabled'])
            self.assertEqual(row['smooth_strength_percent'], 0)

    def test_saved_geometry_and_maximized_state_are_local_and_validated(self):
        window = self.window
        window.settings.setValue('ui/main_geometry', 'corrupt')
        window.settings.setValue('ui/main_normal_geometry', [4000, 2000, 9000, 9000])
        window.settings.setValue('ui/main_maximized', False)
        window._restore_window_geometry()
        window.show(); APP.processEvents()
        rect = window.geometry()
        available = window.screen().availableGeometry()
        self.assertTrue(available.contains(rect))
        self.assertFalse(window.isMaximized())
        window._save_window_geometry()
        saved = window.settings.value('ui/main_normal_geometry')
        window._restore_window_geometry(); APP.processEvents()
        self.assertEqual(list(window.geometry().getRect()), saved)
        window.showMaximized(); APP.processEvents(); window._save_window_geometry()
        window.showNormal(); window._restore_window_geometry(); APP.processEvents()
        self.assertTrue(window.isMaximized())


if __name__ == '__main__':
    unittest.main()
