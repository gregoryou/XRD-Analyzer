"""Usability regressions of the main window (offscreen Qt; settings/home redirected to a temp folder)."""
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

os.environ['QT_QPA_PLATFORM'] = 'offscreen'
import numpy as np
from PyQt6.QtCore import Qt
from test_plot_creator_persistence import xrd, APP, isolate_settings


class UsabilityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='xrd-usability-')
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        isolate_settings(self.base / 'settings')
        home = patch.object(Path, 'home', return_value=self.base)
        home.start(); self.addCleanup(home.stop)
        self.messages = {}
        for method in ('information', 'warning', 'critical', 'question'):
            mock = patch.object(xrd.QMessageBox, method)
            self.messages[method] = mock.start(); self.addCleanup(mock.stop)
        self.messages['question'].return_value = xrd.QMessageBox.StandardButton.No
        self.window = xrd.MainWindow()

    def tearDown(self):
        self.window.project_path = None; self.window.project_dirty = False
        self.window.db.conn.close(); self.window.deleteLater(); APP.processEvents()

    def add_dataset(self, two_thetas=(9.5, 10.2, 100.0)):
        x = np.linspace(5, 120, 1000)
        y = np.ones_like(x)
        ds = xrd.Dataset(self.base / 'peaks.dat', 'Próbka testowa', '', x, y)
        ds.peaks = [xrd.PeakResult(index=int(np.argmin(abs(x - t))), two_theta=t, d_A=xrd.bragg_d(t, ds.wavelength_A),
                                   height=h, rel_intensity=h, fwhm_deg=0.1) for t, h in zip(two_thetas, (50.0, 9.0, 100.0))]
        self.window.datasets = [ds]
        self.window.sample_combo.addItem(ds.name)
        self.window.update_table()
        return ds

    def test_menu_bar_has_groups_and_shortcuts(self):
        titles = [action.text() for action in self.window.menuBar().actions()]
        for expected in ("&Plik", "&Analiza", "Analiza naukowa", "&Wykres", "&Bazy faz", "Pomo&c"):
            self.assertIn(expected, titles)
        shortcuts = {}
        for top in self.window.menuBar().actions():
            for action in top.menu().actions():
                if not action.isSeparator():
                    shortcuts[action.text()] = action.shortcut().toString()
        self.assertEqual(shortcuts["Zapisz projekt"], "Ctrl+S")
        self.assertEqual(shortcuts["Otwórz projekt…"], "Ctrl+O")
        self.assertEqual(shortcuts["Analizuj bieżący"], "F5")
        self.assertIn("Profile / Scherrer / WH / sin²ψ / Wilkens…", shortcuts)

    def test_numeric_fields_commit_after_editing_not_per_keystroke(self):
        for name in ('wavelength', 'prominence', 'min_distance', 'min_snr', 'max_peaks', 'phase_tol', 'label_top_n'):
            with self.subTest(field=name):
                self.assertFalse(getattr(self.window, name).keyboardTracking())

    def test_peak_table_sorts_numbers_numerically(self):
        self.add_dataset()
        self.window.table.sortItems(2, Qt.SortOrder.AscendingOrder)
        self.assertEqual([self.window.table.item(r, 2).text() for r in range(3)], ["9.5000", "10.2000", "100.0000"])
        self.window.table.sortItems(4, Qt.SortOrder.DescendingOrder)
        self.assertEqual([self.window.table.item(r, 4).text() for r in range(3)], ["100.00", "50.00", "9.00"])

    def test_csv_file_is_loaded_with_correct_intensities(self):
        csv_path = self.base / 'scan.csv'
        csv_path.write_text("2theta,intensity\n" + "\n".join(f"{35 + 0.02 * i:.2f},{1000 + i}" for i in range(50)), encoding='utf-8')
        with patch.object(xrd.QFileDialog, 'getOpenFileNames', return_value=([str(csv_path)], '')):
            self.window.open_files()
        self.assertEqual(len(self.window.datasets), 1)
        self.assertAlmostEqual(float(self.window.datasets[0].y[0]), 1000.0)
        self.assertEqual(self.window.statusBar().currentMessage()[:7], "Wczytan")

    def test_export_writes_real_numbers_and_appends_extension_from_filter(self):
        try:
            import openpyxl
        except ImportError:
            self.skipTest("openpyxl")
        self.add_dataset()
        target = self.base / 'peaks'
        with patch.object(xrd.QFileDialog, 'getSaveFileName', return_value=(str(target), 'Excel (*.xlsx)')):
            self.window.export_table()
        sheet = openpyxl.load_workbook(self.base / 'peaks.xlsx')['XRD peaks']
        self.assertIsInstance(sheet.cell(2, 3).value, float)

    def test_export_failure_is_reported_instead_of_raised(self):
        self.add_dataset()
        missing = self.base / 'no-such-folder' / 'peaks.csv'
        with patch.object(xrd.QFileDialog, 'getSaveFileName', return_value=(str(missing), 'CSV (*.csv)')):
            self.window.export_table()
        self.assertTrue(self.messages['critical'].called)

    def test_wavelength_typing_does_not_invalidate_results(self):
        ds = self.add_dataset()
        ds.metadata['scientific_analysis'] = {'results': {'keep': True}}
        self.window.wavelength.lineEdit().setText("1")  # intermediate text while typing
        self.assertIn('scientific_analysis', ds.metadata)
        self.assertAlmostEqual(ds.wavelength_A, xrd.DEFAULT_WAVELENGTH_A)


if __name__ == '__main__':
    unittest.main()
