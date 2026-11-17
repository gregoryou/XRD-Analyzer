"""Regressions for wavelength-dependent results and angular FWHM."""
import math
import unittest
from copy import deepcopy
from types import SimpleNamespace
from unittest.mock import patch

from test_numeric_smoke import app, np, Path, Value, detector_settings


DEPENDENT_KEYS = (
    'phase_fits', 'structural_family_fits', 'manual_cod_peak_assignments',
    'ai_phase_ranking', 'ai_user_confirmed_candidate', 'publication_phase_aliases',
)


def analyzed_dataset():
    peak = app.PeakResult(1, 60.0, app.bragg_d(60.0, 1.5406), 100.0, 100.0, 0.2,
                          phase='old phase', hkl='111', delta_2theta=0.01,
                          custom_label='keep this label', phase_source='cod')
    dataset = app.Dataset(Path('synthetic.dat'), 'synthetic', '',
                          np.array([59., 60., 61.]), np.array([0., 100., 0.]), peaks=[peak])
    for key in DEPENDENT_KEYS:
        dataset.metadata[key] = [{'phase': 'old phase'}]
        dataset.metadata_sources[key] = 'test source'
    dataset.metadata['allowed_elements'] = 'Fe Ni'
    return dataset


def profile_dialog(dataset, profile):
    return SimpleNamespace(dataset=dataset, all_datasets=[dataset],
        _selected_profile=lambda: profile, fill_missing_cb=Value(False),
        _load_dataset_fields=lambda: None, parent=lambda: None)


class WavelengthRegressionTests(unittest.TestCase):
    def assert_updated(self, dataset):
        self.assertEqual(dataset.wavelength_A, 2.2897)
        self.assertAlmostEqual(dataset.peaks[0].d_A, 2.2897, places=12)
        self.assertEqual(dataset.peaks[0].two_theta, 60.0)
        self.assertEqual(dataset.peaks[0].height, 100.0)
        self.assertEqual(dataset.peaks[0].fwhm_deg, 0.2)
        self.assertEqual(dataset.peaks[0].custom_label, 'keep this label')
        self.assertEqual(dataset.peaks[0].phase, 'unassigned')
        self.assertEqual(dataset.peaks[0].hkl, '')
        self.assertEqual(dataset.peaks[0].phase_source, '')
        self.assertTrue(math.isnan(dataset.peaks[0].delta_2theta))
        for key in DEPENDENT_KEYS:
            self.assertFalse(dataset.metadata.get(key), key)
            self.assertNotIn(key, dataset.metadata_sources)

    def test_profile_recalculates_and_invalidates_dependent_results(self):
        dataset = analyzed_dataset()
        app.apply_instrument_profile(dataset, {'name': 'test', 'wavelength_A': 2.2897}, False)
        self.assert_updated(dataset)
        self.assertEqual(dataset.metadata['allowed_elements'], 'Fe Ni')

    def test_same_wavelength_preserves_results(self):
        dataset = analyzed_dataset()
        before = deepcopy(dataset.peaks[0].__dict__)
        dependent = {key: deepcopy(dataset.metadata[key]) for key in DEPENDENT_KEYS}
        app.apply_instrument_profile(dataset, {'name': 'test', 'wavelength_A': 1.5406}, False)
        self.assertEqual(dataset.peaks[0].__dict__, before)
        self.assertEqual({key: dataset.metadata[key] for key in DEPENDENT_KEYS}, dependent)

    def test_fill_missing_preserves_verified_wavelength(self):
        dataset = analyzed_dataset()
        dataset.metadata_sources['wavelength_A'] = 'file'
        app.apply_instrument_profile(dataset, {'name': 'test', 'wavelength_A': 2.2897}, True)
        self.assertEqual(dataset.wavelength_A, 1.5406)
        self.assertEqual(dataset.peaks[0].phase, 'old phase')

    def test_manual_metadata_path(self):
        dataset = analyzed_dataset()
        dialog = SimpleNamespace(dataset=dataset,
            FIELD_SPECS=app.InstrumentMeasurementDialog.FIELD_SPECS,
            _widget_value=lambda key: 2.2897 if key == 'wavelength_A' else None,
            _load_dataset_fields=lambda: None, parent=lambda: None)
        with patch.object(app.QMessageBox, 'information') as message:
            app.InstrumentMeasurementDialog.save_current_metadata(dialog)
        self.assert_updated(dataset)
        self.assertIn('ponown', message.call_args.args[2].lower())

    def test_main_window_metadata_path(self):
        dataset = analyzed_dataset()
        window = SimpleNamespace(current_dataset=lambda: dataset,
            sample_name=SimpleNamespace(text=lambda: 'synthetic'),
            material_type=SimpleNamespace(text=lambda: ''),
            allowed_elements_edit=SimpleNamespace(text=lambda: 'Fe Ni'),
            wavelength=Value(2.2897), sample_combo=SimpleNamespace(currentIndex=lambda: -1),
            update_table=lambda: None, redraw=lambda: None, update_summary=lambda: None,
            mark_project_dirty=lambda: None)
        with patch.object(app.QMessageBox, 'information') as message:
            app.MainWindow.update_metadata(window)
        self.assert_updated(dataset)
        message.assert_called_once()

    def test_profile_dialog_notifies_only_when_results_invalidated(self):
        dataset = analyzed_dataset()
        dialog = profile_dialog(dataset, {'name': 'test', 'wavelength_A': 2.2897})
        with patch.object(app.QMessageBox, 'information') as message:
            app.InstrumentMeasurementDialog.apply_profile_current(dialog)
            message.assert_called_once()
        self.assert_updated(dataset)
        with patch.object(app.QMessageBox, 'information') as message:
            app.InstrumentMeasurementDialog.apply_profile_current(dialog)
            message.assert_not_called()

    def test_batch_profile_path(self):
        dataset = analyzed_dataset()
        dialog = profile_dialog(dataset, {'name': 'test', 'wavelength_A': 2.2897})
        with patch.object(app.QMessageBox, 'information') as message:
            app.InstrumentMeasurementDialog.apply_profile_all(dialog)
        self.assert_updated(dataset)
        self.assertIn('ponown', message.call_args.args[2].lower())


class AngularWidthRegressionTests(unittest.TestCase):
    def test_nonuniform_gaussian_width(self):
        x = np.concatenate((np.arange(10., 20., .01), np.arange(20., 30.0001, .1)))
        y = 100.0 * np.exp(-0.5 * ((x - 25.0) / 0.2) ** 2)
        dataset = app.Dataset(Path('synthetic.dat'), 'synthetic', '', x, y)
        app.MainWindow._analyze_dataset_peaks_only(detector_settings(), dataset)
        self.assertEqual(len(dataset.peaks), 1)
        self.assertAlmostEqual(dataset.peaks[0].fwhm_deg,
                               2.0 * math.sqrt(2.0 * math.log(2.0)) * 0.2, delta=0.006)

    def test_width_filter_uses_degrees_on_nonuniform_grid(self):
        x = np.concatenate((np.arange(10., 20., .01), np.arange(20., 30.0001, .1)))
        y = 100.0 * np.exp(-0.5 * ((x - 25.0) / 0.2) ** 2)
        dataset = app.Dataset(Path('synthetic.dat'), 'synthetic', '', x, y)
        settings = detector_settings()
        settings.min_fwhm = Value(0.1)
        app.MainWindow._analyze_dataset_peaks_only(settings, dataset)
        self.assertEqual(len(dataset.peaks), 1)

    def test_invalid_axes_rejected_before_changing_results(self):
        for x in ([3., 2., 1.], [1., 1., 2.], [1., math.nan, 2.]):
            with self.subTest(x=x):
                dataset = analyzed_dataset()
                dataset.x = np.array(x)
                with self.assertRaises(ValueError):
                    app.MainWindow._analyze_dataset_peaks_only(detector_settings(), dataset)
                self.assertEqual(dataset.peaks[0].phase, 'old phase')


if __name__ == '__main__':
    unittest.main()
