"""Silicon-standard calibration: model fitting, application to peaks, bookkeeping (no window)."""
import json
import math
import os
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, os.environ.get('XRD_TEST_APP_DIR', str(Path(__file__).resolve().parents[1])))
import XRD_Analyzer_GUI as app
import xrd_calibration as cal
import xrd_rd_import as rd
from xrd_profile_analysis import coupled_profile

LAM1, LAM2 = 1.540598, 1.544426
SI_FILE = Path(__file__).resolve().parent / 'data' / 'Si_polycrystal_PW3064.rd'


def synthetic_standard(zero=0.0, displacement=0.0, seed=1, step=0.0334, peak_height=900.0, a=cal.SI_LATTICE_A):
    """Si scan whose reflections are shifted by zero + displacement*cos(theta); Poisson noise."""
    rng = np.random.default_rng(seed)
    x = np.arange(24.0, 116.0, step)
    clean = np.full_like(x, 100.0)
    for hkl, theory in cal.si_reflections(LAM1, a, 25.0, 114.0):
        shifted = theory + zero + displacement * math.cos(math.radians(theory / 2))
        clean = clean + coupled_profile(x, shifted, peak_height * 3.0, 0.0, 0.10 + 0.0007 * (theory - 28), background_reference_deg=shifted)
    return x, rng.poisson(clean).astype(float)


def synthetic_sample_peaks(delta_fn, angles=(43.6, 50.8, 74.7, 90.6, 96.1)):
    """Peak list as the detector would produce it on the instrument with error delta_fn."""
    peaks = []
    for true in angles:
        observed = true + float(delta_fn(true))
        peaks.append(app.PeakResult(0, observed, LAM1 / (2 * math.sin(math.radians(observed / 2))), 1000.0, 100.0, 0.3, snr=50.0))
    return peaks


class SiliconGeometryTests(unittest.TestCase):
    def test_reflection_list(self):
        rows = dict(cal.si_reflections(LAM1, cal.SI_LATTICE_A, 20.0, 120.0))
        self.assertAlmostEqual(rows['111'], 28.442, places=2)
        self.assertAlmostEqual(rows['220'], 47.302, places=2)
        self.assertAlmostEqual(rows['400'], 69.129, places=2)
        for forbidden in ('200', '222', '420'):
            self.assertNotIn(forbidden, rows)
        for expected in ('311', '331', '422', '511', '440', '531'):
            self.assertIn(expected, rows)
        self.assertEqual(sorted(rows.values()), list(rows.values()))

    def test_allowed_rule(self):
        self.assertTrue(all(cal.si_reflection_allowed(*h) for h in ((1, 1, 1), (2, 2, 0), (4, 0, 0), (3, 3, 1), (5, 1, 1))))
        self.assertFalse(any(cal.si_reflection_allowed(*h) for h in ((1, 0, 0), (2, 0, 0), (1, 1, 0), (2, 2, 2), (4, 2, 0))))


class RealStandardTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.scan = rd.parse_philips_rd(SI_FILE)
        cls.calibration = cal.calibrate_from_standard(cls.scan.x, cls.scan.y, cls.scan.wavelength_A,
                                                      cls.scan.wavelength_kalpha2_A, source_name='Si PW3064')

    def test_most_reflections_are_usable(self):
        used = [r for r in self.calibration.reflections if r.used]
        self.assertGreaterEqual(len(used), 7)
        self.assertEqual(used[0].hkl, '111')

    def test_correction_is_small_smooth_and_explains_the_data(self):
        axis = self.calibration.axis
        self.assertLess(axis.rms_deg, 0.003)
        for r in self.calibration.reflections:
            if r.used:
                self.assertLess(abs(r.delta), 0.02)                         # the instrument is well aligned
                self.assertLess(abs(float(axis.corrected(r.two_theta_obs)) - r.two_theta_theory), 0.004)

    def test_instrument_profile_widens_with_angle(self):
        ipf = self.calibration.ipf
        self.assertIsNotNone(ipf)
        self.assertGreater(len(ipf['two_theta_deg']), 6)
        total = np.asarray(ipf['gaussian_fwhm_deg']) + np.asarray(ipf['lorentzian_fwhm_deg'])
        self.assertGreater(total[-1], total[0])
        self.assertTrue(0.05 < total[0] < 0.25)

    def test_notes_explain_the_assumptions(self):
        text = ' '.join(self.calibration.notes)
        self.assertIn('niecertyfikowany', text)
        self.assertIn('5.43102', text)


class AxisModelTests(unittest.TestCase):
    def reflections(self, zero, displacement, seed):
        x, y = synthetic_standard(zero, displacement, seed)
        return cal.fit_standard_reflections(x, y, LAM1, LAM2)

    def test_pure_zero_shift_is_recovered_and_the_simple_model_is_chosen(self):
        chosen = []
        for seed in range(6):
            rows = self.reflections(0.08, 0.0, seed)
            used = [r for r in rows if r.used]
            self.assertGreaterEqual(len(used), 8)
            axis = cal.fit_axis_correction([r.two_theta_obs for r in used], [r.two_theta_theory for r in used], [r.u_two_theta for r in used])
            self.assertAlmostEqual(float(axis.delta(60.0)), 0.08, delta=0.01)
            chosen.append(axis.model)
        self.assertGreaterEqual(chosen.count('zero'), 4, chosen)

    def test_displacement_is_recovered_and_detected(self):
        rows = self.reflections(0.03, -0.12, 11)
        used = [r for r in rows if r.used]
        axis = cal.fit_axis_correction([r.two_theta_obs for r in used], [r.two_theta_theory for r in used], [r.u_two_theta for r in used])
        self.assertIn('cos_theta', axis.terms)
        for tt in (30.0, 60.0, 100.0):
            truth = 0.03 - 0.12 * math.cos(math.radians(tt / 2))
            self.assertAlmostEqual(float(axis.delta(tt)), truth, delta=0.012)

    def test_a_sample_peak_list_is_corrected_back_to_the_truth(self):
        zero, disp = 0.05, -0.08
        rows = self.reflections(zero, disp, 12)
        calibration = cal.calibrate_from_reflections(rows, LAM1, LAM2)
        truth = lambda tt: zero + disp * math.cos(math.radians(tt / 2))
        peaks = synthetic_sample_peaks(truth)
        ds = app.Dataset(Path('s'), 'sample', '', np.arange(20.0, 120.0, 0.05), np.ones(2000), LAM1)
        ds.peaks = peaks
        cal.install_calibration([ds], ds, calibration)                 # ds is its own "standard": nothing applied
        ds.metadata.pop(cal.ROLE_KEY)
        record = cal.apply_calibration_to_fresh_peaks(ds)
        self.assertIsNotNone(record)
        for p, true in zip(ds.peaks, (43.6, 50.8, 74.7, 90.6, 96.1)):
            self.assertAlmostEqual(p.two_theta, true, delta=0.012)
            self.assertAlmostEqual(p.d_A, LAM1 / (2 * math.sin(math.radians(true / 2))), delta=2e-4)

    def test_model_selection_rules_and_errors(self):
        obs = np.array([28.4, 47.3, 56.1, 69.1]); theory = np.array([28.44, 47.30, 56.12, 69.13])
        axis = cal.fit_axis_correction(obs, theory, None, 'zero')
        self.assertEqual(axis.terms, ['zero'])
        with self.assertRaises(ValueError):
            cal.fit_axis_correction(obs[:2], theory[:2])
        with self.assertRaises(ValueError):
            cal.fit_axis_correction(obs, theory, None, 'zero+cos+sin2')      # needs 5 reflections
        with self.assertRaises(ValueError):
            cal.fit_axis_correction(obs, theory, None, 'nonsense')
        self.assertEqual(cal.fit_axis_correction(obs, theory).model in cal.MODELS, True)

    def test_too_few_reflections_in_the_scan(self):
        x, y = synthetic_standard(0.0, 0.0, 3)
        narrow = (x > 25) & (x < 50)
        with self.assertRaises(ValueError) as caught:
            cal.calibrate_from_standard(x[narrow], y[narrow], LAM1, LAM2)
        self.assertIn('co najmniej 3', str(caught.exception))

    def test_axis_correction_object(self):
        axis = cal.AxisCorrection(['zero', 'cos_theta'], [0.02, -0.05], [0.003, 0.01], [[1e-5, 0], [0, 1e-4]], 0.001, 1.0, 8, 28.0, 107.0, 'zero+cos')
        tt = np.array([30.0, 60.0, 100.0])
        np.testing.assert_allclose(axis.delta(tt), 0.02 - 0.05 * np.cos(np.radians(tt / 2)))
        np.testing.assert_allclose(axis.raw_from_corrected(axis.corrected(tt)), tt, atol=1e-9)
        self.assertTrue(np.all(axis.u_delta(tt) > 0))
        self.assertEqual(axis.covers(np.array([20.0, 60.0, 109.0, 125.0])).tolist(), [False, True, True, False])
        self.assertIsInstance(axis.delta(60.0), float)
        again = cal.AxisCorrection.from_dict(json.loads(json.dumps(axis.to_dict(), allow_nan=False)))
        np.testing.assert_allclose(again.params, axis.params)
        self.assertIn('Δ(2θ)', axis.describe())


class BookkeepingTests(unittest.TestCase):
    def setUp(self):
        x, y = synthetic_standard(0.06, 0.0, 21)
        self.calibration = cal.calibrate_from_standard(x, y, LAM1, LAM2, model='zero', source_name='Si test')
        self.axis_x = np.arange(20.0, 120.0, 0.05)

    def sample(self, name='probka', wavelength=LAM1):
        ds = app.Dataset(Path(name), name, '', self.axis_x, np.ones(len(self.axis_x)), wavelength)
        ds.peaks = synthetic_sample_peaks(lambda tt: 0.06)
        return ds

    def test_install_apply_undo_and_idempotence(self):
        std, a, b = self.sample('std'), self.sample('a'), self.sample('b', wavelength=1.78897)
        result = cal.install_calibration([std, a, b], std, self.calibration)
        self.assertEqual(std.metadata[cal.ROLE_KEY], 'standard')
        self.assertEqual(len(result['applied']), 1)
        self.assertIn('inna długość fali', result['skipped'][0])
        raw = [p.two_theta for p in a.peaks]
        cal.reapply_calibration(a)
        once = [p.two_theta for p in a.peaks]
        cal.reapply_calibration(a)
        twice = [p.two_theta for p in a.peaks]
        np.testing.assert_allclose(once, twice, atol=1e-12)
        self.assertTrue(all(abs((r - o) - 0.06) < 0.005 for r, o in zip(raw, once)))
        self.assertAlmostEqual(once[0], 43.6, delta=0.01)
        self.assertEqual([p.two_theta for p in std.peaks], raw)                # the standard itself is never shifted
        self.assertEqual([p.two_theta for p in b.peaks], raw)
        self.assertIn('skorygowane', cal.calibration_summary(a))
        self.assertIn('NIE została zastosowana', cal.calibration_summary(b))
        self.assertIn('Wzorzec kalibracyjny', cal.calibration_summary(std))
        removed = cal.remove_calibration([std, a, b])
        self.assertGreaterEqual(removed, 3)
        np.testing.assert_allclose([p.two_theta for p in a.peaks], raw, atol=1e-9)
        for ds in (std, a, b):
            self.assertFalse({cal.CAL_KEY, cal.ROLE_KEY, cal.APPLIED_KEY} & set(ds.metadata))
        self.assertEqual(cal.calibration_summary(a), '')

    def test_replacing_a_calibration_does_not_stack_corrections(self):
        std, a = self.sample('std'), self.sample('a')
        cal.install_calibration([std, a], std, self.calibration)
        cal.reapply_calibration(a)
        first = [p.two_theta for p in a.peaks]
        x, y = synthetic_standard(0.02, 0.0, 22)
        newer = cal.calibrate_from_standard(x, y, LAM1, LAM2, model='zero', source_name='Si newer')
        cal.install_calibration([std, a], std, newer)
        cal.reapply_calibration(a)
        second = [p.two_theta for p in a.peaks]
        raw = [p.two_theta + 0.0 for p in self.sample('fresh').peaks]
        self.assertTrue(all(abs((r - s) - 0.02) < 0.006 for r, s in zip(raw, second)), (raw, second))
        self.assertNotEqual(first, second)

    def test_fresh_analysis_applies_once_and_reports_skipped_calibration(self):
        std, a, b = self.sample('std'), self.sample('a'), self.sample('b', wavelength=1.78897)
        cal.install_calibration([std, a, b], std, self.calibration)
        for ds in (a, b, std):
            raw = [p.two_theta for p in ds.peaks]
            cal.apply_calibration_to_fresh_peaks(ds)
            cal.apply_calibration_to_fresh_peaks(self.refresh(ds, raw))     # a second analysis starts again from raw peaks
        self.assertAlmostEqual(a.peaks[0].two_theta, 43.6, delta=0.01)
        self.assertEqual(b.peaks[0].two_theta, std.peaks[0].two_theta)
        self.assertIn('inna długość fali', b.metadata['calibration_skipped'])

    @staticmethod
    def refresh(ds, raw):
        for p, value in zip(ds.peaks, raw):
            p.two_theta = value
            p.d_A = ds.wavelength_A / (2 * math.sin(math.radians(value / 2)))
        return ds

    def test_calibrated_axis(self):
        std, a = self.sample('std'), self.sample('a')
        cal.install_calibration([std, a], std, self.calibration)
        axis, used = cal.calibrated_axis(a)
        self.assertIsNotNone(used)
        self.assertTrue(np.all(np.diff(axis) > 0))
        self.assertAlmostEqual(float(np.mean(a.x - axis)), 0.06, delta=0.005)
        raw_axis, none = cal.calibrated_axis(std)
        self.assertIsNone(none)
        np.testing.assert_array_equal(raw_axis, std.x)

    def test_damaged_calibration_data_is_reported_not_raised(self):
        ds = self.sample()
        ds.metadata[cal.CAL_KEY] = {'axis': {'terms': ['zero']}}
        calibration, reason = cal.calibration_status(ds)
        self.assertIsNone(calibration)
        self.assertIn('uszkodzone', reason)
        self.assertIsNone(cal.apply_calibration_to_fresh_peaks(ds))

    def test_failed_reflection_fits_survive_json_and_the_project_format(self):
        failed = cal.StandardReflection('511', 94.95, math.nan, math.nan, math.nan, math.nan, math.nan, None, None, None, False, 'nieudane')
        text = json.dumps(failed.to_dict(), allow_nan=False)
        again = cal.StandardReflection.from_dict(json.loads(text))
        self.assertTrue(math.isnan(again.two_theta_obs))
        self.assertFalse(again.used)
        calibration = self.calibration
        calibration.reflections.append(failed)
        json.dumps(calibration.to_dict(), allow_nan=False)
        self.assertEqual(len(cal.StandardCalibration.from_dict(json.loads(json.dumps(calibration.to_dict()))).reflections),
                         len(calibration.reflections))

    def test_instrument_profile_needs_two_reflections(self):
        rows = [cal.StandardReflection('111', 28.44, 28.43, 0.001, 0.1, 0.0, 0.1, None, None, 0.99, True)]
        self.assertIsNone(cal.build_instrument_profile(rows, 'x', LAM1))
        rows.append(cal.StandardReflection('220', 47.3, 47.29, 0.001, 0.12, 0.0, 0.12, None, None, 0.99, True))
        model = cal.build_instrument_profile(rows, 'x', LAM1)
        self.assertEqual(len(model.two_theta_deg), 2)


class ProgramIntegrationTests(unittest.TestCase):
    """The calibration as seen through the program's own analysis and wavelength code."""

    @staticmethod
    def window():
        from types import SimpleNamespace as NS
        value = lambda v: NS(value=lambda: v)
        return NS(baseline_cb=NS(isChecked=lambda: False), prominence=value(3.0), min_snr=value(8.0), min_rel_height=value(5.0),
                  min_fwhm=value(0.1), min_distance=value(0.25), max_peaks=value(50), phase_tol=value(0.35))

    @staticmethod
    def pattern(shift=0.0, seed=0):
        rng = np.random.default_rng(seed)
        x = np.arange(30.0, 100.0, 0.05)
        clean = np.full_like(x, 300.0)
        for centre in (43.6, 50.8, 74.7, 90.6):
            top = coupled_profile(x, centre + shift, 1.0, 0.12, 0.08, background_reference_deg=centre + shift)
            clean = clean + 4000 * top / top.max()
        return x, rng.poisson(clean).astype(float)

    def calibration(self, model='zero'):
        x, y = synthetic_standard(0.06, 0.0, 31)
        return cal.calibrate_from_standard(x, y, LAM1, LAM2, model=model, source_name='Si')

    def test_analysis_corrects_samples_but_not_the_standard(self):
        x, y = self.pattern(shift=0.06)
        sample = app.Dataset(Path('s'), 'sample', '', x, y, LAM1)
        standard = app.Dataset(Path('si'), 'standard', '', x, y.copy(), LAM1)
        cal.install_calibration([sample, standard], standard, self.calibration())
        uncalibrated = app.Dataset(Path('u'), 'uncalibrated', '', x, y.copy(), LAM1)
        for ds in (sample, standard, uncalibrated):
            app.MainWindow._analyze_dataset_peaks_only(self.window(), ds)
        # The peak list holds the maximum of the Kalpha1+Kalpha2 envelope (a few hundredths of a degree above
        # Kalpha1); the calibration removes the instrument error (0.06 deg here) and nothing else.
        raw = [p.two_theta for p in uncalibrated.peaks]
        self.assertGreaterEqual(len(raw), 4)
        self.assertEqual(len(raw), len(sample.peaks))
        for corrected, before in zip((p.two_theta for p in sample.peaks), raw):
            self.assertAlmostEqual(before - corrected, 0.06, delta=0.006)
        for standard_position, before in zip((p.two_theta for p in standard.peaks), raw):
            self.assertEqual(standard_position, before)                # the standard itself is never shifted
        for p in sample.peaks:
            self.assertAlmostEqual(p.d_A, LAM1 / (2 * math.sin(math.radians(p.two_theta / 2))), places=9)
        self.assertEqual(sample.metadata[cal.APPLIED_KEY]['n_peaks'], len(sample.peaks))
        self.assertNotIn(cal.APPLIED_KEY, standard.metadata)

    def test_reanalysis_does_not_apply_the_correction_twice(self):
        x, y = self.pattern(shift=0.06, seed=2)
        sample = app.Dataset(Path('s'), 'sample', '', x, y, LAM1)
        standard = app.Dataset(Path('si'), 'standard', '', x, y.copy(), LAM1)
        cal.install_calibration([sample, standard], standard, self.calibration())
        results = []
        for _ in range(3):
            app.MainWindow._analyze_dataset_peaks_only(self.window(), sample)
            results.append([p.two_theta for p in sample.peaks])
        np.testing.assert_allclose(results[0], results[1], atol=1e-12)
        np.testing.assert_allclose(results[1], results[2], atol=1e-12)

    def test_changing_the_wavelength_withdraws_the_correction(self):
        x, y = self.pattern(shift=0.06, seed=3)
        sample = app.Dataset(Path('s'), 'sample', '', x, y, LAM1)
        standard = app.Dataset(Path('si'), 'standard', '', x, y.copy(), LAM1)
        cal.install_calibration([sample, standard], standard, self.calibration())
        app.MainWindow._analyze_dataset_peaks_only(self.window(), sample)
        corrected = [p.two_theta for p in sample.peaks]
        app.set_dataset_wavelength(sample, 1.78897, 'user')
        raw = [p.two_theta for p in sample.peaks]
        self.assertTrue(all(abs((r - c) - 0.06) < 0.01 for r, c in zip(raw, corrected)))     # raw positions are back
        self.assertNotIn(cal.APPLIED_KEY, sample.metadata)
        self.assertIn('inna długość fali', sample.metadata['calibration_skipped'])
        for p in sample.peaks:
            self.assertAlmostEqual(p.d_A, 1.78897 / (2 * math.sin(math.radians(p.two_theta / 2))), places=9)
        app.set_dataset_wavelength(sample, LAM1, 'user')                                   # and it comes back
        self.assertIn(cal.APPLIED_KEY, sample.metadata)

    def test_project_save_format_round_trip(self):
        x, y = self.pattern(shift=0.06, seed=4)
        sample = app.Dataset(Path('s'), 'sample', '', x, y, LAM1)
        standard = app.Dataset(Path('si'), 'standard', '', x, y.copy(), LAM1)
        cal.install_calibration([sample, standard], standard, self.calibration())
        app.MainWindow._analyze_dataset_peaks_only(self.window(), sample)
        stored = json.loads(json.dumps(sample.metadata, allow_nan=False, default=str))
        restored = app.Dataset(Path('s'), 'sample', '', x, y, LAM1)
        restored.metadata.update(stored)
        calibration, reason = cal.calibration_status(restored)
        self.assertEqual(reason, '')
        np.testing.assert_allclose(calibration.axis.params, self.calibration().axis.params, atol=0.01)
        self.assertIn('skorygowane', cal.calibration_summary(restored))


if __name__ == '__main__':
    unittest.main()
