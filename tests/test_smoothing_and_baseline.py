"""Peak-preserving display smoothing and step-independent AsLS baseline (numerical, no window)."""
import os
import sys
import unittest
from pathlib import Path

import numpy as np
from scipy.signal import savgol_filter

sys.path.insert(0, os.environ.get('XRD_TEST_APP_DIR', str(Path(__file__).resolve().parents[1])))
import XRD_Analyzer_GUI as app
from xrd_plot_support import (DEFAULT_SMOOTH_METHOD, SMOOTH_METHODS, display_smoothing_settings,
                              peak_preserving_filter, smooth_for_display)
from xrd_profile_analysis import coupled_profile

CENTRES = (44.0, 51.5, 75.2, 91.5)


def pattern(step, g_fwhm, seed, centres=CENTRES, noisy=True):
    """Background + hump + reflections; returns x, noisy y, clean y."""
    rng = np.random.default_rng(seed)
    x = np.arange(30.0, 100.0, step)
    clean = 30 + 0.3 * (x - 30) + 120 * np.exp(-0.5 * ((x - 45) / 8.0) ** 2)
    for centre in centres:
        clean = clean + coupled_profile(x, centre, 2000.0, g_fwhm, 0.6 * g_fwhm, background_reference_deg=centre)
    return x, (rng.poisson(clean).astype(float) if noisy else clean.copy()), clean


def kept(smoothed, clean, x, centres=CENTRES):
    return float(np.mean([smoothed[abs(x - c) < 0.6].max() / clean[abs(x - c) < 0.6].max() for c in centres]))


def background_noise(smoothed, clean, x, centres=CENTRES):
    far = np.ones_like(x, dtype=bool)
    for c in centres:
        far &= abs(x - c) > 3
    return float(np.std((smoothed - clean)[far]))


class PeakPreservingSmoothingTests(unittest.TestCase):
    def test_default_method_is_still_the_plain_filter(self):
        x, y, _ = pattern(0.02, 0.12, 1)
        np.testing.assert_array_equal(smooth_for_display(x, y, 100), smooth_for_display(x, y, 100, method="savgol"))
        np.testing.assert_allclose(smooth_for_display(x, y, 100), savgol_filter(y, 13, 2, mode="interp"))

    def test_plain_filter_flattens_narrow_peaks_and_new_method_does_not(self):
        for step in (0.01, 0.02, 0.05):
            x, y, clean = pattern(step, 0.06, 2)
            plain = smooth_for_display(x, y, 100)
            better = smooth_for_display(x, y, 100, method="peak_preserving")
            with self.subTest(step=step):
                self.assertLess(kept(plain, clean, x), 0.85)       # documents the defect being fixed
                self.assertGreater(kept(better, clean, x), 0.97)
                self.assertLess(np.sqrt(np.mean((better - clean) ** 2)), 0.2 * np.sqrt(np.mean((plain - clean) ** 2)))

    def test_background_noise_reduction_is_unchanged(self):
        for step in (0.01, 0.02):
            x, y, clean = pattern(step, 0.12, 3)
            ratio = (background_noise(smooth_for_display(x, y, 100, method="peak_preserving"), clean, x)
                     / background_noise(smooth_for_display(x, y, 100), clean, x))
            with self.subTest(step=step):
                self.assertTrue(0.9 < ratio < 1.1, ratio)
                self.assertLess(background_noise(smooth_for_display(x, y, 100, method="peak_preserving"), clean, x),
                                0.6 * background_noise(y, clean, x))

    def test_broad_peaks_are_treated_like_the_plain_filter(self):
        x, y, clean = pattern(0.02, 0.6, 4)
        plain = smooth_for_display(x, y, 100)
        better = smooth_for_display(x, y, 100, method="peak_preserving")
        self.assertLess(np.sqrt(np.mean((better - plain) ** 2)), 0.1 * np.std(y - clean))

    def test_narrow_and_broad_peaks_in_one_pattern(self):
        x = np.arange(10.0, 80.0, 0.02)
        rng = np.random.default_rng(5)
        clean = (30 + coupled_profile(x, 30.0, 3000.0, 0.05, 0.03, background_reference_deg=30.0)
                 + coupled_profile(x, 60.0, 3000.0, 0.8, 0.5, background_reference_deg=60.0))
        y = rng.poisson(clean).astype(float)
        plain = smooth_for_display(x, y, 100)
        better = smooth_for_display(x, y, 100, method="peak_preserving")
        narrow, broad = abs(x - 30) < 0.5, abs(x - 60) < 2
        self.assertGreater(better[narrow].max() / clean[narrow].max(), 0.9)
        self.assertLess(plain[narrow].max() / clean[narrow].max(), 0.8)
        # the broad reflection keeps the full noise reduction of the plain filter
        self.assertLess(np.sqrt(np.mean((better - clean)[broad] ** 2)), 1.1 * np.sqrt(np.mean((plain - clean)[broad] ** 2)))

    def test_peak_position_is_preserved(self):
        x, y, clean = pattern(0.01, 0.06, 6, centres=(60.0,))
        smoothed = smooth_for_display(x, y, 100, method="peak_preserving")
        window = abs(x - 60.0) < 0.5
        self.assertLessEqual(abs(x[window][np.argmax(smoothed[window])] - x[window][np.argmax(clean[window])]), 0.02)

    def test_noise_free_data_are_not_distorted(self):
        x, _, clean = pattern(0.02, 0.06, 7, noisy=False)
        np.testing.assert_array_equal(smooth_for_display(x, clean, 100, method="peak_preserving"), clean)

    def test_strength_is_an_exact_mixture_and_zero_is_identity(self):
        x, y, _ = pattern(0.02, 0.06, 8)
        full = smooth_for_display(x, y, 100, method="peak_preserving")
        np.testing.assert_allclose(smooth_for_display(x, y, 30, method="peak_preserving"), 0.7 * y + 0.3 * full)
        np.testing.assert_array_equal(smooth_for_display(x, y, 0, method="peak_preserving"), y)

    def test_result_stays_in_data_range_and_inputs_are_untouched(self):
        x, y, _ = pattern(0.02, 0.06, 9)
        x0, y0 = x.copy(), y.copy()
        out = smooth_for_display(x, y, 100, method="peak_preserving")
        np.testing.assert_array_equal(x, x0); np.testing.assert_array_equal(y, y0)
        self.assertGreaterEqual(out.min(), y.min()); self.assertLessEqual(out.max(), y.max())
        np.testing.assert_array_equal(out, smooth_for_display(x, y, 100, method="peak_preserving"))  # deterministic

    def test_irregular_axis_and_tiny_arrays(self):
        x, y, _ = pattern(0.02, 0.06, 10)
        jitter = np.sort(x + np.random.default_rng(1).uniform(-0.004, 0.004, len(x)))
        out = smooth_for_display(jitter, y, 100, method="peak_preserving")
        self.assertEqual(out.shape, y.shape); self.assertTrue(np.all(np.isfinite(out)))
        for n in (0, 1, 4, 5, 6, 9, 14, 40):
            xs = np.arange(n) * 0.02
            ys = np.random.default_rng(n).poisson(50, n).astype(float)
            result = smooth_for_display(xs, ys, 100, method="peak_preserving")
            self.assertEqual(result.shape, ys.shape); self.assertTrue(np.all(np.isfinite(result)))
        self.assertIsNone(peak_preserving_filter(np.ones(4), 0.02))

    def test_unusual_data_do_not_break_the_filter(self):
        x = np.arange(30.0, 60.0, 0.02)
        for values in (np.full(len(x), 5.0), np.zeros(len(x)), np.random.default_rng(2).poisson(0.3, len(x)).astype(float)):
            out = smooth_for_display(x, values, 100, method="peak_preserving")
            self.assertEqual(out.shape, values.shape); self.assertTrue(np.all(np.isfinite(out)))
        spike = np.random.default_rng(3).poisson(5, len(x)).astype(float); spike[700] = 4000.0
        out = smooth_for_display(x, spike, 100, method="peak_preserving")
        self.assertTrue(np.all(np.isfinite(out))); self.assertLessEqual(out.max(), 4000.0)

    def test_display_settings_methods(self):
        self.assertEqual(DEFAULT_SMOOTH_METHOD, "peak_preserving")
        self.assertEqual(display_smoothing_settings({})["smooth_method"], "peak_preserving")
        for method in SMOOTH_METHODS:
            self.assertEqual(display_smoothing_settings({"smooth_method": method, "smooth_enabled": True})["smooth_method"], method)
        stored = display_smoothing_settings({"smooth_method": "savgol", "smooth_enabled": True})
        self.assertTrue(stored["smooth_enabled"])                 # figures saved before 20.3.2 keep their look
        invalid = display_smoothing_settings({"smooth_method": "bogus", "smooth_enabled": True})
        self.assertFalse(invalid["smooth_enabled"]); self.assertEqual(invalid["smooth_method"], "peak_preserving")
        nested = display_smoothing_settings({"display_smoothing": {"enabled": True, "method": "savgol", "strength_percent": 40}})
        self.assertEqual((nested["smooth_enabled"], nested["smooth_method"], nested["smooth_strength_percent"]), (True, "savgol", 40.0))


class BaselineStepIndependenceTests(unittest.TestCase):
    @staticmethod
    def scan(step, seed=3):
        rng = np.random.default_rng(seed)
        x = np.arange(20.0, 100.0, step)
        background = 60 + 0.5 * (x - 30) + 300 * np.exp(-0.5 * ((x - 45) / 5.0) ** 2)
        sharp = coupled_profile(x, 60.0, 2000.0, 0.15, 0.09, background_reference_deg=60.0)
        broad = coupled_profile(x, 80.0, 2000.0, 1.5, 0.9, background_reference_deg=80.0)
        return x, rng.poisson(background + sharp + broad).astype(float), background, sharp, broad

    def test_without_axis_the_legacy_result_is_unchanged(self):
        x, y, *_ = self.scan(0.02)
        np.testing.assert_array_equal(app.baseline_asls(y), app.baseline_asls(y, x=None))

    def test_reference_step_gives_the_same_baseline_as_before(self):
        # 0.05 deg is the step of the real 316L / CuNi / Monel scans for which lam=2e6 was tuned
        # (a 0.02 deg reference lifted the baseline under their peaks and removed a reflection).
        x, y, *_ = self.scan(0.05)
        np.testing.assert_allclose(app.baseline_asls(y), app.baseline_asls(y, x=x), rtol=1e-6)

    def test_fine_steps_no_longer_absorb_broad_reflections(self):
        for step in (0.005, 0.01):
            x, y, background, sharp, broad = self.scan(step)
            near_broad = abs(x - 80) < 0.6
            legacy = (y - app.baseline_asls(y))[near_broad].max() / broad[near_broad].max()
            scaled = (y - app.baseline_asls(y, x=x))[near_broad].max() / broad[near_broad].max()
            with self.subTest(step=step):
                self.assertGreater(scaled, 0.8)
                self.assertGreater(scaled - legacy, 0.2)

    def test_baseline_error_is_similar_for_every_step(self):
        errors = []
        for step in (0.005, 0.01, 0.02, 0.05, 0.1):
            x, y, background, *_ = self.scan(step)
            errors.append(float(np.sqrt(np.mean((app.baseline_asls(y, x=x) - background) ** 2))))
        self.assertLess(max(errors) / min(errors), 1.25, errors)

    def test_sharp_reflections_keep_their_height_at_every_step(self):
        for step in (0.005, 0.02, 0.1):
            x, y, background, sharp, broad = self.scan(step)
            near = abs(x - 60) < 0.3
            self.assertGreater((y - app.baseline_asls(y, x=x))[near].max() / sharp[near].max(), 0.95)

    def test_unusable_axis_falls_back_to_the_legacy_stiffness(self):
        x, y, *_ = self.scan(0.02)
        for bad in (x[:-1], np.zeros(len(y)), x[::-1]):
            np.testing.assert_array_equal(app.baseline_asls(y, x=bad), app.baseline_asls(y))


if __name__ == '__main__':
    unittest.main()
