"""Width-aware peak acceptance (no window needed)."""
import math
import os
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace as NS

import numpy as np
from scipy.signal import find_peaks

sys.path.insert(0, os.environ.get('XRD_TEST_APP_DIR', str(Path(__file__).resolve().parents[1])))
import XRD_Analyzer_GUI as app
from xrd_detection_support import find_candidate_peaks, matched_significance, noise_sigma_from_differences
from xrd_profile_analysis import coupled_profile


def old_rule(sig, sigma, scale, pct=3.0, snr=8.0, rel_h=5.0):
    idx, props = find_peaks(sig, prominence=max(pct / 100 * scale, snr * sigma), height=rel_h / 100 * scale)
    return idx


def new_rule(sig, snr=8.0, pct=3.0, rel_h=5.0):
    sigma = noise_sigma_from_differences(sig)
    return find_candidate_peaks(sig, sigma, float(sig.max()), pct, snr, rel_h)[0], sigma


def _shape(x, centre, top):
    profile = coupled_profile(x, centre, 1.0, 0.30, 0.18, background_reference_deg=centre)
    return top * profile / profile.max()


def _flank_clean(weak_top):
    """Main reflection 14000 counts high, weak broad one weak_top counts high, 0.8-1.05 deg away."""
    x = np.arange(35.0, 55.0, 0.05)
    return x, _shape(x, 43.7, 14000.0) + _shape(x, 44.75, weak_top)


def _weak_height_for_prominence(target):
    """Amplitude of the weak reflection that gives the requested noise-free prominence (bisection)."""
    low, high = 10.0, 30000.0
    for _ in range(40):
        mid = 0.5 * (low + high)
        x, clean = _flank_clean(mid)
        idx, props = find_peaks(clean, prominence=0.0)
        near = [i for i in idx if abs(x[i] - 44.75) < 0.3]
        prominence = max((props['prominences'][list(idx).index(i)] for i in near), default=0.0)
        low, high = (mid, high) if prominence < target else (low, mid)
    return 0.5 * (low + high)


def flank_pattern(seed, prominence=760.0):
    """Strong reflection with a weak, broad one on its flank; background 12000, Poisson noise (sigma ~ 118)."""
    x, clean = _flank_clean(_weak_height_for_prominence(prominence))
    rng = np.random.default_rng(seed)
    return x, rng.poisson(12000 + clean).astype(float) - 12000.0, clean


class CandidateSelectionTests(unittest.TestCase):
    def test_weak_broad_peak_next_to_a_strong_one_is_found_where_the_old_rule_misses_it(self):
        found_new = found_old = 0
        trials = 40
        for seed in range(trials):
            x, sig, clean = flank_pattern(seed)
            # the weak reflection must be a real local maximum for the test to be meaningful
            self.assertTrue(clean[np.argmin(abs(x - 44.75))] > clean[np.argmin(abs(x - 44.45))])
            self.assertTrue(clean[np.argmin(abs(x - 44.75))] > clean[np.argmin(abs(x - 45.05))])
            idx_new, sigma = new_rule(sig)
            idx_old = old_rule(sig, sigma, float(sig.max()))
            found_new += bool(np.any(abs(x[idx_new] - 44.75) < 0.2))
            found_old += bool(np.any(abs(x[idx_old] - 44.75) < 0.2))
        # 6.4 sigma prominence, the case measured on the real 316L base-material scan (742 counts, limit 931)
        self.assertEqual(found_new, trials)
        self.assertLessEqual(found_old, 0.5 * trials)

    def test_nothing_the_old_rule_detects_is_lost(self):
        rng = np.random.default_rng(1)
        for seed in range(25):
            rng = np.random.default_rng(100 + seed)
            x = np.arange(20.0, 120.0, 0.05)
            clean = 3000 + 0.5 * x * 10
            for centre, height in zip(rng.uniform(25, 115, 6), rng.uniform(300, 30000, 6)):
                clean = clean + coupled_profile(x, centre, float(height), 0.25, 0.15, background_reference_deg=centre)
            sig = rng.poisson(clean).astype(float) - 3000 - 5 * x
            idx_new, sigma = new_rule(sig)
            old = set(old_rule(sig, sigma, float(sig.max())).tolist())
            self.assertTrue(old <= set(idx_new.tolist()), seed)

    def test_pure_noise_gives_almost_no_peaks(self):
        rng = np.random.default_rng(0)
        false_peaks = 0
        for level in (200, 2000, 12000):
            for _ in range(60):
                sig = rng.poisson(level, 2000).astype(float) - level
                false_peaks += len(new_rule(sig)[0])
        self.assertLessEqual(false_peaks, 3)

    def test_periodic_ripple_and_broad_hump_are_not_peaks(self):
        rng = np.random.default_rng(2)
        x = np.arange(20.0, 60.0, 0.05)
        ripple = 250 * np.cos(2 * np.pi * np.arange(len(x)) / 5.0)
        hump = 800 * np.exp(-0.5 * ((x - 30) / 3.0) ** 2)
        strong = coupled_profile(x, 50.0, 20000.0, 0.3, 0.2, background_reference_deg=50.0)
        sig = rng.poisson(10000 + strong).astype(float) - 10000 + ripple + hump
        idx, _ = new_rule(sig)
        positions = x[idx]
        self.assertTrue(np.any(abs(positions - 50.0) < 0.3))
        self.assertFalse(np.any(positions < 45.0), positions)

    def test_low_snr_setting_reduces_to_the_original_rule(self):
        x, sig, _ = flank_pattern(3)
        sigma = noise_sigma_from_differences(sig)
        for snr in (3.0, 4.0, 5.0):
            idx, _ = find_candidate_peaks(sig, sigma, float(sig.max()), 3.0, snr, 5.0)
            np.testing.assert_array_equal(idx, old_rule(sig, sigma, float(sig.max()), snr=snr))

    def test_matched_significance_separates_signal_from_noise(self):
        rng = np.random.default_rng(3)
        noise = rng.normal(0, 100, 3000)
        peaks, _ = find_peaks(noise, prominence=300)
        self.assertLess(matched_significance(noise, 100.0, peaks).max(), 9.0)
        bump = noise + 900 * np.exp(-0.5 * ((np.arange(3000) - 1500) / 4.0) ** 2)
        self.assertGreater(matched_significance(bump, 100.0, [1500])[0], 15.0)
        self.assertEqual(len(matched_significance(bump, 100.0, [])), 0)

    def test_noise_estimate_ignores_the_peaks(self):
        rng = np.random.default_rng(4)
        x = np.arange(20.0, 100.0, 0.05)
        sig = rng.normal(0, 50, len(x)) + coupled_profile(x, 60.0, 50000.0, 0.3, 0.2, background_reference_deg=60.0)
        self.assertAlmostEqual(noise_sigma_from_differences(sig), 50.0, delta=4.0)

    def test_analysis_method_uses_the_new_rule(self):
        x, sig, _ = flank_pattern(5)
        ds = app.Dataset(Path('x'), 'flank', '', x, sig + 600.0, 1.540598)   # baseline switched off below
        V = lambda v: NS(value=lambda: v)
        window = NS(baseline_cb=NS(isChecked=lambda: False), prominence=V(3.0), min_snr=V(8.0), min_rel_height=V(5.0),
                    min_fwhm=V(0.1), min_distance=V(0.25), max_peaks=V(50), phase_tol=V(0.35))
        app.MainWindow._analyze_dataset_peaks_only(window, ds)
        positions = [p.two_theta for p in ds.peaks]
        self.assertTrue(any(abs(p - 43.7) < 0.2 for p in positions))
        self.assertTrue(any(abs(p - 44.75) < 0.2 for p in positions), positions)


if __name__ == '__main__':
    unittest.main()
