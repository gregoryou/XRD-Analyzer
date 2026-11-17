"""Small deterministic scientific smoke tests; no MainWindow or user files."""
import math
import os
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, os.environ.get('XRD_TEST_APP_DIR', str(Path(__file__).resolve().parents[1])))
import XRD_Analyzer_GUI as app


class Value:
    def __init__(self, value):
        self._value = value

    def value(self):
        return self._value

    def isChecked(self):
        return bool(self._value)


def detector_settings():
    return SimpleNamespace(
        baseline_cb=Value(False), min_distance=Value(0.01),
        prominence=Value(1.0), min_snr=Value(0.0),
        min_rel_height=Value(1.0), min_fwhm=Value(0.0),
        max_peaks=Value(10), phase_tol=Value(0.35),
    )


class NumericSmokeTests(unittest.TestCase):
    def test_bragg_known_angle(self):
        self.assertAlmostEqual(app.bragg_d(60.0, 1.5), 1.5, places=12)
        self.assertAlmostEqual(app.d_to_2theta(1.5, 1.5), 60.0, places=12)

    def test_bragg_roundtrip_positive_physical_domain(self):
        for wavelength in (0.5, 1.5406, 2.2897):
            for angle in (1.0, 20.0, 60.0, 120.0, 179.0):
                with self.subTest(wavelength=wavelength, angle=angle):
                    self.assertAlmostEqual(
                        app.d_to_2theta(app.bragg_d(angle, wavelength), wavelength),
                        angle, places=9,
                    )

    def test_basic_lattice_spacings(self):
        self.assertAlmostEqual(app.cubic_d(4.0, 2, 0, 0), 2.0)
        self.assertAlmostEqual(app.hex_d(3.0, 5.0, 0, 0, 2), 2.5)
        self.assertAlmostEqual(app.tetragonal_d(3.0, 5.0, 0, 0, 2), 2.5)
        self.assertAlmostEqual(app.hex_d(3.0, 5.0, 1, 0, 0), 3.0 * math.sqrt(3) / 2)

    def test_uniform_gaussian_peak_width(self):
        x = np.linspace(20.0, 30.0, 2001)
        sigma = 0.2
        y = 100.0 * np.exp(-0.5 * ((x - 25.0) / sigma) ** 2)
        dataset = app.Dataset(Path('synthetic.dat'), 'synthetic', '', x, y)
        app.MainWindow._analyze_dataset_peaks_only(detector_settings(), dataset)
        self.assertEqual(len(dataset.peaks), 1)
        self.assertAlmostEqual(dataset.peaks[0].two_theta, 25.0, places=10)
        self.assertAlmostEqual(dataset.peaks[0].fwhm_deg, 2.0 * math.sqrt(2.0 * math.log(2.0)) * sigma, delta=0.0001)


if __name__ == '__main__':
    unittest.main()
