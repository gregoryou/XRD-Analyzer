"""Provisional ferrite label for a lone, strong, otherwise unexplained peak (internal database only)."""
import math
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

sys.path.insert(0, os.environ.get('XRD_TEST_APP_DIR', str(Path(__file__).resolve().parents[1])))
import XRD_Analyzer_GUI as app

LAM = 1.540598


def two_theta(a, n):
    return 2 * math.degrees(math.asin(LAM * math.sqrt(n) / (2 * a)))


def make_dataset(extra_peaks=((44.52, 3779.0, 33.0),), a=3.592, ferrite_200=0.0):
    x = np.arange(20.0, 120.0, 0.05)
    y = np.full_like(x, 5000.0)
    peaks = []
    for n, height in ((3, 14000.0), (4, 3300.0), (8, 1500.0), (11, 1900.0), (12, 1000.0)):
        tt = two_theta(a, n)
        i = int(np.argmin(abs(x - tt)))
        peaks.append(app.PeakResult(i, float(tt), app.bragg_d(tt, LAM), height, 100 * height / 14000.0, 0.4, snr=height / 110.0))
    for tt, height, snr in extra_peaks:
        i = int(np.argmin(abs(x - tt)))
        peaks.append(app.PeakResult(i, float(tt), app.bragg_d(tt, LAM), height, 100 * height / 14000.0, 0.5, snr=snr))
    if ferrite_200:
        y += ferrite_200 * np.exp(-0.5 * ((x - 64.8) / 0.3) ** 2)
    ds = app.Dataset(Path('316L.udf'), 'material 316L', '', x, y, LAM)
    ds.baseline = np.full_like(x, 5000.0)
    ds.peaks = sorted(peaks, key=lambda p: p.two_theta)
    return ds


class TentativeFerriteTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='xrd-tentative-')
        self.addCleanup(self.temp.cleanup)
        home = patch.object(Path, 'home', return_value=Path(self.temp.name))
        home.start(); self.addCleanup(home.stop)
        self.db = app.PhaseDatabase()
        self.addCleanup(self.db.conn.close)

    def assigned(self, ds, elements=('Cr', 'Fe', 'Ni')):
        app.assign_phases(ds, self.db, 0.35, set(elements), min_phase_matches=2, min_pattern_rel_intensity=2.0, max_pattern_peaks=30)
        return ds

    def peak_at(self, ds, tt):
        return min(ds.peaks, key=lambda p: abs(p.two_theta - tt))

    def test_regular_assignment_leaves_the_lone_peak_unassigned(self):
        ds = self.assigned(make_dataset())
        self.assertEqual(self.peak_at(ds, 44.52).phase, 'unassigned')        # the author's multi-reflection rule stands
        self.assertTrue(self.peak_at(ds, two_theta(3.592, 3)).phase.lower().startswith('austenite'))

    def test_lone_strong_peak_gets_a_provisional_ferrite_label(self):
        ds = self.assigned(make_dataset())
        records = app.assign_tentative_single_reflections(ds, self.db, 0.35, {'Cr', 'Fe', 'Ni'})
        peak = self.peak_at(ds, 44.52)
        self.assertEqual(len(records), 1)
        self.assertIn('δ ferrite', peak.phase)                 # austenite is present: the delta-ferrite entry is preferred
        self.assertEqual((peak.hkl, peak.phase_source), ('110', app.TENTATIVE_SOURCE))
        self.assertLess(peak.delta_2theta, 0.1)
        self.assertEqual(ds.metadata['tentative_phases'][0]['hkl'], '110')
        self.assertIn('200', ds.metadata['tentative_phases'][0]['note'])

    def test_austenite_peaks_are_never_changed(self):
        ds = self.assigned(make_dataset())
        before = [(p.phase, p.hkl, p.phase_source) for p in ds.peaks if p.phase != 'unassigned']
        app.assign_tentative_single_reflections(ds, self.db, 0.35, {'Cr', 'Fe', 'Ni'})
        after = [(p.phase, p.hkl, p.phase_source) for p in ds.peaks if p.phase_source != app.TENTATIVE_SOURCE and p.phase != 'unassigned']
        self.assertEqual(before, after)

    def test_labels_show_that_the_assignment_is_provisional(self):
        ds = self.assigned(make_dataset())
        app.assign_tentative_single_reflections(ds, self.db, 0.35, {'Cr', 'Fe', 'Ni'})
        peak = self.peak_at(ds, 44.52)
        short = app.format_peak_label(ds, peak, app.PEAK_LABEL_STYLE_SHORT)
        full = app.format_peak_label(ds, peak, app.PEAK_LABEL_STYLE_FULL)
        self.assertIn('?', short)
        self.assertIn(app.TENTATIVE_SUFFIX, full)
        self.assertEqual(app._hkl_for_display(peak, app.HKL_DISPLAY_INTERNAL), '(110)')

    def test_weak_noisy_or_distant_peaks_are_not_labelled(self):
        for peak in ((44.52, 3779.0, 6.0), (41.0, 3779.0, 33.0), (47.0, 3779.0, 33.0)):
            with self.subTest(peak=peak):
                ds = self.assigned(make_dataset(extra_peaks=(peak,)))
                self.assertEqual(app.assign_tentative_single_reflections(ds, self.db, 0.35, {'Cr', 'Fe', 'Ni'}), [])
                self.assertEqual(self.peak_at(ds, peak[0]).phase_source, '')

    def test_requires_a_compatible_chemistry_and_steel_phases(self):
        ds = self.assigned(make_dataset(), elements=('Cu', 'Ni'))
        self.assertEqual(app.assign_tentative_single_reflections(ds, self.db, 0.35, {'Cu', 'Ni'}), [])
        self.assertEqual(app.assign_tentative_single_reflections(ds, self.db, 0.35, set()), [])
        titanium = self.assigned(make_dataset(extra_peaks=((40.2, 3779.0, 33.0),)), elements=('Ti', 'Al', 'V'))
        self.assertEqual(app.assign_tentative_single_reflections(titanium, self.db, 0.35, {'Ti', 'Al', 'V'}), [])

    def test_second_reflection_is_reported_but_not_required(self):
        absent = self.assigned(make_dataset())
        app.assign_tentative_single_reflections(absent, self.db, 0.35, {'Cr', 'Fe', 'Ni'})
        self.assertIn('brak potwierdzenia', absent.metadata['tentative_phases'][0]['note'])
        present = self.assigned(make_dataset(ferrite_200=700.0))
        app.assign_tentative_single_reflections(present, self.db, 0.35, {'Cr', 'Fe', 'Ni'})
        self.assertIn('zgodne', present.metadata['tentative_phases'][0]['note'])

    def test_repeatable_and_cleared_with_the_assignments(self):
        ds = self.assigned(make_dataset())
        first = app.assign_tentative_single_reflections(ds, self.db, 0.35, {'Cr', 'Fe', 'Ni'})
        again = app.assign_tentative_single_reflections(ds, self.db, 0.35, {'Cr', 'Fe', 'Ni'})
        self.assertEqual(len(first), 1)
        self.assertEqual(again, [])                        # the peak is no longer unassigned
        app.clear_phase_assignments(ds)
        self.assertNotIn('tentative_phases', ds.metadata)
        self.assertTrue(all(p.phase == 'unassigned' and not p.phase_source for p in ds.peaks))


if __name__ == '__main__':
    unittest.main()
