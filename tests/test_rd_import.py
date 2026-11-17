"""Philips .rd (V3RD) reader.

Most tests use SYNTHETIC files built from the layout seen in the one real file available
(810-byte header, uint16 counts, float64 start/step/end); they verify the validation logic and the
error handling. Optional tests can read a separate experimental silicon scan when the user supplies
tests/data/Si_polycrystal_PW3064.rd. The research measurement is intentionally NOT distributed
in this repository. Compatibility with other instruments or V5RD is NOT established.
"""
import os
import struct
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, os.environ.get('XRD_TEST_APP_DIR', str(Path(__file__).resolve().parents[1])))
import xrd_rd_import as rd
import XRD_Analyzer_GUI as app


def poisson_scan(n=2000, seed=1):
    x = np.linspace(20.0, 100.0, n)
    clean = 400 + 300 * np.exp(-0.5 * ((x - 44.0) / 0.12) ** 2) + 150 * np.exp(-0.5 * ((x - 51.0) / 0.15) ** 2)
    return np.random.default_rng(seed).poisson(clean).astype(np.uint16)


def build_v3(counts, start=20.0, step=0.04, axis_offset=214, axis_kind='<f4', order='sde', convention=0,
             wavelengths=(1.540598, 1.544426), wl_offset=94, extra_blocks=(), trailer=b'', magic=b'V3RD'):
    header = bytearray(810)
    header[0:4] = magic
    header[8:40] = b'Si standard NIST 640  '.ljust(32, b' ')
    n = len(counts)
    size = np.dtype(axis_kind).itemsize
    end = start + (n - convention) * step if convention else start + (n - 1) * step
    values = {'s': start, 'd': step, 'e': end}
    for slot, key in enumerate(order):
        chunk = np.array([values[key]], dtype=axis_kind).tobytes()
        header[axis_offset + slot * size:axis_offset + (slot + 1) * size] = chunk
    if wavelengths:
        for slot, value in enumerate(wavelengths):
            struct.pack_into('<f', header, wl_offset + 4 * slot, value)
    for offset, (s, d, e) in extra_blocks:
        struct.pack_into('<fff', header, offset, s, d, e)
    return bytes(header) + counts.astype('<u2').tobytes() + trailer


class RdReaderTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='xrd-rd-')
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)

    def write(self, payload, name='scan.rd'):
        path = self.base / name
        path.write_bytes(payload)
        return path

    def test_valid_v3_file_is_read_with_axis_counts_and_wavelength(self):
        counts = poisson_scan()
        scan = rd.parse_philips_rd(self.write(build_v3(counts)))
        self.assertEqual(scan.version, '3')
        np.testing.assert_array_equal(scan.y, counts.astype(float))
        np.testing.assert_allclose(scan.x, 20.0 + 0.04 * np.arange(len(counts)), atol=1e-4)
        self.assertAlmostEqual(scan.wavelength_A, 1.540598, places=5)
        self.assertAlmostEqual(scan.wavelength_kalpha2_A, 1.544426, places=5)
        self.assertEqual(scan.anode, 'Cu')
        self.assertIn('Si standard', scan.header_text)
        self.assertTrue(any('heurystycznie' in w for w in scan.warnings))

    def test_axis_block_may_be_float64_in_any_order_and_any_position(self):
        counts = poisson_scan(1500, seed=2)
        for kind, order, offset in (('<f8', 'sde', 300), ('<f4', 'esd', 531), ('<f4', 'dse', 77), ('<f8', 'des', 401)):
            with self.subTest(kind=kind, order=order, offset=offset):
                scan = rd.parse_philips_rd(self.write(build_v3(counts, start=15.2, step=0.0167, axis_offset=offset,
                                                                  axis_kind=kind, order=order)))
                self.assertEqual(len(scan.x), len(counts))
                self.assertAlmostEqual(scan.x[0], 15.2, places=4)
                self.assertAlmostEqual(scan.x[1] - scan.x[0], 0.0167, places=5)

    def test_end_stored_as_end_of_last_step_is_read_with_one_point_less(self):
        # With end = start + N*step the same numbers also read as "N+1 points, end = last point".
        # The reader cannot tell them apart; it takes the longer reading only if the file has the data.
        counts = poisson_scan(1200, seed=3)
        scan = rd.parse_philips_rd(self.write(build_v3(counts, convention=1)))
        self.assertIn(len(scan.x), (len(counts) - 1, len(counts)))
        np.testing.assert_allclose(scan.x[:5], 20.0 + 0.04 * np.arange(5), atol=1e-4)

    def test_a_few_trailing_values_are_ignored_and_reported(self):
        counts = poisson_scan(1800, seed=4)
        scan = rd.parse_philips_rd(self.write(build_v3(counts, trailer=struct.pack('<4H', 1, 2, 3, 4))))
        self.assertEqual(len(scan.y), len(counts))
        self.assertTrue(any('Pominięto 4' in w for w in scan.warnings))

    def test_duplicate_identical_blocks_are_not_ambiguous(self):
        counts = poisson_scan(1000, seed=5)
        n = len(counts)
        block = (20.0, 0.04, 20.0 + (n - 1) * 0.04)
        scan = rd.parse_philips_rd(self.write(build_v3(counts, extra_blocks=((600, block),))))
        self.assertEqual(len(scan.x), n)

    def test_conflicting_blocks_are_refused(self):
        counts = poisson_scan(1000, seed=6)
        n = len(counts)
        other = (25.0, 0.03, 25.0 + (n - 1) * 0.03)
        with self.assertRaises(ValueError) as caught:
            rd.parse_philips_rd(self.write(build_v3(counts, extra_blocks=((600, other),))))
        self.assertIn('kilka różnych', str(caught.exception))

    def test_missing_axis_block_is_refused_with_instructions(self):
        counts = poisson_scan(900, seed=7)
        payload = bytearray(build_v3(counts))
        payload[214:226] = bytes(12)
        with self.assertRaises(ValueError) as caught:
            rd.parse_philips_rd(self.write(bytes(payload)))
        self.assertIn('Oś kątowa nie została odgadnięta', str(caught.exception))
        self.assertIn('xrd_rd_import.py', str(caught.exception))

    def test_header_describing_more_points_than_the_file_holds_is_anchored_at_the_end_with_a_warning(self):
        counts = poisson_scan(1000, seed=8)
        payload = bytearray(build_v3(counts))
        struct.pack_into('<f', payload, 222, 95.0)         # end angle 95 deg implies 1876 points, the file has 1000
        scan = rd.parse_philips_rd(self.write(bytes(payload)))
        self.assertEqual(len(scan.x), 1000)
        self.assertAlmostEqual(scan.x[-1], 95.0, places=4)
        self.assertEqual(scan.layout['anchor'], 'end')
        self.assertTrue(any('zakotwiczona' in w for w in scan.warnings))

    def test_axis_that_matches_no_whole_number_of_points_is_refused(self):
        counts = poisson_scan(1000, seed=8)
        payload = bytearray(build_v3(counts))
        struct.pack_into('<f', payload, 222, 95.013)       # not a whole number of steps from the start
        with self.assertRaises(ValueError):
            rd.parse_philips_rd(self.write(bytes(payload)))

    def test_trimmed_scan_matches_the_layout_of_the_real_file(self):
        full = poisson_scan(1500, seed=14)
        n_skip = 281
        payload = build_v3(full, start=15.046, step=0.0334225372, axis_kind='<f8', order='dse', axis_offset=214)
        trimmed = payload[:810] + full[n_skip:].astype('<u2').tobytes()
        scan = rd.parse_philips_rd(self.write(trimmed))
        n = len(full) - n_skip
        self.assertEqual(len(scan.x), n)
        self.assertAlmostEqual(scan.x[-1], 15.046 + 0.0334225372 * (len(full) - 1), places=4)
        self.assertAlmostEqual(scan.x[0], 15.046 + 0.0334225372 * n_skip, places=4)
        self.assertEqual(scan.layout['missing_points'], n_skip)

    def test_two_different_trimmed_blocks_are_refused(self):
        full = poisson_scan(1500, seed=15)
        payload = bytearray(build_v3(full, start=15.0, step=0.04, axis_offset=214, axis_kind='<f8', order='dse'))
        other = np.array([25.0, 0.03, 25.0 + 0.03 * 1899], dtype='<f8').tobytes()
        payload[600:624] = other
        trimmed = bytes(payload[:810]) + full[300:].astype('<u2').tobytes()
        with self.assertRaises(ValueError) as caught:
            rd.parse_philips_rd(self.write(trimmed))
        self.assertIn('kilka', str(caught.exception))

    def test_real_silicon_standard_file(self):
        real = Path(__file__).resolve().parent / 'data' / 'Si_polycrystal_PW3064.rd'
        if not real.is_file():
            self.skipTest('Optional experimental silicon scan not distributed in open-source repository')
        scan = rd.parse_philips_rd(real)
        self.assertEqual((len(scan.x), scan.version, scan.anode), (2711, '3', 'Cu'))
        self.assertAlmostEqual(scan.wavelength_A, 1.540598, places=6)
        self.assertAlmostEqual(scan.wavelength_kalpha2_A, 1.544426, places=6)
        self.assertAlmostEqual(scan.x[0], 24.4377, places=3)
        self.assertAlmostEqual(scan.x[-1], 115.0128, places=3)
        self.assertEqual(scan.layout['missing_points'], 281)
        self.assertIn('PW3050/60', scan.header_text)
        # the axis is only believable if the silicon reflections come out where they belong
        from xrd_profile_analysis import fit_cubic_lattice, fit_peak_profile
        positions, hkls = [], []
        for hkl, nominal in (('111', 28.442), ('220', 47.302), ('311', 56.122), ('400', 69.129), ('331', 76.375),
                             ('422', 88.029), ('511', 94.951), ('440', 106.706)):
            window = abs(scan.x - nominal) < 0.7
            fit = fit_peak_profile(scan.x[window], scan.y[window], center_hint=nominal, wavelength_A=1.540598,
                                   wavelength_kalpha2_A=1.544426)
            self.assertLess(abs(fit.two_theta_kalpha1_deg - nominal), 0.03, hkl)
            positions.append(fit.two_theta_kalpha1_deg); hkls.append(hkl)
        lattice = fit_cubic_lattice(np.array(positions), hkls, 1.540598)
        self.assertLess(abs(lattice.a_A - 5.43102), 0.001)         # silicon, a = 5.43102 angstrom

    def test_real_silicon_standard_through_load_dataset(self):
        source = Path(__file__).resolve().parent / 'data' / 'Si_polycrystal_PW3064.rd'
        if not source.is_file():
            self.skipTest('Optional experimental silicon scan not distributed in open-source repository')
        ds = app.load_dataset(source)
        self.assertEqual(len(ds.x), 2711)
        self.assertEqual(ds.metadata['radiation'], 'Cu Kα')
        self.assertIn('zakotwiczona', ds.metadata['notes'])

    def test_data_that_do_not_look_like_counts_are_refused(self):
        noise = np.random.default_rng(9).integers(0, 65535, 1500).astype(np.uint16)
        with self.assertRaises(ValueError) as caught:
            rd.parse_philips_rd(self.write(build_v3(noise)))
        self.assertIn('nie wyglądają na zliczenia', str(caught.exception))

    def test_missing_wavelength_is_reported_not_invented(self):
        scan = rd.parse_philips_rd(self.write(build_v3(poisson_scan(1000, seed=10), wavelengths=())))
        self.assertIsNone(scan.wavelength_A)
        self.assertTrue(any('domyślną' in w for w in scan.warnings))

    def test_other_anodes_are_recognised(self):
        scan = rd.parse_philips_rd(self.write(build_v3(poisson_scan(1000, seed=11), wavelengths=(1.788965, 1.792850))))
        self.assertEqual(scan.anode, 'Co')
        self.assertAlmostEqual(scan.wavelength_A, 1.788965, places=5)

    def test_v5_and_foreign_files_are_refused_with_clear_messages(self):
        with self.assertRaises(ValueError) as v5:
            rd.parse_philips_rd(self.write(build_v3(poisson_scan(500), magic=b'V5RD'), 'v5.rd'))
        self.assertIn('V5RD', str(v5.exception))
        with self.assertRaises(ValueError) as foreign:
            rd.parse_philips_rd(self.write(b'RAW1.01' + bytes(2000), 'other.rd'))
        self.assertIn('To nie jest plik Philips RD', str(foreign.exception))
        with self.assertRaises(ValueError):
            rd.parse_philips_rd(self.write(b'V3RD' + bytes(100), 'short.rd'))
        with self.assertRaises(ValueError):
            rd.parse_philips_rd(self.write(b'V3', 'tiny.rd'))

    def test_plausibility_check(self):
        self.assertTrue(rd.count_plausibility(poisson_scan(500).astype(float))[0])
        self.assertFalse(rd.count_plausibility(-np.ones(100))[0])
        self.assertFalse(rd.count_plausibility(np.full(100, 7.0))[0])
        self.assertFalse(rd.count_plausibility(np.r_[np.arange(50.0), [np.nan] * 50])[0])
        self.assertTrue(rd.count_plausibility(np.random.default_rng(1).poisson(3, 200).astype(float))[0])   # low counts

    def test_diagnostics_text(self):
        good = rd.describe_philips_rd(self.write(build_v3(poisson_scan(800, seed=12))))
        self.assertIn('ODCZYT POWIÓDŁ SIĘ', good)
        bad = rd.describe_philips_rd(self.write(build_v3(poisson_scan(800, seed=12), magic=b'V5RD'), 'v5.rd'))
        self.assertIn('ODCZYT NIE POWIÓDŁ SIĘ', bad)
        self.assertIn('Początek pliku (hex)', bad)

    def test_load_dataset_integration(self):
        counts = poisson_scan(1600, seed=13)
        ds = app.load_dataset(self.write(build_v3(counts, start=18.5, step=0.0167), 'Si_standard.rd'))
        self.assertEqual(len(ds.x), len(counts))
        np.testing.assert_allclose(ds.y, counts.astype(float))
        self.assertAlmostEqual(ds.wavelength_A, 1.540598, places=5)
        self.assertIn('Philips RD v3', ds.metadata['source_format'])
        self.assertEqual(ds.metadata['radiation'], 'Cu Kα')
        self.assertEqual(ds.metadata_sources['wavelength_A'], 'file (heurystyka)')
        self.assertEqual(ds.material_type, '')                      # never invented from the file name
        bad = self.base / 'broken.rd'
        bad.write_bytes(b'V3RD' + bytes(2000))
        with self.assertRaises(ValueError):
            app.load_dataset(bad)


if __name__ == '__main__':
    unittest.main()
