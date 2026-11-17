"""ASCII/CSV import, number parsing and display helpers (no window, no user files)."""
import os
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, os.environ.get('XRD_TEST_APP_DIR', str(Path(__file__).resolve().parents[1])))
import XRD_Analyzer_GUI as app
import xrd_analysis_dialog as dialog

ANGLES = np.round(np.arange(35.0, 35.0 + 0.02 * 8, 0.02), 2)
COUNTS = [1000 + 37 * i for i in range(8)]


def _real_qt_items():
    try:
        return isinstance(app.QTableWidgetItem("1").text(), str)
    except Exception:
        return False


class AsciiImportTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='xrd-import-')
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)

    def write(self, name, lines, encoding='utf-8'):
        path = self.base / name
        path.write_text("\n".join(lines), encoding=encoding)
        return path

    def assert_loaded(self, path, expected_y=COUNTS, expected_x=ANGLES):
        ds = app.parse_csv_generic(path)
        np.testing.assert_allclose(ds.y, expected_y)
        np.testing.assert_allclose(ds.x, expected_x)
        return ds

    def test_comma_separated_csv_keeps_intensities(self):
        # Regression: "35.00,1000" used to be read as x=35.00, y=0.1.
        self.assert_loaded(self.write('a.csv', ["2theta,intensity"] + [f"{x:.2f},{y}" for x, y in zip(ANGLES, COUNTS)]))

    def test_comma_separated_csv_with_decimal_intensities_and_extra_columns(self):
        lines = ["2theta,counts,esd"] + [f"{x:.2f},{y}.5,{5 + i}" for i, (x, y) in enumerate(zip(ANGLES, COUNTS))]
        self.assert_loaded(self.write('b.csv', lines), expected_y=np.array(COUNTS) + 0.5)

    def test_comma_followed_by_space(self):
        self.assert_loaded(self.write('c.csv', ["2theta, intensity"] + [f"{x:.2f}, {y}" for x, y in zip(ANGLES, COUNTS)]))

    def test_semicolon_with_decimal_comma(self):
        lines = ["2theta;intensity"] + [f"{str(x).replace('.', ',')};{y}" for x, y in zip(ANGLES, COUNTS)]
        self.assert_loaded(self.write('d.csv', lines))

    def test_whitespace_with_decimal_comma(self):
        self.assert_loaded(self.write('e.dat', [f"{str(x).replace('.', ',')} {y}" for x, y in zip(ANGLES, COUNTS)]))

    def test_tab_space_and_utf8_bom(self):
        self.assert_loaded(self.write('f.txt', [f"{x:.2f}\t{y}" for x, y in zip(ANGLES, COUNTS)]))
        self.assert_loaded(self.write('g.asc', [f"{x:.2f} {y}" for x, y in zip(ANGLES, COUNTS)]))
        bom = self.base / 'h.csv'
        bom.write_bytes(b'\xef\xbb\xbf' + "\r\n".join(["x,y"] + [f"{x:.2f},{y}" for x, y in zip(ANGLES, COUNTS)]).encode())
        self.assert_loaded(bom)

    def test_header_lines_with_numbers_are_not_data(self):
        lines = ["Step size: 0.02 Time: 1.0", "Scan 1 of 2 range 5 90"] + [f"{x:.2f} {y}" for x, y in zip(ANGLES, COUNTS)]
        ds = self.assert_loaded(self.write('i.dat', lines))
        self.assertEqual(len(ds.x), 8)

    def test_descending_axis_is_reversed_and_noted(self):
        lines = [f"{x:.2f},{y}" for x, y in zip(ANGLES[::-1], COUNTS[::-1])]
        ds = self.assert_loaded(self.write('j.csv', lines))
        self.assertIn('odwrócono', ds.metadata['notes'])

    def test_unsorted_axis_is_not_silently_changed(self):
        lines = [f"{x:.2f},{y}" for x, y in zip(ANGLES, COUNTS)]
        lines[2], lines[5] = lines[5], lines[2]
        ds = app.parse_csv_generic(self.write('k.csv', lines))
        self.assertFalse(np.all(np.diff(ds.x) > 0))

    def test_too_few_rows_is_a_readable_error(self):
        with self.assertRaises(ValueError) as caught:
            app.parse_csv_generic(self.write('l.csv', ["a,b", "hello world", "1,2"]))
        self.assertIn('co najmniej 5 wierszy', str(caught.exception))

    def test_parse_xy_line(self):
        for line, expected in [("35.00,1000", [35.0, 1000.0]), ("35,02 1000", [35.02, 1000.0]),
                               ("2theta,intensity", []), ("# comment 1 2", []), ("", []),
                               ("35.0;1e3;5", [35.0, 1000.0, 5.0]), ("nan 5", []), ("35.00, 1000", [35.0, 1000.0])]:
            with self.subTest(line=line):
                self.assertEqual(app.parse_xy_line(line), expected)


class HelperTests(unittest.TestCase):
    def test_number_or_text(self):
        self.assertEqual(app._number_or_text("12,5"), 12.5)
        self.assertEqual(app._number_or_text("abc"), "abc")
        self.assertEqual(app._number_or_text(""), "")

    def test_safe_filename(self):
        self.assertEqual(app._safe_filename("316L/CuNi: weld?"), "316L_CuNi_ weld")
        self.assertEqual(app._safe_filename("///"), "XRD")

    def test_display_helpers_use_polish_placeholders(self):
        self.assertEqual(app._display_text(None), "nie ustawiono")
        self.assertEqual(app._display_text("not set"), "nie ustawiono")
        self.assertEqual(app._display_number("", "kV"), "nie ustawiono")
        self.assertEqual(app._display_number("40", "kV"), "40 kV")
        self.assertEqual(app._display_number("40 kV", "kV"), "40 kV")  # non-numeric text must not raise

    @unittest.skipUnless(_real_qt_items(), "wymaga prawdziwego PyQt6")
    def test_numeric_sort_item(self):
        a, b, blank = app.NumericSortItem("9.5"), app.NumericSortItem("10.2"), app.NumericSortItem("")
        self.assertTrue(a < b)
        self.assertFalse(b < a)
        self.assertTrue(b < blank)
        self.assertFalse(blank < b)


class DialogHelperTests(unittest.TestCase):
    def test_parse_number_accepts_decimal_comma(self):
        self.assertEqual(dialog.parse_number("43,5", "2θ"), 43.5)
        self.assertEqual(dialog.parse_number(" 1.25 ", "x"), 1.25)
        self.assertIsNone(dialog.parse_number("", "x", allow_empty=True))

    def test_parse_number_errors_are_readable(self):
        with self.assertRaises(ValueError) as empty:
            dialog.parse_number("", "C_hkl dla refleksu 111")
        self.assertIn("Uzupełnij pole: C_hkl dla refleksu 111", str(empty.exception))
        with self.assertRaises(ValueError) as text:
            dialog.parse_number("abc", "Od [°]")
        self.assertIn("nie jest liczbą", str(text.exception))
        with self.assertRaises(ValueError):
            dialog.parse_number("nan", "x")

    def test_format_with_uncertainty(self):
        self.assertEqual(dialog.format_with_uncertainty(0.2318882, 1.796663e-05), "0.231888 ± 0.000018")
        self.assertEqual(dialog.format_with_uncertainty(3.58765, 0.0012), "3.5877 ± 0.0012")
        self.assertEqual(dialog.format_with_uncertainty(43.4, 12.3), "43 ± 12")
        self.assertEqual(dialog.format_with_uncertainty(1.5, None), "1.5")
        self.assertEqual(dialog.format_with_uncertainty(None, 0.1), "—")
        self.assertNotIn("e-", dialog.format_with_uncertainty(0.5, 1e-9))

    def test_format_value_lists_and_digits(self):
        self.assertEqual(dialog.format_value(0.2318882, 5), "0.23189")
        self.assertEqual(dialog.format_value([220.123456789, 280.0], 5), "[220.12; 280]")
        self.assertEqual(dialog.format_value(None), "—")
        self.assertEqual(dialog.format_value(float('nan')), "—")

    def test_read_csv_records_semicolon_decimal_comma_and_header_delimiter(self):
        with tempfile.TemporaryDirectory(prefix='xrd-csv-') as folder:
            path = Path(folder) / 'series.csv'
            path.write_text("psi_deg; phi_deg; two_theta_deg; u_two_theta_deg; hkl; phase\n"
                            "0;0;43,5012;0,001;311;FCC test\n"
                            "30;0;43,4800;0,001;311;FCC test\n", encoding='utf-8')
            rows = dialog.read_csv_records(path)
            self.assertEqual(rows[0]['two_theta_deg'], "43.5012")
            self.assertEqual(rows[1]['u_two_theta_deg'], "0.001")
            self.assertEqual(rows[1]['phase'], "FCC test")
            comma = Path(folder) / 'comma.csv'
            comma.write_text("psi_deg,phi_deg\n0,0\n", encoding='utf-8')
            self.assertEqual(dialog.read_csv_records(comma)[0]['psi_deg'], "0")
            single = Path(folder) / 'single.csv'
            single.write_text("a\n1\n", encoding='utf-8')
            with self.assertRaises(ValueError):
                dialog.read_csv_records(single)


if __name__ == '__main__':
    unittest.main()
