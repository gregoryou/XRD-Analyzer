import unittest
from test_numeric_smoke import app, np


class GemmiCODIntegrationTests(unittest.TestCase):
    def test_fcc_cod_cell_reflections_and_systematic_absences(self):
        self.assertTrue(app.ensure_gemmi(),app.GEMMI_IMPORT_ERROR)
        row=(1,'Cu','copper','','','F m -3 m',225,3.6,3.6,3.6,90,90,90,2026,'synthetic',1)
        refs=app.CODDatabase._row_to_reflections(row,1.540598,20,100)
        angles=[value[0] for value in refs]
        self.assertGreater(len(angles),3)
        expected=np.rad2deg(2*np.arcsin(1.540598*np.sqrt(3)/(2*3.6)))
        self.assertAlmostEqual(angles[0],expected,places=10)
        forbidden=np.rad2deg(2*np.arcsin(1.540598/(2*3.6)))
        self.assertFalse(any(abs(angle-forbidden)<1e-5 for angle in angles))


if __name__=='__main__':unittest.main()
