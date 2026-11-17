import unittest
from test_numeric_smoke import app, np, Path, Value, detector_settings


class AngularDistanceTests(unittest.TestCase):
    def run_detector(self, x, centers, distance):
        y = sum(height*np.exp(-.5*((x-center)/.025)**2) for center,height in centers)
        dataset = app.Dataset(Path('angles.dat'), 'angles', '', x, y)
        settings = detector_settings(); settings.min_distance = Value(distance)
        app.MainWindow._analyze_dataset_peaks_only(settings, dataset)
        return [p.two_theta for p in dataset.peaks]

    def test_uniform_separated(self):
        x = np.linspace(20,40,4001)
        np.testing.assert_allclose(self.run_detector(x, [(35,100),(35.6,80)], .5), [35,35.6])

    def test_uniform_close_keeps_stronger(self):
        x = np.linspace(20,40,4001)
        np.testing.assert_allclose(self.run_detector(x, [(35,80),(35.3,100)], .5), [35.3])

    def test_dense_background_sparse_peaks_keeps_separated(self):
        x = np.r_[np.arange(20,30,.005),np.arange(30,40,.05)]
        np.testing.assert_allclose(self.run_detector(x, [(35,100),(35.6,80)], .5), [35,35.6], atol=1e-8)

    def test_sparse_background_dense_peaks_suppresses_close(self):
        x = np.r_[np.arange(1,34,.05),np.arange(34,35,.005),np.arange(35,80,.05)]
        np.testing.assert_allclose(self.run_detector(x, [(34.3,80),(34.6,100)], .5), [34.6], atol=1e-8)

    def test_boundary_ties_and_width_semantics(self):
        x=np.array([20.,20.5,21.])
        np.testing.assert_array_equal(app.angular_distance_selection([0,1,2],x,[5,5,5],.5),[0,1,2])
        p=app.PeakResult(1,30,2,1,100,.2)
        self.assertEqual(p.width_half_prominence_deg,.2)
        rec=app.MainWindow._serialize_peak(None,p)
        self.assertEqual(rec['width_definition'],'scipy_peak_width_half_prominence')
        self.assertEqual(rec['width_half_prominence_deg'],rec['fwhm_deg'])


if __name__ == '__main__': unittest.main()
