"""Analytical and synthetic regressions for the independent profile module."""
import json
import math
from pathlib import Path
import sys
import unittest

import numpy as np
from scipy.special import voigt_profile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from xrd_profile_analysis import (
    CU_KALPHA1_A, CU_KALPHA2_A, InstrumentBroadeningModel, PeakFitResult,
    bragg_d, correct_instrument_broadening, coupled_profile,
    evaluate_peak_profile, fit_cubic_lattice, fit_peak_profile,
    fit_williamson_hall, kalpha2_position, pseudo_voigt, scherrer_size,
    tch_fwhm_eta,
)


class ProfileAnalysisTests(unittest.TestCase):
    def test_gaussian_and_lorentz_limits_and_half_maximum(self):
        for g, l, expected_eta in ((.2, 0, 0), (0, .2, 1)):
            h, eta = tch_fwhm_eta(g, l)
            self.assertAlmostEqual(h, .2, places=12)
            self.assertAlmostEqual(eta, expected_eta, places=12)
            profile = pseudo_voigt(np.array([-.1, 0, .1]), 0, h, eta)
            self.assertAlmostEqual(profile[0]/profile[1], .5, places=12)
            self.assertAlmostEqual(profile[2]/profile[1], .5, places=12)

    def test_tch_width_agrees_with_independent_exact_voigt(self):
        # scipy.special.voigt_profile convolves the profiles, independent of TCH.
        for g, l in ((.18, .04), (.08, .15), (.12, .12)):
            x = np.linspace(0, 1, 100001)
            y = voigt_profile(x, g/(2*np.sqrt(2*np.log(2))), l/2)
            half_width = np.interp(y[0]/2, y[::-1], x[::-1])
            self.assertLess(abs(tch_fwhm_eta(g, l)[0]/(2*half_width)-1), .003)

    def test_pseudo_voigt_normalization(self):
        x = np.linspace(-1000, 1000, 500001)
        self.assertAlmostEqual(np.trapezoid(pseudo_voigt(x, 0, .2, .4), x), 1, places=4)

    def test_coupled_doublet_bragg_geometry(self):
        center = 91.2
        second = kalpha2_position(center)
        d1 = CU_KALPHA1_A/(2*np.sin(np.deg2rad(center/2)))
        d2 = CU_KALPHA2_A/(2*np.sin(np.deg2rad(second/2)))
        self.assertAlmostEqual(d1, d2, places=13)
        self.assertGreater(second, center)
        with self.assertRaises(ValueError):
            kalpha2_position(179)

    def test_doublet_fit_irregular_grid_background_covariance_and_raw_unchanged(self):
        rng = np.random.default_rng(1206)
        x = np.sort(rng.uniform(42, 45, 850))
        exact = coupled_profile(x, 43.4, 500, .17, .065, 71, 4, 43.5)
        y = exact+rng.normal(0, 1.5, len(x))
        original_x, original_y = x.copy(), y.copy()
        fit = fit_peak_profile(x, y, 43.35, sigma_y=np.full(len(x), 1.5), hkl="111")
        self.assertTrue(fit.fit_success)
        self.assertEqual(fit.status, "CONDITIONAL")
        self.assertAlmostEqual(fit.two_theta_kalpha1_deg, 43.4, delta=.001)
        self.assertAlmostEqual(fit.gaussian_fwhm_deg, .17, delta=.002)
        self.assertAlmostEqual(fit.lorentzian_fwhm_deg, .065, delta=.002)
        self.assertAlmostEqual(fit.fwhm_obs_deg, tch_fwhm_eta(.17, .065)[0], delta=.001)
        self.assertGreater(fit.u_two_theta_deg, 0)
        self.assertGreater(fit.u_fwhm_obs_deg, 0)
        self.assertAlmostEqual(fit.reduced_chi_squared, 1, delta=.2)
        self.assertGreater(fit.r_squared, .999)
        self.assertTrue(np.array_equal(x, original_x))
        self.assertTrue(np.array_equal(y, original_y))
        saved = json.loads(json.dumps(fit.to_dict(), allow_nan=False))
        restored = PeakFitResult(**saved)
        np.testing.assert_allclose(evaluate_peak_profile(x, restored), fit.fitted_y, rtol=1e-13)

    def test_single_gaussian_fit_without_doublet(self):
        # Generate independently, without the pseudo_voigt/coupled_profile helper.
        x = np.linspace(35, 37, 501)
        y = 30+600*np.exp(-4*np.log(2)*((x-36.1)/.23)**2)
        fit = fit_peak_profile(x, y, 36.1, kalpha2_ratio=0)
        self.assertAlmostEqual(fit.two_theta_kalpha1_deg, 36.1, delta=1e-5)
        self.assertAlmostEqual(fit.fwhm_obs_deg, .23, delta=1e-4)
        self.assertLess(fit.eta, .001)

    def test_single_lorentzian_fit_without_doublet(self):
        x = np.linspace(35, 37, 501)
        y = 30+600/(1+4*((x-36.1)/.23)**2)
        fit = fit_peak_profile(x, y, 36.1, kalpha2_ratio=0)
        self.assertAlmostEqual(fit.two_theta_kalpha1_deg, 36.1, delta=1e-5)
        self.assertAlmostEqual(fit.fwhm_obs_deg, .23, delta=1e-4)
        self.assertGreater(fit.eta, .999)

    def test_invalid_profiles_rejected(self):
        x = np.linspace(40, 42, 20)
        for bad_x, bad_y in ((x[::-1], np.ones(20)), (x, np.full(20, np.nan)),
                             (x[:8], np.ones(8)), (np.ones(20)*40, np.arange(20))):
            with self.assertRaises(ValueError):
                fit_peak_profile(bad_x, bad_y)
        with self.assertRaises(ValueError):
            fit_peak_profile(x, np.zeros(20))

    def test_bragg_uncertainties_units(self):
        d, ud = bragg_d(60, 1.5, .02, .0001)
        self.assertAlmostEqual(d, 1.5, places=13)
        expected = math.hypot(.0001, 1.5*math.cos(math.pi/6)*math.radians(.02))
        self.assertAlmostEqual(ud, expected, places=13)

    def _fitted_observed(self):
        x = np.linspace(42, 45, 750)
        gobs = math.hypot(.2, .08)
        lobs = .06+.035
        y = coupled_profile(x, 43.4, 800, gobs, lobs, 24, 1, 43.5)
        return fit_peak_profile(x, y, 43.4, sigma_y=1)

    def test_instrument_component_correction_and_uncertainty(self):
        fit = self._fitted_observed()
        ipf = InstrumentBroadeningModel([40, 50], [.08, .08], [.035, .035],
                                       [.002, .002], [.001, .001], standard="Synthetic")
        correction = correct_instrument_broadening(fit, ipf)
        self.assertEqual(correction.status, "CONDITIONAL")
        self.assertAlmostEqual(correction.gaussian_corrected_deg, .2, places=6)
        self.assertAlmostEqual(correction.lorentzian_corrected_deg, .06, places=6)
        self.assertAlmostEqual(correction.fwhm_corrected_deg, tch_fwhm_eta(.2, .06)[0], places=6)
        self.assertGreater(correction.u_fwhm_corrected_deg, 0)
        self.assertTrue(correction.uncertainty_complete)

    def test_instrument_no_extrapolation_or_negative_deconvolution(self):
        fit = self._fitted_observed()
        ipf = InstrumentBroadeningModel([44, 50], [.08, .08], [.035, .035])
        result = correct_instrument_broadening(fit, ipf)
        self.assertEqual(result.status, "INSUFFICIENT_DATA")
        self.assertIsNone(result.fwhm_corrected_deg)
        too_broad = InstrumentBroadeningModel([40, 50], [.4, .4], [.035, .035])
        result = correct_instrument_broadening(fit, too_broad)
        self.assertEqual(result.status, "REJECTED")
        self.assertIsNone(result.fwhm_corrected_deg)
        with self.assertRaises(ValueError):
            ipf.at(51)

    def test_instrument_interpolation_errors_and_input_validation(self):
        ipf = InstrumentBroadeningModel([30, 70], [.1, .2], [.02, .06], [.01, .02], [.02, .02])
        g, l, ug, ul = ipf.at(50)
        self.assertAlmostEqual(g, .15)
        self.assertAlmostEqual(l, .04)
        self.assertAlmostEqual(ug, math.sqrt(.25*.01**2+.25*.02**2))
        self.assertAlmostEqual(ul, .02/math.sqrt(2))
        for angles, widths in (([40, 40], [.1, .1]), ([40, 50], [-.1, .1])):
            with self.assertRaises(ValueError):
                InstrumentBroadeningModel(angles, widths, [.1, .1])

    def test_scherrer_known_size_and_angstrom_to_nm(self):
        D = 75
        angle = 62
        beta = math.degrees(.9*.1540598/(D*math.cos(math.radians(angle/2))))
        result = scherrer_size(angle, beta, u_beta_deg=.01*beta, u_two_theta_deg=.002)
        self.assertAlmostEqual(result.D_nm, D, places=10)
        expected = D*math.hypot(.01, math.tan(math.radians(angle/2))*math.radians(.002)/2)
        self.assertAlmostEqual(result.u_D_nm, expected, places=10)
        unresolved = scherrer_size(angle, beta, u_beta_deg=beta)
        self.assertIsNone(unresolved.D_nm)
        self.assertEqual(unresolved.status, "INSUFFICIENT_DATA")
        with self.assertRaises(ValueError):
            scherrer_size(angle, -1)

    def test_wh_weighted_known_size_strain_and_leave_one_out(self):
        angles = np.array([35, 48, 65, 82, 108.])
        theta = np.deg2rad(angles/2)
        D, epsilon = 80, .0015
        y = .9*.1540598/D+4*epsilon*np.sin(theta)
        beta = np.rad2deg(y/np.cos(theta))
        result = fit_williamson_hall(angles, beta, sigma_beta_deg=np.full(5, .0004))
        self.assertAlmostEqual(result.epsilon, epsilon, places=12)
        self.assertAlmostEqual(result.D_nm, D, places=10)
        self.assertGreater(result.u_epsilon, 0)
        self.assertGreater(result.u_D_nm, 0)
        self.assertAlmostEqual(sum(result.leverage), 2, places=10)
        self.assertEqual(len(result.leave_one_out), 5)
        self.assertEqual(result.status, "CONDITIONAL")
        self.assertFalse(result.valid_for_quantitative_interpretation)
        for item in result.leave_one_out:
            self.assertAlmostEqual(item["D_nm"], D, places=9)
        json.dumps(result.to_dict(), allow_nan=False)

    def test_wh_negative_intercept_and_insufficient_reflections(self):
        angles = np.array([40, 60, 80, 100.])
        theta = np.deg2rad(angles/2)
        beta = np.rad2deg((-.0004+.003*4*np.sin(theta))/np.cos(theta))
        fit = fit_williamson_hall(angles, beta, sigma_beta_deg=.001)
        self.assertLess(fit.intercept, 0)
        self.assertIsNone(fit.D_nm)
        self.assertEqual(fit.status, "REJECTED")
        for bad_angles, bad_beta in ((angles[:2], beta[:2]), (np.array([40, 40, 40]), beta[:3]),
                                    (angles, -beta)):
            with self.assertRaises(ValueError):
                fit_williamson_hall(bad_angles, bad_beta)

    def test_wh_negative_strain_and_unresolved_intercept(self):
        angles = np.array([40, 60, 80, 100.])
        theta = np.deg2rad(angles/2)
        beta = np.rad2deg((.01-.0003*4*np.sin(theta))/np.cos(theta))
        fit = fit_williamson_hall(angles, beta, sigma_beta_deg=.001)
        self.assertLess(fit.epsilon, 0)
        self.assertEqual(fit.status, "REJECTED")
        beta = np.rad2deg((.000001+.003*4*np.sin(theta))/np.cos(theta))
        fit = fit_williamson_hall(angles, beta, sigma_beta_deg=.01)
        self.assertIsNone(fit.D_nm)

    def test_wh_weighting_reduces_known_bad_point_influence(self):
        angles = np.array([35, 48, 65, 82, 108.])
        theta = np.deg2rad(angles/2)
        beta = np.rad2deg((.002+.0015*4*np.sin(theta))/np.cos(theta))
        beta[-1] += .25
        unweighted = fit_williamson_hall(angles, beta)
        weighted = fit_williamson_hall(angles, beta, sigma_beta_deg=[.001]*4+[1])
        self.assertLess(abs(weighted.epsilon-.0015), abs(unweighted.epsilon-.0015)/100)

    def test_lattice_recovers_cubic_a_and_zero_shift(self):
        hkls = [(1, 1, 1), (2, 0, 0), (2, 2, 0), (3, 1, 1), (2, 2, 2)]
        norms = np.sqrt([sum(v*v for v in hkl) for hkl in hkls])
        angles = np.rad2deg(2*np.arcsin(CU_KALPHA1_A*norms/(2*3.6012)))+.075
        result = fit_cubic_lattice(angles, hkls, sigma_2theta_deg=np.full(5, .005))
        self.assertAlmostEqual(result.a_A, 3.6012, places=10)
        self.assertAlmostEqual(result.zero_shift_deg, .075, places=9)
        self.assertGreater(result.u_a_A, 0)
        self.assertGreater(result.u_zero_shift_deg, 0)
        self.assertLess(result.rms_2theta_deg, 1e-8)
        self.assertEqual(result.status, "CONDITIONAL")
        no_shift = fit_cubic_lattice(angles-.075, ["111", "200", "220", "311", "222"],
                                    fit_zero_shift=False, sigma_2theta_deg=.005)
        self.assertAlmostEqual(no_shift.a_A, 3.6012, places=10)
        self.assertEqual(no_shift.zero_shift_deg, 0)

    def test_lattice_rejects_zero_or_duplicate_families(self):
        with self.assertRaises(ValueError):
            fit_cubic_lattice([40, 60, 80], ["000", "200", "220"])
        with self.assertRaises(ValueError):
            fit_cubic_lattice([40, 41, 42], ["111", "111", "111"])
        with self.assertRaises(ValueError):
            fit_cubic_lattice([40, 60], ["111", "200"])


if __name__ == "__main__":
    unittest.main()
