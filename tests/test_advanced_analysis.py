"""Numerical regressions for sin²psi and constrained Wilkens whole-profile."""
import copy
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from xrd_advanced_analysis import (
    VALIDATED, CONDITIONAL, INSUFFICIENT_DATA, REJECTED,
    fit_sin2psi, bragg_d_with_uncertainty, wilkens_f_star,
    wilkens_fourier, lognormal_sphere_fourier,
    simulate_wilkens_profile, fit_wilkens_profiles, export_whole_profile_request,
)


class StressTests(unittest.TestCase):
    def setUp(self):
        self.options=dict(xec_half_s2_per_MPa=6e-6,d0_A=2.,
            xec_s1_per_MPa=-1.2e-6,u_half_s2_per_MPa=1e-7,u_d0_A=1e-4,
            xec_hkl=[2,1,1],xec_phase="alpha",xec_source="synthetic independent XEC",
            assumptions_confirmed=True)
        self.rows=[]
        for psi in [0,-15,15,-25,25,-35,35,-45,45,-60,60,-70,70]:
            self.rows.append(dict(psi_deg=psi,phi_deg=20,phase="alpha",hkl=[2,1,1],
                d_A=2*(1-1.2e-6*180+6e-6*325*np.sin(np.deg2rad(psi))**2),u_d_A=2e-5))

    def test_known_tensile_stress_and_confidence(self):
        r=fit_sin2psi(self.rows,**self.options)
        self.assertEqual(r["status"],VALIDATED,r["warnings"])
        self.assertAlmostEqual(r["stress_MPa"],325,places=8)
        self.assertLess(r["stress_ci95_MPa"][0],325)
        self.assertGreater(r["stress_ci95_MPa"][1],325)
        self.assertEqual(r["independent_psi_count"],7)
        self.assertEqual(len(r["psi_splitting"]),6)
        self.assertGreater(r["u_stress_MPa"],0)
        json.dumps(r,allow_nan=False)

    def test_compressive_stress_sign(self):
        for row in self.rows:
            row["d_A"]=2*(1-6e-6*200*np.sin(np.deg2rad(row["psi_deg"]))**2)
        r=fit_sin2psi(self.rows,**self.options)
        self.assertAlmostEqual(r["stress_MPa"],-200,places=8)

    def test_d0_not_replaced_by_intercept(self):
        r=fit_sin2psi(self.rows,**self.options)
        self.assertNotEqual(r["d0_A"],r["intercept_A"])
        self.assertAlmostEqual(r["stress_MPa"],325,places=8)

    def test_missing_inputs_and_repeated_tilts(self):
        self.assertEqual(fit_sin2psi(self.rows)["status"],INSUFFICIENT_DATA)
        self.assertEqual(fit_sin2psi(self.rows[:2],**self.options)["status"],INSUFFICIENT_DATA)
        repeat=[dict(self.rows[0]),dict(self.rows[1]),dict(self.rows[2])]*3
        self.assertEqual(fit_sin2psi(repeat,**self.options)["status"],INSUFFICIENT_DATA)

    def test_mixed_phase_hkl_phi_and_invalid_hkl_rejected(self):
        for key,value in [("phase","beta"),("hkl",[1,1,1]),("phi_deg",200),("hkl",[0,0,0]),("hkl",[1.5,1,1])]:
            rows=copy.deepcopy(self.rows)
            rows[0][key]=value
            r=fit_sin2psi(rows,**self.options)
            self.assertEqual(r["status"],REJECTED,(key,r))
            self.assertIsNone(r["stress_MPa"])

    def test_splitting_and_curvature_reject_plane_stress(self):
        rows=copy.deepcopy(self.rows)
        for row in rows:
            row["d_A"]+=.001*np.sin(2*np.deg2rad(row["psi_deg"]))
        r=fit_sin2psi(rows,**self.options)
        self.assertEqual(r["status"],REJECTED)
        self.assertTrue(any(p["significant"] for p in r["psi_splitting"]))
        self.assertIsNone(r["stress_MPa"])
        rows=copy.deepcopy(self.rows)
        for row in rows:
            row["d_A"]+=.005*np.sin(np.deg2rad(row["psi_deg"]))**4
        self.assertEqual(fit_sin2psi(rows,**self.options)["status"],REJECTED)

    def test_uncertainty_increases_with_xec_uncertainty(self):
        a=fit_sin2psi(self.rows,**self.options)
        options=dict(self.options,u_half_s2_per_MPa=1e-6)
        b=fit_sin2psi(self.rows,**options)
        self.assertGreater(b["u_stress_MPa"],a["u_stress_MPa"])

    def test_angle_to_d_and_uncertainty(self):
        d,u=bragg_d_with_uncertainty(60,.01,1.54)
        self.assertAlmostEqual(d,1.54,places=12)
        step=1e-4
        numeric=(bragg_d_with_uncertainty(60+step,0,1.54)[0]-
                 bragg_d_with_uncertainty(60-step,0,1.54)[0])/(2*step)
        self.assertAlmostEqual(u,abs(numeric)*.01,places=11)
        rows=copy.deepcopy(self.rows)
        for row in rows:
            row["two_theta_deg"]=np.rad2deg(2*np.arcsin(1.54/(2*row.pop("d_A"))))
            row.pop("u_d_A")
            row["u_two_theta_deg"]=.001
            row["wavelength_A"]=1.54
        self.assertAlmostEqual(fit_sin2psi(rows,**self.options)["stress_MPa"],325,places=7)

    def test_missing_uncertainties_conditional_and_mixed_blocked(self):
        rows=copy.deepcopy(self.rows)
        for row in rows: row.pop("u_d_A")
        self.assertEqual(fit_sin2psi(rows,**self.options)["status"],CONDITIONAL)
        rows[0]["u_d_A"]=2e-5
        self.assertEqual(fit_sin2psi(rows,**self.options)["status"],INSUFFICIENT_DATA)


class WilkensKernelTests(unittest.TestCase):
    def test_exact_kernel_against_independent_quadrature(self):
        # Independent quadrature of the arcsin integral; no dilogarithm used.
        for z in [.001,.01,.1,.5,.999]:
            integ=quad(lambda y:np.arcsin(y)/y,0,z,epsabs=1e-12)[0]
            expected=(-np.log(z)+1.75-np.log(2)+256/(45*np.pi*z)
                +2/np.pi*(1-1/(4*z*z))*integ
                -(769/(180*z)+41*z/90+z**3/45)*np.sqrt(1-z*z)/np.pi
                -(11/(12*z*z)+3.5+z*z/3)*np.arcsin(z)/np.pi+z*z/6)
            self.assertAlmostEqual(float(wilkens_f_star(z)),expected,places=9)
        self.assertAlmostEqual(float(wilkens_f_star(1)),1.179209446261134,places=12)
        self.assertAlmostEqual(float(wilkens_f_star(2)),.7041880564639003,places=12)

    def test_continuity_positivity_and_monotonicity(self):
        z=np.geomspace(1e-7,1e5,2001)
        v=wilkens_f_star(z)
        self.assertTrue(np.all(v>0))
        self.assertTrue(np.all(np.diff(v)<0))
        self.assertLess(abs(float(wilkens_f_star(1-1e-7)-wilkens_f_star(1+1e-7))),2e-7)
        left=(wilkens_f_star(1)-wilkens_f_star(1-1e-5))/1e-5
        right=(wilkens_f_star(1+1e-5)-wilkens_f_star(1))/1e-5
        self.assertLess(abs(float(left-right)),2e-5)

    def test_rho_si_and_re_convention(self):
        L=1e-4
        rho,re,b,g,c=1e14,100,.25,5,.2
        a=float(wilkens_fourier(L,rho,re,b,g,c))
        expected=np.exp(-np.pi/2*(rho*1e-18)*b*b*g*g*c*L*L*np.log(re/L))
        self.assertAlmostEqual(a,expected,places=13)
        self.assertEqual(float(wilkens_fourier(0,rho,re,b,g,c)),1)
        np.testing.assert_array_equal(wilkens_fourier([0,10,100],0,re,b,g,c),[1,1,1])

    def test_spherical_size_limit_and_lognormal_mixture(self):
        L=np.array([0,10,25,50,100.])
        np.testing.assert_allclose(lognormal_sphere_fourier(L,50,0),[1,.704,.3125,0,0])
        # Numerically integrate independent number-lognormal spheres, volume weighted.
        median,s=55.,.35
        denominator=median**3*np.exp(4.5*s*s)
        for length in [0,10,50,100]:
            def integrand(z):
                diameter=median*np.exp(s*z)
                ratio=length/diameter
                sphere=max(0.,1-1.5*ratio+.5*ratio**3) if ratio<=1 else 0
                return sphere*diameter**3*np.exp(-z*z/2)/np.sqrt(2*np.pi)
            expected=quad(integrand,-10,10,epsabs=1e-6,points=[np.log(max(length,1e-30)/median)/s] if length else None)[0]/denominator
            self.assertAlmostEqual(float(lognormal_sphere_fourier(length,median,s)),expected,places=8)

    def test_invalid_kernel_inputs(self):
        for z in [0,-1,np.nan,np.inf]:
            with self.assertRaises(ValueError): wilkens_f_star(z)
        with self.assertRaises(ValueError): wilkens_fourier([1],-1,100,.25,5,.2)
        with self.assertRaises(ValueError): lognormal_sphere_fourier([1],50,-.1)


class WholeProfileTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.profiles=[]
        cls.truth=dict(rho_m2=4e14,Re_nm=110,median_nm=55,sigma_ln=.3,burgers_vector_nm=.25)
        for center,hkl,c in [(44,[1,1,1],.2),(64,[2,0,0],.3),(82,[2,2,0],.15),(110,[3,1,1],.24)]:
            p=dict(x_deg=np.linspace(center-1.5,center+1.5,301),hkl=hkl,
                contrast_factor=c,center_deg=center,instrument_gaussian_fwhm_deg=.04,
                instrument_lorentzian_fwhm_deg=.01)
            p["y"]=simulate_wilkens_profile(p,**cls.truth,area=500,background=30,
                                             fourier_points=16385)
            p["sigma_y"]=np.ones(301)
            cls.profiles.append(p)
        cls.options=dict(burgers_vector_nm=.25,phase="synthetic",structure="FCC",multistarts=2,
                         instrument_source="synthetic IPF",contrast_source="specified synthetic C",slip_system="synthetic")

    def test_synthetic_recovery_with_different_fourier_grid(self):
        before=[p["y"].copy() for p in self.profiles]
        r=fit_wilkens_profiles(self.profiles,**self.options)
        self.assertEqual(r["status"],CONDITIONAL,r["warnings"])
        for key in ("rho_m2","Re_nm","median_nm","sigma_ln"):
            self.assertLess(abs(r[key]/self.truth[key]-1),.003,key)
        self.assertLess(r["rho_ci95_m2"][0],self.truth["rho_m2"])
        self.assertGreater(r["rho_ci95_m2"][1],self.truth["rho_m2"])
        self.assertEqual(len(r["profiles"]),4)
        self.assertEqual(len(r["multistart_results"]),2)
        self.assertAlmostEqual(r["M"],r["Re_nm"]*1e-9*np.sqrt(r["rho_m2"]),places=12)
        for p,y in zip(self.profiles,before): np.testing.assert_array_equal(p["y"],y)
        json.dumps(r,allow_nan=False)

    def test_noisy_recovery(self):
        rng=np.random.default_rng(482)
        profiles=copy.deepcopy(self.profiles)
        for p in profiles:
            p["sigma_y"]=np.full_like(p["y"],3.)
            p["y"]+=rng.normal(0,3,len(p["y"]))
        r=fit_wilkens_profiles(profiles,**self.options)
        self.assertEqual(r["status"],CONDITIONAL,r["warnings"])
        self.assertLess(abs(r["rho_m2"]/self.truth["rho_m2"]-1),.03)

    def test_doublet_recovery(self):
        profiles=copy.deepcopy(self.profiles)
        for p in profiles:
            p["y"]=simulate_wilkens_profile(p,**self.truth,area=500,background=30,
                kalpha2_ratio=.5,fourier_points=16385)
        r=fit_wilkens_profiles(profiles,**self.options,kalpha2_ratio=.5)
        self.assertEqual(r["status"],CONDITIONAL,r["warnings"])
        self.assertLess(abs(r["rho_m2"]/self.truth["rho_m2"]-1),.003)

    def test_forward_fft_against_cosine_quadrature(self):
        p=self.profiles[0]
        wavelength=.1540598
        center=p["center_deg"]
        g=2*np.sin(np.deg2rad(center/2))/wavelength
        conv=np.cos(np.deg2rad(center/2))*np.pi/(180*wavelength)
        for index in [150,160,200]:
            x=p["x_deg"][index]
            q=2*np.sin(np.deg2rad(x/2))/wavelength-g
            def integrand(length):
                ipf=np.exp(-2*np.pi**2*(.04*conv/np.sqrt(8*np.log(2)))**2*length**2-np.pi*.01*conv*length)
                return float(lognormal_sphere_fourier(length,55,.3)*ipf*
                    wilkens_fourier(length,4e14,110,.25,g,.2))
            expected=2*quad(integrand,0,2048,weight="cos",wvar=2*np.pi*q,epsabs=1e-7,limit=300)[0]
            jac=np.cos(np.deg2rad(x/2))*np.pi/(180*wavelength)
            observed=(p["y"][index]-30)/500/jac
            self.assertLess(abs(observed-expected)/expected,.001)

    def test_missing_data_rejected_inputs_and_no_fake_density(self):
        r=fit_wilkens_profiles([],phase="x",structure="FCC")
        self.assertEqual(r["status"],INSUFFICIENT_DATA)
        self.assertIsNone(r["rho_m2"])
        for key,value in [("overlap",True),("contrast_factor",-1),("x_deg",np.array([1,2]))]:
            profiles=copy.deepcopy(self.profiles)
            profiles[0][key]=value
            r=fit_wilkens_profiles(profiles,**self.options)
            self.assertEqual(r["status"],REJECTED,(key,r))
            self.assertIsNone(r["rho_m2"])
        profiles=copy.deepcopy(self.profiles)
        del profiles[0]["instrument_gaussian_fwhm_deg"]
        self.assertEqual(fit_wilkens_profiles(profiles,**self.options)["status"],INSUFFICIENT_DATA)
        profiles=copy.deepcopy(self.profiles)
        del profiles[0]["sigma_y"]
        self.assertEqual(fit_wilkens_profiles(profiles,**self.options)["status"],INSUFFICIENT_DATA)

    def test_export_portable_json(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"whole_profile.json"
            export_whole_profile_request(path,self.profiles,**self.options)
            data=json.loads(path.read_text(encoding="utf-8"))
            self.assertEqual(data["schema"],"xrd-analyzer-whole-profile/1")
            self.assertEqual(len(data["profiles"]),4)
            self.assertEqual(data["settings"]["burgers_vector_nm"],.25)


if __name__ == "__main__":
    unittest.main()
