"""Residual stress and constrained whole-profile dislocation analysis.

Pure numerical module: no Qt, file execution or mutation of measured arrays.
See SCIENTIFIC_ADVANCED.md for conventions, sources and model limitations.
"""
from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from scipy.optimize import least_squares
from scipy.special import erfc, spence
from scipy.stats import chi2, norm, t

VALIDATED = "VALIDATED"
CONDITIONAL = "CONDITIONAL"
INSUFFICIENT_DATA = "INSUFFICIENT_DATA"
REJECTED = "REJECTED"


def _result(status, method, warnings=None, **values):
    return dict(status=status, method=method, warnings=list(warnings or []), **values)


def _number(value, name, positive=False, nonnegative=False):
    v = float(value)
    if not np.isfinite(v) or (positive and v <= 0) or (nonnegative and v < 0):
        raise ValueError(f"Nieprawidłowe {name}.")
    return v


def _hkl(value):
    if isinstance(value, (list, tuple)) and len(value) == 3:
        if any(not np.isfinite(float(i)) or int(i) != float(i) for i in value):
            raise ValueError("Indeksy hkl muszą być całkowite.")
        vals = value
    else:
        text = str(value).strip()
        if re.search(r"[^\d\s,+\-()\[\]{}]", text):
            raise ValueError("Nieprawidłowy zapis indeksów hkl.")
        vals = re.findall(r"-?\d+", text)
        if len(vals) == 1 and len(vals[0]) == 3 and vals[0].isdigit():
            vals = list(vals[0])
    if len(vals) != 3:
        raise ValueError("hkl musi zawierać trzy indeksy, np. [2, 1, 1].")
    hkl = tuple(int(i) for i in vals)
    if hkl == (0, 0, 0):
        raise ValueError("hkl=(0,0,0) jest niedopuszczalne.")
    return hkl


def bragg_d_with_uncertainty(two_theta_deg, u_two_theta_deg, wavelength_A,
                             u_wavelength_A=0.0):
    """Independent peak-position and wavelength standard uncertainties."""
    tt = _number(two_theta_deg, "2theta")
    if not 0 < tt < 180:
        raise ValueError("2theta musi należeć do (0,180) stopni.")
    ut = _number(u_two_theta_deg, "u(2theta)", nonnegative=True)
    w = _number(wavelength_A, "lambda", positive=True)
    uw = _number(u_wavelength_A, "u(lambda)", nonnegative=True)
    theta = np.deg2rad(tt / 2)
    d = w / (2 * np.sin(theta))
    u = np.hypot(uw / (2 * np.sin(theta)),
                 w * np.cos(theta) / (4 * np.sin(theta) ** 2) * np.deg2rad(ut))
    return float(d), float(u)


def fit_sin2psi(measurements: Sequence[dict], *, xec_half_s2_per_MPa=None,
                d0_A=None, xec_s1_per_MPa=None, u_half_s2_per_MPa=0.0,
                u_d0_A=0.0, xec_hkl="", xec_phase="", xec_source="",
                assumptions_confirmed=False):
    """Fit d = a + b sin²(psi), sigma_phi = b/(d0*(S2/2)).

    At least three distinct sin²psi are required. d0 is independently measured;
    it is never substituted by the fitted intercept. Same phase, hkl and phi
    (modulo 360°, preserving signed-psi semantics) are required. psi-splitting
    and lack of linearity reject the plane-stress result, preserving diagnostics.
    """
    method = "sin2psi_plane_stress"
    empty = dict(stress_MPa=None, u_stress_MPa=None, stress_ci95_MPa=None)
    missing = []
    if xec_half_s2_per_MPa is None:
        missing.append("Brak rentgenowskiej stałej S2/2 [1/MPa] dla badanego hkl.")
    if d0_A is None:
        missing.append("Brak niezależnej odległości beznaprężeniowej d0 [Å].")
    if len(measurements) < 3:
        missing.append("Wymagane co najmniej trzy niezależne wartości sin²ψ.")
    if missing:
        return _result(INSUFFICIENT_DATA, method, missing, **empty)
    try:
        q = _number(xec_half_s2_per_MPa, "S2/2", positive=True)
        uq = _number(u_half_s2_per_MPa, "u(S2/2)", nonnegative=True)
        d0 = _number(d0_A, "d0", positive=True)
        ud0 = _number(u_d0_A, "u(d0)", nonnegative=True)
        if xec_s1_per_MPa is not None:
            _number(xec_s1_per_MPa, "S1")
        ds, uds, psis, phis, hkls, phases = [], [], [], [], [], []
        for m in measurements:
            hkls.append(_hkl(m.get("hkl", "")))
            phases.append(str(m.get("phase", "")).strip())
            psi = _number(m["psi_deg"], "psi")
            if not -90 <= psi <= 90:
                raise ValueError("psi poza zakresem [-90,90] stopni.")
            psis.append(psi)
            phis.append(_number(m["phi_deg"], "phi") % 360)
            if "d_A" in m:
                d = _number(m["d_A"], "d", positive=True)
                u = _number(m.get("u_d_A", 0), "u(d)", nonnegative=True)
            else:
                d, u = bragg_d_with_uncertainty(m["two_theta_deg"],
                    m.get("u_two_theta_deg", 0), m["wavelength_A"],
                    m.get("u_wavelength_A", 0))
            ds.append(d)
            uds.append(u)
        if len(set(hkls)) != 1 or len(set(phases)) != 1 or not phases[0]:
            raise ValueError("Wszystkie pomiary muszą dotyczyć tego samego hkl i nazwanej fazy.")
        if any(abs((v - phis[0] + 180) % 360 - 180) > 1e-6 for v in phis):
            raise ValueError("Analiza wymaga jednego azymutu phi; nie mieszaj serii phi.")
        if xec_hkl and _hkl(xec_hkl) != hkls[0]:
            raise ValueError("XEC przypisano do innego hkl.")
        if xec_phase and str(xec_phase).strip() != phases[0]:
            raise ValueError("XEC przypisano do innej fazy.")
    except (ValueError, KeyError, TypeError, OverflowError) as e:
        return _result(REJECTED, method, [str(e)], **empty)
    psi = np.asarray(psis)
    x = np.sin(np.deg2rad(psi)) ** 2
    y, u = np.asarray(ds), np.asarray(uds)
    n_unique = len(np.unique(np.round(x, 10)))
    if n_unique < 3 or np.ptp(x) < 1e-6:
        return _result(INSUFFICIENT_DATA, method,
            ["Wymagane co najmniej trzy różne wartości sin²ψ i niezerowy zakres."],
            independent_psi_count=n_unique, **empty)
    if np.any(u == 0) and np.any(u > 0):
        return _result(INSUFFICIENT_DATA, method,
            ["Niepewności muszą być podane dla wszystkich punktów albo dla żadnego."], **empty)
    warnings = []
    weighted = bool(np.all(u > 0))
    if not weighted:
        warnings.append("Brak niepewności pozycji: regresja nieważona, niepewność z residuów.")
    if not xec_source or not xec_hkl or not xec_phase:
        warnings.append("Niepełna identyfikacja źródła, fazy lub hkl stałych XEC.")
    if not assumptions_confirmed:
        warnings.append("Należy potwierdzić płaski stan naprężenia, brak gradientu d0 i adekwatność XEC.")
    if uq == 0 or ud0 == 0:
        warnings.append("Nie uwzględniono pełnej niepewności XEC/d0 (co najmniej jedna wynosi zero).")
    if n_unique < 5 or np.ptp(x) < 0.3:
        warnings.append("Mało niezależnych nachyleń lub wąski zakres sin²ψ; zwiększ serię pomiarową.")
    design = np.column_stack((np.ones(len(x)), x))
    scale = u if weighted else np.ones(len(x))
    aw, yw = design / scale[:, None], y / scale
    coef, _, rank, _ = np.linalg.lstsq(aw, yw, rcond=None)
    if rank < 2:
        return _result(REJECTED, method, ["Osobliwa regresja sin²ψ."], **empty)
    fitted = design @ coef
    residual = y - fitted
    dof = len(y) - 2
    redchi = float(np.sum((residual / scale) ** 2) / dof)
    cov = np.linalg.inv(aw.T @ aw)
    # Absolute experimental uncertainties never shrink for an unusually good fit.
    cov *= max(1.0, redchi) if weighted else redchi
    slope = float(coef[1])
    stress = slope / (q * d0)
    fit_var = float(cov[1, 1]) / (q * d0) ** 2
    external_var = (stress * uq / q) ** 2 + (stress * ud0 / d0) ** 2
    us = float(np.sqrt(max(0, fit_var + external_var)))
    # Conservative t coverage for regression plus independent normal XEC/d0 terms.
    ci_half = float(np.sqrt((t.ppf(.975, dof) ** 2) * fit_var + norm.ppf(.975) ** 2 * external_var))
    h = np.einsum("ij,jk,ik->i", aw, np.linalg.inv(aw.T @ aw), aw)
    diag_scale = u if weighted else np.full(len(y), np.sqrt(max(redchi, np.finfo(float).eps)))
    standardized = residual / (diag_scale * np.sqrt(np.maximum(1 - h, 1e-12)))
    split = []
    for p in sorted(set(np.round(np.abs(psi), 6)) - {0.0}):
        ip, im = np.isclose(psi, p, atol=1e-6), np.isclose(psi, -p, atol=1e-6)
        if np.any(ip) and np.any(im):
            if weighted:
                wp, wm = 1/u[ip]**2, 1/u[im]**2
                delta = float(np.average(y[ip], weights=wp) - np.average(y[im], weights=wm))
                uncertainty = float(np.sqrt(1/wp.sum() + 1/wm.sum()))
            else:
                delta = float(y[ip].mean()-y[im].mean())
                uncertainty = float(np.sqrt(max(redchi, 0)*(1/ip.sum()+1/im.sum())))
            z = delta / uncertainty if uncertainty > 0 else (0.0 if abs(delta) < 1e-12 else math.inf)
            split.append(dict(abs_psi_deg=float(p), delta_d_A=delta,
                              u_delta_d_A=uncertainty, z_score=float(z) if np.isfinite(z) else None,
                              significant=bool(abs(z) > 3)))
    rejected = False
    if not split:
        warnings.append("Brak par +ψ/-ψ: nie można ocenić rozszczepienia ψ.")
    if any(s["significant"] for s in split):
        warnings.append("Wykryto istotne rozszczepienie +ψ/-ψ (>3u); prosty model płaskiego stanu odrzucono.")
        rejected = True
    p_fit = float(chi2.sf(redchi*dof, dof)) if weighted else None
    if p_fit is not None and p_fit < .001:
        warnings.append("Liniowy model nie opisuje danych w granicach niepewności (p<0.001).")
        rejected = True
    if np.max(np.abs((y-d0)/d0)) > .02:
        warnings.append("Odkształcenie przekracza 2%; sprawdź d0, hkl i liniowy model sprężysty.")
        rejected = True
    loo = []
    for val in np.unique(np.round(x, 10)):
        keep = np.round(x, 10) != val
        if len(np.unique(np.round(x[keep], 10))) >= 2:
            c = np.linalg.lstsq(aw[keep], yw[keep], rcond=None)[0]
            loo.append(float(c[1] / (q*d0)))
    if loo and max(abs(v-stress) for v in loo) > max(3*us, abs(stress)*.2, 1.0):
        warnings.append("Duża wrażliwość naprężenia na pominięcie jednego niezależnego ψ.")
    den = float(np.sum((y - np.mean(y))**2))
    status = REJECTED if rejected else (CONDITIONAL if warnings else VALIDATED)
    return _result(status, method, warnings,
        stress_MPa=None if rejected else float(stress),
        u_stress_MPa=None if rejected else us,
        stress_ci95_MPa=None if rejected else [float(stress-ci_half), float(stress+ci_half)],
        trial_stress_MPa=float(stress),
        independent_psi_count=n_unique, n_measurements=len(y),
        phase=phases[0], hkl=list(hkls[0]), phi_deg=phis[0],
        psi_deg=psis, sin2psi=x.tolist(), d_A=y.tolist(), u_d_A=u.tolist(),
        fitted_d_A=fitted.tolist(), residuals_A=residual.tolist(),
        slope_A=slope, intercept_A=float(coef[0]), covariance_A2=cov.tolist(),
        weighted=weighted, reduced_chi_squared=redchi if weighted else None,
        lack_of_fit_p=p_fit, rmse_A=float(np.sqrt(np.mean(residual**2))),
        r_squared=float(1 - np.sum(residual**2)/den) if den > 0 else None,
        leverage=h.tolist(), standardized_residuals=standardized.tolist(),
        psi_splitting=split, leave_one_psi_out_stress_MPa=loo,
        xec_half_s2_per_MPa=q, u_half_s2_per_MPa=uq, xec_s1_per_MPa=xec_s1_per_MPa,
        d0_A=d0, u_d0_A=ud0, xec_source=xec_source,
        uncertainty_assumptions="Independent d errors, XEC and d0; common wavelength/systematic errors require covariance-aware analysis.")


def wilkens_f_star(eta):
    """Exact Wilkens f*(eta), ITC H eqs. 3.6.35–37; eta > 0.

    Small-eta series is used only below 1e-3 to avoid cancellation, where
    omitted terms are below double precision. This is not a full-range fit.
    """
    original = np.asarray(eta, dtype=float)
    if np.any(~np.isfinite(original)) or np.any(original <= 0):
        raise ValueError("eta musi być skończone i dodatnie.")
    z = original.reshape(-1)
    out = np.empty_like(z)
    lo, hi = z < 1e-3, z >= 1
    a = z[lo]
    out[lo] = -np.log(a)+1.75-np.log(2)+a*a/6-32*a**3/(225*np.pi)
    a = z[hi]
    out[hi] = 256/(45*np.pi*a)-(11/24+np.log(2*a)/4)/a**2
    mid = ~(lo | hi)
    a = z[mid]
    asin = np.arcsin(a)
    root = np.sqrt(1-a*a)
    # scipy.spence(z) = Li2(1-z).
    integral = np.log(2*a)*asin + .5*np.imag(spence(2*a*a-2j*a*root))
    out[mid] = (-np.log(a)+1.75-np.log(2)+256/(45*np.pi*a)
        +2/np.pi*(1-1/(4*a*a))*integral
        -(769/(180*a)+41*a/90+a**3/45)*root/np.pi
        -(11/(12*a*a)+3.5+a*a/3)*asin/np.pi+a*a/6)
    return out.reshape(original.shape)


def wilkens_fourier(length_nm, rho_m2, Re_nm, burgers_vector_nm, g_nm_inv, contrast_factor):
    """A_dislocation(L), with physical Re convention f ~ log(Re/L).

    eta = exp(7/4)/2 * |L|/Re, M = Re[m]*sqrt(rho[m^-2]).
    """
    L = np.abs(np.asarray(length_nm, dtype=float))
    if np.any(~np.isfinite(L)):
        raise ValueError("Nieprawidłowe długości Fouriera.")
    rho = _number(rho_m2, "rho", nonnegative=True) * 1e-18
    Re = _number(Re_nm, "Re", positive=True)
    b = _number(burgers_vector_nm, "b", positive=True)
    g = _number(g_nm_inv, "g", positive=True)
    c = _number(contrast_factor, "C_hkl", positive=True)
    out = np.ones_like(L)
    mask = L > 0
    f = wilkens_f_star((np.exp(1.75)/2)*L[mask]/Re)
    out[mask] = np.exp(-.5*np.pi*rho*b*b*g*g*c*L[mask]**2*f)
    return out


def lognormal_sphere_fourier(length_nm, median_nm, sigma_ln):
    """Volume-weighted sphere autocorrelation averaged over number lognormal D."""
    L = np.abs(np.asarray(length_nm, dtype=float))
    if np.any(~np.isfinite(L)):
        raise ValueError("Nieprawidłowe długości Fouriera.")
    m = _number(median_nm, "mediana D", positive=True)
    s = _number(sigma_ln, "sigma_ln", nonnegative=True)
    if s < 1e-6:
        z = L/m
        return np.where(z < 1, 1-1.5*z+.5*z**3, 0.0)
    def tail(k):
        arg = (np.log(np.maximum(L, np.finfo(float).tiny)/m)-k*s*s)/(np.sqrt(2)*s)
        return m**k*np.exp(.5*k*k*s*s)*.5*erfc(arg)
    a = (tail(3)-1.5*L*tail(2)+.5*L**3*tail(0))/(m**3*np.exp(4.5*s*s))
    return np.clip(a, 0, 1)


def _wp_setup(profiles, wavelength_A, kalpha2_wavelength_A, kalpha2_ratio,
              fourier_step_nm, fourier_points):
    wave = _number(wavelength_A, "lambda1", positive=True)*.1
    r = _number(kalpha2_ratio, "Kalpha2/Kalpha1", nonnegative=True)
    wave2 = _number(kalpha2_wavelength_A, "lambda2", positive=True)*.1 if r else wave
    dl = _number(fourier_step_nm, "krok Fouriera", positive=True)
    n = int(fourier_points)
    if n < 1025 or n > 131073:
        raise ValueError("Liczba punktów Fouriera musi mieścić się w [1025,131073].")
    L = np.arange(n)*dl
    qgrid = np.arange(n)/(2*(n-1)*dl)
    prepared = []
    for p in profiles:
        x = np.asarray(p["x_deg"], float)
        y = np.asarray(p["y"], float)
        if x.ndim != 1 or len(x) < 31 or x.shape != y.shape or not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
            raise ValueError("Każdy pełny profil wymaga >=31 skończonych par x,y.")
        if np.ptp(y) <= np.finfo(float).eps * max(1.0, float(np.max(np.abs(y)))):
            raise ValueError("Profil stały nie zawiera informacji o poszerzeniu refleksu.")
        if not np.all(np.diff(x) > 0) or x[0] <= 0 or x[-1] >= 180:
            raise ValueError("Osie 2theta muszą rosnąć i należeć do (0,180).")
        center = _number(p["center_deg"], "pozycja Kalpha1")
        if not x[0] < center < x[-1]:
            raise ValueError("Centrum refleksu musi leżeć wewnątrz całego profilu.")
        if p.get("overlap", False):
            raise ValueError("Ten model wymaga rozdzielonych refleksów; nakładanie zgłoszone w profilu.")
        g = 2*np.sin(np.deg2rad(center/2))/wave
        q = 2*np.sin(np.deg2rad(x/2))/wave-g
        q2 = 2*np.sin(np.deg2rad(x/2))/wave2-g
        if max(np.max(np.abs(q)), np.max(np.abs(q2))) >= qgrid[-1]*.9:
            raise ValueError("Za duży zakres q dla kroku Fouriera; zmniejsz fourier_step_nm.")
        hg = _number(p["instrument_gaussian_fwhm_deg"], "IPF Gauss", nonnegative=True)
        hl = _number(p["instrument_lorentzian_fwhm_deg"], "IPF Lorentz", nonnegative=True)
        conv = np.cos(np.deg2rad(center/2))*np.pi/(180*wave)
        ipf = np.exp(-2*np.pi**2*(hg*conv/np.sqrt(8*np.log(2)))**2*L**2-np.pi*hl*conv*L)
        c = _number(p["contrast_factor"], "czynnik kontrastu", positive=True)
        if c > 1:
            raise ValueError("C_hkl > 1 wymaga odrębnej weryfikacji konwencji.")
        u = np.asarray(p.get("sigma_y", np.ones_like(y)), float)
        if u.shape != y.shape or np.any(~np.isfinite(u)) or np.any(u <= 0):
            raise ValueError("sigma_y musi zawierać dodatnią niepewność każdego punktu.")
        jac1 = np.cos(np.deg2rad(x/2))*np.pi/(180*wave)
        jac2 = np.cos(np.deg2rad(x/2))*np.pi/(180*wave2)
        prepared.append(dict(x=x,y=y,u=u,q=q,q2=q2,g=g,c=c,ipf=ipf,
            jac1=jac1,jac2=jac2,center=center,hkl=_hkl(p["hkl"]), sigma_known="sigma_y" in p))
    return L,qgrid,prepared,r


def _wp_shape(prep, L, qgrid, params, b, ratio):
    rho, re, med, sd = params
    a = (lognormal_sphere_fourier(L,med,sd)*prep["ipf"]*
         wilkens_fourier(L,rho,re,b,prep["g"],prep["c"]))
    n = 2*(len(L)-1)
    spectrum = np.fft.irfft(a,n=n)[:len(L)]*n*(L[1]-L[0])
    shape = np.interp(np.abs(prep["q"]),qgrid,spectrum)*prep["jac1"]
    if ratio:
        shape += ratio*np.interp(np.abs(prep["q2"]),qgrid,spectrum)*prep["jac2"]
    return shape, float(a[-1])


def simulate_wilkens_profile(profile: dict, *, rho_m2, Re_nm, median_nm,
                sigma_ln, burgers_vector_nm, wavelength_A=1.540598,
                kalpha2_wavelength_A=1.544426, kalpha2_ratio=0.0,
                area=1.0, background=0.0, background_slope=0.0,
                fourier_step_nm=.25, fourier_points=8193):
    """Forward model for verification/planning; returns intensity, no data mutation."""
    p = dict(profile)
    # _wp_setup rejects flat measured profiles; use a non-flat placeholder only
    # in this forward-model entry point, where intensity is not an observation.
    p["y"] = np.arange(len(p["x_deg"]), dtype=float)
    L,qgrid,prep,r = _wp_setup([p],wavelength_A,kalpha2_wavelength_A,kalpha2_ratio,
                               fourier_step_nm,fourier_points)
    shape,_ = _wp_shape(prep[0],L,qgrid,(rho_m2,Re_nm,median_nm,sigma_ln),burgers_vector_nm,r)
    return area*shape+background+background_slope*(prep[0]["x"]-prep[0]["center"])


def fit_wilkens_profiles(profiles: Sequence[dict], *, burgers_vector_nm=None,
                phase="", structure="", wavelength_A=1.540598,
                kalpha2_wavelength_A=1.544426, kalpha2_ratio=0.0,
                instrument_source="", contrast_source="", slip_system="",
                initial=None, multistarts=4, fourier_step_nm=.25,
                fourier_points=8193, max_nfev=500):
    """Experimental constrained whole-profile fit, NOT the external CMWP engine.

    Shared rho, Re, number-lognormal median D and sigma_ln; fixed peak centers,
    fixed per-reflection contrast/IPF; variable-projected area and linear
    backgrounds. Supports isolated reflections and an explicit Kalpha doublet.
    Successful fits always remain CONDITIONAL until external validation.
    """
    method="Wilkens_constrained_whole_profile_experimental"
    empty=dict(rho_m2=None,rho_ci95_m2=None)
    missing=[]
    if len(profiles)<3:
        missing.append("Wymagane co najmniej trzy całe profile różnych refleksów jednej fazy.")
    if burgers_vector_nm is None:
        missing.append("Brak długości wektora Burgersa b [nm].")
    if not phase or not structure:
        missing.append("Brak fazy lub typu struktury.")
    if any("sigma_y" in p for p in profiles) and not all("sigma_y" in p for p in profiles):
        missing.append("Podaj sigma_y dla wszystkich profili albo dla żadnego.")
    for i,p in enumerate(profiles):
        absent=[k for k in ("x_deg","y","hkl","center_deg","contrast_factor",
            "instrument_gaussian_fwhm_deg","instrument_lorentzian_fwhm_deg") if k not in p]
        if absent:
            missing.append(f"Profil {i+1}: brak {', '.join(absent)}.")
        if p.get("phase",phase)!=phase:
            missing.append(f"Profil {i+1}: niezgodna faza.")
    if missing:
        return _result(INSUFFICIENT_DATA,method,missing,**empty)
    try:
        b=_number(burgers_vector_nm,"b",positive=True)
        L,qgrid,prep,ratio=_wp_setup(profiles,wavelength_A,kalpha2_wavelength_A,kalpha2_ratio,
                                    fourier_step_nm,fourier_points)
        if len(set(p["hkl"] for p in prep)) < 3:
            raise ValueError("Wymagane co najmniej trzy różne hkl; powtórzenia nie wystarczą.")
        if all(p["ipf"][-1] == 1 for p in prep):
            raise ValueError("Zerowa IPF we wszystkich profilach wymaga niezależnej kalibracji; brak korekcji.")
        if not 1 <= int(multistarts) <= 16:
            raise ValueError("multistarts musi mieścić się w [1,16].")
    except (ValueError,KeyError,TypeError,OverflowError) as e:
        return _result(REJECTED,method,[str(e)],**empty)
    starts=[(1e14,100,50,.25),(1e15,30,100,.4),(1e13,500,25,.15),(5e14,200,80,.55)]
    if initial:
        try:
            starts[0]=tuple(_number(initial[k],k,positive=True) for k in ("rho_m2","Re_nm","median_nm","sigma_ln"))
        except (KeyError,ValueError,TypeError) as e:
            return _result(REJECTED,method,[f"Nieprawidłowy punkt startowy: {e}"],**empty)
    bounds=np.log(np.array([[1e10,1,2,.025],[1e17,3000,1000,1.0]]))
    weighted=all(p["sigma_known"] for p in prep)
    def prediction(params):
        predicted,nuisance,tails=[],[],[]
        for p in prep:
            shape,tail=_wp_shape(p,L,qgrid,params,b,ratio)
            A=np.column_stack((shape,np.ones(len(shape)),p["x"]-p["center"]))
            coeff=np.linalg.lstsq(A/p["u"][:,None],p["y"]/p["u"],rcond=None)[0]
            # Positive peak area is a physical constraint. If active, refit background only.
            if coeff[0]<0:
                coeff[0]=0
                coeff[1:]=np.linalg.lstsq(A[:,1:]/p["u"][:,None],p["y"]/p["u"],rcond=None)[0]
            predicted.append(A@coeff)
            nuisance.append(coeff)
            tails.append(tail)
        return predicted,nuisance,tails
    def residual(logparams):
        predicted,_,_=prediction(np.exp(logparams))
        return np.concatenate([(v-p["y"])/p["u"] for v,p in zip(predicted,prep)])
    fits=[]
    try:
        for i in range(int(multistarts)):
            s=np.log(starts[i%len(starts)])
            s=np.clip(s,bounds[0]+1e-6,bounds[1]-1e-6)
            fit=least_squares(residual,s,bounds=bounds,method="trf",max_nfev=int(max_nfev),
                              ftol=1e-9,xtol=1e-9,gtol=1e-8)
            fits.append(fit)
        best=min(fits,key=lambda f:float(np.sum(f.fun**2)))
    except (ValueError,np.linalg.LinAlgError,FloatingPointError) as e:
        return _result(REJECTED,method,[f"Błąd dopasowania: {e}"],**empty)
    params=np.exp(best.x)
    predicted,nuisance,tails=prediction(params)
    dof=sum(len(p["x"]) for p in prep)-4-3*len(prep)
    redchi=float(np.sum(best.fun**2)/dof)
    sv=np.linalg.svd(best.jac,compute_uv=False)
    condition=float(sv[0]/sv[-1]) if sv[-1]>0 else math.inf
    cov=np.linalg.pinv(best.jac.T@best.jac)*(max(1.,redchi) if weighted else redchi)
    se=np.sqrt(np.maximum(np.diag(cov),0))
    corr=cov/np.maximum(np.outer(se,se),np.finfo(float).tiny)
    with np.errstate(over="ignore", invalid="ignore"):
        ci=np.exp(np.column_stack((best.x-t.ppf(.975,dof)*se,best.x+t.ppf(.975,dof)*se)))
    warnings=["Eksperymentalny model ograniczony Wilkensa; nie jest pełną implementacją programu CMWP.",
              "Założono sferyczne domeny lognormalne, znane centra/kontrasty/IPF i niezależne rozdzielone refleksy.",
              "Przed interpretacją bezwzględnej rho wymagana niezależna walidacja standardem/CMWP oraz analiza wrażliwości tła, IPF i zakresu."]
    rejected=False
    if not weighted:
        warnings.append("Brak kompletnych sigma_y: statystyka opiera się na residuach; CI nie obejmuje błędów systematycznych.")
    if any(np.any(p["y"] < 0) for p in prep):
        warnings.append("Ujemne intensywności: sprawdź korekcję tła i model niepewności; nie są to surowe zliczenia Poissona.")
    if not instrument_source or not contrast_source or not slip_system:
        warnings.append("Niepełna identyfikacja wzorca IPF, źródła kontrastów lub układu poślizgu.")
    if not best.success:
        rejected=True
        warnings.append("Optymalizator nie uzyskał zbieżności.")
    if condition>1e7 or np.any(se>1.5) or np.any(~np.isfinite(ci)):
        rejected=True
        warnings.append("Parametry nieidentyfikowalne lub niepewności zbyt duże; nie raportuje się rho.")
    if np.any(best.x-bounds[0]<.01) or np.any(bounds[1]-best.x<.01):
        rejected=True
        warnings.append("Wynik osiągnął granicę parametrów; nie raportuje się rho.")
    if max(tails)>1e-7:
        rejected=True
        warnings.append("Niedostateczny zakres długości Fouriera; zwiększ fourier_points.")
    # Independently refine both Fourier extent and spacing at the fitted point.
    # This detects truncation/interpolation error, not model misspecification.
    refined_n = 4 * (int(fourier_points) - 1) + 1
    numerical_change = None
    if refined_n <= 131073:
        fine_L, fine_q, fine_prep, _ = _wp_setup(profiles, wavelength_A,
            kalpha2_wavelength_A, kalpha2_ratio, fourier_step_nm / 2, refined_n)
        changes = []
        for p, fp in zip(prep, fine_prep):
            coarse, _ = _wp_shape(p, L, qgrid, params, b, ratio)
            fine, _ = _wp_shape(fp, fine_L, fine_q, params, b, ratio)
            changes.append(float(np.max(np.abs(coarse - fine)) / max(np.max(fine), 1e-30)))
        numerical_change = max(changes)
        if numerical_change > .002:
            rejected = True
            warnings.append("Niewystarczająca zbieżność całkowania Fouriera (>0.2% wysokości); zagęść siatkę.")
    else:
        warnings.append("Pomijano podwojone sprawdzenie siatki Fouriera ze względu na limit pamięci.")
    # Guard against insufficient angular resolution and Fourier interpolation.
    oversampling=[]
    for p in prep:
        shape,_=_wp_shape(p,L,qgrid,params,b,ratio)
        half=np.flatnonzero(shape>.5*np.max(shape))
        oversampling.append(int(len(half)))
    if min(oversampling)<6:
        rejected=True
        warnings.append("Za mało punktów w szerokości linii; wymagany drobniejszy krok pomiarowy.")
    if weighted and chi2.sf(redchi*dof,dof)<.001:
        rejected=True
        warnings.append("Model nie opisuje intensywności w granicach sigma_y (p<0.001).")
    trial=[dict(success=bool(f.success),chi_squared=float(np.sum(f.fun**2)),
                rho_m2=float(np.exp(f.x[0])),Re_nm=float(np.exp(f.x[1]))) for f in fits]
    comparable=[v for v in trial if v["chi_squared"]<=float(np.sum(best.fun**2))*1.01+1e-8]
    if len(comparable)>1 and max(v["rho_m2"] for v in comparable)/min(v["rho_m2"] for v in comparable)>1.5:
        rejected=True
        warnings.append("Różne starty dają podobny błąd i różne rho; niejednoznaczność rozwiązania.")
    maxcorr=float(np.max(np.abs(corr-np.eye(4))))
    if maxcorr>.98:
        warnings.append("Silna korelacja parametrów (>0.98); CI lokalne wymagają bootstrap/profile likelihood.")
    diagnostics=[]
    for p,v,c in zip(prep,predicted,nuisance):
        diagnostics.append(dict(hkl=list(p["hkl"]),x_deg=p["x"].tolist(),
            observed_y=p["y"].tolist(),fitted_y=v.tolist(),residuals_y=(p["y"]-v).tolist(),
            area=float(c[0]),background=float(c[1]),background_slope=float(c[2])))
    return _result(REJECTED if rejected else CONDITIONAL,method,warnings,
        rho_m2=None if rejected else float(params[0]),
        rho_ci95_m2=None if rejected else ci[0].tolist(),trial_rho_m2=float(params[0]),
        Re_nm=float(params[1]),median_nm=float(params[2]),sigma_ln=float(params[3]),
        M=float(params[1]*1e-9*np.sqrt(params[0])),
        cutoff_convention="eta=exp(7/4)*L/(2*Re), f~ln(Re/L), M=Re*sqrt(rho)",
        parameter_order=["ln(rho_m2)","ln(Re_nm)","ln(median_nm)","ln(sigma_ln)"],
        parameter_ci95=ci.tolist() if np.all(np.isfinite(ci)) else None,covariance_log_parameters=cov.tolist(),
        parameter_correlation=corr.tolist(),jacobian_condition=condition if np.isfinite(condition) else None,
        reduced_chi_squared=redchi if weighted else None,weighted=weighted,
        residual_sum_squares=float(sum(np.sum((p["y"]-v)**2) for p,v in zip(prep,predicted))),
        multistart_results=trial,profiles=diagnostics,phase=phase,structure=structure,
        burgers_vector_nm=b,instrument_source=instrument_source,contrast_source=contrast_source,
        slip_system=slip_system,numerical_fourier_tail=max(tails),
        numerical_profile_relative_change=numerical_change,
        uncertainty_scope="Local covariance, fixed b/contrast/IPF/centers; excludes model/systematic uncertainties.")


def export_whole_profile_request(path, profiles, **settings):
    """Portable transparent input contract, not an official CMWP file format."""
    payload=dict(schema="xrd-analyzer-whole-profile/1",engine="experimental_wilkens",
                 profiles=profiles,settings=settings)
    def plain(v):
        if isinstance(v,np.ndarray): return v.tolist()
        if isinstance(v,np.generic): return v.item()
        raise TypeError(f"Nieobsługiwany typ {type(v).__name__}")
    Path(path).write_text(json.dumps(payload,ensure_ascii=False,indent=2,default=plain,allow_nan=False),encoding="utf-8")
    return str(path)
