"""Numerical line-profile analysis; no GUI, smoothing or implicit calibration.

Angles and component FWHM are degrees 2theta, wavelengths/d are angstrom,
domain sizes are nm. Public results are dataclasses with JSON-safe to_dict().
All material interpretations remain CONDITIONAL pending experimental validation.
See SCIENTIFIC_PROFILE.md for equations, assumptions and uncertainty limits.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
import math
import re
from typing import Optional

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import norm, t as student_t

CU_KALPHA1_A = 1.540598
CU_KALPHA2_A = 1.544426


def _json_safe(value):
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (tuple, list, np.ndarray)):
        return [_json_safe(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


class SerializableResult:
    def to_dict(self):
        return _json_safe(asdict(self))


def _vector(value, name, n=None, minimum=1):
    result = np.asarray(value, dtype=float)
    if result.ndim != 1 or result.size < minimum or not np.all(np.isfinite(result)):
        raise ValueError(f"{name}: wymagana skończona tablica 1D (min. {minimum} punktów).")
    if n is not None and result.size != n:
        raise ValueError(f"{name}: niezgodna liczba punktów.")
    return result


def _positive(value, name, allow_zero=False):
    value = float(value)
    if not math.isfinite(value) or (value < 0 if allow_zero else value <= 0):
        raise ValueError(f"{name}: wymagana wartość {'nieujemna' if allow_zero else 'dodatnia'}.")
    return value


def _angles(value, name="2θ", n=None, minimum=1):
    a = _vector(value, name, n, minimum)
    if np.any((a <= 0) | (a >= 180)):
        raise ValueError(f"{name}: kąty muszą leżeć w (0, 180) stopni.")
    return a


def _sigma(value, n, name):
    if value is None:
        return None
    v = np.full(n, float(value)) if np.ndim(value) == 0 else _vector(value, name, n)
    if not np.all(np.isfinite(v)) or np.any(v <= 0):
        raise ValueError(f"{name}: niepewności muszą być dodatnie i skończone.")
    return v


def _ci(value, uncertainty, factor=1.959963984540054):
    return None if uncertainty is None else (value-factor*uncertainty, value+factor*uncertainty)


def _r2(y, residuals, weights=None):
    w = np.ones(len(y)) if weights is None else weights
    total = float(np.sum(w*(y-np.average(y, weights=w))**2))
    return None if total <= np.finfo(float).eps else 1-float(np.sum(w*residuals**2))/total


def _covariance(jacobian, residuals, absolute_sigma):
    """Scale columns before rank/condition checks, preserving original units."""
    scale = np.linalg.norm(jacobian, axis=0)
    if np.any(scale == 0):
        return None, math.inf
    scaled = jacobian/scale
    _, singular, vt = np.linalg.svd(scaled, full_matrices=False)
    condition = float(singular[0]/singular[-1]) if singular[-1] > 0 else math.inf
    if singular[-1] <= singular[0]*1e-10:
        return None, condition
    cov = (vt.T/singular**2) @ vt
    cov = cov/scale[:, None]/scale[None, :]
    if not absolute_sigma:
        dof = jacobian.shape[0]-jacobian.shape[1]
        if dof <= 0:
            return None, condition
        cov *= float(np.dot(residuals, residuals))/dof
    return cov, condition


def _gradient(fun, point):
    p = np.asarray(point, dtype=float)
    out = np.empty(p.size)
    for i in range(p.size):
        step = max(abs(p[i])*1e-5, 1e-7)
        a, b = p.copy(), p.copy()
        a[i] += step
        if p[i] >= 0 and p[i]-step < 0:
            out[i] = (fun(a)-fun(p))/step
        else:
            b[i] -= step
            out[i] = (fun(a)-fun(b))/(2*step)
    return out


def tch_fwhm_eta(gaussian_fwhm_deg, lorentzian_fwhm_deg):
    """TCH approximation from convolution-component widths (G, L)."""
    g = _positive(gaussian_fwhm_deg, "G", True)
    l = _positive(lorentzian_fwhm_deg, "L", True)
    if g == 0 and l == 0:
        return 0.0, 0.0
    h = (g**5 + 2.69269*g**4*l + 2.42843*g**3*l**2
         + 4.47163*g**2*l**3 + 0.07842*g*l**4 + l**5)**0.2
    q = l/h
    eta = float(np.clip(1.36603*q - 0.47719*q*q + 0.11116*q**3, 0, 1))
    return float(h), eta


def pseudo_voigt(x, center_deg, fwhm_deg, eta):
    """Unit-area pseudo-Voigt; eta is Lorentzian area fraction."""
    h = _positive(fwhm_deg, "FWHM")
    if not math.isfinite(eta) or not 0 <= eta <= 1:
        raise ValueError("eta musi należeć do [0, 1].")
    z = (np.asarray(x, dtype=float)-center_deg)/h
    gaussian = math.sqrt(4*math.log(2)/math.pi)/h*np.exp(-4*math.log(2)*z*z)
    lorentzian = 2/(math.pi*h)/(1+4*z*z)
    return (1-eta)*gaussian + eta*lorentzian


def kalpha2_position(two_theta_kalpha1_deg, wavelength_A=CU_KALPHA1_A,
                     wavelength_kalpha2_A=CU_KALPHA2_A):
    a = float(_angles([two_theta_kalpha1_deg])[0])
    wave1 = _positive(wavelength_A, "λ1")
    wave2 = _positive(wavelength_kalpha2_A, "λ2")
    arg = wave2/wave1*math.sin(math.radians(a/2))
    if arg >= 1:
        raise ValueError("Kα2 nie spełnia warunku Bragga dla tego refleksu.")
    return math.degrees(2*math.asin(arg))


def bragg_d(two_theta_deg, wavelength_A=CU_KALPHA1_A, u_two_theta_deg=None,
            u_wavelength_A=0.0):
    """Return (d_A, u_d_A); independent wavelength/position uncertainty."""
    theta = math.radians(float(_angles([two_theta_deg])[0])/2)
    wave = _positive(wavelength_A, "λ")
    uw = _positive(u_wavelength_A, "u(λ)", True)
    d = wave/(2*math.sin(theta))
    if u_two_theta_deg is None:
        return d, None
    ut = math.radians(_positive(u_two_theta_deg, "u(2θ)", True))
    ud = math.hypot(uw/(2*math.sin(theta)), wave*math.cos(theta)/(4*math.sin(theta)**2)*ut)
    return d, ud


def coupled_profile(x, two_theta_kalpha1_deg, amplitude, gaussian_fwhm_deg,
                    lorentzian_fwhm_deg, background_0=0.0, background_1=0.0,
                    background_reference_deg=None, wavelength_A=CU_KALPHA1_A,
                    wavelength_kalpha2_A=CU_KALPHA2_A, kalpha2_ratio=0.5):
    """Unit-area components; amplitude is integrated Kα1 intensity (intensity·deg)."""
    a = _positive(amplitude, "Pole Kα1", True)
    ratio = _positive(kalpha2_ratio, "Kα2/Kα1", True)
    h, eta = tch_fwhm_eta(gaussian_fwhm_deg, lorentzian_fwhm_deg)
    x = np.asarray(x, dtype=float)
    ref = two_theta_kalpha1_deg if background_reference_deg is None else background_reference_deg
    y = background_0+background_1*(x-ref)+a*pseudo_voigt(x, two_theta_kalpha1_deg, h, eta)
    if ratio > 0:
        center2 = kalpha2_position(two_theta_kalpha1_deg, wavelength_A, wavelength_kalpha2_A)
        y = y+ratio*a*pseudo_voigt(x, center2, h, eta)
    return y


@dataclass
class PeakFitResult(SerializableResult):
    hkl: str
    two_theta_kalpha1_deg: float
    u_two_theta_deg: Optional[float]
    two_theta_ci95_deg: Optional[tuple]
    d_A: float
    u_d_A: Optional[float]
    fwhm_obs_deg: float
    u_fwhm_obs_deg: Optional[float]
    fwhm_ci95_deg: Optional[tuple]
    gaussian_fwhm_deg: float
    lorentzian_fwhm_deg: float
    eta: float
    amplitude: float
    background_0: float
    background_1: float
    background_reference_deg: float
    wavelength_A: float
    wavelength_kalpha2_A: float
    kalpha2_ratio: float
    x_deg: list
    fitted_y: list
    residuals: list
    covariance: Optional[list]
    parameter_names: list
    rmse: float
    r_squared: Optional[float]
    reduced_chi_squared: Optional[float]
    jacobian_condition: float
    n_points: int
    fit_success: bool
    status: str
    warnings: list
    width_definition: str = "TCH single Kalpha1-component FWHM, degrees 2theta"
    uncertainty_method: str = "linearized covariance"

    @property
    def warning(self):
        return " ".join(self.warnings)


def fit_peak_profile(x, y, center_hint=None, *, wavelength_A=CU_KALPHA1_A,
                     wavelength_kalpha2_A=CU_KALPHA2_A, kalpha2_ratio=0.5,
                     sigma_y=None, hkl="", u_wavelength_A=0.0):
    """Fit one isolated reflection/doublet in an explicitly selected raw-data ROI.

    Six parameters: center, Kα1 area, G, L, linear-background intercept/slope.
    sigma_y supplied: absolute known 1σ weights; absent: residual-estimated OLS.
    No robust loss, smoothing, source-ratio refinement or automatic overlap removal.
    """
    x = _angles(x, minimum=12)
    y = _vector(y, "intensywność", len(x))
    if np.any(np.diff(x) <= 0):
        raise ValueError("Oś 2θ musi być ściśle rosnąca, bez duplikatów.")
    wave = _positive(wavelength_A, "λ1")
    wave2 = _positive(wavelength_kalpha2_A, "λ2")
    ratio = _positive(kalpha2_ratio, "Kα2/Kα1", True)
    if ratio > 1:
        raise ValueError("Kα2/Kα1 musi należeć do [0, 1]; sprawdź oznaczenie linii.")
    sigma = _sigma(sigma_y, len(x), "u(I)")
    divisor = np.ones(len(x)) if sigma is None else sigma
    ref, span = float(np.mean(x)), float(np.ptp(x))
    edge = max(2, len(x)//10)
    edge_idx = np.r_[0:edge, len(x)-edge:len(x)]
    b1, b0 = np.polyfit(x[edge_idx]-ref, y[edge_idx], 1)
    signal = np.maximum(y-(b0+b1*(x-ref)), 0)
    if np.max(signal) <= np.finfo(float).eps*max(float(np.max(np.abs(y))), 1):
        raise ValueError("Brak dodatniego refleksu ponad tłem w wybranym zakresie.")
    center = float(x[np.argmax(signal)]) if center_hint is None else float(center_hint)
    upper_center = float(x[-1])
    if ratio > 0 and wave2 > wave:
        upper_center = min(upper_center, math.degrees(2*math.asin(wave/wave2))-1e-7)
    if not x[0] < center < upper_center:
        raise ValueError("Przybliżona pozycja refleksu musi leżeć wewnątrz zakresu.")
    above = x[signal >= 0.5*np.max(signal)]
    width = max(float(np.ptp(above)) if len(above) > 1 else span/10, span/len(x)*2)
    width = min(width, span/3)
    area = max(float(np.trapezoid(signal, x))/(1+ratio), np.max(signal)*width/2)
    tiny = max(span*1e-7, 1e-9)
    bounds = ([float(x[0]), 0, tiny, 0, -np.inf, -np.inf],
              [upper_center, np.inf, span, span, np.inf, np.inf])

    def evaluate(p):
        return coupled_profile(x, *p[:4], p[4], p[5], ref, wave, wave2, ratio)

    def residual(p):
        return (evaluate(p)-y)/divisor

    # Two different shape starts reduce a common G/L local-convergence failure.
    candidates = []
    for gfrac, lfrac in ((0.8, 0.3), (0.25, 0.8)):
        initial = [center, area, max(width*gfrac, tiny*2), width*lfrac, b0, b1]
        candidates.append(least_squares(residual, initial, bounds=bounds, x_scale="jac",
                                       max_nfev=3000, ftol=1e-11, xtol=1e-11, gtol=1e-11))
    fit = min(candidates, key=lambda item: item.cost)
    p = fit.x
    predicted, residuals = evaluate(p), y-evaluate(p)
    cov, condition = _covariance(fit.jac, fit.fun, sigma is not None)
    h, eta = tch_fwhm_eta(p[2], p[3])
    warnings = ["Model lokalny symetryczny: zweryfikuj nakładanie refleksów, tło i asymetrię w residuach."]
    if sigma is None:
        warnings.append("Nie podano u(I); niepewności oszacowano z rozrzutu residuów, χ²red nie jest raportowane.")
    active = np.asarray(fit.active_mask) != 0
    if np.any(active[:4]):
        warnings.append("Parametr profilu osiągnął ograniczenie; symetryczne CI mogą być niewiarygodne.")
    if cov is None or condition > 1e6:
        warnings.append("Słaba identyfikowalność parametrów / źle uwarunkowana macierz Jacobiego.")
    ut = math.sqrt(max(cov[0, 0], 0)) if cov is not None else None
    uh = None
    if cov is not None:
        grad = _gradient(lambda v: tch_fwhm_eta(*v)[0], p[2:4])
        uh = math.sqrt(max(float(grad @ cov[2:4, 2:4] @ grad), 0))
    d, ud = bragg_d(p[0], wave, ut, u_wavelength_A)
    if min(p[0]-x[0], x[-1]-p[0]) < 2*h:
        warnings.append("Zakres obejmuje mniej niż 2 FWHM po jednej stronie Kα1; rozszerz okno i porównaj wynik.")
    if ratio > 0 and kalpha2_position(p[0], wave, wave2)+h > x[-1]:
        warnings.append("Kα2 jest częściowo poza zakresem; rozszerz okno dopasowania.")
    if h < 2*float(np.median(np.diff(x))):
        warnings.append("FWHM jest słabo rozdzielone względem kroku pomiarowego.")
    dof = len(x)-len(p)
    factor = norm.ppf(.975) if sigma is not None else student_t.ppf(.975, dof)
    status = "CONDITIONAL" if fit.success and cov is not None else "REJECTED"
    return PeakFitResult(str(hkl), float(p[0]), ut, _ci(p[0], ut, factor), d, ud, h, uh,
                         _ci(h, uh, factor), float(p[2]), float(p[3]), eta, float(p[1]),
                         float(p[4]), float(p[5]), ref, wave, wave2, ratio, x.tolist(),
                         predicted.tolist(), residuals.tolist(), None if cov is None else cov.tolist(),
                         ["two_theta_kalpha1_deg", "amplitude", "gaussian_fwhm_deg",
                          "lorentzian_fwhm_deg", "background_0", "background_1"],
                         float(np.sqrt(np.mean(residuals**2))), _r2(y, residuals),
                         float(np.dot(fit.fun, fit.fun)/dof) if sigma is not None else None,
                         condition, len(x), bool(fit.success), status, warnings)


def evaluate_peak_profile(x, result):
    return coupled_profile(x, result.two_theta_kalpha1_deg, result.amplitude,
                           result.gaussian_fwhm_deg, result.lorentzian_fwhm_deg,
                           result.background_0, result.background_1, result.background_reference_deg,
                           result.wavelength_A, result.wavelength_kalpha2_A, result.kalpha2_ratio)


@dataclass
class InstrumentBroadeningModel(SerializableResult):
    two_theta_deg: list
    gaussian_fwhm_deg: list
    lorentzian_fwhm_deg: list
    u_gaussian_fwhm_deg: Optional[list] = None
    u_lorentzian_fwhm_deg: Optional[list] = None
    name: str = "Profil instrumentalny użytkownika"
    standard: str = ""
    measurement_date: str = ""
    radiation: str = "Cu Kα"
    wavelength_A: float = CU_KALPHA1_A
    optics_id: str = ""
    source_file: str = ""
    model: str = "linear interpolation of TCH Gaussian/Lorentzian components"

    def __post_init__(self):
        x = _angles(self.two_theta_deg, "2θ wzorca", minimum=2)
        if np.any(np.diff(x) <= 0):
            raise ValueError("Kąty wzorca muszą być ściśle rosnące.")
        self.two_theta_deg = x.tolist()
        for key in ("gaussian_fwhm_deg", "lorentzian_fwhm_deg"):
            v = _vector(getattr(self, key), key, len(x))
            if np.any(v < 0):
                raise ValueError("Składowe profilu instrumentalnego nie mogą być ujemne.")
            setattr(self, key, v.tolist())
        if np.any(np.asarray(self.gaussian_fwhm_deg)+self.lorentzian_fwhm_deg <= 0):
            raise ValueError("Wzorzec musi mieć dodatnią szerokość każdego refleksu.")
        for key in ("u_gaussian_fwhm_deg", "u_lorentzian_fwhm_deg"):
            v = getattr(self, key)
            if v is not None:
                v = _vector(v, key, len(x))
                if np.any(v < 0):
                    raise ValueError("Niepewności wzorca nie mogą być ujemne.")
                setattr(self, key, v.tolist())
        _positive(self.wavelength_A, "λ wzorca")

    def at(self, two_theta_deg):
        """Return G, L, uG, uL. Linear interpolation, NEVER extrapolation.

        Independent node errors are combined using squared interpolation weights.
        Calibration systematics/correlations must be assessed separately.
        """
        angle = float(two_theta_deg)
        x = np.asarray(self.two_theta_deg)
        if not math.isfinite(angle) or not x[0] <= angle <= x[-1]:
            raise ValueError("Refleks poza zakresem wzorca; ekstrapolacja jest zabroniona.")
        index = min(max(int(np.searchsorted(x, angle, side="right"))-1, 0), len(x)-2)
        f = (angle-x[index])/(x[index+1]-x[index])
        result = []
        for key in ("gaussian_fwhm_deg", "lorentzian_fwhm_deg"):
            v = getattr(self, key)
            result.append(float((1-f)*v[index]+f*v[index+1]))
        for key in ("u_gaussian_fwhm_deg", "u_lorentzian_fwhm_deg"):
            v = getattr(self, key)
            result.append(None if v is None else float(math.hypot((1-f)*v[index], f*v[index+1])))
        return tuple(result)


@dataclass
class CorrectedPeakResult(SerializableResult):
    hkl: str
    two_theta_deg: float
    fwhm_obs_deg: float
    fwhm_instr_deg: Optional[float]
    fwhm_corrected_deg: Optional[float]
    u_fwhm_corrected_deg: Optional[float]
    fwhm_ci95_deg: Optional[tuple]
    gaussian_corrected_deg: Optional[float]
    lorentzian_corrected_deg: Optional[float]
    status: str
    warnings: list
    uncertainty_complete: bool = False


def correct_instrument_broadening(fit, instrument):
    """G_sample=sqrt(Gobs²-Ginst²), L_sample=Lobs-Linst, followed by TCH.

    Invalid component deconvolution is rejected, never silently clamped.
    Full fit G/L covariance is propagated. Standard-node errors are independent.
    """
    result = CorrectedPeakResult(fit.hkl, fit.two_theta_kalpha1_deg, fit.fwhm_obs_deg,
                                None, None, None, None, None, None, "REJECTED", [])
    if not fit.fit_success or fit.status == "REJECTED":
        result.warnings.append("Nie zaakceptowano dopasowania profilu obserwowanego.")
        return result
    if not math.isclose(fit.wavelength_A, instrument.wavelength_A, rel_tol=1e-6):
        result.warnings.append("Niezgodne długości fali próbki i wzorca instrumentalnego.")
        return result
    try:
        gi, li, ugi, uli = instrument.at(fit.two_theta_kalpha1_deg)
    except ValueError as exc:
        result.status = "INSUFFICIENT_DATA"
        result.warnings.append(str(exc))
        return result
    result.fwhm_instr_deg = tch_fwhm_eta(gi, li)[0]
    go, lo = fit.gaussian_fwhm_deg, fit.lorentzian_fwhm_deg
    g2, ls = go*go-gi*gi, lo-li
    tolerance = 1e-12*max(go*go, gi*gi, 1)
    if g2 < -tolerance or ls < -1e-10:
        result.warnings.append("Składowa obserwowana węższa od instrumentalnej; brak fizycznej dekonwolucji.")
        return result
    gs, ls = math.sqrt(max(g2, 0)), max(ls, 0)
    width = tch_fwhm_eta(gs, ls)[0]
    if width <= 1e-10:
        result.status = "INSUFFICIENT_DATA"
        result.warnings.append("Poszerzenie próbki nie jest rozdzielone od instrumentalnego.")
        return result
    result.gaussian_corrected_deg, result.lorentzian_corrected_deg = gs, ls
    result.fwhm_corrected_deg = width
    result.status = "CONDITIONAL"
    result.warnings.append("Korekcja zakłada zgodną optykę, geometrię, model linii i pomijalne poszerzenie wzorca.")
    if gs > 1e-8 and fit.covariance is not None:
        cov = np.zeros((4, 4))
        cov[:2, :2] = np.asarray(fit.covariance)[2:4, 2:4]
        cov[2, 2], cov[3, 3] = (ugi or 0)**2, (uli or 0)**2
        gh = _gradient(lambda v: tch_fwhm_eta(*v)[0], [gs, ls])
        deriv = np.array([gh[0]*go/gs, gh[1], -gh[0]*gi/gs, -gh[1]])
        uw = math.sqrt(max(float(deriv @ cov @ deriv), 0))
        result.u_fwhm_corrected_deg, result.fwhm_ci95_deg = uw, _ci(width, uw)
        result.uncertainty_complete = ugi is not None and uli is not None
        if width <= 1.96*uw:
            result.warnings.append("Skorygowana szerokość nierozdzielona od zera w 95% CI; rozmiar może być niewyznaczalny.")
    else:
        result.warnings.append("Korekcja blisko granicy rozdzielczości: liniowa propagacja niepewności niedostępna.")
    if not result.uncertainty_complete:
        result.warnings.append("Niepełny budżet niepewności: brak u(G/L) wzorca lub osobliwa korekcja.")
    return result


@dataclass
class ScherrerResult(SerializableResult):
    D_nm: Optional[float]
    u_D_nm: Optional[float]
    D_ci95_nm: Optional[tuple]
    K: float
    beta_deg: float
    status: str
    warnings: list
    interpretation: str = "Warunkowa długość domen koherentnego rozpraszania; nie wielkość ziarna."


def scherrer_size(two_theta_deg, beta_deg, wavelength_A=CU_KALPHA1_A, K=0.9,
                  u_beta_deg=None, u_two_theta_deg=None, u_wavelength_A=0.0, u_K=0.0):
    """Caller MUST provide instrument-corrected single-component FWHM in degrees."""
    theta = math.radians(float(_angles([two_theta_deg])[0])/2)
    beta = _positive(beta_deg, "Skorygowana FWHM")
    wave, k = _positive(wavelength_A, "λ"), _positive(K, "K")
    ub = None if u_beta_deg is None else _positive(u_beta_deg, "u(β)", True)
    ut = None if u_two_theta_deg is None else _positive(u_two_theta_deg, "u(2θ)", True)
    uw, uk = _positive(u_wavelength_A, "u(λ)", True), _positive(u_K, "u(K)", True)
    warnings = ["D jest warunkowe: odkształcenia, defekty, rozkład rozmiarów i anizotropia mogą poszerzać linię."]
    if ub is not None and beta <= 1.96*ub:
        return ScherrerResult(None, None, None, k, beta, "INSUFFICIENT_DATA",
                              warnings+["95% CI szerokości obejmuje zero; D nie jest rozdzielone."])
    size = k*wave*0.1/(math.radians(beta)*math.cos(theta))
    uncertainty = None
    if ub is not None:
        uncertainty = size*math.sqrt((ub/beta)**2 + (uw/wave)**2 + (uk/k)**2 +
                                    (math.tan(theta)*math.radians(ut or 0)/2)**2)
    else:
        warnings.append("Brak u(β); brak niepewności D.")
    if ut is None:
        warnings.append("Nie uwzględniono niepewności położenia refleksu.")
    return ScherrerResult(size, uncertainty, _ci(size, uncertainty), k, beta, "CONDITIONAL", warnings)


@dataclass
class WilliamsonHallResult(SerializableResult):
    epsilon: float
    u_epsilon: float
    epsilon_ci95: tuple
    intercept: float
    u_intercept: float
    intercept_ci95: tuple
    D_nm: Optional[float]
    u_D_nm: Optional[float]
    D_ci95_nm: Optional[tuple]
    r_squared: Optional[float]
    adjusted_r_squared: Optional[float]
    rmse: float
    reduced_chi_squared: Optional[float]
    n_reflections: int
    hkls: list
    x: list
    y: list
    sigma_y: Optional[list]
    residuals: list
    leverage: list
    studentized_residuals: list
    cooks_distance: list
    leave_one_out: list
    covariance: list
    valid_for_quantitative_interpretation: bool
    status: str
    warnings: list
    interpretation: str = "Pozorne mikroodkształcenie (apparent microstrain), model izotropowy FWHM-WH."


def fit_williamson_hall(two_theta_deg, beta_deg, wavelength_A=CU_KALPHA1_A, K=.9,
                        sigma_beta_deg=None, hkls=None):
    """WH y=beta(rad)*cos(theta), x=4*sin(theta); slope=epsilon.

    beta MUST be instrument-corrected FWHM of the same phase. At least three
    distinct reflections. X/position errors and inter-reflection covariances are
    not modelled; supplied beta errors are absolute independent standard errors.
    """
    angles = _angles(two_theta_deg, minimum=3)
    beta = _vector(beta_deg, "β", len(angles))
    if np.any(beta <= 0) or len(np.unique(np.round(angles, 8))) < 3:
        raise ValueError("WH wymaga co najmniej trzech różnych refleksów z dodatnią skorygowaną FWHM.")
    if hkls is not None and len(hkls) != len(angles):
        raise ValueError("Niezgodna liczba hkl.")
    if hkls is not None:
        labelled = [str(v).strip() for v in hkls if str(v).strip()]
        if len(set(labelled)) != len(labelled):
            raise ValueError("WH wymaga niezależnych refleksów; powtórzono oznaczenie hkl.")
    wave, k = _positive(wavelength_A, "λ"), _positive(K, "K")
    sigma_beta = _sigma(sigma_beta_deg, len(angles), "u(β)")
    theta = np.deg2rad(angles/2)
    x, y = 4*np.sin(theta), np.deg2rad(beta)*np.cos(theta)
    sy = None if sigma_beta is None else np.deg2rad(sigma_beta)*np.cos(theta)
    divisor = np.ones(len(x)) if sy is None else sy
    design = np.column_stack([np.ones(len(x)), x])
    a, b = design/divisor[:, None], y/divisor
    params, _, rank, _ = np.linalg.lstsq(a, b, rcond=None)
    if rank != 2:
        raise ValueError("Osobliwa regresja WH.")
    residuals = y-design@params
    cov, condition = _covariance(a, residuals/divisor, sy is not None)
    if cov is None or condition > 1e6:
        raise ValueError("Refleksy mają za mały zakres kątowy do stabilnej regresji WH.")
    intercept, epsilon = map(float, params)
    ui, ue = np.sqrt(np.maximum(np.diag(cov), 0))
    dof = len(x)-2
    factor = float(norm.ppf(.975) if sy is not None else student_t.ppf(.975, dof))
    ici, eci = _ci(intercept, ui, factor), _ci(epsilon, ue, factor)
    warnings = ["Interpretacja warunkowa: ta sama faza/hkl, izotropowe poszerzenie i brak istotnego overlapu.",
                "Regresja uwzględnia u(β), ale pomija niepewności kątów i korelacje wspólnego wzorca."]
    if sy is None:
        warnings.append("Brak u(β): użyto regresji nieważonej i niepewności z residuów.")
    size = usize = size_ci = None
    if intercept <= 0:
        warnings.append("Wyraz wolny WH ≤ 0: D_WH niewyznaczalne fizycznie dla tego dopasowania.")
    elif ici[0] <= 0:
        warnings.append("95% CI wyrazu wolnego obejmuje zero: D_WH nie jest statystycznie rozdzielone.")
    else:
        size = k*wave*.1/intercept
        usize = size*ui/intercept
        size_ci = (k*wave*.1/ici[1], k*wave*.1/ici[0])
    if epsilon < 0:
        warnings.append("Ujemne nachylenie: model nie daje fizycznego dodatniego mikroodkształcenia.")
    elif eci[0] <= 0:
        warnings.append("95% CI ε obejmuje zero; mikroodkształcenie nierozdzielone statystyczcznie.")
    hat = np.sum((a @ np.linalg.inv(a.T@a))*a, axis=1)
    whitened = residuals/divisor
    mse = float(np.dot(whitened, whitened)/dof)
    scale = 1.0 if sy is not None else math.sqrt(max(mse, 0))
    studentized = whitened/np.maximum(scale*np.sqrt(np.maximum(1-hat, 1e-15)), 1e-15)
    cooks = studentized**2*hat/(2*np.maximum(1-hat, 1e-15))
    if np.any(np.abs(studentized) > 3):
        warnings.append("Co najmniej jedno residuum standaryzowane przekracza 3; sprawdź refleksy odstające.")
    if np.any(cooks > 4/len(x)):
        warnings.append("Co najmniej jeden refleks silnie wpływa na regresję (Cook > 4/n).")
    loo = []
    for i in range(len(x)):
        keep = np.arange(len(x)) != i
        pi = np.linalg.lstsq(a[keep], b[keep], rcond=None)[0]
        loo.append({"omitted_index": i, "intercept": float(pi[0]), "epsilon": float(pi[1]),
                    "D_nm": float(k*wave*.1/pi[0]) if pi[0] > 0 else None})
    if any(item["intercept"] <= 0 for item in loo) and intercept > 0:
        warnings.append("Pominięcie jednego refleksu zmienia znak wyrazu wolnego; niestabilny D_WH.")
    r2 = _r2(y, residuals, 1/divisor**2)
    if r2 is not None and r2 < .8:
        warnings.append("Słaba liniowość WH (R² < 0.8); próg diagnostyczny, nie poziom ufności.")
    physical = epsilon >= 0 and intercept > 0 and ici[0] > 0
    # Never label an unvalidated experimental/material model VALIDATED from R² alone.
    return WilliamsonHallResult(epsilon, float(ue), eci, intercept, float(ui), ici, size,
                                None if usize is None else float(usize), size_ci, r2,
                                None if r2 is None else 1-(1-r2)*(len(x)-1)/(len(x)-2),
                                float(np.sqrt(np.mean(residuals**2))), mse if sy is not None else None,
                                len(x), list(map(str, hkls)) if hkls is not None else [""]*len(x),
                                x.tolist(), y.tolist(), None if sy is None else sy.tolist(),
                                residuals.tolist(), hat.tolist(), studentized.tolist(), cooks.tolist(),
                                loo, cov.tolist(), False, "CONDITIONAL" if physical else "REJECTED", warnings)


def _parse_hkl(value):
    if isinstance(value, str):
        text = value.strip().strip("()[]")
        if re.fullmatch(r"\d{3}", text):
            values = [int(v) for v in text]
        else:
            values = [int(v) for v in re.findall(r"[-+]?\d+", text)]
    else:
        values = list(value)
    if len(values) != 3 or any(float(v) != int(v) for v in values):
        raise ValueError("hkl musi zawierać trzy liczby całkowite.")
    if sum(int(v)**2 for v in values) == 0:
        raise ValueError("hkl=(0,0,0) nie jest refleksem.")
    return tuple(int(v) for v in values)


@dataclass
class LatticeFitResult(SerializableResult):
    a_A: float
    u_a_A: Optional[float]
    a_ci95_A: Optional[tuple]
    zero_shift_deg: float
    u_zero_shift_deg: Optional[float]
    rms_2theta_deg: float
    chi2_red: Optional[float]
    residuals_deg: list
    predicted_2theta_deg: list
    used_hkls: list
    covariance: Optional[list]
    a_zero_shift_correlation: Optional[float]
    n_reflections: int
    fit_success: bool
    status: str
    warnings: list


def fit_cubic_lattice(two_theta_deg, hkls, wavelength_A=CU_KALPHA1_A,
                      fit_zero_shift=True, sigma_2theta_deg=None):
    """Measured 2theta = 2 asin(lambda*sqrt(h²+k²+l²)/(2a)) + zero shift."""
    nparam = 2 if fit_zero_shift else 1
    angles = _angles(two_theta_deg, minimum=nparam+1)
    if len(hkls) != len(angles):
        raise ValueError("Niezgodna liczba kątów i hkl.")
    hkl = [_parse_hkl(v) for v in hkls]
    norms = np.sqrt([sum(c*c for c in v) for v in hkl])
    if len(np.unique(norms)) < nparam+1:
        raise ValueError("Za mało niezależnych rodzin refleksów do dopasowania a / zero shift.")
    wave = _positive(wavelength_A, "λ")
    sigma = _sigma(sigma_2theta_deg, len(angles), "u(2θ)")
    divisor = np.ones(len(angles)) if sigma is None else sigma
    a0 = float(np.median(wave*norms/(2*np.sin(np.deg2rad(angles/2)))))
    lower = float(wave*max(norms)/2*(1+1e-10))

    def predict(params):
        return np.rad2deg(2*np.arcsin(wave*norms/(2*params[0])))+(params[1] if fit_zero_shift else 0)

    fit = least_squares(lambda p: (predict(p)-angles)/divisor,
                        [max(a0, lower*1.001), 0] if fit_zero_shift else [max(a0, lower*1.001)],
                        bounds=([lower, -5], [max(100, a0*10), 5]) if fit_zero_shift else ([lower], [max(100, a0*10)]),
                        x_scale="jac", max_nfev=2000, xtol=1e-13, ftol=1e-13, gtol=1e-13)
    predicted = predict(fit.x)
    residuals = angles-predicted
    cov, condition = _covariance(fit.jac, fit.fun, sigma is not None)
    ua = math.sqrt(max(cov[0, 0], 0)) if cov is not None else None
    uz = math.sqrt(max(cov[1, 1], 0)) if cov is not None and fit_zero_shift else None
    correlation = None
    warnings = ["Pozorny parametr sieci; różnice mogą wynikać ze składu, odkształceń i geometrii pomiaru.",
                "Niepewność λ i systematyczne błędy pozycji nie są zawarte w u(a)."]
    if fit_zero_shift and ua and uz:
        correlation = float(cov[0, 1]/(ua*uz))
        if abs(correlation) > .95:
            warnings.append("Silna korelacja a i zero shift; rozszerz zakres kątowy / sprawdź wzorzec.")
    if sigma is None:
        warnings.append("Brak u(2θ); niepewność z residuów, bez χ²red.")
    if np.any(fit.active_mask):
        warnings.append("Parametr osiągnął ograniczenie dopasowania.")
    if cov is None or condition > 1e6:
        warnings.append("Parametry słabo identyfikowalne.")
    dof = len(angles)-nparam
    factor = norm.ppf(.975) if sigma is not None else student_t.ppf(.975, dof)
    return LatticeFitResult(float(fit.x[0]), ua, _ci(fit.x[0], ua, factor),
                            float(fit.x[1]) if fit_zero_shift else 0.0, uz,
                            float(np.sqrt(np.mean(residuals**2))),
                            float(np.dot(fit.fun, fit.fun)/dof) if sigma is not None else None,
                            residuals.tolist(), predicted.tolist(), [str(v) for v in hkl],
                            None if cov is None else cov.tolist(), correlation, len(angles),
                            bool(fit.success), "CONDITIONAL" if fit.success and cov is not None else "REJECTED", warnings)
