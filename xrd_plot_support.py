"""Presentation-only smoothing and defensive restoration of optional UI state.

No function in this module writes to Dataset or performs scientific peak analysis.
"""
from __future__ import annotations

import math
import warnings
from copy import deepcopy

import numpy as np
from matplotlib.colors import is_color_like
from scipy.ndimage import maximum_filter1d, median_filter, uniform_filter1d
from scipy.signal import find_peaks, peak_widths, savgol_coeffs, savgol_filter

# Display smoothing methods. "savgol" is the original plain Savitzky-Golay filter and is kept so
# that figures saved with it look exactly as before; "peak_preserving" is the default for new rows.
SMOOTH_METHODS = ("peak_preserving", "savgol")
DEFAULT_SMOOTH_METHOD = "peak_preserving"
SMOOTH_METHOD_LABELS = {
    "peak_preserving": "Z zachowaniem pików (zalecane)",
    "savgol": "Klasyczny Savitzky–Golay",
}
PEAK_KEEP_FRACTION = 0.99   # share of peak height a display filter must preserve
_FWHM_TO_SIGMA = 2.3548200450309493


def finite_number(value, default, lower=None, upper=None):
    try:
        number = float(value)
        if not math.isfinite(number):
            return default
    except (TypeError, ValueError, OverflowError):
        return default
    if lower is not None:
        number = max(lower, number)
    if upper is not None:
        number = min(upper, number)
    return number


def optional_bool(value, default=False):
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and value in (0, 1):
        return bool(value)
    if isinstance(value, str) and value.lower() in ("true", "false", "0", "1"):
        return value.lower() in ("true", "1")
    return default


def display_smoothing_settings(row):
    """Current flat schema, with migration of the report's nested draft schema."""
    row = row if isinstance(row, dict) else {}
    nested = row.get("display_smoothing", {})
    nested = nested if isinstance(nested, dict) else {}
    method = row.get("smooth_method", nested.get("method", DEFAULT_SMOOTH_METHOD))
    valid_method = method in SMOOTH_METHODS
    return {
        "smooth_enabled": optional_bool(row.get("smooth_enabled", nested.get("enabled", False))) and valid_method,
        "smooth_strength_percent": finite_number(row.get("smooth_strength_percent", row.get("smooth_percent", nested.get("strength_percent", 0.0))), 0.0, 0.0, 100.0),
        "smooth_method": method if valid_method else DEFAULT_SMOOTH_METHOD,
        "smooth_window_deg": finite_number(row.get("smooth_window_deg", nested.get("window_deg", 0.25)), 0.25, 0.001, 5.0),
        "smooth_polyorder": int(finite_number(row.get("smooth_polyorder", nested.get("polyorder", 2)), 2, 1, 4)),
    }


def display_smoothing_window(n, step_deg, window_deg=0.25, polyorder=2):
    """Odd window nearest above round(window_deg / step_deg), at most n.

    The minimum odd window is strictly larger than polyorder + 1. Fewer than
    five points are deliberately left unchanged; 0 means no applicable filter.
    """
    if n < 5 or not math.isfinite(step_deg) or step_deg <= 0:
        return 0
    window = max(polyorder + 2, int(math.floor(window_deg / step_deg + 0.5)))
    window += (window % 2 == 0)
    window = min(window, n if n % 2 else n - 1)
    return window if window > polyorder else 0


def _odd(n):
    n = int(n)
    return n if n % 2 else n + 1


def _profile_retention(width_pts, window, order, eta=0.5):
    """Height kept by SG(window, order) for a pseudo-Voigt (eta = Lorentzian share) of FWHM width_pts."""
    k = np.arange(window) - window // 2
    gauss = np.exp(-0.5 * (k / (width_pts / _FWHM_TO_SIGMA)) ** 2)
    lorentz = 1.0 / (1.0 + (2.0 * k / width_pts) ** 2)
    return float(np.dot(savgol_coeffs(window, order), (1.0 - eta) * gauss + eta * lorentz))


def _widest_safe_window(width_pts, max_window, keep, order=4):
    """Largest odd window whose order-4 SG keeps >= keep of a peak of this width; 0 if none does."""
    best = 0
    for window in range(order + 3, max_window + 1, 2):
        if _profile_retention(width_pts, window, order) < keep:
            break
        best = window
    return best


def _local_noise(values, span):
    """Robust local noise level (std) from second differences; follows Poisson-like growth with counts."""
    n = len(values)
    d2 = np.empty(n)
    d2[1:-1] = values[:-2] - 2.0 * values[1:-1] + values[2:]
    d2[0], d2[-1] = d2[1], d2[-2]
    size = min(_odd(max(51, span)), n if n % 2 else n - 1)
    sigma = 1.4826 * median_filter(np.abs(d2), size=size, mode="nearest") / math.sqrt(6.0)
    return np.maximum(sigma, 1e-9 * max(float(np.ptp(values)), 1e-30))


def peak_preserving_filter(values, step_deg, window_deg=0.25, polyorder=2, keep=PEAK_KEEP_FRACTION):
    """Smooth like SG(window_deg, polyorder), but keep narrow peaks at their measured height.

    The strong SG filter flattens peaks that are only a few samples wide. Peaks are located on a
    lightly smoothed curve; where the strong filter would remove more than (1 - keep) of a peak's
    height, the neighbourhood of that peak is replaced by the widest order-4 SG filter that still
    keeps >= keep of the height (or by the measured values for peaks only ~3 samples wide). The
    transition is weighted by where the two filters actually disagree, so the background keeps the
    full noise reduction and broad peaks are treated exactly like the plain filter.
    Data without measurable noise are returned unchanged. Result stays inside [min(y), max(y)].
    Returns None when fewer than 5 points make filtering inapplicable.
    """
    values = np.asarray(values, dtype=float)
    n = len(values)
    window = _odd(max(polyorder + 2, int(math.floor(window_deg / step_deg + 0.5))))
    window = min(window, n if n % 2 else n - 1)
    if n < 5 or window <= polyorder:
        return None
    strong = savgol_filter(values, window, polyorder, mode="interp")
    sigma = _local_noise(values, 8 * window)
    if float(np.median(sigma)) <= 1e-7 * max(float(np.ptp(values)), 1e-30):
        return values.copy()                      # noise-free data: any smoothing would only distort
    light = savgol_filter(values, 5, 2, mode="interp")
    light_gain = float(np.linalg.norm(savgol_coeffs(5, 2)))
    span = 25 * window + 1
    candidates, props = find_peaks(light, prominence=(None, None), wlen=span)
    if len(candidates) == 0:
        return strong
    peaks = candidates[props["prominences"] >= 8.0 * light_gain * sigma[candidates]]
    if len(peaks) == 0:
        return strong
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        widths = peak_widths(light, peaks, rel_height=0.5, wlen=span)[0]
    result = strong.copy()
    fine_cache = {}
    max_fine = min(n if n % 2 else n - 1, 4 * window + 1)

    def fine(size):
        if size not in fine_cache:
            length = max(size, window)
            c_strong = np.zeros(length); c_fine = np.zeros(length)
            a = (length - window) // 2
            c_strong[a:a + window] = savgol_coeffs(window, polyorder)
            if size == 1:
                filtered = values
                c_fine[length // 2] = 1.0
            else:
                filtered = savgol_filter(values, size, 4, mode="interp")
                b = (length - size) // 2
                c_fine[b:b + size] = savgol_coeffs(size, 4)
            fine_cache[size] = (filtered, float(np.linalg.norm(c_fine - c_strong)))
        return fine_cache[size]

    for centre, width in zip(peaks, widths):
        width = float(width)
        if width <= 0 or _profile_retention(width, window, polyorder) >= keep:
            continue                              # the strong filter does not harm this peak
        size = _widest_safe_window(width, max_fine, keep) or 1
        filtered, gain = fine(size)
        half = int(math.ceil(1.5 * max(window, size) + 2.5 * width))
        lo, hi = max(0, centre - half), min(n, centre + half + 1)
        sigma_diff = gain * sigma[lo:hi]
        difference = np.abs(filtered[lo:hi] - strong[lo:hi])
        weight = np.clip((difference - 3.0 * sigma_diff) / (3.0 * sigma_diff), 0.0, 1.0)
        weight = maximum_filter1d(weight, size=window, mode="nearest")
        weight = np.clip(uniform_filter1d(weight, size=max(3, size), mode="nearest"), 0.0, 1.0)
        taper = np.ones(hi - lo)
        edge = min(window, (hi - lo) // 2)
        if edge > 0:
            ramp = np.linspace(0.0, 1.0, edge, endpoint=False)
            taper[:edge] = ramp
            taper[-edge:] = ramp[::-1]
        weight *= taper
        result[lo:hi] = (1.0 - weight) * result[lo:hi] + weight * filtered[lo:hi]
    return np.clip(result, float(values.min()), float(values.max()))


def smooth_for_display(x, y, strength_percent, window_deg=0.25, polyorder=2, method="savgol"):
    """Return an independent y display array, aligned with the original x.

    alpha = strength_percent / 100; display = (1-alpha)*y + alpha*F(y), where F is
    method="savgol": plain Savitzky-Golay (nominal 0.25 degree window, quadratic by default), or
    method="peak_preserving": the same noise reduction in the background, but narrow peaks keep
    their height (see peak_preserving_filter).
    For irregular x, interpolate to a uniform grid at approximately median dx,
    apply SG, and interpolate the filtered values back to the original x.
    The auxiliary grid has at most 4*N points to avoid extreme-gap allocation.
    Endpoints use scipy's mode='interp'. No clipping, no data mutation.
    Nonfinite/nonmonotonic/short data are returned unchanged on a copy.
    """
    x = np.asarray(x, dtype=float)
    original = np.array(y, dtype=float, copy=True)
    settings = display_smoothing_settings({"smooth_strength_percent": strength_percent,
                                           "smooth_window_deg": window_deg,
                                           "smooth_polyorder": polyorder})
    alpha = settings["smooth_strength_percent"] / 100.0
    if alpha == 0 or x.ndim != 1 or original.ndim != 1 or len(x) != len(original) or len(x) < 5:
        return original
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(original)):
        return original
    differences = np.diff(x)
    if not np.all(differences > 0):
        return original
    step = float(np.median(differences))
    uniform = np.allclose(differences, step, rtol=1e-6, atol=1e-12)
    if uniform:
        grid, values = x, original
    else:
        estimated = (float(x[-1]) - float(x[0])) / step
        count = min(4 * len(x), max(5, int(min(estimated, 4 * len(x))) + 1))
        grid = np.linspace(x[0], x[-1], count)
        values = np.interp(grid, x, original)
        step = float(grid[1] - grid[0])
    window = display_smoothing_window(len(grid), step, settings["smooth_window_deg"], settings["smooth_polyorder"])
    if window == 0:
        return original
    filtered = None
    if method == "peak_preserving":
        try:
            filtered = peak_preserving_filter(values, step, settings["smooth_window_deg"], settings["smooth_polyorder"])
        except (ValueError, np.linalg.LinAlgError):
            filtered = None                       # fall back to the plain filter instead of a broken line
    if filtered is None:
        filtered = savgol_filter(values, window, settings["smooth_polyorder"], mode="interp")
    if not uniform:
        filtered = np.interp(x, grid, filtered)
    return (1.0 - alpha) * original + alpha * filtered


def sanitize_plot_state(state):
    """Repair each optional scalar separately, retaining all valid siblings."""
    result = deepcopy(state) if isinstance(state, dict) else {}
    numeric = {
        "mode": (0, 0, 2), "data_mode": (0, 0, 3),
        "vertical_spacing": (1.15, 0, 20), "line_width": (1.2, 0.2, 5),
        "legend_font_size": (8, 5, 24), "peak_font_size": (8, 5, 24),
        "title_font_size": (12, 6, 36), "axis_font_size": (10, 6, 28),
        "tick_font_size": (9, 5, 24), "sample_font_size": (9, 5, 24),
        "width_in": (11, 4, 30), "height_in": (7.5, 3, 30),
        "group_break_after": (1, 0, 1000000),
    }
    for name, (default, low, high) in numeric.items():
        if name in result:
            result[name] = finite_number(result[name], default, low, high)
    for name, default in {"keep_peak_height_with_spacing": True, "legend": True,
                          "sample_names": True, "grid": True, "peak_markers": False,
                          "peak_labels": False, "phase_marker_strip": False,
                          "group_separator": False}.items():
        if name in result:
            result[name] = optional_bool(result[name], default)
    positions = result.get("label_positions", {})
    result["label_positions"] = {}
    if isinstance(positions, dict):
        for key, value in positions.items():
            if not isinstance(value, dict) or value.get("coord") not in ("data", "axes_y"):
                continue
            x = finite_number(value.get("x"), None)
            y = finite_number(value.get("y"), None)
            if x is not None and y is not None:
                result["label_positions"][str(key)] = dict(value, x=x, y=y)
    manual = result.get("manual_legend_loc")
    if not isinstance(manual, (list, tuple)) or len(manual) != 2 or any(finite_number(v, None) is None for v in manual):
        result["manual_legend_loc"] = None
    separators = result.get("separators", [])
    result["separators"] = []
    if isinstance(separators, list):
        for item in separators:
            if not isinstance(item, dict):
                continue
            sep = dict(item)
            sep["position"] = int(finite_number(sep.get("position", 0), 0, 0, 1000000))
            sep["offset_pct"] = finite_number(sep.get("offset_pct", 0), 0, -48, 48)
            sep["color"] = valid_color(sep.get("color"), "#808080")
            sep["labels_visible"] = optional_bool(sep.get("labels_visible", True), True)
            result["separators"].append(sep)
    return result


def valid_color(value, default):
    return value if isinstance(value, str) and is_color_like(value) else default


def safe_window_rect(saved_rect, screens, default_size=(1500, 900), minimum_size=(900, 600)):
    """Clamp a normal-window (x,y,w,h) into one available screen (pure logic).

    Screens use Qt logical coordinates, including negative monitor origins.
    Prefer greatest overlap, otherwise the closest screen centre. Invalid saved
    geometry falls back to a centred 90%-of-screen normal window.
    """
    valid_screens = []
    for screen in screens:
        if isinstance(screen, (list, tuple)) and len(screen) == 4:
            values = [finite_number(v, None) for v in screen]
            if None not in values and values[2] > 0 and values[3] > 0:
                valid_screens.append(tuple(int(v) for v in values))
    if not valid_screens:
        valid_screens = [(0, 0, 1280, 720)]
    saved = None
    if isinstance(saved_rect, (list, tuple)) and len(saved_rect) == 4:
        values = [finite_number(v, None) for v in saved_rect]
        if None not in values and values[2] > 0 and values[3] > 0:
            saved = tuple(int(v) for v in values)
    screen = valid_screens[0]
    if saved:
        x, y, w, h = saved
        def rank(candidate):
            sx, sy, sw, sh = candidate
            overlap = max(0, min(x + w, sx + sw) - max(x, sx)) * max(0, min(y + h, sy + sh) - max(y, sy))
            distance = (x + w / 2 - sx - sw / 2) ** 2 + (y + h / 2 - sy - sh / 2) ** 2
            return overlap, -distance
        screen = max(valid_screens, key=rank)
    sx, sy, sw, sh = screen
    if saved is None:
        w, h = min(default_size[0], int(sw * 0.9)), min(default_size[1], int(sh * 0.9))
        x, y = sx + (sw - w) // 2, sy + (sh - h) // 2
    else:
        x, y, w, h = saved
    w = min(sw, max(min(minimum_size[0], sw), w))
    h = min(sh, max(min(minimum_size[1], sh), h))
    x = max(sx, min(x, sx + sw - w))
    y = max(sy, min(y, sy + sh - h))
    return int(x), int(y), int(w), int(h)
