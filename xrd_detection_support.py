"""Peak candidate selection (pure numpy/scipy, no GUI).

The original rule demanded a prominence of at least ``min_snr`` times the point-to-point noise. That
is far too strict for broad reflections: a peak ten samples wide is much easier to see than its
single-point prominence suggests, and a reflection sitting on the flank of a strong one (a ferrite
110 next to an austenite 111) was lost at 6.4 sigma although its real significance is far higher.

Rule used now: a candidate is accepted when EITHER it meets the original criterion (so nothing
that used to be detected disappears) OR all of the following hold: prominence >= min(5, min_snr)
sigma, width-aware significance prominence * sqrt(max(1, 0.752 * width_in_samples)) / sigma >=
min_snr, and a multi-scale matched-filter significance >= min_snr (this last test removes ripples
and broad background residuals that are not white noise). On pure noise this admits about one false
peak per 200 scans of 2000 points.
"""
from __future__ import annotations

import math
import warnings
from typing import Tuple

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks, peak_prominences, peak_widths

WIDE_PEAK_FLOOR_SIGMA = 5.0
MATCHED_SCALES = (1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0)   # Gaussian sigma in samples


def noise_sigma_from_differences(sig: np.ndarray) -> float:
    """Robust white-noise estimate from first differences (insensitive to slow background)."""
    dif = np.diff(np.asarray(sig, dtype=float))
    if not len(dif):
        return 1e-12
    med = float(np.nanmedian(dif))
    mad = float(np.nanmedian(np.abs(dif - med)))
    return max(1e-12, 1.4826 * mad / math.sqrt(2.0))


def matched_significance(sig: np.ndarray, sigma: float, candidates) -> np.ndarray:
    """Best significance over Gaussian scales: smoothed prominence / noise of the smoothed curve."""
    sig = np.asarray(sig, dtype=float)
    cand = np.asarray(candidates, dtype=int)
    n = len(sig)
    best = np.zeros(len(cand))
    if not len(cand):
        return best
    for scale in MATCHED_SCALES:
        smooth = gaussian_filter1d(sig, scale, mode="nearest")
        gain = 1.0 / math.sqrt(2.0 * math.sqrt(math.pi) * scale)   # std of smoothed white noise / sigma
        radius = max(2, int(math.ceil(1.5 * scale)))
        near = np.empty(len(cand), dtype=int)
        for k, centre in enumerate(cand):
            lo, hi = max(0, centre - radius), min(n, centre + radius + 1)
            near[k] = lo + int(np.argmax(smooth[lo:hi]))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            prominence = peak_prominences(smooth, near)[0]
        best = np.maximum(best, prominence / (max(sigma, 1e-12) * gain))
    return best


def find_candidate_peaks(sig, noise_sigma: float, scale: float, prominence_pct: float, min_snr: float,
                         min_rel_height_pct: float) -> Tuple[np.ndarray, np.ndarray]:
    """Return (indices, prominences) of accepted peaks of the baseline-subtracted signal."""
    sig = np.asarray(sig, dtype=float)
    prom_abs = max(prominence_pct / 100.0 * scale, min_snr * noise_sigma)
    height_abs = min_rel_height_pct / 100.0 * scale
    floor_sigma = min(WIDE_PEAK_FLOOR_SIGMA, float(min_snr))
    loose = min(prom_abs, max(prominence_pct / 100.0 * scale, floor_sigma * noise_sigma))
    idx, props = find_peaks(sig, prominence=loose, height=height_abs)
    if not len(idx):
        return idx, np.zeros(0)
    prominences = np.asarray(props["prominences"], dtype=float)
    strong = prominences >= prom_abs
    if strong.all():
        return idx, prominences
    widths_pts = peak_widths(sig, idx, rel_height=0.5)[0]
    boost = np.sqrt(np.maximum(1.0, 0.752 * widths_pts))
    extra = (~strong) & (prominences >= floor_sigma * noise_sigma) & (prominences * boost / noise_sigma >= min_snr)
    if extra.any():
        extra[extra] = matched_significance(sig, noise_sigma, idx[extra]) >= min_snr
    keep = strong | extra
    return idx[keep], prominences[keep]
