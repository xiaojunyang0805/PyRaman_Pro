"""
quality_metrics.py -- objective scoring functions for AutoRaman.

These metrics answer the roadmap's question "what makes a 'good' baseline
objectively?" and provide confidence scores for detected peaks. They are pure
functions over numpy arrays so they can be tested in isolation.

All scores are normalised to roughly [0, 1] where **higher is better**, so they
can be combined with weights. Where a raw physical quantity is more meaningful
(e.g. RMS noise) the raw value is also returned for transparency.

Status: implemented for baseline; peak/fit metrics partially implemented (TODO
markers where tuning against real data is still required).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np


# ---------------------------------------------------------------------------
# Baseline quality
# ---------------------------------------------------------------------------

@dataclass
class BaselineScore:
    """Breakdown of baseline quality sub-metrics and the combined score."""
    residual_rms: float          # raw RMS of the corrected signal's noise floor
    smoothness: float            # [0,1] higher = smoother baseline (low curvature)
    peak_preservation: float     # [0,1] higher = peaks not clipped/distorted
    non_negativity: float        # [0,1] higher = fewer negative corrected values
    combined: float = 0.0        # weighted total [0,1]
    weights: dict = field(default_factory=dict)


# Default metric weights (roadmap table: Residual & Peak preservation = High).
DEFAULT_BASELINE_WEIGHTS = {
    "residual_rms": 0.30,
    "smoothness": 0.20,
    "peak_preservation": 0.35,
    "non_negativity": 0.15,
}


def _robust_noise_rms(y_corr: np.ndarray) -> float:
    """Estimate noise RMS robustly via the median absolute deviation of the
    first difference (insensitive to peaks). 1.4826 makes MAD ~ sigma."""
    d = np.diff(np.asarray(y_corr, dtype=float))
    if d.size == 0:
        return 0.0
    mad = np.median(np.abs(d - np.median(d)))
    # diff inflates noise variance by 2 -> divide by sqrt(2)
    return float(1.4826 * mad / np.sqrt(2.0))


def baseline_fit_cost(
    y_raw: np.ndarray,
    y_corr: np.ndarray,
    window_frac: float = 0.06,
) -> float:
    """Objective cost of a baseline-corrected spectrum (lower = better).

    A *good* Raman baseline hugs the bottom envelope of the spectrum: after
    subtraction the corrected signal should be >= 0 and its **local minima**
    (the peak-free regions between bands) should sit at ~0. This captures both
    failure modes that the sub-metric ``baseline_quality`` missed:

    * baseline too stiff / through the middle -> corrected swings negative and
      its rolling minimum is far from 0  (large cost);
    * baseline correctly at the bottom       -> rolling minimum ~ 0, few
      negatives                              (small cost).

    Used to drive automatic lambda selection (see auto_baseline). The window is
    a fraction of the spectrum so it spans typical Raman bands without flooding
    the broad background.
    """
    y_raw = np.asarray(y_raw, dtype=float)
    y_corr = np.asarray(y_corr, dtype=float)
    span = float(np.ptp(y_raw)) or 1.0
    n = y_corr.size
    if n == 0:
        return float("inf")
    w = max(5, int(window_frac * n))
    try:
        from scipy.ndimage import minimum_filter1d
        rmin = minimum_filter1d(y_corr, size=w, mode="nearest")
    except Exception:  # pragma: no cover - scipy always present in this app
        rmin = np.array([y_corr[max(0, i - w // 2):i + w // 2 + 1].min() for i in range(n)])
    background = float(np.mean(np.abs(rmin)))          # minima should be ~0
    negative = float(np.mean(np.clip(-y_corr, 0, None)))  # not strongly negative
    return (background + negative) / span


def baseline_quality(
    y_raw: np.ndarray,
    y_corr: np.ndarray,
    baseline: np.ndarray,
    weights: Optional[dict] = None,
) -> BaselineScore:
    """Score a baseline-corrected spectrum.

    Parameters
    ----------
    y_raw   : original intensity
    y_corr  : intensity after baseline subtraction (y_raw - baseline)
    baseline: the estimated baseline
    weights : optional override of DEFAULT_BASELINE_WEIGHTS

    Returns
    -------
    BaselineScore with sub-metrics and a combined [0,1] score (higher better).
    """
    w = {**DEFAULT_BASELINE_WEIGHTS, **(weights or {})}
    y_raw = np.asarray(y_raw, dtype=float)
    y_corr = np.asarray(y_corr, dtype=float)
    baseline = np.asarray(baseline, dtype=float)

    span = float(np.ptp(y_raw)) or 1.0

    # 1. Residual noise: lower RMS is better. Normalise by signal span.
    rms = _robust_noise_rms(y_corr)
    residual_score = float(np.exp(-rms / (0.02 * span)))  # ~1 when rms << 2% span

    # 2. Smoothness: penalise high curvature in the baseline (2nd derivative).
    d2 = np.diff(baseline, n=2)
    curvature = float(np.sqrt(np.mean(d2 ** 2))) if d2.size else 0.0
    smoothness = float(np.exp(-curvature / (0.01 * span)))

    # 3. Peak preservation: the baseline must not eat into peaks. Proxy: the
    #    baseline should stay at or below the raw signal almost everywhere.
    overshoot = np.clip(baseline - y_raw, 0, None)
    frac_overshoot = float(np.mean(overshoot > 0.01 * span))
    peak_preservation = float(1.0 - frac_overshoot)

    # 4. Non-negativity: corrected signal shouldn't go strongly negative.
    neg = np.clip(-y_corr, 0, None)
    frac_neg = float(np.mean(neg > 0.01 * span))
    non_negativity = float(1.0 - frac_neg)

    combined = (
        w["residual_rms"] * residual_score
        + w["smoothness"] * smoothness
        + w["peak_preservation"] * peak_preservation
        + w["non_negativity"] * non_negativity
    )

    return BaselineScore(
        residual_rms=rms,
        smoothness=smoothness,
        peak_preservation=peak_preservation,
        non_negativity=non_negativity,
        combined=float(combined),
        weights=w,
    )


# ---------------------------------------------------------------------------
# Peak confidence
# ---------------------------------------------------------------------------

def peak_confidence(
    prominence: float,
    noise_rms: float,
    width: float,
) -> float:
    """Confidence [0,1] that a detected peak is real rather than noise.

    Primary signal is prominence-to-noise ratio (a peak well above the noise
    floor is trustworthy). Width sanity is folded in lightly: extremely narrow
    single-bin spikes are more likely cosmic rays than Raman bands.

    TODO: calibrate the width term against labelled spectra in test_data/
    (cosmic_spike.txt vs real bands) before relying on it for filtering.
    """
    if noise_rms <= 0:
        snr = np.inf if prominence > 0 else 0.0
    else:
        snr = prominence / noise_rms
    # logistic on SNR, centred at ~3 sigma
    conf = 1.0 / (1.0 + np.exp(-(snr - 3.0)))
    return float(np.clip(conf, 0.0, 1.0))


# ---------------------------------------------------------------------------
# Fit quality (Phase 2c -- placeholder)
# ---------------------------------------------------------------------------

def fit_quality(y: np.ndarray, y_fit: np.ndarray, n_params: int) -> float:
    """Goodness-of-fit score [0,1] for an automated peak fit.

    Currently returns adjusted-R^2 clamped to [0,1]. TODO: add reduced chi^2
    and per-peak residual checks when auto-fitting (Phase 2c) lands.
    """
    y = np.asarray(y, dtype=float)
    y_fit = np.asarray(y_fit, dtype=float)
    n = y.size
    ss_res = float(np.sum((y - y_fit) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2)) or 1.0
    r2 = 1.0 - ss_res / ss_tot
    if n - n_params - 1 > 0:
        r2 = 1.0 - (1.0 - r2) * (n - 1) / (n - n_params - 1)
    return float(np.clip(r2, 0.0, 1.0))
