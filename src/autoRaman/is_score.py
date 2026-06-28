"""
is_score.py -- IS-Score (Integrity Spectrum Score) for AutoRaman.

A ground-truth-free, 0..1 metric that judges the quality of a Raman baseline
correction directly from ``(raw spectrum, baseline-corrected spectrum, axis)``.
A value near 1 means a good fit; near 0 means poor (over- or under-fitting).

This is a **clean-room re-implementation** of the method described in:

    S. Innocente, A. Visentin, S. Andersson-Engels, K. Komolibus, R. Gautam,
    "Automated Baseline Correction Evaluation Score for Raman Spectroscopy",
    ACS Omega 2026, 11, 25057-25068. DOI: 10.1021/acsomega.5c09870

Only the *method* (the published equations) is reproduced here, from the paper;
no source code from the authors' repository was used. The module depends only
on numpy / scipy so it stays GUI-free and headless-testable, consistent with
the rest of the ``autoRaman`` package.

The score is the sum of five penalty blocks subtracted from 1
(IS = clip(1 - sum(penalties), 0, 1)):

    1. Single peak/dip   -- baseline must not cut into peaks or float off dips.
    2. Band region (RRP) -- extend (1) across each band's full region.
    3. AUC / over-fit    -- compare against a deliberately over-fitted baseline.
    4. Mean ratio        -- dips left above the baseline signal underfitting.
    5. Intensity         -- baseline rising above the spectrum is penalized.

Where the paper leaves a constant or sub-step underspecified, a physically
motivated, bounded choice is made and noted in a comment. The goal -- as for
ranking and parameter optimisation -- is the correct *ordering* of baselines
(good > over/under-fit), not byte-identical reproduction of the authors' scores.

The block weights (:data:`PENALTY_WEIGHTS`) are calibration knobs.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np
from scipy.signal import find_peaks, savgol_filter

from .quality_metrics import _robust_noise_rms

# np.trapz was renamed to np.trapezoid in NumPy 2.0 (trapz removed in 2.x).
_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))


# ---------------------------------------------------------------------------
# Result type and tunable weights
# ---------------------------------------------------------------------------

@dataclass
class ISScoreResult:
    """IS-Score plus the per-block penalty breakdown (for transparency / GUI)."""
    score: float
    single_band: float = 0.0
    band_region: float = 0.0
    auc: float = 0.0
    mean_ratio: float = 0.0
    intensity: float = 0.0
    n_peaks: int = 0
    n_dips: int = 0
    note: str = ""
    weights: dict = field(default_factory=dict)

    def breakdown(self) -> str:
        return (
            "IS-Score {:.3f}  [single {:.3f}, region {:.3f}, auc {:.3f}, "
            "mean-ratio {:.3f}, intensity {:.3f}]  ({} peaks / {} dips)".format(
                self.score, self.single_band, self.band_region, self.auc,
                self.mean_ratio, self.intensity, self.n_peaks, self.n_dips)
        )


# Per-block weights. Defaults calibrated against the synthetic over/good/under
# baselines in autoRaman.tests so a good correction lands ~0.8-0.9 with clear
# separation from over- and under-fits. Tune here, not in the block code.
PENALTY_WEIGHTS = {
    "single_band": 0.30,   # peak-eating (overfit) -- strongest signal
    "band_region": 0.20,   # region-wide overfit
    "auc": 0.15,           # global over-correction (negative area)
    "mean_ratio": 0.25,    # leftover background at dips (underfit)
    "intensity": 0.10,     # baseline above spectrum
}

# A band must clear this prominence (fraction of the [0,1]-normalised span) to
# count -- mirrors the paper's prominence gate that ignores insignificant bumps.
_MIN_PROMINENCE = 0.02

# Band counted only if it recurs in >= this many denoised iterations (paper: >=2).
_MIN_DENOISED_AGREEMENT = 2

# Positions in raw vs denoised match within this many cm^-1 (paper: 5).
_POS_TOL_CM = 5.0


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def _odd(n: int) -> int:
    n = int(n)
    return n if n % 2 == 1 else n + 1


def _smooth(y: np.ndarray, win: int, poly: int = 3) -> np.ndarray:
    """Savitzky-Golay smoothing with safe window/order handling."""
    n = y.size
    if n < 5:
        return y.copy()
    win = _odd(win)
    win = min(win, n if n % 2 == 1 else n - 1)
    if win <= poly:
        win = _odd(poly + 2)
        if win > n:
            return y.copy()
    return savgol_filter(y, win, poly)


def _normalize(a: np.ndarray, lo: float, hi: float) -> np.ndarray:
    span = hi - lo
    if span <= 0:
        return np.zeros_like(a)
    return (a - lo) / span


@dataclass
class _Band:
    apex: int          # sample index of the band's apex (peak or dip)
    left: int          # left edge sample index
    right: int         # right edge sample index
    is_peak: bool      # True = peak, False = dip


def _detect_bands(s: np.ndarray, axis: np.ndarray, invert: bool) -> List[_Band]:
    """Robustly detect peaks (or dips, if ``invert``) on a [0,1] signal.

    Robustness (paper "Bands and Edges Detection"): a band is kept only if it
    recurs across multiple Savitzky-Golay-denoised versions of the signal,
    within +/-5 cm^-1 of the raw detection. Edges are the prominence contour
    bases returned by scipy.find_peaks (the local minima bounding the band).
    """
    sig = -s if invert else s
    span = float(np.ptp(sig)) or 1.0
    min_prom = _MIN_PROMINENCE * span

    idx, props = find_peaks(sig, prominence=min_prom)
    if idx.size == 0:
        return []
    left_bases = props["left_bases"]
    right_bases = props["right_bases"]

    denoised_positions = []
    for win in (7, 11, 15):
        di, _ = find_peaks(_smooth(sig, win), prominence=min_prom)
        denoised_positions.append(axis[di] if di.size else np.empty(0))

    bands: List[_Band] = []
    for k, i in enumerate(idx):
        pos = axis[i]
        agree = sum(
            bool(dp.size and np.any(np.abs(dp - pos) <= _POS_TOL_CM))
            for dp in denoised_positions
        )
        if agree >= _MIN_DENOISED_AGREEMENT:
            bands.append(_Band(int(i), int(left_bases[k]), int(right_bases[k]),
                               not invert))
    return bands


def _joint_prominence(I: np.ndarray, b: _Band) -> float:
    """Prominence of a band measured in the (joint-normalised) frame ``I``.

    For a peak: apex height above the higher of its two edges.
    For a dip:  apex depth below the lower of its two edges.
    """
    apex, lo, hi = I[b.apex], I[b.left], I[b.right]
    if b.is_peak:
        return float(apex - max(lo, hi))
    return float(min(lo, hi) - apex)


# ---------------------------------------------------------------------------
# Penalty block 1 -- single peak / dip
# ---------------------------------------------------------------------------

def _penalty_single_band(I: np.ndarray, B: np.ndarray,
                         peaks: List[_Band], dips: List[_Band]) -> float:
    """Eqs 2-6: flag bands where the baseline cuts into a peak (overfit) or
    sits wrong at a dip (under/over-fit), then aggregate adaptively."""
    total = len(peaks) + len(dips)
    if total == 0:
        return 0.0

    penalized = 0
    for p in peaks:
        prom = _joint_prominence(I, p)
        if prom <= 0:
            continue
        # overfit: baseline rises above (peak - 75% prominence) -> eats the peak.
        if B[p.apex] > I[p.apex] - 0.75 * prom:
            penalized += 1
    for d in dips:
        prom = _joint_prominence(I, d)
        if prom <= 0:
            continue
        # underfit: baseline well below the valley (leftover background);
        # overfit: baseline well above the valley (over-correction). 50% gate.
        if abs(I[d.apex] - B[d.apex]) > 0.5 * prom:
            penalized += 1

    if penalized == 0:
        return 0.0

    ratio = penalized / total                          # eq 2
    beta = max(1.0, total * (1.0 - ratio))             # eq 3
    w = total / (total + beta)                         # eq 4
    log_term = np.log2(1.5 + ratio)                    # eq 5
    return float(w * log_term)                         # eq 6


# ---------------------------------------------------------------------------
# Penalty block 2 -- band region (Raman Region Prominence)
# ---------------------------------------------------------------------------

def _penalty_band_region(I: np.ndarray, B: np.ndarray,
                         peaks: List[_Band], dips: List[_Band]) -> float:
    """Eqs 7-15: extend the single-band idea across each peak's full region
    using the Raman Region Prominence (a baseline-scaled pseudo-prominence).

    The RRP at the apex distance ``x`` and true prominence ``y`` scales the
    region's intensity-to-baseline distance into a pseudo-prominence field
    ``rrp = (y/x) * |I-B|``. ``rrp_dist = |rrp - (I-B)| = |I-B| * |y/x - 1|``
    is then ~0 when the baseline sits a peak-consistent distance below the
    signal across the whole region, and grows when the baseline cuts into the
    band (overfit, ``x`` small -> ``y/x`` large) or floats off it.

    Applied to peaks only: at dips ``x = I-B -> 0`` makes ``y/x`` blow up, and
    dip mis-fit is already captured by the single-band and mean-ratio blocks.
    The paper's 75th-percentile filtering keeps the penalty robust to a few
    outlier points.
    """
    if not peaks:
        return 0.0

    per_region = []
    for b in peaks:
        lo, hi = b.left, b.right
        if hi - lo < 2:
            continue
        x = I[b.apex] - B[b.apex]                       # eq 7 (apex distance)
        y = _joint_prominence(I, b)                     # actual prominence
        if x <= 1e-6 or y <= 0:
            continue
        diff = I[lo:hi + 1] - B[lo:hi + 1]              # eq 9 (signed)
        rrp = (y / x) * np.abs(diff)                    # eqs 10-11
        rrp_dist = np.abs(rrp - np.abs(diff))           # eq 12 = |diff|*|y/x-1|
        if rrp_dist.size:
            thr = np.percentile(rrp_dist, 75)           # keep robust subset
            kept = rrp_dist[rrp_dist <= thr]
            if kept.size:
                per_region.append(float(np.mean(kept)))
    if not per_region:
        return 0.0
    # mean region deviation, already in [0,1]-normalised intensity units.
    return float(np.clip(np.mean(per_region), 0.0, 1.0))


# ---------------------------------------------------------------------------
# Penalty block 3 -- AUC vs an artificial over-fitted baseline
# ---------------------------------------------------------------------------

def _penalty_auc(I: np.ndarray, B: np.ndarray, axis: np.ndarray) -> float:
    """Global, area-based overfitting penalty (paper AUC section, robust form).

    The paper compares the baseline-corrected AUC against that of a
    deliberately over-fitted reference baseline. Reproducing that reference
    faithfully (multi-scale Savitzky-Golay + spline + Gaussian bumps) proved
    numerically fragile, so this uses the directly observable consequence of
    over-fitting instead: when a baseline cuts into bands, the corrected
    spectrum develops *negative area* (the baseline rose above the signal and
    removed real intensity). The penalty is that negative area as a fraction of
    the total corrected area -- ~0 for a good fit, growing with over-correction,
    and bounded in [0,1].
    """
    corr = I - B
    neg_area = _trapz(np.clip(-corr, 0, None), axis)
    pos_area = _trapz(np.clip(corr, 0, None), axis)
    denom = neg_area + pos_area
    if denom <= 0:
        return 0.0
    return float(np.clip(neg_area / denom, 0.0, 1.0))


# ---------------------------------------------------------------------------
# Penalty block 4 -- mean ratio (dips above baseline => underfit)
# ---------------------------------------------------------------------------

def _penalty_mean_ratio(I: np.ndarray, B: np.ndarray) -> float:
    """Eqs 22-23: across several denoised versions, measure how much spectral
    dip content is left *above* the baseline (leftover background = underfit).
    """
    noise = max(_robust_noise_rms(I - B), 1e-9)
    span = float(np.ptp(I)) or 1.0
    fracs = []
    for win in (8, 16, 32, 40):
        sm = _smooth(I, win, poly=4)
        di, _ = find_peaks(-sm, prominence=0.005 * span)
        if di.size == 0:
            continue
        # dips lying clearly above the baseline -> uncorrected background.
        above = (I[di] - B[di]) > noise
        fracs.append(float(np.mean(above)))
    if not fracs:
        return 0.0
    return float(np.clip(np.mean(fracs), 0.0, 1.0))


# ---------------------------------------------------------------------------
# Penalty block 5 -- intensity (baseline above spectrum)
# ---------------------------------------------------------------------------

def _penalty_intensity(I: np.ndarray, B: np.ndarray) -> float:
    """Noise-thresholded count of points where the baseline exceeds the
    spectrum, normalised by the number of points (paper "Intensity")."""
    noise = max(_robust_noise_rms(I - B), 1e-9)
    return float(np.mean(B > I + noise))


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def is_score(
    raw: np.ndarray,
    corrected: np.ndarray,
    axis: np.ndarray,
    *,
    weights: Optional[dict] = None,
    return_breakdown: bool = False,
):
    """Compute the IS-Score for one baseline correction.

    Parameters
    ----------
    raw        : measured Raman intensity I(v).
    corrected  : baseline-corrected intensity (raw - baseline).
    axis       : Raman shift / wavenumber axis (same length).
    weights    : optional override of :data:`PENALTY_WEIGHTS`.
    return_breakdown : return an :class:`ISScoreResult` instead of a float.

    Returns
    -------
    float in [0, 1] (higher = better), or :class:`ISScoreResult` if
    ``return_breakdown`` is True. Never raises on degenerate input -- it
    returns a defined, neutral-to-low score and a ``note``.
    """
    raw = np.asarray(raw, dtype=float).flatten()
    corrected = np.asarray(corrected, dtype=float).flatten()
    axis = np.asarray(axis, dtype=float).flatten()
    w = {**PENALTY_WEIGHTS, **(weights or {})}

    def _result(score, note="", **kw):
        res = ISScoreResult(score=float(np.clip(score, 0.0, 1.0)),
                            note=note, weights=w, **kw)
        return res if return_breakdown else res.score

    n = raw.size
    if n < 5 or corrected.size != n or axis.size != n:
        return _result(0.0, note="degenerate input (too short / length mismatch)")

    baseline = raw - corrected

    # --- normalisations -----------------------------------------------------
    s_lo, s_hi = float(np.min(raw)), float(np.max(raw))
    s = _normalize(raw, s_lo, s_hi)                    # spectrum-only (detection)
    if s_hi - s_lo <= 0:
        return _result(0.0, note="flat spectrum")

    g_lo = float(min(np.min(raw), np.min(baseline)))  # joint (evaluation)
    g_hi = float(max(np.max(raw), np.max(baseline)))
    I = _normalize(raw, g_lo, g_hi)
    B = _normalize(baseline, g_lo, g_hi)

    # --- feature detection (on spectrum-only normalised signal) -------------
    peaks = _detect_bands(s, axis, invert=False)
    dips = _detect_bands(s, axis, invert=True)

    # --- penalty blocks -----------------------------------------------------
    p1 = _penalty_single_band(I, B, peaks, dips)
    p2 = _penalty_band_region(I, B, peaks, dips)
    p3 = _penalty_auc(I, B, axis)
    p4 = _penalty_mean_ratio(I, B)
    p5 = _penalty_intensity(I, B)

    total = (w["single_band"] * p1 + w["band_region"] * p2 + w["auc"] * p3
             + w["mean_ratio"] * p4 + w["intensity"] * p5)
    score = 1.0 - total

    note = "" if peaks else "no peaks detected -- score less reliable"
    return _result(score, note=note,
                   single_band=p1, band_region=p2, auc=p3,
                   mean_ratio=p4, intensity=p5,
                   n_peaks=len(peaks), n_dips=len(dips))
