"""
auto_peaks.py -- automated peak detection for AutoRaman.

Strategy (roadmap Phase 2b):
    * Primary detector: scipy.signal.find_peaks with prominence/width gating
      derived from an estimated noise floor (robust, no manual threshold).
    * Alternative detectors (CWT, 2nd-derivative) are scaffolded for the
      challenging cases (shoulders, overlapping bands, low SNR) the roadmap
      calls out.
    * Each peak carries a confidence score (quality_metrics.peak_confidence)
      and the parameters used, so the user can carry the proposal into the
      manual peak-fitting dialog and adjust.

Operates on a *baseline-corrected* spectrum (run AutoBaseline first, or pass an
already-corrected y). Like the rest of the package it is GUI-free.

Status: SCAFFOLD. find_peaks path is functional; CWT and 2nd-derivative paths
are stubbed with clear TODOs.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

import numpy as np
from scipy.signal import find_peaks, peak_widths

from .quality_metrics import peak_confidence, _robust_noise_rms


@dataclass
class Peak:
    """A detected peak. Mirrors what the manual fit dialog seeds a band with."""
    position: float        # x (wavenumber) at the apex
    index: int             # sample index of the apex
    height: float          # corrected intensity at the apex
    prominence: float      # prominence above local baseline
    fwhm: float            # full width at half maximum (in x units)
    confidence: float      # [0,1] from peak_confidence
    method: str = "find_peaks"


class AutoPeaks:
    """Detect peaks in a (baseline-corrected) spectrum."""

    def __init__(
        self,
        min_confidence: float = 0.5,
        prominence_sigma: float = 3.0,
    ):
        """
        Parameters
        ----------
        min_confidence  : drop peaks below this confidence (default 0.5).
        prominence_sigma: required prominence as a multiple of noise RMS
                          (default 3 sigma) -- the core noise-adaptive gate.
        """
        self.min_confidence = min_confidence
        self.prominence_sigma = prominence_sigma

    # -- primary detector ---------------------------------------------------
    def detect(
        self,
        x: np.ndarray,
        y_corr: np.ndarray,
    ) -> List[Peak]:
        """Detect peaks via prominence-gated find_peaks.

        The prominence threshold is set from the spectrum's own noise floor, so
        no manual threshold is required -- this is the key "no clicking on
        peaks" automation.
        """
        x = np.asarray(x, dtype=float)
        y = np.asarray(y_corr, dtype=float)

        noise = _robust_noise_rms(y)
        min_prom = self.prominence_sigma * noise if noise > 0 else None

        idx, props = find_peaks(y, prominence=min_prom)
        if idx.size == 0:
            return []

        # widths at half-prominence -> convert sample units to x units.
        widths_samp, _, _, _ = peak_widths(y, idx, rel_height=0.5)
        dx = float(np.median(np.abs(np.diff(x)))) if x.size > 1 else 1.0

        peaks: List[Peak] = []
        for k, i in enumerate(idx):
            prom = float(props["prominences"][k])
            fwhm = float(widths_samp[k] * dx)
            conf = peak_confidence(prominence=prom, noise_rms=noise, width=fwhm)
            if conf < self.min_confidence:
                continue
            peaks.append(
                Peak(
                    position=float(x[i]),
                    index=int(i),
                    height=float(y[i]),
                    prominence=prom,
                    fwhm=fwhm,
                    confidence=conf,
                )
            )
        peaks.sort(key=lambda p: p.position)
        return peaks

    # -- alternative detectors (scaffold) -----------------------------------
    def detect_cwt(self, x, y_corr, widths=None) -> List[Peak]:  # pragma: no cover
        """Continuous Wavelet Transform detector (scipy.find_peaks_cwt).

        TODO (Phase 2b): noise-robust multi-scale detection for low-SNR spectra
        (see test_data/low_snr_noisy.txt). Needs width-range selection and a
        confidence mapping from CWT ridge strength.
        """
        raise NotImplementedError("CWT detector is a Phase 2b stub.")

    def detect_second_derivative(self, x, y_corr) -> List[Peak]:  # pragma: no cover
        """Second-derivative detector for resolving shoulders/overlaps.

        TODO (Phase 2b): smooth (Savitzky-Golay), take -y'', find zero-crossings
        / negative minima to surface shoulder peaks find_peaks misses.
        """
        raise NotImplementedError("Second-derivative detector is a Phase 2b stub.")


class PeakValidator:
    """Filter false positives (cosmic spikes, noise) from a peak list.

    TODO (Phase 2b): use width and neighbour-symmetry to reject single-bin
    cosmic spikes (test_data/cosmic_spike.txt) distinct from real bands, beyond
    the basic confidence gate already applied in AutoPeaks.detect.
    """

    def filter(self, peaks: List[Peak]) -> List[Peak]:  # pragma: no cover
        raise NotImplementedError("PeakValidator is a Phase 2b stub.")
