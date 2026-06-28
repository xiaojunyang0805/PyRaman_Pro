"""
auto_analyzer.py -- one-call orchestrator for AutoRaman (roadmap Phase 2c).

Ties the pieces together for the "Auto Analyze" action:

    raw spectrum
        -> AutoBaseline.evaluate  (pick best baseline, keep runners-up)
        -> AutoPeaks.detect       (find peaks on the corrected signal)
        -> AnalysisReport         (everything the UI needs to show + let the
                                   user accept/override)

Crucially this returns a *report of proposals*, not a mutated spectrum. The
PyRamanGUI layer presents the report and the user accepts, tweaks, or ignores
it -- the manual baseline/fit dialogs remain the source of truth. Automation is
a fast first draft, not a replacement.

Status: SCAFFOLD. Baseline + peak orchestration works; auto-fitting (Phase 2c)
is left as a documented extension point.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np

from .auto_baseline import AutoBaseline, BaselineProposal
from .auto_peaks import AutoPeaks, Peak


@dataclass
class AnalysisReport:
    """Human- and UI-facing summary of an automated analysis run."""
    baseline: Optional[BaselineProposal]
    baseline_alternatives: List[BaselineProposal] = field(default_factory=list)
    peaks: List[Peak] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return self.baseline is not None and self.baseline.score is not None

    def summary(self) -> str:
        """Short text summary suitable for a status bar / log pane."""
        if not self.ok:
            return "AutoRaman: baseline correction failed -- use the manual dialog."
        b = self.baseline
        lines = [
            f"Baseline: {b.method}  (IS-Score {b.score.combined:.2f}, "
            f"noise RMS {b.score.residual_rms:.3g})",
            f"  params: {b.params}",
            f"Peaks detected: {len(self.peaks)}",
        ]
        for p in self.peaks[:12]:
            lines.append(
                f"  {p.position:8.1f} cm^-1  h={p.height:.3g}  "
                f"FWHM={p.fwhm:.1f}  conf={p.confidence:.2f}"
            )
        if len(self.peaks) > 12:
            lines.append(f"  ... and {len(self.peaks) - 12} more")
        for w in self.warnings:
            lines.append(f"  ! {w}")
        return "\n".join(lines)


class AutoAnalyzer:
    """High-level entry point: full automated baseline + peak analysis."""

    def __init__(
        self,
        auto_baseline: Optional[AutoBaseline] = None,
        auto_peaks: Optional[AutoPeaks] = None,
    ):
        self.auto_baseline = auto_baseline or AutoBaseline()
        self.auto_peaks = auto_peaks or AutoPeaks()

    def analyze(
        self,
        x: np.ndarray,
        y: np.ndarray,
        keep_alternatives: int = 3,
        optimize: bool = False,
    ) -> AnalysisReport:
        """Run baseline selection then peak detection; return a proposal report.

        Does not modify inputs and does not apply anything -- the caller (GUI)
        decides. ``keep_alternatives`` controls how many runner-up baselines to
        surface so the user can switch if the top pick looks wrong.
        ``optimize=True`` tunes each method's key parameter per-spectrum to
        maximise the IS-Score (slower, higher quality).
        """
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        warnings: List[str] = []

        proposals = self.auto_baseline.evaluate(x, y, optimize=optimize)
        valid = [p for p in proposals if p.score is not None]
        if not valid:
            return AnalysisReport(
                baseline=None,
                warnings=["All baseline methods failed; fall back to manual."],
            )

        best = valid[0]
        alternatives = valid[1 : 1 + keep_alternatives]

        # Low-confidence guard: tell the user when to distrust the automation.
        if best.score.combined < 0.5:
            warnings.append(
                "Low baseline confidence -- review manually "
                "(possible heavy fluorescence or unusual shape)."
            )

        peaks = self.auto_peaks.detect(x, best.y_corr)
        if not peaks:
            warnings.append("No peaks passed the confidence gate.")

        return AnalysisReport(
            baseline=best,
            baseline_alternatives=alternatives,
            peaks=peaks,
            warnings=warnings,
        )

    # Extension point -------------------------------------------------------
    def auto_fit(self, x, y_corr, peaks):  # pragma: no cover
        """Auto peak fitting (Phase 2c).

        TODO: seed FitFunctions (Lorentz/Gauss/Voigt) from detected peak
        position/height/FWHM, fit, and score with quality_metrics.fit_quality.
        Will reuse peakFitting.FitFunctions from the existing app.
        """
        raise NotImplementedError("Auto-fitting is a Phase 2c extension point.")
