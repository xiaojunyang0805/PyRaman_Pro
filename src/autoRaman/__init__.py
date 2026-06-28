"""
AutoRaman - Automated Raman spectral analysis for PyRamanGUI (Track A).

This package adds *automated* baseline correction and peak detection on top of
PyRamanGUI. It is deliberately built as additive automation:

    Automation proposes; the user disposes.

Every result this module produces (a chosen baseline method, a set of detected
peaks) carries the underlying method name and parameters, so the user can take
the suggestion straight into PyRamanGUI's existing *manual* dialogs and adjust
it. The manual workflow is never removed or bypassed -- automated processing can
be limited in special cases (heavy fluorescence, overlapping shoulders, unusual
line shapes), so the manual path must always remain available.

Design notes
------------
* This package is GUI-free (numpy / scipy / pybaselines only). It does NOT
  import PyQt5, so it can be unit-tested headlessly against the spectra in
  ``test_data/``. The PyQt layer (Plot.py / a future "Auto Analyze" action)
  calls into these functions and renders the proposals.
* It mirrors the algorithm set in ``analysisMethods.BaselineCorrectionMethods``
  but calls pybaselines / scipy directly to keep the analysis core decoupled
  from the Qt application.

Module map (see doc/PROJECT_ROADMAP.md -> Phase 2)
-------------------------------------------------
    quality_metrics.py  -- scoring functions (baseline, peak, fit quality)
    auto_baseline.py    -- evaluate baseline methods, rank, pick best
    auto_peaks.py       -- detect peaks (scipy find_peaks / CWT / 2nd deriv)
    auto_analyzer.py    -- orchestrator: one-call full auto analysis + report

Status: SCAFFOLD. Interfaces are stable; several internals are marked TODO.
"""

from .auto_baseline import (
    AutoBaseline,
    BaselineProposal,
)
from .auto_peaks import (
    AutoPeaks,
    Peak,
)
from .auto_analyzer import (
    AutoAnalyzer,
    AnalysisReport,
)
from .is_score import (
    is_score,
    ISScoreResult,
)

__all__ = [
    "AutoBaseline",
    "BaselineProposal",
    "AutoPeaks",
    "Peak",
    "AutoAnalyzer",
    "AnalysisReport",
    "is_score",
    "ISScoreResult",
]

__version__ = "0.2.0"
