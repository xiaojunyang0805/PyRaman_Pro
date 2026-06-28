# AutoRaman module

Automated baseline correction and peak detection for PyRamanGUI (**Track A**,
roadmap Phase 2). Status: **v0.2.0** — baseline selection is driven by the
**IS-Score** (a published, ground-truth-free baseline-quality metric) and the
per-spectrum `ParameterOptimizer` is implemented; peak detection internals are
still marked `TODO` for accuracy tuning.

### IS-Score (baseline-quality metric)

`is_score.py` is a **clean-room re-implementation** of the Integrity Spectrum
Score from:

> S. Innocente, A. Visentin, S. Andersson-Engels, K. Komolibus, R. Gautam,
> *Automated Baseline Correction Evaluation Score for Raman Spectroscopy*,
> ACS Omega 2026, 11, 25057–25068. DOI: 10.1021/acsomega.5c09870

It scores a baseline correction on `[0, 1]` (higher = better) from only
`(raw, corrected, axis)` — no ground truth, no manual inspection. Only the
*method* (published equations) is reproduced, from the paper; no upstream source
code is used (the authors' repo has no license). Depends only on numpy/scipy.

```python
from autoRaman import is_score
score = is_score(raw, raw - baseline, axis)            # float in [0, 1]
res = is_score(raw, raw - baseline, axis, return_breakdown=True)  # ISScoreResult
print(res.breakdown())
```

`AutoBaseline` ranks candidate methods by IS-Score, and `ParameterOptimizer`
tunes each method's key knob (lambda / poly order / half-window) to maximise it
(`AutoBaseline.evaluate(x, y, optimize=True)`).

## Core principle: automation is additive

Automation **proposes; the user disposes.** Every result carries the underlying
`method` name and `params`, so a proposal can be replayed in PyRamanGUI's
existing **manual** baseline/fit dialogs and hand-tuned. The manual workflow is
never removed — automated processing can be limited in special cases (heavy
fluorescence, overlapping shoulders, unusual line shapes), so the manual path
must always remain available.

`AutoAnalyzer.analyze()` returns a **report of proposals**; it does **not**
mutate the spectrum or auto-apply anything. The GUI layer decides.

## Layout

| File | Role | Status |
|---|---|---|
| `is_score.py` | **IS-Score** baseline-quality metric (`is_score`, `ISScoreResult`) | Done — drives ranking + optimisation |
| `quality_metrics.py` | Objective scoring: `baseline_quality`, `peak_confidence`, `fit_quality` (legacy heuristics; `_robust_noise_rms` reused by IS-Score) | Baseline done; peak/fit partial |
| `auto_baseline.py` | `AutoBaseline` runs a method panel, ranks by **IS-Score**; `BaselineProposal`; `ParameterOptimizer` (IS-Score-maximising search) | Done |
| `auto_peaks.py` | `AutoPeaks` (prominence-gated `find_peaks`); `Peak`; CWT / 2nd-deriv / `PeakValidator` (stubs) | Primary detector works; alternatives TODO |
| `auto_analyzer.py` | `AutoAnalyzer` orchestrator; `AnalysisReport` | Baseline+peak orchestration works; auto-fit TODO |
| `tests/test_smoke.py` | End-to-end smoke test over `test_data/` | Passing |

## Design

- **GUI-free.** Depends only on numpy / scipy / pybaselines (no PyQt), so it is
  headless-testable. The Qt layer (a future "Auto Analyze" action in `Plot.py`)
  imports and calls it.
- **Decoupled from the app core.** The baseline method panel mirrors a subset of
  `analysisMethods.BaselineCorrectionMethods` but calls pybaselines/rampy
  directly, keeping the analysis core independent of the Qt application.

## Usage

```python
import numpy as np
from autoRaman import AutoAnalyzer

x, y = np.loadtxt("spectrum.txt", unpack=True)
report = AutoAnalyzer().analyze(x, y)

print(report.summary())
if report.ok:
    best = report.baseline          # BaselineProposal: .method, .params, .y_corr, .baseline, .score
    for p in report.peaks:          # list[Peak]: .position, .height, .fwhm, .confidence
        ...
    # Hand best.method + best.params to the manual dialog to fine-tune.
```

## Run the smoke test

```bash
cd PyRamanGUI/src
python -m autoRaman.tests.test_smoke
```

## Next steps (roadmap Phase 2)

- **2a** — *done.* `ParameterOptimizer` does a coarse-to-fine search over the
  key knob (`lam`, `poly_order`, `half_window`) maximising the IS-Score, which
  also replaced the home-grown `baseline_quality` heuristic as the ranking
  metric. (Optional future work: Bayesian search; calibrate `PENALTY_WEIGHTS`
  against expert-labelled spectra.)
- **2b** — calibrate `peak_confidence` and implement `PeakValidator` to reject
  cosmic spikes (`test_data/cosmic_spike.txt`); add CWT + 2nd-derivative
  detectors for low-SNR and shoulder peaks. *(Current scaffold over-detects —
  the confidence gate is intentionally loose until this work lands.)*
- **2c** — `AutoAnalyzer.auto_fit`: seed `peakFitting.FitFunctions` from detected
  peaks, fit, score with `fit_quality`.
- **GUI** — wire an "Auto Analyze" action in `Plot.py` that shows the report and
  lets the user accept / switch baseline / drop into the manual dialog.
```
