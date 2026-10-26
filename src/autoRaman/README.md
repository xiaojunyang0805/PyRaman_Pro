# autoRaman — analysis core of PyRaman Pro

`autoRaman` is the GUI-free analysis package behind PyRaman Pro's **Auto Analyze**,
the live **IS-Score** readout, **Optimize (IS-Score)** and peak detection. It
depends only on numpy, scipy and pybaselines (no PyQt), so it can be used from
scripts and tested headlessly.

Automation **proposes; the user decides.** Every result carries the method name
and parameters used, and `AutoAnalyzer.analyze()` returns a report of proposals
without modifying the spectrum. In the app the proposal is copied into the
normal controls, where it can be fine-tuned or replaced by the manual anchor
baseline.

## IS-Score (baseline-quality metric)

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

On strongly structured backgrounds (e.g. polystyrene on glass) the IS-Score
maximum can favour a baseline that leaves broad humps in; see
[docs/Polystyrene_demo.md](../../docs/Polystyrene_demo.md).

## How the automatic analysis works

1. **Baseline** — `AutoBaseline` runs a panel of methods: AsLS, arPLS, airPLS,
   derpsALS, polynomial and rolling ball. For the penalised-least-squares
   methods λ is chosen per spectrum from a grid (10¹ … 10⁷) by minimising a
   fit cost; the candidates are then **ranked by IS-Score**.
2. **Peaks** — `AutoPeaks` finds peaks on the corrected spectrum with
   `scipy.signal.find_peaks`, gated by an estimated noise floor.
3. **Report** — `AutoAnalyzer` returns an `AnalysisReport` with the best
   baseline, the runners-up and the peaks.

`ParameterOptimizer` (behind the app's **Optimize** button) tunes one method's
main parameter to maximise the IS-Score directly.

## Layout

| File | Role | Status |
|---|---|---|
| `is_score.py` | IS-Score (`is_score`, `ISScoreResult`) | Done |
| `auto_baseline.py` | `AutoBaseline`, `BaselineProposal`, `ParameterOptimizer` | Done |
| `auto_peaks.py` | `AutoPeaks`, `Peak`; CWT / 2nd-derivative detectors and `PeakValidator` | `find_peaks` detector works; the others are not implemented yet |
| `auto_analyzer.py` | `AutoAnalyzer`, `AnalysisReport` | Baseline + peaks work; automatic fitting not implemented yet |
| `quality_metrics.py` | fit cost, noise estimate, peak confidence, fit quality | Used by the modules above |
| `tests/` | IS-Score unit tests and an end-to-end smoke test over `test_data/` | Passing |

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
```

## Tests

```bash
cd src
python -m pytest autoRaman/tests -q
```

## Possible future work

- Peak detection for difficult cases: CWT and 2nd-derivative detectors for
  shoulders and low-SNR spectra, and `PeakValidator` to reject cosmic spikes.
  The current detector errs on the side of finding too many peaks.
- Automatic peak fitting from the detected peaks (`AutoAnalyzer.auto_fit`).
- Calibrating the IS-Score penalty weights against expert-labelled spectra.
