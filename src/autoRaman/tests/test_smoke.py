"""
Smoke test for the AutoRaman scaffold.

Runs the full auto pipeline against the synthetic spectra in test_data/ and
checks the contract holds (baseline chosen, peaks found on a clean spectrum,
nothing crashes on the hard cases). This is a scaffold-level test, not a
correctness/accuracy suite -- accuracy tuning comes with Phase 2a/2b.

Run from PyRamanGUI/src/:
    python -m autoRaman.tests.test_smoke
or with pytest:
    pytest autoRaman/tests/test_smoke.py
"""

from __future__ import annotations

import os
import sys

import numpy as np

# Make 'autoRaman' importable when run directly (python autoRaman/tests/...).
_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = os.path.abspath(os.path.join(_HERE, "..", ".."))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from autoRaman import AutoAnalyzer, AutoBaseline, AutoPeaks  # noqa: E402

# test_data/ lives at the repo root: src/ -> ../test_data
_TEST_DATA = os.path.abspath(os.path.join(_SRC, "..", "test_data"))


def _load(name: str):
    path = os.path.join(_TEST_DATA, name)
    arr = np.loadtxt(path)
    return arr[:, 0], arr[:, 1]


def _run_one(name: str):
    x, y = _load(name)
    report = AutoAnalyzer().analyze(x, y)
    print(f"\n=== {name} ===")
    print(report.summary())
    return report


def main():
    cases = [
        "clean_high_snr.txt",
        "high_fluorescence.txt",
        "low_snr_noisy.txt",
        "cosmic_spike.txt",
        "complex_spectrum.txt",
    ]
    available = [c for c in cases if os.path.exists(os.path.join(_TEST_DATA, c))]
    if not available:
        print(f"!! No test spectra found under {_TEST_DATA}")
        return 1

    failures = []
    for name in available:
        try:
            report = _run_one(name)
            # Contract: clean spectrum must yield a baseline and >=1 peak.
            if name == "clean_high_snr.txt":
                assert report.ok, "clean spectrum: no baseline chosen"
                assert len(report.peaks) >= 1, "clean spectrum: no peaks found"
            # Every successful proposal now carries an IS-Score in [0,1] and the
            # report ranks by it (best first).
            if report.ok:
                top = report.baseline
                assert top.is_score is not None, "no IS-Score on best proposal"
                assert 0.0 <= top.is_score <= 1.0, f"IS-Score out of range: {top.is_score}"
                for alt in report.baseline_alternatives:
                    assert alt.is_score <= top.is_score + 1e-9, "alternatives not ranked"
        except Exception as exc:  # noqa: BLE001
            failures.append((name, repr(exc)))
            print(f"   FAILED: {exc}")

    # The IS-Score-driven ParameterOptimizer must run end-to-end on a real
    # spectrum and never produce a worse best score than the fast path.
    try:
        x, y = _load("high_fluorescence.txt")
        ab = AutoBaseline()
        fast = ab.best(x, y)
        opt = ab.best(x, y, optimize=True)
        print(f"\noptimize check: fast IS={fast.is_score:.3f}  "
              f"optimized IS={opt.is_score:.3f}")
        assert opt.is_score >= fast.is_score - 1e-6, "optimization regressed the score"
    except FileNotFoundError:
        pass
    except Exception as exc:  # noqa: BLE001
        failures.append(("optimize", repr(exc)))
        print(f"   FAILED (optimize): {exc}")

    print("\n" + "=" * 50)
    if failures:
        print(f"SMOKE TEST FAILED ({len(failures)} case(s)):")
        for name, err in failures:
            print(f"  - {name}: {err}")
        return 1
    print(f"SMOKE TEST PASSED ({len(available)} spectra processed).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
