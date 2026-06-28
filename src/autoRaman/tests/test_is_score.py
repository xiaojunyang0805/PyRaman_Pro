"""
Tests for the IS-Score baseline-quality metric (autoRaman.is_score).

These are *relative* / ordering tests, not byte-for-byte reproduction of the
authors' published values (we cannot, without their blood-cell data set and the
exact thresholds). The contract verified here is what ranking and parameter
optimisation rely on:

    * the score is always a defined float in [0, 1];
    * a good baseline scores higher than a deliberate over-fit or under-fit;
    * degenerate inputs never raise;
    * on real spectra, a mid-strength lambda beats both extremes.

Run from PyRamanGUI/src/:
    python -m autoRaman.tests.test_is_score
or with pytest:
    pytest autoRaman/tests/test_is_score.py
"""

from __future__ import annotations

import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = os.path.abspath(os.path.join(_HERE, "..", ".."))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from autoRaman.is_score import is_score, ISScoreResult  # noqa: E402

_TEST_DATA = os.path.abspath(os.path.join(_SRC, "..", "test_data"))


# ---------------------------------------------------------------------------
# Synthetic fixtures
# ---------------------------------------------------------------------------

def _gauss(x, c, a, w):
    return a * np.exp(-0.5 * ((x - c) / w) ** 2)


def _make_case():
    """A fluorescence-laden spectrum with a known background and peaks, plus a
    good / over-fit / under-fit baseline correction of it."""
    rng = np.random.default_rng(0)
    x = np.linspace(200, 2000, 1000)
    bg = 800 * np.exp(-(x - 200) / 900) + 200
    peaks = (_gauss(x, 600, 400, 12) + _gauss(x, 950, 600, 10)
             + _gauss(x, 1350, 300, 15) + _gauss(x, 1600, 500, 11))
    raw = bg + peaks + rng.normal(0, 8, x.size)

    good = raw - bg                         # remove the true background
    overfit = raw - (bg + 0.6 * peaks)      # baseline eats 60% of the peaks
    underfit_bl = np.linspace(bg[0] * 0.4, bg[-1] * 0.4, x.size)
    underfit = raw - underfit_bl            # baseline far below -> leftover bg
    return x, raw, good, overfit, underfit


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_score_in_range():
    x, raw, good, overfit, underfit = _make_case()
    for corr in (good, overfit, underfit):
        s = is_score(raw, corr, x)
        assert isinstance(s, float)
        assert 0.0 <= s <= 1.0


def test_good_beats_overfit_and_underfit():
    x, raw, good, overfit, underfit = _make_case()
    s_good = is_score(raw, good, x)
    s_over = is_score(raw, overfit, x)
    s_under = is_score(raw, underfit, x)
    assert s_good > s_over, (s_good, s_over)
    assert s_good > s_under, (s_good, s_under)
    # a good fit should land clearly in the upper range
    assert s_good > 0.8, s_good


def test_breakdown_shape():
    x, raw, good, *_ = _make_case()
    res = is_score(raw, good, x, return_breakdown=True)
    assert isinstance(res, ISScoreResult)
    assert 0.0 <= res.score <= 1.0
    assert res.n_peaks >= 1
    # every block penalty is a finite, non-negative number
    for p in (res.single_band, res.band_region, res.auc,
              res.mean_ratio, res.intensity):
        assert np.isfinite(p) and p >= 0.0
    assert "IS-Score" in res.breakdown()


def test_degenerate_inputs_do_not_raise():
    # flat spectrum
    x = np.linspace(0, 100, 50)
    assert 0.0 <= is_score(np.ones_like(x), np.zeros_like(x), x) <= 1.0
    # too short
    assert 0.0 <= is_score(np.array([1.0, 2.0]), np.array([0.0, 0.0]),
                           np.array([0.0, 1.0])) <= 1.0
    # length mismatch
    assert 0.0 <= is_score(np.ones(10), np.ones(8), np.arange(10)) <= 1.0


def test_determinism():
    x, raw, good, *_ = _make_case()
    assert is_score(raw, good, x) == is_score(raw, good, x)


def test_real_data_lambda_ordering():
    """On a fluorescent spectrum, a mid-strength arPLS lambda should beat both
    the over-fit (tiny lambda) and under-fit (huge lambda) extremes."""
    fn = os.path.join(_TEST_DATA, "high_fluorescence.txt")
    if not os.path.exists(fn):
        return  # test data not present -> skip silently
    try:
        from pybaselines import whittaker
    except Exception:  # pragma: no cover - pybaselines optional in some envs
        return
    a = np.loadtxt(fn)
    x, y = a[:, 0], a[:, 1]
    scores = {}
    for lam in (1e2, 1e6, 1e8):
        bl, _ = whittaker.arpls(y, lam=lam)
        scores[lam] = is_score(y, y - bl, x)
    assert scores[1e6] > scores[1e2], scores
    assert scores[1e6] > scores[1e8], scores


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failures = []
    for t in tests:
        try:
            t()
            print(f"  PASS  {t.__name__}")
        except Exception as exc:  # noqa: BLE001
            failures.append((t.__name__, repr(exc)))
            print(f"  FAIL  {t.__name__}: {exc}")
    print("=" * 50)
    if failures:
        print(f"IS-SCORE TESTS FAILED ({len(failures)}/{len(tests)})")
        return 1
    print(f"IS-SCORE TESTS PASSED ({len(tests)} tests)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
