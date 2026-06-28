"""
auto_baseline.py -- automated baseline method selection for AutoRaman.

Strategy (roadmap Phase 2a):
    1. Run a panel of baseline methods over the spectrum.
    2. Score each result with quality_metrics.baseline_quality().
    3. Rank and return proposals; the caller decides what to apply.

The output is a list of :class:`BaselineProposal`, each carrying the method
name and the exact parameters used. That makes the proposal *portable into the
manual workflow*: a user can open PyRamanGUI's baseline dialog, pick the same
method, and start from the suggested parameters -- then hand-tune for the
special cases where automation falls short.

This module calls pybaselines / rampy directly (no PyQt) so it stays
headless-testable. The method panel intentionally mirrors a subset of
``analysisMethods.BaselineCorrectionMethods`` so behaviour matches the GUI.

Status: SCAFFOLD with a working evaluate/rank loop. Parameter optimisation
(grid / Bayesian search per method) is stubbed -- see ParameterOptimizer TODO.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

from .quality_metrics import baseline_quality, baseline_fit_cost, BaselineScore
from .is_score import is_score as _is_score


# ---------------------------------------------------------------------------
# Result type
# ---------------------------------------------------------------------------

@dataclass
class BaselineProposal:
    """One scored baseline candidate.

    Attributes mirror what the manual dialog needs so the proposal can be
    replayed by hand: ``method`` + ``params`` identify it exactly.
    """
    method: str
    params: dict
    y_corr: np.ndarray = field(repr=False)
    baseline: np.ndarray = field(repr=False)
    score: BaselineScore = None
    error: Optional[str] = None       # populated if the method failed
    fit_cost: Optional[float] = None  # bottom-hugging cost (lower = better)
    is_score: Optional[float] = None  # IS-Score [0,1] (higher = better)

    @property
    def combined(self) -> float:
        return self.score.combined if self.score else -1.0


# ---------------------------------------------------------------------------
# Default method panel (subset of BaselineCorrectionMethods, GUI-free)
# ---------------------------------------------------------------------------

def _asls(x, y, p=0.001, lam=1e7):
    from pybaselines import whittaker
    bl, _ = whittaker.asls(y, lam=lam, p=p)
    return y - bl, bl


def _arpls(x, y, lam=1e7):
    from pybaselines import whittaker
    bl, _ = whittaker.arpls(y, lam=lam)
    return y - bl, bl


def _airpls(x, y, lam=1e7):
    from pybaselines import whittaker
    bl, _ = whittaker.airpls(y, lam=lam)
    return y - bl, bl


def _derpsalsa(x, y, lam=1e6, p=0.01):
    from pybaselines import whittaker
    bl, _ = whittaker.derpsalsa(y, lam=lam, p=p)
    return y - bl, bl


def _poly(x, y, poly_order=3):
    from pybaselines import polynomial
    bl, _ = polynomial.poly(y, x_data=x, poly_order=int(poly_order))
    return y - bl, bl


def _rolling_ball(x, y, half_window=50):
    from pybaselines import morphological
    bl, _ = morphological.rolling_ball(y, half_window=int(half_window))
    return y - bl, bl


# name -> (callable(x, y, **params), default params).
# Keep keys identical to BaselineCorrectionMethods so the GUI can map them.
DEFAULT_METHOD_PANEL: Dict[str, Tuple[Callable, dict]] = {
    "Asymmetric Least Square": (_asls, {"p": 0.001, "lam": 1e7}),
    "Asymmetrically Reweighted Penalized Least Squares": (_arpls, {"lam": 1e7}),
    "Adaptive Iteratively Reweighted Penalized Least Squares": (_airpls, {"lam": 1e7}),
    "Derivative Peak-Screening Asymmetric Least Square": (_derpsalsa, {"lam": 1e6, "p": 0.01}),
    "Polynomial (without regions)": (_poly, {"poly_order": 3}),
    "Rolling Ball": (_rolling_ball, {"half_window": 50}),
}


# Lambda search grid for the Whittaker family. Spans 1e1..1e7: the LOW end
# (1e1-1e3) is needed for spectra with a broad, structured/wavy background where
# a flexible baseline must follow the humps to reveal sharp Raman peaks (e.g.
# polystyrene beads on glass); the high end (1e6-1e7) suits smooth fluorescence.
# The right value is picked per-spectrum by the chosen objective.
_WHITTAKER_LAM_GRID = [10.0 ** e for e in (
    1.0, 1.5, 2.0, 2.5, 3.0, 3.5,
    4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0, 6.5, 7.0)]


def _optimize_lambda(fn, x, y, base_params):
    """Pick the lambda in the grid that minimises baseline_fit_cost.

    Returns (y_corr, baseline, best_lambda, best_cost) or (None, None, None, inf).
    """
    best = (None, None, None, float("inf"))
    for lam in _WHITTAKER_LAM_GRID:
        try:
            params = dict(base_params)
            params["lam"] = lam
            y_corr, baseline = fn(x, y, **params)
            y_corr = np.asarray(y_corr, dtype=float).flatten()
            baseline = np.asarray(baseline, dtype=float).flatten()
            cost = baseline_fit_cost(y, y_corr)
            if cost < best[3]:
                best = (y_corr, baseline, lam, cost)
        except Exception:  # noqa: BLE001 - skip a failing lambda
            continue
    return best


class AutoBaseline:
    """Evaluate a panel of baseline methods and rank them by quality score."""

    def __init__(
        self,
        method_panel: Optional[Dict[str, Tuple[Callable, dict]]] = None,
        weights: Optional[dict] = None,
    ):
        self.method_panel = method_panel or DEFAULT_METHOD_PANEL
        self.weights = weights

    def evaluate(
        self,
        x: np.ndarray,
        y: np.ndarray,
        optimize: bool = False,
    ) -> List[BaselineProposal]:
        """Run every method in the panel and return proposals sorted best-first.

        Proposals are ranked by **IS-Score** (Integrity Spectrum Score, [0,1],
        higher = better) -- a ground-truth-free metric that penalises both over-
        and under-fitting (see ``is_score``). ``baseline_fit_cost`` is retained
        as a tie-breaker.

        With ``optimize=True`` each method's key parameter (lambda / poly order /
        half-window) is tuned per-spectrum to maximise the IS-Score via
        :class:`ParameterOptimizer`; the default fast path keeps the original
        lambda-cost search for the Whittaker family.

        Failures are captured (not raised) so one bad method never sinks the
        whole auto-analysis; failed proposals sort to the bottom.
        """
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        proposals: List[BaselineProposal] = []
        optimizer = ParameterOptimizer() if optimize else None

        for name, (fn, params) in self.method_panel.items():
            try:
                if optimizer is not None:
                    y_corr, baseline, used_params, isc = optimizer.optimize(fn, x, y, params)
                    if y_corr is None:
                        raise RuntimeError("optimisation produced no valid baseline")
                    cost = baseline_fit_cost(y, y_corr)
                elif "lam" in params:
                    # search lambda for the Whittaker family (was a fixed default)
                    y_corr, baseline, lam, cost = _optimize_lambda(fn, x, y, params)
                    if y_corr is None:
                        raise RuntimeError("all lambda candidates failed")
                    used_params = dict(params)
                    used_params["lam"] = lam
                    isc = _is_score(y, y_corr, x)
                else:
                    y_corr, baseline = fn(x, y, **params)
                    y_corr = np.asarray(y_corr, dtype=float).flatten()
                    baseline = np.asarray(baseline, dtype=float).flatten()
                    cost = baseline_fit_cost(y, y_corr)
                    used_params = dict(params)
                    isc = _is_score(y, y_corr, x)

                score = baseline_quality(y, y_corr, baseline, self.weights)
                # The IS-Score is the headline confidence the GUI displays and
                # the value proposals are ranked by.
                score.combined = float(isc)
                proposals.append(
                    BaselineProposal(name, used_params, y_corr, baseline, score,
                                     fit_cost=cost, is_score=float(isc))
                )
            except Exception as exc:  # noqa: BLE001 - report, don't crash
                proposals.append(
                    BaselineProposal(
                        name, dict(params),
                        np.empty(0), np.empty(0),
                        score=None, error=f"{type(exc).__name__}: {exc}",
                        fit_cost=float("inf"), is_score=None,
                    )
                )

        # rank by IS-Score (desc); fit_cost breaks ties; failures sink to bottom.
        proposals.sort(
            key=lambda p: (
                p.is_score if p.is_score is not None else -1.0,
                -(p.fit_cost if p.fit_cost is not None else float("inf")),
            ),
            reverse=True,
        )
        return proposals

    def best(self, x: np.ndarray, y: np.ndarray,
             optimize: bool = False) -> Optional[BaselineProposal]:
        """Return the single highest-scoring proposal (or None if all failed)."""
        proposals = self.evaluate(x, y, optimize=optimize)
        top = proposals[0] if proposals else None
        return top if (top and top.score) else None


# Per-knob search grids for the optimiser. The Whittaker lambda grid is the
# module-level _WHITTAKER_LAM_GRID; these cover the other method families.
_POLY_ORDER_GRID = [2, 3, 4, 5, 6, 7, 8]
_HALF_WINDOW_GRID = [20, 35, 50, 75, 100, 150, 250]


class ParameterOptimizer:
    """Per-method parameter search that maximises the IS-Score.

    For each method the single most influential knob is tuned per-spectrum:

    * Whittaker family (asls / arpls / airpls / derpsalsa) -> ``lam`` (with a
      coarse-to-fine refine around the best decade);
    * polynomial -> ``poly_order``;
    * rolling ball -> ``half_window``.

    The objective is the IS-Score (Innocente et al., ACS Omega 2026) computed on
    the corrected spectrum -- the same metric AutoBaseline ranks by -- so the
    optimiser and the ranker agree on what "good" means. Returns
    ``(y_corr, baseline, used_params, is_score)`` for the best candidate, or
    ``(None, None, base_params, -1.0)`` if every candidate failed.
    """

    def optimize(self, fn: Callable, x: np.ndarray, y: np.ndarray,
                 base_params: dict):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)

        if "lam" in base_params:
            return self._search(fn, x, y, base_params, "lam",
                                _WHITTAKER_LAM_GRID, refine=True)
        if "poly_order" in base_params:
            return self._search(fn, x, y, base_params, "poly_order",
                                _POLY_ORDER_GRID)
        if "half_window" in base_params:
            return self._search(fn, x, y, base_params, "half_window",
                                _HALF_WINDOW_GRID)
        # no tunable knob -> evaluate once at the defaults
        return self._search(fn, x, y, base_params, None, [None])

    # -- internals ----------------------------------------------------------
    def _eval_one(self, fn, x, y, params):
        """Run one parameter set; return (y_corr, baseline, is_score) or None."""
        try:
            y_corr, baseline = fn(x, y, **params)
            y_corr = np.asarray(y_corr, dtype=float).flatten()
            baseline = np.asarray(baseline, dtype=float).flatten()
            return y_corr, baseline, float(_is_score(y, y_corr, x))
        except Exception:  # noqa: BLE001 - skip a failing candidate
            return None

    def _search(self, fn, x, y, base_params, knob, grid, refine=False):
        best = (None, None, dict(base_params), -1.0)  # y_corr, bl, params, score
        for val in grid:
            params = dict(base_params)
            if knob is not None:
                params[knob] = val
            res = self._eval_one(fn, x, y, params)
            if res and res[2] > best[3]:
                best = (res[0], res[1], params, res[2])

        # coarse-to-fine: refine around the best lambda decade (Whittaker only).
        if refine and best[0] is not None and knob == "lam":
            center = best[2][knob]
            fine = [center * f for f in (0.33, 0.5, 0.66, 1.5, 2.0, 3.0)]
            for val in fine:
                params = dict(best[2])
                params[knob] = val
                res = self._eval_one(fn, x, y, params)
                if res and res[2] > best[3]:
                    best = (res[0], res[1], params, res[2])

        return best
