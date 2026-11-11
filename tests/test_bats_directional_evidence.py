"""Directional evidence and candidate-onset helpers (bats twins).

The evidence math is pure numpy, but the functions live in
``notebooks/goldshtein_bats_posthoc.py``, whose module-level imports
(arviz/joblib/tqdm) are not guaranteed in every CI dependency set, so
the whole module skips when those are absent.
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("arviz")
pytest.importorskip("joblib")
pytest.importorskip("tqdm")
pytest.importorskip("matplotlib")

_NOTEBOOKS = Path(__file__).resolve().parents[1] / "notebooks"
_spec = importlib.util.spec_from_file_location(
    "goldshtein_bats_posthoc_under_test",
    _NOTEBOOKS / "goldshtein_bats_posthoc.py",
)
ph = importlib.util.module_from_spec(_spec)
# exec_module needs the module registered first (a @dataclass in the
# module resolves its namespace via sys.modules).
sys.modules[_spec.name] = ph
_spec.loader.exec_module(ph)


def _synthetic_draws(n_obs=30, n_draws=50, seed=0):
    rng = np.random.default_rng(seed)
    mu = rng.uniform(-np.pi, np.pi, size=(n_draws, n_obs))
    kappa = 10.0 ** rng.uniform(-2, 2.5, size=(n_draws, n_obs))
    w = rng.uniform(1e-6, 1 - 1e-6, size=(n_draws, n_obs))
    y = rng.uniform(-np.pi, np.pi, size=n_obs)
    return y, mu, kappa, w


def test_evidence_log_ratio_identity():
    """r/(1-r) divided by the weight odds equals exp(logB) per draw."""
    y, mu, kappa, w = _synthetic_draws()
    logB = ph.evidence_log_ratio(y[None, :], mu, kappa)
    logit_w = np.log(w) - np.log1p(-w)
    r = ph.sigmoid(logit_w + logB)
    lhs = (r / (1 - r)) / (w / (1 - w))
    assert np.max(np.abs(np.log(lhs) - logB)) < 1e-10


def test_evidence_independent_of_weight():
    """E does not use the mixing weight; the raw probability r does."""
    y, mu, kappa, _ = _synthetic_draws()
    e = ph.sigmoid(ph.evidence_log_ratio(y[None, :], mu, kappa))
    logB = ph.evidence_log_ratio(y[None, :], mu, kappa)
    r_a = ph.sigmoid(np.log(0.2) - np.log1p(-0.2) + logB)
    r_b = ph.sigmoid(np.log(0.8) - np.log1p(-0.8) + logB)
    assert np.array_equal(e, ph.sigmoid(logB))
    assert np.max(np.abs(r_a - r_b)) > 1e-3


def test_evidence_half_at_zero_concentration():
    y, mu = _synthetic_draws()[:2]
    kappa0 = np.zeros_like(mu)
    assert np.max(np.abs(ph.evidence_log_ratio(y, mu, kappa0))) < 1e-12


def test_evidence_finite_at_extreme_kappa():
    kx = np.array([1e-8, 1e-3, 1.0, 50.0, 700.0, 701.0, 1e4])
    out = ph.evidence_log_ratio(np.array([0.0]), np.array([0.0]), kx)
    assert np.isfinite(out).all()


def _onset_frame(support, gaps=None, theta=None):
    n = len(support)
    df = pd.DataFrame({
        "tag": ["a"] * n,
        "_step": np.arange(n),
        "evi_logb_gt0": support,
        "gap_before": gaps if gaps is not None else np.zeros(n),
        "HeadingRadians": theta if theta is not None else np.zeros(n),
    })
    return df


def _runs(res):
    return [(int(r["onset_step"]), int(r["end_step"]), int(r["run_len"]),
             bool(r["at_boundary"])) for _, r in res.iterrows()]


def test_candidate_onset_gap_reset_and_boundary():
    # steps: 0..2 supported, 3 not, 4..7 supported (gap before 4), 8-9
    # not, 10..15 supported.
    support = np.array([1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1],
                       dtype=float) * 0.9
    gaps = np.zeros(16)
    gaps[4] = 3.0
    res = ph.candidate_onset_steps(_onset_frame(support, gaps),
                                   min_run_steps=3)
    assert _runs(res) == [(0, 2, 3, True), (4, 7, 4, False), (10, 15, 6, False)]


def test_candidate_onset_missing_support_breaks_runs():
    support = np.array([1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1],
                       dtype=float) * 0.9
    support[4:] = np.nan
    res = ph.candidate_onset_steps(_onset_frame(support), min_run_steps=3)
    assert [(a, b) for a, b, _, _ in _runs(res)] == [(0, 2)]


def test_candidate_onset_gap_splits_run():
    support = np.ones(16) * 0.9
    df = _onset_frame(support)
    df.loc[[3], "gap_before"] = 3.0
    res = ph.candidate_onset_steps(df, min_run_steps=2)
    assert [(a, b) for a, b, _, _ in _runs(res)] == [(0, 2), (3, 15)]


def test_candidate_onset_toward_goal_fraction():
    support = np.ones(16) * 0.9
    theta = np.zeros(16)  # all headings exactly at the destination
    res = ph.candidate_onset_steps(_onset_frame(support, None, theta),
                                   min_run_steps=3)
    assert len(res) == 1
    row = res.iloc[0]
    assert row["toward_goal_frac"] == 1.0
    assert row["mean_support"] == pytest.approx(0.9)
    assert bool(row["at_boundary"]) is True
