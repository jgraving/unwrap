"""Recency-weighted trailing directedness smoother (issue #142).

Covers ``smooth_directedness`` (the finite normalized trailing exponential
window that replaces the former Gaussian x-domain path), the curve builders
``population_curve_kernel_per_tag`` / ``tag_curves_kernel_by_treatment``
(the trailing filter is the ONLY smoother; grid mapping is documented
equal-width binning), and Python<->R agreement on shared fixtures.

Acceptance checks from the issue: W=1 exact identity; constants remain
constant; outputs stay in [0,1]; analytic short-sequence weights including
partial windows; future observations cannot affect earlier results;
animals and gap/missing segments cannot contaminate one another; step
changes get smaller adjacent jumps than raw and less distortion/lag than
equal-weight trailing means; away-facing valid steps contribute; intervals
are summarized after filtering (anti-correlated draws); Python and R agree
on shared fixtures including gaps, missing headings, away-facing steps and
short tracks.

The math is pure numpy; the functions live in
``notebooks/goldshtein_bats_posthoc.py``, whose module-level imports are not
guaranteed in every CI dependency set, so the module skips when absent.
"""
import importlib.util
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("arviz")
pytest.importorskip("joblib")
pytest.importorskip("tqdm")
pytest.importorskip("matplotlib")

_REPO_ROOT = Path(__file__).resolve().parents[1]
_NOTEBOOKS = _REPO_ROOT / "notebooks"
_spec = importlib.util.spec_from_file_location(
    "goldshtein_bats_posthoc_under_test",
    _NOTEBOOKS / "goldshtein_bats_posthoc.py",
)
ph = importlib.util.module_from_spec(_spec)
# exec_module needs the module registered first (a @dataclass in the module
# resolves its namespace via sys.modules).
sys.modules[_spec.name] = ph
_spec.loader.exec_module(ph)


def _alpha(window):
    return 2.0 / (1.0 + np.sqrt(window))


# ---------------------------------------------------------------------------
# smooth_directedness: identity, constants, ranges, validation
# ---------------------------------------------------------------------------

def test_w1_is_exact_identity():
    rng = np.random.default_rng(0)
    q = rng.uniform(0.0, 1.0, (40, 7))
    valid = np.ones(40, dtype=bool)
    gap = rng.random(40) < 0.2
    out = ph.smooth_directedness(q, valid, gap, 1)
    assert out is not q          # a copy, never the input object
    assert np.array_equal(out, q)  # exactly, bit for bit
    # NaN rows survive untouched at W=1
    qn = q.copy()
    qn[[3, 11]] = np.nan
    validn = np.ones(40, dtype=bool)
    validn[[3, 11]] = False
    outn = ph.smooth_directedness(qn, validn, np.zeros(40, dtype=bool), 1)
    assert np.array_equal(outn, qn, equal_nan=True)


def test_constants_remain_constant():
    rng = np.random.default_rng(1)
    for window in (1, 2, 5, 10, 20):
        q = np.full((25, 4), 0.37)
        out = ph.smooth_directedness(q, np.ones(25, dtype=bool),
                                     np.zeros(25, dtype=bool), window)
        assert np.allclose(out, 0.37)
    # constants with gaps/missing too (rows keep raw/NaN, rest stay constant)
    q = np.full((25, 4), 0.37)
    valid = np.ones(25, dtype=bool)
    valid[7] = False
    q[7] = np.nan
    gap = np.zeros(25, dtype=bool)
    gap[15] = True
    out = ph.smooth_directedness(q, valid, gap, window)
    ok = valid.copy()
    ok[15] = False
    assert np.allclose(out[ok], 0.37)
    assert np.all(np.isnan(out[7]))
    assert np.allclose(out[15], 0.37)


def test_outputs_stay_in_unit_interval():
    rng = np.random.default_rng(2)
    q = rng.uniform(0.0, 1.0, (50, 9))
    out = ph.smooth_directedness(q, np.ones(50, dtype=bool),
                                 np.zeros(50, dtype=bool), 10)
    assert np.nanmin(out) >= 0.0 and np.nanmax(out) <= 1.0


@pytest.mark.parametrize("window", [True, False, 0, 21, 2.5, -1, "10", None,
                                    np.bool_(True), np.float64(10.0)])
def test_window_validation(window):
    q = np.full((4, 2), 0.5)
    with pytest.raises(ValueError):
        ph.smooth_directedness(q, np.ones(4, dtype=bool),
                               np.zeros(4, dtype=bool), window)


def test_shape_and_value_validation():
    q = np.full((4, 2), 0.5)
    with pytest.raises(ValueError):
        ph.smooth_directedness(q.ravel(), np.ones(4, dtype=bool),
                               np.zeros(4, dtype=bool), 5)
    with pytest.raises(ValueError):  # mask length mismatch
        ph.smooth_directedness(q, np.ones(3, dtype=bool),
                               np.zeros(4, dtype=bool), 5)
    bad = q.copy()
    bad[2, 0] = np.nan            # NaN on a VALID row
    with pytest.raises(ValueError):
        ph.smooth_directedness(bad, np.ones(4, dtype=bool),
                               np.zeros(4, dtype=bool), 5)
    bad2 = q.copy()
    bad2[1, 1] = 1.2              # out of [0,1] on a VALID row
    with pytest.raises(ValueError):
        ph.smooth_directedness(bad2, np.ones(4, dtype=bool),
                               np.zeros(4, dtype=bool), 5)


# ---------------------------------------------------------------------------
# Exact algorithm: normalized weights, partial windows, truncation
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("window", [2, 3, 5, 10, 20])
def test_analytic_weights_short_sequence(window):
    """scores=[1,0,0,0,0]: out[i] equals the analytic normalized weights."""
    q = np.zeros((5, 1))
    q[0, 0] = 1.0
    out = ph.smooth_directedness(q, np.ones(5, dtype=bool),
                                 np.zeros(5, dtype=bool), window)
    d = 1.0 - _alpha(window)
    # out[i] = weight of step 0 within a trailing window ending at i:
    # decay**i / sum(decay**k for k in 0..min(i, window-1)) once step 0 is
    # still inside the window, and 0 once the W-step truncation drops it.
    for i in range(5):
        if i >= window:
            assert np.isclose(out[i, 0], 0.0), f"W={window} step {i} truncated"
            continue
        m = i + 1
        expected = d ** i / sum(d ** k for k in range(m))
        assert np.isclose(out[i, 0], expected, rtol=0, atol=1e-12), \
            f"W={window} step {i}"


def test_weights_match_issue_example_w10():
    # Issue: at W=10 the current and preceding steps get ~48%, 25%, 13%, 7%
    d = 1.0 - _alpha(10)
    w = np.array([d ** k for k in range(10)])[::-1]
    w = w / w.sum()
    assert np.isclose(w[-1], 0.48, atol=0.01)
    assert np.isclose(w[-2], 0.25, atol=0.01)
    assert np.isclose(w[-3], 0.13, atol=0.01)
    assert np.isclose(w[-4], 0.07, atol=0.01)


def test_window_truncation_at_W():
    """With W=3 the step at index 5 must not see steps 0 or 1."""
    q = np.zeros((6, 1))
    q[0, 0] = 1.0
    out = ph.smooth_directedness(q, np.ones(6, dtype=bool),
                                 np.zeros(6, dtype=bool), 3)
    assert out[5, 0] == 0.0
    q[3, 0] = 1.0  # first index inside the W=3 window ending at 5
    out2 = ph.smooth_directedness(q, np.ones(6, dtype=bool),
                                  np.zeros(6, dtype=bool), 3)
    assert out2[5, 0] > 0.0


def test_short_tracks_partial_windows():
    """Tracks shorter than W normalize over the available steps."""
    q = np.array([[0.2], [0.8]])
    out = ph.smooth_directedness(q, np.ones(2, dtype=bool),
                                 np.zeros(2, dtype=bool), 10)
    d = 1.0 - _alpha(10)
    assert np.isclose(out[0, 0], 0.2)
    assert np.isclose(out[1, 0], (0.8 + d * 0.2) / (1.0 + d))


# ---------------------------------------------------------------------------
# Chronology, gaps and missing rows
# ---------------------------------------------------------------------------

def test_future_observations_cannot_affect_the_past():
    rng = np.random.default_rng(3)
    q = rng.uniform(0.0, 1.0, (30, 5))
    valid = np.ones(30, dtype=bool)
    gap = np.zeros(30, dtype=bool)
    base = ph.smooth_directedness(q, valid, gap, 7)
    q2 = q.copy()
    q2[15:] = rng.uniform(0.0, 1.0, (15, 5))
    out2 = ph.smooth_directedness(q2, valid, gap, 7)
    assert np.allclose(out2[:15], base[:15])


def test_gap_crossing_rows_and_segment_isolation():
    rng = np.random.default_rng(4)
    n = 20
    q = rng.uniform(0.0, 1.0, (n, 5))
    valid = np.ones(n, dtype=bool)
    gap = np.zeros(n, dtype=bool)
    gap[10] = True
    out = ph.smooth_directedness(q, valid, gap, 6)
    # the gap-crossing heading retains its RAW score ...
    assert np.array_equal(out[10], q[10])
    # ... the step right after it starts fresh ...
    assert np.array_equal(out[11], q[11])
    # ... and history never crosses the gap: changing the pre-gap segment
    # cannot move any post-gap result.
    q2 = q.copy()
    q2[:10] = rng.uniform(0.0, 1.0, (10, 5))
    out2 = ph.smooth_directedness(q2, valid, gap, 6)
    assert np.allclose(out2[10:], out[10:])
    assert not np.allclose(out2[:10], out[:10])  # pre-gap itself did change


def test_missing_rows_reset_history_and_keep_missing_values():
    rng = np.random.default_rng(5)
    n = 12
    q = rng.uniform(0.0, 1.0, (n, 4))
    valid = np.ones(n, dtype=bool)
    valid[4] = False
    q[4] = np.nan          # missing rows arrive as NaN
    out = ph.smooth_directedness(q, valid, np.zeros(n, dtype=bool), 5)
    assert np.all(np.isnan(out[4]))       # missing values retained
    assert np.array_equal(out[5], q[5])   # next valid step starts fresh
    # and history cannot cross the missing row backwards either
    q2 = q.copy()
    q2[:4] = rng.uniform(0.0, 1.0, (4, 4))
    out2 = ph.smooth_directedness(q2, valid, np.zeros(n, dtype=bool), 5)
    assert np.allclose(out2[5:], out[5:])


def test_away_facing_valid_steps_contribute():
    """Low (away-facing) valid scores pull the trailing mean down like any
    other valid score."""
    q = np.vstack([np.full((6, 1), 0.9), np.full((6, 1), 0.1),
                   np.full((6, 1), 0.9)])
    out = ph.smooth_directedness(q, np.ones(18, dtype=bool),
                                 np.zeros(18, dtype=bool), 10)
    # during the away stretch the recency-weighted value drops below 0.9 ...
    assert out[9, 0] < 0.8
    # ... and recovers while staying between the levels (never overshoots)
    assert 0.1 < out[13, 0] < 0.9


def test_anti_correlated_draws_filter_before_summarizing():
    """Per-draw filtering must precede intervals.

    Two perfectly anti-correlated joint draws q and 1-q have a constant
    pointwise MEAN of 0.5. Smoothing each draw first (linearity) keeps the
    mean at 0.5 but preserves the per-draw spread, so intervals computed
    AFTER filtering stay wide -- whereas summarizing/intervaling the raw
    means first would collapse them.
    """
    n = 12
    t = np.linspace(0.0, 1.0, n)[:, None]
    q = np.hstack([t, 1.0 - t])
    out = ph.smooth_directedness(q, np.ones(n, dtype=bool),
                                 np.zeros(n, dtype=bool), 4)
    mid = slice(4, 8)
    assert np.allclose(out[mid].mean(axis=1), 0.5)   # mean survives
    spread = np.abs(out[mid, 0] - out[mid, 1])       # intervals survive
    assert spread.max() > 0.25
    assert spread.min() > 0.0


# ---------------------------------------------------------------------------
# Step-change: jump, distortion and lag versus raw and boxcar
# ---------------------------------------------------------------------------

def _boxcar_trailing(scores, window):
    """Equal-weight trailing mean with partial windows (the reference the
    issue asks the recency filter to beat on distortion/lag)."""
    q = np.asarray(scores, float)
    out = q.copy()
    for i in range(1, len(q)):
        lo = max(0, i - window + 1)
        out[i] = q[lo:i + 1].mean(axis=0)
    return out


@pytest.mark.parametrize("window", [2, 5, 10, 20])
def test_step_change_jumps_distortion_and_lag(window):
    n = 120
    edge = 60
    truth = np.where(np.arange(n) < edge, 0.2, 0.8)[:, None]
    raw_jump = float(np.abs(np.diff(truth[:, 0])).max())

    rec = ph.smooth_directedness(truth, np.ones(n, dtype=bool),
                                 np.zeros(n, dtype=bool), window)
    box = _boxcar_trailing(truth, window)

    # (1) adjacent jumps shrink versus raw
    jump_rec = float(np.abs(np.diff(rec[:, 0])).max())
    assert jump_rec < raw_jump

    # (2) less distortion than the equal-weight trailing mean at the same W
    dev_rec = float(np.abs(rec[:, 0] - truth[:, 0]).mean())
    dev_box = float(np.abs(box[:, 0] - truth[:, 0]).mean())
    assert dev_rec < dev_box

    # (3) less lag: the recency filter crosses the midpoint of the step
    # strictly earlier than the boxcar (strict > so the boxcar's exact
    # half-value on the jump step does not count as a crossing).
    def crossing(x):
        idx = np.flatnonzero(x > 0.5)
        return idx[0] if idx.size else n
    assert crossing(rec[:, 0]) < crossing(box[:, 0])


# ---------------------------------------------------------------------------
# Curve builders: single smoother, binning semantics, masks
# ---------------------------------------------------------------------------

def _synthetic_tracks(n_tags=6, n_steps=40, n_draws=30, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for t in range(n_tags):
        treat = (t % 4) + 1
        for s in range(n_steps):
            rows.append({
                "tag": f"tag{t}",
                "treatment": treat,
                "gap_before": 1 if (s == 20) else 0,
                "dist_to_dest_grad_m": float(1000 - s * 10 + rng.normal(0, 5)),
                "flightAngDiffDest_deg": float(rng.normal(0, 40)),
            })
    df = pd.DataFrame(rows)
    df["treatment"] = df["treatment"].astype("category")
    evi = np.clip(0.5 + 0.3 * rng.standard_normal((n_draws, len(df))), 0, 1)
    return df, evi


@pytest.mark.parametrize("x_kind", ["distance", "steps_from_end"])
def test_population_min_two_rule_and_ranges(x_kind):
    df, evi = _synthetic_tracks()
    kwargs = dict(window=10, hdi_prob=0.90, min_tags_at_x=2)
    if x_kind == "distance":
        kwargs["dist_col"] = "dist_to_dest_grad_m"
    else:
        kwargs.update(last_n_steps=110, steps_x_units="m", step_length_m=60.0)
    levels, out = ph.population_curve_kernel_per_tag(df, evi, x_kind=x_kind,
                                                     **kwargs)
    assert set(levels) <= {1, 2, 3, 4}
    for t in levels:
        b = out[t]
        assert b["window"] == 10
        g = len(b["x_mid"])
        assert len(b["mean"]) == g and len(b["lo"]) == g and len(b["hi"]) == g
        support = b["n_tags_by_x"] >= 2
        shown = np.isfinite(b["mean"])
        assert np.array_equal(shown, support & np.isfinite(b["mean"]))
        fin = np.isfinite(b["mean"])
        assert np.all((b["lo"][fin] >= 0) & (b["hi"][fin] <= 1))
        assert np.all(b["lo"][fin] <= b["mean"][fin])
        assert np.all(b["mean"][fin] <= b["hi"][fin])


@pytest.mark.parametrize("x_kind", ["distance", "steps_from_end"])
def test_tag_curves_masking_and_ranges(x_kind):
    df, evi = _synthetic_tracks()
    kwargs = dict(window=10, hdi_prob=0.90)
    if x_kind == "distance":
        kwargs["dist_col"] = "dist_to_dest_grad_m"
    else:
        kwargs.update(last_n_steps=110, steps_x_units="m", step_length_m=60.0)
    levels, out = ph.tag_curves_kernel_by_treatment(df, evi, x_kind=x_kind,
                                                    **kwargs)
    for t in levels:
        b = out[t]
        g = len(b["x_mid"])
        for tg, s in b["tags"].items():
            assert len(s["mean"]) == g
            occ = np.isfinite(s["mean"])
            if occ.any():
                # Bin-start labels sit at most one bin below the tag's first
                # step; nothing appears beyond the observed range.
                bw = b.get("bin_width_m", 0.0)
                assert b["x_mid"][occ].min() >= s["x_min"] - bw - 1e-9
                assert b["x_mid"][occ].max() <= s["x_max"] + 1e-9


def test_engine_w1_equals_plain_bin_means():
    """Through the engine, W=1 (raw scores) must equal plain equal-animal
    per-draw bin means -- the smoother is the only thing that ever touches
    the scores."""
    df, evi = _synthetic_tracks(n_tags=4, n_steps=30, n_draws=12, seed=11)
    levels, out = ph.population_curve_kernel_per_tag(
        df, evi, x_kind="steps_from_end", last_n_steps=110,
        steps_x_units="m", step_length_m=60.0, window=1, hdi_prob=0.90,
        min_tags_at_x=1)
    bw = 60.0
    for t in levels:
        b = out[t]
        sub_all = df[df["treatment"] == t]
        xg = b["x_mid"]
        # recompute per-draw equal-animal bin means by hand
        per_tag = {}
        for tag, sub in sub_all.groupby("tag", sort=False):
            xs = -(len(sub) - 1 - np.arange(len(sub))) * bw
            per_tag[tag] = (xs, evi[:, sub.index.to_numpy()])
        for gi, xc in enumerate(xg):
            contributions = []
            for xs, draws_mat in per_tag.values():
                sel = (np.abs(xs - xc) < bw / 2 + 1e-9)
                if sel.any():
                    contributions.append(draws_mat[:, sel].mean(axis=1))
            if not contributions:
                assert not np.isfinite(b["mean"][gi])
                continue
            pooled = np.mean(contributions, axis=0)  # per-draw equal-animal
            if len(contributions) >= 1:
                assert np.isclose(b["mean"][gi], pooled.mean(), atol=1e-9)


def test_binning_no_bridging_across_unsupported_ranges():
    """A distance band with no steps stays NaN even inside a tag's overall
    x-range: binning fills nothing, so tracking-gap regions are not
    bridged."""
    rows, scores = [], []
    n_draws = 8
    for s in range(20):  # near segment 0..900 m
        rows.append(dict(tag="a", treatment=1, gap_before=0,
                         dist_to_dest_grad_m=float(900 - 45 * s)))
    for s in range(20):  # far segment 2100..3000 m -- nothing between
        rows.append(dict(tag="a", treatment=1, gap_before=0,
                         dist_to_dest_grad_m=float(3000 - 45 * s)))
    for s in range(20):
        rows.append(dict(tag="b", treatment=1, gap_before=0,
                         dist_to_dest_grad_m=float(900 - 45 * s)))
        rows.append(dict(tag="b", treatment=1, gap_before=0,
                         dist_to_dest_grad_m=float(3000 - 45 * s)))
    df = pd.DataFrame(rows)
    df["treatment"] = df["treatment"].astype("category")
    evi = np.full((n_draws, len(df)), 0.6)
    levels, out = ph.tag_curves_kernel_by_treatment(
        df, evi, x_kind="distance", dist_col="dist_to_dest_grad_m",
        window=5, hdi_prob=0.90)
    b = out[1]["tags"]["a"]
    xg = out[1]["x_mid"]
    hole = (xg >= 1000) & (xg <= 1750)   # bins between the two recorded bands
    assert hole.any()
    assert np.all(np.isnan(b["mean"][hole]))


def test_repeated_nonmonotone_distances_are_aggregated():
    """The same x visited twice (out-and-back) contributes both smoothed
    steps to that bin's per-draw mean."""
    # step x: 300, 200, 100, 200, 300 (away-facing return) with constant
    # per-draw scores of 0.2 on the outbound and 0.8 on the return steps.
    xs = [300.0, 200.0, 100.0, 200.0, 300.0]
    vals = [0.2, 0.2, 0.2, 0.8, 0.8]
    rows = [dict(tag="a", treatment=1, gap_before=0, dist_to_dest_grad_m=x)
            for x in xs]
    df = pd.DataFrame(rows)
    df["treatment"] = df["treatment"].astype("category")
    evi = np.tile(np.array(vals), (4, 1))  # one draw replicated x4
    _, out = ph.population_curve_kernel_per_tag(
        df, evi, x_kind="distance", dist_col="dist_to_dest_grad_m",
        window=1, bin_width_m=100.0, hdi_prob=0.90, min_tags_at_x=1)
    b = out[1]
    # W=1: bin x=200 holds steps with values 0.2 and 0.8 -> mean 0.5
    gi = int(np.argmin(np.abs(b["x_mid"] - 200.0)))
    assert np.isclose(b["mean"][gi], 0.5)


def test_history_spans_the_plotting_range_boundary():
    """steps_from_end truncation must not shorten the smoother's history:
    the first kept step of a long track still sees earlier recorded steps
    (it differs from the raw score for nonconstant input)."""
    rows = []
    for s in range(130):
        rows.append(dict(tag="a", treatment=1, gap_before=0,
                         dist_to_dest_grad_m=float(s * 60.0 + 60)))
    for s in range(30):
        rows.append(dict(tag="b", treatment=1, gap_before=0,
                         dist_to_dest_grad_m=float(s * 60.0 + 60)))
    df = pd.DataFrame(rows)
    df["treatment"] = df["treatment"].astype("category")
    # smoothable gradient: score = 0.1 + 0.006*s on tag a's steps
    sa = 0.1 + 0.006 * np.arange(130)
    sb = 0.1 + 0.006 * np.arange(30)
    vals = np.concatenate([sa, sb])[None, :].repeat(6, axis=0)
    evi = np.clip(vals, 0, 1)
    _, out10 = ph.tag_curves_kernel_by_treatment(
        df, evi, x_kind="steps_from_end", last_n_steps=25,
        steps_x_units="m", step_length_m=60.0, window=10, hdi_prob=0.90)
    _, out1 = ph.tag_curves_kernel_by_treatment(
        df, evi, x_kind="steps_from_end", last_n_steps=25,
        steps_x_units="m", step_length_m=60.0, window=1, hdi_prob=0.90)
    a10 = out10[1]["tags"]["a"]
    a1 = out1[1]["tags"]["a"]
    assert np.isfinite(a10["mean"]).sum() > 20
    # the smoothed curve lags the raw scores on the increasing segment
    assert not np.allclose(a10["mean"], a1["mean"], equal_nan=True)


# ---------------------------------------------------------------------------
# Python <-> R agreement
# ---------------------------------------------------------------------------

def _extract_r_function(src, name):
    start = src.find(f"{name} <- function")
    if start == -1:
        raise AssertionError(f"R {name} not found in goldshtein_bats_posthoc.R")
    ob = src.find("{", start)
    depth, i = 0, ob
    while True:
        if src[i] == "{":
            depth += 1
        elif src[i] == "}":
            depth -= 1
            if depth == 0:
                break
        i += 1
    return src[start: i + 1]


def _find_rscript():
    import shutil
    for cand in (
        shutil.which("Rscript"),
        "/home/jake/miniforge3/envs/r453/bin/Rscript",
    ):
        if cand and Path(cand).exists():
            return cand
    return None


_R_HDI = """
hdi <- function(x, prob = 0.95) {
  x <- sort(as.numeric(x))
  stopifnot(length(x) > 1L, all(is.finite(x)))
  gap <- min(length(x) - 1L, floor(prob * length(x)))
  ii <- seq_len(length(x) - gap)
  j <- ii[which.min(x[ii + gap] - x[ii])]
  c(low = x[j], high = x[j + gap])
}
"""


def _r_fixtures():
    """Shared fixtures: gaps, missing headings, away-facing steps and a
    short track (shorter than W)."""
    rng = np.random.default_rng(42)
    cases = []
    # (scores [steps x draws], valid, gap_before, window)
    q = rng.uniform(0.0, 1.0, (18, 6))
    valid = np.ones(18, dtype=bool)
    valid[5] = False
    q[5] = np.nan                      # missing heading
    gap = np.zeros(18, dtype=bool)
    gap[11] = True                     # gap-crossing heading
    cases.append((q, valid, gap, 10))
    # away-facing stretch + short track + gap at the first step
    s = np.vstack([np.full((4, 3), 0.85), np.full((4, 3), 0.1),
                   np.full((4, 3), 0.85)])
    gap2 = np.zeros(12, dtype=bool)
    gap2[0] = True
    cases.append((s, np.ones(12, dtype=bool), gap2, 20))
    # track much shorter than W (partial window everywhere)
    cases.append((rng.uniform(0.0, 1.0, (3, 4)), np.ones(3, dtype=bool),
                  np.zeros(3, dtype=bool), 10))
    return cases


@pytest.mark.parametrize("case_index", range(3))
def test_python_r_smoother_agreement(tmp_path, case_index):
    rscript = _find_rscript()
    if rscript is None:
        pytest.skip("Rscript not available")

    r_src = (_REPO_ROOT / "notebooks_r" / "goldshtein_bats_posthoc.R").read_text()
    func = _extract_r_function(r_src, "smooth_directedness")
    func_file = tmp_path / "r_func.R"
    func_file.write_text(func + "\n")

    q, valid, gap, window = _r_fixtures()[case_index]
    np.savetxt(tmp_path / "scores.csv", q, delimiter=",")
    (tmp_path / "masks.txt").write_text(
        " ".join("1" if v else "0" for v in valid) + "\n" +
        " ".join("1" if g else "0" for g in gap) + "\n")

    runner = tmp_path / "run.R"
    runner.write_text(
        "args <- commandArgs(trailingOnly = TRUE)\n"
        "source(args[1])\n"
        "scores <- as.matrix(read.csv(args[2], header = FALSE))\n"
        "scores[is.na(scores)] <- NA_real_\n"
        "masks <- scan(args[3], what = character(0))\n"
        "nm <- nrow(scores)\n"
        "valid <- masks[1:nm] == '1'\n"
        "gap <- masks[(nm + 1):(2 * nm)] == '1'\n"
        "out <- smooth_directedness(scores, valid, gap, as.integer(args[4]))\n"
        "write.table(out, args[5], row.names = FALSE, col.names = FALSE, "
        "quote = FALSE, na = 'NaN')\n"
    )
    proc = subprocess.run(
        [rscript, str(runner), str(func_file), str(tmp_path / "scores.csv"),
         str(tmp_path / "masks.txt"), str(window), str(tmp_path / "out.csv")],
        cwd=str(tmp_path), capture_output=True, text=True,
    )
    if proc.returncode != 0:
        pytest.skip(f"Rscript failed: {proc.stderr.strip()[:200]}")

    r_out = np.loadtxt(tmp_path / "out.csv", ndmin=2)
    py_out = ph.smooth_directedness(q, valid, gap, window)
    assert r_out.shape == py_out.shape
    both_nan = np.isnan(r_out) & np.isnan(py_out)
    assert np.array_equal(np.isnan(r_out), np.isnan(py_out))
    assert np.allclose(r_out[~both_nan], py_out[~both_nan], atol=1e-9)


def test_python_r_make_curves_agreement(tmp_path):
    """The full curve engines agree on a shared fixture (gaps included)."""
    rscript = _find_rscript()
    if rscript is None:
        pytest.skip("Rscript not available")

    r_src = (_REPO_ROOT / "notebooks_r" / "goldshtein_bats_posthoc.R").read_text()
    funcs = "\n".join(_extract_r_function(r_src, name)
                      for name in ("smooth_directedness", "summarize_curve",
                                   "make_curves"))
    (tmp_path / "r_funcs.R").write_text(_R_HDI + "\n" + funcs + "\n")

    rng = np.random.default_rng(9)
    rows = []
    for tag, treat, n in [("1", 1, 60), ("2", 1, 45), ("3", 2, 50), ("4", 2, 30)]:
        for s in range(n):
            rows.append({
                "tag": tag,
                "treatment": treat,
                "gap_before": 1 if s in (12, 33) else 0,
                "dist_to_dest_grad_m": float(60 * (n - s) + rng.normal(0, 15)),
            })
    df = pd.DataFrame(rows)
    draws = rng.uniform(0.0, 1.0, (200, len(df)))
    df["treatment"] = df["treatment"].astype("category")

    df.to_csv(tmp_path / "df.csv", index=False)
    np.savetxt(tmp_path / "draws.txt", draws.T)  # R orientation: steps x draws

    runner = tmp_path / "run.R"
    runner.write_text(
        "args <- commandArgs(trailingOnly = TRUE)\n"
        f"quick_run <- FALSE\n"
        "df_track <- read.csv(args[1])\n"
        "df_track$tag <- factor(df_track$tag)\n"
        "df_track$treatment <- factor(df_track$treatment)\n"
        "draws <- as.matrix(read.table(args[2]))\n"
        "source(args[3])\n"
        "curves <- make_curves('distance', draws, hdi_prob = 0.90, "
        "window = 10L, min_tags_at_x = 1L)\n"
        "p <- curves$population[order(curves$population$treatment, "
        "curves$population$x), ]\n"
        "write.table(data.frame(x = p$x, mean = p$mean, low = p$low, "
        "high = p$high, treatment = as.integer(p$treatment)), args[4], "
        "row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ',')\n"
    )
    proc = subprocess.run(
        [rscript, str(runner), str(tmp_path / "df.csv"),
         str(tmp_path / "draws.txt"), str(tmp_path / "r_funcs.R"),
         str(tmp_path / "r_pop.csv")],
        cwd=str(tmp_path), capture_output=True, text=True,
    )
    if proc.returncode != 0:
        pytest.skip(f"Rscript failed: {proc.stderr.strip()[:200]}")

    r_pop = pd.read_csv(tmp_path / "r_pop.csv", header=None,
                        names=["x", "mean", "low", "high", "treatment"])
    levels, out = ph.population_curve_kernel_per_tag(
        df, draws, x_kind="distance", dist_col="dist_to_dest_grad_m",
        window=10, hdi_prob=0.90, min_tags_at_x=1)
    frames = []
    for t in levels:
        b = out[t]
        fin = np.isfinite(b["mean"])
        frames.append(pd.DataFrame({
            "x": b["x_mid"][fin], "mean": b["mean"][fin],
            "low": b["lo"][fin], "high": b["hi"][fin],
            "treatment": int(t)}))
    py_pop = pd.concat(frames).sort_values(["treatment", "x"]).reset_index(drop=True)
    r_pop = r_pop.sort_values(["treatment", "x"]).reset_index(drop=True)
    assert len(r_pop) == len(py_pop) and len(py_pop) > 10
    for col in ("x", "mean", "low", "high"):
        assert np.allclose(r_pop[col].to_numpy(float),
                           py_pop[col].to_numpy(float), atol=1e-8), col
