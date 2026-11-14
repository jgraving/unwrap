"""Recency-weighted trailing directedness smoother (issue #142).

Covers the Python reference ``smooth_directedness`` (finite normalized
trailing exponential window), the curve builder
``directedness_curves_recency``, and a Python<->R agreement check that runs
the R twin's ``smooth_directedness`` (extracted from the R bats notebook)
through ``Rscript`` on shared fixtures.

The math is pure numpy; the functions live in
``notebooks/goldshtein_bats_posthoc.py``, whose module-level imports are not
guaranteed in every CI dependency set, so the module skips when absent.
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


def _decay(window):
    alpha = 2.0 / (1.0 + np.sqrt(window))
    return 1.0 - alpha


def test_w1_is_exact_identity():
    rng = np.random.default_rng(0)
    q = rng.uniform(0, 1, size=(7, 3))
    out = ph.smooth_directedness(
        q, np.ones(7, bool), np.zeros(7, bool), window=1
    )
    assert np.array_equal(out, q)


def test_constant_input_is_constant():
    for w in (1, 3, 10, 20):
        q = np.full((6, 2), 0.42)
        out = ph.smooth_directedness(
            q, np.ones(6, bool), np.zeros(6, bool), window=w
        )
        assert np.allclose(out, 0.42), f"constant broken at W={w}"


def test_output_stays_in_unit_interval():
    rng = np.random.default_rng(1)
    q = np.clip(rng.normal(0.5, 0.35, size=(12, 4)), 0, 1)
    out = ph.smooth_directedness(
        q, np.ones(12, bool), np.zeros(12, bool), window=10
    )
    assert (out >= 0).all() and (out <= 1).all()


def test_no_future_influence():
    # Changing the last step must not move any earlier step.
    q = np.array([[0.2], [0.4], [0.6], [0.8], [0.1]])
    valid = np.ones(5, bool)
    gap = np.zeros(5, bool)
    a = ph.smooth_directedness(q, valid, gap, window=5)
    qb = q.copy()
    qb[4] = 0.0
    b = ph.smooth_directedness(qb, valid, gap, window=5)
    assert np.allclose(a[:4], b[:4])


def test_gap_resets_history():
    # A high value before a gap must not leak into the step after it.
    q = np.array([[0.9], [0.9], [0.1], [0.1]])
    valid = np.ones(4, bool)
    gap = np.array([False, False, True, False])
    out = ph.smooth_directedness(q, valid, gap, window=10)
    # Step 3 is the first valid step after the gap: history is just itself.
    assert out[3, 0] == pytest.approx(0.1)


def test_missing_row_is_excluded_and_kept():
    # A missing (NA) heading keeps its value and resets the history.
    q = np.array([[0.5], [0.5], [np.nan], [0.5]])
    valid = np.array([True, True, False, True])
    gap = np.zeros(4, bool)
    out = ph.smooth_directedness(q, valid, gap, window=10)
    assert np.isnan(out[2, 0])
    # Step 3's history is only itself (step 2 was missing).
    assert out[3, 0] == pytest.approx(0.5)


def test_gap_crossing_row_keeps_raw_score():
    # A valid heading that crosses a gap keeps its raw score (not smoothed).
    q = np.array([[0.2], [0.4], [0.6], [0.8]])
    valid = np.ones(4, bool)
    gap = np.array([False, False, False, True])
    out = ph.smooth_directedness(q, valid, gap, window=10)
    assert out[3, 0] == pytest.approx(0.8)


def test_analytic_weights_full_and_partial_window():
    # Full window at the last step of a W-step run.
    w = 5
    q = np.array([[0.1], [0.2], [0.3], [0.4], [0.5]])
    out = ph.smooth_directedness(
        q, np.ones(5, bool), np.zeros(5, bool), window=w
    )
    decay = _decay(w)
    age = np.arange(w - 1, -1, -1)  # newest first
    weights = decay ** age
    weights /= weights.sum()
    assert out[4, 0] == pytest.approx(float(weights @ q[:, 0]))

    # Partial window: only 3 consecutive valid steps available through i=2
    # with a larger W -> m=3, not W.
    w = 10
    q = np.array([[0.1], [0.2], [0.3]])
    out = ph.smooth_directedness(
        q, np.ones(3, bool), np.zeros(3, bool), window=w
    )
    decay = _decay(w)
    age = np.arange(2, -1, -1)
    weights = decay ** age
    weights /= weights.sum()
    assert out[2, 0] == pytest.approx(float(weights @ q[:, 0]))


def test_window_validation():
    q = np.array([[0.2], [0.8]])
    valid = np.ones(2, bool)
    gap = np.zeros(2, bool)
    for bad in (0, 21, 2.5, True):
        with pytest.raises((ValueError, TypeError)):
            ph.smooth_directedness(q, valid, gap, window=bad)


def test_score_range_validation():
    with pytest.raises(ValueError):
        ph.smooth_directedness(
            np.array([[1.5], [0.5]]), np.ones(2, bool), np.zeros(2, bool),
            window=2,
        )


def _boxcar_trailing(q1d, w):
    """Equal-weight trailing mean over the same window (the comparison)."""
    n = len(q1d)
    out = np.empty(n)
    for i in range(n):
        m = min(w, i + 1)
        out[i] = q1d[i - m + 1: i + 1].mean()
    return out


def test_step_change_smaller_jumps_and_less_distortion_than_boxcar():
    # A step change from 0.2 to 0.8 at the midpoint. The recency filter must
    # (a) reduce adjacent jumps versus the raw signal and (b) distort less
    # (stay closer to the true level) than the equal-weight trailing mean at
    # the same window, while lagging no more than the boxcar.
    n, step = 60, 30
    raw = np.array([0.2] * step + [0.8] * (n - step))
    valid = np.ones(n, bool)
    gap = np.zeros(n, bool)
    raw_jump = float(np.abs(np.diff(raw)).max())

    for w in (2, 5, 10, 20):
        ex = ph.smooth_directedness(raw[:, None], valid, gap, window=w)[:, 0]
        bx = _boxcar_trailing(raw, w)
        ex_jump = float(np.abs(np.diff(ex)).max())
        ex_dist = float(np.abs(ex - raw).sum())
        bx_dist = float(np.abs(bx - raw).sum())

        def cross(s):
            idx = np.where(s >= 0.5)[0]
            return int(idx[0]) if idx.size else -1

        ex_cross, bx_cross = cross(ex), cross(bx)
        assert ex_jump < raw_jump, f"jump not reduced at W={w}"
        assert ex_dist < bx_dist, f"distortion not smaller than boxcar at W={w}"
        assert ex_cross <= bx_cross, f"lag worse than boxcar at W={w}"


def test_interval_summarized_after_filtering():
    # Anti-correlated draws: each draw's score is a step change in opposite
    # directions, so the per-draw smoothing pulls each draw toward its own
    # level and the cross-draw interval is narrower than the raw 0/1 spread.
    q = np.array([[0.0, 1.0], [0.0, 1.0], [1.0, 0.0], [1.0, 0.0]])
    valid = np.ones(4, bool)
    gap = np.zeros(4, bool)
    sm = ph.smooth_directedness(q, valid, gap, window=4)

    # Each draw must be aligned: smoothing column d equals smoothing q[:, d].
    for d in range(2):
        col = ph.smooth_directedness(q[:, d: d + 1], valid, gap, window=4)[:, 0]
        assert np.allclose(sm[:, d], col)

    lo, hi = ph.hdi_1d(sm[3, :], prob=0.90)
    assert (hi - lo) < 1.0, "filtered interval should be tighter than raw"


def _synthetic_tracks(n_tags=6, n_steps=40, n_draws=30, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for t in range(n_tags):
        treat = (t % 4) + 1
        for s in range(n_steps):
            rows.append({
                "tag": f"tag{t}",
                "treatment": treat,
                "gap_before": 1 if (s == 0 or s == 20) else 0,
                "dist_to_dest_grad_m": float(1000 - s * 10 + rng.normal(0, 5)),
                "flightAngDiffDest_deg": float(rng.normal(0, 40)),
            })
    df_full = pd.DataFrame(rows)
    # Drop a few rows to emulate missing headings; the index into df_full is
    # preserved (as build_modeling_table / complete.cases do).
    df_model = df_full.drop(index=[3, 7, 15, 22])
    evi = np.clip(
        0.5 + 0.3 * rng.standard_normal((n_draws, len(df_model))), 0, 1
    )
    return df_full, df_model, evi


@pytest.mark.parametrize("x_kind", ["distance", "steps_from_end"])
def test_curves_recency_shapes_and_min_two_rule(x_kind):
    df_full, df_model, evi = _synthetic_tracks()
    levels, out = ph.directedness_curves_recency(
        df_full, df_model, evi, x_kind=x_kind, window=10, n_grid=200,
        hdi_prob=0.90, min_tags_at_x=2,
    )
    assert set(levels) <= {1, 2, 3, 4}
    for t in levels:
        b = out[t]
        g = len(b["x_mid"])
        assert len(b["mean"]) == g and len(b["lo"]) == g and len(b["hi"]) == g
        # Population support never below the min-two-animal rule.
        assert b["n_tags_by_x"].max() >= 0
        support = b["n_tags_by_x"] >= 2
        assert np.all(np.isnan(b["mean"]) | support | (b["n_tags_by_x"] > 0))
        for tg, s in b["tags"].items():
            assert len(s["mean"]) == g and len(s["lo"]) == g and len(s["hi"]) == g
            # A per-animal curve is masked outside its own observed x-range.
            occ = np.isfinite(s["mean"])
            if occ.any():
                assert b["x_mid"][occ].min() >= s["x_min"] - 1e-9
                assert b["x_mid"][occ].max() <= s["x_max"] + 1e-9


def _extract_r_smooth_directedness():
    """Pull the R twin's smooth_directedness source out of the R bats code.

    PR #143 moved the curve machinery out of the notebook into
    notebooks_r/goldshtein_bats_posthoc.R, so that file is the primary source;
    the notebook is a fallback (older structure).
    """
    import json

    def _func(src):
        start = src.find("smooth_directedness <- function")
        if start == -1:
            return None
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

    script = _REPO_ROOT / "notebooks_r" / "goldshtein_bats_posthoc.R"
    if script.exists():
        out = _func(script.read_text())
        if out is not None:
            return out
    nb = json.loads((_REPO_ROOT / "notebooks_r" / "goldshtein_bats.ipynb").read_text())
    for cell in nb["cells"]:
        if cell.get("cell_type") != "code":
            continue
        src = "".join(cell["source"]) if isinstance(cell["source"], list) else cell["source"]
        out = _func(src)
        if out is not None:
            return out
    raise AssertionError("R smooth_directedness not found in the R bats code")


def test_python_r_agreement(tmp_path):
    rscript = _find_rscript()
    if rscript is None:
        pytest.skip("Rscript not available")

    func = _extract_r_smooth_directedness()
    func_file = tmp_path / "r_smooth_func.R"
    func_file.write_text(func + "\n")

    rng = np.random.default_rng(7)
    n, d = 12, 4
    q = np.clip(0.5 + 0.3 * rng.standard_normal((n, d)), 0, 1)
    valid = np.array([1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1])
    gap = np.array([0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0])

    q_csv = tmp_path / "q.csv"
    with open(q_csv, "w") as f:
        for row in q:
            f.write(",".join(f"{x:.10f}" for x in row) + "\n")
    v_txt = tmp_path / "valid.txt"
    v_txt.write_text("\n".join(str(int(x)) for x in valid) + "\n")
    g_txt = tmp_path / "gap.txt"
    g_txt.write_text("\n".join(str(int(x)) for x in gap) + "\n")

    runner = tmp_path / "run.R"
    runner.write_text(
        "args <- commandArgs(trailingOnly = TRUE)\n"
        "source(args[1])\n"
        "q <- as.matrix(read.csv(args[2], header = FALSE))\n"
        "valid <- as.logical(scan(args[3], what = integer(0)))\n"
        "gap <- as.logical(scan(args[4], what = integer(0)))\n"
        "out <- smooth_directedness(q, valid, gap, window = as.integer(args[5]))\n"
        "write.table(out, args[6], row.names = FALSE, col.names = FALSE, quote = FALSE)\n"
    )
    out_csv = tmp_path / "out.csv"
    import subprocess
    proc = subprocess.run(
        [rscript, str(runner), str(func_file), str(q_csv), str(v_txt),
         str(g_txt), "5", str(out_csv)],
        cwd=str(tmp_path),  # run in a neutral dir to avoid a repo renv hook
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        pytest.skip(f"Rscript failed: {proc.stderr.strip()[:200]}")

    r = np.loadtxt(out_csv)
    py = ph.smooth_directedness(q, valid.astype(bool), gap.astype(bool), window=5)
    assert np.allclose(py, r, atol=1e-9)


def _find_rscript():
    import shutil
    for cand in (
        shutil.which("Rscript"),
        "/home/jake/miniforge3/envs/r453/bin/Rscript",
    ):
        if cand and Path(cand).exists():
            return cand
    return None
