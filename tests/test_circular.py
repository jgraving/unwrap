"""Vectorised circular-statistics primitives shared by the research notebooks."""
import numpy as np
import pytest
from scipy.integrate import simpson

from unwrap.circular import (
    A_of_kappa,
    bootstrap_track_condition,
    hierarchical_bootstrap_group,
    circ_dist,
    circ_hdi_centered,
    circ_mean,
    circ_xerr_from_interval,
    circ_mod,
    circular_mean,
    linear_hdi,
    resultant_length_angles,
    sorted_groups,
    stack_posterior_1d,
    von_mises_kde,
    vonmises_kappa_from_R,
    vonmises_pdf,
    wrap_0_2pi,
)


# ---------------------------------------------------------------------------
# Wrapping
# ---------------------------------------------------------------------------

def test_wrap_pi_bounds():
    x = np.array([-4 * np.pi, -np.pi, 0, 0.5, np.pi, 2 * np.pi, 3 * np.pi, 10.0])
    w = circ_mod(x)
    assert (w >= -np.pi).all() and (w <= np.pi).all()
    # +pi maps to the closed -pi end
    assert circ_mod(np.pi) == pytest.approx(-np.pi)


def test_wrap_0_2pi():
    w = wrap_0_2pi(np.array([0, 1, 2 * np.pi, -1, 3 * np.pi]))
    assert (w >= 0).all() and (w <= 2 * np.pi).all()


def test_circ_mod_matches_notebook_identity():
    # The notebooks defined  (x + pi) % (2 pi) - pi  -- identical.
    x = np.linspace(-10 * np.pi, 10 * np.pi, 2000)
    assert np.allclose(circ_mod(x), (x + np.pi) % (2 * np.pi) - np.pi)


def test_circ_dist_symmetry_and_fold():
    a = np.array([0.0, 0.5, 3.0])
    b = np.array([0.1, 1.5, 3.2])
    assert np.allclose(circ_dist(a, b), circ_dist(b, a))
    # a raw gap > pi folds back to the short arc
    assert circ_dist(0.0, 3.2) == pytest.approx(2 * np.pi - 3.2)


# ---------------------------------------------------------------------------
# Moments
# ---------------------------------------------------------------------------

def test_circ_mean_periodic_identity():
    assert circ_mean(np.array([0.1, 0.1])) == pytest.approx(0.1)
    # symmetric pair about 0 -> 0 (mod 2pi)
    v = circ_mean(np.array([0.2, -0.2]))
    assert min(v, 2 * np.pi - v) < 1e-9


def test_circ_mean_winds_around():
    # angles clustered across the -pi/pi seam still average correctly (no seam bias)
    a = np.array([-3.0, -2.9, -2.8, -2.7])
    assert circ_mod(circ_mean(a)) == pytest.approx(-2.85, abs=0.01)


def test_circular_mean_negative_pi_range():
    m = circular_mean(np.array([3.0, 3.05]))
    assert -np.pi <= m <= np.pi
    assert m == pytest.approx(3.02, abs=0.05)


def test_resultant_length_limits():
    assert resultant_length_angles(np.array([0.0, 0.0, 0.01])) == pytest.approx(1.0, abs=1e-3)
    assert resultant_length_angles(np.array([0, 0, np.pi, np.pi])) == pytest.approx(0.0, abs=1e-12)


# ---------------------------------------------------------------------------
# Densities / kappa
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("kappa", [0.1, 1.0, 5.0])
def test_vonmises_pdf_integrates_to_one(kappa):
    t = np.linspace(-np.pi, np.pi, 40001)
    assert simpson(vonmises_pdf(t, 0.4, kappa), x=t) == pytest.approx(1.0, rel=1e-4)


def test_von_mises_kde_positive_and_bimodal_support():
    rng = np.random.default_rng(0)
    data = rng.normal(0.0, 0.1, 200)
    pts = np.linspace(-np.pi, np.pi, 721)
    d = von_mises_kde(data, pts, kappa=8.0)
    assert (d > 0).all()
    assert d.argmin() != d.argmax()  # has a bump and a valley


def test_A_of_kappa_limits_and_stability():
    assert A_of_kappa(np.array([0.0]))[0] == pytest.approx(0.0)
    assert A_of_kappa(np.array([1.0]))[0] == pytest.approx(0.44638996, rel=1e-6)
    # scaled bessel keeps large kappa finite (naive i1/i0 would overflow)
    assert np.isfinite(A_of_kappa(np.array([1e6]))).all()


@pytest.mark.parametrize("R", [0.2, 0.5, 0.9, 0.99])
def test_kappa_from_R_roundtrip(R):
    k = vonmises_kappa_from_R(np.array([R]))[0]
    assert A_of_kappa(np.array([k]))[0] == pytest.approx(R, rel=1e-6)


# ---------------------------------------------------------------------------
# Credible intervals
# ---------------------------------------------------------------------------

def test_linear_hdi_covers_center():
    s = np.linspace(0, 1, 101)
    lo, hi = linear_hdi(s, 0.95)
    assert lo <= 0.5 <= hi


def test_linear_hdi_empty():
    lo, hi = linear_hdi(np.array([]), 0.95)
    assert np.isnan(lo) and np.isnan(hi)


def test_circ_hdi_centered_contiguous():
    rng = np.random.default_rng(1)
    m, lo, hi = circ_hdi_centered(rng.normal(1.0, 0.1, 5000), 0.95)
    assert np.isfinite(m) and np.isfinite(lo) and np.isfinite(hi)
    # interval is narrow around the centre in wrapped difference
    assert abs(circ_mod(m - lo)) < 0.5 and abs(circ_mod(hi - m)) < 0.5


def test_circ_xerr_from_interval():
    m, lo, hi = circ_hdi_centered(np.random.default_rng(2).normal(0.0, 0.2, 5000), 0.95)
    wlo, whi = circ_xerr_from_interval(m, lo, hi)
    assert wlo > 0 and whi > 0
    # half-widths are the wrapped angular distance, in degrees
    assert wlo == pytest.approx(abs(np.degrees(circ_mod(m - lo))))
    assert whi == pytest.approx(abs(np.degrees(circ_mod(hi - m))))


# ---------------------------------------------------------------------------
# Bootstrap / posterior helpers
# ---------------------------------------------------------------------------

def test_bootstrap_track_condition_shape_and_range():
    rng = np.random.default_rng(3)
    y = rng.normal(0.5, 0.2, 50)
    mu, R = bootstrap_track_condition(y, B=200, seed=0)
    assert mu.shape == (200,) and R.shape == (200,)
    assert (R >= 0).all() and (R <= 1).all()
    # near-deterministic input -> high mean alignment
    assert R.mean() > 0.9


def test_bootstrap_track_condition_empty():
    mu, R = bootstrap_track_condition(np.array([]), B=10, seed=0)
    assert mu.size == 0 and R.size == 0


def test_hierarchical_bootstrap_group_keys():
    import pandas as pd

    rng = np.random.default_rng(4)
    rows = []
    for t in range(6):
        for c in (0, 1):
            for _ in range(40):
                rows.append(("T0" + str(t), c, rng.normal(0.3 if c == 0 else -0.3, 0.2)))
    df = pd.DataFrame(rows, columns=["Track", "PostStimulus", "HeadingRadians"])
    out = hierarchical_bootstrap_group(df, B=50, seed=0)
    assert set(out) >= {"mu_pre", "mu_post", "R_pre", "R_post", "tracks"}
    assert len(out["tracks"]) == 6
    assert (np.isnan(out["R_pre"]).sum() == 0)
    assert (out["R_pre"] <= 1).all()


def test_stack_posterior_1d_missing_raises():
    import xarray as xr

    post = xr.DataArray(np.ones((2, 3)), dims=("chain", "draw"), name="a").to_dataset()
    with pytest.raises(KeyError):
        stack_posterior_1d(post, "nope")


def test_sorted_groups_numeric_then_lexical():
    assert sorted_groups(["3", "1", "2", "10"]) == ["1", "2", "3", "10"]
    assert sorted_groups(["b", "a", "c"]) == ["a", "b", "c"]
