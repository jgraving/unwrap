"""
Vectorised circular-statistics utilities shared by the research notebooks.

Every helper here operates on angles in radians unless a name says otherwise
(the ``*_deg`` forms return degrees). Two wrapping conventions appear in the
original analyses and both are kept, so notebooks can move verbatim onto
these definitions without changing behaviour:

* :func:`wrap_pi` / :func:`circ_mod`  -- wrap to ``(-pi, pi]``
* :func:`wrap_0_2pi`                  -- wrap to ``[0, 2pi)``

Headings that the notebooks report as "radians" are already wrapped by
tracking software to ``(-pi, pi)``; the helpers below accept either range.

The Bessel-ratio helpers build on :mod:`unwrap.bessel`, which uses the
exponentially-scaled ``i0e`` / ``i1e`` so that
``A_of_kappa(large kappa)`` does not overflow the way the naive
``scipy.special.i1/kappa`` / ``i0`` implementations do.
"""

from __future__ import annotations

import numpy as np
from scipy import stats as scs

from unwrap.bessel import bessel_ratio

__all__ = [
    "wrap_pi",
    "circ_mod",
    "wrap_0_2pi",
    "circ_mean",
    "circular_mean",
    "resultant_length_angles",
    "circ_dist",
    "circ_pdf",
    "vonmises_pdf",
    "von_mises_kde",
    "A_of_kappa",
    "vonmises_kappa_from_R",
    "linear_hdi",
    "linear_hdi_mean",
    "circ_hdi_centered",
    "circ_hdi",
    "circ_xerr_from_interval",
    "circular_summary",
    "stack_posterior_1d",
    "bootstrap_track_condition",
    "hierarchical_bootstrap_group",
    "sorted_groups",
]

TWO_PI = 2.0 * np.pi


# ---------------------------------------------------------------------------
# Wrapping
# ---------------------------------------------------------------------------

def wrap_pi(x):
    """Wrap angles to ``(-pi, pi]`` (alias of :func:`circ_mod`)."""
    return (np.asarray(x) + np.pi) % TWO_PI - np.pi


def circ_mod(x):
    """Wrap angles to ``(-pi, pi]``. Vectorised, array-safe."""
    return wrap_pi(x)


def wrap_0_2pi(x):
    """Wrap angles to ``[0, 2pi)``."""
    return np.asarray(x) % TWO_PI


def circ_dist(a, b, axis=None):
    """Circular distance ``|a - b|`` wrapped to ``[0, pi]``."""
    d = np.abs(circ_mod(np.asarray(a, float) - np.asarray(b, float)))
    return d


# ---------------------------------------------------------------------------
# Moments
# ---------------------------------------------------------------------------

def circ_mean(a, axis=None):
    """Circular mean in ``[0, 2pi)`` (vector mean of unit phasors).

    Note the original zebra-analysis convention: result in ``[0, 2pi)``,
    so callers that mix it with ``(-pi, pi]`` angles wrap explicitly.
    """
    a = np.asarray(a, float)
    s = np.sin(a).mean(axis=axis)
    c = np.cos(a).mean(axis=axis)
    return np.arctan2(s, c) % TWO_PI


def circular_mean(theta):
    """Circular mean in ``(-pi, pi]`` (unit-phasor convention).

    The bats-analysis flavour of the mean; differs from :func:`circ_mean`
    only in the output range, which is what the bats notebooks plot.
    """
    z = np.exp(1j * np.asarray(theta, float))
    return float(circ_mod(np.angle(z.mean())))


def resultant_length_angles(a, axis=-1):
    """Mean resultant length ``R`` in ``[0, 1]`` (alignment strength)."""
    a = np.asarray(a, float)
    return np.sqrt(np.mean(np.cos(a), axis=axis) ** 2 + np.mean(np.sin(a), axis=axis) ** 2)


# ---------------------------------------------------------------------------
# Densities / kappa<->R
# ---------------------------------------------------------------------------

def vonmises_pdf(theta, mu, kappa):
    """Von Mises pdf ``exp(kappa cos(theta-mu)) / (2pi I0(kappa))``."""
    theta = np.asarray(theta, float)
    return np.exp(kappa * np.cos(theta - mu)) / (TWO_PI * np.i0(kappa))


def von_mises_kde(data: np.ndarray, points: np.ndarray, kappa: float) -> np.ndarray:
    """Von Mises kernel density estimate (density per radian).

    Each ``data`` point contributes a von Mises kernel of concentration
    ``kappa``; the returned array has the shape of ``points``.
    """
    data = circ_mod(np.asarray(data, float))
    points = np.asarray(points, float)
    c = 1.0 / (TWO_PI * np.i0(kappa))
    cos_diffs = np.cos(points[:, None] - data[None, :])
    return c * np.exp(kappa * cos_diffs).mean(axis=1)


def A_of_kappa(kappa):
    """Mean resultant length ``A(kappa) = I1(kappa)/I0(kappa)`` (vectorised).

    Wraps :func:`unwrap.bessel.bessel_ratio`, which uses the scaled
    ``i0e``/``i1e`` so large ``kappa`` stays finite (the naive
    ``scipy.special.i1/i0`` ratios saturate or overflow past ~700).
    """
    return np.vectorize(bessel_ratio, otypes=[float])(np.asarray(kappa, float))


def vonmises_kappa_from_R(R):
    """Invert :func:`A_of_kappa`: concentration for a given mean resultant length."""
    from unwrap.bessel import inverse_bessel_ratio

    return np.vectorize(inverse_bessel_ratio, otypes=[float])(np.asarray(R, float))


# ---------------------------------------------------------------------------
# Credible intervals
# ---------------------------------------------------------------------------

def linear_hdi(samples, prob=0.95):
    """Shortest Highest Density Interval for LINEAR samples.

    Returns ``(lo, hi)``; ``(nan, nan)`` for empty input.
    """
    s = np.sort(np.asarray(samples, float))
    n = s.size
    if n == 0:
        return np.nan, np.nan
    k = int(np.ceil(prob * n))
    spans = s[k - 1:] - s[:n - k + 1]
    j = int(np.argmin(spans))
    return float(s[j]), float(s[j + k - 1])


def linear_hdi_mean(x, hdi_prob=0.95):
    """``(mean, hdi_lo, hdi_hi)`` for linear samples via ``arviz.hdi``."""
    import arviz as az

    x = np.asarray(x, float)
    h = az.hdi(x, prob=hdi_prob, skipna=True)
    return float(np.nanmean(x)), float(h[0]), float(h[1])


def circ_hdi(samples, prob=0.95):
    """Circular HDI ``(lo, hi)`` in ``(-pi, pi]`` (arviz's circular hdi).

    arviz 1.x removed the top-level ``circular_hdi``; the same computation
    is ``az.hdi(..., circular=True)``.
    """
    import arviz as az

    h = az.hdi(np.asarray(samples, float), prob=prob, circular=True)
    return float(h[0]), float(h[1])


def circ_hdi_centered(samples, prob=0.95):
    """Circular HDI as ``circ_mean + linear_HDI(wrapped deviations)``.

    The zebra-analysis flavour: centre on the circular mean and take a
    shortest HDI on the wrapped deviations, unwrapping the endpoints
    relative to the centre so the interval is contiguous in angle.
    """
    a = np.asarray(samples, float)
    m = circ_mean(a)
    d = circ_mod(a - m)
    lo_d, hi_d = linear_hdi(d, prob)
    lo = float(circ_mod(m + lo_d))
    hi = float(circ_mod(m + hi_d))
    return m, lo, hi


def circ_xerr_from_interval(m, lo, hi):
    """Half-widths ``(w_lo, w_hi)`` in DEGREES of a wrapped (lo, hi) about mean m.

    For ``ax.errorbar`` on an axis in degrees: the lower width extends from the
    point to ``lo`` (clockwise), the upper width to ``hi`` (counter-clockwise),
    each measured as the minimal wrapped angular distance.
    """
    w_lo = abs(float(circ_mod(np.asarray(m) - np.asarray(lo))))
    w_hi = abs(float(circ_mod(np.asarray(hi) - np.asarray(m))))
    return np.degrees(w_lo), np.degrees(w_hi)


def _circular_mcse(samples, batches=5):
    """Batched Monte Carlo SE of a circular mean (arviz 0.23 ``_mc_error``).

    The trace is split into ``batches`` equally-sized blocks; the SE is the
    circular standard deviation of the block means divided by
    ``sqrt(batches)`` (identical to arviz 0.23's non-numba path).
    """
    a = np.asarray(samples, float).ravel()
    if a.size == 0:
        return np.nan
    if batches == 1:
        return float(scs.circstd(a, high=np.pi, low=-np.pi)) / np.sqrt(a.size)
    b = np.resize(a, (batches, int(a.size / batches)))
    means = scs.circmean(b, high=np.pi, low=-np.pi, axis=1)
    return float(scs.circstd(means, high=np.pi, low=-np.pi)) / np.sqrt(batches)


def circular_summary(idata, var_names=None, circ_var_names=(), *, hdi_prob=0.94,
                     skipna=False, group="posterior"):
    """``az.summary`` with the circular columns arviz 1.x dropped.

    arviz 0.23's ``az.summary(..., circ_var_names=[...])`` replaced the
    first five statistics of the named variables - mean, sd, HDI
    lower/upper, ``mcse_mean`` - with their circular equivalents (scipy
    ``circmean``/``circstd`` on ``(-pi, pi]``, the circular HDI, and the
    batched circular MCSE) while ``mcse_sd``, the effective sample sizes
    and ``r_hat`` stayed linear. arviz 1.x removed both the kwarg and the
    top-level ``circular_hdi``. This reproduces that behaviour exactly,
    including the 0.23 frame: the base comes from
    ``az.summary(kind="all", ci_kind="hdi")`` (``hdi_prob`` was renamed
    ``ci_prob`` and 1.x's default CI is the ETI, not the HDI) with the
    columns renamed to the 0.23 layout
    (``mean, sd, hdi_<p>%/hdi_<100-p>%/mcse_mean, mcse_sd, ess_bulk,
    ess_tail, r_hat``) and the circular columns recomputed for the named
    rows (rounded the way 0.23's auto rounding did).
    """
    import arviz as az
    from arviz_base.labels import IdxLabeller

    base = az.summary(idata, var_names=var_names, ci_prob=hdi_prob,
                      ci_kind="hdi", kind="all", skipna=skipna, group=group)
    if not circ_var_names:
        return base

    pct_lo = f"{100 * (1 - hdi_prob) / 2:g}"
    pct_hi = f"{100 * (1 + hdi_prob) / 2:g}"
    lb_col = next(c for c in base.columns if c.endswith("_lb"))
    ub_col = next(c for c in base.columns if c.endswith("_ub"))
    out = base.rename(columns={lb_col: f"hdi_{pct_lo}%", ub_col: f"hdi_{pct_hi}%"})
    out = out[["mean", "sd", f"hdi_{pct_lo}%", f"hdi_{pct_hi}%",
               "mcse_mean", "mcse_sd", "ess_bulk", "ess_tail", "r_hat"]]

    labeller = IdxLabeller()
    nan_policy = "omit" if skipna else "propagate"
    ds = getattr(idata, group)
    for var in circ_var_names:
        for name, sel, isel, values in az.xarray_var_iter(ds[var],
                                                          skip_dims=("chain", "draw")):
            label = labeller.make_label_flat(name, sel, isel)
            if label not in out.index:
                continue
            a = np.asarray(values, float).ravel()
            h_lo, h_hi = circ_hdi(a, hdi_prob)
            out.loc[label, "mean"] = round(float(scs.circmean(a, high=np.pi, low=-np.pi,
                                                              nan_policy=nan_policy)), 3)
            out.loc[label, "sd"] = round(float(scs.circstd(a, high=np.pi, low=-np.pi,
                                                           nan_policy=nan_policy)), 3)
            out.loc[label, f"hdi_{pct_lo}%"] = round(float(h_lo), 3)
            out.loc[label, f"hdi_{pct_hi}%"] = round(float(h_hi), 3)
            out.loc[label, "mcse_mean"] = round(_circular_mcse(a), 3)
    return out


# ---------------------------------------------------------------------------
# Posterior bookkeeping + bootstrap
# ---------------------------------------------------------------------------

def stack_posterior_1d(post, name):
    """Flatten an arviz posterior variable to a 1-D ``(chain*draw)`` array.

    ``post`` is the ``idata.posterior`` xarray; missing variables raise
    ``KeyError`` so a typo cannot silently zero out a panel.
    """
    if name not in post.data_vars:
        raise KeyError(f"Missing posterior variable '{name}'.")
    return post[name].stack(sample=("chain", "draw")).values.astype(float)


def bootstrap_track_condition(y, B=1200, seed=0):
    """Nonparametric bootstrap of ``(mu, R)`` for one track's heading vector.

    Returns two ``(B,)`` arrays of bootstrapped circular means and mean
    resultant lengths; empty ``y`` yields empty arrays.
    """
    rng = np.random.default_rng(seed)
    y = np.asarray(y, float)
    if y.size == 0:
        return np.array([]), np.array([])
    z = np.exp(1j * y)
    n = z.size

    mu = np.empty(B, float)
    R = np.empty(B, float)

    for b in range(B):
        idx = rng.integers(0, n, size=n)
        zbar = z[idx].mean()
        mu[b] = circ_mod(np.angle(zbar))
        R[b] = np.abs(zbar)

    return mu, R


def hierarchical_bootstrap_group(df, B=1500, seed=123):
    """Two-level (individual -> step) bootstrap of group-level summaries.

    Expects a long DataFrame with columns ``Track`` (id), ``PostStimulus``
    (0/1 condition) and ``HeadingRadians``. Each replicate resamples
    individuals with replacement AND resamples the heading steps inside each
    chosen individual, then reports per condition:

    * group mean direction (angle of the mean of individual mean phasors), and
    * group alignment ``R`` (mean resultant length of the individual mean
      directions, i.e. alignment between individuals).

    Individuals whose resampled step vector averages to zero are skipped.
    Returns a dict with ``mu_pre, mu_post, R_pre, R_post`` (each ``B``) and
    the ``tracks`` ordering used.
    """
    rng = np.random.default_rng(seed)
    df = df.copy()
    df["Track"] = df["Track"].astype(str)

    tracks = sorted_groups(df["Track"])
    nT = len(tracks)

    unit = {}
    post_unit = {}
    for t in tracks:
        y_pre = df.loc[(df["Track"] == t) & (df["PostStimulus"] == 0), "HeadingRadians"].dropna().to_numpy(float)
        y_post = df.loc[(df["Track"] == t) & (df["PostStimulus"] == 1), "HeadingRadians"].dropna().to_numpy(float)
        unit[t] = np.exp(1j * y_pre) if y_pre.size else None
        post_unit[t] = np.exp(1j * y_post) if y_post.size else None

    mu_pre_s = np.empty(B, float)
    mu_post_s = np.empty(B, float)
    R_pre_s = np.empty(B, float)
    R_post_s = np.empty(B, float)

    for b in range(B):
        samp_tracks = rng.choice(tracks, size=nT, replace=True)

        u_pre = []
        u_post = []

        for t in samp_tracks:
            z = unit[t]
            if z is not None:
                idx = rng.integers(0, z.size, size=z.size)
                zbar = z[idx].mean()
                if np.abs(zbar) > 0:
                    u_pre.append(zbar / np.abs(zbar))

            z = post_unit[t]
            if z is not None:
                idx = rng.integers(0, z.size, size=z.size)
                zbar = z[idx].mean()
                if np.abs(zbar) > 0:
                    u_post.append(zbar / np.abs(zbar))

        if len(u_pre) > 0:
            zG = np.mean(u_pre)
            mu_pre_s[b] = circ_mod(np.angle(zG))
            R_pre_s[b] = np.abs(zG)
        else:
            mu_pre_s[b] = np.nan
            R_pre_s[b] = np.nan

        if len(u_post) > 0:
            zG = np.mean(u_post)
            mu_post_s[b] = circ_mod(np.angle(zG))
            R_post_s[b] = np.abs(zG)
        else:
            mu_post_s[b] = np.nan
            R_post_s[b] = np.nan

    return {
        "mu_pre": mu_pre_s,
        "mu_post": mu_post_s,
        "R_pre": R_pre_s,
        "R_post": R_post_s,
        "tracks": tracks,
    }


def sorted_groups(x):
    """Unique group ids, numerically sorted where possible, else lexicographic."""
    import pandas as pd

    vals = pd.Series(x).dropna().astype(str).unique().tolist()
    try:
        return sorted(vals, key=lambda v: float(v))
    except (TypeError, ValueError):
        return sorted(vals)
