"""Post-hoc computations for notebooks/zebra_escape.ipynb.

Three entry points, each producing everything one figure cell displays (the
matplotlib "how" lives in zebra_plots.py, which calls into this module so
the notebook cells stay one call):

- :func:`posterior_directions`: per-draw per-track mean directions and
  concentrations (design grid, linear predictors, then the links);
- :func:`summary_stats`: the group/individual quantities behind the 2x2
  posterior summary figure;
- :func:`ppc_bootstrap_stats`: data-based (hierarchical bootstrap) and
  model-based (posterior) summaries behind the PPC figure, including the
  per-observation response parameters read straight from the PyMC graph.
"""

import numpy as np
import pandas as pd
import pymc as pm
from scipy.special import i0, i1
from scipy.stats import circmean, vonmises


# ===========================================================================
# Circular helpers (shared by the entry points below)
# ===========================================================================
def _softplus(x):
    return np.log1p(np.exp(-np.abs(x))) + np.maximum(x, 0.0)


def _circ_mod(x):
    # wrap to (-pi, pi]
    return (x + np.pi) % (2 * np.pi) - np.pi


def _wrap_0_2pi(x):
    return x % (2 * np.pi)


def _circ_mean(a, axis=None):
    a = np.asarray(a)
    s = np.sin(a).mean(axis=axis)
    c = np.cos(a).mean(axis=axis)
    return (np.arctan2(s, c)) % (2 * np.pi)


def _circ_mean_centered(theta):
    theta = np.asarray(theta, float)
    z = np.exp(1j * theta).mean()
    return float(_circ_mod(np.angle(z)))


def _circ_hdi(samples, prob=0.95):
    """Smallest arc on [0, 2pi) containing `prob` mass."""
    a = np.sort(samples % (2 * np.pi))
    n = len(a)
    k = int(np.ceil(prob * n))
    a_ext = np.concatenate([a, a + 2 * np.pi])
    spans = a_ext[k - 1:n + k - 1] - a
    j = np.argmin(spans)
    low = a[j]
    high = (a_ext[j + k - 1]) % (2 * np.pi)
    return low, high


def _as_samples(da):
    """Flatten xarray (chain, draw, *levels) into (S, *levels)."""
    vals = da.values
    if vals.ndim == 2:  # (chain, draw)
        c, d = vals.shape
        return vals.reshape(c * d)
    elif vals.ndim == 1:
        return vals
    elif vals.ndim == 3:  # (chain, draw, level)
        c, d, L = vals.shape
        return vals.reshape(c * d, L)
    else:
        raise ValueError(f"Unexpected shape {vals.shape}")


def _circ_mean_series(a):
    a = np.asarray(a)
    if a.size == 0:
        return np.nan
    return _circ_mean(_wrap_0_2pi(a))


def _arc_width(a, b):
    return (b - a) % (2 * np.pi)


def _linear_hdi(samples, prob=0.95):
    s = np.sort(np.asarray(samples, float))
    n = s.size
    if n == 0:
        return np.nan, np.nan
    k = int(np.ceil(prob * n))
    spans = s[k - 1:] - s[:n - k + 1]
    j = int(np.argmin(spans))
    return float(s[j]), float(s[j + k - 1])


def _resultant_length(a, axis=-1):
    a = np.asarray(a)
    return np.sqrt(np.mean(np.cos(a), axis=axis) ** 2 +
                   np.mean(np.sin(a), axis=axis) ** 2)


def _circ_hdi_centered(samples, prob=0.95):
    a = np.asarray(samples, float)
    m = _circ_mean_centered(a)
    d = _circ_mod(a - m)
    lo_d, hi_d = _linear_hdi(d, prob)
    lo = float(_circ_mod(m + lo_d))
    hi = float(_circ_mod(m + hi_d))
    return m, lo, hi


def _circ_xerr_from_interval(m, lo, hi):
    w_lo = abs(float(_circ_mod(m - lo)))
    w_hi = abs(float(_circ_mod(hi - m)))
    return np.degrees(w_lo), np.degrees(w_hi)


def _A_of_kappa(kappa):
    kappa = np.asarray(kappa, float)
    return i1(kappa) / i0(kappa)


def _stack_1d(post, name):
    if name not in post.data_vars:
        raise KeyError(f"Missing posterior variable '{name}'.")
    return post[name].stack(sample=("chain", "draw")).values.astype(float)


# ===========================================================================
# Entry point 1: per-draw per-track directions and concentrations
# ===========================================================================
def posterior_directions(idata, df_blocks):
    """Per-draw per-track mean directions and concentrations from the fit.

    The fixed and random effect nodes are located in the posterior by
    pattern matching, so the extraction is robust to ``bambi``'s internal
    naming. Hyperparameter nodes (the random-effect ``*_kappa`` and
    ``*_sigma`` scales) are excluded by name suffix: in the current
    ``bambi`` they sit next to the random-effect value nodes in the
    posterior and share their factor dimension, which a pattern match on
    ``|Track`` alone would pick up incorrectly.

    The Track x Condition design grid is evaluated per posterior draw and
    the circular-modulo / softplus links are applied, recovering the
    per-track direction and concentration in each condition.

    Returns
    -------
    dict
        ``mu``, ``kappa`` : (S, 2*nT) arrays - pre rows first, post after;
        ``tracks``, ``nT``, ``idx_pre``, ``idx_post`` : grid bookkeeping;
        ``grid`` : the (2*nT, 2) design table.
    """
    post = idata.posterior  # assumes `idata` already exists

    # fixed effects (mu)
    mu_intercept_name = "Intercept"
    mu_post_name = "PostStimulus" if "PostStimulus" in post.data_vars else None

    # group effect for mu (random intercept per Track)
    mu_g_name = None
    for v in post.data_vars:
        if ("|Track" in v) and (not v.startswith("kappa_")) and ("PostStimulus|Track" not in v) \
                and (not v.endswith(("_kappa", "_sigma"))) and any(cn not in ("chain", "draw") for cn in post[v].coords):
            mu_g_name = v
            break

    # random slope for PostStimulus by Track for mu
    mu_g_post_name = None
    for v in post.data_vars:
        if ("PostStimulus|Track" in v) and (not v.startswith("kappa_")) \
                and (not v.endswith(("_kappa", "_sigma"))) and any(cn not in ("chain", "draw") for cn in post[v].coords):
            mu_g_post_name = v
            break

    # fixed effects (kappa)
    k_intercept_name = "kappa_Intercept" if "kappa_Intercept" in post.data_vars else None
    k_post_name = "kappa_PostStimulus" if "kappa_PostStimulus" in post.data_vars else None

    # group effect for kappa (random intercept)
    k_g_name = None
    for v in post.data_vars:
        if v.startswith("kappa_") and ("|Track" in v) and ("PostStimulus|Track" not in v) \
                and (not v.endswith(("_kappa", "_sigma"))) and any(cn not in ("chain", "draw") for cn in post[v].coords):
            k_g_name = v
            break

    # random slope for PostStimulus by Track for kappa
    k_g_post_name = None
    for v in post.data_vars:
        if v.startswith("kappa_") and ("PostStimulus|Track" in v) \
                and (not v.endswith(("_kappa", "_sigma"))) and any(cn not in ("chain", "draw") for cn in post[v].coords):
            k_g_post_name = v
            break

    # sanity checks
    if mu_intercept_name not in post:
        raise RuntimeError("Could not find mu Intercept in posterior.")
    if (k_intercept_name is None) or (k_intercept_name not in post):
        raise RuntimeError("Could not find kappa_Intercept in posterior.")

    # Tracks and grid: each Track x {pre (0), post (1)}
    tracks = np.sort(df_blocks["Track"].unique())
    nT = len(tracks)

    grid = pd.DataFrame({
        "PostStimulus": np.repeat([0, 1], nT),
        "Track":        np.tile(tracks, 2),
    })

    idx_pre  = np.arange(0, nT)
    idx_post = np.arange(nT, 2 * nT)

    # Build linear predictors eta_mu, eta_kappa per posterior sample
    S = post.sizes["chain"] * post.sizes["draw"]
    eta_mu   = np.zeros((S, 2 * nT))
    eta_kapp = np.zeros((S, 2 * nT))

    # intercepts
    eta_mu   += _as_samples(post[mu_intercept_name])[:, None]
    eta_kapp += _as_samples(post[k_intercept_name])[:, None]

    # fixed effect PostStimulus
    if mu_post_name is not None:
        b = _as_samples(post[mu_post_name])  # (S,)
        eta_mu[:, idx_post] += b[:, None]
    if k_post_name is not None:
        bk = _as_samples(post[k_post_name])
        eta_kapp[:, idx_post] += bk[:, None]

    # random intercepts per Track for mu
    if mu_g_name is not None:
        mu_g = _as_samples(post[mu_g_name])  # (S, L)
        coord_name = [cn for cn in post[mu_g_name].coords if cn not in ("chain", "draw")][0]
        levels = np.array(post[mu_g_name][coord_name].values).astype(str)
        level_idx = {lev: i for i, lev in enumerate(levels)}
        miss = [t for t in tracks if t not in level_idx]
        if miss:
            raise RuntimeError(f"Tracks missing in posterior RE for mu: {miss}")
        for j, t in enumerate(tracks):
            eta_mu[:, j]       += mu_g[:, level_idx[t]]   # pre row
            eta_mu[:, nT + j]  += mu_g[:, level_idx[t]]   # post row

    # random slope PostStimulus|Track for mu (post rows only)
    if mu_g_post_name is not None:
        mu_g_post = _as_samples(post[mu_g_post_name])  # (S, L)
        coord_name_mu_post = [cn for cn in post[mu_g_post_name].coords
                              if cn not in ("chain", "draw")][0]
        levels_mu_post = np.array(post[mu_g_post_name][coord_name_mu_post].values).astype(str)
        level_idx_mu_post = {lev: i for i, lev in enumerate(levels_mu_post)}
        miss_mu_post = [t for t in tracks if t not in level_idx_mu_post]
        if miss_mu_post:
            raise RuntimeError(f"Tracks missing in mu PostStimulus|Track RE: {miss_mu_post}")
        for j, t in enumerate(tracks):
            eta_mu[:, nT + j] += mu_g_post[:, level_idx_mu_post[t]]

    # random intercepts per Track for kappa
    if k_g_name is not None:
        kg = _as_samples(post[k_g_name])  # (S, L)
        coord_name_k = [cn for cn in post[k_g_name].coords if cn not in ("chain", "draw")][0]
        levels_k = np.array(post[k_g_name][coord_name_k].values).astype(str)
        level_idx_k = {lev: i for i, lev in enumerate(levels_k)}
        missk = [t for t in tracks if t not in level_idx_k]
        if missk:
            raise RuntimeError(f"Tracks missing in posterior RE for kappa: {missk}")
        for j, t in enumerate(tracks):
            eta_kapp[:, j]      += kg[:, level_idx_k[t]]
            eta_kapp[:, nT + j] += kg[:, level_idx_k[t]]

    # random slope kappa_PostStimulus|Track (post rows only)
    if k_g_post_name is not None:
        kg_post = _as_samples(post[k_g_post_name])  # (S, L)
        coord_name_kg_post = [cn for cn in post[k_g_post_name].coords
                              if cn not in ("chain", "draw")][0]
        levels_kg_post = np.array(post[k_g_post_name][coord_name_kg_post].values).astype(str)
        level_idx_kg_post = {lev: i for i, lev in enumerate(levels_kg_post)}
        misskg_post = [t for t in tracks if t not in level_idx_kg_post]
        if misskg_post:
            raise RuntimeError(f"Tracks missing in kappa PostStimulus|Track RE: {misskg_post}")
        for j, t in enumerate(tracks):
            eta_kapp[:, nT + j] += kg_post[:, level_idx_kg_post[t]]

    # Apply links
    mu    = _circ_mod(eta_mu)     # (-pi, pi]
    kappa = _softplus(eta_kapp)   # > 0

    return {
        "mu": mu,
        "kappa": kappa,
        "tracks": tracks,
        "nT": nT,
        "idx_pre": idx_pre,
        "idx_post": idx_post,
        "grid": grid,
    }


# ===========================================================================
# Entry point 2: quantities for the 2x2 posterior summary figure
# ===========================================================================
def summary_stats(idata, df_blocks):
    """Group-level and individual-level posterior quantities.

    From the per-draw per-track directions and concentrations
    (:func:`posterior_directions`), computes the group-level mean-direction
    draws and their KDE radii, the group alignment (resultant length of the
    per-track mean directions) with 95% HDIs, the individual mean-direction
    KDE radii and empirical per-track means, and the individual
    concentration (kappa -> mean resultant length via A(kappa)) means and
    95% HDIs. Everything the two-by-two figure draws comes from this dict.
    """
    post = posterior_directions(idata, df_blocks)
    mu = post["mu"]
    kappa = post["kappa"]
    tracks = post["tracks"]
    nT = post["nT"]
    idx_pre = post["idx_pre"]
    idx_post = post["idx_post"]
    S = mu.shape[0]

    # GROUP-LEVEL mean direction per draw (the "population level")
    # mu is (-pi, pi]; wrap before averaging
    pre_mu_draw  = _circ_mean(_wrap_0_2pi(mu[:, idx_pre]),  axis=1)   # (S,)
    post_mu_draw = _circ_mean(_wrap_0_2pi(mu[:, idx_post]), axis=1)   # (S,)

    mu_group_pre  = _circ_mod(pre_mu_draw)   # (-pi, pi]
    mu_group_post = _circ_mod(post_mu_draw)

    # Group-level KDEs (von Mises kernel density over posterior draws)
    theta = np.linspace(-np.pi, np.pi, 200)
    kappa_bw_group = 500

    kde_mu_group_pre = np.mean(
        vonmises.pdf(theta[None, :], kappa=kappa_bw_group, loc=mu_group_pre[:, None]),
        axis=0
    )
    kde_mu_group_post = np.mean(
        vonmises.pdf(theta[None, :], kappa=kappa_bw_group, loc=mu_group_post[:, None]),
        axis=0
    )

    # Group alignment = R of individual mean directions (posterior, per draw)
    mu_ind_pre  = _circ_mod(mu[:, idx_pre])    # (S, nT)
    mu_ind_post = _circ_mod(mu[:, idx_post])   # (S, nT)

    R_group_pre_draw  = _resultant_length(mu_ind_pre,  axis=1)  # (S,)
    R_group_post_draw = _resultant_length(mu_ind_post, axis=1)  # (S,)

    R_group_pre_mean,  R_group_pre_lo,  R_group_pre_hi  = (
        float(np.mean(R_group_pre_draw)), *_linear_hdi(R_group_pre_draw, prob=0.95))
    R_group_post_mean, R_group_post_lo, R_group_post_hi = (
        float(np.mean(R_group_post_draw)), *_linear_hdi(R_group_post_draw, prob=0.95))

    # INDIVIDUAL-LEVEL concentration: convert kappa -> mean resultant length
    kappa_pre  = kappa[:, idx_pre]   # (S, nT)
    kappa_post = kappa[:, idx_post]  # (S, nT)

    R_ind_pre_samples  = i1(kappa_pre)  / i0(kappa_pre)   # (S, nT)
    R_ind_post_samples = i1(kappa_post) / i0(kappa_post)  # (S, nT)

    R_ind_pre_mean  = R_ind_pre_samples.mean(axis=0)
    R_ind_post_mean = R_ind_post_samples.mean(axis=0)

    R_ind_pre_hdi  = np.array([_linear_hdi(R_ind_pre_samples[:, j],  prob=0.95) for j in range(nT)])
    R_ind_post_hdi = np.array([_linear_hdi(R_ind_post_samples[:, j], prob=0.95) for j in range(nT)])

    # Empirical individual means (dots on the individual polar panel)
    emp_mu_pre = []
    emp_mu_post = []
    for t in tracks:
        d_t = df_blocks[df_blocks["Track"] == t]
        h_pre  = d_t.loc[d_t["PostStimulus"] == 0, "HeadingRadians"].dropna().values
        h_post = d_t.loc[d_t["PostStimulus"] == 1, "HeadingRadians"].dropna().values
        m_pre  = np.nan if h_pre.size  == 0 else circmean(h_pre,  high=np.pi, low=-np.pi)
        m_post = np.nan if h_post.size == 0 else circmean(h_post, high=np.pi, low=-np.pi)
        emp_mu_pre.append(m_pre)
        emp_mu_post.append(m_post)
    emp_mu_pre = np.array(emp_mu_pre)
    emp_mu_post = np.array(emp_mu_post)

    # INDIVIDUAL-LEVEL KDEs for mu (thin curves per individual)
    kappa_bw_ind = 300
    kde_mu_ind_pre  = np.zeros((nT, theta.size))
    kde_mu_ind_post = np.zeros((nT, theta.size))
    for j in range(nT):
        kde_mu_ind_pre[j] = np.mean(
            vonmises.pdf(theta[None, :], kappa=kappa_bw_ind, loc=mu_ind_pre[:, j][:, None]),
            axis=0
        )
        kde_mu_ind_post[j] = np.mean(
            vonmises.pdf(theta[None, :], kappa=kappa_bw_ind, loc=mu_ind_post[:, j][:, None]),
            axis=0
        )

    # Shared scaling for polar densities (1 = baseline, 2 = max)
    max_kde = max(
        kde_mu_group_pre.max(), kde_mu_group_post.max(),
        kde_mu_ind_pre.max(), kde_mu_ind_post.max()
    )

    def scale_kde_to_unit(kde):
        return 1.0 + (kde / (max_kde + 1e-12))

    r_group_pre  = scale_kde_to_unit(kde_mu_group_pre)
    r_group_post = scale_kde_to_unit(kde_mu_group_post)
    r_ind_pre    = scale_kde_to_unit(kde_mu_ind_pre)
    r_ind_post   = scale_kde_to_unit(kde_mu_ind_post)

    # Group-level mean directions (single scalar each) for the polar dots.
    mu_group_pre_mean = _circ_mod(_circ_mean(_wrap_0_2pi(mu_group_pre)))
    mu_group_post_mean = _circ_mod(_circ_mean(_wrap_0_2pi(mu_group_post)))

    return {
        "nT": nT,
        "theta": theta,
        "r_group_pre": r_group_pre,
        "r_group_post": r_group_post,
        "r_ind_pre": r_ind_pre,
        "r_ind_post": r_ind_post,
        "mu_group_pre_mean": mu_group_pre_mean,
        "mu_group_post_mean": mu_group_post_mean,
        "R_group_pre_mean": R_group_pre_mean,
        "R_group_pre_lo": R_group_pre_lo,
        "R_group_pre_hi": R_group_pre_hi,
        "R_group_post_mean": R_group_post_mean,
        "R_group_post_lo": R_group_post_lo,
        "R_group_post_hi": R_group_post_hi,
        "emp_mu_pre": emp_mu_pre,
        "emp_mu_post": emp_mu_post,
        "R_ind_pre_mean": R_ind_pre_mean,
        "R_ind_post_mean": R_ind_post_mean,
        "R_ind_pre_hdi": R_ind_pre_hdi,
        "R_ind_post_hdi": R_ind_post_hdi,
    }


# ===========================================================================
# Entry point 3: bootstrap (data) vs posterior (model) for the PPC figure
# ===========================================================================
def _sorted_tracks(x):
    vals = pd.Series(x).dropna().astype(str).unique().tolist()
    try:
        return sorted(vals, key=lambda v: float(v))
    except Exception:
        return sorted(vals)


def _hierarchical_bootstrap_group(df, B=1500, seed=123):
    """Hierarchical bootstrap (Individuals -> steps) for group-level stats.

    Returns:
      mu_pre, mu_post : group mean directions (angles from individual mean directions)
      R_pre, R_post   : alignment of individual mean directions within condition (Pre/Post)
    """
    rng = np.random.default_rng(seed)
    df = df.copy()
    df["Track"] = df["Track"].astype(str)

    tracks = _sorted_tracks(df["Track"])
    nT = len(tracks)

    pre = {}
    post = {}
    for t in tracks:
        y_pre  = df.loc[(df["Track"] == t) & (df["PostStimulus"] == 0), "HeadingRadians"].dropna().to_numpy(float)
        y_post = df.loc[(df["Track"] == t) & (df["PostStimulus"] == 1), "HeadingRadians"].dropna().to_numpy(float)
        pre[t]  = np.exp(1j * y_pre)  if y_pre.size  else None
        post[t] = np.exp(1j * y_post) if y_post.size else None

    mu_pre_s  = np.empty(B, float)
    mu_post_s = np.empty(B, float)
    R_pre_s   = np.empty(B, float)
    R_post_s  = np.empty(B, float)

    for b in range(B):
        samp_tracks = rng.choice(tracks, size=nT, replace=True)

        u_pre  = []
        u_post = []

        for t in samp_tracks:
            z = pre[t]
            if z is not None:
                idx = rng.integers(0, z.size, size=z.size)
                zbar = z[idx].mean()
                if np.abs(zbar) > 0:
                    u_pre.append(zbar / np.abs(zbar))

            z = post[t]
            if z is not None:
                idx = rng.integers(0, z.size, size=z.size)
                zbar = z[idx].mean()
                if np.abs(zbar) > 0:
                    u_post.append(zbar / np.abs(zbar))

        if len(u_pre) > 0:
            zG = np.mean(u_pre)
            mu_pre_s[b] = _circ_mod(np.angle(zG))
            R_pre_s[b]  = np.abs(zG)
        else:
            mu_pre_s[b] = np.nan
            R_pre_s[b]  = np.nan

        if len(u_post) > 0:
            zG = np.mean(u_post)
            mu_post_s[b] = _circ_mod(np.angle(zG))
            R_post_s[b]  = np.abs(zG)
        else:
            mu_post_s[b] = np.nan
            R_post_s[b]  = np.nan

    return {
        "mu_pre": mu_pre_s,
        "mu_post": mu_post_s,
        "R_pre": R_pre_s,
        "R_post": R_post_s,
        "tracks": tracks,
    }


def _bootstrap_track_condition(y, B=1200, seed=0):
    rng = np.random.default_rng(seed)
    y = np.asarray(y, float)
    if y.size == 0:
        return np.array([]), np.array([])
    z = np.exp(1j * y)
    n = z.size

    mu = np.empty(B, float)
    R  = np.empty(B, float)

    for b in range(B):
        idx = rng.integers(0, n, size=n)
        zbar = z[idx].mean()
        mu[b] = _circ_mod(np.angle(zbar))
        R[b]  = np.abs(zbar)

    return mu, R


def ppc_bootstrap_stats(idata, model, df_blocks):
    """Data-based (hierarchical bootstrap) and model-based (posterior)
    summaries for the PPC figure.

    Empirical side: a hierarchical bootstrap (resample tracks, resample
    windows within tracks) at the group level and per track-condition.
    Posterior side: per-observation response parameters (mean direction and
    concentration) read straight from the PyMC graph, aggregated by
    Track x Condition, and compared against the bootstraps. All 95% HDI
    summaries and the per-track tables returned come from this dict.
    """
    # ----------------------------
    # Data
    # ----------------------------
    df = df_blocks.copy()
    df["Track"] = (
        df["Track"].astype(str)
        .str.replace(r"^observation088-", "", regex=True)
    )
    tracks = _sorted_tracks(df["Track"])
    nT = len(tracks)

    # ----------------------------
    # Empirical bootstraps
    # ----------------------------
    B_GROUP = 2000
    B_INDIV = 1500

    bootG = _hierarchical_bootstrap_group(df, B=B_GROUP, seed=123)

    boot_mu_pre  = {}
    boot_mu_post = {}
    boot_R_pre   = {}
    boot_R_post  = {}

    for j, t in enumerate(tracks):
        y_pre  = df.loc[(df["Track"] == t) & (df["PostStimulus"] == 0), "HeadingRadians"].dropna().to_numpy(float)
        y_post = df.loc[(df["Track"] == t) & (df["PostStimulus"] == 1), "HeadingRadians"].dropna().to_numpy(float)

        mu_b, R_b = _bootstrap_track_condition(y_pre,  B=B_INDIV, seed=1000 + j)
        boot_mu_pre[t] = mu_b
        boot_R_pre[t]  = R_b

        mu_b, R_b = _bootstrap_track_condition(y_post, B=B_INDIV, seed=2000 + j)
        boot_mu_post[t] = mu_b
        boot_R_post[t]  = R_b

    # ----------------------------
    # Posterior stats for comparison
    # ----------------------------
    post = idata.posterior

    # fixed effects mu for group mean direction panel
    mu_intercept = _stack_1d(post, "Intercept")
    mu_post_fx   = _stack_1d(post, "PostStimulus") if "PostStimulus" in post.data_vars else np.zeros_like(mu_intercept)
    mu_pop_pre_draw  = _circ_mod(mu_intercept)
    mu_pop_post_draw = _circ_mod(mu_intercept + mu_post_fx)

    # Track-level posterior mu and kappa via response_params -> aggregate by Track x Condition
    # N.B.: bambi bug, tracked upstream at bambinos/bambi#1084:
    # Model.predict(kind="response_params") applies the vonmises *mean* link,
    # angle(exp(1j x)), to every component, so its "kappa" is the wrapped
    # linear predictor instead of softplus(linear predictor). Temporary fix
    # until the official bambi fix ships: read the response parameters straight
    # from the PyMC graph with compute_deterministics, which uses each
    # component's actual link.
    #
    # Official (currently broken for the kappa component):
    # preds = model.predict(idata, kind="response_params", inplace=False)
    # mu_da    = preds.posterior["mu"]
    # kappa_da = preds.posterior["kappa"]

    params = pm.compute_deterministics(
        idata.posterior,
        var_names=["mu", "kappa"],
        model=model.backend.model,
        progressbar=False,
    )
    mu_da, kappa_da = params["mu"], params["kappa"]
    assert np.isfinite(mu_da.values).all()
    assert np.isfinite(kappa_da.values).all()
    assert (kappa_da.values > 0).all()
    print("response params from PyMC graph OK (kappa all positive)")

    obs_dim = [d for d in mu_da.dims if d not in ("chain", "draw")][0]

    mu_samps    = mu_da.stack(sample=("chain", "draw")).transpose("sample", obs_dim).values
    kappa_samps = kappa_da.stack(sample=("chain", "draw")).transpose("sample", obs_dim).values
    S_total, n_obs = mu_samps.shape

    track_codes = pd.Categorical(df["Track"], categories=tracks).codes
    cond_codes  = df["PostStimulus"].to_numpy(int)
    group = track_codes + nT * cond_codes
    G = 2 * nT

    cnt = np.bincount(group, minlength=G).astype(float)
    cnt[cnt == 0] = np.nan

    z_mu = np.exp(1j * mu_samps)
    sum_re = np.vstack([np.bincount(group, weights=np.real(z_mu[s]), minlength=G) for s in range(S_total)])
    sum_im = np.vstack([np.bincount(group, weights=np.imag(z_mu[s]), minlength=G) for s in range(S_total)])
    mean_re = sum_re / cnt[None, :]
    mean_im = sum_im / cnt[None, :]
    mu_group_draw = _circ_mod(np.arctan2(mean_im, mean_re))  # (S, G)

    R_obsparam = _A_of_kappa(kappa_samps)
    sum_R = np.vstack([np.bincount(group, weights=R_obsparam[s], minlength=G) for s in range(S_total)])
    R_group_draw = sum_R / cnt[None, :]  # (S, G)

    mu_pre_draw  = mu_group_draw[:, 0:nT]
    mu_post_draw = mu_group_draw[:, nT:2 * nT]
    R_pre_draw   = R_group_draw[:, 0:nT]
    R_post_draw  = R_group_draw[:, nT:2 * nT]

    # Group alignment (posterior): alignment of individual mean directions per condition
    R_group_pre_draw  = _resultant_length(mu_pre_draw,  axis=1)
    R_group_post_draw = _resultant_length(mu_post_draw, axis=1)

    # ----------------------------
    # Summaries: 95% HDIs
    # ----------------------------
    # Group mean direction: bootstrap vs posterior fixed effects
    muG_pre_m_b,  muG_pre_lo_b,  muG_pre_hi_b  = _circ_hdi_centered(bootG["mu_pre"][~np.isnan(bootG["mu_pre"])],   prob=0.95)
    muG_post_m_b, muG_post_lo_b, muG_post_hi_b = _circ_hdi_centered(bootG["mu_post"][~np.isnan(bootG["mu_post"])], prob=0.95)
    muG_pre_m_p,  muG_pre_lo_p,  muG_pre_hi_p  = _circ_hdi_centered(mu_pop_pre_draw,  prob=0.95)
    muG_post_m_p, muG_post_lo_p, muG_post_hi_p = _circ_hdi_centered(mu_pop_post_draw, prob=0.95)

    # Group alignment: bootstrap Pre/Post and posterior Pre/Post
    R_pre_b = bootG["R_pre"][~np.isnan(bootG["R_pre"])]
    R_post_b = bootG["R_post"][~np.isnan(bootG["R_post"])]

    R_pre_m_b = float(np.nanmean(R_pre_b)) if R_pre_b.size else np.nan
    R_pre_lo_b, R_pre_hi_b = _linear_hdi(R_pre_b, prob=0.95) if R_pre_b.size else (np.nan, np.nan)

    R_post_m_b = float(np.nanmean(R_post_b)) if R_post_b.size else np.nan
    R_post_lo_b, R_post_hi_b = _linear_hdi(R_post_b, prob=0.95) if R_post_b.size else (np.nan, np.nan)

    R_pre_m_p = float(np.mean(R_group_pre_draw))
    R_pre_lo_p, R_pre_hi_p = _linear_hdi(R_group_pre_draw, prob=0.95)

    R_post_m_p = float(np.mean(R_group_post_draw))
    R_post_lo_p, R_post_hi_p = _linear_hdi(R_group_post_draw, prob=0.95)

    # Individual-level tables (mu and R)
    track_mu_rows = []
    track_R_rows  = []
    for j, t in enumerate(tracks):
        # bootstrap mu
        mu_b = boot_mu_pre[t]
        if mu_b.size:
            m, lo, hi = _circ_hdi_centered(mu_b, prob=0.95)
            track_mu_rows.append(("Pre", t, m, lo, hi, "bootstrap"))
        mu_b = boot_mu_post[t]
        if mu_b.size:
            m, lo, hi = _circ_hdi_centered(mu_b, prob=0.95)
            track_mu_rows.append(("Post", t, m, lo, hi, "bootstrap"))

        # posterior mu
        m, lo, hi = _circ_hdi_centered(mu_pre_draw[:, j], prob=0.95)
        track_mu_rows.append(("Pre", t, m, lo, hi, "posterior"))
        m, lo, hi = _circ_hdi_centered(mu_post_draw[:, j], prob=0.95)
        track_mu_rows.append(("Post", t, m, lo, hi, "posterior"))

        # bootstrap R
        Rb = boot_R_pre[t]
        if Rb.size:
            track_R_rows.append(("Pre", t, float(np.mean(Rb)), *_linear_hdi(Rb, prob=0.95), "bootstrap"))
        Rb = boot_R_post[t]
        if Rb.size:
            track_R_rows.append(("Post", t, float(np.mean(Rb)), *_linear_hdi(Rb, prob=0.95), "bootstrap"))

        # posterior R
        Rp = R_pre_draw[:, j]
        track_R_rows.append(("Pre", t, float(np.mean(Rp)), *_linear_hdi(Rp, prob=0.95), "posterior"))
        Rp = R_post_draw[:, j]
        track_R_rows.append(("Post", t, float(np.mean(Rp)), *_linear_hdi(Rp, prob=0.95), "posterior"))

    mu_tbl = pd.DataFrame(track_mu_rows, columns=["Condition", "Track", "mean", "low", "high", "Source"])
    R_tbl  = pd.DataFrame(track_R_rows,  columns=["Condition", "Track", "mean", "low", "high", "Source"])

    return {
        "tracks": tracks,
        "nT": nT,
        "muG_pre_m_b": muG_pre_m_b,  "muG_pre_lo_b": muG_pre_lo_b,  "muG_pre_hi_b": muG_pre_hi_b,
        "muG_post_m_b": muG_post_m_b, "muG_post_lo_b": muG_post_lo_b, "muG_post_hi_b": muG_post_hi_b,
        "muG_pre_m_p": muG_pre_m_p,  "muG_pre_lo_p": muG_pre_lo_p,  "muG_pre_hi_p": muG_pre_hi_p,
        "muG_post_m_p": muG_post_m_p, "muG_post_lo_p": muG_post_lo_p, "muG_post_hi_p": muG_post_hi_p,
        "R_pre_m_b": R_pre_m_b,  "R_pre_lo_b": R_pre_lo_b,  "R_pre_hi_b": R_pre_hi_b,
        "R_post_m_b": R_post_m_b, "R_post_lo_b": R_post_lo_b, "R_post_hi_b": R_post_hi_b,
        "R_pre_m_p": R_pre_m_p,  "R_pre_lo_p": R_pre_lo_p,  "R_pre_hi_p": R_pre_hi_p,
        "R_post_m_p": R_post_m_p, "R_post_lo_p": R_post_lo_p, "R_post_hi_p": R_post_hi_p,
        "mu_tbl": mu_tbl,
        "R_tbl": R_tbl,
    }
