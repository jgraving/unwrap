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
from scipy.stats import circmean

from unwrap.inference import graph_free_rvs
from unwrap.link import softplus as _softplus
from unwrap.circular import (
    A_of_kappa as _A_of_kappa,
    bootstrap_track_condition as _bootstrap_track_condition,
    circ_mean as _circ_mean,
    circ_mod as _circ_mod,
    circular_mean as _circ_mean_centered,
    hierarchical_bootstrap_group as _hierarchical_bootstrap_group,
    linear_hdi as _linear_hdi,
    resultant_length_angles as _resultant_length,
    sorted_groups as _sorted_tracks,
    stack_posterior_1d as _stack_1d,
    vonmises_logpdf,
    vonmises_pdf,
    wrap_0_2pi as _wrap_0_2pi,
)


# ===========================================================================
# Circular helpers (shared by the entry points below)
# ===========================================================================


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


def _circ_hdi_centered(samples, prob=0.95):
    a = np.asarray(samples, float)
    m = _circ_mean_centered(a)
    d = _circ_mod(a - m)
    lo_d, hi_d = _linear_hdi(d, prob)
    lo = float(_circ_mod(m + lo_d))
    hi = float(_circ_mod(m + hi_d))
    return m, lo, hi


def posterior_directions(idata, df_blocks):
    """Per-draw per-track mean directions and concentrations from the fit.

    The fixed and random effect nodes are located in the posterior by
    pattern matching, so the extraction is robust to ``bambi``'s internal
    naming. Auxiliary nodes are excluded by name suffix: the random-effect
    ``*_kappa`` / ``*_sigma`` hyperparameter scales sit next to the
    random-effect value nodes in the posterior and share their factor
    dimension, and the ``*_offset`` factors of the non-centred kappa
    random effects (the stored effect value is ``sigma * offset``) — a
    pattern match on ``|Track`` alone would pick any of these up
    incorrectly.

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
                and (not v.endswith(("_kappa", "_sigma", "_offset"))) and any(cn not in ("chain", "draw") for cn in post[v].coords):
            mu_g_name = v
            break

    # random slope for PostStimulus by Track for mu
    mu_g_post_name = None
    for v in post.data_vars:
        if ("PostStimulus|Track" in v) and (not v.startswith("kappa_")) \
                and (not v.endswith(("_kappa", "_sigma", "_offset"))) and any(cn not in ("chain", "draw") for cn in post[v].coords):
            mu_g_post_name = v
            break

    # fixed effects (kappa)
    k_intercept_name = "kappa_Intercept" if "kappa_Intercept" in post.data_vars else None
    k_post_name = "kappa_PostStimulus" if "kappa_PostStimulus" in post.data_vars else None

    # group effect for kappa (random intercept)
    k_g_name = None
    for v in post.data_vars:
        if v.startswith("kappa_") and ("|Track" in v) and ("PostStimulus|Track" not in v) \
                and (not v.endswith(("_kappa", "_sigma", "_offset"))) and any(cn not in ("chain", "draw") for cn in post[v].coords):
            k_g_name = v
            break

    # random slope for PostStimulus by Track for kappa
    k_g_post_name = None
    for v in post.data_vars:
        if v.startswith("kappa_") and ("PostStimulus|Track" in v) \
                and (not v.endswith(("_kappa", "_sigma", "_offset"))) and any(cn not in ("chain", "draw") for cn in post[v].coords):
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
        vonmises_pdf(theta[None, :], mu_group_pre[:, None], kappa_bw_group),
        axis=0
    )
    kde_mu_group_post = np.mean(
        vonmises_pdf(theta[None, :], mu_group_post[:, None], kappa_bw_group),
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

    R_ind_pre_samples  = _A_of_kappa(kappa_pre)   # (S, nT)
    R_ind_post_samples = _A_of_kappa(kappa_post)  # (S, nT)

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
            vonmises_pdf(theta[None, :], mu_ind_pre[:, j][:, None], kappa_bw_ind),
            axis=0
        )
        kde_mu_ind_post[j] = np.mean(
            vonmises_pdf(theta[None, :], mu_ind_post[:, j][:, None], kappa_bw_ind),
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


def ppc_bootstrap_stats(idata, model, df_blocks):
    """Data-based (hierarchical bootstrap) and model-based (posterior)
    summaries for the PPC figure.

    Empirical side: a hierarchical bootstrap (resample tracks, resample
    windows within tracks) at the group level and per track-condition.
    Posterior side: per-observation response parameters (mean direction and
    concentration) read straight from the PyMC graph (via
    :func:`graph_free_rvs`, which reconstructs the non-centered kappa
    ``*_offset`` free RVs that bambi drops from the stored posterior),
    cross-checked against the per-observation log-likelihood, then
    aggregated by Track x Condition and compared against the bootstraps.
    All 95% HDI summaries and the per-track tables returned come from this
    dict.
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
        graph_free_rvs(idata, model),
        var_names=["mu", "kappa"],
        model=model.backend.model,
        progressbar=False,
    )
    mu_da, kappa_da = params["mu"], params["kappa"]
    assert np.isfinite(mu_da.values).all()
    assert np.isfinite(kappa_da.values).all()
    assert (kappa_da.values > 0).all()

    # Cross-check the graph parameters against the per-observation
    # log-likelihood stored in the fit result: they must agree to
    # floating-point precision.
    _obs_ll = idata.log_likelihood["HeadingRadians"]
    _expected = vonmises_logpdf(
        df_blocks["HeadingRadians"].to_numpy(),
        mu_da.transpose(*_obs_ll.dims).values,
        kappa_da.transpose(*_obs_ll.dims).values,
    )
    np.testing.assert_allclose(_expected, _obs_ll.values, rtol=1e-10, atol=1e-10)
    print("response params from PyMC graph OK; log-likelihood cross-check passed")

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
