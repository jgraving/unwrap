"""Post-hoc computations for notebooks/beetles.ipynb.

The change-in-dispersion ($\\kappa_\\mu$) analysis and the supplemental
posterior-predictive check are heavy post-hoc passes over the fitted model's
posterior: per-individual mean-direction estimates with von Mises fits of
their concentration (section "Extract individual $\\mu$ estimates"), the same
fit repeated for every MCMC draw (section "Extract $\\kappa_\\mu$ estimates
for each draw"), and the bootstrap (data) vs posterior (model) PPC statistics
(bootstrap resampling, per-draw per-group directions and concentrations from
the posterior, and a cross-check of the graph response parameters against the
log-likelihood bambi stored during fitting). The notebook calls these
functions as one-liners; figure construction is in ``beetles_plots.py``.

The computations are transcribed verbatim from the notebook's former inline
cells, including their upstream quirks (e.g. the intercept's failed-draw
fallback in :func:`kappamu_per_draw` writing into the light slot).
"""

import re

import numpy as np
import pandas as pd
import pymc as pm
from scipy.special import i0, i1
from scipy.stats import circmean
from scipy.stats import vonmises, directional_stats, vonmises as _vonmises

from unwrap.bessel import inverse_bessel_ratio
from unwrap.link.numpy import circular_modulo


# ---------------------------------------------------------------------------
# Circular-statistics helpers.
# ---------------------------------------------------------------------------
def circ_mod(x):
    # wrap to (-pi, pi]
    return (x + np.pi) % (2 * np.pi) - np.pi


def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', str(s))]


def circ_mean(theta):
    theta = np.asarray(theta, float)
    z = np.exp(1j * theta).mean()
    return float(circ_mod(np.angle(z)))


def linear_hdi(x, prob=0.95):
    x = np.sort(np.asarray(x, float))
    n = x.size
    if n == 0: return np.nan, np.nan
    k = int(np.ceil(prob * n))
    spans = x[k-1:] - x[:n-k+1]
    j = int(np.argmin(spans))
    return float(x[j]), float(x[j + k - 1])


def circ_hdi_centered(samples, prob=0.95):
    a = np.asarray(samples, float)
    if a.size == 0: return np.nan, np.nan, np.nan
    m = circ_mean(a)
    d = circ_mod(a - m)
    lo_d, hi_d = linear_hdi(d, prob)
    lo = float(circ_mod(m + lo_d))
    hi = float(circ_mod(m + hi_d))
    return m, lo, hi


def A_of_kappa(kappa):
    kappa = np.asarray(kappa, float)
    return i1(kappa) / i0(kappa)


def resultant_length_angles(a, axis=-1):
    a = np.asarray(a, float)
    return np.sqrt(np.mean(np.cos(a), axis=axis)**2 + np.mean(np.sin(a), axis=axis)**2)


def _stack_1d(post, name):
    if name not in post.data_vars:
        raise KeyError(f"Missing posterior variable '{name}'.")
    return post[name].stack(sample=("chain", "draw")).values.astype(float)


def rho_circular(angles):
    """
    Calculates the mean resultant length using sine and cosine decomposition.

    Parameters:
    angles (array-like): Input angles in radians.

    Returns:
    float: Mean resultant length (R-bar) between 0 and 1.
    """
    angles = np.asarray(angles)

    # 1. Decompose angles into x and y components on the unit circle
    x_coords = np.cos(angles)
    y_coords = np.sin(angles)

    # 2. Calculate the mean of the components
    C = np.mean(x_coords)
    S = np.mean(y_coords)

    # 3. Calculate the magnitude of the mean vector (R-bar)
    # This uses the Pythagorean theorem: sqrt(C^2 + S^2)
    r_bar = np.sqrt(C**2 + S**2)

    return r_bar


# ---------------------------------------------------------------------------
# Section 1: individual mu estimates + kappa_mu (was cell "Extract individual
# mu estimates and calculate mean kappa_mu").
# ---------------------------------------------------------------------------
def extract_mu_effects(idata_treat, df):
    """Fixed and random effects of the mean direction, plus the individual
    IDs, from the treatment fit's posterior (verbatim extraction)."""
    # --- 1.1. Extract Fixed Effects (Common-level) for mu ---
    mu_intercept_fixed = idata_treat.posterior["Intercept"].values  # Note: No 'mu_' prefix if it's the default component
    mu_light_fixed = idata_treat.posterior["Light"].values

    # --- 1.2. Extract Random Effects (Group-specific) for mu ---
    mu_intercept_random = idata_treat.posterior["1|Individual"].values
    mu_light_random = idata_treat.posterior["Light|Individual"].values

    # --- 1.3. Get the list of individual IDs ---
    individual_ids = df['Individual'].unique()
    n_individuals = len(individual_ids)

    return {
        "mu_intercept_fixed": mu_intercept_fixed,
        "mu_light_fixed": mu_light_fixed,
        "mu_intercept_random": mu_intercept_random,
        "mu_light_random": mu_light_random,
        "individual_ids": individual_ids,
        "n_individuals": n_individuals,
    }


def kappamu_individual(idata_treat, df, *, verbose=True):
    """Post-hoc estimate of the second-level concentration $\\kappa_\\mu$:
    the concentration of the (per-individual) distribution of individual
    mean-direction estimates, off and on.

    Returns a dict with the per-individual estimates and the four fitted
    (mixture) parameters. Prints the two ``--- Kappa_mu ... ---`` summary
    lines (as the original cell did).
    """
    e = extract_mu_effects(idata_treat, df)
    mu_intercept_fixed = e["mu_intercept_fixed"]
    mu_light_fixed = e["mu_light_fixed"]
    mu_intercept_random = e["mu_intercept_random"]
    mu_light_random = e["mu_light_random"]

    # Sum the fixed and random effects for each draw and individual
    mu_intercept_draws = circular_modulo(mu_intercept_fixed[..., np.newaxis] + mu_intercept_random)

    # Calculate the mean posterior direction (often the circular mean or just the arithmetic mean of the posterior)
    mu_intercept_individual_estimates = circular_modulo(circmean(mu_intercept_draws, axis=(0, 1)))

    # for the light on condition
    mu_light_draws = circular_modulo(
        mu_intercept_fixed[..., np.newaxis] + mu_light_fixed[..., np.newaxis] +
        mu_intercept_random + mu_light_random
    )

    # Calculate the mean posterior direction (point estimate)
    mu_light_individual_estimates = circular_modulo(circmean(mu_light_draws, axis=(0, 1)))

    # Fit the von Mises distribution to the set of individual mu estimates
    # fscale=1 fixes the scale parameter to 1
    kappa_mu_intercept_fit, mu_mu_intercept_fit, scale_fit = vonmises.fit(
        mu_intercept_individual_estimates, fscale=1
    )
    # N.B. an alternative method would be:
    kappa_mu_intercept_fit = inverse_bessel_ratio(rho_circular(mu_intercept_individual_estimates))
    if verbose:
        print(f"\n--- Kappa_mu from the Distribution of Individual Mu (Intercept) ---")
        print(f"The Kappa_mu (concentration) of the individual mu estimates: {kappa_mu_intercept_fit:.4f}")
    # The estimated mu_mu_intercept_fit is the overall mean of the individual means.

    # Fit the von Mises distribution to the set of individual mu estimates
    kappa_mu_light_fit, mu_mu_light_fit, scale_fit = vonmises.fit(
        mu_light_individual_estimates, fscale=1
    )
    # N.B. an alternative method would be:
    # kappa_mu_light_fit = inverse_bessel_ratio(rho_circular(mu_light_individual_estimates))

    if verbose:
        print(f"\n--- Kappa_mu from the Distribution of Individual Mu (Light) ---")
        print(f"The Kappa_mu (concentration) of the individual mu estimates: {kappa_mu_light_fit:.4f}")

    return {
        "mu_intercept_individual_estimates": mu_intercept_individual_estimates,
        "mu_light_individual_estimates": mu_light_individual_estimates,
        "kappa_mu_intercept_fit": float(kappa_mu_intercept_fit),
        "mu_mu_intercept_fit": mu_mu_intercept_fit,
        "kappa_mu_light_fit": float(kappa_mu_light_fit),
        "mu_mu_light_fit": mu_mu_light_fit,
    }


# ---------------------------------------------------------------------------
# Section 2: kappa_mu for each draw (was cell "Extract kappa_mu estimates for
# each draw").
# ---------------------------------------------------------------------------
def kappamu_per_draw(idata_treat, df, *, verbose=True):
    """Per-draw second-level $\\kappa_\\mu$ posterior (a von Mises fit of the
    per-individual mean directions repeated for every MCMC draw).

    Returns ``{"id_data": <flat ArviZ dict>, "kappamu_intercept_posterior",
    "kappamu_light_posterior", "kappamu_delta_posterior"}``. Prints the
    progress and the two summary header lines (as the original cell did);
    the ``az.summary`` display is left to the caller.
    """
    e = extract_mu_effects(idata_treat, df)
    mu_intercept_fixed = e["mu_intercept_fixed"]
    mu_light_fixed = e["mu_light_fixed"]
    mu_intercept_random = e["mu_intercept_random"]
    mu_light_random = e["mu_light_random"]

    # --- 2.2. Determine Iteration Count and Individual Count ---
    n_chain, n_draws, n_individuals = mu_intercept_random.shape
    total_draws = n_chain * n_draws

    # Initialize an array to store the second-level kappa for each draw
    kappamu_intercept_posterior = np.empty(total_draws)
    kappamu_light_posterior = np.empty(total_draws)
    kappamu_delta_posterior = np.empty(total_draws)

    # --- 2.3. Reshape all arrays for easier iteration (Combine chains and draws) ---
    # Shape is now (total_draws, n_individuals)
    mu_i_random_flat = mu_intercept_random.reshape(total_draws, n_individuals)
    mu_L_random_flat = mu_light_random.reshape(total_draws, n_individuals)
    mu_I_fixed_flat = mu_intercept_fixed.reshape(total_draws, 1)
    mu_L_fixed_flat = mu_light_fixed.reshape(total_draws, 1)

    # --- 2.4. Loop through every draw and fit the von Mises distribution ---
    if verbose:
        print(f"Calculating second-level kappamu for {total_draws} draws...")

    for i in range(total_draws):

        # 1. Calculate the individual-level mu for the current draw (Light condition)
        # The sum: fixed_I + fixed_L + random_I + random_L
        mu_intercept_individual_draw = (
            mu_I_fixed_flat[i] +
            mu_i_random_flat[i]
        )
        mu_light_individual_draw = (
            mu_I_fixed_flat[i] + mu_L_fixed_flat[i] +
            mu_i_random_flat[i] + mu_L_random_flat[i]
        )

        # 2. Wrap the angles to the preferred [-pi, pi) range
        mu_intercept_wrapped = circular_modulo(mu_intercept_individual_draw)
        mu_light_wrapped = circular_modulo(mu_light_individual_draw)

        # 3. Fit the von Mises distribution to this set of wrapped individual means
        # fscale=1 fixes the scale
        try:
            kappamu_intercept_fit, _, _ = vonmises.fit(mu_intercept_wrapped, fscale=1)
            # equivalently, but using mean vector length rather than MLE fit
            # kappamu_intercept_fit = inverse_bessel_ratio(directional_stats([mu_intercept_wrapped]).mean_resultant_length) #some potential for bias (rho == 1.0)
            kappamu_intercept_posterior[i] = kappamu_intercept_fit
        except RuntimeError:
            # Handle cases where the MLE optimizer fails (e.g., if all mu's are identical)
            # Assigning a high kappa value, or skipping the draw might be necessary.
            # Here we assign a NaN, which will be filtered later.
            # kappamu_light_posterior[i] = np.nan
            kappamu_light_posterior[i] = 100#high enough
        try:
            kappamu_light_fit, _, _ = vonmises.fit(mu_light_wrapped, fscale=1)
            # kappamu_light_fit = inverse_bessel_ratio(directional_stats([mu_light_wrapped]).mean_resultant_length)
            kappamu_light_posterior[i] = kappamu_light_fit
        except RuntimeError:
            # Handle cases where the MLE optimizer fails (e.g., if all mu's are identical)
            # Assigning a high kappa value, or skipping the draw might be necessary.
            # Here we assign a NaN, which will be filtered later.
            # kappamu_light_posterior[i] = np.nan
            kappamu_light_posterior[i] = 100#high enough

        kappamu_delta_posterior[i] = kappamu_light_posterior[i] -kappamu_intercept_posterior[i] #change in kappa_mu

    # --- 2.5. Summarize the Posterior Distribution of second-level kappa ---
    # Filter out any NaNs if the optimizer failed on certain draws
    kappamu_intercept_posterior = kappamu_intercept_posterior[~np.isnan(kappamu_intercept_posterior)]
    kappamu_light_posterior = kappamu_light_posterior[~np.isnan(kappamu_light_posterior)]

    # Convert the array to an ArviZ InferenceData object for easy summary
    # We need to reshape it back to (n_chain, n_draws) for ArviZ, though it's flat now.
    # For simplicity, we'll use az.summary on the flat array directly.
    # (1, n) so arviz 1.x can read the flat arrays under its default
    # sample_dims = ("chain", "draw").
    kappamu_posterior_idata = {"kappa_mu_intercept": kappamu_intercept_posterior[None, :],
                               "kappa_mu_light": kappamu_light_posterior[None, :],
                               "kappa_mu_delta": kappamu_delta_posterior[None, :],
                               }

    if verbose:
        print("\n--- Summary of Kappa_mu (Light Condition) ---")
        print("This represents the posterior distribution of the concentration of individual mu's.")

    return {
        "id_data": kappamu_posterior_idata,
        "kappamu_intercept_posterior": kappamu_intercept_posterior,
        "kappamu_light_posterior": kappamu_light_posterior,
        "kappamu_delta_posterior": kappamu_delta_posterior,
    }


# ---------------------------------------------------------------------------
# PyMC-graph evaluation helpers for the non-centered fit.
# ---------------------------------------------------------------------------
def graph_free_rvs(idata, model):
    """The fit result with a posterior group containing all of the model's
    free RVs, i.e. evaluable by the PyMC graph.

    The kappa random effects are fitted non-centered: in the PyMC graph the
    free RVs are then ``*_sigma`` (scale) and ``*_offset`` (factor), but
    bambi stores the transformed centered effects and the ``*_sigma`` scales
    and drops the ``*_offset`` factors. ``pm.compute_deterministics`` and
    ``pm.compute_log_likelihood`` both evaluate the graph against every free
    RV, so the missing factors are reconstructed: the graph computes
    ``effect = sigma * offset`` exactly (sigma scalar per draw, offset one
    factor per group), hence ``offset = effect / sigma`` (0 where sigma is
    0).
    """
    posterior = idata.posterior
    missing = [
        rv.name
        for rv in model.backend.model.free_RVs
        if rv.name not in posterior.data_vars
    ]
    if not missing:
        return idata
    aug = idata.copy(deep=False)
    ds = posterior.to_dataset()
    for name in missing:
        sigma = ds[name.replace("_offset", "_sigma")]
        effect = ds[name[: -len("_offset")]]
        # The graph computes effect = sigma * offset (sigma scalar per
        # draw, offset one factor per group), so the offset carries the
        # group dimension of the centered effect.
        e = np.asarray(effect)
        s = np.asarray(sigma).reshape(np.asarray(sigma).shape + (1,) * (e.ndim - np.asarray(sigma).ndim))
        offset = np.divide(e, s, out=np.zeros(e.shape), where=s != 0)
        ds[name] = (tuple(effect.dims), offset)
    aug["posterior"] = ds
    return aug


# ---------------------------------------------------------------------------
# Section 3: supplemental PPC statistics (bootstrap data vs posterior model).
# ---------------------------------------------------------------------------
def hierarchical_bootstrap_group(df, indivs, B=1500, seed=123):
    rng = np.random.default_rng(seed)
    nT = len(indivs)
    pre, post = {}, {}
    for t in indivs:
        y_off = df.loc[(df["Individual"] == t) & (df["Light"] == 0), "HeadingRadians"].dropna().to_numpy(float)
        y_on  = df.loc[(df["Individual"] == t) & (df["Light"] == 1), "HeadingRadians"].dropna().to_numpy(float)
        pre[t]  = np.exp(1j * y_off) if y_off.size else None
        post[t] = np.exp(1j * y_on)  if y_on.size  else None

    mu_off_s, mu_on_s = np.empty(B, float), np.empty(B, float)
    R_off_s, R_on_s   = np.empty(B, float), np.empty(B, float)

    for b in range(B):
        samp_indivs = rng.choice(indivs, size=nT, replace=True)
        u_off, u_on = [], []
        for t in samp_indivs:
            for d, u_list in [(pre[t], u_off), (post[t], u_on)]:
                if d is not None:
                    idx = rng.integers(0, d.size, size=d.size)
                    zbar = d[idx].mean()
                    if np.abs(zbar) > 0: u_list.append(zbar / np.abs(zbar))

        for s_arr, r_arr, u_list in [(mu_off_s, R_off_s, u_off), (mu_on_s, R_on_s, u_on)]:
            if len(u_list) > 0:
                zG = np.mean(u_list)
                s_arr[b], r_arr[b] = circ_mod(np.angle(zG)), np.abs(zG)
            else:
                s_arr[b], r_arr[b] = np.nan, np.nan

    return {"mu_off": mu_off_s, "mu_on": mu_on_s, "R_off": R_off_s, "R_on": R_on_s}


def bootstrap_individual_condition(y, B=1200, seed=0):
    rng = np.random.default_rng(seed)
    y = np.asarray(y, float)
    if y.size == 0: return np.array([]), np.array([])
    z = np.exp(1j * y)
    n = z.size
    mu, R = np.empty(B, float), np.empty(B, float)
    for b in range(B):
        idx = rng.integers(0, n, size=n)
        zbar = z[idx].mean()
        mu[b], R[b] = circ_mod(np.angle(zbar)), np.abs(zbar)
    return mu, R


def ppc_bootstrap_stats(idata_treat, model_treat, df, *, verbose=True):
    """Bootstrap (data) vs posterior (model) PPC statistics for the
    treatment fit: the nonparametric hierarchical group bootstrap, the
    within-condition per-individual bootstraps, and the per-draw per-group
    directions and concentrations from the posterior (graph response
    parameters, cross-checked against the stored log-likelihood).

    Returns a dict of every array ``beetles_plots.plot_ppc`` needs. The
    original cell also computed 95% HDI tables (``mu_tbl``/``R_tbl``) that
    were never displayed or used - the figure recomputes them - so they are
    not returned here.
    """
    # ----------------------------
    # Data Processing
    # ----------------------------
    # df = df_blocks.copy()  # Commented out per request
    df = df.copy()
    df["Individual"] = df["Individual"].astype(str).str.replace(r"^observation088-", "", regex=True)
    # Natural sort individuals (1, 2, ..., 10)
    indivs = sorted(df["Individual"].unique(), key=natural_sort_key)
    nT = len(indivs)

    # ----------------------------
    # Bootstrap routines
    # ----------------------------
    bootG = hierarchical_bootstrap_group(df, indivs, B=2000, seed=123)
    boot_mu_off, boot_mu_on, boot_R_off, boot_R_on = {}, {}, {}, {}

    for j, t in enumerate(indivs):
        y_off = df.loc[(df["Individual"] == t) & (df["Light"] == 0), "HeadingRadians"].dropna().to_numpy(float)
        y_on  = df.loc[(df["Individual"] == t) & (df["Light"] == 1), "HeadingRadians"].dropna().to_numpy(float)

        m_b, r_b = bootstrap_individual_condition(y_off, B=1500, seed=1000 + j)
        boot_mu_off[t], boot_R_off[t] = m_b, r_b
        m_b, r_b = bootstrap_individual_condition(y_on,  B=1500, seed=2000 + j)
        boot_mu_on[t], boot_R_on[t] = m_b, r_b

    # ----------------------------
    # Posterior Stats (idata_treat and model_treat)
    # ----------------------------
    post = idata_treat.posterior
    mu_intercept = _stack_1d(post, "Intercept")
    mu_light_fx  = _stack_1d(post, "Light")
    mu_light_fx  = _stack_1d(post, "Light") if "Light" in post.data_vars else np.zeros_like(mu_intercept)
    mu_pop_off_draw = circ_mod(mu_intercept)
    mu_pop_on_draw  = circ_mod(mu_intercept + mu_light_fx)

    # N.B.: bambi 0.17.2 bug, tracked upstream at bambinos/bambi#1084:
    # Model.predict(kind="response_params") applies the vonmises *mean* link,
    # angle(exp(1j x)), to every component, so its "kappa" is the wrapped
    # linear predictor instead of softplus(linear predictor). Temporary fix
    # until the official bambi fix ships: read the response parameters
    # straight from the PyMC graph with compute_deterministics, which uses
    # each component's actual link.
    #
    # Official (currently broken for the kappa component):
    # preds = model_treat.predict(idata_treat, kind="response_params", inplace=False)
    # mu_da, kappa_da = preds.posterior["mu"], preds.posterior["kappa"]

    _model = model_treat.backend.model
    params = pm.compute_deterministics(
        graph_free_rvs(idata_treat, model_treat),
        var_names=["mu", "kappa"],
        model=_model,
        progressbar=False,
    )
    mu_da, kappa_da = params["mu"], params["kappa"]
    assert np.isfinite(mu_da.values).all()
    assert np.isfinite(kappa_da.values).all()
    assert (kappa_da.values > 0).all()

    # Cross-check the graph parameters against the log-likelihood bambi stored
    # during fitting: they must agree to floating-point precision.
    _obs_ll = idata_treat.log_likelihood["HeadingRadians"]
    _expected = _vonmises.logpdf(
        df["HeadingRadians"].to_numpy(),
        kappa_da.transpose(*_obs_ll.dims).values,
        loc=mu_da.transpose(*_obs_ll.dims).values,
    )
    np.testing.assert_allclose(_expected, _obs_ll.values, rtol=1e-10, atol=1e-10)
    if verbose:
        print("response params from PyMC graph OK; log-likelihood cross-check passed")

    obs_dim = [d for d in mu_da.dims if d not in ("chain", "draw")][0]
    mu_samps = mu_da.stack(sample=("chain", "draw")).transpose("sample", obs_dim).values
    kappa_samps = kappa_da.stack(sample=("chain", "draw")).transpose("sample", obs_dim).values

    indiv_codes = pd.Categorical(df["Individual"], categories=indivs).codes
    light_codes = df["Light"].to_numpy(int)
    group = indiv_codes + nT * light_codes
    G = 2 * nT
    cnt = np.bincount(group, minlength=G).astype(float); cnt[cnt == 0] = np.nan

    z_mu = np.exp(1j * mu_samps)
    sum_re = np.vstack([np.bincount(group, weights=np.real(z_mu[s]), minlength=G) for s in range(len(mu_samps))])
    sum_im = np.vstack([np.bincount(group, weights=np.imag(z_mu[s]), minlength=G) for s in range(len(mu_samps))])
    mu_group_draw = circ_mod(np.arctan2(sum_im / cnt, sum_re / cnt))

    R_obsparam = A_of_kappa(kappa_samps)
    sum_R = np.vstack([np.bincount(group, weights=R_obsparam[s], minlength=G) for s in range(len(kappa_samps))])
    R_group_draw = sum_R / cnt

    mu_off_draw, mu_on_draw = mu_group_draw[:, 0:nT], mu_group_draw[:, nT:2*nT]
    R_off_draw, R_on_draw   = R_group_draw[:, 0:nT], R_group_draw[:, nT:2*nT]
    R_group_off_draw = resultant_length_angles(mu_off_draw, axis=1)
    R_group_on_draw  = resultant_length_angles(mu_on_draw,  axis=1)

    return {
        "indivs": list(indivs),
        "bootG": bootG,
        "boot_mu_off": boot_mu_off, "boot_mu_on": boot_mu_on,
        "boot_R_off": boot_R_off,   "boot_R_on": boot_R_on,
        "mu_pop_off_draw": mu_pop_off_draw, "mu_pop_on_draw": mu_pop_on_draw,
        "mu_off_draw": mu_off_draw, "mu_on_draw": mu_on_draw,
        "R_off_draw": R_off_draw,   "R_on_draw": R_on_draw,
        "R_group_off_draw": R_group_off_draw, "R_group_on_draw": R_group_on_draw,
    }
