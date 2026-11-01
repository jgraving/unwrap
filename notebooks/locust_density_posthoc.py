"""Post-hoc computations for notebooks/locust_density_order.ipynb.

Entry point per figure/evaluation block in the notebook (the matplotlib
"how" lives in ``locust_density_plots.py``, which calls into this module
so the notebook cells stay one call):

- :func:`screening_table`: the per-animal mean step vectors (angle and
  magnitude) and the one-observation-per-animal modelling table (wrapped
  response ``y``, log2-scaled ``t_density``, normalised ``t_order``) that
  the screening model and the 5-model PSIS-LOO sweep are fitted to;
- :func:`fit_table`: the step-level modelling table of the final
  hierarchical fit (categorical ``animal_id``, wrapped ``y``, scaled
  predictors);
- :func:`attach_log_likelihood`: the per-observation von Mises
  log-likelihood of the final fit, read off the PyMC graph (see the note
  below) and stored in the ArviZ-standard ``log_likelihood`` group;
- :func:`model_sweep`: the five candidate mu/kappa-matched
  specifications, seeded fits with pointwise PSIS-LOO, ranked by ΔELPD
  vs the best (order-inclusion flag included);
- :func:`alignment_stats`: posterior predictive directional alignment
  ``A(kappa) * cos(mu)`` on the order/density grids (screening figure:
  uncentered predictors, 68 % HDI, grid of 50; final figure:
  training-mean-centered, 90 % HDI; the PPC-bootstrap figure: 95 % HDI,
  grid of 80);
- :func:`elpd0_stats`: the posterior predictive density at the march
  direction 0 deg, summarised in log space;
- :func:`polar_pdf_stats`: per-cell posterior predictive von Mises PDFs
  (median + 90 % HDI band + circular 90 % interval of mu) for the polar
  small-multiples;
- :func:`coefficient_stats`: the fixed-effect posteriors with 95 % HDIs
  (circular for the mu intercept, linear for the rest);
- :func:`data_cell_bootstrap` / :func:`model_cell_ppc`: the hierarchical
  bootstrap of the alignment statistic (resample animals, then steps
  within animals) for each observed order-density cell, on the data vs
  on the individual-level posterior predictive;
- :func:`individual_ppc`: per-animal alignment, within-animal bootstrap
  (data) vs individual-level posterior predictive (model).

Log-likelihood of the final fit
-------------------------------
``Model.predict(kind="response_params")`` applies the mu *circular-modulo*
link to every component (upstream bambi bug, bambinos/bambi#1084), so its
"kappa" is a wrapped angle and can be negative. The response parameters
are therefore read straight from the PyMC graph with
:func:`pymc.compute_deterministics`, which uses each component's actual
link. One consequence of the RE construction: with the
``transform=Softplus()`` sigma prior the sampled free RV is the log-space
draw (NUTS reports it as ``kappa_1|animal_id_sigma_log_exp``; the
posterior stores it under the friendly name ``kappa_1|animal_id_sigma``),
and ``compute_deterministics`` selects posterior variables by exactly the
model's free-RV names — but the ``kappa_1|animal_id_offset`` factors are
only kept when ``Model.fit(omit_offsets=False)`` is used. The fit in the
notebook therefore passes ``omit_offsets=False`` (as the original
analysis did); without it the graph cannot be evaluated.

The per-observation log-likelihood is then assembled by hand — the von
Mises logpdf in closed form from the response parameters — and stored in
the ``log_likelihood`` group as variable ``y`` (the response name),
which is what ``az.loo(idata, var_name="y")`` reads.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import arviz as az
import bambi as bmb
import pymc as pm
import xarray as xr
from joblib import Parallel, delayed
from scipy.special import i0, i0e, i1

from unwrap.link import BAMBI_LINKS, circular_modulo, inverse_softplus
from unwrap.prior import CircularUniformPrior

# ----------------------------------------------------------------------
# Transform helpers (the original analysis' definitions, kept verbatim).
# ----------------------------------------------------------------------
DENS_MIN = 1.0
DENS_MAX = 64.0


def normalize(x):
    return (np.asarray(x, float) * 2.0) - 1.0


def log2_scale(x, min_size=DENS_MIN, max_size=DENS_MAX):
    x = np.asarray(x, float)
    return np.log2(x / min_size) / np.log2(max_size / min_size)


def _softplus_np(x):
    """Numerically stable softplus for numpy arrays (the link evaluated
    on posterior draws)."""
    x = np.asarray(x)
    return np.where(x > 20, x, np.log1p(np.exp(x)))


def _circular_mean(theta):
    z = np.exp(1j * np.asarray(theta, float))
    return float(circular_modulo(np.angle(z.mean())))


def _A_of_kappa(kappa):
    return i1(kappa) / i0(kappa)


# ----------------------------------------------------------------------
# Posterior linear predictors (shared by the figure entries).
# ----------------------------------------------------------------------
def _training_means(df):
    """bambi-style training means for the centered predictors (each animal
    experienced exactly one condition, so the means run over the rows)."""
    return {
        "t_order": float(np.mean(normalize(df["order"].astype(float)))),
        "t_density": float(np.mean(normalize(log2_scale(df["density"].astype(float))))),
    }


def _posterior_linpred(idata, prefix, X_df):
    """Stacked linear predictor (S, n) for the mu (prefix="") or kappa
    (prefix="kappa_") submodel, from the posterior draws."""
    post = idata.posterior
    if prefix:
        coefs = [v for v in post.data_vars if v.startswith(prefix)]
        key_intercept = f"{prefix}Intercept"
    else:
        coefs = [v for v in post.data_vars if (not v.startswith("kappa_"))]
        key_intercept = "Intercept"

    if key_intercept not in post.data_vars:
        raise ValueError(f"Intercept for prefix '{prefix}' not found in posterior.")

    def flat(name):
        return post[name].stack(sample=("chain", "draw")).values

    n = len(X_df)
    eta = flat(key_intercept)[:, None] * np.ones((1, n), dtype=float)

    for coef in coefs:
        if coef == key_intercept:
            continue
        col = coef.replace(prefix, "")
        if col in X_df.columns:
            eta += flat(coef)[:, None] * X_df[col].to_numpy()[None, :]

    return eta


def _make_X(order_vals, density_vals, idata, train_means=None):
    """Step-level design for the prediction grids. With ``train_means``
    the predictors are centered the way the fit saw them (and the
    interaction is built from the centered columns); without it the raw
    scaled columns are used (the screening fit)."""
    X = pd.DataFrame({"order": order_vals, "density": density_vals})
    X["t_order"] = normalize(X["order"].astype(float))
    X["t_density"] = normalize(log2_scale(X["density"].astype(float)))

    if train_means is not None:
        X["t_order"] = X["t_order"] - train_means["t_order"]
        X["t_density"] = X["t_density"] - train_means["t_density"]

    post = idata.posterior
    has_interaction = (
        "t_density:t_order" in post.data_vars
        or "kappa_t_density:t_order" in post.data_vars
    )
    if has_interaction:
        X["t_density:t_order"] = X["t_density"] * X["t_order"]
    return X


def _levels(df):
    """Descending level vectors and the prediction grids (descending so
    legends/rows run largest -> smallest, as in the original figures)."""
    dens_levels = np.sort(df["density"].unique())[::-1]
    order_levels = np.sort(df["order"].unique())[::-1]
    return dens_levels, order_levels


def _keep_cols(X):
    return ["t_order", "t_density"] + (
        ["t_density:t_order"] if "t_density:t_order" in X.columns else []
    )


# ----------------------------------------------------------------------
# Entry point 1: the one-observation-per-animal screening table
# ----------------------------------------------------------------------
def screening_table(merged_df):
    """Per-animal mean step vectors and the screening modelling table.

    Each animal was observed under exactly one order-density condition,
    so one row per (order, density, ``fname``): the mean (dX, dY) step
    vector, its angle/magnitude, the wrapped response ``y`` and the
    scaled predictors ``t_density`` (log2, 1-64 m^-2 -> [0, 1] then
    mapped to [-1, 1]) and ``t_order`` (0-1 -> [-1, 1]).
    """
    angles_df = (
        merged_df
        .copy()
        .groupby(["order", "density", "fname"], as_index=False)[["dX", "dY"]]
        .mean()
        .rename(columns={"dX": "mean_dX", "dY": "mean_dY"})
    )

    # Angle (no normalisation), magnitude for reference
    angles_df["angle_rad"] = np.arctan2(angles_df["mean_dY"], angles_df["mean_dX"])
    angles_df["angle_deg"] = np.degrees(angles_df["angle_rad"])
    angles_df["magnitude"] = np.hypot(angles_df["mean_dX"], angles_df["mean_dY"])
    angles_df = angles_df.sort_values(["order", "density", "fname"], kind="stable").reset_index(drop=True)

    df = angles_df.copy()

    # Make sure y is wrapped to (-pi, pi]
    df["y"] = ((df["angle_rad"] + np.pi) % (2 * np.pi)) - np.pi
    df["t_density"] = normalize(log2_scale(df["density"].astype(float), DENS_MIN, DENS_MAX))
    df["t_order"] = normalize(df["order"].astype(float))

    # Sanity checks
    assert df["t_density"].between(-1, 1).all()
    assert df["t_order"].between(-1, 1).all()

    return df


# ----------------------------------------------------------------------
# Entry point 2: the step-level modelling table of the final fit
# ----------------------------------------------------------------------
def fit_table(merged_df):
    """The step-level modelling table: categorical ``animal_id``, wrapped
    ``y`` and the scaled predictors (the dead per-animal subsample branch
    of the original cell is dropped: ``n_samples_per_id`` was None)."""
    df = merged_df.copy()

    # Ensure categorical id
    df["animal_id"] = df["animal_id"].astype("category")

    # Wrap y to (-pi, pi]
    df["y"] = ((df["direction"] + np.pi) % (2 * np.pi)) - np.pi
    df["t_density"] = normalize(log2_scale(df["density"].astype(float), DENS_MIN, DENS_MAX))
    df["t_order"] = normalize(df["order"].astype(float))

    print(f"num_samples: {df.shape[0]}")
    return df


# ----------------------------------------------------------------------
# Entry point 3: per-observation log-likelihood of the final fit
# ----------------------------------------------------------------------
def attach_log_likelihood(idata, model, df):
    """Attach the per-observation von Mises log-likelihood to the final
    fit (see the module note on the bambi #1084 quirk and the RE value
    variable naming). Returns the (mutated) fit result."""
    params = pm.compute_deterministics(
        idata.posterior,
        var_names=["mu", "kappa"],
        model=model.backend.model,
        progressbar=False,
    )
    mu = params["mu"]      # dims: (chain, draw, __obs__ or similar)
    kappa = params["kappa"]
    assert (kappa.values > 0).all(), "kappa response params must be positive"
    print("response params from PyMC graph OK (kappa all positive)")

    # Find the observation dimension used by the graph
    obs_dim = [d for d in mu.dims if d not in ("chain", "draw")][0]

    # Observed y aligned to that observation dimension
    y = xr.DataArray(
        df["y"].to_numpy(),
        dims=[obs_dim],
        coords={obs_dim: mu.coords[obs_dim]},
    )

    # von Mises logpdf: k*cos(y-mu) - log(2*pi) - log(I0(k))
    # stable log(I0(k)) using i0e: I0(k) = i0e(k)*exp(k)
    logI0 = np.log(i0e(kappa)) + kappa
    loglik = kappa * np.cos(y - mu) - np.log(2 * np.pi) - logI0

    # Store in the ArviZ-standard location (variable named after the response)
    idata["log_likelihood"] = xr.Dataset({"y": loglik})
    return idata


# ----------------------------------------------------------------------
# Entry point 4: the five-candidate model sweep (PSIS-LOO)
# ----------------------------------------------------------------------
#: Candidate specifications (both submodels carry the same terms).
MODELS = [
    ("Intercept-only", []),
    ("Density", ["t_density"]),
    ("Order", ["t_order"]),
    ("Density + Order", ["t_density", "t_order"]),
    ("Density × Order", ["t_density", "t_order", "t_density:t_order"]),
]

MU_INTERCEPT_PRIOR = CircularUniformPrior()
MU_SLOPE_PRIOR = bmb.Prior("Normal", mu=0.0, sigma=0.5 * np.pi)
K_INTERCEPT_PRIOR = bmb.Prior("Normal", mu=inverse_softplus(2.0), sigma=1.0)
K_SLOPE_PRIOR = bmb.Prior("Normal", mu=0.0, sigma=1.0)


def _rhs_from_terms(terms):
    return "1" if not terms else ("1 + " + " + ".join(terms))


def _build_formulas(terms):
    mu_str = f"y ~ {_rhs_from_terms(terms)}"
    k_str = f"kappa ~ {_rhs_from_terms(terms)}"
    return mu_str, k_str, bmb.Formula(mu_str, k_str)


def _build_priors(terms):
    pri = {}
    pri["Intercept"] = MU_INTERCEPT_PRIOR
    for term in (terms or []):
        pri[term] = MU_SLOPE_PRIOR

    pri["kappa"] = {"Intercept": K_INTERCEPT_PRIOR}
    for term in (terms or []):
        pri["kappa"][term] = K_SLOPE_PRIOR
    return pri


def _elpd_diff_from_pointwise(model_loo_i, base_loo_i):
    diff_i = (np.asarray(model_loo_i) - np.asarray(base_loo_i)).reshape(-1)
    n = diff_i.size
    elpd_diff = float(diff_i.sum())
    se_diff = float(np.sqrt(n * np.var(diff_i, ddof=1))) if n > 1 else float("nan")
    return elpd_diff, se_diff


def model_sweep(screen_df, draws=1000, seed=123):
    """Fit the five candidate specifications to the screening table, rank
    them by pointwise ΔELPD (PSIS-LOO) versus the best, and return the
    sorted results table (with the order-inclusion flag the figures
    use). Prints a per-model fit header and divergence count."""
    meta_rows = []
    loo_i_by_name = {}

    for pretty_name, terms in MODELS:
        print(f"\n=== Fitting {pretty_name} ===")

        mu_str, k_str, formula = _build_formulas(terms)
        priors = _build_priors(terms)

        mdl = bmb.Model(
            formula=formula,
            data=screen_df,
            family="vonmises",
            priors=priors,
            link={"mu": BAMBI_LINKS["circular_modulo"], "kappa": BAMBI_LINKS["softplus"]},
        )

        fit_idata = mdl.fit(
            draws=draws,
            chains=4,
            target_accept=0.9,
            random_seed=seed,
            progressbar=False,
        )
        # The PyMC 6 / bambi 0.21 replacement for
        # idata_kwargs={"log_likelihood": True}.
        fit_idata = pm.compute_log_likelihood(
            fit_idata, model=mdl.backend.model, progressbar=False
        )
        print("Divergences:", int(fit_idata.sample_stats.diverging.sum()))

        loo = az.loo(fit_idata, pointwise=True)

        loo_i = np.asarray(loo.elpd_i).reshape(-1)
        loo_i_by_name[pretty_name] = loo_i

        meta_rows.append({
            "pretty_name": pretty_name,
            "terms": terms,
            "mu_formula": mu_str,
            "kappa_formula": k_str,
            "elpd": float(loo.elpd),
            "se_elpd": float(loo.se),
            "p_eff": float(loo.p),
            "pareto_k_max": float(np.max(np.asarray(loo.pareto_k))),
        })

    meta_df = pd.DataFrame(meta_rows)

    # Manual pointwise ΔELPD vs BEST (highest elpd):
    # ΔELPD = model - best (best = 0, others <= 0)
    best_name = meta_df.loc[meta_df["elpd"].idxmax(), "pretty_name"]
    base_loo_i = loo_i_by_name[best_name]

    meta_df["best_model"] = best_name

    deltas = []
    dse = []
    for name in meta_df["pretty_name"]:
        d, se = _elpd_diff_from_pointwise(loo_i_by_name[name], base_loo_i)
        deltas.append(d)
        dse.append(se)

    meta_df["delta_elpd"] = deltas
    meta_df["se_delta_elpd"] = dse

    # Order-inclusion flag for the figure groupings
    meta_df["has_order"] = meta_df["terms"].apply(
        lambda t: any("t_order" in s for s in t)
    )

    # Sort best-first (delta 0 at top)
    meta_df = meta_df.sort_values("delta_elpd", ascending=False).reset_index(drop=True)

    print(f"\nBest model: {best_name}")
    return meta_df


# ----------------------------------------------------------------------
# Entry point 5: posterior predictive directional alignment
# ----------------------------------------------------------------------
def alignment_stats(idata, df, hdi_prob=0.9, centered=True, grid_size=50):
    """Posterior predictive alignment ``A(kappa) * cos(mu)`` on the
    order/density grids: one curve per density level (panel i, over order)
    and one per order level (panel ii, over density), with posterior
    mean and HDI. ``centered=False`` reproduces the screening figure
    (uncentered predictors); the final-fit figures center on the
    training means. Returns the predictions plus grid bookkeeping."""
    dens_levels, order_levels = _levels(df)
    train_means = _training_means(df) if centered else None

    order_grid = np.linspace(order_levels.min(), order_levels.max(), grid_size)
    density_grid = np.geomspace(max(dens_levels.min(), 1e-6), dens_levels.max(), grid_size)

    def _at(X):
        eta_mu = _posterior_linpred(idata, prefix="", X_df=X)
        eta_kappa = _posterior_linpred(idata, prefix="kappa_", X_df=X)

        mu = circular_modulo(eta_mu)
        kappa = _softplus_np(eta_kappa)

        align = _A_of_kappa(kappa) * np.cos(mu)

        # az.hdi treats the LAST axis of a raw 2D array as the samples
        # axis, so hand it (n_points, S).
        hdi = az.hdi(align.T, prob=hdi_prob, skipna=True)
        return {"mean": align.mean(axis=0), "low": hdi[:, 0], "high": hdi[:, 1]}

    pred_i = {}
    for d in dens_levels:
        X = _make_X(order_grid, np.full_like(order_grid, d, dtype=float),
                    idata, train_means)
        pred_i[d] = _at(X[_keep_cols(X)])

    pred_ii = {}
    for o in order_levels:
        X = _make_X(np.full_like(density_grid, o, dtype=float), density_grid,
                    idata, train_means)
        pred_ii[o] = _at(X[_keep_cols(X)])

    return {
        "pred_i": pred_i,
        "pred_ii": pred_ii,
        "order_grid": order_grid,
        "density_grid": density_grid,
        "dens_levels": dens_levels,
        "order_levels": order_levels,
    }


# ----------------------------------------------------------------------
# Entry point 6: posterior predictive density at the march direction
# ----------------------------------------------------------------------
def elpd0_stats(idata, df, hdi_prob=0.9, grid_size=50):
    """Posterior predictive density of a step heading exactly at the march
    direction (0 deg): for each grid point the von Mises logpdf at 0 under
    each posterior draw, summarised in log space (mean/HDI) and
    exponentiated for plotting."""
    dens_levels, order_levels = _levels(df)
    train_means = _training_means(df)

    order_grid = np.linspace(order_levels.min(), order_levels.max(), grid_size)
    density_grid = np.geomspace(max(dens_levels.min(), 1e-6), dens_levels.max(), grid_size)

    def _at(X):
        eta_mu = _posterior_linpred(idata, prefix="", X_df=X)
        eta_kappa = _posterior_linpred(idata, prefix="kappa_", X_df=X)

        mu = circular_modulo(eta_mu)
        kappa = _softplus_np(eta_kappa)

        # log predictive density at y=0 for each posterior draw
        logp0 = vonmises_logpdf(0.0, kappa, loc=mu)  # (S, n_points)

        mean_log = np.mean(logp0, axis=0)
        # az.hdi treats the LAST axis of a raw 2D array as the samples axis.
        hdi_log = az.hdi(logp0.T, prob=hdi_prob, skipna=True)

        return {
            "mean": np.exp(mean_log),
            "low": np.exp(hdi_log[:, 0]),
            "high": np.exp(hdi_log[:, 1]),
        }

    pred_i = {}
    for d in dens_levels:
        X = _make_X(order_grid, np.full_like(order_grid, d, dtype=float),
                    idata, train_means)
        pred_i[d] = _at(X[_keep_cols(X)])

    pred_ii = {}
    for o in order_levels:
        X = _make_X(np.full_like(density_grid, o, dtype=float), density_grid,
                    idata, train_means)
        pred_ii[o] = _at(X[_keep_cols(X)])

    return {
        "pred_i": pred_i,
        "pred_ii": pred_ii,
        "order_grid": order_grid,
        "density_grid": density_grid,
        "dens_levels": dens_levels,
        "order_levels": order_levels,
    }


def _vonmises_logpdf_scalar(y, kappa, loc):
    """von Mises logpdf with stable log I0, numpy broadcasting."""
    logI0 = np.log(i0e(kappa)) + kappa
    return kappa * np.cos(y - loc) - np.log(2 * np.pi) - logI0


def vonmises_logpdf(y, kappa, loc):
    return _vonmises_logpdf_scalar(y, np.asarray(kappa, float), np.asarray(loc, float))


# ----------------------------------------------------------------------
# Entry point 7: per-cell posterior predictive von Mises PDFs (polar)
# ----------------------------------------------------------------------
def polar_pdf_stats(idata, df, hdi_prob=0.90, max_draws=1200, seed=123,
                    theta=None):
    """Per order-density cell: the posterior predictive von Mises PDF on a
    900-point circle (median + HDI band), the circular mean of the mu
    draws and its 90 % interval drawn as an arc on the inner ring.
    Posterior draws are subsampled (seeded) to ``max_draws`` cells
    sequentially, exactly as the original cell did."""
    dens_levels, order_levels = _levels(df)
    train_means = _training_means(df)
    if theta is None:
        theta = np.linspace(-np.pi, np.pi, 900, endpoint=False)

    rng = np.random.default_rng(seed)

    summ = {}
    peak_max = 0.0
    for o in order_levels:
        for d in dens_levels:
            X = _make_X(np.array([o], dtype=float), np.array([d], dtype=float),
                        idata, train_means)
            X = X[_keep_cols(X)]

            eta_mu = _posterior_linpred(idata, prefix="", X_df=X)
            eta_kappa = _posterior_linpred(idata, prefix="kappa_", X_df=X)

            mu = circular_modulo(eta_mu[:, 0])
            kappa = _softplus_np(eta_kappa[:, 0])

            mu, kappa = _subsample_draws(mu, kappa, max_draws, rng)

            pdfs = vonmises_pdf(theta[None, :], mu[:, None], kappa[:, None])

            med = np.median(pdfs, axis=0)
            # az.hdi treats the LAST axis of a raw 2D array as the samples axis.
            hdi = az.hdi(pdfs.T, prob=hdi_prob, skipna=True)
            lo, hi = hdi[:, 0], hdi[:, 1]

            mu0, mu_lo, mu_hi = _mu_interval_central(mu, prob=hdi_prob)
            summ[(float(o), float(d))] = dict(
                med=med, lo=lo, hi=hi, mu0=mu0, mu_lo=mu_lo, mu_hi=mu_hi
            )
            peak_max = max(peak_max, hi.max())

    return {
        "summ": summ,
        "theta": theta,
        "order_levels": order_levels,
        "dens_levels": dens_levels,
        "peak_max": peak_max,
    }


def vonmises_pdf(theta, mu, kappa):
    return (
        np.exp(kappa * np.cos(theta - mu))
        / (2 * np.pi * np.i0(kappa))
    )


def _subsample_draws(mu, kappa, max_draws, rng):
    S = mu.shape[0]
    if S <= max_draws:
        return mu, kappa
    idx = rng.choice(S, size=max_draws, replace=False)
    return mu[idx], kappa[idx]


def _mu_interval_central(mu_draws, prob=0.9):
    """Practical circular interval: unwrap around the circular mean, then
    take the central interval in the unwrapped coordinates."""
    mu0 = np.angle(np.mean(np.exp(1j * mu_draws)))
    d = np.angle(np.exp(1j * (mu_draws - mu0)))  # deltas in (-pi, pi]
    a = (1 - prob) / 2 * 100
    b = (1 + prob) / 2 * 100
    lo_d, hi_d = np.percentile(d, [a, b])
    lo = circular_modulo(mu0 + lo_d)
    hi = circular_modulo(mu0 + hi_d)
    return mu0, lo, hi


# ----------------------------------------------------------------------
# Entry point 8: coefficient posteriors (fixed effects, mu vs kappa)
# ----------------------------------------------------------------------
#: (pretty name, mu variable, kappa variable)
TERMS = [
    ("Intercept", "Intercept", "kappa_Intercept"),
    ("Density", "t_density", "kappa_t_density"),
    ("Order", "t_order", "kappa_t_order"),
    ("Density × Order", "t_density:t_order", "kappa_t_density:t_order"),
]


def circular_hdi(theta, hdi_prob=0.95):
    theta = np.asarray(theta, float)
    mu_c = _circular_mean(theta)
    centered = circular_modulo(theta - mu_c)
    h = az.hdi(centered, prob=hdi_prob, skipna=True)
    low = float(circular_modulo(float(h[0]) + mu_c))
    high = float(circular_modulo(float(h[1]) + mu_c))
    return low, high, mu_c


def _linear_hdi_mean(x, hdi_prob=0.95):
    x = np.asarray(x, float)
    h = az.hdi(x, prob=hdi_prob, skipna=True)
    return float(x.mean()), float(h[0]), float(h[1])


def _circ_dist(a, b):
    return float(circular_modulo(b - a))


def coefficient_stats(idata, hdi_prob=0.95):
    """Fixed-effect posteriors (mu and kappa, pre-link scale) with HDIs:
    circular for the mu intercept, linear mean+HDI otherwise; returns the
    per-term point estimates and (asymmetric) error-bar widths in the
    order of ``TERMS``."""
    post = idata.posterior

    rows = []
    for pretty, mu_name, k_name in TERMS:
        if mu_name not in post.data_vars:
            raise KeyError(f"Missing '{mu_name}' in idata.posterior.")
        if k_name not in post.data_vars:
            raise KeyError(f"Missing '{k_name}' in idata.posterior.")
        mu_draws = post[mu_name].stack(sample=("chain", "draw")).values.astype(float)
        k_draws = post[k_name].stack(sample=("chain", "draw")).values.astype(float)

        if mu_name == "Intercept":
            mu_lo, mu_hi, mu_m = circular_hdi(mu_draws, hdi_prob=hdi_prob)
            rows.append(("mu", pretty, float(mu_m), float(mu_lo), float(mu_hi)))
        else:
            mu_m, mu_lo, mu_hi = _linear_hdi_mean(mu_draws, hdi_prob)
            rows.append(("mu", pretty, float(mu_m), float(mu_lo), float(mu_hi)))

        k_m, k_lo, k_hi = _linear_hdi_mean(k_draws, hdi_prob)
        rows.append(("kappa", pretty, float(k_m), float(k_lo), float(k_hi)))

    terms_order = [t[0] for t in TERMS]

    mu_stats = {term: (mean, lo, hi) for (sub, term, mean, lo, hi) in rows if sub == "mu"}
    k_stats = {term: (mean, lo, hi) for (sub, term, mean, lo, hi) in rows if sub == "kappa"}

    mu_mean = np.array([mu_stats[t][0] for t in terms_order], float)
    mu_lo = np.array([mu_stats[t][1] for t in terms_order], float)
    mu_hi = np.array([mu_stats[t][2] for t in terms_order], float)

    k_mean = np.array([k_stats[t][0] for t in terms_order], float)
    k_lo = np.array([k_stats[t][1] for t in terms_order], float)
    k_hi = np.array([k_stats[t][2] for t in terms_order], float)

    mu_xerr = np.zeros((2, len(terms_order)), float)
    for i, term in enumerate(terms_order):
        if term == "Intercept":
            mu_xerr[:, i] = [
                abs(_circ_dist(mu_lo[i], mu_mean[i])),
                abs(_circ_dist(mu_mean[i], mu_hi[i])),
            ]
        else:
            mu_xerr[:, i] = [mu_mean[i] - mu_lo[i], mu_hi[i] - mu_mean[i]]

    k_xerr = np.vstack([k_mean - k_lo, k_hi - k_mean])

    return {
        "terms_order": terms_order,
        "mu_mean": mu_mean,
        "mu_xerr": mu_xerr,
        "k_mean": k_mean,
        "k_xerr": k_xerr,
    }


# ----------------------------------------------------------------------
# Entry point 9: empirical hierarchical bootstrap per treatment cell
# ----------------------------------------------------------------------
def data_cell_bootstrap(df, B=1000, ci=0.95, seed=123, n_jobs=-1, verbose=0):
    """Hierarchical bootstrap of the alignment statistic per observed
    order-density cell: resample animals (with replacement), then resample
    steps within each resampled animal; the statistic is
    ``Re(mean(exp(i*y)))`` over the resampled steps, averaged over the
    resampled animals. Returns the per-cell mean and HDI."""
    rng_master = np.random.default_rng(seed)

    records = []
    for aid, g in df.groupby("animal_id", observed=True):
        dens = float(g["density"].iloc[0])
        ordv = float(g["order"].iloc[0])
        u = np.exp(1j * g["y"].to_numpy())
        records.append({"animal_id": aid, "density": dens, "order": ordv, "u": u})
    animal_tbl = pd.DataFrame(records)

    cell_groups = {}
    for (dens, ordv), sub in animal_tbl.groupby(["density", "order"], observed=True):
        cell_groups[(float(dens), float(ordv))] = sub["animal_id"].to_numpy()

    u_by_animal = dict(zip(animal_tbl["animal_id"], animal_tbl["u"]))

    def _hdi_1d(x, prob):
        h = az.hdi(x, prob=prob, skipna=True)
        h = np.asarray(h)
        if h.ndim == 2:
            return float(h[0, 0]), float(h[0, 1])
        return float(h[0]), float(h[1])

    def _bootstrap_one_cell(cell_key, animals, cell_seed):
        rng = np.random.default_rng(cell_seed)
        animals = np.asarray(animals)
        nA = len(animals)
        boots = np.empty(B, dtype=float)

        for b in range(B):
            samp_animals = animals[rng.integers(0, nA, size=nA)]
            z_sum = 0.0 + 0.0j

            for aid in samp_animals:
                u = u_by_animal[aid]
                idx = rng.integers(0, len(u), size=len(u))
                z_sum += u[idx].mean()

            boots[b] = np.real(z_sum / nA)

        low, high = _hdi_1d(boots, ci)

        return {
            "density": cell_key[0],
            "order": cell_key[1],
            "n_animals": int(nA),
            "emp_mean": float(boots.mean()),
            "emp_low": float(low),
            "emp_high": float(high),
        }

    cells = list(cell_groups.items())
    cell_seeds = rng_master.integers(0, 2**32 - 1, size=len(cells), dtype=np.uint32)

    out = Parallel(n_jobs=n_jobs, prefer="threads", verbose=verbose)(
        delayed(_bootstrap_one_cell)(cell_key, animals, int(cell_seeds[i]))
        for i, (cell_key, animals) in enumerate(cells)
    )

    return pd.DataFrame(out)


# ----------------------------------------------------------------------
# Entry point 10: model-side hierarchical bootstrap (cell-level PPC)
# ----------------------------------------------------------------------
def model_cell_ppc(df, idata, B=600, ci=0.95, seed=2025, n_jobs=-1,
                   bootstrap_steps=True):
    """Posterior predictive version of :func:`data_cell_bootstrap`: for
    each bootstrap draw, simulate the same number of steps per animal from
    the individual-level likelihood (fixed effects + that animal's random
    effect) and apply the identical hierarchical bootstrap to the
    simulated data. Returns the per-cell model mean and HDI."""
    rng_master = np.random.default_rng(seed)
    post = idata.posterior

    animal_dim = None
    for d in post.dims:
        if "animal_id" in d:
            animal_dim = d
            break
    if animal_dim is None:
        raise KeyError(f"Could not find an animal dimension containing 'animal_id'. "
                       f"Found dims: {list(post.dims)}")

    def stack_samples_1d(name, default_zero=False):
        if name not in post.data_vars:
            if default_zero:
                S = int(post.sizes["chain"] * post.sizes["draw"])
                return np.zeros(S, dtype=float)
            raise KeyError(f"Missing posterior variable: {name}")
        return post[name].stack(sample=("chain", "draw")).values.astype(float)

    def stack_samples_re(name):
        if name not in post.data_vars:
            raise KeyError(f"Missing posterior random effect: {name}")
        da = post[name].stack(sample=("chain", "draw")).transpose("sample", animal_dim)
        return da.values.astype(float)

    b0 = stack_samples_1d("Intercept")
    b_td = stack_samples_1d("t_density", default_zero=True)
    b_to = stack_samples_1d("t_order", default_zero=True)
    b_int = stack_samples_1d("t_density:t_order", default_zero=True)

    k0 = stack_samples_1d("kappa_Intercept")
    k_td = stack_samples_1d("kappa_t_density", default_zero=True)
    k_to = stack_samples_1d("kappa_t_order", default_zero=True)
    k_int = stack_samples_1d("kappa_t_density:t_order", default_zero=True)

    if "1|animal_id" in post.data_vars:
        mu_re = stack_samples_re("1|animal_id")
    else:
        # pooled-mu spec: direction has no per-animal effect
        S_tot = int(post.sizes["chain"] * post.sizes["draw"])
        mu_re = np.zeros((S_tot, int(post.sizes[animal_dim])), dtype=float)
    kappa_re = stack_samples_re("kappa_1|animal_id")

    S_total = b0.shape[0]
    post_animal_ids = post.coords[animal_dim].values.astype(str)
    post_index = {aid: i for i, aid in enumerate(post_animal_ids)}

    df_local = df.copy()
    df_local["animal_id"] = df_local["animal_id"].astype(str)

    animal_tbl = (
        df_local.groupby("animal_id", observed=True)
        .agg(density=("density", "first"),
             order=("order", "first"),
             n_steps=("y", "size"))
        .reset_index()
    )

    train_means = _training_means(df_local)
    animal_tbl["t_order"] = normalize(animal_tbl["order"].astype(float)) - train_means["t_order"]
    animal_tbl["t_density"] = normalize(log2_scale(animal_tbl["density"].astype(float))) - train_means["t_density"]
    animal_tbl["t_int"] = animal_tbl["t_order"] * animal_tbl["t_density"]

    animal_tbl["post_idx"] = animal_tbl["animal_id"].map(post_index)
    if animal_tbl["post_idx"].isna().any():
        bad = animal_tbl.loc[animal_tbl["post_idx"].isna(), "animal_id"].tolist()[:10]
        raise ValueError(
            "Some df animal_id levels not found in posterior coords. "
            f"Examples: {bad}\nPosterior coord example: {post_animal_ids[:10].tolist()}"
        )
    animal_tbl["post_idx"] = animal_tbl["post_idx"].astype(int)

    dens_arr = animal_tbl["density"].to_numpy(float)
    order_arr = animal_tbl["order"].to_numpy(float)
    nsteps_arr = animal_tbl["n_steps"].to_numpy(int)
    torder_arr = animal_tbl["t_order"].to_numpy(float)
    tdens_arr = animal_tbl["t_density"].to_numpy(float)
    tint_arr = animal_tbl["t_int"].to_numpy(float)
    postidx_arr = animal_tbl["post_idx"].to_numpy(int)

    cell_to_rows = {}
    for i in range(len(animal_tbl)):
        key = (float(dens_arr[i]), float(order_arr[i]))
        cell_to_rows.setdefault(key, []).append(i)

    def _hdi_1d(x, prob):
        h = az.hdi(x, prob=prob, skipna=True)
        h = np.asarray(h)
        if h.ndim == 2:
            return float(h[0, 0]), float(h[0, 1])
        return float(h[0]), float(h[1])

    def _one_cell(cell_key, row_idx_list, cell_seed):
        rng = np.random.default_rng(cell_seed)
        rows = np.asarray(row_idx_list, dtype=int)
        nA = len(rows)
        boots = np.empty(B, dtype=float)

        for b in range(B):
            samp_rows = rows[rng.integers(0, nA, size=nA)]
            s = int(rng.integers(0, S_total))

            z_sum = 0.0 + 0.0j
            for rr in samp_rows:
                j = int(postidx_arr[rr])
                n_steps = int(nsteps_arr[rr])

                eta_mu = (
                    b0[s]
                    + b_td[s] * tdens_arr[rr]
                    + b_to[s] * torder_arr[rr]
                    + b_int[s] * tint_arr[rr]
                    + mu_re[s, j]
                )
                eta_kappa = (
                    k0[s]
                    + k_td[s] * tdens_arr[rr]
                    + k_to[s] * torder_arr[rr]
                    + k_int[s] * tint_arr[rr]
                    + kappa_re[s, j]
                )

                mu = circular_modulo(eta_mu)
                kappa = _softplus_np(eta_kappa)

                y_rep = rng.vonmises(mu, kappa, size=n_steps)
                u = np.exp(1j * y_rep)

                if bootstrap_steps:
                    idx = rng.integers(0, n_steps, size=n_steps)
                    z_sum += u[idx].mean()
                else:
                    z_sum += u.mean()

            boots[b] = float(np.real(z_sum / nA))

        low, high = _hdi_1d(boots, ci)
        return {
            "density": cell_key[0],
            "order": cell_key[1],
            "n_animals": int(nA),
            "model_mean": float(boots.mean()),
            "model_low": float(low),
            "model_high": float(high),
        }

    cells = list(cell_to_rows.items())
    cell_seeds = rng_master.integers(0, 2**32 - 1, size=len(cells), dtype=np.uint32)

    out = Parallel(n_jobs=n_jobs, prefer="threads", verbose=0)(
        delayed(_one_cell)(cell_key, row_idx_list, int(cell_seeds[i]))
        for i, (cell_key, row_idx_list) in enumerate(cells)
    )

    return pd.DataFrame(out)


# ----------------------------------------------------------------------
# Entry point 11: per-animal PPC (data bootstrap vs model PPC)
# ----------------------------------------------------------------------
def individual_ppc(df, idata, b_boot=800, s_ppc=400, boot_seed=123,
                   ppc_seed=2025, hdp=0.95):
    """Per animal, within one order-density condition: the within-animal
    bootstrap of the alignment statistic (data) against the individual
    level posterior predictive (model: simulate the same number of steps
    from that animal's own mu/kappa). Returns a one-row-per-animal table;
    the 6x6 panel grid is drawn by the plots module."""
    post = idata.posterior

    animal_dim = None
    for d in post.dims:
        if "animal_id" in d:
            animal_dim = d
            break
    if animal_dim is None:
        raise KeyError(f"Could not find an animal dimension containing 'animal_id'. "
                       f"Found dims: {list(post.dims)}")

    def stack_samples_1d(name):
        return post[name].stack(sample=("chain", "draw")).values

    def stack_samples_re(name):
        da = post[name].stack(sample=("chain", "draw")).transpose("sample", animal_dim)
        return da.values

    # fixed effects (mu)
    b0 = stack_samples_1d("Intercept").astype(float)
    b_td = stack_samples_1d("t_density").astype(float)
    b_to = stack_samples_1d("t_order").astype(float)
    b_int = stack_samples_1d("t_density:t_order").astype(float)

    # fixed effects (kappa pre-softplus)
    k0 = stack_samples_1d("kappa_Intercept").astype(float)
    k_td = stack_samples_1d("kappa_t_density").astype(float)
    k_to = stack_samples_1d("kappa_t_order").astype(float)
    k_int = stack_samples_1d("kappa_t_density:t_order").astype(float)

    # random effects (pooled-mu spec: direction has no per-animal effect)
    if "1|animal_id" in post.data_vars:
        mu_re = stack_samples_re("1|animal_id")
    else:
        mu_re = np.zeros((b0.shape[0], int(post.sizes[animal_dim])), dtype=float)
    kappa_re = stack_samples_re("kappa_1|animal_id")

    S_total = b0.shape[0]
    n_animals_post = mu_re.shape[1]
    print(f"[OK] posterior draws S={S_total}, animals={n_animals_post}")
    print(f"[OK] mu_re shape={mu_re.shape} kappa_re shape={kappa_re.shape}")

    df_local = df.copy()
    df_local["animal_id"] = df_local["animal_id"].astype(str)

    y_by_animal = {aid: g["y"].to_numpy() for aid, g in df_local.groupby("animal_id", observed=True)}

    animal_tbl = (
        df_local.groupby("animal_id", observed=True)
        .agg(density=("density", "first"), order=("order", "first"), n_steps=("y", "size"))
        .reset_index()
    )

    train_means = _training_means(df_local)
    animal_tbl["t_order"] = normalize(animal_tbl["order"].astype(float)) - train_means["t_order"]
    animal_tbl["t_density"] = normalize(log2_scale(animal_tbl["density"].astype(float))) - train_means["t_density"]
    animal_tbl["t_int"] = animal_tbl["t_order"] * animal_tbl["t_density"]

    # map animal_id -> posterior index using the coords on the animal dim
    post_animal_ids = post.coords[animal_dim].values.astype(str)
    post_index = {aid: i for i, aid in enumerate(post_animal_ids)}

    animal_tbl["post_idx"] = animal_tbl["animal_id"].map(post_index)
    if animal_tbl["post_idx"].isna().any():
        bad = animal_tbl.loc[animal_tbl["post_idx"].isna(), "animal_id"].tolist()[:10]
        raise ValueError(
            "Some df animal_id levels not found in posterior coords. "
            f"Examples: {bad}\nPosterior coord example: {post_animal_ids[:10].tolist()}"
        )
    animal_tbl["post_idx"] = animal_tbl["post_idx"].astype(int)

    def alignment_stat(y):
        y = np.asarray(y, dtype=float)
        return float(np.real(np.exp(1j * y).mean()))

    def bootstrap_alignment(y, B, seed):
        rng = np.random.default_rng(seed)
        u = np.exp(1j * np.asarray(y, dtype=float))
        n = len(u)
        boots = np.empty(B, dtype=float)
        for b in range(B):
            idx = rng.integers(0, n, size=n)
            boots[b] = float(np.real(u[idx].mean()))
        return boots

    rng_master = np.random.default_rng(ppc_seed)
    out_rows = []

    for row in animal_tbl.itertuples(index=False):
        aid = row.animal_id
        j = int(row.post_idx)
        n_steps = int(row.n_steps)
        y_obs = y_by_animal[aid]

        # ---- DATA bootstrap ----
        boots = bootstrap_alignment(y_obs, B=b_boot, seed=boot_seed + j)
        data_mean = float(boots.mean())
        data_low, data_high = az.hdi(boots, prob=hdp, skipna=True)
        data_low, data_high = float(data_low), float(data_high)

        # ---- MODEL PPC ----
        draw_idx = rng_master.choice(S_total, size=min(s_ppc, S_total), replace=False)

        eta_mu = (
            b0[draw_idx]
            + b_td[draw_idx] * float(row.t_density)
            + b_to[draw_idx] * float(row.t_order)
            + b_int[draw_idx] * float(row.t_int)
            + mu_re[draw_idx, j]
        )

        eta_kappa = (
            k0[draw_idx]
            + k_td[draw_idx] * float(row.t_density)
            + k_to[draw_idx] * float(row.t_order)
            + k_int[draw_idx] * float(row.t_int)
            + kappa_re[draw_idx, j]
        )

        mu = circular_modulo(eta_mu)
        kappa = _softplus_np(eta_kappa)

        rng = np.random.default_rng(ppc_seed + 10_000 + j)
        reps = np.empty(len(draw_idx), dtype=float)
        for i in range(len(draw_idx)):
            y_rep = rng.vonmises(mu[i], kappa[i], size=n_steps)
            reps[i] = alignment_stat(y_rep)

        model_mean = float(reps.mean())
        model_low, model_high = az.hdi(reps, prob=hdp, skipna=True)

        out_rows.append({
            "animal_id": aid,
            "density": float(row.density),
            "order": float(row.order),
            "data_mean": data_mean,
            "data_low": data_low,
            "data_high": data_high,
            "model_mean": model_mean,
            "model_low": float(model_low),
            "model_high": float(model_high),
        })

    return pd.DataFrame(out_rows)
