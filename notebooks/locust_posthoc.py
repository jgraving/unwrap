"""Post-hoc computations for notebooks/locust_direction_mixture_model.ipynb.

Entry point per figure in the notebook (the matplotlib "how" lives in
locust_direction_plots.py, which calls into this module so the notebook
cells stay one call):

- :func:`loo_elpd_rows`: PSIS-LOO ELPD and SE for the bimodal and unimodal fits;
- :func:`elpd_by_individual`: pointwise LOO contributions aggregated per animal;
- :func:`population_density` / :func:`bootstrap_histogram_stats` /
  :func:`bootstrap_von_mises_kde_stats` (combined by
  :func:`population_validation_stats`): the three population-level estimates
  of the directional density behind the validation figures;
- :func:`individual_ppc_stats`: per-individual posterior predictive bands and
  bootstrap KDE bands of the 52-panel PPC grid;
- :func:`hypothesis_density_stats` / :func:`hypothesis_elpd_stats`: the
  population mixture's posterior predictive density and ELPD at the
  hypothesized directions (forward vs perpendicular);
- :func:`hypothesis_individual_stats` /
  :func:`hypothesis_elpd_individual_stats`: the same tests per animal, plus a
  pooled row;
- :func:`population_alignment_stats` / :func:`individual_alignment_stats`:
  weight-averaged directional-alignment scores around 0 and the axial
  +/-pi/2.

Animal indexing: every entry point indexes individuals densely from 0 - the
notebook's data cell renumbers ``df["animal_id"]`` so that holds.
"""

from concurrent.futures import ProcessPoolExecutor, as_completed
import os

import arviz as az
import numpy as np
import pandas as pd
import pymc as pm
import pytensor.tensor as pt
from joblib import Parallel, delayed
from numba import jit
from scipy.special import i0
from scipy.stats import vonmises

# Magic constants from the original analysis (kept verbatim).
N_BINS = 16           # histogram bins over the full circle
N_BOOTSTRAPS = 1000   # population bootstrap reps
N_POINTS = 360        # population density / KDE evaluation grid
KAPPA_KDE_POP = 30.0  # fixed KDE kernel concentration (population KDE)
KAPPA_KDE = 15.0      # fixed KDE kernel concentration (individual KDE)
N_BOOT_KDE = 1000     # bootstrap reps per individual
BASE_SEED = 123       # individual KDE seeds: BASE_SEED + 10_000 + idx
BOOTSTRAP_SEED = 123  # population bootstrap seed (deterministic pipeline)
N_ANGLES = 360        # individual PPC evaluation grid (endpoint=False)
EPS = 1e-12
Z = 1.96
PERP_COMBINE = "avg"  # combine p(+/-pi/2) by averaging (log scale)
USE_MEDIAN = True     # True: median (+HDI) on the log scale; False: mean
HDI_PROB = 0.95
P_UNIFORM = 1.0 / (2.0 * np.pi)

N_JOBS = max(1, (os.cpu_count() or 2) - 1)


# ===========================================================================
# Shared helpers
# ===========================================================================
def wrap_pm_pi(x):
    """Wrap angles to (-pi, pi]."""
    x = np.asarray(x, float)
    return (x + np.pi) % (2 * np.pi) - np.pi


def _norm_w(w):
    """Clip and renormalize mixture weights (guard numeric drift)."""
    w = np.asarray(w, float)
    w = np.clip(w, EPS, 1.0)
    return w / np.clip(w.sum(axis=-1, keepdims=True), EPS, None)


def _combine_perp(logp_L, logp_R):
    """Combine the perpendicular-direction log densities.

    ``PERP_COMBINE = "avg"``: log(0.5 p(pi/2) + 0.5 p(-pi/2)) on the log
    scale (``"sum"`` would keep the unnormalized total mass)."""
    if PERP_COMBINE == "sum":
        return np.logaddexp(logp_L, logp_R)
    elif PERP_COMBINE == "avg":
        return np.logaddexp(logp_L, logp_R) - np.log(2.0)
    else:
        raise ValueError("PERP_COMBINE must be 'avg' or 'sum'")


def _mixture_logpdf(theta, mu, kappa, w):
    """Mixture log-density at scalar `theta`.

    ``mu``, ``kappa``, ``w``: (S, I, 2) arrays (posterior draws, animals,
    components); returns log p(theta) of shape (S, I)."""
    mu = wrap_pm_pi(mu)
    kappa = np.asarray(kappa, float)
    w = _norm_w(w)
    logc0 = np.log(w[..., 0]) + vonmises.logpdf(theta, kappa[..., 0], loc=mu[..., 0])
    logc1 = np.log(w[..., 1]) + vonmises.logpdf(theta, kappa[..., 1], loc=mu[..., 1])
    return np.logaddexp(logc0, logc1)


def _von_mises_kde(data, points, kappa):
    """Von Mises kernel density estimate of `data` at `points` (kernel
    concentration fixed at `kappa`). Returns an array of len(points)."""
    data = wrap_pm_pi(np.asarray(data, float))
    c = 1.0 / (2.0 * np.pi * i0(kappa))
    cos_diffs = np.cos(points[:, None] - data[None, :])
    return c * np.exp(kappa * cos_diffs).mean(axis=1)


# ===========================================================================
# Entry point 1: PSIS-LOO ELPD for the two fitted models
# ===========================================================================
def loo_elpd_rows(models):
    """Per-model PSIS-LOO expected log predictive density.

    ``models`` maps a display name to its fit result (bimodal first,
    unimodal second - the figure fixes that order), and returns
    ``[(name, elpd, se), ...]`` in the same order.

    ArviZ 1.3 reports the estimate as ``loo.elpd`` / ``loo.se`` (the
    pre-1.0 ``elpd_loo`` key is gone); the scale kwarg was removed with the
    log scale now the only output.
    """
    rows = []
    for name, idata in models.items():
        loo = az.loo(idata, pointwise=True)
        rows.append((name, float(loo.elpd), float(loo.se)))
    return rows


# ===========================================================================
# Entry point 2: pointwise LOO aggregated per individual
# ===========================================================================
def _pointwise_loo_contrib(idata):
    """Pointwise ELPD (and its SE) per observation for one fit result."""
    loo = az.loo(idata, pointwise=True)
    elpd_i = np.asarray(loo.elpd_i)
    # ArviZ no longer reports a per-observation SE (pre-1.0 arviz never did
    # either; the original notebook's defensive fallback was zeros), so the
    # aggregated per-individual error bars are zero-length, as in the
    # original analysis.
    se_i = np.zeros_like(elpd_i, dtype=float)
    return elpd_i.astype(float), se_i.astype(float)


def _aggregate_by_individual(values_i, se_i, individuals, n_ind):
    """Sum within-individual and combine SEs as sqrt(sum(se_i^2))."""
    val_ind = np.zeros(n_ind, dtype=float)
    se_ind = np.zeros(n_ind, dtype=float)
    for k in range(n_ind):
        mask = (individuals == k)
        val_ind[k] = values_i[mask].sum()
        se_ind[k] = np.sqrt(np.sum(se_i[mask] ** 2))
    return val_ind, se_ind


def elpd_by_individual(models, df):
    """Per-individual stratified ELPD contributions for both models, plus
    the pointwise difference (bimodal - unimodal) aggregated per animal.

    ``models`` is the same name->fit dictionary as :func:`loo_elpd_rows`
    (both entries are required); ``df`` provides the observation->animal
    mapping (``df["animal_id"]`` must be dense 0..n-1)."""
    individuals = np.asarray(df["animal_id"].to_numpy(), int)
    n_ind = individuals.max() + 1

    elpd_i_bi, se_i_bi = _pointwise_loo_contrib(models["Bimodal"])
    elpd_i_uni, se_i_uni = _pointwise_loo_contrib(models["Unimodal"])

    if (elpd_i_bi.shape[0] != individuals.shape[0]
            or elpd_i_uni.shape[0] != individuals.shape[0]):
        raise ValueError(
            f"Mismatch: loo_i lengths (bimodal={elpd_i_bi.shape[0]}, "
            f"unimodal={elpd_i_uni.shape[0]}) "
            f"but individuals={individuals.shape[0]}"
        )

    elpd_bi, se_bi = _aggregate_by_individual(
        elpd_i_bi, se_i_bi, individuals, n_ind)
    elpd_uni, se_uni = _aggregate_by_individual(
        elpd_i_uni, se_i_uni, individuals, n_ind)

    # Difference from the POINTWISE contributions (approx SE for the
    # difference assumes independence across models).
    delta_i = elpd_i_bi - elpd_i_uni
    se_delta_i = np.sqrt(se_i_bi ** 2 + se_i_uni ** 2)
    delta, se_delta = _aggregate_by_individual(
        delta_i, se_delta_i, individuals, n_ind)

    return {
        "elpd_bi": elpd_bi, "elpd_uni": elpd_uni,
        "se_bi": se_bi, "se_uni": se_uni,
        "delta": delta, "se_delta": se_delta,
    }


# ===========================================================================
# Entry points 3-6: population-level directional density estimates
# ===========================================================================
def population_density(trace_bimodal):
    """Posterior predictive density of the population mixture on a 360-point
    grid (-pi..pi): the median and 2.5/97.5 percentiles across posterior
    draws of (mu_pop, kappa_pop, weights_pop), via a compiled logp.

    Returns a dict with keys ``angles``, ``density_median``,
    ``density_lower``, ``density_upper``."""
    posterior = trace_bimodal["posterior"]
    mu_p = posterior["mu_pop"].values.reshape(-1, 2)
    kappa_p = posterior["kappa_pop"].values.reshape(-1, 2)
    weights_p = posterior["weights_pop"].values.reshape(-1, 2)

    angles_var = pt.vector("angles")
    components = pm.VonMises.dist(mu=mu_p[:, None], kappa=kappa_p[:, None])
    mixture = pm.Mixture.dist(w=weights_p[:, None], comp_dists=components)
    logp = pm.logp(mixture, angles_var[None, :])
    logp_fn = pm.compile([angles_var], logp)

    angle_values = np.linspace(-np.pi, np.pi, N_POINTS).astype(np.float64)
    logp_values = logp_fn(angle_values)

    logp_lower, logp_median, logp_upper = np.percentile(
        logp_values, [2.5, 50.0, 97.5], axis=0)

    return {
        "angles": angle_values,
        "density_lower": np.exp(logp_lower),
        "density_median": np.exp(logp_median),
        "density_upper": np.exp(logp_upper),
    }


def _hierarchical_bootstrap(data, labels, compute_func, n_bootstraps,
                            n_jobs=-1, random_state=None):
    """Two-stage hierarchical bootstrap:
      (1) Bootstrap labels (clusters) with replacement.
      (2) For each selected label, bootstrap its observations with
          replacement.
      (3) Aggregate per-label results (mean along axis=0).

    Each iteration draws from its own seed (the seeds are pre-drawn from
    ``random_state``), so the result is reproducible and independent of the
    number of parallel workers.
    """
    data = np.asarray(data, float)
    labels = np.asarray(labels)

    rng = np.random.default_rng(random_state)

    unique_labels, inverse = np.unique(labels, return_inverse=True)
    label_to_indices = [np.flatnonzero(inverse == i)
                        for i in range(len(unique_labels))]
    n_labels = len(unique_labels)

    def single_bootstrap(seed):
        r = np.random.default_rng(seed)
        sampled_label_ids = r.integers(
            0, n_labels, size=n_labels, endpoint=False)
        per_label_results = []
        for lid in sampled_label_ids:
            idx = label_to_indices[lid]
            boot_idx = idx[r.integers(0, len(idx), size=len(idx), endpoint=False)]
            per_label_results.append(compute_func(data[boot_idx]))
        per_label_results = np.array(per_label_results)
        # duplicates are naturally weighted by multiplicity
        return np.mean(per_label_results, axis=0)

    seeds = rng.integers(0, np.iinfo(np.int32).max, size=n_bootstraps, endpoint=False)

    bootstrap_results = Parallel(n_jobs=n_jobs)(
        delayed(single_bootstrap)(int(s)) for s in seeds)

    return np.asarray(bootstrap_results)


def _bootstrap_summary_stats(bootstrap_samples, axis=0, lower=2.5, upper=97.5):
    """Median and percentile bounds across bootstrap iterations."""
    if lower >= upper:
        raise ValueError("`lower` must be less than `upper`.")
    lower_bounds, median, upper_bounds = np.percentile(
        bootstrap_samples, [lower, 50.0, upper], axis=axis)
    errors = np.array([median - lower_bounds, upper_bounds - median])
    return {
        "median": median,
        "lower_bounds": lower_bounds,
        "upper_bounds": upper_bounds,
        "errors": errors,
    }


@jit(nopython=True)
def _histogram(data_sample, bins):
    """Histogram counts normalized to density (numba; the rightmost bin
    edge closes the interval)."""
    n_bins = len(bins) - 1
    hist = np.zeros(n_bins)

    for value in data_sample:
        for i in range(n_bins):
            if bins[i] <= value < bins[i + 1]:
                hist[i] += 1
                break
        if value == bins[-1]:
            hist[-1] += 1

    bin_widths = np.diff(bins)
    hist = hist / (np.sum(hist) * bin_widths)
    return hist


def bootstrap_histogram_stats(df, n_bins=N_BINS, n_bootstraps=N_BOOTSTRAPS):
    """Hierarchical-bootstrap circular histogram of the observed headings
    (animals resampled with replacement, then their steps). Returns bin
    centers and the median / 2.5/97.5 percentiles of the bootstrap
    distribution of mean histograms, as ``hist_bin_centers``,
    ``hist_median``, ``hist_lower_bounds``, ``hist_upper_bounds``,
    ``hist_errors``."""
    data = wrap_pm_pi(df["direction"].to_numpy(float))
    labels = df["animal_id"].to_numpy()

    bins = np.linspace(-np.pi, np.pi, n_bins + 1)
    bin_centers = (bins[:-1] + bins[1:]) / 2

    def compute_histogram(data_sample):
        return _histogram(data_sample, bins=bins)

    bootstrap_histograms = _hierarchical_bootstrap(
        data, labels, compute_histogram, n_bootstraps,
        random_state=BOOTSTRAP_SEED)

    summary = _bootstrap_summary_stats(bootstrap_histograms)
    return {
        "hist_bin_centers": bin_centers,
        "hist_median": summary["median"],
        "hist_lower_bounds": summary["lower_bounds"],
        "hist_upper_bounds": summary["upper_bounds"],
        "hist_errors": summary["errors"],
    }


def bootstrap_von_mises_kde_stats(df, kappa=KAPPA_KDE_POP,
                                  n_bootstraps=N_BOOTSTRAPS,
                                  n_points=N_POINTS):
    """Bootstrapped von Mises KDE of the observed headings on the population
    grid (animals resampled with replacement, then their steps). Returns the
    grid ``kde_angles`` and the median / 2.5/97.5 percentiles of the
    bootstrap KDE distribution as ``kde_median``, ``kde_lower_bounds``,
    ``kde_upper_bounds``."""
    data = wrap_pm_pi(df["direction"].to_numpy(float))
    labels = df["animal_id"].to_numpy()

    angles = np.linspace(-np.pi, np.pi, n_points)

    def compute_kde(bootstrapped_sample):
        return _von_mises_kde(bootstrapped_sample, angles, kappa)

    bootstrap_kdes = _hierarchical_bootstrap(
        data, labels, compute_kde, n_bootstraps, random_state=BOOTSTRAP_SEED)

    summary = _bootstrap_summary_stats(bootstrap_kdes)
    return {
        "kde_angles": angles,
        "kde_median": summary["median"],
        "kde_lower_bounds": summary["lower_bounds"],
        "kde_upper_bounds": summary["upper_bounds"],
    }


def population_validation_stats(trace_bimodal, df):
    """All three population-level directional density estimates (the
    bimodal model's posterior predictive density plus the bootstrap
    histogram and von Mises KDE baselines) in one dict, for the cartesian
    and polar validation figures."""
    return {
        **population_density(trace_bimodal),
        **bootstrap_histogram_stats(df),
        **bootstrap_von_mises_kde_stats(df),
    }


# ===========================================================================
# Entry point 7: per-individual posterior predictive check
# ===========================================================================
def _bootstrap_kde_summary(data, points, kappa, n_boot, seed):
    """Non-hierarchical bootstrap within one individual: median and
    2.5/97.5 percentile bands for the von Mises KDE at `points`."""
    data = wrap_pm_pi(np.asarray(data, float))
    n = data.size
    if n == 0:
        A = points.size
        return (np.full(A, np.nan), np.full(A, np.nan), np.full(A, np.nan))

    rng = np.random.default_rng(seed)
    boot = np.empty((n_boot, points.size), dtype=np.float32)
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        boot[b] = _von_mises_kde(data[idx], points, kappa).astype(np.float32)

    med = np.quantile(boot, 0.5, axis=0)
    lo = np.quantile(boot, 0.025, axis=0)
    hi = np.quantile(boot, 0.975, axis=0)
    return med.astype(float), lo.astype(float), hi.astype(float)


def _kde_worker(idx, data_1d, points, kappa, n_boot, base_seed):
    """Picklable worker for the parallel per-individual bootstrap KDEs."""
    med, lo, hi = _bootstrap_kde_summary(
        data_1d, points, kappa, n_boot=n_boot, seed=base_seed + 10_000 + idx)
    return idx, med, lo, hi


def individual_ppc_stats(trace_bimodal, df):
    """52-panel individual-level PPC arrays.

    Model side: per-individual posterior mean mixture density plus a
    pointwise 95% credible band across draws. Data side: per-individual
    non-hierarchical bootstrap von Mises KDE (median + 95% band, computed
    in parallel processes, seeds BASE_SEED + 10_000 + idx). Also returns
    the per-individual observed directions, the posterior-mean component
    parameters (for the component arrows) and the global density scale of
    the grid. Assumes a 2-component bimodal fit with posterior arrays of
    shape (chain, draw, n_individuals, 2) and dense ``df["animal_id"]``.
    """
    posterior = trace_bimodal["posterior"]
    mu_all = posterior["mu"].values            # (chain, draw, I, 2)
    kappa_all = posterior["kappa"].values      # (chain, draw, I, 2)
    weights_all = posterior["weights"].values  # (chain, draw, I, 2)

    n_chain, n_draw, n_individuals, n_components = mu_all.shape
    assert n_components == 2, "This code assumes 2 mixture components."

    n_samples = n_chain * n_draw
    mu_flat = mu_all.reshape(n_samples, n_individuals, 2)
    kappa_flat = kappa_all.reshape(n_samples, n_individuals, 2)
    weights_flat = weights_all.reshape(n_samples, n_individuals, 2)

    # Posterior means (for the mean density and the component arrows).
    mu_mean = wrap_pm_pi(mu_flat.mean(axis=0))               # (I, 2)
    kappa_mean = np.asarray(kappa_flat.mean(axis=0), float)  # (I, 2)
    weights_mean = np.asarray(weights_flat.mean(axis=0), float)
    weights_mean = weights_mean / weights_mean.sum(axis=1, keepdims=True)

    directions = wrap_pm_pi(df["direction"].to_numpy(float))
    individuals = np.asarray(df["animal_id"].to_numpy(), int)

    per_ind_data = [directions[individuals == idx]
                    for idx in range(n_individuals)]

    angles = np.linspace(-np.pi, np.pi, N_ANGLES, endpoint=False)
    A = angles.size

    # (A) Model: pointwise 95% credible band + mean PDF per individual.
    pdf_lower = np.zeros((n_individuals, A), dtype=float)
    pdf_upper = np.zeros((n_individuals, A), dtype=float)
    pdf_mean_vals = np.zeros((n_individuals, A), dtype=float)

    max_pdf_global = 0.0
    theta = angles[None, :].astype(float)  # (1, A)

    for idx in range(n_individuals):
        mu_draws = mu_flat[:, idx, :]          # (S, 2)
        kappa_draws = kappa_flat[:, idx, :]    # (S, 2)
        weight_draws = weights_flat[:, idx, :]  # (S, 2)

        mu0 = wrap_pm_pi(mu_draws[:, 0])[:, None]
        mu1 = wrap_pm_pi(mu_draws[:, 1])[:, None]
        k0 = np.asarray(kappa_draws[:, 0], float)[:, None]
        k1 = np.asarray(kappa_draws[:, 1], float)[:, None]

        w0 = np.clip(np.asarray(weight_draws[:, 0], float)[:, None], EPS, 1.0)
        w1 = np.clip(np.asarray(weight_draws[:, 1], float)[:, None], EPS, 1.0)
        ws = w0 + w1
        w0 = w0 / ws
        w1 = w1 / ws

        # Vectorized mixture PDF over (S, A).
        log_c1 = np.log(w0) + vonmises.logpdf(theta, k0, loc=mu0)
        log_c2 = np.log(w1) + vonmises.logpdf(theta, k1, loc=mu1)
        log_pdf_samples = np.logaddexp(log_c1, log_c2).astype(np.float32)
        pdf_samples = np.exp(log_pdf_samples).astype(np.float32)

        pdf_lower[idx] = np.quantile(pdf_samples, 0.025, axis=0)
        pdf_upper[idx] = np.quantile(pdf_samples, 0.975, axis=0)

        # Mean posterior PDF using the mean parameters.
        mu_ind = mu_mean[idx]
        kappa_ind = kappa_mean[idx]
        weights_ind = np.clip(weights_mean[idx], EPS, 1.0)
        weights_ind = weights_ind / weights_ind.sum()

        log_c1_mean = np.log(weights_ind[0]) + vonmises.logpdf(
            angles, kappa_ind[0], loc=mu_ind[0])
        log_c2_mean = np.log(weights_ind[1]) + vonmises.logpdf(
            angles, kappa_ind[1], loc=mu_ind[1])
        pdf_mean_vals[idx] = np.exp(np.logaddexp(log_c1_mean, log_c2_mean))

        max_pdf_global = max(max_pdf_global, float(np.nanmax(pdf_upper[idx])))

    # (B) Data: bootstrap KDE per individual (parallel).
    kde_med = np.full((n_individuals, A), np.nan, dtype=float)
    kde_lo = np.full((n_individuals, A), np.nan, dtype=float)
    kde_hi = np.full((n_individuals, A), np.nan, dtype=float)

    todo = [idx for idx in range(n_individuals) if per_ind_data[idx].size > 0]
    if len(todo) > 0:
        with ProcessPoolExecutor(max_workers=min(N_JOBS, len(todo))) as ex:
            futs = [
                ex.submit(_kde_worker, idx, per_ind_data[idx], angles,
                          KAPPA_KDE, N_BOOT_KDE, BASE_SEED)
                for idx in todo
            ]
            for f in as_completed(futs):
                idx, med, lo, hi = f.result()
                kde_med[idx] = med
                kde_lo[idx] = lo
                kde_hi[idx] = hi

        max_pdf_global = max(max_pdf_global, float(np.nanmax(kde_hi)))

    return {
        "angles": angles,
        "per_ind_data": per_ind_data,
        "pdf_lower": pdf_lower, "pdf_upper": pdf_upper,
        "pdf_mean_vals": pdf_mean_vals,
        "kde_med": kde_med, "kde_lo": kde_lo, "kde_hi": kde_hi,
        "mu_mean": mu_mean, "kappa_mean": kappa_mean,
        "weights_mean": weights_mean,
        "max_pdf_global": max_pdf_global,
        "n_individuals": n_individuals,
    }


# ===========================================================================
# Entry points 8-11: hypothesis tests at the theory-predicted directions
# ===========================================================================
def _population_mixture_logp(trace_bimodal):
    """Per-draw log predictive density of the population mixture at the
    hypothesized directions (forward 0, perpendicular +/ -pi/2)."""
    posterior = trace_bimodal["posterior"]
    components = pm.VonMises.dist(
        mu=posterior["mu_pop"].values.reshape(-1, 2),
        kappa=posterior["kappa_pop"].values.reshape(-1, 2),
    )
    directions_mixture = pm.Mixture.dist(
        w=posterior["weights_pop"].values.reshape(-1, 2),
        comp_dists=components,
    )
    logp_forward = pm.logp(directions_mixture, value=0).eval().ravel()
    logp_L = pm.logp(directions_mixture, value=-np.pi / 2).eval().ravel()
    logp_R = pm.logp(directions_mixture, value=np.pi / 2).eval().ravel()
    logp_perp = _combine_perp(logp_L, logp_R)
    return logp_forward, logp_perp


def hypothesis_density_stats(trace_bimodal):
    """Population posterior PREDICTIVE DENSITY at the hypothesized
    directions: median (or mean) and 95% HDI summarized on the log scale,
    then exponentiated back to the density scale. Returns ``centers``
    (2,), ``xerr`` (2, 2) and the uniform baseline ``P_UNIFORM``."""
    logp_forward, logp_perp = _population_mixture_logp(trace_bimodal)

    vals_log = [np.asarray(logp_forward), np.asarray(logp_perp)]
    centers_log = np.array(
        [np.median(v) if USE_MEDIAN else v.mean() for v in vals_log],
        dtype=float)
    hdis_log = np.vstack([az.hdi(v, prob=HDI_PROB) for v in vals_log]).astype(float)

    centers = np.exp(centers_log)
    hdis = np.exp(hdis_log)
    xerr = np.vstack([centers - hdis[:, 0], hdis[:, 1] - centers])

    return {"centers": centers, "xerr": xerr, "p_uniform": P_UNIFORM}


def hypothesis_elpd_stats(trace_bimodal):
    """Population ELPD (mean of the per-draw log predictive densities) at
    the hypothesized directions with +/- Z standard errors, on the log
    scale. Returns ``elpds`` (2,), ``xerr`` (2,) and the log-uniform
    reference ``logp_uniform``."""
    logp_forward, logp_perp = _population_mixture_logp(trace_bimodal)

    vals_log = [np.asarray(logp_forward), np.asarray(logp_perp)]
    elpds = np.array([v.mean() for v in vals_log], dtype=float)
    ses = np.array([v.std(ddof=1) / np.sqrt(v.size) for v in vals_log],
                   dtype=float)
    xerr = Z * ses

    return {
        "elpds": elpds,
        "xerr": xerr,
        "logp_uniform": np.log(P_UNIFORM),
    }


def _individual_hypothesis_arrays(trace_bimodal):
    """Per-draw per-individual log predictive densities at the three
    hypothesized directions; returns logp_fwd / logp_L / logp_R of shape
    (S, I) and n_ind."""
    posterior = trace_bimodal["posterior"]
    mu_all = posterior["mu"].values
    kappa_all = posterior["kappa"].values
    w_all = posterior["weights"].values

    n_chain, n_draw, n_ind, n_comp = mu_all.shape
    assert n_comp == 2, f"Expected 2 components, got {n_comp}"

    S = n_chain * n_draw
    mu = mu_all.reshape(S, n_ind, 2)
    kappa = kappa_all.reshape(S, n_ind, 2)
    w = w_all.reshape(S, n_ind, 2)

    logp_fwd = _mixture_logpdf(0.0, mu, kappa, w)
    logp_L = _mixture_logpdf(-np.pi / 2.0, mu, kappa, w)
    logp_R = _mixture_logpdf(+np.pi / 2.0, mu, kappa, w)
    return logp_fwd, logp_L, logp_R, n_ind


def _center_and_hdi_log(x_log, use_median=True, hdi_prob=HDI_PROB):
    center = np.median(x_log) if use_median else float(np.mean(x_log))
    low, high = az.hdi(x_log, prob=hdi_prob)
    return float(center), float(low), float(high)


def hypothesis_individual_stats(trace_bimodal):
    """Per-individual predictive DENSITY HDIs at the hypothesized
    directions (median + 95% HDI on the log scale, exponentiated back),
    plus P(VE > OF) per individual, with a pooled row (the expectation
    across individuals per draw) appended. Assumed 2-component bimodal
    fit; dense ``df["animal_id"]`` indexing (n_ind from the fit)."""
    logp_fwd, logp_L, logp_R, n_ind = _individual_hypothesis_arrays(trace_bimodal)
    logp_perp = _combine_perp(logp_L, logp_R)

    p_VE_gt_OF = (logp_perp > logp_fwd).mean(axis=0)  # (I,)

    center_f_log = np.empty(n_ind, float)
    center_p_log = np.empty(n_ind, float)
    hdi_f_log = np.empty((n_ind, 2), float)
    hdi_p_log = np.empty((n_ind, 2), float)

    for i in range(n_ind):
        c, lo, hi = _center_and_hdi_log(
            logp_fwd[:, i], use_median=USE_MEDIAN, hdi_prob=HDI_PROB)
        center_f_log[i] = c
        hdi_f_log[i] = (lo, hi)
        c, lo, hi = _center_and_hdi_log(
            logp_perp[:, i], use_median=USE_MEDIAN, hdi_prob=HDI_PROB)
        center_p_log[i] = c
        hdi_p_log[i] = (lo, hi)

    # Pooled row: expectation across individuals per draw (on the log
    # scale).
    logp_fwd_pool_draw = logp_fwd.mean(axis=1)
    logp_perp_pool_draw = logp_perp.mean(axis=1)

    c, lo, hi = _center_and_hdi_log(
        logp_fwd_pool_draw, use_median=USE_MEDIAN, hdi_prob=HDI_PROB)
    center_f_pool_log = c
    hdi_f_pool_log = (lo, hi)
    c, lo, hi = _center_and_hdi_log(
        logp_perp_pool_draw, use_median=USE_MEDIAN, hdi_prob=HDI_PROB)
    center_p_pool_log = c
    hdi_p_pool_log = (lo, hi)

    p_VE_gt_OF_pool = (logp_perp_pool_draw > logp_fwd_pool_draw).mean()

    # Append pooled and back-transform to the density scale.
    center_f_log_all = np.concatenate([center_f_log, [center_f_pool_log]])
    center_p_log_all = np.concatenate([center_p_log, [center_p_pool_log]])
    hdi_f_log_all = np.vstack([hdi_f_log, np.array(hdi_f_pool_log, float)])
    hdi_p_log_all = np.vstack([hdi_p_log, np.array(hdi_p_pool_log, float)])
    p_VE_gt_OF_all = np.concatenate([p_VE_gt_OF, [p_VE_gt_OF_pool]])

    return {
        "center_f_all": np.exp(center_f_log_all),
        "center_p_all": np.exp(center_p_log_all),
        "hdi_f_all": np.exp(hdi_f_log_all),
        "hdi_p_all": np.exp(hdi_p_log_all),
        "p_VE_gt_OF_all": p_VE_gt_OF_all,
        "n_ind": n_ind,
        "p_uniform": P_UNIFORM,
    }


def hypothesis_elpd_individual_stats(trace_bimodal):
    """Per-individual ELPD (mean of the per-draw log predictive
    densities) at the hypothesized directions with +/- Z standard errors,
    plus P(VE > OF) per individual, with a pooled row appended (the
    expectation across individuals, NOT the sum)."""
    logp_fwd, logp_L, logp_R, n_ind = _individual_hypothesis_arrays(
        trace_bimodal)
    logp_perp = _combine_perp(logp_L, logp_R)

    S = logp_fwd.shape[0]

    elpd_f = logp_fwd.mean(axis=0)
    elpd_p = logp_perp.mean(axis=0)
    se_f = logp_fwd.std(axis=0, ddof=1) / np.sqrt(S)
    se_p = logp_perp.std(axis=0, ddof=1) / np.sqrt(S)
    p_VE_gt_OF = (logp_perp > logp_fwd).mean(axis=0)

    # Pooled row (expectation across individuals per draw).
    logp_fwd_pool_draw = logp_fwd.mean(axis=1)
    logp_perp_pool_draw = logp_perp.mean(axis=1)

    elpd_f_pool = logp_fwd_pool_draw.mean()
    elpd_p_pool = logp_perp_pool_draw.mean()
    se_f_pool = logp_fwd_pool_draw.std(ddof=1) / np.sqrt(S)
    se_p_pool = logp_perp_pool_draw.std(ddof=1) / np.sqrt(S)
    p_VE_gt_OF_pool = (logp_perp_pool_draw > logp_fwd_pool_draw).mean()

    return {
        "elpd_f_all": np.concatenate([elpd_f, [elpd_f_pool]]),
        "elpd_p_all": np.concatenate([elpd_p, [elpd_p_pool]]),
        "se_f_all": np.concatenate([se_f, [se_f_pool]]),
        "se_p_all": np.concatenate([se_p, [se_p_pool]]),
        "p_gt_all": np.concatenate([p_VE_gt_OF, [p_VE_gt_OF_pool]]),
        "n_ind": n_ind,
        "logp_uniform": np.log(P_UNIFORM),
    }


# ===========================================================================
# Entry points 12-13: directional-alignment scores
# ===========================================================================
def population_alignment_stats(trace_bimodal):
    """Population posterior MEAN alignment (weights only; no kappa):
    forward alignment (weight-averaged cosine of the component means
    around 0) and axial perpendicular alignment (around +/-pi/2), with
    means and 95% HDIs. Returns ``means`` (2,), ``xerr`` (2, 2)."""
    posterior = trace_bimodal["posterior"]
    mu = posterior["mu_pop"].values.reshape(-1, 2)
    w = posterior["weights_pop"].values.reshape(-1, 2)

    phi_forward = 0.0
    align_forward = np.sum(w * np.cos(mu - phi_forward), axis=1)
    phi_perp = np.pi / 2
    align_perp_axial = np.sum(w * np.cos(2 * (mu - phi_perp)), axis=1)

    vals = [np.asarray(align_forward), np.asarray(align_perp_axial)]
    means = np.array([v.mean() for v in vals], dtype=float)
    hdis = np.vstack([az.hdi(v, prob=HDI_PROB) for v in vals]).astype(float)
    xerr = np.vstack([means - hdis[:, 0], hdis[:, 1] - means])

    return {"means": means, "xerr": xerr}


def individual_alignment_stats(trace_bimodal):
    """Per-individual posterior MEAN alignment scores (weights only; no
    kappa): OF = sum_k w_k cos(mu_k - 0), VE = sum_k w_k cos(2 (mu_k -
    pi/2)) (axial). Returns per-individual means and 95% HDIs plus
    P(VE > OF) per individual."""
    posterior = trace_bimodal["posterior"]
    mu_all = posterior["mu"].values
    w_all = posterior["weights"].values

    n_chain, n_draw, n_ind, n_comp = mu_all.shape
    assert n_comp == 2, f"Expected 2 components, got {n_comp}"

    S = n_chain * n_draw
    mu = wrap_pm_pi(mu_all.reshape(S, n_ind, 2))
    w = _norm_w(w_all.reshape(S, n_ind, 2))

    phi_forward = 0.0
    phi_perp = np.pi / 2

    align_forward = np.sum(w * np.cos(mu - phi_forward), axis=-1)
    # axial (+/- pi/2 treated as one direction)
    align_perp_axial = np.sum(w * np.cos(2.0 * (mu - phi_perp)), axis=-1)

    p_VE_gt_OF = (align_perp_axial > align_forward).mean(axis=0)  # (I,)

    mean_f = align_forward.mean(axis=0)
    mean_p = align_perp_axial.mean(axis=0)

    hdi_f = np.vstack([az.hdi(align_forward[:, i], prob=HDI_PROB)
                       for i in range(n_ind)])
    hdi_p = np.vstack([az.hdi(align_perp_axial[:, i], prob=HDI_PROB)
                       for i in range(n_ind)])

    return {
        "mean_f": mean_f, "mean_p": mean_p,
        "hdi_f": hdi_f, "hdi_p": hdi_p,
        "p_VE_gt_OF": p_VE_gt_OF,
        "n_ind": n_ind,
    }
