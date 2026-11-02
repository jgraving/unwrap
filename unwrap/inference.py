"""Inference helpers shared by the fitted example notebooks.

Two families that the study scripts each reimplemented:

* **Graph evaluation for Bambi fits.** :func:`graph_free_rvs`,
  :func:`graph_response_params` and :func:`crosscheck_vonmises_loglike`
  work around the upstream Bambi bug in which
  ``Model.predict(kind="response_params")`` applies the vonmises family's
  angular *mean* link to every response component, so the ``kappa`` linear
  predictor comes back wrapped into an angle instead of
  ``softplus(eta)`` (bambinos/bambi#1084). The parameters are read from
  the PyMC graph instead, which uses each component's actual link.

* **Pointwise PSIS-LOO bookkeeping.** :func:`loo_pointwise`,
  :func:`se_of_sum`, :func:`elpd_diff_from_pointwise`,
  :func:`pareto_k_report` and :func:`aggregate_by_group` wrap ArviZ 1.x
  pointwise LOO objects (``elpd_i`` / ``pareto_k``; the pre-1.0 ``loo_i``
  key is gone) and provide the standard-error, delta-ELPD, Pareto-k and
  per-group aggregation arithmetic used by every model comparison.
"""
import numpy as np
import arviz as az
import pymc as pm
from scipy.stats import vonmises

__all__ = [
    "graph_free_rvs",
    "graph_response_params",
    "crosscheck_vonmises_loglike",
    "loo_pointwise",
    "se_of_sum",
    "elpd_diff_from_pointwise",
    "pareto_k_report",
    "aggregate_by_group",
]


# ---------------------------------------------------------------------------
# Graph evaluation for Bambi fits (bambinos/bambi#1084 workaround)
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
        effect = ds[name[:-len("_offset")]]
        # The graph computes effect = sigma * offset (sigma scalar per
        # draw, offset one factor per group), so the offset carries the
        # group dimension of the centered effect.
        e = np.asarray(effect)
        s = np.asarray(sigma).reshape(np.asarray(sigma).shape + (1,) * (e.ndim - np.asarray(sigma).ndim))
        offset = np.divide(e, s, out=np.zeros(e.shape), where=s != 0)
        ds[name] = (tuple(effect.dims), offset)
    aug["posterior"] = ds
    return aug


def graph_response_params(idata, model, var_names, *, progressbar=False):
    """Response parameters (e.g. ``["mu", "kappa"]``) read from the PyMC
    graph of a Bambi fit.

    Evaluates the named deterministics over the posterior of
    :func:`graph_free_rvs`, so every component is passed through its
    *actual* link — the safe alternative to
    ``model.predict(kind="response_params")``, which applies the vonmises
    mean link to all components (bambinos/bambi#1084). Returns the
    deterministics as ``{name: DataArray}`` with ``(chain, draw, ...obs)``
    dimensions and asserts that every value is finite.
    """
    params = pm.compute_deterministics(
        graph_free_rvs(idata, model),
        var_names=list(var_names),
        model=model.backend.model,
        progressbar=progressbar,
    )
    for name, da in params.items():
        if not np.isfinite(da.values).all():
            raise AssertionError(f"non-finite values in graph response parameter {name!r}")
    return params


def _to_obs_values(da, obs_ll):
    """Values of a ``(chain, draw, *obs)`` parameter aligned with the
    rows of the stored log-likelihood.

    Bambi names the observation dimension identically in the posterior
    and the log-likelihood group, so a transpose lines things up when the
    order differs elsewhere; if the names differ at all (e.g. a pure PyMC
    model), the positional order is the same ``(chain, draw, *obs)`` and
    the values need no transpose.
    """
    if tuple(da.dims) == tuple(obs_ll.dims):
        return da.transpose(*obs_ll.dims).values
    if set(da.dims) == set(obs_ll.dims):
        return da.transpose(*obs_ll.dims).values
    return da.values


def crosscheck_vonmises_loglike(idata, mu, kappa, theta, *, response,
                                rtol=1e-10, atol=1e-10):
    """Assert graph response parameters against bambi's stored log-likelihood.

    ``mu`` and ``kappa`` are ``(chain, draw, n_obs)`` parameters (e.g. from
    :func:`graph_response_params`); ``theta`` are the observed angles in
    the same row order; ``response`` is the name of the observed variable
    inside ``idata.log_likelihood``. The von Mises log-density of the
    observations under the graph parameters must equal the log-likelihood
    stored during fitting to floating-point precision, so the cross-check
    proves both that the graph evaluation used the right links and that
    the posterior and the stored log-likelihood agree.

    Returns the expected log-likelihood array (post-assertion).
    """
    obs_ll = idata.log_likelihood[response]
    expected = vonmises.logpdf(
        np.asarray(theta, float),
        kappa=_to_obs_values(kappa, obs_ll),
        loc=_to_obs_values(mu, obs_ll),
    )
    np.testing.assert_allclose(expected, obs_ll.values, rtol=rtol, atol=atol)
    return expected


# ---------------------------------------------------------------------------
# Pointwise PSIS-LOO bookkeeping (ArviZ 1.x)
# ---------------------------------------------------------------------------
def loo_pointwise(idata, *, var_name=None):
    """PSIS-LOO with the pointwise result unpacked.

    Returns ``(elpd_i, pareto_k, loo)``: the per-observation ELPD and
    Pareto-k as flat arrays, and the full ArviZ LOO object (total ``elpd``,
    ``se``, ...) for callers that want more. ArviZ 1.x reports the
    pointwise ELPD as ``elpd_i``; the pre-1.0 name was ``loo_i``.
    """
    loo = az.loo(idata, pointwise=True, var_name=var_name)
    loo_i = np.asarray(loo.elpd_i).reshape(-1)
    k = np.asarray(loo.pareto_k).reshape(-1)
    return loo_i, k, loo


def se_of_sum(values):
    """Rough SE for a SUM of pointwise values: ``sqrt(n * var(v))``.

    Non-finite entries are dropped; ``nan`` when fewer than two remain.
    """
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    n = v.size
    if n <= 1:
        return np.nan
    return float(np.sqrt(n * np.var(v, ddof=1)))


def elpd_diff_from_pointwise(model_elpd_i, base_elpd_i):
    """Delta-ELPD between two pointwise LOO results and its SE.

    ``model_elpd_i - base_elpd_i`` summed over observations (positive: the
    model is preferred), with the :func:`se_of_sum` standard error.
    """
    diff_i = (np.asarray(model_elpd_i) - np.asarray(base_elpd_i)).reshape(-1)
    return float(diff_i.sum()), se_of_sum(diff_i)


def pareto_k_report(pareto_k, threshold=0.7):
    """Pareto-k diagnostics for a pointwise LOO result.

    Returns ``{"n", "n_gt", "max_k", "mean_excess", "max_excess"}`` where
    the excess quantities are measured above ``threshold`` (0.0 when no
    k exceeds it). Non-finite k values are dropped.
    """
    k = np.asarray(pareto_k, float)
    k = k[np.isfinite(k)]
    if k.size == 0:
        return {"n": 0, "n_gt": 0, "max_k": np.nan, "mean_excess": np.nan, "max_excess": np.nan}

    excess = k - threshold
    gt = excess > 0
    n_gt = int(np.sum(gt))

    if n_gt == 0:
        return {
            "n": int(k.size),
            "n_gt": 0,
            "max_k": float(np.max(k)),
            "mean_excess": 0.0,
            "max_excess": 0.0,
        }

    return {
        "n": int(k.size),
        "n_gt": n_gt,
        "max_k": float(np.max(k)),
        "mean_excess": float(np.mean(excess[gt])),
        "max_excess": float(np.max(excess[gt])),
    }


def aggregate_by_group(values_i, se_i, group, n_groups):
    """Sum pointwise values within groups and combine their SEs.

    ``group`` maps each observation to a dense code in ``0..n_groups-1``;
    returns ``(group_vals, group_ses)`` with the SEs combined as
    ``sqrt(sum(se_i^2))``.
    """
    values_i = np.asarray(values_i, float)
    se_i = np.asarray(se_i, float)
    group = np.asarray(group, int)
    val = np.zeros(n_groups, dtype=float)
    se = np.zeros(n_groups, dtype=float)
    for g in range(n_groups):
        mask = group == g
        val[g] = values_i[mask].sum()
        se[g] = float(np.sqrt(np.sum(se_i[mask] ** 2)))
    return val, se
