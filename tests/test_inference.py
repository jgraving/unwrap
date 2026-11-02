"""Tests for ``unwrap.inference`` (graph evaluation + pointwise LOO helpers)."""
import numpy as np
import pymc as pm
import xarray as xr
from arviz import InferenceData
from pytest import approx, raises
from scipy.stats import vonmises

from unwrap.inference import (
    aggregate_by_group,
    crosscheck_vonmises_loglike,
    elpd_diff_from_pointwise,
    graph_free_rvs,
    graph_response_params,
    loo_pointwise,
    pareto_k_report,
    se_of_sum,
)


# ---------------------------------------------------------------------------
# se_of_sum / elpd_diff_from_pointwise / pareto_k_report / aggregate_by_group
# ---------------------------------------------------------------------------
def test_se_of_sum_known_value():
    # sqrt(3 * var([1, 2, 3], ddof=1)) = sqrt(3 * 1) = sqrt(3)
    assert se_of_sum([1.0, 2.0, 3.0]) == approx(np.sqrt(3.0))


def test_se_of_sum_constant_is_zero():
    assert se_of_sum([1.0, 1.0, 1.0, 1.0]) == approx(0.0)


def test_se_of_sum_single_value_is_nan():
    assert np.isnan(se_of_sum([1.0]))
    assert np.isnan(se_of_sum([]))


def test_se_of_sum_drops_nonfinite():
    assert se_of_sum([1.0, np.nan, 3.0]) == approx(se_of_sum([1.0, 3.0]))


def test_elpd_diff_from_pointwise():
    model_i = np.array([0.5, -0.25, 1.0])
    base_i = np.array([0.25, 0.25, 0.0])
    diff_i = model_i - base_i  # [0.25, -0.5, 1.0]
    diff, se = elpd_diff_from_pointwise(model_i, base_i)
    assert diff == approx(0.75)
    assert se == approx(se_of_sum(diff_i))


def test_elpd_diff_flattens_2d_pointwise():
    diff, se = elpd_diff_from_pointwise(np.zeros((2, 3)), np.ones((2, 3)))
    assert diff == approx(-6.0)
    assert se == approx(0.0)  # a constant difference has zero variance


def test_pareto_k_report_all_below_threshold():
    rep = pareto_k_report([0.1, 0.4, 0.69])
    assert rep["n"] == 3
    assert rep["n_gt"] == 0
    assert rep["max_k"] == approx(0.69)
    assert rep["mean_excess"] == 0.0
    assert rep["max_excess"] == 0.0


def test_pareto_k_report_with_excess():
    rep = pareto_k_report([0.1, 0.8, 0.95], threshold=0.7)
    assert rep["n"] == 3
    assert rep["n_gt"] == 2
    assert rep["max_k"] == approx(0.95)
    assert rep["mean_excess"] == approx((0.1 + 0.25) / 2)
    assert rep["max_excess"] == approx(0.25)


def test_pareto_k_report_empty():
    rep = pareto_k_report([np.nan, np.inf])
    assert rep["n"] == 0
    assert rep["n_gt"] == 0
    assert np.isnan(rep["max_k"])


def test_aggregate_by_group():
    values = [1.0, 2.0, 3.0, 4.0]
    se = [0.5, 0.5, 0.5, 0.5]
    group = [0, 1, 0, 1]
    val, ses = aggregate_by_group(values, se, group, n_groups=2)
    assert val[0] == approx(4.0)
    assert val[1] == approx(6.0)
    assert ses[0] == approx(np.sqrt(0.5**2 + 0.5**2))


# ---------------------------------------------------------------------------
# loo_pointwise (synthetic InferenceData - az.loo needs posterior + log_likelihood)
# ---------------------------------------------------------------------------
def _make_idata(nchains=2, ndraws=200, nobs=50, seed=0):
    rng = np.random.default_rng(seed)
    posterior = xr.Dataset(
        {"x": (("chain", "draw", "obs_id"), rng.normal(size=(nchains, ndraws, 4)))}
    )
    log_likelihood = xr.Dataset(
        {"y": (("chain", "draw", "obs"), rng.normal(size=(nchains, ndraws, nobs)))}
    )
    sample_stats = xr.Dataset(
        {
            "diverging": (("chain", "draw"), np.zeros((nchains, ndraws), dtype=bool)),
            "energy": (("chain", "draw"), rng.normal(size=(nchains, ndraws))),
        }
    )
    return InferenceData.from_dict(
        {"posterior": posterior, "log_likelihood": log_likelihood, "sample_stats": sample_stats}
    )


def test_loo_pointwise_shapes_and_values():
    import arviz as az

    idata = _make_idata()
    loo_i, k, loo = loo_pointwise(idata)
    assert loo_i.shape == (50,)
    assert k.shape == (50,)
    assert np.isfinite(loo_i).all()
    # must agree with a direct az.loo call
    ref = az.loo(idata, pointwise=True)
    np.testing.assert_array_equal(loo_i, np.asarray(ref.elpd_i).reshape(-1))
    assert np.isfinite(loo.elpd)


def test_loo_pointwise_var_name_selects_the_log_likelihood():
    idata = _make_idata()
    loo_i_all, _, _ = loo_pointwise(idata)
    loo_i_sel, _, _ = loo_pointwise(idata, var_name="y")
    np.testing.assert_array_equal(loo_i_sel, loo_i_all)


# ---------------------------------------------------------------------------
# graph_free_rvs (stub model - only the free-RV names and the posterior are used)
# ---------------------------------------------------------------------------
class _FreeRV:
    def __init__(self, name):
        self.name = name


class _StubModel:
    def __init__(self, names):
        self._backend = type("Backend", (), {})()
        self._backend.model = type("PMModel", (), {})()
        self._backend.model.free_RVs = [_FreeRV(n) for n in names]

    @property
    def backend(self):
        return self._backend


def _make_bambi_like_idata(effect="kappa_1|Group", sigma_name="kappa_1|Group_sigma"):
    rng = np.random.default_rng(1)
    effect_vals = rng.normal(size=(2, 100, 3))
    sigma_vals = np.abs(rng.normal(size=(2, 100))) + 0.1
    posterior = xr.Dataset(
        {
            "kappa_Intercept": (("chain", "draw"), rng.normal(size=(2, 100))),
            effect: (("chain", "draw", "group"), effect_vals),
            sigma_name: (("chain", "draw"), sigma_vals),
        }
    )
    return InferenceData.from_dict({"posterior": posterior})


def test_graph_free_rvs_reconstructs_missing_offsets():
    idata = _make_bambi_like_idata()
    model = _StubModel(["kappa_Intercept", "kappa_1|Group_offset", "kappa_1|Group_sigma"])
    aug = graph_free_rvs(idata, model)
    ds = aug.posterior.to_dataset()
    assert "kappa_1|Group_offset" in ds.data_vars
    # offset = effect / sigma with sigma broadcast over the group dim
    expected = idata.posterior["kappa_1|Group"].values / idata.posterior[
        "kappa_1|Group_sigma"
    ].values[..., None]
    np.testing.assert_array_equal(ds["kappa_1|Group_offset"].values, expected)
    # the original idata is untouched
    assert "kappa_1|Group_offset" not in idata.posterior.to_dataset().data_vars


def test_graph_free_rvs_noop_when_all_stored():
    idata = _make_bambi_like_idata()
    model = _StubModel(["kappa_Intercept", "kappa_1|Group", "kappa_1|Group_sigma"])
    assert graph_free_rvs(idata, model) is idata


def test_graph_free_rvs_zero_sigma_gives_zero_offset():
    idata = _make_bambi_like_idata()
    ds = idata.posterior.to_dataset()
    ds["kappa_1|Group_sigma"] = (ds["kappa_1|Group_sigma"].dims,
                                 np.zeros_like(ds["kappa_1|Group_sigma"].values, dtype=float))
    idata2 = InferenceData.from_dict({"posterior": ds})
    model = _StubModel(["kappa_Intercept", "kappa_1|Group_offset", "kappa_1|Group_sigma"])
    aug = graph_free_rvs(idata2, model)
    got = aug.posterior.to_dataset()["kappa_1|Group_offset"].values
    # divide only where s != 0, with the out-array pre-filled with zeros
    assert (got == 0.0).all()


# ---------------------------------------------------------------------------
# graph_response_params + crosscheck_vonmises_loglike (real small PyMC graph, no sampling)
# ---------------------------------------------------------------------------
def _small_vonmises_setup():
    """Two observed angles under a (mu(x), kappa(exp(x))) von Mises density."""
    with pm.Model() as pm_model:
        x = pm.Normal("x", 0.0, 1.0, shape=2)
        pm.Deterministic("mu", x)
        pm.Deterministic("kappa", pm.math.exp(x))
        pm.VonMises("y", mu=x, kappa=pm.math.exp(x), observed=_THETA)
    # all free RVs are stored -> graph_free_rvs is a no-op on this idata
    rng = np.random.default_rng(7)
    x_vals = rng.normal(size=(4, 25, 2))
    posterior = xr.Dataset({"x": (("chain", "draw", "obs_dim"), x_vals)})
    ll = vonmises.logpdf(_THETA, kappa=np.exp(x_vals), loc=x_vals)
    log_likelihood = xr.Dataset({"y": (("chain", "draw", "obs_dim"), ll)})
    sample_stats = xr.Dataset(
        {"diverging": (("chain", "draw"), np.zeros((4, 25), dtype=bool))}
    )
    idata = InferenceData.from_dict(
        {"posterior": posterior, "log_likelihood": log_likelihood, "sample_stats": sample_stats}
    )
    # bambi-specific wrapper: .backend.model is the PyMC model
    wrapper = type("ModelWrapper", (), {})()
    wrapper.backend = type("Backend", (), {"model": pm_model})()
    return idata, wrapper


_THETA = np.array([0.5, -1.25])


def test_graph_response_params_reads_actual_links():
    idata, wrapper = _small_vonmises_setup()
    params = graph_response_params(idata, wrapper, ["mu", "kappa"])
    x = idata.posterior["x"].values
    np.testing.assert_array_equal(params["mu"].values, x)
    np.testing.assert_allclose(params["kappa"].values, np.exp(x))


def test_graph_response_params_rejects_nonfinite():
    idata, wrapper = _small_vonmises_setup()
    ds = idata.posterior.to_dataset()
    ds["x"] = (ds["x"].dims, np.full(ds["x"].shape, np.nan))
    idata2 = InferenceData.from_dict(
        {"posterior": ds, "log_likelihood": idata.log_likelihood.to_dataset(),
         "sample_stats": idata.sample_stats.to_dataset()}
    )
    with raises(AssertionError, match="non-finite"):
        graph_response_params(idata2, wrapper, ["mu", "kappa"])


def test_crosscheck_vonmises_loglike_passes_for_matching_fit():
    idata, wrapper = _small_vonmises_setup()
    params = graph_response_params(idata, wrapper, ["mu", "kappa"])
    got = crosscheck_vonmises_loglike(
        idata, params["mu"], params["kappa"], _THETA, response="y"
    )
    np.testing.assert_allclose(got, idata.log_likelihood["y"].values)


def test_crosscheck_vonmises_loglike_fails_for_wrong_stored_likelihood():
    idata, wrapper = _small_vonmises_setup()
    params = graph_response_params(idata, wrapper, ["mu", "kappa"])
    ds = idata.log_likelihood.to_dataset()
    ds["y"] = (ds["y"].dims, ds["y"].values - 1.0)  # shift the stored ll
    idata2 = InferenceData.from_dict(
        {"posterior": idata.posterior.to_dataset(), "log_likelihood": ds,
         "sample_stats": idata.sample_stats.to_dataset()}
    )
    with raises(AssertionError):
        crosscheck_vonmises_loglike(
            idata2, params["mu"], params["kappa"], _THETA, response="y"
        )
