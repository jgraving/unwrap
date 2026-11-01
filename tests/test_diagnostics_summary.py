"""Mixed-variable diagnostics and ArviZ-compatible circular summary tables."""

import arviz as az
import numpy as np
import pandas as pd
import pytest
import xarray as xr
from numpy.testing import assert_allclose, assert_array_equal

from unwrap.diagnostics import rhat, rhat_components, summary


DIAGNOSTICS = ("mcse_mean", "mcse_sd", "ess_bulk", "ess_tail", "r_hat")
PROJECTED_DIAGNOSTICS = tuple(f"{metric}_{axis}" for axis in ("sin", "cos") for metric in DIAGNOSTICS)
RAW_ANGLE_DIAGNOSTICS = ("mcse_sd", "ess_bulk", "ess_tail")


@pytest.fixture
def mixed_posterior():
    rng = np.random.default_rng(907)
    direction = rng.normal(2.8, 0.45, (4, 400, 2))
    direction += 2 * np.pi * np.arange(4)[:, None, None]
    return xr.Dataset(
        {
            "direction": (("chain", "draw", "animal"), direction),
            "speed": (("chain", "draw", "animal"), rng.lognormal(0, 0.3, (4, 400, 2))),
            "kappa": (("chain", "draw"), rng.gamma(3, 1, (4, 400))),
        },
        coords={"chain": [2, 4, 6, 8], "draw": np.arange(400), "animal": ["ant", "bee"]},
    )


def assert_summary_equal(actual, expected):
    if isinstance(expected, xr.Dataset):
        xr.testing.assert_identical(actual, expected)
    else:
        pd.testing.assert_frame_equal(actual, expected)


@pytest.mark.parametrize("method", ["rank", "split", "folded", "z_scale", "identity"])
@pytest.mark.parametrize("circ_method", ["circular", "sin", "cos"])
def test_mixed_rhat_uses_requested_methods_per_variable(mixed_posterior, method, circ_method):
    actual = rhat(
        mixed_posterior, circ_var_names="direction", method=method, circ_method=circ_method
    )
    ordinary = az.rhat(mixed_posterior, method=method)
    xr.testing.assert_identical(actual[["speed", "kappa"]], ordinary[["speed", "kappa"]])
    if circ_method == "circular":
        expected = rhat(mixed_posterior, var_names="direction", method="circular")
    else:
        expected = rhat_components(mixed_posterior, var_names="direction")[circ_method]
    xr.testing.assert_identical(actual.direction, expected.direction)
    assert list(actual.data_vars) == list(mixed_posterior.data_vars)


@pytest.mark.parametrize("names", [None, []])
def test_no_circular_names_preserve_rhat_exactly(mixed_posterior, names):
    xr.testing.assert_identical(
        rhat(mixed_posterior, circ_var_names=names), az.rhat(mixed_posterior)
    )


def test_empty_circular_names_preserve_numpy_and_explicit_circular_methods(mixed_posterior):
    draws = mixed_posterior.kappa.values
    assert_array_equal(rhat(draws, circ_var_names=[]), az.rhat(draws))
    xr.testing.assert_identical(
        rhat(mixed_posterior, method="circular", circ_var_names=[]),
        rhat(mixed_posterior, method="circular"),
    )


@pytest.mark.parametrize("container", ["dataset", "dict", "dataarray", "inferencedata"])
def test_mixed_rhat_accepts_arviz_labeled_containers(mixed_posterior, container):
    data = mixed_posterior
    if container == "dict":
        data = {name: value.values for name, value in data.data_vars.items()}
    elif container == "dataarray":
        data = data.direction
    elif container == "inferencedata":
        data = az.InferenceData(posterior=data, prior=data * 0)
    actual = rhat(data, circ_var_names=["direction"])
    expected = rhat(data, var_names="direction", method="circular")
    xr.testing.assert_allclose(actual.direction, expected.direction)


def test_mixed_rhat_respects_selection_and_named_dimensions(mixed_posterior):
    transposed = mixed_posterior.transpose("animal", "draw", "chain", missing_dims="ignore")
    included = rhat(transposed, var_names=["direction", "kappa"], circ_var_names=["direction"])
    excluded = rhat(transposed, var_names=["~speed"], circ_var_names=["direction"])
    xr.testing.assert_identical(included, excluded)
    assert_array_equal(included.animal, ["ant", "bee"])
    xr.testing.assert_identical(
        rhat(mixed_posterior, var_names=["~direction"], circ_var_names=["direction"]),
        az.rhat(mixed_posterior, var_names=["~direction"]),
    )


def test_mixed_rhat_rejects_ambiguous_or_unknown_circular_names(mixed_posterior):
    with pytest.raises(KeyError, match="missing"):
        rhat(mixed_posterior, circ_var_names=["missing"], var_names=["speed"])
    with pytest.raises(KeyError, match="dir"):
        rhat(mixed_posterior, circ_var_names=["dir"])
    with pytest.raises(ValueError):
        rhat(mixed_posterior, circ_var_names=["direction"], method="circular")
    with pytest.raises(ValueError):
        rhat(mixed_posterior.kappa.values, circ_var_names=["kappa"])
    with pytest.raises((TypeError, ValueError)):
        rhat(mixed_posterior, circ_var_names=["direction"], circ_method="missing")


def test_mixed_rhat_ignores_methods_for_unselected_variable_types(mixed_posterior):
    xr.testing.assert_identical(
        rhat(
            mixed_posterior,
            circ_var_names=["direction"],
            var_names=["speed"],
            circ_method="unused",
        ),
        az.rhat(mixed_posterior, var_names=["speed"]),
    )
    xr.testing.assert_identical(
        rhat(
            mixed_posterior,
            circ_var_names=["direction"],
            var_names=["direction"],
            method="unused",
        ),
        rhat(mixed_posterior, var_names=["direction"], method="circular"),
    )


def test_mixed_wrappers_preserve_input_and_rhat_metadata(mixed_posterior):
    mixed_posterior.attrs["model"] = "direction and speed"
    mixed_posterior.direction.attrs["units"] = "radians"
    mixed_posterior.animal.attrs["description"] = "individual identifier"
    original = mixed_posterior.copy(deep=True)
    actual = rhat(mixed_posterior, circ_var_names=["direction"])
    reference = az.rhat(mixed_posterior)
    assert actual.attrs == reference.attrs
    assert actual.direction.attrs == reference.direction.attrs
    xr.testing.assert_identical(actual.animal, reference.animal)
    summary(mixed_posterior, circ_var_names=["direction"])
    xr.testing.assert_identical(mixed_posterior, original)


def test_mixed_rhat_retains_lazy_dask_results(mixed_posterior):
    pytest.importorskip("dask.array")
    lazy = mixed_posterior.chunk({"chain": -1, "draw": -1, "animal": 1})
    actual = rhat(
        lazy,
        circ_var_names="direction",
        dask_kwargs={"dask": "parallelized", "output_dtypes": [float]},
    )
    assert all(hasattr(value.data, "dask") for value in actual.data_vars.values())
    xr.testing.assert_allclose(actual.compute(), rhat(mixed_posterior, circ_var_names="direction"))


@pytest.mark.parametrize("fmt", ["wide", "long", "xarray"])
@pytest.mark.parametrize("kind", ["all", "stats", "diagnostics"])
@pytest.mark.parametrize("stat_focus", ["mean", "median"])
def test_summary_without_circular_names_is_exact_arviz(mixed_posterior, fmt, kind, stat_focus):
    kwargs = dict(fmt=fmt, kind=kind, stat_focus=stat_focus, coords={"animal": ["bee"]})
    assert_summary_equal(summary(mixed_posterior, **kwargs), az.summary(mixed_posterior, **kwargs))


def test_summary_empty_circular_selection_is_exact_arviz(mixed_posterior):
    assert_summary_equal(
        summary(mixed_posterior, circ_var_names=[], circ_method="sin"),
        az.summary(mixed_posterior, circ_var_names=[]),
    )


@pytest.mark.parametrize("circ_method", ["circular", "sin", "cos"])
def test_summary_reports_selected_rhat_and_preserves_arviz_statistics(mixed_posterior, circ_method):
    actual = summary(
        mixed_posterior, circ_var_names=["direction"], circ_method=circ_method, fmt="xarray"
    )
    base = az.summary(mixed_posterior, circ_var_names=["direction"], fmt="xarray")
    for variable in ("speed", "kappa"):
        xr.testing.assert_identical(actual[variable].sel(metric=base.metric), base[variable])
        assert np.isnan(actual[variable].sel(metric=list(PROJECTED_DIAGNOSTICS))).all()
    assert_allclose(
        actual.direction.sel(metric=["mean", "sd", "hdi_3%", "hdi_97%", "mcse_mean"]),
        base.direction.sel(metric=["mean", "sd", "hdi_3%", "hdi_97%", "mcse_mean"]),
    )
    expected_rhat = rhat(mixed_posterior, circ_var_names="direction", circ_method=circ_method)
    assert_allclose(actual.direction.sel(metric="r_hat"), expected_rhat.direction)
    assert np.isnan(actual.direction.sel(metric=list(RAW_ANGLE_DIAGNOSTICS))).all()
    assert set(actual.metric.values) == set(base.metric.values) | set(PROJECTED_DIAGNOSTICS)


def test_summary_projected_diagnostics_are_arviz_values(mixed_posterior):
    actual = summary(mixed_posterior, circ_var_names="direction", fmt="xarray")
    for axis, transform in (("sin", np.sin), ("cos", np.cos)):
        expected = az.summary(transform(mixed_posterior[["direction"]]), kind="diagnostics", fmt="xarray")
        for metric in DIAGNOSTICS:
            assert_allclose(
                actual.direction.sel(metric=f"{metric}_{axis}"),
                expected.direction.sel(metric=metric),
            )


def test_summary_diagnostics_kind_has_same_circular_mcse_and_rhat(mixed_posterior):
    complete = summary(mixed_posterior, circ_var_names="direction", fmt="xarray")
    diagnostics = summary(
        mixed_posterior, circ_var_names="direction", kind="diagnostics", fmt="xarray"
    )
    xr.testing.assert_identical(complete.sel(metric=diagnostics.metric), diagnostics)


def test_summary_optional_projected_columns_do_not_restore_raw_angle_diagnostics(mixed_posterior):
    actual = summary(
        mixed_posterior,
        circ_var_names="direction",
        include_circular_diagnostics=False,
        fmt="xarray",
    )
    assert not set(PROJECTED_DIAGNOSTICS) & set(actual.metric.values)
    assert np.isnan(actual.direction.sel(metric=list(RAW_ANGLE_DIAGNOSTICS))).all()
    assert_allclose(
        actual.direction.sel(metric="r_hat"),
        rhat(mixed_posterior, var_names="direction", method="circular").direction,
    )


@pytest.mark.parametrize("fmt", ["wide", "long", "xarray"])
def test_summary_stats_kind_matches_arviz_without_diagnostic_columns(mixed_posterior, fmt):
    kwargs = dict(circ_var_names=["direction"], kind="stats", fmt=fmt, round_to="none")
    assert_summary_equal(summary(mixed_posterior, **kwargs), az.summary(mixed_posterior, **kwargs))


def test_summary_extend_false_keeps_only_custom_statistics(mixed_posterior):
    kwargs = dict(
        circ_var_names=["direction"],
        stat_funcs={"minimum": np.min, "maximum": np.max},
        extend=False,
        round_to="none",
    )
    assert_summary_equal(summary(mixed_posterior, **kwargs), az.summary(mixed_posterior, **kwargs))


def test_summary_custom_statistic_is_retained_with_default_statistics(mixed_posterior):
    actual = summary(
        mixed_posterior, circ_var_names="direction", stat_funcs={"minimum": np.min}, fmt="xarray"
    )
    assert_allclose(actual.direction.sel(metric="minimum"), mixed_posterior.direction.min(("chain", "draw")))
    assert_allclose(actual.speed.sel(metric="minimum"), mixed_posterior.speed.min(("chain", "draw")))


@pytest.mark.parametrize("name", ["r_hat", "r_hat_sin"])
def test_summary_rejects_custom_statistics_that_would_overwrite_diagnostics(mixed_posterior, name):
    with pytest.raises(ValueError, match="Custom statistic names"):
        summary(
            mixed_posterior,
            circ_var_names=["direction"],
            stat_funcs={name: lambda draws: 1.234567},
        )


def test_summary_custom_projection_names_keep_arviz_rounding_when_projections_disabled(
    mixed_posterior,
):
    stat_funcs = {
        "r_hat_sin": lambda draws: 1.234567,
        "ess_bulk_cos": lambda draws: 12.345678,
    }
    actual = summary(
        mixed_posterior,
        circ_var_names=["direction"],
        stat_funcs=stat_funcs,
        include_circular_diagnostics=False,
    )
    expected = az.summary(mixed_posterior, circ_var_names=["direction"], stat_funcs=stat_funcs)
    pd.testing.assert_frame_equal(actual[list(stat_funcs)], expected[list(stat_funcs)])


def test_summary_named_coords_filtering_and_transposed_dimensions(mixed_posterior):
    transposed = mixed_posterior.transpose("animal", "draw", "chain", missing_dims="ignore")
    actual = summary(
        transposed,
        var_names="direct",
        filter_vars="like",
        circ_var_names=["direction"],
        coords={"animal": ["bee"]},
        fmt="xarray",
    )
    expected = summary(
        mixed_posterior.sel(animal=["bee"]),
        var_names=["direction"],
        circ_var_names=["direction"],
        fmt="xarray",
    )
    xr.testing.assert_identical(actual, expected)
    assert_array_equal(actual.animal, ["bee"])
    assert list(actual.data_vars) == ["direction"]


@pytest.mark.parametrize("fmt", ["wide", "long", "xarray"])
def test_summary_scalar_coordinate_selection_matches_preselected_data(mixed_posterior, fmt):
    kwargs = dict(circ_var_names=["direction"], fmt=fmt)
    actual = summary(mixed_posterior, coords={"animal": "bee"}, **kwargs)
    expected = summary(mixed_posterior.sel(animal="bee"), **kwargs)
    assert_summary_equal(actual, expected)


def test_summary_group_and_prior_fallback_follow_arviz(mixed_posterior):
    prior = mixed_posterior + 0.3
    data = az.InferenceData(posterior=mixed_posterior, prior=prior)
    expected = summary(prior, circ_var_names=["direction"], fmt="xarray")
    xr.testing.assert_identical(
        summary(data, group="prior", circ_var_names=["direction"], fmt="xarray"), expected
    )
    xr.testing.assert_identical(
        summary(az.InferenceData(prior=prior), circ_var_names=["direction"], fmt="xarray"), expected
    )


def test_summary_selected_out_circular_variables_are_ignored(mixed_posterior):
    for kwargs in ({"var_names": ["~direction"]}, {"var_names": "^speed$", "filter_vars": "regex"}):
        actual = summary(mixed_posterior, circ_var_names=["direction"], **kwargs)
        expected = az.summary(mixed_posterior, **kwargs)
        assert_summary_equal(actual, expected)
    actual = summary(
        mixed_posterior, circ_var_names=["direction"], var_names=["speed"], stat_focus="median"
    )
    assert_summary_equal(actual, az.summary(mixed_posterior, var_names=["speed"], stat_focus="median"))


@pytest.mark.parametrize("kind", ["all", "stats", "diagnostics"])
def test_summary_unknown_circular_name_is_an_error_even_if_excluded(mixed_posterior, kind):
    with pytest.raises(KeyError, match="missing"):
        summary(mixed_posterior, circ_var_names=["missing"], var_names=["speed"], kind=kind)


def test_summary_median_rejects_active_circular_variables(mixed_posterior):
    with pytest.raises(TypeError, match="Circular stats not supported"):
        summary(mixed_posterior, circ_var_names=["direction"], stat_focus="median")


def test_summary_formats_represent_the_same_values(mixed_posterior):
    kwargs = dict(circ_var_names=["direction"], round_to="none")
    labeled = summary(mixed_posterior, fmt="xarray", **kwargs)
    wide = summary(mixed_posterior, fmt="wide", **kwargs)
    long = summary(mixed_posterior, fmt="long", **kwargs)
    assert list(wide.columns) == list(labeled.metric.values)
    assert list(wide.index) == ["direction[ant]", "direction[bee]", "speed[ant]", "speed[bee]", "kappa"]
    expected_long = labeled.to_dataframe().reset_index().set_index("metric")
    expected_long.index = list(expected_long.index)
    pd.testing.assert_frame_equal(long, expected_long)
    for variable in ("direction", "speed"):
        for animal in ("ant", "bee"):
            assert_allclose(wide.loc[f"{variable}[{animal}]"], labeled[variable].sel(animal=animal))
    assert_allclose(wide.loc["kappa"], labeled.kappa)


def test_summary_custom_labeller_does_not_break_circular_row_assignment(mixed_posterior):
    labeller = az.labels.MapLabeller(var_name_map={"direction": "heading", "speed": "velocity"})
    actual = summary(
        mixed_posterior, circ_var_names="direction", labeller=labeller, round_to="none"
    )
    assert "heading[ant]" in actual.index
    assert "velocity[ant]" in actual.index
    assert np.isnan(actual.loc["heading[ant]", "ess_bulk"])
    assert np.isfinite(actual.loc["velocity[ant]", "ess_bulk"])
    assert_allclose(
        actual.loc["heading[ant]", "r_hat"],
        rhat(mixed_posterior, circ_var_names="direction").direction.sel(animal="ant"),
    )


@pytest.mark.parametrize("fmt", ["wide", "long", "xarray"])
def test_summary_explicit_rounding_happens_after_diagnostics_are_added(mixed_posterior, fmt):
    unrounded = summary(mixed_posterior, circ_var_names="direction", fmt=fmt, round_to="none")
    rounded = summary(mixed_posterior, circ_var_names="direction", fmt=fmt, round_to=4)
    assert_summary_equal(rounded, unrounded.round(4))


def test_summary_default_wide_rounding_uses_metric_families(mixed_posterior):
    unrounded = summary(mixed_posterior, circ_var_names="direction", round_to="none")
    expected = unrounded.round(
        {name: 2 if name.startswith("r_hat") else 0 if name.startswith("ess_") else 3 for name in unrounded.columns}
    )
    assert_summary_equal(summary(mixed_posterior, circ_var_names="direction"), expected)


@pytest.mark.parametrize("invalid", [np.nan, np.inf, -np.inf])
def test_summary_nonfinite_circular_values_propagate_nan(mixed_posterior, invalid):
    mixed_posterior.direction.values[0, 0, 0] = invalid
    result = summary(mixed_posterior, circ_var_names="direction", fmt="xarray")
    assert np.isnan(result.direction.sel(animal="ant", metric="r_hat"))
    assert np.isnan(result.direction.sel(animal="ant", metric=["r_hat_sin", "r_hat_cos"])).all()
    assert np.isfinite(result.speed.sel(metric="r_hat")).all()
