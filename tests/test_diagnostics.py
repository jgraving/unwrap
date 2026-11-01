"""ArviZ compatibility and mathematical behavior of circular convergence checks."""

import logging

import arviz as az
import numpy as np
import pytest
import xarray as xr
from numpy.testing import assert_allclose, assert_array_equal
from scipy.linalg import eigvalsh

from unwrap.diagnostics import rhat, rhat_components


@pytest.fixture
def angular_draws():
    return np.random.default_rng(103).normal(0.4, 0.45, (4, 1200))


@pytest.fixture
def posterior():
    rng = np.random.default_rng(381)
    return xr.Dataset(
        {
            "direction": (("chain", "draw", "animal"), rng.normal(0.4, 0.4, (4, 600, 3))),
            "offset": (("chain", "draw"), rng.normal(0, 0.3, (4, 600))),
        },
        coords={"chain": [2, 4, 6, 8], "draw": np.arange(600), "animal": ["ant", "bee", "moth"]},
    )


@pytest.mark.parametrize("method", ["rank", "split", "folded", "z_scale", "identity"])
def test_existing_arviz_methods_are_unchanged(angular_draws, posterior, method):
    assert_array_equal(rhat(angular_draws, method=method), az.rhat(angular_draws, method=method))
    xr.testing.assert_identical(
        rhat(posterior, var_names=["direction"], method=method),
        az.rhat(posterior, var_names=["direction"], method=method),
    )


def test_default_is_arviz_rank(angular_draws):
    assert_array_equal(rhat(angular_draws), az.rhat(angular_draws))


def test_invalid_method_raises_arviz_error(angular_draws):
    with pytest.raises(TypeError, match="R-hat method.*not found"):
        rhat(angular_draws, method="not-a-method")


@pytest.mark.parametrize("container", ["dataset", "dataarray", "inferencedata", "dict"])
def test_circular_labeled_inputs_match_individual_arrays(posterior, container):
    if container == "dataarray":
        data = posterior.direction
    elif container == "inferencedata":
        # An unrelated group must not affect the posterior diagnostic.
        data = az.InferenceData(posterior=posterior, prior=posterior * 0)
    elif container == "dict":
        data = {"direction": posterior.direction.values, "offset": posterior.offset.values}
    else:
        data = posterior

    result = rhat(data, method="circular")
    assert isinstance(result, xr.Dataset)
    assert result.direction.ndim == 1
    assert_allclose(
        result.direction.values,
        [rhat(posterior.direction.values[:, :, i], method="circular") for i in range(3)],
    )
    if container != "dict":
        assert_array_equal(result.animal, posterior.animal)
    if container != "dataarray":
        assert_allclose(result.offset, rhat(posterior.offset.values, method="circular"))


def test_unnamed_dataarray_uses_arviz_default_name(posterior):
    result = rhat(posterior.offset.rename(None), method="circular")
    assert list(result.data_vars) == ["x"]
    assert_allclose(result.x, rhat(posterior.offset.values, method="circular"))


def test_variable_selection_and_dimension_order(posterior):
    expected = rhat(posterior, method="circular")["direction"]
    transposed = posterior.transpose("animal", "draw", "chain", missing_dims="ignore")
    included = rhat(transposed, var_names="direction", method="circular")
    excluded = rhat(transposed, var_names=["~offset"], method="circular")
    xr.testing.assert_identical(included, excluded)
    xr.testing.assert_identical(included.direction, expected)
    with pytest.raises(KeyError, match="missing"):
        rhat(posterior, var_names=["missing"], method="circular")


def test_selection_happens_before_diagnosing_other_variables(posterior, caplog):
    data = posterior.assign(broken=posterior.offset * np.nan)
    with caplog.at_level(logging.WARNING):
        result = rhat(data, var_names="direction", method="circular")
    assert np.isfinite(result.direction).all()
    assert not caplog.records


def test_circular_dask_matches_eager(posterior):
    pytest.importorskip("dask.array")
    lazy = posterior.chunk({"chain": -1, "draw": -1, "animal": 1})
    result = rhat(
        lazy,
        method="circular",
        dask_kwargs={"dask": "parallelized", "output_dtypes": [float]},
    )
    assert hasattr(result.direction.data, "dask")
    xr.testing.assert_allclose(result.compute(), rhat(posterior, method="circular"))


def test_standard_arviz_dask_behavior_is_delegated(posterior):
    pytest.importorskip("dask.array")
    lazy = posterior.chunk({"chain": -1, "draw": -1, "animal": 1})
    kwargs = {"dask": "parallelized", "output_dtypes": [float]}
    actual = rhat(lazy, dask_kwargs=kwargs)
    expected = az.rhat(lazy, dask_kwargs=kwargs)
    xr.testing.assert_allclose(actual.compute(), expected.compute())


@pytest.mark.parametrize("shape", [(), (8,), (1, 100), (2, 3), (2, 0), (0, 10)])
def test_insufficient_chains_or_draws_return_nan_and_log(shape, caplog):
    with caplog.at_level(logging.WARNING):
        result = rhat(np.zeros(shape), method="circular")
    assert np.isnan(result)
    assert "shape" in caplog.text.lower()


def test_numpy_parameter_dimensions_require_labeled_input():
    with pytest.raises(TypeError, match="uni-dimensional ndarray"):
        rhat(np.zeros((4, 100, 3)), method="circular")


@pytest.mark.parametrize("invalid", [np.nan, np.inf, -np.inf])
def test_nonfinite_angles_return_nan_and_log(angular_draws, invalid, caplog):
    angular_draws[0, 0] = invalid
    with caplog.at_level(logging.WARNING):
        result = rhat(angular_draws, method="circular")
    assert np.isnan(result)
    assert caplog.records


@pytest.mark.parametrize("reference", [0.0, 0.4, np.pi])
def test_constant_equivalent_directions_return_nan(reference):
    turns = np.array([-2, 0, 1, 3])[:, None]
    values = np.broadcast_to(reference + turns * 2 * np.pi, (4, 100)).copy()
    assert np.isnan(rhat(values, method="circular"))


def test_different_constant_directions_return_infinity():
    values = np.broadcast_to(np.array([0.0, 0.0, 1.0, 1.0])[:, None], (4, 100)).copy()
    assert np.isposinf(rhat(values, method="circular"))


def test_constant_directions_with_draw_specific_period_offsets():
    rng = np.random.default_rng(91)
    turns = 2 * np.pi * rng.integers(-5, 6, size=(4, 100))
    assert np.isnan(rhat(0.4 + turns, method="circular"))
    different = np.array([0.4, 0.4, 1.2, 1.2])[:, None] + turns
    assert np.isposinf(rhat(different, method="circular"))


def test_large_common_origin_does_not_erase_distinct_fixed_directions():
    values = np.broadcast_to((1e16 + np.array([0.0, 2.0, 4.0, 6.0]))[:, None], (4, 100))
    assert np.isposinf(rhat(values, method="circular"))


def test_large_distant_constant_is_not_assumed_period_equivalent():
    values = np.broadcast_to(np.array([0.0, 0.0, 1e16, 1e16])[:, None], (4, 100))
    assert np.isposinf(rhat(values, method="circular"))


def test_rank_deficient_within_covariance_returns_nan_and_logs(caplog):
    # Two exact directions have a one-dimensional chord: W cannot be inverted.
    values = np.tile([-np.pi / 4, np.pi / 4], (4, 100))
    with caplog.at_level(logging.WARNING):
        result = rhat(values, method="circular")
    assert np.isnan(result)
    assert "singular" in caplog.text.lower() or "rank" in caplog.text.lower()


def test_period_and_rotation_invariance(angular_draws):
    baseline = rhat(angular_draws, method="circular")
    rng = np.random.default_rng(381)
    turns = rng.integers(-5, 6, size=angular_draws.shape)
    for rotation in (0.0, 0.9, np.pi, -2.6):
        assert_allclose(
            rhat(angular_draws + rotation + 2 * np.pi * turns, method="circular"),
            baseline,
            rtol=1e-10,
            atol=1e-12,
        )


def test_whole_turn_offsets_do_not_create_between_chain_disagreement(angular_draws):
    offset = angular_draws + 2 * np.pi * np.arange(4)[:, None]
    assert az.rhat(offset) > 2
    assert rhat(offset, method="circular") < 1.01
    assert_allclose(rhat(offset, method="circular"), rhat(angular_draws, method="circular"))


def test_iid_same_direction_has_near_unit_rhat(angular_draws):
    assert 0.99 < rhat(angular_draws, method="circular") < 1.01


def test_detects_different_mean_directions():
    rng = np.random.default_rng(103)
    means = np.array([0, 0, 1.5, 1.5])[:, None]
    draws = rng.normal(means, 0.12, (4, 1200))
    assert rhat(draws, method="circular") > 2


def test_detects_different_concentrations_at_same_mean():
    rng = np.random.default_rng(103)
    scales = np.array([0.08, 0.08, 1.2, 1.2])[:, None]
    draws = rng.normal(0, scales, (4, 2400))
    assert rhat(draws, method="circular") > 1.1


def test_split_detects_nonstationarity_with_matching_whole_chain_means():
    rng = np.random.default_rng(103)
    draws = np.concatenate(
        [rng.normal(0, 0.12, (4, 600)), rng.normal(1.5, 0.12, (4, 600))], axis=1
    )
    assert az.rhat(draws, method="identity") < 1.01
    assert rhat(draws, method="circular") > 2


def test_odd_length_split_drops_the_middle_draw(angular_draws):
    odd = np.insert(angular_draws, angular_draws.shape[1] // 2, 1.7, axis=1)
    assert_allclose(rhat(odd, method="circular"), rhat(angular_draws, method="circular"))


def test_matches_independent_multivariate_covariance_reference(angular_draws):
    """Reference uses generic covariance/eigensolver, independent of circular formulas."""
    half = angular_draws.shape[1] // 2
    split = np.concatenate([angular_draws[:, :half], angular_draws[:, -half:]])
    vectors = np.stack([np.cos(split), np.sin(split)], axis=-1)
    chains, draws, _ = vectors.shape
    within = np.mean([np.cov(chain, rowvar=False, ddof=1) for chain in vectors], axis=0)
    mean_covariance = np.cov(vectors.mean(axis=1), rowvar=False, ddof=1)
    reference = np.sqrt(
        (draws - 1) / draws + (chains + 1) / chains * eigvalsh(mean_covariance, within)[-1]
    )
    assert_allclose(rhat(angular_draws, method="circular"), reference, rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize("scale", [1e-8, 1e-16, 1e-40, 1e-100, 1e-150])
@pytest.mark.parametrize("shifted", [False, True])
def test_tiny_angles_match_stable_small_angle_limit(scale, shifted):
    rng = np.random.default_rng(672)
    standardized = rng.normal(size=(4, 1200))
    if shifted:
        standardized += np.array([0, 0, 100, 100])[:, None]

    # As theta -> 0, joint (cos(theta)-1, sin(theta)) is, up to
    # invertible coordinate scaling, (z**2, z). Compute this reference with
    # O(1) inputs so it cannot lose tiny angular deviations to wrap roundoff.
    half = standardized.shape[1] // 2
    split = np.concatenate([standardized[:, :half], standardized[:, -half:]])
    features = np.stack([split**2, split], axis=-1)
    chains, draws, _ = features.shape
    within = np.mean([np.cov(chain, rowvar=False, ddof=1) for chain in features], axis=0)
    mean_covariance = np.cov(features.mean(axis=1), rowvar=False, ddof=1)
    reference = np.sqrt(
        (draws - 1) / draws + (chains + 1) / chains * eigvalsh(mean_covariance, within)[-1]
    )
    actual = rhat(scale * standardized, method="circular")
    assert_allclose(actual, reference, rtol=1e-6, atol=1e-10)
    if shifted:
        assert actual > 10
    else:
        assert actual < 1.01


@pytest.mark.parametrize("scale", [1e-78, 1e-150])
def test_extreme_disagreement_does_not_overflow_to_apparent_agreement(scale):
    rng = np.random.default_rng(823)
    standardized = rng.normal(size=(2, 1000))
    draws = np.concatenate([scale * standardized, np.full((2, 1000), np.pi)])

    # Independent small-angle limit: the two noisy chains have features
    # (-delta**2 / 2, delta), while the other two remain fixed at the antipode.
    # After splitting, S_rr -> m/(m-1), and R * scale**2 approaches this value.
    deviations = standardized - standardized[0, 0]
    split = np.concatenate([deviations[:, :500], deviations[:, 500:]])
    features = np.stack([-0.5 * split**2, split], axis=-1)
    within = np.mean([np.cov(chain, rowvar=False) for chain in features], axis=0) / 2
    reference = np.sqrt(9 / 7 * np.linalg.inv(within)[0, 0])
    with np.errstate(over="raise", invalid="raise"):
        actual = rhat(draws, method="circular")
    assert np.isfinite(actual)
    assert actual > 1e100
    assert_allclose(actual * scale**2, reference, rtol=1e-12)


@pytest.mark.parametrize("per_draw", [False, True])
def test_float32_period_equivalent_constants_respect_input_precision(per_draw):
    rng = np.random.default_rng(73)
    turns = rng.integers(-5, 6, (4, 100)) if per_draw else np.arange(4)[:, None]
    values = np.broadcast_to(0.4 + 2 * np.pi * turns, (4, 100)).astype(np.float32)
    assert np.isnan(rhat(values, method="circular"))
    distinct = values + np.array([0, 0, 0.8, 0.8], dtype=np.float32)[:, None]
    assert np.isinf(rhat(distinct, method="circular"))


def test_large_chain_count_preserves_representable_extreme_rhat():
    standardized = np.random.default_rng(5).normal(size=(8192, 8))
    scaled_results = []
    for scale in (2e-150, 2e-155):
        draws = scale * standardized
        draws[-1] = np.pi
        actual = rhat(draws, method="circular")
        assert np.isfinite(actual)
        assert actual > 1e290
        # The small-angle covariance ratio grows as scale**-4, hence R as
        # scale**-2. Apply the scale in stages to avoid underflow in the test.
        scaled_results.append(actual * scale * scale)
    assert_allclose(scaled_results[0], scaled_results[1], rtol=1e-12)


def test_first_moment_diagnostic_can_miss_different_antipodal_mode_pairs():
    """Document the limitation instead of treating near-unit R-hat as proof."""
    rng = np.random.default_rng(132)
    chains = []
    for rotation in (0.0, 0.0, np.pi / 2, np.pi / 2):
        halves = []
        for _ in range(2):
            jitter = rng.normal(0, 0.06, 600)
            half = np.concatenate([jitter, jitter + np.pi]) + rotation
            rng.shuffle(half)
            halves.append(half)
        chains.append(np.concatenate(halves))
    draws = np.asarray(chains)
    # Both first moments vanish. Higher angular moments distinguish the chains.
    second_moment = np.cos(2 * draws).mean(axis=1)
    assert np.all(second_moment[:2] > 0.98)
    assert np.all(second_moment[2:] < -0.98)
    assert rhat(draws, method="circular") < 1.01


def test_component_rhats_delegate_to_arviz_rank(angular_draws):
    result = rhat_components(angular_draws)
    assert set(result) == {"sin", "cos"}
    for name, transform in (("sin", np.sin), ("cos", np.cos)):
        assert np.isscalar(result[name])
        assert_array_equal(result[name], az.rhat(transform(angular_draws), method="rank"))


def test_component_projection_preserves_float32_arviz_behavior(angular_draws):
    draws = angular_draws.astype(np.float32)
    result = rhat_components(draws)
    for name, transform in (("sin", np.sin), ("cos", np.cos)):
        assert_array_equal(result[name], az.rhat(transform(draws), method="rank"))


@pytest.mark.parametrize("container", ["dataset", "dataarray", "inferencedata", "dict"])
def test_component_labeled_inputs_match_arviz(posterior, container):
    if container == "dataarray":
        data = posterior.direction
    elif container == "inferencedata":
        data = az.InferenceData(posterior=posterior, prior=posterior * 0)
    elif container == "dict":
        data = {"direction": posterior.direction.values, "offset": posterior.offset.values}
    else:
        data = posterior
    result = rhat_components(data)
    selected = az.extract(data, combined=False, keep_dataset=True)
    for name, transform in (("sin", np.sin), ("cos", np.cos)):
        assert isinstance(result[name], xr.Dataset)
        expected = az.rhat(transform(selected), method="rank")
        if container == "dict":
            # ArviZ stamps each independent dict conversion with created_at.
            xr.testing.assert_equal(result[name], expected)
        else:
            xr.testing.assert_identical(result[name], expected)


def test_component_selection_transposed_dimensions_and_coords(posterior):
    transposed = posterior.transpose("animal", "draw", "chain", missing_dims="ignore")
    included = rhat_components(transposed, var_names="direction")
    excluded = rhat_components(transposed, var_names=["~offset"])
    for name, transform in (("sin", np.sin), ("cos", np.cos)):
        xr.testing.assert_identical(included[name], excluded[name])
        xr.testing.assert_identical(included[name], az.rhat(transform(posterior[["direction"]])))
        assert_array_equal(included[name].animal, posterior.animal)


def test_component_transforms_only_selected_variables(posterior, caplog):
    # Variable selection must happen before trigonometric transforms. A caller
    # can retain non-angular metadata and concentrations in the same Dataset.
    mixed = posterior.assign(kappa=posterior.offset * np.nan, label=posterior.offset.astype(str))
    with caplog.at_level(logging.WARNING):
        result = rhat_components(mixed, var_names="direction")
    for component in result.values():
        assert list(component.data_vars) == ["direction"]
        assert np.isfinite(component.direction).all()
    assert not caplog.records


def test_component_full_turn_offsets_leave_diagnostics_unchanged(angular_draws):
    rng = np.random.default_rng(487)
    offsets = 2 * np.pi * rng.integers(-4, 5, angular_draws.shape)
    original = rhat_components(angular_draws)
    shifted = rhat_components(angular_draws + offsets)
    for name, transform in (("sin", np.sin), ("cos", np.cos)):
        assert_array_equal(shifted[name], az.rhat(transform(angular_draws + offsets)))
        # Trigonometric roundoff can resolve exact ties in folded ranks.
        assert_allclose(shifted[name], original[name], rtol=1e-5, atol=1e-12)


@pytest.mark.parametrize("invalid", [np.nan, np.inf, -np.inf])
def test_component_nonfinite_angles_return_nan(angular_draws, invalid):
    angular_draws[0, 0] = invalid
    with np.errstate(invalid="ignore"):
        result = rhat_components(angular_draws)
    assert all(np.isnan(value) for value in result.values())


def test_constant_projection_keeps_arviz_nan_behavior():
    # Cosine is constant even though the underlying direction varies.
    values = np.tile([-np.pi / 4, np.pi / 4], (4, 600))
    with np.errstate(invalid="ignore", divide="ignore"):
        result = rhat_components(values)
    assert np.isfinite(result["sin"])
    assert np.isnan(result["cos"])


def test_component_dask_matches_eager(posterior):
    pytest.importorskip("dask.array")
    lazy = posterior.chunk({"chain": -1, "draw": -1, "animal": 1})
    result = rhat_components(
        lazy,
        var_names="direction",
        dask_kwargs={"dask": "parallelized", "output_dtypes": [float]},
    )
    eager = rhat_components(posterior, var_names="direction")
    for name in result:
        assert hasattr(result[name].direction.data, "dask")
        xr.testing.assert_allclose(result[name].compute(), eager[name])


def test_component_marginals_can_agree_while_joint_angular_patterns_differ():
    """Equal sine/cosine marginals do not determine their joint distribution."""
    rng = np.random.default_rng(311)
    chains = []
    for center in (np.pi / 4, np.pi / 4, 3 * np.pi / 4, 3 * np.pi / 4):
        halves = []
        for _ in range(2):
            jitter = rng.uniform(0.005, 0.06, 300)
            half = np.concatenate([jitter, -jitter, jitter + np.pi, -jitter + np.pi]) + center
            rng.shuffle(half)
            halves.append(half)
        chains.append(np.concatenate(halves))
    values = np.asarray(chains)
    result = rhat_components(values)
    assert max(result.values()) < 1.01
    second_moment = np.sin(2 * values).mean(axis=1)
    assert np.all(second_moment[:2] > 0.99)
    assert np.all(second_moment[2:] < -0.99)


def test_averaging_component_rhats_can_hide_a_projection_disagreement():
    """A scalar average is not an appropriate replacement for reporting both."""
    rng = np.random.default_rng(841)
    noise = rng.normal(0, 0.1, 8000)
    noise -= noise.mean()
    # Reflection about pi/2 preserves the sine marginal but shifts the cosine
    # marginal. The shared finite pool isolates this representation effect.
    values = np.array(
        [rng.permutation(np.pi / 2 + sign * (noise + 0.016)) for sign in (1, 1, -1, -1)]
    )
    result = rhat_components(values)
    assert result["cos"] > 1.01
    assert np.mean(list(result.values())) < 1.01
