"""Regression tests for full-mixture heading densities and goal-sector probabilities."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import xarray as xr
from scipy.integrate import quad, trapezoid
from scipy.stats import vonmises

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "notebooks"))
import goldshtein_bats_layout_data as layout
import goldshtein_bats_posthoc as posthoc


# The first state alone contains the 90% HDI. A 95% interval must include both.
STATE_COUNTS = [74, 6]


def expected_summary(draws):
    return [np.average(draws, weights=STATE_COUNTS), draws[0], draws[0]]


@pytest.fixture
def predictions(monkeypatch):
    # Two joint states with positively correlated directed weight and concentration.
    # The control mean bearing straddles the ±180° seam.
    state = np.repeat([0, 1], STATE_COUNTS)
    mu = np.deg2rad(np.array([[178, -179, -3, 178], [-178, -173, 3, -178]])[state])
    kappa = np.array([[1, 3, 9, 1], [9, 11, 1, 9]])[state]
    omega = np.array([[0.2, 0.3, 0.8, 0.2], [0.8, 0.9, 0.2, 0.8]])[state]
    monkeypatch.setattr(layout, "_population_parameters", lambda *_: (mu, kappa, omega))
    levels = [1, 2, 3, 4]
    data = layout._mixture_predictions(None, levels, dict.fromkeys(levels, 1), posthoc)
    return data, mu[[0, -1]], kappa[[0, -1]], omega[[0, -1]]


def test_draw_weighted_component_densities_and_intervals(predictions):
    data, mu, kappa, omega = predictions
    angle = data["density_bearings"]
    uniform = 1 / (2 * np.pi)
    goal_draws = []
    for i, level in enumerate([1, 2, 3, 4]):
        directed = vonmises.pdf(angle, kappa[:, i, None], loc=mu[:, i, None])
        weight = omega[:, i, None]
        mixture = weight * directed + (1 - weight) * uniform
        curve = data["density_profiles"][level]
        np.testing.assert_allclose(curve["mean"], np.average(mixture, axis=0, weights=STATE_COUNTS))
        np.testing.assert_allclose(curve["lo"], mixture[0])
        np.testing.assert_allclose(curve["hi"], mixture[0])

        assert trapezoid(curve["mean"], angle) == pytest.approx(1)
        goal = mixture[:, angle == 0].ravel()
        goal_draws.append(goal)
        np.testing.assert_allclose(data["density_at_zero_by_treatment"][level], expected_summary(goal))
    for level, draws in enumerate(goal_draws[1:], start=2):
        difference = draws - goal_draws[0]
        np.testing.assert_allclose(data["full_contrasts"][level], expected_summary(difference), atol=1e-14)


def sector_integral(mu, kappa, weight):
    """Independent quadrature of both components over the specified ±45° sector."""
    return quad(lambda x: weight * vonmises.pdf(x, kappa, loc=mu) + (1 - weight) / (2 * np.pi),
                -np.pi / 4, np.pi / 4)[0]


def test_population_sector_contrasts_are_paired_in_each_draw(predictions):
    data, mu, kappa, omega = predictions
    assert layout.PANEL_C_HDI_PROB == 0.90
    assert data["goal_window_half_width_degrees"] == 45
    group_draws = []
    for i, level in enumerate([1, 2, 3, 4]):
        draws = np.array([sector_integral(m, k, w)
                          for m, k, w in zip(mu[:, i], kappa[:, i], omega[:, i], strict=True)])
        group_draws.append(draws)
        np.testing.assert_allclose(data["goal_window_groups"][level], expected_summary(draws), atol=1e-14)
    for level, draws in enumerate(group_draws[1:], start=2):
        delta = draws - group_draws[0]
        np.testing.assert_allclose(data["goal_window_contrasts"][level], expected_summary(delta), atol=1e-14)


def test_each_fitted_bat_is_counted_once_in_its_observed_treatment():
    # A and B share a treatment but differ in their individual effects. Each is
    # observed a different number of times, and labels appear in non-model order.
    tag_index = np.array([1, 0, 1, 2, 0, 1])
    tags = np.array(["a", "b", "c"])
    groups = [1, 1, 3]
    observed = pd.DataFrame({"tag": tags[tag_index], "treatment": np.array(groups)[tag_index]},
                            index=[29, 13, 7, 31, 9, 3])
    mu = np.array([[0.1, 1.4, 2.9], [-0.5, 0.1, -2.9]])
    kappa = np.array([[9., 2., 5.], [2., 9., 4.]])
    weight = np.array([[0.9, 0.2, 0.], [0.3, 0.8, 0.]])
    state = np.repeat([0, 1], STATE_COUNTS)
    directed_share = weight[state][:, tag_index]
    posterior = xr.Dataset({
        "mu": (("chain", "draw", "mu_dim_0"), mu[state][:, tag_index][None]),
        "kappa": (("chain", "draw", "kappa_dim_0"), kappa[state][:, tag_index][None]),
        "p_goal": (("chain", "draw", "p_goal_dim_0", "p_goal_dim_1"),
                   np.stack([directed_share, 1 - directed_share], axis=-1)[None]),
    })
    rows = layout._goal_window_by_bat(posterior, observed, posthoc)
    assert [row["tag"] for row in rows] == list(tags)
    assert [row["treatment"] for row in rows] == groups
    for i, row in enumerate(rows):
        draws = np.array([sector_integral(m, k, w)
                          for m, k, w in zip(mu[:, i], kappa[:, i], weight[:, i], strict=True)])
        np.testing.assert_allclose([row["mean"], row["lo"], row["hi"]], expected_summary(draws))
    assert rows[2]["mean"] == pytest.approx(0.25)
    observed.loc[7, "treatment"] = 4
    with pytest.raises(ValueError, match="more than one treatment"):
        layout._goal_window_by_bat(posterior, observed, posthoc)
