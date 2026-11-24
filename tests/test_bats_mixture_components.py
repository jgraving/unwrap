"""Regression tests for population density profiles and individual goal-window probabilities."""
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


@pytest.fixture
def predictions(monkeypatch):
    # Two joint states with positively correlated directed weight and concentration.
    # The control mean bearing straddles the ±180° seam.
    state = np.arange(80) % 2
    mu = np.deg2rad(np.array([[178, -179, -3, 178], [-178, -173, 3, -178]])[state])
    kappa = np.array([[1, 3, 9, 1], [9, 11, 1, 9]])[state]
    omega = np.array([[0.2, 0.3, 0.8, 0.2], [0.8, 0.9, 0.2, 0.8]])[state]
    monkeypatch.setattr(layout, "_population_parameters", lambda *_: (mu, kappa, omega))
    levels = [1, 2, 3, 4]
    data = layout._mixture_predictions(None, levels, dict.fromkeys(levels, 1), posthoc)
    return data, mu[:2], kappa[:2], omega[:2]


def test_draw_weighted_component_densities_and_intervals(predictions):
    data, mu, kappa, omega = predictions
    angle = data["density_bearings"]
    uniform = 1 / (2 * np.pi)
    for i, level in enumerate([1, 2, 3, 4]):
        directed = vonmises.pdf(angle, kappa[:, i, None], loc=mu[:, i, None])
        weight = omega[:, i, None]
        mixture = weight * directed + (1 - weight) * uniform
        np.testing.assert_allclose(data["directed_density_profiles"][level]["mean"],
                                   directed.mean(axis=0))
        curve = data["density_profiles"][level]
        np.testing.assert_allclose(curve["mean"], mixture.mean(axis=0))
        np.testing.assert_allclose(curve["lo"], mixture.min(axis=0))
        np.testing.assert_allclose(curve["hi"], mixture.max(axis=0))

        assert trapezoid(curve["mean"], angle) == pytest.approx(1)
        assert trapezoid(data["directed_density_profiles"][level]["mean"], angle) == pytest.approx(1)


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
    directed_share = weight[:, tag_index]
    posterior = xr.Dataset({
        "mu": (("chain", "draw", "mu_dim_0"), mu[:, tag_index][None]),
        "kappa": (("chain", "draw", "kappa_dim_0"), kappa[:, tag_index][None]),
        "p_goal": (("chain", "draw", "p_goal_dim_0", "p_goal_dim_1"),
                   np.stack([directed_share, 1 - directed_share], axis=-1)[None]),
    })
    rows = layout._goal_window_by_bat(posterior, observed, posthoc)
    assert [row["tag"] for row in rows] == list(tags)
    assert [row["treatment"] for row in rows] == groups
    for i, row in enumerate(rows):
        draws = []
        for m, k, w in zip(mu[:, i], kappa[:, i], weight[:, i], strict=True):
            probability, _ = quad(lambda x: w * vonmises.pdf(x, k, loc=m) + (1 - w) / (2 * np.pi),
                                  -np.pi / 4, np.pi / 4)
            draws.append(probability)
        assert row["probability"] == pytest.approx(np.mean(draws), abs=1e-12)
    assert rows[2]["probability"] == pytest.approx(0.25)
    observed.loc[7, "treatment"] = 4
    with pytest.raises(ValueError, match="more than one treatment"):
        layout._goal_window_by_bat(posterior, observed, posthoc)
