"""Regression tests for the publication figure's paired mixture summaries."""
import sys
from pathlib import Path

import numpy as np
import pytest
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
    np.testing.assert_array_equal(data["undirected_density_profile"], np.full_like(angle, uniform))
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

        contributions = data["component_contributions"][level]
        np.testing.assert_allclose(contributions["directed"], (weight * directed).mean(axis=0))
        np.testing.assert_allclose(contributions["undirected"], (1 - weight.mean()) * uniform)
        np.testing.assert_allclose(contributions["directed"] + contributions["undirected"],
                                   data["density_profiles"][level]["mean"])
        assert not np.allclose(contributions["directed"], weight.mean() * directed.mean(axis=0))


def test_paired_contrasts_and_circular_intervals(predictions):
    data, mu, kappa, omega = predictions
    contrasts = data["component_contrasts"]
    np.testing.assert_allclose(contrasts[2]["omega"], [0.1, 0.1, 0.1])
    np.testing.assert_allclose(contrasts[2]["kappa"], [2, 2, 2])
    np.testing.assert_allclose(contrasts[2]["mu_deg"], [4, 3, 5])
    mean, start, end = contrasts[3]["mu_deg"]
    np.testing.assert_allclose([abs(mean), start, end], [180, 179, -179])
    np.testing.assert_allclose(data["component_parameters"][1]["omega"], [0.5, 0.2, 0.8])
    np.testing.assert_allclose(data["component_parameters"][1]["kappa"], [5, 1, 9])
    mean, start, end = data["component_parameters"][1]["mu_deg"]
    np.testing.assert_allclose([abs(mean), start, end], [180, 178, -178])
    full_at_goal = omega * vonmises.pdf(0.0, kappa, loc=mu) + (1 - omega) / (2 * np.pi)
    for i, level in enumerate([2, 3, 4], start=1):
        paired = full_at_goal[:, i] - full_at_goal[:, 0]
        np.testing.assert_allclose(data["full_contrasts"][level],
                                   [paired.mean(), paired.min(), paired.max()], atol=1e-15)
