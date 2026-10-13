"""Compare custom likelihoods with independent SciPy/PyMC references."""
import numpy as np
import pymc as pm
import pytest
from numpy.testing import assert_allclose
from scipy.special import logsumexp
from scipy.stats import vonmises

from unwrap.distributions import CircularUniform, VonMisesMixture, VonMisesSS, vmss_logp


@pytest.mark.parametrize("mu,kappa", [(0.0, 0.1), (2.9, 2.0), (-2.8, 20.0)])
def test_sufficient_statistics_match_raw_observations(mu, kappa):
    mu, kappa = np.float64(mu), np.float64(kappa)
    angles = np.array([[-3.0, 3.1, 2.8], [-0.5, 0.0, 0.4]])
    vectors = np.exp(1j * angles).mean(axis=1)
    theta, resultant = np.angle(vectors), np.abs(vectors)
    counts = np.full(2, angles.shape[1])
    expected = vonmises.logpdf(angles, kappa, loc=mu).sum(axis=1)
    assert_allclose(vmss_logp(theta, mu, kappa, resultant, counts).eval(), expected)
    with pm.Model() as model:
        VonMisesSS("angles", mu=mu, kappa=kappa, R=resultant, n=counts, observed=theta)
    assert_allclose(model.compile_logp()({}), expected.sum())


def test_sufficient_statistics_reject_angles_outside_support():
    assert np.isneginf(vmss_logp(np.pi + 0.1, 0, 1, 0.5, 3).eval())


def test_mixture_matches_weighted_densities():
    weights = np.array([0.25, 0.75])
    means = np.array([-1.0, 1.2])
    kappas = np.array([2.0, 5.0])
    angles = np.linspace(-np.pi, np.pi, 9)
    expected = logsumexp(np.log(weights) + vonmises.logpdf(
        angles[:, None], kappas, loc=means), axis=-1)
    mixture = VonMisesMixture.dist(w=weights, mu=means, kappa=kappas)
    assert_allclose(pm.logp(mixture, angles).eval(), expected)


def test_circular_uniform_has_constant_log_density():
    with pm.Model():
        prior = CircularUniform("direction")
    assert_allclose(pm.logp(prior, np.linspace(-np.pi, np.pi, 9)).eval(),
                    -np.log(2 * np.pi), atol=1e-12)
