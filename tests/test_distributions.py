"""Compare custom likelihoods with independent SciPy/PyMC references."""
import numpy as np
import pymc as pm
from numpy.testing import assert_allclose
from scipy.special import logsumexp
from scipy.stats import vonmises

from unwrap.distributions import CircularUniform, VonMisesMixture


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
