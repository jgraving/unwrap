"""Bimodal inputs retain the same weighted von Mises density."""
import numpy as np
import pymc as pm
import pytest
import pytensor
import pytensor.tensor as pt
from numpy.testing import assert_allclose
from scipy.special import logsumexp
from scipy.stats import vonmises

from unwrap.distributions import BimodalVonMises


def reference(y, mu, kappa, p):
    weights = np.stack([p, 1 - np.asarray(p)], axis=-1)
    return logsumexp(np.log(weights) + vonmises.logpdf(
        np.asarray(y)[..., None], kappa, loc=mu), axis=-1)


@pytest.mark.parametrize('convert', [np.asarray, list, pt.as_tensor_variable])
def test_component_input_types(convert):
    mu, kappa = convert([-1., 1.2]), convert([2., 5.])
    y = np.linspace(-np.pi, np.pi, 7)
    rv = BimodalVonMises.dist(mu, kappa, .25)
    assert_allclose(pm.logp(rv, y).eval(), reference(y, [-1., 1.2], [2., 5.], .25))


def test_symbolic_unknown_shape_and_gradient():
    mu, kappa = pt.vector('mu'), pt.vector('kappa')
    rv = BimodalVonMises.dist(mu, kappa, .25)
    lp = pm.logp(rv, .3)
    f = pytensor.function([mu, kappa], [lp, pt.grad(lp, mu)])
    means, concentrations = np.array([-1., 1.2]), np.array([2., 5.])
    actual, grad = f(means, concentrations)
    assert_allclose(actual, reference(.3, means, concentrations, .25))
    h = 1e-6
    expected_grad = [(reference(.3, means + h * e, concentrations, .25)
                      - reference(.3, means - h * e, concentrations, .25)) / (2*h)
                     for e in np.eye(2)]
    assert_allclose(grad, expected_grad, atol=1e-7)
    # A one-component axis must not silently broadcast to two components.
    for bad_mu, bad_kappa in [([0.], concentrations), (means, [2.]),
                              ([0., 1., 2.], concentrations)]:
        with pytest.raises((AssertionError, ValueError)):
            f(bad_mu, bad_kappa)


@pytest.mark.parametrize('mu,kappa', [(0., [2., 3.]), ([0.], [2., 3.]),
                                      ([0., 1.], 2.), ([0., 1.], [2., 3., 4.])])
def test_invalid_known_component_shapes(mu, kappa):
    with pytest.raises(ValueError, match='last dimension of size 2'):
        BimodalVonMises.dist(mu, kappa, .5)


def test_batched_parameters_and_probability_list():
    mu = np.array([[0., 1.], [.5, -1.], [-.4, 2.]])
    kappa = np.array([[2., 3.], [1., 4.], [3., 2.]])
    p = [.2, .5, .8]
    y = np.array([.1, .2, .3])
    rv = BimodalVonMises.dist(mu, kappa, p)
    assert_allclose(pm.logp(rv, y).eval(), reference(y, mu, kappa, p))
    draws = pm.draw(rv, draws=4, random_seed=42)
    assert draws.shape == (4, 3)
    assert np.isfinite(draws).all()


def test_named_distribution_with_random_parameters():
    with pm.Model(coords={'observation': range(3)}) as model:
        mu = pm.Normal('mu', shape=2)
        kappa = pm.HalfNormal('kappa', shape=2)
        p = pm.Beta('p', 2, 2)
        BimodalVonMises('y', mu, kappa, p, observed=[.1, .2, .3], dims='observation')
    point = model.initial_point()
    assert np.isfinite(model.compile_logp()(point))
    assert np.isfinite(model.compile_dlogp()(point)).all()
