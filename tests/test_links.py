"""Mathematical invariants and agreement between numerical backends."""
import numpy as np
import pytest
import pytensor.tensor as pt
from numpy.testing import assert_allclose

from unwrap.link import numpy as np_links
from unwrap.link import pymc as pt_links
from unwrap.link.bambi import LINKS


@pytest.mark.parametrize("name", ["modulo", "circular_modulo", "atan2"])
def test_wrapping_preserves_direction_and_periodicity(name):
    angles = np.array([-20.1, -np.pi, -0.1, 0, np.pi, 19.3])
    wrap = getattr(np_links, name)
    wrapped = wrap(angles)
    assert_allclose(np.exp(1j * wrapped), np.exp(1j * angles), atol=1e-13)
    assert_allclose(np.exp(1j * wrap(angles + 8 * np.pi)),
                    np.exp(1j * wrapped), atol=1e-13)
    lower, upper = (0, 2 * np.pi) if name == "modulo" else (-np.pi, np.pi)
    assert np.all((wrapped >= lower) & (wrapped <= upper))


def test_circular_modulo_boundary_convention():
    assert_allclose(np_links.circular_modulo(np.array([-np.pi, np.pi, 3 * np.pi])), -np.pi)


@pytest.mark.parametrize("name", sorted(np_links.LINKS))
def test_backend_and_bambi_link_agreement(name):
    eta = np.array([-6.0, -0.2, 0.0, 0.8, 6.0])
    expected = np_links.LINKS[name].linkinv(eta)
    assert_allclose(pt_links.LINKS[name].linkinv(pt.as_tensor_variable(eta)).eval(), expected, atol=1e-12)
    assert_allclose(LINKS[name].linkinv(eta), expected, atol=1e-12)
    assert_allclose(LINKS[name].linkinv_backend(pt.as_tensor_variable(eta)).eval(), expected, atol=1e-12)


@pytest.mark.parametrize("name", ["softplus", "reciprocal_softplus"])
def test_positive_link_roundtrip(name):
    eta = np.linspace(-10, 10, 21)
    for backend in (np_links, pt_links):
        link = backend.LINKS[name]
        positive = link.linkinv(eta)
        result = link.link(positive)
        if backend is pt_links:
            positive, result = positive.eval(), result.eval()
        assert np.all(positive > 0)
        assert_allclose(result, eta, atol=1e-10)


@pytest.mark.parametrize("backend", [np_links, pt_links])
def test_softplus_extreme_inputs_stay_finite(backend):
    eta = np.array([-100.0, 0.0, 100.0, 1000.0])
    mu = backend.softplus(eta)
    if backend is pt_links:
        mu = mu.eval()
    assert np.isfinite(mu).all()
    assert np.all(mu > 0)
    assert_allclose(mu, np.logaddexp(0, eta), rtol=1e-12)
    recovered = backend.inverse_softplus(mu)
    if backend is pt_links:
        recovered = recovered.eval()
    assert_allclose(recovered, eta, atol=1e-12)


def test_symbolic_inverse_softplus_gradient_at_large_input():
    mu = pt.dscalar("mu")
    derivative = pt.grad(pt_links.inverse_softplus(mu), mu)
    assert_allclose(derivative.eval({mu: 1000.0}), 1.0, atol=1e-12)
