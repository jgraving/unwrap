"""Concentration/resultant-length conversion on its documented domain."""
import numpy as np
import pytest
from unwrap.bessel import bessel_ratio, inverse_bessel_ratio


@pytest.mark.parametrize("kappa", [0, 0.1, 0.5, 1, 2, 10, 50, 100])
def test_concentration_roundtrip(kappa):
    resultant = bessel_ratio(kappa)
    assert 0 <= resultant < 1
    assert inverse_bessel_ratio(resultant) == pytest.approx(kappa, rel=1e-6, abs=1e-8)


def test_resultant_limits():
    assert bessel_ratio(0) == 0
    assert bessel_ratio(np.inf) == 1
    assert inverse_bessel_ratio(0) == 0


@pytest.mark.parametrize("resultant", [-0.1, 1, 1.1, np.nan, np.inf])
def test_invalid_resultant_rejected(resultant):
    with pytest.raises(ValueError, match=r"R must be in"):
        inverse_bessel_ratio(resultant)


@pytest.mark.parametrize("kappa", [1000.0, 1e4, 1e6])
def test_high_concentration_is_finite_and_invertible(kappa):
    resultant = bessel_ratio(kappa)
    assert np.isfinite(resultant) and 0 < resultant < 1
    # Large-kappa expansion: 1 - 1/(2*kappa) - 1/(8*kappa**2).
    expected = 1 - 0.5 / kappa - 0.125 / kappa**2
    assert resultant == pytest.approx(expected, abs=2e-10)
    assert inverse_bessel_ratio(resultant) == pytest.approx(kappa, rel=1e-7)


def test_small_positive_resultant_is_not_rounded_to_zero():
    # I1(k)/I0(k) ~ k/2 near zero.
    assert inverse_bessel_ratio(1e-12) == pytest.approx(2e-12, rel=1e-8, abs=0)


@pytest.mark.parametrize("kappa", [-1, -np.inf, np.nan])
def test_invalid_concentration_rejected(kappa):
    with pytest.raises(ValueError, match="kappa"):
        bessel_ratio(kappa)


@pytest.mark.parametrize("kwargs", [
    {"tolerance": 0}, {"tolerance": -1}, {"tolerance": np.nan},
    {"tolerance": np.inf}, {"max_iterations": 0}, {"max_iterations": 1.5},
])
def test_invalid_solver_controls_rejected(kwargs):
    with pytest.raises(ValueError):
        inverse_bessel_ratio(0.5, **kwargs)


def test_solver_reports_nonconvergence():
    with pytest.raises(RuntimeError, match="converge"):
        inverse_bessel_ratio(0.7, max_iterations=1)
