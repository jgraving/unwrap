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
