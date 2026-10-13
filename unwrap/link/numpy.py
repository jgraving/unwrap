"""NumPy circular and positive transforms for regression.

Angles are radians. In ``LINKS``, ``link`` maps a response parameter to
predictor scale and ``linkinv`` maps a predictor to response-parameter scale.
Circular wrappers and the rectifier are many-to-one: using them in both
slots is a convention, not a claim that they are invertible.

Inputs are real scalars or NumPy arrays; operations are elementwise. The corresponding
NumPy and PyTensor implementations share names and aliases. See the package
module docstring for backend selection and Bambi usage.
"""
import numpy as np
from unwrap.link.utils import LinksContainer

PI = np.pi
EPS = np.finfo(float).eps


def inverse_softplus(mu):
    """Map a positive value to the predictor scale: ``log(expm1(mu))``.

    ``mu`` must be positive for a finite real result. This is the forward
    ``link`` paired with ``softplus``. Zero maps to negative infinity;
    negative inputs are outside the real-valued domain."""
    mu = np.asarray(mu)
    return mu + np.log(-np.expm1(-mu))


def softplus(eta):
    """Map a real predictor to a positive value: ``log(1 + exp(eta))``.

    This is ``LINKS["softplus"].linkinv``, useful for concentration ``kappa``.
    Floating-point underflow can produce zero at sufficiently negative inputs."""
    return np.logaddexp(0, eta)


def reciprocal_softplus(eta):
    """Map a real predictor to ``1 / softplus(eta)``.

    This decreasing positive transform is paired with
    ``inverse_reciprocal_softplus``. It is not the inverse of softplus.
    Extreme inputs may exceed the floating-point range."""
    return np.exp(-np.log(softplus(eta)))


def inverse_reciprocal_softplus(mu):
    """Map a positive value to ``inverse_softplus(1 / mu)``.

    This is the forward link for ``reciprocal_softplus``; ``mu > 0``."""
    return inverse_softplus(np.exp(-np.log(mu)))


def modulo(eta):
    """Wrap radians to the half-open interval ``[0, 2*pi)``.

    Directions differing by an integer number of turns map to the same value.
    The mapping is discontinuous at multiples of ``2*pi`` and is not invertible."""
    return eta % (2 * PI)


def circular_modulo(eta):
    """Wrap radians to the half-open interval ``[-pi, pi)``.

    In particular, both ``-pi`` and ``pi`` map to ``-pi``. The mapping is
    discontinuous at odd multiples of ``pi`` and is not invertible.
    ``shifted_modulo`` and ``rotated_modulo`` are exact function aliases;
    they do not accept a shift or rotation parameter."""
    return ((eta + PI) % (2 * PI)) - PI

shifted_modulo = circular_modulo
rotated_modulo = circular_modulo

def rectifier(eta):
    """Floor values at double-precision machine epsilon.

    Returns ``max(eps, eta)``, unlike a zero-floor ReLU. This many-to-one
    map is used in both link slots; applying it cannot recover clipped inputs."""
    return np.maximum(EPS, eta)


def atan2(eta):
    """Represent a direction as ``atan2(sin(eta), cos(eta))`` in radians.

    Results lie in ``[-pi, pi]``. The sign of a boundary value can depend on
    floating-point rounding, unlike the fixed half-open modulo convention.
    This unary angle wrapper is not NumPy's two-argument ``arctan2``."""
    return np.arctan2(np.sin(eta), np.cos(eta))


LINKS = {
    "softplus": LinksContainer(inverse_softplus, softplus),
    "logexpm1": LinksContainer(inverse_softplus, softplus),  # Alias for softplus
    "log1pexp": LinksContainer(inverse_softplus, softplus),  # Alias for softplus
    "reciprocal_softplus": LinksContainer(inverse_reciprocal_softplus, reciprocal_softplus),
    "inv_softplus": LinksContainer(inverse_reciprocal_softplus, reciprocal_softplus),  # Alias for reciprocal_softplus
    "modulo": LinksContainer(modulo, modulo),
    "mod": LinksContainer(modulo, modulo),  # Alias for modulo
    "shifted_modulo": LinksContainer(shifted_modulo, shifted_modulo),
    "shift_mod": LinksContainer(shifted_modulo, shifted_modulo),  # Alias for circmod
    "rotated_modulo": LinksContainer(rotated_modulo, rotated_modulo),
    "rot_mod": LinksContainer(rotated_modulo, rotated_modulo),  # Alias for circmod
    "circular_modulo": LinksContainer(circular_modulo, circular_modulo),
    "circ_mod": LinksContainer(circular_modulo, circular_modulo),  # Alias for circmod
    "rectifier": LinksContainer(rectifier, rectifier),
    "relu": LinksContainer(rectifier, rectifier),  # Alias for rectifier
    "atan2": LinksContainer(atan2, atan2),
}
