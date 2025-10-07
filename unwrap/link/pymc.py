import pytensor.tensor as pt
import pymc as pm
import numpy as np
from unwrap.link.utils import LinksContainer

PI = np.pi
EPS = np.finfo(float).eps


def inverse_softplus(mu):
    return pt.log(pt.expm1(mu))


def softplus(eta):
    return pm.math.log1pexp(eta)


def reciprocal_softplus(eta):
    return pt.exp(-pt.log(softplus(eta)))


def inverse_reciprocal_softplus(mu):
    return inverse_softplus(pt.exp(-pt.log(mu)))


def modulo(eta):
    return eta % (2 * PI)


def rotated_modulo(eta):
    return ((eta + PI) % (2 * PI)) - PI

shifted_modulo = rotated_modulo

def rectifier(eta):
    return pt.maximum(EPS, eta)


def atan2(eta):
    return pt.arctan2(pt.sin(eta), pt.cos(eta))


LINKS = {
    "softplus": LinksContainer(inverse_softplus, softplus),
    "logexpm1": LinksContainer(inverse_softplus, softplus),  # Alias for softplus
    "log1pexp": LinksContainer(inverse_softplus, softplus),  # Alias for softplus
    "reciprocal_softplus": LinksContainer(
        inverse_reciprocal_softplus, reciprocal_softplus
    ),
    "inv_softplus": LinksContainer(
        inverse_reciprocal_softplus, reciprocal_softplus
    ),  # Short alias
    "modulo": LinksContainer(modulo, modulo),
    "mod": LinksContainer(modulo, modulo),  # Alias for modulo
    "shifted_modulo": LinksContainer(rotated_modulo, rotated_modulo),
    "shift_mod": LinksContainer(
        rotated_modulo, rotated_modulo
    ),  # Alias for rotated_modulo
    "rotated_modulo": LinksContainer(rotated_modulo, rotated_modulo),
    "rot_mod": LinksContainer(
        rotated_modulo, rotated_modulo
    ),  # Alias for rotated_modulo
    "rectifier": LinksContainer(rectifier, rectifier),
    "relu": LinksContainer(rectifier, rectifier),  # Alias for rectifier
    "atan2": LinksContainer(atan2, atan2),
}
