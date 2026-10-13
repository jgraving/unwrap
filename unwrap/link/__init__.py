"""Links for circular means and positive concentration parameters.

The functions imported from ``unwrap.link`` use NumPy. For symbolic model
expressions use ``unwrap.link.pymc``; for Bambi use ``unwrap.link.bambi.LINKS``.
Importing this package currently loads all three backends.

Registry conventions
--------------------
``NUMPY_LINKS`` and ``PYMC_LINKS`` map names to ``LinksContainer(link, linkinv)``.
``BAMBI_LINKS`` maps the same names to Bambi Link objects; each combines NumPy
forward/inverse functions with a symbolic ``linkinv_backend`` for fitting.

Names and aliases
-----------------
* ``softplus``: aliases ``logexpm1``, ``log1pexp``.
* ``reciprocal_softplus``: alias ``inv_softplus`` (not inverse_softplus).
* ``modulo``: alias ``mod``; range [0, 2*pi).
* ``circular_modulo``: aliases ``circ_mod``, ``shifted_modulo``, ``shift_mod``,
  ``rotated_modulo``, ``rot_mod``; range [-pi, pi).
* ``rectifier``: alias ``relu``; floor at float64 machine epsilon.
* ``atan2``: unary angle wrapper with range [-pi, pi].

Aliases here are registry keys; not every key is a standalone NumPy function.
For the positive transforms, ``link`` takes a positive response parameter and
``linkinv`` takes an unconstrained predictor. Circular wrappers and rectifier
use the same function in both slots and have no mathematical inverse.

Examples
--------
>>> import numpy as np
>>> circular_modulo(np.array([0.0, 2 * np.pi]))
array([0., 0.])

Pass Bambi link objects to its existing von Mises family::

    from unwrap.link.bambi import LINKS
    links = {"mu": LINKS["circular_modulo"], "kappa": LINKS["softplus"]}
    # bmb.Model(..., family="vonmises", link=links)
"""

# Default: NumPy backend functions
from .numpy import (
    inverse_softplus,
    softplus,
    reciprocal_softplus,
    inverse_reciprocal_softplus,
    modulo,
    shifted_modulo,
    rotated_modulo,
    circular_modulo,
    rectifier,
    atan2,
    LINKS as NUMPY_LINKS,
)

# Backend-specific LINKS
from .bambi import LINKS as BAMBI_LINKS
from .pymc import LINKS as PYMC_LINKS

__all__ = [
    # Default NumPy link functions
    "inverse_softplus",
    "softplus",
    "reciprocal_softplus",
    "inverse_reciprocal_softplus",
    "modulo",
    "shifted_modulo",
    "rotated_modulo",
    "circular_modulo",
    "rectifier",
    "atan2",
    # Backend-specific LINKS dicts
    "NUMPY_LINKS",
    "BAMBI_LINKS",
    "PYMC_LINKS",
]
