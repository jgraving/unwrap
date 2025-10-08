# unwrap/link/__init__.py

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
