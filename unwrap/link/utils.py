"""Small shared types used by the numerical link registries."""
from collections import namedtuple

LinksContainer = namedtuple("LinksContainer", ["link", "linkinv"])

LinksContainer.__doc__ = "Pair of parameter-to-predictor and predictor-to-parameter functions."
LinksContainer.link.__doc__ = "Forward link: response parameter to linear predictor."
LinksContainer.linkinv.__doc__ = "Inverse link: linear predictor to response parameter."

# Function names are shared; implementations remain backend-specific.
_TRANSFORMS = {
    "softplus": ("inverse_softplus", "softplus"),
    "reciprocal_softplus": ("inverse_reciprocal_softplus", "reciprocal_softplus"),
    "modulo": ("modulo", "modulo"),
    "circular_modulo": ("circular_modulo", "circular_modulo"),
    "rectifier": ("rectifier", "rectifier"),
    "atan2": ("atan2", "atan2"),
}

# Preserve the established registry names and insertion order.
_NAMES = {
    "softplus": "softplus",
    "logexpm1": "softplus",
    "log1pexp": "softplus",
    "reciprocal_softplus": "reciprocal_softplus",
    "inv_softplus": "reciprocal_softplus",
    "modulo": "modulo",
    "mod": "modulo",
    "shifted_modulo": "circular_modulo",
    "shift_mod": "circular_modulo",
    "rotated_modulo": "circular_modulo",
    "rot_mod": "circular_modulo",
    "circular_modulo": "circular_modulo",
    "circ_mod": "circular_modulo",
    "rectifier": "rectifier",
    "relu": "rectifier",
    "atan2": "atan2",
}


def _build_links(namespace):
    """Bind the shared registry to one backend's named transform functions."""
    return {
        name: LinksContainer(*(namespace[fn] for fn in _TRANSFORMS[canonical]))
        for name, canonical in _NAMES.items()
    }
