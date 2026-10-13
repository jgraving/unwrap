"""Small shared types used by the numerical link registries."""
from collections import namedtuple

LinksContainer = namedtuple("LinksContainer", ["link", "linkinv"])

LinksContainer.__doc__ = "Pair of parameter-to-predictor and predictor-to-parameter functions."
LinksContainer.link.__doc__ = "Forward link: response parameter to linear predictor."
LinksContainer.linkinv.__doc__ = "Inverse link: linear predictor to response parameter."
