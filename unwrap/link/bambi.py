"""Bambi adapters combining numerical and symbolic link implementations.

``LINKS[name]`` is a ``bambi.Link``. Each registry key is also exported as a
module attribute (for example ``circular_modulo``); these attributes are Link
objects, not the functions with the same names in the other backends.

Alias keys keep their own Bambi ``name`` while sharing transform behaviour.
"""
import bambi as bmb
from unwrap.link.numpy import LINKS as NP_LINKS
from unwrap.link.pymc import LINKS as PYMC_LINKS

common_keys = NP_LINKS.keys() & PYMC_LINKS.keys()

LINKS = {
    name: bmb.Link(
        name=name,
        link=NP_LINKS[name].link,
        linkinv=NP_LINKS[name].linkinv,
        linkinv_backend=PYMC_LINKS[name].linkinv,
    )
    for name in common_keys
}

# Export each link as a module-level variable
globals().update(LINKS)

# Optional: define __all__ for clean imports
__all__ = list(LINKS.keys()) + ["LINKS"]
