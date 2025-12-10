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
