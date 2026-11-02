"""Colab bootstrap helpers for the example notebooks.

The "Open in Colab" button on a notebook page starts the notebook
straight from GitHub, so the Colab runtime has no repository checkout
and ``unwrap`` is not installed. The first code cell of each example
notebook bootstraps the package from its git URL and then calls
:func:`install_colab` to finish the environment.
"""

import subprocess
import sys
from collections.abc import Sequence
from importlib.metadata import PackageNotFoundError, requires

__all__ = ["install_colab", "notebook_extra_requirements"]

_EXTRA = "notebook"


def _normalize_marker(marker: str) -> str:
    return marker.replace('"', "'").replace(" ", "")


def notebook_extra_requirements() -> list[str]:
    """Return the installed distribution's ``notebook`` extra as plain
    requirement strings (no ``extra ==`` marker).
    """
    try:
        dist_requires = requires("unwrap") or []
    except PackageNotFoundError:
        return []
    return [
        requirement.partition(";")[0].strip()
        for requirement in dist_requires
        if _normalize_marker(requirement.partition(";")[2]) == f"extra=='{_EXTRA}'"
    ]


def install_colab(extras: Sequence[str] = ()) -> None:
    """Install the ``[notebook]`` stack into the running environment.

    Intended for the first code cell of the example notebooks. Outside
    Google Colab this is a no-op (local environments follow the Setup
    section), so the cell is safe to run anywhere.

    ``extras`` adds per-notebook requirements: the bats notebook passes
    the legacy ``.xls``/map pieces that sit on top of ``[notebook]``
    (``pandas[excel]``, ``osmnx``, ``xlrd``; the geospatial stack itself
    is part of the extra).
    """
    try:
        import google.colab  # noqa: F401
    except ImportError:
        return
    requirements = notebook_extra_requirements() + list(extras)
    if not requirements:
        # unwrap is importable but has no distribution metadata (run
        # from a source checkout): the bootstrap installed the package
        # itself, so there is nothing further to add.
        return
    subprocess.check_call(
        [sys.executable, "-m", "pip", "install", "-q", *requirements]
    )
