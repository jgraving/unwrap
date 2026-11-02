"""Colab bootstrap helpers for the example notebooks.

The "Open in Colab" button on a notebook page starts the notebook
straight from GitHub, so the Colab runtime has no repository checkout
and ``unwrap`` is not installed. The first code cell of each example
notebook bootstraps the package from its git URL and then calls
:func:`install_colab` to finish the environment.
"""

import subprocess
import sys
import urllib.request
from collections.abc import Sequence
from importlib.metadata import PackageNotFoundError, requires
from pathlib import Path

__all__ = ["install_colab", "notebook_extra_requirements", "scripts_dir"]

_EXTRA = "notebook"
_REPO = "jgraving/unwrap"
_DEFAULT_REF = "main"


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


def scripts_dir() -> Path:
    """Where the example notebooks' sibling scripts belong.

    Every example notebook inserts
    ``os.path.dirname(os.path.dirname(unwrap.__file__)) / "notebooks"``
    into ``sys.path`` and imports its ``<study>_preprocess`` /
    ``<study>_posthoc`` / ``<study>_plots`` scripts from there: in a
    checkout that is the repository's ``notebooks/`` folder, in an
    installed environment the matching folder next to the package.
    """
    return Path(__file__).resolve().parent.parent / "notebooks"


def _fetch_script(name: str, ref: str, directory: Path) -> None:
    url = f"https://raw.githubusercontent.com/{_REPO}/{ref}/notebooks/{name}"
    with urllib.request.urlopen(url) as response:  # noqa: S310
        (directory / name).write_bytes(response.read())


def install_colab(extras: Sequence[str] = (),
                  scripts: Sequence[str] = (),
                  ref: str = _DEFAULT_REF) -> None:
    """Install the ``[notebook]`` stack into the running environment.

    Intended for the first code cell of the example notebooks. Outside
    Google Colab this is a no-op (local environments follow the Setup
    section), so the cell is safe to run anywhere.

    ``extras`` adds per-notebook requirements: the bats notebook passes
    the legacy ``.xls``/map pieces that sit on top of ``[notebook]``
    (``pandas[excel]``, ``osmnx``, ``xlrd``; the geospatial stack itself
    is part of the extra).

    ``scripts`` fetches the notebook's sibling scripts (same names as in
    the repository's ``notebooks/`` folder) from ``ref`` — the same ref
    the cell's pip bootstrap installed the package from — into
    :func:`scripts_dir`, the directory the notebook's own ``sys.path``
    line targets. Colab only received the ``.ipynb`` itself, so without
    this the ``import beetles_preprocess``-style cells would fail.
    """
    try:
        import google.colab  # noqa: F401
    except ImportError:
        return
    requirements = notebook_extra_requirements() + list(extras)
    if requirements:
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", "-q", *requirements]
        )
    if scripts:
        directory = scripts_dir()
        directory.mkdir(parents=True, exist_ok=True)
        for name in scripts:
            _fetch_script(name, ref, directory)
