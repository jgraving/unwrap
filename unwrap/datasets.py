"""
Locust dataset access.

The raw datasets are hosted in a separate repository (``jgraving/unwrap-data``)
so that this repo stays small. This module fetches the parts you need and
returns the same DataFrames the notebooks always saw.

Layout inside ``unwrap-data``:

* ``directional/*.parquet``       -- one file per locust from the VR
  directional experiment (Sayin et al. 2025). Columns: ``X, Y, fname``.
* ``density_order/XY_part_*.parquet``
                                 -- raw tracking parts from the
  density/order experiment. Columns: ``X, Y, fname, loss, order,
  density, groups``. The parts are size-capped, not logical units;
  read them all.

Both datasets are already snappy-compressed parquet files (the same
data the original ``.h5`` tables or tracking-export dumps contained,
just in a format git can carry).

Usage
-----
``load_directional()`` returns a single per-animal DataFrame
(columns ``X, Y, fname, animal_id``).

``load_density_order()`` returns the concatenated parts as one
DataFrame with the experiment's columns (``X, Y, fname, loss, order,
density, groups``).

Either way call :func:`ensure_clone` first (or pass
``data_root=`` to skip fetching) so the files are local.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

DEFAULT_REPO = "jgraving/unwrap-data"
DEFAULT_DATA_ROOT = Path(os.path.expanduser("~/.cache/unwrap-data"))


@dataclass
class UnwrapDataRepo:
    """Handle to a local copy of the ``unwrap-data`` repository."""

    root: Path
    shallow: bool = True

    @property
    def directional_dir(self) -> Path:
        return self.root / "directional"

    @property
    def density_order_dir(self) -> Path:
        return self.root / "density_order"

    def _git(self, *args: str, check: bool = True) -> str:
        out = subprocess.run(
            ["git", *args],
            cwd=self.root,
            capture_output=True,
            text=True,
            check=check,
        )
        return out.stdout


def ensure_clone(
    repo: str = DEFAULT_REPO,
    data_root: Path | None = None,
    *,
    shallow: bool = True,
    force: bool = False,
) -> UnwrapDataRepo:
    """Ensure ``data_root/unwrap-data`` is a local clone of ``repo``.

    Runs ``git clone --depth 1`` on first call, ``git pull``
    afterwards (shallow clones pull with ``--depth 1``).

    Parameters
    ----------
    repo:
        GitHub slug, e.g. ``jgraving/unwrap-data``.
    data_root:
        Parent directory the clone lives in.
        Defaults to ``~/.cache/unwrap-data``.
    shallow:
        Keep the clone to one commit (much faster, much smaller).
    force:
        Delete any existing clone before fetching.
    """
    data_root = Path(data_root) if data_root is not None else DEFAULT_DATA_ROOT
    root = data_root / "unwrap-data"

    if force and root.exists():
        shutil.rmtree(root)
    if not root.exists():
        data_root.mkdir(parents=True, exist_ok=True)
        shallow_flag = "--depth 1" if shallow else ""
        subprocess.run(
            ["git", "clone", shallow_flag, f"https://github.com/{repo}.git", root],
            check=True,
        )
        return UnwrapDataRepo(root=root, shallow=shallow)
    if (root / ".git").is_dir():
        if shallow:
            root._git("fetch", "--depth", "1", "origin", "main", check=False)
            root._git("checkout", "main", check=False)
            root._git("reset", "--hard", "origin/main", check=False)
        else:
            root._git("pull", "--ff-only", check=False)
        return UnwrapDataRepo(root=root, shallow=shallow)
    raise RuntimeError(
        f"{root} exists but is not a git checkout. "
        "Pass force=True to wipe and reclone."
    )


def _require_dir(p: Path) -> None:
    if not p.is_dir():
        raise FileNotFoundError(
            f"expected {p} to exist. Run ensure_clone() first, or pass "
            "data_root= to a directory containing the dataset subdirectory."
        )


def _resolve_dir(kind: str, data_root: Path | None) -> Path:
    """Resolve the parquet parts directory for ``kind`` ("directional" or
    "density_order"). ``data_root`` may be None (clone/fetch the default
    repo), an unwrap-data repo root (containing the dataset subdir), or the
    dataset parts directory itself (a local checkout)."""
    subdir = "directional" if kind == "directional" else "density_order"
    if data_root is None:
        repo = ensure_clone()
        d = repo.directional_dir if kind == "directional" else repo.density_order_dir
        _require_dir(d)
        return d
    root = Path(data_root)
    if (root / subdir).is_dir():
        d = root / subdir
    elif list(root.glob("*.parquet")) or list(root.glob("XY_part_*.parquet")):
        d = root
    else:
        d = root / subdir
    _require_dir(d)
    return d


def load_directional(
    data_root: Path | None = None,
    *,
    subset: list[str] | None = None,
) -> pd.DataFrame:
    """Load all per-animal directional trajectories into one DataFrame.

    Adds a 0-based integer ``animal_id`` column, in sorted filename
    order, matching the mapping used in the original slm notebook.
    Also adds a ``source_table`` column that records which original
    source file the animal came from (``XY_raw_50cm_`` or ``XY_raw_8cm_``).

    Parameters
    ----------
    data_root:
        Directory that contains ``directional/`` (or the ``directional/``
        directory itself; pass an actual dir holding ``*.parquet`` to
        use it directly). Defaults to the cloned ``unwrap-data`` repo
        (``~/.cache/unwrap-data/unwrap-data``).
    subset:
        If given, only load animals whose ``fname`` is in this list.
        Each file holds one animal; the match is made on the trailing
        ``YYYYMMDD_HHMM`` of the filename so we never open a file
        just to peek at its columns.
    """
    base = _resolve_dir("directional", data_root)
    files = sorted(base.glob("*.parquet"))
    if subset is not None:
        # fname is a trailing "YYYYMMDD_HHMM"; the filename is
        # "<prefix><fname>.parquet". Match on the tail to avoid reading
        # each (one-animal) file just to peek at ``fname``.
        import re

        keep = set(subset)
        tail = re.compile(r"(\d{8}_\d{4})")
        files = [f for f in files if tail.search(f.name).group(1) in keep]
        if not files:
            raise LookupError(f"no directional files matched subset={subset!r}")
    frames = []
    for f in files:
        part = pd.read_parquet(f)
        # Recover which original source table this animal came from
        # (XY_raw_50cm_.h5 vs XY_raw_8cm_.h5) from the filename prefix,
        # so downstream results stay comparable to the original nb.
        part["source_table"] = (
            "XY_raw_50cm_" if f.name.startswith("XY_raw_50cm_") else "XY_raw_8cm_"
        )
        frames.append(part)
    df = pd.concat(frames, ignore_index=True)
    df["animal_id"] = pd.factorize(df["fname"], sort=False)[0].astype("int32")
    return df


def load_density_order(
    data_root: Path | None = None,
    *,
    columns_only: list[str] | None = None,
) -> pd.DataFrame:
    """Load all density/order tracking parts into one DataFrame.

    Invalid / partial files (corrupted footer) are skipped with a printed
    warning, matching the original notebook's behaviour.

    Parameters
    ----------
    data_root:
        Directory that contains ``density_order/`` (or the parts
        directory itself). Defaults to the cloned ``unwrap-data`` repo
        (``~/.cache/unwrap-data/unwrap-data``).
    columns_only:
        If given, only these columns are read (memory saver).
    """
    base = _resolve_dir("density_order", data_root)
    files = sorted(base.glob("*.parquet"))
    import pyarrow.parquet as pq

    good, bad = [], []
    for f in files:
        try:
            pq.ParquetFile(f)
            good.append(f)
        except Exception as e:  # noqa: BLE001
            bad.append((f.name, str(e)))
    if bad:
        print("Skipping invalid/corrupted density_order files:")
        for name, why in bad:
            print(f"  - {name} -> {why}")
    if not good:
        raise RuntimeError(f"no valid parquet parts under {base}")
    df = pd.concat(
        [pd.read_parquet(f, columns=columns_only) for f in good],
        ignore_index=True,
    )
    return df
