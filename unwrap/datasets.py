"""
Dataset access.

The raw datasets are hosted in a separate repository (``jgraving/unwrap-data``)
so that this repo stays small. This module fetches the parts you need and
returns the same DataFrames the notebooks always saw. It serves four
experiments: the two locust studies (:func:`load_directional`,
:func:`load_density_order`), the beetle security-light trials
(:func:`load_beetle_security_light`), and the zebra escape bearings
(:func:`load_zebra_bearings`, plus :func:`load_zebra_bearings_raw` for the
full per-frame raw source).

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

# dataset kind -> subdirectory of the unwrap-data repo
_SUBDIRS = {
    "directional": "directional",
    "density_order": "density_order",
    "beetle_security_light": "beetle_security_light",
    "zebra_bearings": "zebra_bearings",
}


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
    """Resolve the parquet parts directory for ``kind``. ``data_root`` may be
    None (clone/fetch the default repo), an unwrap-data repo root (containing
    the dataset subdir), or the dataset directory itself (a local checkout)."""
    subdir = _SUBDIRS[kind]
    if data_root is None:
        d = ensure_clone().root / subdir
        _require_dir(d)
        return d
    root = Path(data_root)
    d = root / subdir if (root / subdir).is_dir() else root
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


def _resolve_file(kind: str, filename: str, data_root: Path | None) -> Path:
    """Resolve a single dataset file. ``data_root`` may be None (cloned
    repo), a repo root, or the dataset directory itself."""
    base = _resolve_dir(kind, data_root)
    f = (base / filename).expanduser()
    if not f.is_file():
        raise FileNotFoundError(f"expected {f} to exist.")
    return f


def load_beetle_security_light(
    data_root: Path | None = None,
    *,
    sep: str = "\t",
) -> pd.DataFrame:
    """Load the 2019 security-light beetle experiment.

    Returns the raw per-trial table (``Experiment``, ``Beetle``, ``Trial``,
    ``Heading`` in degrees). The ``beetles.ipynb`` notebook applies its own
    Light/HeadingRadians transforms on top, so this returns the data as-is.
    """
    f = _resolve_file(
        "beetle_security_light", "SecurityLightBergsig.txt", data_root
    )
    return pd.read_csv(f, sep=sep)


def load_zebra_bearings(
    data_root: Path | None = None,
    *,
    columns_only: list[str] | None = None,
) -> pd.DataFrame:
    """Load zebra escape-bearings block data.

    Returns columns ``Track, Observation, Condition, HeadingRadians,
    PostStimulus``. ``Condition`` is the ordered categorical
    {pre-approach, post-approach} exactly as the ``zebra_escape.ipynb``
    notebook expects.
    """
    f = _resolve_file("zebra_bearings", "zebra_bearings_blocks.parquet", data_root)
    df = pd.read_parquet(f, columns=columns_only)
    if "Condition" in df.columns:
        df["Condition"] = pd.Categorical(
            df["Condition"],
            categories=["pre-approach", "post-approach"],
            ordered=True,
        )
    return df


def load_zebra_bearings_raw(
    data_root: Path | None = None,
    *,
    sort: bool = True,
) -> pd.DataFrame:
    """Load the raw per-frame zebra escape-bearing observations.

    The raw frames are split into one small parquet file per tracked animal
    (``zebra_bearings/raw/<track>.parquet``) so each stays within git's
    per-file limit. Every row carries a ``row_index`` naming its row in the
    original ``zebra_bearings.csv``; with ``sort=True`` the reassembled frame
    is value-for-value the original CSV (verified lossless).

    Returns ``observation``, ``track_number``, ``frame``,
    ``compass_bearing``, ``pre_scare``, ``row_index`` -- 4 zebras, 59 tracks,
    3,897,536 frames. This is the *raw* source; :func:`load_zebra_bearings`
    (the preprocessed 300-frame block view the notebook fits) is a few
    hundred times smaller and is what ``zebra_escape.ipynb`` actually uses.
    Load this only when you need to reproduce that preprocessing.
    """
    base = _resolve_dir("zebra_bearings", data_root)
    raw = base / "raw"
    _require_dir(raw)
    files = sorted(raw.glob("*.parquet"))
    if not files:
        raise FileNotFoundError(f"no raw zebra files under {raw}")
    df = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
    if sort:
        df = df.sort_values("row_index", kind="mergesort").reset_index(drop=True)
    return df
