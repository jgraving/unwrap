"""
Dataset access.

The raw datasets are hosted in a separate repository (``jgraving/unwrap-data``)
so that this repo stays small. This module fetches the parts you need and
returns the same DataFrames the notebooks always saw. It serves five
experiments: the two locust studies (:func:`load_directional`,
:func:`load_density_order`), the beetle security-light trials
(:func:`load_beetle_security_light`), the zebra escape bearings
(:func:`load_zebra_bearings`, plus :func:`load_zebra_bearings_raw` for the
full per-frame raw source and :func:`load_zebra_postures` /
:func:`load_zebra_mosaic` for the escape-map figure's posture positions
and aerial-imagery backdrop), and the bat homing-tracking tracks
(:func:`load_bats_tracks`, :func:`load_bats_info`,
:func:`load_bats_lake_outline`).

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
``data_root=`` to skip fetching) so the files are local. A clone is
made once and then reused as-is: subsequent loader calls cost only
local file reads, with no re-download and no network round-trip. Pass
``update=True`` to :func:`ensure_clone` to pull newly published data.
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
    "goldshtein_bats": "goldshtein_bats",
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


def _git_stream(root: Path, *args: str, check: bool = True) -> None:
    """Run a git command inside a local checkout, with stderr inherited so
    ``--progress`` reporting (download percentage, speed) is visible in
    notebooks and terminals instead of being swallowed."""
    subprocess.run(["git", *args], cwd=root, check=check)


def ensure_clone(
    repo: str = DEFAULT_REPO,
    data_root: Path | None = None,
    *,
    shallow: bool = True,
    force: bool = False,
    update: bool = False,
) -> UnwrapDataRepo:
    """Ensure ``data_root/unwrap-data`` is a local clone of ``repo``.

    The clone is made **once** (``git clone --depth 1``, with git's progress
    meter streamed to the notebook) and is then reused as-is: subsequent
    ``load_*()`` calls resolve to local files with no network round-trip and
    no re-download. Use ``update=True`` to pull in newly published data
    (fetch + hard reset, also shown with a progress meter); without it an
    existing checkout is trusted in place.

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
    update:
        Re-fetch the latest ``main`` and hard-reset to it (progress meter
        shown). Default ``False`` means "use what's already here".
    """
    data_root = Path(data_root) if data_root is not None else DEFAULT_DATA_ROOT
    root = data_root / "unwrap-data"

    if force and root.exists():
        shutil.rmtree(root)
    if not root.exists():
        data_root.mkdir(parents=True, exist_ok=True)
        shallow_flags = ["--depth", "1"] if shallow else []
        subprocess.run(
            ["git", "clone", *shallow_flags, "--progress",
             f"https://github.com/{repo}.git", root],
            check=True,
        )
        return UnwrapDataRepo(root=root, shallow=shallow)
    if not (root / ".git").is_dir():
        raise RuntimeError(
            f"{root} exists but is not a git checkout. "
            "Pass force=True to wipe and reclone."
        )
    if update:
        if shallow:
            _git_stream(root, "fetch", "--depth", "1", "--progress", "origin", "main")
            _git_stream(root, "checkout", "main")
            _git_stream(root, "reset", "--hard", "origin/main")
        else:
            _git_stream(root, "pull", "--progress", "--ff-only")
    return UnwrapDataRepo(root=root, shallow=shallow)


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
    tracks: list[str] | None = None,
) -> pd.DataFrame:
    """Load the raw per-frame zebra escape-bearing observations.

    The raw frames are split into one small parquet file per tracked animal
    (``zebra_bearings/raw/<track>.parquet``) so each stays within git's
    per-file limit. Every row carries a ``row_index`` naming its row in the
    original ``zebra_bearings.csv``; with ``sort=True`` the reassembled frame
    follows the original CSV order. Pass full track IDs such as
    ``observation088-5`` in ``tracks`` to load only those files.

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
    files = (sorted(raw.glob("*.parquet")) if tracks is None
             else [raw / f"{track}.parquet" for track in tracks])
    if not files:
        raise FileNotFoundError(f"no raw zebra files under {raw}")
    df = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
    if sort:
        df = df.sort_values("row_index", kind="mergesort").reset_index(drop=True)
    return df


def load_zebra_postures(
    data_root: Path | None = None,
    *,
    observation: str = "observation088",
    tracks: list[int] | None = None,
) -> pd.DataFrame:
    """Load zebra posture keypoint positions (the input to the escape-map figure).

    The original ``<observation>-postures.npy`` array
    ``(tracks, frames, keypoints, 2)`` in UTM 37N meters was split per
    tracked animal into wide float64 parquets
    (``zebra_bearings/postures/<observation>-<track>.parquet``) so each
    file stays within git's per-file limit; the split was verified
    element-for-element against the npy (``tools/prepare_zebra_map_data.py``).

    Returns ``frame`` (int, original frame index), ``track_number`` (str)
    and one column per keypoint coordinate (``kp0_x, kp0_y, ..., kp9_y``,
    float64, NaN where the keypoint was not detected). Keypoint 0 is the
    animal's centre (its trajectory), and the per-frame body bearings from
    keypoints 8 (rump) and 3 (shoulder) are already computed in
    :func:`load_zebra_bearings_raw`.

    ``tracks`` selects track numbers (the part after the dash in
    ``track_number``); default loads every track of the observation.
    """
    base = _resolve_dir("zebra_bearings", data_root)
    pos = base / "postures"
    _require_dir(pos)
    if tracks is None:
        files = sorted(pos.glob(f"{observation}-*.parquet"))
    else:
        files = [pos / f"{observation}-{int(t)}.parquet" for t in tracks]
    files = [f for f in files if f.exists()]
    if not files:
        raise FileNotFoundError(f"no posture files under {pos} for {observation}")
    return pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)


def load_zebra_mosaic(
    data_root: Path | None = None,
    *,
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
    max_side_px: int | None = 4000,
) -> tuple:
    """Load the aerial orthomosaic behind the escape-map figure.

    The published map panels crop the original 38645x44856 px orthomosaic
    (RGBA uint8, 0.02707 m/px, UTM zone 37N) to a single native-resolution
    window, stored losslessly as PNG row-strips plus a JSON manifest
    (``zebra_bearings/mosaic/``; see ``tools/prepare_zebra_map_data.py``).

    Returns ``(img, extent)`` where ``img`` is a ``(H, W, 4)`` uint8 RGBA
    array and ``extent`` is ``(left, right, bottom, top)`` in UTM 37N
    meters, ready for ``ax.imshow(img, extent=extent, origin="upper")``.

    ``xlim``/``ylim`` (meters) return only the requested window; tiles and
    columns outside it are never read. ``max_side_px`` downsamples the
    returned image (nearest-neighbor) so it stays a sane size when
    embedded in the saved PDF; pass ``None`` for the full 0.027 m/px crop.
    """
    import json

    import numpy as np
    from PIL import Image

    base = _resolve_dir("zebra_bearings", data_root)
    mos = base / "mosaic"
    _require_dir(mos)
    manifest = json.loads((mos / "manifest.json").read_text())
    left, right = manifest["extent"]["left"], manifest["extent"]["right"]
    top, bottom = manifest["extent"]["top"], manifest["extent"]["bottom"]
    ps = manifest["pixel_size_m"]
    n_rows, n_cols = manifest["n_rows"], manifest["n_cols"]

    # global pixel window (row 0 = top = y == top; cols increase with x)
    c0, c1 = 0, n_cols
    if xlim is not None:
        c0 = max(0, int(np.floor((min(xlim) - left) / ps)))
        c1 = min(n_cols, int(np.ceil((max(xlim) - left) / ps)))
    r0, r1 = 0, n_rows
    if ylim is not None:
        r0 = max(0, int(np.floor((top - max(ylim)) / ps)))
        r1 = min(n_rows, int(np.ceil((top - min(ylim)) / ps)))
    if c0 >= c1 or r0 >= r1:
        raise ValueError("requested window lies outside the mosaic crop")

    strips = []
    for tile in manifest["tiles"]:
        t0, t1 = tile["row0"], tile["row0"] + tile["nrows"]
        if t1 <= r0 or t0 >= r1:
            continue
        img = Image.open(mos / tile["file"]).convert("RGBA")
        a = np.asarray(img)
        del img
        if t0 < r0 or t0 + a.shape[0] > r1:  # keep only intersecting rows
            lo, hi = max(r0 - t0, 0), min(r1 - t0, a.shape[0])
            a = a[lo:hi]
        strips.append(a[:, c0:c1])
    out = np.concatenate(strips, axis=0)

    extent = (left + c0 * ps, left + c1 * ps, top - r1 * ps, top - r0 * ps)
    if max_side_px is not None and max(out.shape[:2]) > max_side_px:
        scale = max_side_px / max(out.shape[:2])
        w_new = max(1, int(round(out.shape[1] * scale)))
        h_new = max(1, int(round(out.shape[0] * scale)))
        img = Image.fromarray(out, mode="RGBA")
        out = np.asarray(img.resize((w_new, h_new), Image.NEAREST))
    return out, extent


def _read_bats_excel(filename: str, data_root: Path | None, *,
                     columns_only: list[str] | None) -> pd.DataFrame:
    """Read a bats ``.xls`` (legacy BIFF) from the bats dataset subdir.

    These are the original tracking-export/workbook files, kept verbatim in
    ``unwrap-data``. They are read with the ``calamine`` engine; declare
    ``python-calamine`` (pandas >= 3 uses it by default, and it reads the
    legacy ``.xls`` format that ``xlrd``/``openpyxl`` do not).
    """
    f = _resolve_file("goldshtein_bats", filename, data_root)
    df = pd.read_excel(f, engine="calamine")
    if columns_only is not None:
        df = df[columns_only]
    return df


def load_bats_tracks(
    data_root: Path | None = None,
    *,
    columns_only: list[str] | None = None,
) -> pd.DataFrame:
    """Load the raw post-localization bat homing tracks (Goldshtein et al.).

    One row per tracking fix; each track is tagged with its destination
    (see the ``tag`` id and ``distDest`` / ``flightAngDiffDest`` columns).
    28,350 rows x 17 cols: ``tag, time, x, y, lon, lat, alt, si, distance,
    speed, flightAngle, turningAngle, distDeparture, distDest,
    flightAngDiffDest, afterLocalization, illumination``. This is the input
    table ``goldshtein_bats.ipynb`` starts from; the notebook does its own
    coordinate transform to web-mercator on top.
    """
    return _read_bats_excel("dataAll_trackDest.xls", data_root,
                            columns_only=columns_only)


def load_bats_info(
    data_root: Path | None = None,
    *,
    columns_only: list[str] | None = None,
) -> pd.DataFrame:
    """Load the per-bat metadata table for the homing experiments.

    67 bats x 56 cols (``tag, year, month, treatment, translocationPoint,
    treatmentName, flyHome, flyAgamon, ...``). Join to
    :func:`load_bats_tracks` on ``tag``.
    """
    return _read_bats_excel("all_bats_info.xls", data_root,
                            columns_only=columns_only)


def load_bats_lake_outline(
    data_root: Path | None = None,
) -> "object":
    """Load the Agamon (Hula) lake shoreline as a GeoDataFrame.

    The polygon used as the natural-homing destination in the harmonic-field
    target-bearing construction. Requires ``geopandas`` (imported lazily so
    the rest of this module works without it). Returns a
    ``geopandas.GeoDataFrame`` with a single polygon (EPSG:4326).
    """
    import geopandas as gpd

    f = _resolve_file("goldshtein_bats", "agamon_hula_lake_outline.geojson",
                      data_root)
    return gpd.read_file(f)


def bats_lake_outline_path(
    data_root: Path | None = None,
) -> Path:
    """Path to the Agamon (Hula) lake outline ``.geojson`` in ``unwrap-data``.

    For cells that ``gpd.read_file(LAKE_GEOJSON)`` a local ``LAKE_GEOJSON``
    path variable rather than taking a returned GeoDataFrame. Resolves the
    file (cloning ``unwrap-data`` on first use) so the path is valid on this
    machine.
    """
    return _resolve_file(
        "goldshtein_bats", "agamon_hula_lake_outline.geojson", data_root
    )
