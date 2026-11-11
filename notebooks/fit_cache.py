"""Development fit cache for the Python notebooks.

When ``UNWRAP_PY_FIT_CACHE=1`` a notebook's PyMC / Bambi fits are reused
across local runs instead of resampled: before fitting, a cell asks
:func:`unwrap_load_fit_cache` for a saved snapshot; when none is valid it
fits as usual and stores the result with :func:`unwrap_save_fit_cache`.

This mirrors the R convention in ``R/unwrap_functions.R``
(``unwrap_load_fit_cache`` / ``unwrap_save_fit_cache``). A snapshot is
reused only when the name, seed, sampling budget (tune/draws), and a
data checksum all match, so a stale cache can never sneak into a run.
``name`` is caller-provided and must identify the MODEL as well as its
budget: a model change (priors/formula/family) under an unchanged name
would silently reuse the old snapshot, so bump the name (e.g. a ``_v2``
suffix) whenever the model changes.

With ``UNWRAP_PY_FIT_CACHE`` unset (the default, and how committed
notebooks are validated) these helpers are silent no-ops and nothing is
read from or written to disk.

Storage
-------
Fits are ``arviz.InferenceData`` objects, which in arviz 1.x are
``xarray.DataTree``. They are saved with ``InferenceData.to_netcdf`` and
loaded with ``arviz.from_netcdf`` (NetCDF4), the xarray-native, portable
format - no pickling. The cache metadata is carried in the file's global
attributes so a single ``.nc`` file holds both the draws and the validity
check.
"""

from __future__ import annotations

import hashlib
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Enable / path
# ---------------------------------------------------------------------------

_ENV_FLAG = "UNWRAP_PY_FIT_CACHE"
_ENV_DIR = "UNWRAP_PY_FIT_CACHE_DIR"


def unwrap_fit_cache_enabled() -> bool:
    """True only when ``UNWRAP_PY_FIT_CACHE=1`` is set."""
    return os.environ.get(_ENV_FLAG) == "1"


def _unwrap_data_root() -> Path | None:
    """Best-effort path to the local ``unwrap-data`` checkout, or None."""
    try:
        from unwrap.datasets import DEFAULT_DATA_ROOT
    except Exception:  # pragma: no cover - unwrap not importable
        return None
    root = DEFAULT_DATA_ROOT / "unwrap-data"
    return root if root.is_dir() else None


def unwrap_fit_cache_path(name: str) -> Path | None:
    """Path of the snapshot for ``name`` (None when the cache is disabled).

    The directory is ``UNWRAP_PY_FIT_CACHE_DIR`` when set; otherwise the
    ``fit_cache/`` directory of the local ``unwrap-data`` checkout (mirroring
    where the R fit snapshots live); otherwise ``~/.unwrap_fit_cache``.
    """
    if not unwrap_fit_cache_enabled():
        return None
    override = os.environ.get(_ENV_DIR, "").strip()
    if override:
        directory = Path(override)
    else:
        root = _unwrap_data_root()
        directory = (root / "fit_cache") if root is not None else (
            Path.home() / ".unwrap_fit_cache")
    return directory / f"{name}.nc"


# ---------------------------------------------------------------------------
# Data checksum
# ---------------------------------------------------------------------------

def data_checksum(response, *, group_arrays=()) -> str:
    """A short, stable checksum identifying the data a fit was built from.

    ``response`` is the model's heading/response values (a 1-D array);
    ``group_arrays`` is an optional iterable of ``(name, array)`` for any
    grouping columns (e.g. animal / tag id). The checksum combines the row
    count, a hash of the response, and a hash of each grouping array, so a
    change to the data - different values, a different number of rows, or a
    re-ordered / re-labelled grouping - invalidates the cache. It is
    deliberately a coarse fingerprint (not a byte-for-byte hash of the whole
    frame): cheap to compute and stable across pandas / numpy versions.
    """
    response = np.asarray(response)
    parts = [f"n:{len(response)}", f"resp:{response.shape}:{_hash(response)}"]
    for name, values in group_arrays:
        values = np.asarray(values)
        parts.append(f"{name}:{values.shape}:{_hash(values)}")
    payload = "|".join(parts)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]


def _hash(values) -> str:
    return hashlib.sha256(np.asarray(values).tobytes()).hexdigest()


# ---------------------------------------------------------------------------
# Metadata (carried in the NetCDF global attributes)
# ---------------------------------------------------------------------------

def _set_meta(idata, *, name, seed, tune, draws, checksum) -> None:
    idata.attrs["unwrap_cache_name"] = str(name)
    idata.attrs["unwrap_cache_seed"] = int(seed)
    idata.attrs["unwrap_cache_tune"] = int(tune)
    idata.attrs["unwrap_cache_draws"] = int(draws)
    idata.attrs["unwrap_cache_checksum"] = str(checksum)
    idata.attrs["unwrap_cache_saved_utc"] = datetime.now(
        timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _meta_matches(attrs, *, name, seed, tune, draws, checksum) -> bool:
    def _num(key):
        value = attrs.get(f"unwrap_cache_{key}")
        try:
            return int(value)
        except (TypeError, ValueError):
            return None

    return (
        attrs.get("unwrap_cache_name") == str(name)
        and _num("seed") == int(seed)
        and _num("tune") == int(tune)
        and _num("draws") == int(draws)
        and attrs.get("unwrap_cache_checksum") == str(checksum)
    )


# ---------------------------------------------------------------------------
# Load / save
# ---------------------------------------------------------------------------

def unwrap_load_fit_cache(name, seed, tune, draws, checksum):
    """Return a valid cached ``InferenceData`` for ``name``, or None.

    A snapshot is returned only when it exists and its stored name, seed,
    tune, draws, and data checksum all match the arguments; otherwise None
    is returned (and the caller refits). Silent no-op (returns None) when
    the cache is disabled.
    """
    path = unwrap_fit_cache_path(name)
    if path is None or not path.exists():
        return None
    try:
        import arviz as az

        idata = az.from_netcdf(path)
    except Exception:
        message("Unwrap fit cache: '%s' is unreadable; refitting." % name)
        return None
    if not _meta_matches(idata.attrs, name=name, seed=seed, tune=tune,
                         draws=draws, checksum=checksum):
        message("Unwrap fit cache: '%s' is stale; refitting." % name)
        return None
    message("Unwrap fit cache: reusing '%s' (seed %s, tune %s, draws %s)."
            % (name, seed, tune, draws))
    return idata


def unwrap_save_fit_cache(name, idata, seed, tune, draws, checksum):
    """Store ``idata`` under ``name`` (silent no-op when disabled)."""
    path = unwrap_fit_cache_path(name)
    if path is None:
        return None
    path.parent.mkdir(parents=True, exist_ok=True)
    _set_meta(idata, name=name, seed=seed, tune=tune, draws=draws,
              checksum=checksum)
    # Write to a temp file in the same directory, then rename, so a
    # crashed bake never leaves a half-written snapshot behind.
    fd, tmp = tempfile.mkstemp(suffix=".tmp", dir=str(path.parent))
    os.close(fd)
    try:
        idata.to_netcdf(tmp)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
    return None


def message(text: str) -> None:
    """Print a cache notice (kept separate so callers stay quiet by default)."""
    print(text)
