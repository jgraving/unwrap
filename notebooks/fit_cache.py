"""Load-or-save fit cache for the Python notebooks (issue #130, bats part).

Environment control
-------------------
``UNWRAP_PY_FIT_CACHE``
    When unset or not a truthy value (``1``/``true``/``yes``/``on``), this
    module is a complete no-op: :func:`load_or_fit` just calls ``fit_fn``
    and performs no file I/O. This is the default, so every existing
    notebook run behaves exactly as before.

``UNWRAP_PY_FIT_CACHE_DIR``
    Directory holding the snapshots (default: ``<repo>/fit_cache``,
    mirroring the R fit caches pushed to the unwrap-data repo).

Snapshots
---------
A fit named ``name`` is stored as two files:

* ``<name>.nc``         the xarray ``DataTree`` in the netCDF format the
                        array package natively saves and loads
                        (``DataTree.to_netcdf`` / ``xarray.open_datatree``,
                        a netCDF-4 file): the whole object, including
                        ``log_likelihood`` when it was requested; in slim
                        snapshots the variables listed under
                        ``extra["stripped"]`` are absent and recomputed
                        from the stored draws on load
* ``<name>.meta.json``  sha256 of the netCDF file plus validation
                       metadata: the posterior variable list and shapes,
                       and the caller-supplied fit parameters (seed, tune,
                       draws, model id, library versions).

A snapshot is reused only when the sha256 matches and every caller-supplied
metadata value equals the stored one; any mismatch refits and overwrites
the snapshot (loudly, on the console).
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import xarray

__all__ = ["CacheMismatch", "enabled", "cache_dir", "load_or_fit"]

ENV_ENABLE = "UNWRAP_PY_FIT_CACHE"
ENV_DIR = "UNWRAP_PY_FIT_CACHE_DIR"

#: <repo>/fit_cache (this file lives in notebooks/)
DEFAULT_DIR = Path(__file__).resolve().parent.parent / "fit_cache"

_TRUTHY = {"1", "true", "yes", "on"}


class CacheMismatch(RuntimeError):
    """A snapshot exists but fails sha256 or metadata validation."""


def enabled() -> bool:
    """True when the caller opted into the cache via the environment."""
    return str(os.environ.get(ENV_ENABLE, "")).strip().lower() in _TRUTHY


def cache_dir() -> Path:
    override = os.environ.get(ENV_DIR, "").strip()
    return Path(override).expanduser() if override else DEFAULT_DIR


def _snapshot_files(name: str):
    root = cache_dir() / name
    return root.with_suffix(".nc"), root.with_suffix(".meta.json")


def _digest(blob: bytes) -> str:
    return hashlib.sha256(blob).hexdigest()


def _describe(idata) -> dict:
    post = idata.posterior
    return {
        "n_posterior_vars": len(post.data_vars),
        "posterior_vars": sorted(map(str, post.data_vars)),
        "posterior_shapes": {
            str(k): list(v.shape) for k, v in post.items()
        },
    }


def load_snapshot(name: str, *, expect: dict | None = None):
    """Load the cached fit ``name``, validating sha256 and metadata.

    Raises ``FileNotFoundError`` (nothing saved) or ``CacheMismatch``
    (corrupted or stale snapshot).
    """
    snap, meta_path = _snapshot_files(name)
    if not (snap.exists() and meta_path.exists()):
        raise FileNotFoundError(
            f"no saved fit snapshot for '{name}' in {cache_dir()}"
        )
    info = json.loads(meta_path.read_text())
    blob = snap.read_bytes()
    if info.get("sha256") != _digest(blob):
        raise CacheMismatch(f"snapshot '{name}': netCDF file sha256 does not match the metadata")
    if expect:
        stored = info.get("extra", {})
        bad = {k: (stored.get(k), v) for k, v in expect.items() if stored.get(k) != v}
        if bad:
            raise CacheMismatch(f"snapshot '{name}': fit metadata mismatch {bad}")
    obj = xarray.open_datatree(snap).load()
    if _describe(obj) != {k: info[k] for k in
                          ("n_posterior_vars", "posterior_vars", "posterior_shapes")}:
        raise CacheMismatch(f"snapshot '{name}': posterior layout does not match the metadata")
    return obj


def save_snapshot(name: str, idata, *, meta: dict | None = None) -> None:
    snap, meta_path = _snapshot_files(name)
    snap.parent.mkdir(parents=True, exist_ok=True)
    idata.to_netcdf(snap)
    blob = snap.read_bytes()
    info = _describe(idata)
    info["sha256"] = _digest(blob)
    info["extra"] = dict(meta or {})
    meta_path.write_text(json.dumps(info, indent=2, sort_keys=True) + "\n")


def slim_snapshot(idata, *, drop_groups=(), drop_posterior=()):
    """Return a copy of ``idata`` without the named large groups and/or
    posterior variables, for small fit snapshots (issue #130): the bats
    vMU snapshot strips the per-observation deterministics (``mu``, ``kappa``,
    ``p_goal``) and the ``log_likelihood`` group, which are cheaply
    recomputed from the stored latent draws when a snapshot is loaded
    (notebook cell 53 / posthoc.restore_per_obs_deterministics).

    Only strips what is actually present; other groups and variables pass
    through untouched.
    """
    drop_groups = {g for g in drop_groups if g in idata.children}
    drop_posterior = set(drop_posterior)

    dt = idata
    if drop_groups:
        dt = dt.drop_nodes(sorted(drop_groups))
    if drop_posterior:
        present = [v for v in dt.posterior.data_vars if v in drop_posterior]
        if present:
            post_ds = dt.posterior.to_dataset().drop_vars(present)
            dt["posterior"] = xarray.DataTree(name="posterior", dataset=post_ds)
    return dt


def load_or_fit(name: str, fit_fn, *, meta: dict | None = None,
                pick=None, label: str | None = None):
    """Reuse a saved fit, or run ``fit_fn()`` and save its snapshot.

    With the cache disabled (the default) this returns ``fit_fn()``
    unchanged and touches no files. With it enabled:

    * a valid snapshot is loaded and returned (``fit_fn`` is not called);
    * a missing or mismatched snapshot refits via ``fit_fn`` and overwrites
      the snapshot, printing exactly what happened.

    ``meta`` (fit parameters such as seed/tune/draws/model id) is validated
    on load and stored on save. ``pick`` maps ``fit_fn``'s return value to
    the ``az.InferenceData`` that gets snapshotted (identity by default);
    the full ``fit_fn`` return value is returned to the caller on refits.
    """
    if not enabled():
        return fit_fn()

    lbl = label or name
    pick = pick or (lambda r: r)
    try:
        obj = load_snapshot(name, expect=meta)
        print(f"FIT-CACHE hit: {lbl} <- {(cache_dir() / name).with_suffix('.nc')}")
        return obj
    except FileNotFoundError:
        print(f"FIT-CACHE miss: {lbl}; fitting now and saving "
              f"{(cache_dir() / name).with_suffix('.nc')}")
    except CacheMismatch as exc:
        print(f"FIT-CACHE mismatch ({exc}); refitting and overwriting the snapshot")

    res = fit_fn()
    save_snapshot(name, pick(res), meta=meta)
    print(f"FIT-CACHE saved: {lbl}")
    return res
