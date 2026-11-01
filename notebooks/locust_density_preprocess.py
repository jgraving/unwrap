"""Preprocessing for notebooks/locust_density_order.ipynb.

The raw per-animal tracking parts (``unwrap.data/density_order``) are
resampled to equal arc length before any modelling: every trajectory is
densely interpolated, then walked at successive points spaced by one
**step length** of 6 body lengths (chosen so consecutive steps are
approximately independent). Each consecutive pair of resampled positions
defines a step heading ``atan2(dY, dX)``.

The arc-length primitives (``rediscretize_trajectory`` and friends) live
in :mod:`unwrap.preprocess` (shared with the locust directional data
set); this module keeps the density-order-specific per-animal pipeline:
the ``order``/``density`` stimulus columns are carried through from the
raw table, and ``animal_id`` numbers the surviving animals densely (one
raw trajectory = one animal, under exactly one order-density condition).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from unwrap.preprocess import STEP_LENGTH, rediscretize_trajectory


def _process_one(x, y, fname, order, density, step_length):
    """Rediscretize one trajectory; returns its step table or None."""
    x_r, y_r = rediscretize_trajectory(x, y, step_length)
    if len(x_r) < 2:
        return None
    dx = np.diff(x_r)
    dy = np.diff(y_r)
    return pd.DataFrame(
        {
            "fname": fname,
            "order": order,
            "density": density,
            "X": x_r[1:],  # location at the end of each step
            "Y": y_r[1:],
            "dX": dx,  # step vectors
            "dY": dy,
            "direction": np.arctan2(dy, dx),  # radians
        }
    )


def rediscretize_density_order(
    df: pd.DataFrame,
    step_length: float = STEP_LENGTH,
    n_jobs: int | None = -1,
) -> pd.DataFrame:
    """Rediscretize the raw density-order trajectories into steps.

    Parameters
    ----------
    df:
        Raw per-animal table (see
        :func:`unwrap.datasets.load_density_order`); must contain
        ``X``, ``Y``, ``fname``, ``order`` and ``density``.
    step_length:
        Target separation between consecutive rediscretized points
        (``BODY_LENGTH * 6`` by default).
    n_jobs:
        Parallel workers for the per-animal rediscretization (``-1`` =
        all cores).

    Returns
    -------
    pd.DataFrame
        One row per step: ``fname``, ``order``, ``density``, ``X``,
        ``Y``, ``dX``, ``dY``, ``direction``, ``animal_id``.
        ``animal_id`` is dense 0..N-1 over the trajectories that
        survived rediscretization (a handful of raw recordings are too
        short or gappy to re-discretize at all).
    """
    if not {"X", "Y", "fname", "order", "density"}.issubset(df.columns):
        raise ValueError(
            f"df must contain X, Y, fname, order, density; got {list(df.columns)}"
        )

    tasks = []
    for fname, grp_raw in df.groupby("fname", sort=False):
        grp = grp_raw.dropna(subset=["X", "Y"])
        if len(grp) < 2:
            continue
        tasks.append(
            (
                grp["X"].to_numpy(dtype=np.float64, copy=False),
                grp["Y"].to_numpy(dtype=np.float64, copy=False),
                fname,
                float(grp["order"].iloc[0]),
                float(grp["density"].iloc[0]),
            )
        )

    results = Parallel(n_jobs=n_jobs)(
        delayed(_process_one)(x, y, fname, order, density, step_length)
        for x, y, fname, order, density in tasks
    )
    merged = pd.concat([r for r in results if r is not None], ignore_index=True)
    if merged.empty:
        raise ValueError("rediscretization produced no data")

    # dense animal ids, in file-scan order (the original numbering; the
    # models use it only as a categorical grouping)
    merged["animal_id"] = merged["fname"].map(
        {name: i for i, name in enumerate(merged["fname"].unique())}
    ).astype("int32")

    return merged[
        ["fname", "order", "density", "X", "Y", "dX", "dY", "direction", "animal_id"]
    ]
