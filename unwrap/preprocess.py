"""
Rediscretization pipeline for the locust directional dataset.

Ported from ``jgraving/sayin_locust_mixture_model/locust_mixture_model.ipynb``
(cells defining ``densely_interpolate``, ``find_next_point``,
``rediscretize_trajectory``, ``process_individual``, and the per-animal
merge). The unwrap directional notebook only loads a baked parquet; this
module restores the raw -> merged preprocessing so the pipeline is
reproducible end to end from ``unwrap-data/directional``.

Constants
---------
``BODY_LENGTH = 0.04`` (locust body length in the experiment units) and
``STEP_LENGTH = 6 * BODY_LENGTH``. Consecutive rediscretized points are
separated by approximately ``STEP_LENGTH``.

Main entry point
----------------
:func:`rediscretize_directional` -- takes the raw per-animal DataFrame
(from :func:`unwrap.datasets.load_directional`) and returns the
same ``merged_df`` the original notebook produced, i.e. columns
``X, Y, dX, dY, fname, direction, animal_id, file_path``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from numba import njit

BODY_LENGTH = 0.04
TOLERANCE = 0.0001  # Tolerance for the difference in step length
STEP_LENGTH = BODY_LENGTH * 6


@njit
def densely_interpolate(x, y, tolerance=TOLERANCE):
    """Linearly interpolate points between given trajectory points."""
    n = len(x)
    counts = np.zeros(n - 1, dtype=np.int64)

    total_points = 0
    for i in range(n - 1):
        segment_length = np.sqrt((x[i + 1] - x[i]) ** 2 + (y[i + 1] - y[i]) ** 2)
        num_points = max(int(segment_length / tolerance), 2)
        counts[i] = num_points
        if i == 0:
            total_points += num_points
        else:
            total_points += num_points - 1

    new_x = np.empty(total_points, dtype=np.float64)
    new_y = np.empty(total_points, dtype=np.float64)

    idx = 0
    for i in range(n - 1):
        num_points = counts[i]
        seg_x = np.linspace(x[i], x[i + 1], num_points)
        seg_y = np.linspace(y[i], y[i + 1], num_points)
        if i == 0:
            for j in range(num_points):
                new_x[idx] = seg_x[j]
                new_y[idx] = seg_y[j]
                idx += 1
        else:
            for j in range(1, num_points):
                new_x[idx] = seg_x[j]
                new_y[idx] = seg_y[j]
                idx += 1

    return new_x, new_y


@njit
def find_next_point(x, y, start_x, start_y, step_length):
    """Find the next point approximately ``step_length`` away from the start."""
    distances = np.sqrt((x - start_x) ** 2 + (y - start_y) ** 2)
    diff = np.abs(distances - step_length)
    index = np.argmin(diff)
    return x[index], y[index], index


def rediscretize_trajectory(x, y, step_length, tolerance=TOLERANCE):
    """Rediscretize a trajectory into points at approximately equal intervals."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)

    dense_x, dense_y = densely_interpolate(x, y, tolerance)
    new_x = [x[0]]
    new_y = [y[0]]

    while True:
        if len(dense_x) == 0:
            break
        next_x, next_y, index = find_next_point(
            dense_x, dense_y, new_x[-1], new_y[-1], step_length
        )

        last_distance = np.sqrt((new_x[-1] - next_x) ** 2 + (new_y[-1] - next_y) ** 2)
        if last_distance < step_length - TOLERANCE:
            break

        new_x.append(next_x)
        new_y.append(next_y)

        dense_x = dense_x[index:]
        dense_y = dense_y[index:]

        if index == 0:  # No progress; prevent infinite loop
            break

    return np.array(new_x), np.array(new_y)


def _process_one(x, y, fname, step_length):
    """Rediscretize one trajectory and compute its directional differences.

    Returns a DataFrame with ``X, Y, dX, dY, fname, direction`` or ``None``
    when the trajectory is too short to rediscretize.
    """
    x_r, y_r = rediscretize_trajectory(x, y, step_length)
    if len(x_r) < 2:
        return None
    dx = np.diff(x_r)
    dy = np.diff(y_r)
    return pd.DataFrame(
        {
            "X": x_r[1:],
            "Y": y_r[1:],
            "dX": dx,
            "dY": dy,
            "fname": fname,
            "direction": np.arctan2(dy, dx),
        }
    )


def rediscretize_directional(
    df: pd.DataFrame,
    step_length: float = STEP_LENGTH,
    n_jobs: int | None = -1,
) -> pd.DataFrame:
    """Rediscretize raw per-animal trajectories.

    Parameters
    ----------
    df:
        Raw directional data (see :func:`unwrap.datasets.load_directional`),
        must contain ``X``, ``Y`` and ``fname`` columns.
    step_length:
        Target separation between consecutive rediscretized points
        (``BODY_LENGTH * 6`` by default).
    n_jobs:
        Parallel workers for the per-animal rediscretization (``-1`` = all
        cores).

    Returns
    -------
    pd.DataFrame
        Columns ``X, Y, dX, dY, fname, direction, animal_id, file_path``.
        ``file_path`` records which of the two original source tables an
        animal came from (``XY_raw_50cm_`` or ``XY_raw_8cm__``), mirroring
        the original notebook's metadata.
    """
    if not {"X", "Y", "fname"}.issubset(df.columns):
        raise ValueError(f"df must contain X, Y, fname; got {list(df.columns)}")

    tasks = [
        (grp["X"].values, grp["Y"].values, name)
        for name, grp in df.groupby("fname", sort=False)
    ]
    results = Parallel(n_jobs=n_jobs)(
        delayed(_process_one)(x, y, fname, step_length)
        for x, y, fname in tasks
    )
    merged = pd.concat([r for r in results if r is not None], ignore_index=True)
    if merged.empty:
        raise ValueError("rediscretization produced no data")

    # animal_id: 0-based, in the same order the input listed animals in.
    order = list(pd.unique(df["fname"]))
    merged["animal_id"] = merged["fname"].map(
        {name: i for i, name in enumerate(order)}
    ).astype("int32")

    # file_path echoes the original source table (50cm / 8cm), carried on
    # the input df as ``source_table`` so results stay comparable to the
    # original notebook's merged_df.
    if "source_table" in df.columns:
        src = df[["fname", "source_table"]].drop_duplicates("fname")
        merged = merged.drop(columns=["file_path"], errors="ignore")
        merged = merged.merge(
            src.rename(columns={"source_table": "file_path"}), on="fname", how="left"
        )
    else:
        merged["file_path"] = "unknown"

    # Keep the same column order as the historical baked
    # ``bimodal_locust.parquet`` (X, Y, dX, dY, direction, animal_id)
    # so downstream code that selects columns by name keeps working.
    return merged[["X", "Y", "dX", "dY", "fname", "direction", "animal_id", "file_path"]].reset_index(drop=True)
