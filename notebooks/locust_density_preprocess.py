"""
Preprocessing for the locust density/order notebook.

Turns the raw per-trial trajectories (``load_density_order``: one block per
``fname`` with ``X, Y, loss, order, density, groups``) into the two frames
the notebook built inline:

* ``merged_df``  -- every trial rediscretized to ``STEP_LENGTH`` = 6 body
  lengths along-track (with a distance-only filter), carrying the trial's
  ``order`` / ``density``, the step vectors ``dX/dY``, the heading
  ``direction`` (rad), and a dense numeric ``animal_id`` per trial
  (``order|density|fname``) -- :func:`rediscretize_density_order`.
* ``angles_df``  -- the per-trial **mean movement vector** (mean ``dX, dY``)
  and its heading ``angle_rad`` / ``angle_deg`` / ``magnitude``, i.e. the
  table the von Mises model cell fits (``y ~ t_density * t_order``) and the
  alignment figures draw -- :func:`mean_vector_angles`.

Everything here is deterministic geometry; no model fitting.
:func:`build_density_table` runs the whole pipeline and returns both frames.
The trajectory functions are lifted verbatim from the notebook's original
cell, so ``merged_df`` is bit-identical to the pre-refactor output.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Rediscretization (verbatim from the notebook)
# ---------------------------------------------------------------------------

try:
    from numba import njit
except Exception:

    def njit(*args, **kwargs):
        def wrap(f):
            return f

        return wrap


BODY_LENGTH = 0.04
TOLERANCE = 1e-4  # max spacing for interpolation step
STEP_LENGTH = BODY_LENGTH * 6


@njit
def densely_interpolate(x, y, tolerance=TOLERANCE):
    n = len(x)
    counts = np.zeros(n - 1, dtype=np.int64)

    total_points = 0
    for i in range(n - 1):
        dx = x[i + 1] - x[i]
        dy = y[i + 1] - y[i]
        segment_length = np.sqrt(dx * dx + dy * dy)
        num_points = max(int(segment_length / tolerance), 2)
        counts[i] = num_points
        if i == 0:
            total_points += num_points
        else:
            total_points += (num_points - 1)

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
            for j in range(1, num_points):  # skip duplicate
                new_x[idx] = seg_x[j]
                new_y[idx] = seg_y[j]
                idx += 1

    return new_x, new_y


@njit
def find_next_point(x, y, start_x, start_y, step_length):
    distances = np.sqrt((x - start_x) ** 2 + (y - start_y) ** 2)
    diff = np.abs(distances - step_length)
    index = np.argmin(diff)
    return x[index], y[index], index


def rediscretize_trajectory(x, y, step_length, tolerance=TOLERANCE):
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if x.size < 2:
        return np.array([x[0]]) if x.size == 1 else np.array([]), np.array([y[0]]) if y.size == 1 else np.array([])

    dense_x, dense_y = densely_interpolate(x, y, tolerance)
    if dense_x.size == 0:
        return np.array([]), np.array([])

    new_x = [x[0]]
    new_y = [y[0]]

    while True:
        if len(dense_x) == 0:
            break
        next_x, next_y, index = find_next_point(dense_x, dense_y, new_x[-1], new_y[-1], step_length)

        last_dx = new_x[-1] - next_x
        last_dy = new_y[-1] - next_y
        last_distance = np.sqrt(last_dx * last_dx + last_dy * last_dy)
        if last_distance < step_length - TOLERANCE:
            break

        new_x.append(next_x)
        new_y.append(next_y)

        dense_x = dense_x[index:]
        dense_y = dense_y[index:]

        if index == 0:  # no progress safeguard
            break

    return np.array(new_x), np.array(new_y)


def process_individual_group(grp, fname, step_length=STEP_LENGTH):
    """
    grp: DataFrame for a single fname from dfXY with columns X,Y,order,density
    """
    x = grp["X"].to_numpy(dtype=np.float64, copy=False)
    y = grp["Y"].to_numpy(dtype=np.float64, copy=False)

    xr, yr = rediscretize_trajectory(x, y, step_length)
    if xr.size < 2:
        return None

    dx = np.diff(xr)
    dy = np.diff(yr)
    direction = np.arctan2(dy, dx)

    o = grp["order"].iloc[0] if "order" in grp.columns else None
    d = grp["density"].iloc[0] if "density" in grp.columns else None

    out = pd.DataFrame(
        {
            "fname": fname,
            "order": o,
            "density": d,
            "X": xr[1:],  # location at the end of each step
            "Y": yr[1:],
            "dX": dx,  # step vectors
            "dY": dy,
            "direction": direction,  # radians
        }
    )
    return out


def _total_distance_m(g_after_dropna):
    dx = g_after_dropna["X"].diff()
    dy = g_after_dropna["Y"].diff()
    return float(np.hypot(dx, dy).fillna(0.0).sum())


def rediscretize_density_order(dfXY):
    """Run the notebook's per-trial rediscretization and return ``merged_df``.

    Mirrors the inline cell exactly: drop NA rows per trial, rediscretize at
    ``STEP_LENGTH``, keep trials that produced >= 2 points, concatenate, then
    assign the dense numeric ``animal_id`` per ``order|density|fname``.
    """
    results = []
    for fname, grp_raw in dfXY.groupby("fname", sort=False):
        # Drop NA rows in X/Y for movement + rediscretization
        grp = grp_raw.dropna(subset=["X", "Y"])
        if len(grp) < 2:
            continue

        # Movement in meters
        dist_m = _total_distance_m(grp)

        # # ---- FILTER: keep only if moved >= MIN_DISTANCE_M and completed at least 80% of trial----
        # if dist_m < 1.0 or grp_raw["loss"].max() > 0.2:
        #     continue

        # proceed to rediscretization
        res = process_individual_group(grp, fname, step_length=STEP_LENGTH)
        if res is not None and not res.empty:
            results.append(res)

    if not results:
        raise RuntimeError(
            "No groups produced rediscretized output after filtering (check dfXY content)."
        )

    merged_df = pd.concat(results, ignore_index=True)

    # Optional: numeric animal_id per (order, density, fname)
    merged_df["unique_id"] = (
        merged_df["order"].astype(str) + "|" +
        merged_df["density"].astype(str) + "|" +
        merged_df["fname"].astype(str)
    )
    id_map = {u: i for i, u in enumerate(merged_df["unique_id"].unique())}
    merged_df["animal_id"] = merged_df["unique_id"].map(id_map).astype(int)
    merged_df = merged_df.drop(columns=["unique_id"])
    return merged_df


def mean_vector_angles(merged_df):
    """Per-trial mean movement vector and its heading -> ``angles_df``.

    The table the von Mises direction model fits (``y ~ t_density * t_order``
    after scaling) and the alignment figures draw.
    """
    # Mean (dX, dY) per individual per condition
    angles_df = (
        merged_df
        .copy()
        .groupby(["order", "density", "fname"], as_index=False)[["dX", "dY"]]
        .mean()
        .rename(columns={"dX": "mean_dX", "dY": "mean_dY"})
    )

    # Angle (no normalization), magnitude for reference
    angles_df["angle_rad"] = np.arctan2(angles_df["mean_dY"], angles_df["mean_dX"])
    angles_df["angle_deg"] = np.degrees(angles_df["angle_rad"])
    angles_df["magnitude"] = np.hypot(angles_df["mean_dX"], angles_df["mean_dY"])

    # Optional: sort nicely
    angles_df = angles_df.sort_values(["order", "density", "fname"], kind="stable").reset_index(drop=True)
    return angles_df


def build_density_table(dfXY):
    """Run the full deterministic pipeline; return ``merged_df`` and ``angles_df``."""
    merged_df = rediscretize_density_order(dfXY)
    angles_df = mean_vector_angles(merged_df)
    return {"merged_df": merged_df, "angles_df": angles_df}
