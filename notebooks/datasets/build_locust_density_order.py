#!/usr/bin/env python3
"""Rebuild the bundled rediscretized locust swarm dataset.

The notebook ``notebooks/locust_density_order.ipynb`` ships with the
preprocessed rediscretized heading steps
(``notebooks/datasets/locust_density_order.parquet``) so it can run from a
fresh checkout. This script reproduces that file from the raw per-animal
trajectory parquet shards (``XY_parquet_parts/``), which are not bundled in
the repository (see the Sayin et al. 2025 supplementary material).

Usage (from the repository root, path to the raw shards as the argument):

    python notebooks/datasets/build_locust_density_order.py PATH/TO/XY_parquet_parts \
        [-o notebooks/datasets/locust_density_order.parquet]

Each shard is one animal's arena positions. The pipeline:

1. concatenates all shards;
2. densely interpolates each trajectory on an arc-length grid;
3. re-samples it at a fixed step length of 6 body lengths (BODY_LENGTH = 0.04
   m, chosen by a lag-1 circular-autocorrelation screen over 1..10 body
   lengths for approximately independent steps);
4. emits one row per step: position, displacement, and heading
   ``atan2(dY, dX)`` in radians, with the stimulus ``order`` and ``density``
   and a zero-based per-animal ``animal_id``.

Output columns: fname, order, density, X, Y, dX, dY, direction, animal_id.
"""
import argparse
import os
import glob as _glob

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

BODY_LENGTH = 0.04
TOLERANCE = 1e-4   # max spacing for interpolation step
STEP_LENGTH = BODY_LENGTH * 6


# Optional: use numba if available
try:
    from numba import njit
except Exception:
    def njit(*args, **kwargs):
        def wrap(f): return f
        return wrap

BODY_LENGTH = 0.04
TOLERANCE = 1e-4   # max spacing for interpolation step
STEP_LENGTH = BODY_LENGTH * 6

@njit
def densely_interpolate(x, y, tolerance=TOLERANCE):
    n = len(x)
    counts = np.zeros(n-1, dtype=np.int64)

    total_points = 0
    for i in range(n - 1):
        dx = x[i+1] - x[i]
        dy = y[i+1] - y[i]
        segment_length = np.sqrt(dx*dx + dy*dy)
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
        seg_x = np.linspace(x[i], x[i+1], num_points)
        seg_y = np.linspace(y[i], y[i+1], num_points)
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
    distances = np.sqrt((x - start_x)**2 + (y - start_y)**2)
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
        last_distance = np.sqrt(last_dx*last_dx + last_dy*last_dy)
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

    out = pd.DataFrame({
        "fname": fname,
        "order": o,
        "density": d,
        "X": xr[1:],              # location at the end of each step
        "Y": yr[1:],
        "dX": dx,                 # step vectors
        "dY": dy,
        "direction": direction    # radians
    })
    return out

# ---- Distance-only filter ----
def _total_distance_m(g_after_dropna):
    dx = g_after_dropna["X"].diff()
    dy = g_after_dropna["Y"].diff()
    return float(np.hypot(dx, dy).fillna(0.0).sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("parquet_dir", help="directory containing the raw trajectory parquet shards")
    ap.add_argument("-o", "--output",
                    default="notebooks/datasets/locust_density_order.parquet")
    args = ap.parse_args()

    files = sorted(_glob.glob(os.path.join(args.parquet_dir, "*.parquet")))
    if not files:
        raise SystemExit(f"No Parquet files found in {args.parquet_dir}")

    good_files, bad_files = [], []
    for f in files:
        try:
            pq.ParquetFile(f)
        except Exception as e:
            bad_files.append((f, str(e)))
            continue
        good_files.append(f)
    for f, e in bad_files:
        print(f"Skipping invalid/corrupted file: {f} -> {e}")

    df = pd.concat([pd.read_parquet(f) for f in good_files], ignore_index=True)
    print(f"loaded {len(df)} rows from {len(good_files)} shards")

    results = []
    for fname, grp_raw in df.groupby("fname", sort=False):
        grp = grp_raw.dropna(subset=["X", "Y"])
        if len(grp) < 2:
            continue
        res = process_individual_group(grp, fname, step_length=STEP_LENGTH)
        if res is not None and not res.empty:
            results.append(res)

    if not results:
        raise SystemExit("No groups produced rediscretized output (check the shards).")

    merged_df = pd.concat(results, ignore_index=True)
    merged_df["unique_id"] = (merged_df["order"].astype(str) + "|"
                              + merged_df["density"].astype(str) + "|"
                              + merged_df["fname"].astype(str))
    id_map = {u: i for i, u in enumerate(merged_df["unique_id"].unique())}
    merged_df["animal_id"] = merged_df["unique_id"].map(id_map).astype(int)
    merged_df = merged_df.drop(columns=["unique_id"])

    merged_df.to_parquet(args.output, index=False)
    print(f"wrote {args.output}: {merged_df.shape[0]} steps from "
          f"{merged_df['animal_id'].nunique()} animals")


if __name__ == "__main__":
    main()
