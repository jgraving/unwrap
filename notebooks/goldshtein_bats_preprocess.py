"""Preprocessing for notebooks/goldshtein_bats.ipynb.

Raw GPS tracks (lat/lon) from the echolocating-bat homing study are
rediscretized to constant 60 m step length along the great-circle arc
(SLERP on unit vectors), with per-tag metadata (treatment and
destination) attached and goal-relative bearings (distance along the
heading and angle error to the destination) computed: for the point
goal (roost) and the polygon goal (Agamon Hula Lake, via a
harmonic/Laplace field solved on a metric grid).

The heavy per-stage implementations (SLERP rediscretization, metadata
merge, destination fields) are the annotated notebook cells moved here
verbatim; the notebook keeps one call per stage plus the descriptive
prints.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import geopandas as gpd
from concurrent.futures import ProcessPoolExecutor, as_completed
from shapely.geometry import Point, LineString
from shapely.prepared import prep
from tqdm.auto import tqdm

from unwrap.datasets import bats_lake_outline_path

# ------------------------- rediscretize -------------------------

# -------------------------
# Fallback spherical helpers
# -------------------------
_R_EARTH = 6371000.0

def _haversine_m(lat1, lon1, lat2, lon2):
    """Vectorized haversine distance in meters."""
    lat1 = np.deg2rad(lat1); lon1 = np.deg2rad(lon1)
    lat2 = np.deg2rad(lat2); lon2 = np.deg2rad(lon2)
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat/2.0)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2.0)**2
    return 2.0 * _R_EARTH * np.arcsin(np.sqrt(a))

def _ll_to_unitvec(lat_deg, lon_deg):
    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)
    cl = np.cos(lat)
    x = cl * np.cos(lon)
    y = cl * np.sin(lon)
    z = np.sin(lat)
    return np.stack([x, y, z], axis=-1)

def _unitvec_to_ll(v):
    v = v / np.linalg.norm(v, axis=-1, keepdims=True)
    x, y, z = v[..., 0], v[..., 1], v[..., 2]
    lat = np.arctan2(z, np.sqrt(x*x + y*y))
    lon = np.arctan2(y, x)
    return np.rad2deg(lat), np.rad2deg(lon)

def _slerp_unitvec(a, b, f):
    """
    Slerp between unit vectors a and b for fractions f in [0,1].
    a,b shape (N,3), f shape (N,)
    """
    a = a / np.linalg.norm(a, axis=1, keepdims=True)
    b = b / np.linalg.norm(b, axis=1, keepdims=True)
    dot = np.sum(a*b, axis=1)
    dot = np.clip(dot, -1.0, 1.0)
    omega = np.arccos(dot)
    so = np.sin(omega)

    out = np.empty_like(a)
    small = so < 1e-12
    out[small] = a[small] + (b[small] - a[small]) * f[small, None]
    k0 = np.sin((1.0 - f[~small]) * omega[~small]) / so[~small]
    k1 = np.sin(f[~small] * omega[~small]) / so[~small]
    out[~small] = k0[:, None]*a[~small] + k1[:, None]*b[~small]
    out = out / np.linalg.norm(out, axis=1, keepdims=True)
    return out


# -------------------------
# Internal: resample one contiguous segment (no gaps inside)
# -------------------------
def _rediscretize_contiguous_latlon_arrays(
    *,
    tag_val,
    lat: np.ndarray,
    lon: np.ndarray,
    t: np.ndarray,
    after: np.ndarray | None,
    spacing_m: float,
    tag_col: str,
    time_col: str,
    lat_col: str,
    lon_col: str,
    after_col: str,
    after_mode: str,
    include_first: bool,
    include_last: bool,
    segment_id: int,
    gap_before_segment: bool,
) -> pd.DataFrame:
    if len(lat) < 2:
        out = pd.DataFrame({
            tag_col: [tag_val],
            time_col: [t[0] if len(t) else np.nan],
            lat_col: [lat[0] if len(lat) else np.nan],
            lon_col: [lon[0] if len(lon) else np.nan],
            "s_along_m": [0.0],
            "segment_id": [segment_id],
            "gap_before": [1 if gap_before_segment else 0],
        })
        if after is not None and len(after):
            out[after_col] = int(after[0])
        return out

    # Prefer geodesic WGS84 if available
    try:
        from pyproj import Geod
        geod = Geod(ellps="WGS84")
        use_pyproj = True
    except Exception:
        geod = None
        use_pyproj = False

    # segment lengths (meters)
    if use_pyproj:
        az12, az21, dist = geod.inv(lon[:-1], lat[:-1], lon[1:], lat[1:])
        seg_len = dist.astype(float)
    else:
        seg_len = _haversine_m(lat[:-1], lon[:-1], lat[1:], lon[1:]).astype(float)

    keep = seg_len > 0
    if not np.any(keep):
        out = pd.DataFrame({
            tag_col: [tag_val],
            time_col: [t[0]],
            lat_col: [lat[0]],
            lon_col: [lon[0]],
            "s_along_m": [0.0],
            "segment_id": [segment_id],
            "gap_before": [1 if gap_before_segment else 0],
        })
        if after is not None:
            out[after_col] = int(after[0])
        return out

    seg_len2 = seg_len.copy()
    seg_len2[~keep] = 0.0
    cum = np.concatenate([[0.0], np.cumsum(seg_len2)])
    total = float(cum[-1])
    if total <= 0:
        out = pd.DataFrame({
            tag_col: [tag_val],
            time_col: [t[0]],
            lat_col: [lat[0]],
            lon_col: [lon[0]],
            "s_along_m": [0.0],
            "segment_id": [segment_id],
            "gap_before": [1 if gap_before_segment else 0],
        })
        if after is not None:
            out[after_col] = int(after[0])
        return out

    # target distances along-track (within THIS segment; reset to 0)
    start = 0.0 if include_first else float(spacing_m)
    targets = np.arange(start, total + 1e-9, float(spacing_m))
    if include_last and (len(targets) == 0 or targets[-1] < total):
        targets = np.concatenate([targets, [total]])

    seg_idx = np.searchsorted(cum, targets, side="right") - 1
    seg_idx = np.clip(seg_idx, 0, len(seg_len2) - 1)

    seg_start = cum[seg_idx]
    L = seg_len2[seg_idx]
    frac = np.zeros_like(targets, dtype=float)
    m = L > 0
    frac[m] = (targets[m] - seg_start[m]) / L[m]
    frac = np.clip(frac, 0.0, 1.0)

    # interpolate lat/lon at each target
    if use_pyproj:
        lon_out = np.empty_like(targets, dtype=float)
        lat_out = np.empty_like(targets, dtype=float)

        order = np.argsort(seg_idx)
        seg_sorted = seg_idx[order]
        frac_sorted = frac[order]

        az12_all = np.zeros_like(seg_len2)
        az12_all[:len(seg_len)] = az12
        az12_all[~keep] = 0.0

        k = 0
        while k < len(seg_sorted):
            i = seg_sorted[k]
            j = k + 1
            while j < len(seg_sorted) and seg_sorted[j] == i:
                j += 1

            d_along = frac_sorted[k:j] * seg_len2[i]
            lo, la, _ = geod.fwd(
                np.full_like(d_along, lon[i]),
                np.full_like(d_along, lat[i]),
                np.full_like(d_along, az12_all[i]),
                d_along
            )
            lon_out[order[k:j]] = lo
            lat_out[order[k:j]] = la
            k = j
    else:
        a = _ll_to_unitvec(lat[seg_idx], lon[seg_idx])
        b = _ll_to_unitvec(lat[seg_idx + 1], lon[seg_idx + 1])
        v = _slerp_unitvec(a, b, frac)
        lat_out, lon_out = _unitvec_to_ll(v)

    # interpolate time linearly by fraction within segment
    t0 = t[seg_idx]
    t1 = t[seg_idx + 1]
    t_out = (1.0 - frac) * t0 + frac * t1

    out = pd.DataFrame({
        tag_col: np.full_like(targets, tag_val, dtype=object),
        time_col: t_out,
        lat_col: lat_out,
        lon_col: lon_out,
        "s_along_m": targets,
        "segment_id": np.full_like(targets, segment_id, dtype=int),
        "gap_before": np.zeros_like(targets, dtype=int),
    })

    # mark the first output point of this segment if it follows a gap
    if gap_before_segment and len(out) > 0:
        out.loc[out.index[0], "gap_before"] = 1

    # assign afterLocalization (binary) at each target
    if after is not None:
        a0 = after[seg_idx]
        a1 = after[seg_idx + 1]
        if after_mode == "nearest":
            out[after_col] = np.where(frac >= 0.5, a1, a0).astype(int)
        elif after_mode == "forward":
            out[after_col] = a1.astype(int)
        else:  # backward
            out[after_col] = a0.astype(int)

    return out


# -------------------------
# Resample ONE tag trajectory (now with gap detection)
# -------------------------
def rediscretize_one_tag_latlon(
    g: pd.DataFrame,
    *,
    spacing_m: float,
    tag_col: str = "tag",
    time_col: str = "time",
    lat_col: str = "lat",
    lon_col: str = "lon",
    after_col: str = "afterLocalization",
    after_mode: str = "nearest",   # {"nearest","forward","backward"}
    include_first: bool = True,
    include_last: bool = False,
    # ---- NEW: gap handling ----
    gap_dist_m: float | None = None,
    gap_factor: float = 10.0,
    gap_time: float | None = None,        # in SAME units as time_col
    gap_time_s: float | None = None,      # convenience: seconds, converted via time_unit
    time_unit: str = "days",              # {"days","seconds"} only used if gap_time_s is given
) -> pd.DataFrame:
    """
    Spatially resample a single trajectory to constant step length spacing_m along-track.

    NEW: Detect "gaps" in the raw trajectory and split into contiguous segments.
         When a gap is detected between raw points i and i+1, we restart rediscretization at i+1
         (s_along_m resets to 0) and mark the first resampled point of the new segment with gap_before=1.

    Output keeps:
      - tag, time, lat, lon
      - afterLocalization (binary) assigned by nearest/forward/backward rule (if present)
      - s_along_m (RESET per contiguous segment)
      - segment_id (0,1,2,... within each tag)
      - gap_before (1 only on the FIRST resampled row after a detected gap)

    Gap detection:
      - distance gap: seg_len_m > gap_dist_m (default gap_factor * spacing_m)
      - time gap: dt > gap_time (or dt > gap_time_s converted using time_unit)
    """
    if spacing_m <= 0:
        raise ValueError("spacing_m must be > 0")
    if after_mode not in {"nearest", "forward", "backward"}:
        raise ValueError("after_mode must be one of: 'nearest', 'forward', 'backward'")
    if time_unit not in {"days", "seconds"}:
        raise ValueError("time_unit must be one of: 'days', 'seconds'")

    # defaults for gap thresholds
    if gap_dist_m is None:
        gap_dist_m = float(gap_factor) * float(spacing_m)

    if gap_time is None and gap_time_s is not None:
        if time_unit == "days":
            gap_time = float(gap_time_s) / 86400.0
        else:
            gap_time = float(gap_time_s)

    keep_cols = [tag_col, time_col, lat_col, lon_col]
    if after_col in g.columns:
        keep_cols.append(after_col)

    g = g[keep_cols].copy()
    g = g.sort_values(time_col) if time_col in g.columns else g.sort_index()

    # filter finite lat/lon/time
    lat_s = pd.to_numeric(g[lat_col], errors="coerce")
    lon_s = pd.to_numeric(g[lon_col], errors="coerce")
    ok = np.isfinite(lat_s.to_numpy()) & np.isfinite(lon_s.to_numpy())

    if time_col in g.columns:
        t_s = pd.to_numeric(g[time_col], errors="coerce")
        ok &= np.isfinite(t_s.to_numpy())

    g = g.loc[ok].copy()
    if len(g) < 2:
        # still add required columns for downstream consistency
        out = g.copy()
        out["s_along_m"] = np.nan
        out["segment_id"] = 0
        out["gap_before"] = 0
        return out

    tag_val = g[tag_col].iloc[0] if tag_col in g.columns else None

    lat = g[lat_col].to_numpy(dtype=float)
    lon = g[lon_col].to_numpy(dtype=float)
    t = g[time_col].to_numpy(dtype=float)

    # afterLocalization array (optional)
    if after_col in g.columns:
        after = pd.to_numeric(g[after_col], errors="coerce").fillna(0).astype(int).to_numpy()
    else:
        after = None

    # Prefer geodesic WGS84 if available for gap detection distances too
    try:
        from pyproj import Geod
        geod = Geod(ellps="WGS84")
        use_pyproj = True
    except Exception:
        geod = None
        use_pyproj = False

    # raw segment distances between successive points (meters)
    if use_pyproj:
        _, _, dist = geod.inv(lon[:-1], lat[:-1], lon[1:], lat[1:])
        seg_len_raw = dist.astype(float)
    else:
        seg_len_raw = _haversine_m(lat[:-1], lon[:-1], lat[1:], lon[1:]).astype(float)

    # raw time deltas
    dt = np.diff(t)

    # detect gaps between i and i+1
    gap_mask = seg_len_raw > float(gap_dist_m)
    if gap_time is not None:
        gap_mask = gap_mask | (dt > float(gap_time))

    # split into contiguous segments using gap_mask
    # gap_mask[k] True means break between raw index k and k+1
    break_idx = np.where(gap_mask)[0]  # indices in [0, n-2]
    starts = np.concatenate([[0], break_idx + 1])
    ends = np.concatenate([break_idx + 1, [len(lat)]])  # end is exclusive

    pieces = []
    for seg_id, (s0, s1) in enumerate(zip(starts, ends)):
        lat_seg = lat[s0:s1]
        lon_seg = lon[s0:s1]
        t_seg = t[s0:s1]
        after_seg = after[s0:s1] if after is not None else None

        # for segment>0, this segment begins after a gap
        gap_before_segment = (seg_id > 0)

        piece = _rediscretize_contiguous_latlon_arrays(
            tag_val=tag_val,
            lat=lat_seg,
            lon=lon_seg,
            t=t_seg,
            after=after_seg,
            spacing_m=float(spacing_m),
            tag_col=tag_col,
            time_col=time_col,
            lat_col=lat_col,
            lon_col=lon_col,
            after_col=after_col,
            after_mode=after_mode,
            include_first=include_first,
            include_last=include_last,
            segment_id=seg_id,
            gap_before_segment=gap_before_segment,
        )
        pieces.append(piece)

    out = pd.concat(pieces, ignore_index=True)

    # keep tag dtype consistent if possible
    try:
        out[tag_col] = out[tag_col].astype(g[tag_col].dtype)
    except Exception:
        pass

    # stable ordering
    out = out.sort_values([tag_col, time_col]).reset_index(drop=True)
    return out


# -------------------------
# Parallel wrapper + tqdm
# -------------------------
def rediscretize_latlon_parallel(
    raw_track_df: pd.DataFrame,
    *,
    spacing_m: float = 15.0,
    tag_col: str = "tag",
    time_col: str = "time",
    lat_col: str = "lat",
    lon_col: str = "lon",
    after_col: str = "afterLocalization",
    after_mode: str = "nearest",
    n_workers: int | None = None,
    include_first: bool = True,
    include_last: bool = False,
    # ---- NEW: gap handling ----
    gap_dist_m: float | None = None,
    gap_factor: float = 10.0,
    gap_time: float | None = None,
    gap_time_s: float | None = None,
    time_unit: str = "days",
) -> pd.DataFrame:
    groups = [(tag, g) for tag, g in raw_track_df.groupby(tag_col, sort=False)]
    results = []

    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        futs = {
            ex.submit(
                rediscretize_one_tag_latlon,
                g,
                spacing_m=float(spacing_m),
                tag_col=tag_col,
                time_col=time_col,
                lat_col=lat_col,
                lon_col=lon_col,
                after_col=after_col,
                after_mode=after_mode,
                include_first=include_first,
                include_last=include_last,
                gap_dist_m=gap_dist_m,
                gap_factor=gap_factor,
                gap_time=gap_time,
                gap_time_s=gap_time_s,
                time_unit=time_unit,
            ): tag
            for tag, g in groups
        }

        for fut in tqdm(as_completed(futs), total=len(futs),
                        desc=f"Rediscretizing lat/lon (Δ={spacing_m:g} m)"):
            tag = futs[fut]
            try:
                results.append(fut.result())
            except Exception as e:
                print(f"[tag {tag}] failed: {e}")

    if not results:
        return raw_track_df.iloc[0:0].copy()

    out = pd.concat(results, ignore_index=True)
    out = out.sort_values([tag_col, time_col]).reset_index(drop=True)
    return out


# -------------------------
# Quick check of step lengths (ignore segment boundaries)
# -------------------------
def check_step_lengths_latlon(
    df_resampled: pd.DataFrame,
    *,
    tag_col="tag",
    time_col="time",
    lat_col="lat",
    lon_col="lon",
    gap_col="gap_before",
):
    rows = []
    for tag, g in df_resampled.groupby(tag_col, sort=False):
        g = g.sort_values(time_col) if time_col in g.columns else g
        lat = g[lat_col].to_numpy(float)
        lon = g[lon_col].to_numpy(float)
        if len(lat) < 2:
            continue

        d = _haversine_m(lat[:-1], lon[:-1], lat[1:], lon[1:])

        # If gap_before exists, drop distances that cross a gap boundary:
        # the distance d[k] is between row k and k+1, so check gap_before at k+1.
        if gap_col in g.columns:
            gb = g[gap_col].to_numpy(int)
            cross = gb[1:] == 1
            d = d.copy()
            d[cross] = np.nan

        rows.append({
            tag_col: tag,
            "n_steps": int(np.sum(np.isfinite(d))),
            "median_step_m": float(np.nanmedian(d)),
            "min_step_m": float(np.nanmin(d)),
            "max_step_m": float(np.nanmax(d)),
        })
    return pd.DataFrame(rows)



def merge_tag_metadata(
    rediscretized_ll_df: pd.DataFrame,
    meta_df: pd.DataFrame,
) -> pd.DataFrame:
    """Attach per-tag metadata (treatment, destination) to the
    rediscretized steps.

    One row per tag is built from ``meta_df`` (duplicates are collapsed:
    ``treatment``/``treatmentName`` first, the ``flyHome``/``flyAgamon``
    destination flags maxed, and ``destination`` derived from the flags:
    Roost / Agamon / Unknown), with warnings if any tag is inconsistent
    (both flags set, or more than one treatment per tag). The result is
    the left merge with the rediscretized table (``many_to_one``), plus
    a combined ``treat_dest`` label ("treatmentName — destination").
    """

    # ---- 1) Build a clean per-tag metadata table ----
    tag_col = "tag"

    need_cols = [tag_col, "treatment", "treatmentName", "flyHome", "flyAgamon"]
    missing = [c for c in need_cols if c not in meta_df.columns]
    if missing:
        raise ValueError(f"meta_df is missing columns: {missing}")

    m = meta_df[need_cols].copy()

    # Coerce destination flags to 0/1
    for c in ["flyHome", "flyAgamon"]:
        m[c] = pd.to_numeric(m[c], errors="coerce").fillna(0).astype(int)

    # Optional: sanity check if both are 1 for any tag
    both = m.groupby(tag_col)[["flyHome", "flyAgamon"]].max()
    bad = both[(both["flyHome"] == 1) & (both["flyAgamon"] == 1)]
    if len(bad):
        print("WARNING: these tags have flyHome==1 AND flyAgamon==1 (check meta_df):")
        print(bad.index.to_list())

    # Create destination label
    def _dest_label(row):
        if row["flyHome"] == 1:
            return "Roost"
        if row["flyAgamon"] == 1:
            return "Agamon"
        return "Unknown"

    # Collapse to one row per tag (robust to duplicates)
    # - treatment / treatmentName should be consistent per tag; we take first but also check
    chk = m.groupby(tag_col).agg(
        n_treat=("treatment", lambda x: x.nunique(dropna=True)),
        n_name=("treatmentName", lambda x: x.nunique(dropna=True)),
    )
    if (chk["n_treat"] > 1).any() or (chk["n_name"] > 1).any():
        print("WARNING: Some tags have multiple treatments/treatmentNames in meta_df.")
        print(chk[(chk["n_treat"] > 1) | (chk["n_name"] > 1)].head(20))

    meta_by_tag = (
        m.groupby(tag_col, as_index=False)
         .agg(
             treatment=("treatment", "first"),
             treatmentName=("treatmentName", "first"),
             flyHome=("flyHome", "max"),
             flyAgamon=("flyAgamon", "max"),
         )
    )

    meta_by_tag["destination"] = meta_by_tag.apply(_dest_label, axis=1)

    # If you want a combined label (handy for plotting facets/legends)
    meta_by_tag["treat_dest"] = (
        meta_by_tag["treatmentName"].astype(str) + " — " + meta_by_tag["destination"].astype(str)
    )

    print("meta_by_tag rows:", len(meta_by_tag), "unique tags:", meta_by_tag[tag_col].nunique())


    # ---- 2) Merge into your rediscretized lat/lon dataframe ----
    # (use validate so pandas enforces many rows per tag in tracks, one row per tag in meta)
    rediscretized_ll_w_meta_df = rediscretized_ll_df.merge(
        meta_by_tag,
        on=tag_col,
        how="left",
        validate="many_to_one",
    )

    # # If you also want it on raw_track_df:
    # raw_track_w_meta_df = raw_track_df.merge(
    #     meta_by_tag,
    #     on=tag_col,
    #     how="left",
    #     validate="many_to_one",
    # )

    print("resampled rows:", len(rediscretized_ll_w_meta_df))
    print("missing treatment rows:",
          rediscretized_ll_w_meta_df["treatment"].isna().sum(),
          " (tags missing:",
          rediscretized_ll_w_meta_df.loc[rediscretized_ll_w_meta_df["treatment"].isna(), tag_col].nunique(),
          ")"
    )
    return rediscretized_ll_w_meta_df

# ============================================================
# Harmonic-field target bearing for a POLYGON goal (LAKE) OR:
# Smoothed-distance (EDT + Gaussian blur) target bearing.
#
# Lake destination ("Agamon"):
#   - dist_to_lake_m: true Euclidean distance to shoreline (polygon boundary)
#   - bearing_to_lake: bearing of the chosen target direction (−∇u)
#   - dist_to_lake_grad_m: distance along the (−∇u) ray to first lake intersection
#
# Home destination:
#   - distance + bearing to roost point (single point target)
#
# Output columns added (without dropping anything else):
#   x_m, y_m
#   dist_to_roost_m, dist_to_lake_m, dist_to_lake_grad_m
#   dist_to_dest_m (Euclidean for lake), dist_to_dest_grad_m (optional alt for lake)
#   dx_step_m, dy_step_m, step_len_m, step_bearing_deg
#   bearing_to_dest_deg
#   flightAngDiffDest_deg, flightAngDiffDest_rad
# ============================================================


# -------------------------
# basic helpers
# -------------------------
def wrap_deg(angle_deg: np.ndarray) -> np.ndarray:
    """Wrap degrees to (-180, 180]."""
    return (np.asarray(angle_deg, float) + 180.0) % 360.0 - 180.0

def bearing_deg_from_xy(dx: np.ndarray, dy: np.ndarray) -> np.ndarray:
    """
    Bearing in degrees FROM North, clockwise (0=N, 90=E).
    Using dx,dy in a metric projected CRS (e.g., UTM).
    """
    dx = np.asarray(dx, float)
    dy = np.asarray(dy, float)
    ang = np.degrees(np.arctan2(dx, dy))  # atan2(x,y) => clockwise from north
    return (ang + 360.0) % 360.0

def make_dest_label(flyHome: pd.Series, flyAgamon: pd.Series) -> pd.Series:
    """Choose destination label per row."""
    flyHome = pd.to_numeric(flyHome, errors="coerce").fillna(0).astype(int)
    flyAgamon = pd.to_numeric(flyAgamon, errors="coerce").fillna(0).astype(int)
    return pd.Series(
        np.select(
            [
                (flyHome == 1) & (flyAgamon == 0),
                (flyAgamon == 1) & (flyHome == 0),
                (flyAgamon == 1) & (flyHome == 1),
            ],
            ["Home", "Agamon", "Ambiguous"],
            default="Unknown",
        ),
        index=flyHome.index,
    )

# -------------------------
# polygon mask on grid (lake interior)
# -------------------------
def polygon_mask_on_grid(lake_poly_m, x_centers, y_centers, chunk=200_000):
    """
    Return boolean mask (ny,nx) where True means the grid point is inside the lake polygon.
    Uses shapely vectorized if available; otherwise chunked contains.
    """
    X, Y = np.meshgrid(x_centers, y_centers)  # (ny,nx)
    xx = X.ravel()
    yy = Y.ravel()

    # Fast path: shapely 2.x vectorized predicates (if present)
    try:
        import shapely
        pts = shapely.points(xx, yy)
        inside = shapely.contains(lake_poly_m, pts)
        return np.asarray(inside, dtype=bool).reshape(Y.shape)
    except Exception:
        pass

    # Older shapely: chunked loop with prepared geometry
    lake_prep = prep(lake_poly_m)
    inside = np.zeros(xx.size, dtype=bool)
    for a in range(0, xx.size, int(chunk)):
        b = min(xx.size, a + int(chunk))
        inside[a:b] = [lake_prep.contains(Point(x, y)) for x, y in zip(xx[a:b], yy[a:b])]
    return inside.reshape(Y.shape)

# -------------------------
# Solve Laplace outside lake, Dirichlet BC: lake=0, outer boundary=1
# -------------------------
def solve_laplace_dirichlet(mask_lake, *, max_iter=20_000, tol=1e-5):
    """
    Solve ∇²u=0 on grid cells where mask_lake=False and not on boundary.
    Dirichlet:
      - u=0 inside lake (mask_lake=True)
      - u=1 on outer boundary (edges)
    Returns u (ny,nx).

    Implementation:
      - Try SciPy sparse CG for speed/robustness
      - Fall back to Jacobi iteration (no extra deps)
    """
    ny, nx = mask_lake.shape

    # Fixed values
    u = np.ones((ny, nx), dtype=float)
    u[mask_lake] = 0.0

    # Outer boundary fixed to 1
    u[0, :] = 1.0
    u[-1, :] = 1.0
    u[:, 0] = 1.0
    u[:, -1] = 1.0

    # Unknowns are "free" cells excluding lake and excluding boundary edges
    free = (~mask_lake).copy()
    free[0, :] = free[-1, :] = free[:, 0] = free[:, -1] = False
    free_flat = np.flatnonzero(free.ravel())

    # ---------- fast path: SciPy sparse ----------
    try:
        import scipy.sparse as sp
        import scipy.sparse.linalg as spla

        idx = -np.ones((ny, nx), dtype=int)
        idx.ravel()[free_flat] = np.arange(free_flat.size)
        n = free_flat.size
        if n == 0:
            return u

        A = sp.lil_matrix((n, n), dtype=float)
        b = np.zeros(n, dtype=float)

        def add_neighbor(eq_i, r, c):
            j = idx[r, c]
            if j >= 0:   # unknown neighbor
                A[eq_i, j] -= 1.0
            else:        # known neighbor => rhs
                b[eq_i] += u[r, c]

        rr, cc = np.where(free)
        for r, c in zip(rr, cc):
            eq = idx[r, c]
            A[eq, eq] = 4.0
            add_neighbor(eq, r - 1, c)
            add_neighbor(eq, r + 1, c)
            add_neighbor(eq, r, c - 1)
            add_neighbor(eq, r, c + 1)

        sol, info = spla.cg(A.tocsr(), b, atol=0.0, rtol=1e-8, maxiter=3000)
        u2 = u.ravel()
        u2[free_flat] = sol
        return u2.reshape((ny, nx))

    except Exception:
        # ---------- fallback: Jacobi ----------
        u_old = u.copy()
        for _ in range(int(max_iter)):
            u_new = u_old.copy()

            avg4 = 0.25 * (
                u_old[:-2, 1:-1] + u_old[2:, 1:-1] +
                u_old[1:-1, :-2] + u_old[1:-1, 2:]
            )
            u_new[1:-1, 1:-1] = avg4

            # enforce Dirichlet constraints
            u_new[mask_lake] = 0.0
            u_new[0, :] = 1.0
            u_new[-1, :] = 1.0
            u_new[:, 0] = 1.0
            u_new[:, -1] = 1.0

            diff = np.nanmax(np.abs(u_new - u_old))
            u_old = u_new
            if diff < tol:
                break

        return u_old

# -------------------------
# Harmonic target vector at points (bearing + unit tx/ty)
# -------------------------
def harmonic_target_vector_to_polygon(
    x: np.ndarray,
    y: np.ndarray,
    lake_poly_m,
    *,
    grid_res_m: float = 20.0,
    domain_margin_m: float | None = None,
    jacobi_tol: float = 1e-5,
):
    """
    Compute harmonic-field target direction at each (x,y):

      - Solve Laplace outside lake: lake=0, outer boundary=1
      - Target vector t(x) = -grad(u)/||grad(u)||
      - Bearing = atan2(t_x, t_y) in degrees from North clockwise

    Returns:
      bearing_deg : (N,)
      tx, ty      : (N,) unit vector components (E,N) of t(x)
      ray_max_m   : scalar length long enough for ray casting across domain
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = x.size

    bearing = np.full(n, np.nan, dtype=float)
    tx_out = np.full(n, np.nan, dtype=float)
    ty_out = np.full(n, np.nan, dtype=float)

    if n == 0:
        return bearing, tx_out, ty_out, np.nan

    # Domain bounds: cover BOTH points and polygon bounds
    pxmin, pymin, pxmax, pymax = lake_poly_m.bounds
    xmin = np.nanmin([np.nanmin(x), pxmin])
    xmax = np.nanmax([np.nanmax(x), pxmax])
    ymin = np.nanmin([np.nanmin(y), pymin])
    ymax = np.nanmax([np.nanmax(y), pymax])

    span = max(xmax - xmin, ymax - ymin)
    if not np.isfinite(span) or span <= 0:
        return bearing, tx_out, ty_out, np.nan

    if domain_margin_m is None:
        domain_margin_m = max(1000.0, 0.10 * span)

    xmin_d = xmin - domain_margin_m
    xmax_d = xmax + domain_margin_m
    ymin_d = ymin - domain_margin_m
    ymax_d = ymax + domain_margin_m

    res = float(grid_res_m)
    if res <= 0:
        raise ValueError("grid_res_m must be > 0")

    nx = int(np.ceil((xmax_d - xmin_d) / res)) + 1
    ny = int(np.ceil((ymax_d - ymin_d) / res)) + 1
    x_centers = xmin_d + np.arange(nx) * res
    y_centers = ymin_d + np.arange(ny) * res

    ray_max_m = float(np.hypot(xmax_d - xmin_d, ymax_d - ymin_d))

    mask_lake = polygon_mask_on_grid(lake_poly_m, x_centers, y_centers)
    u = solve_laplace_dirichlet(mask_lake, tol=jacobi_tol)

    du_dy, du_dx = np.gradient(u, res, res)
    tx = -du_dx
    ty = -du_dy
    norm = np.hypot(tx, ty)
    tx = np.where(norm > 0, tx / norm, np.nan)
    ty = np.where(norm > 0, ty / norm, np.nan)

    # Bilinear sample at points
    fx = (x - x_centers[0]) / res
    fy = (y - y_centers[0]) / res
    valid = (fx >= 0) & (fx <= (nx - 1)) & (fy >= 0) & (fy <= (ny - 1))
    if not np.any(valid):
        return bearing, tx_out, ty_out, ray_max_m

    fxv = fx[valid]
    fyv = fy[valid]
    ix0 = np.floor(fxv).astype(int)
    iy0 = np.floor(fyv).astype(int)
    ix1 = np.clip(ix0 + 1, 0, nx - 1)
    iy1 = np.clip(iy0 + 1, 0, ny - 1)
    wx = fxv - ix0
    wy = fyv - iy0

    def bilinear(F):
        f00 = F[iy0, ix0]
        f10 = F[iy0, ix1]
        f01 = F[iy1, ix0]
        f11 = F[iy1, ix1]
        return (
            (1 - wx) * (1 - wy) * f00 +
            (wx) * (1 - wy) * f10 +
            (1 - wx) * (wy) * f01 +
            (wx) * (wy) * f11
        )

    tx_s = bilinear(tx)
    ty_s = bilinear(ty)

    nn = np.hypot(tx_s, ty_s)
    ok = np.isfinite(nn) & (nn > 0)
    tx_s = np.where(ok, tx_s / nn, np.nan)
    ty_s = np.where(ok, ty_s / nn, np.nan)

    tx_out[valid] = tx_s
    ty_out[valid] = ty_s

    b = np.full(tx_s.size, np.nan, dtype=float)
    b[ok] = bearing_deg_from_xy(tx_s[ok], ty_s[ok])
    bearing[valid] = b

    return bearing, tx_out, ty_out, ray_max_m

# -------------------------
# NEW: Smoothed-distance target vector at points (bearing + unit tx/ty)
# -------------------------
def smoothed_distance_target_vector_to_polygon(
    x: np.ndarray,
    y: np.ndarray,
    lake_poly_m,
    *,
    grid_res_m: float = 20.0,
    domain_margin_m: float | None = None,
    smooth_sigma_pix: float = 3.0,
):
    """
    Compute a smooth "efficiency-like" target direction at each (x,y):

      - Build EDT distance-to-lake on a grid (meters)
      - Gaussian smooth the distance field (sigma_pix in grid cells)
      - Target vector t(x) = -grad(d_smooth)/||grad(d_smooth)||
      - Bearing = atan2(t_x, t_y) in degrees from North clockwise

    Returns:
      bearing_deg : (N,)
      tx, ty      : (N,) unit vector components (E,N) of t(x)
      ray_max_m   : scalar length long enough for ray casting across domain
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = x.size

    bearing = np.full(n, np.nan, dtype=float)
    tx_out = np.full(n, np.nan, dtype=float)
    ty_out = np.full(n, np.nan, dtype=float)

    if n == 0:
        return bearing, tx_out, ty_out, np.nan

    # Domain bounds: cover BOTH points and polygon bounds
    pxmin, pymin, pxmax, pymax = lake_poly_m.bounds
    xmin = np.nanmin([np.nanmin(x), pxmin])
    xmax = np.nanmax([np.nanmax(x), pxmax])
    ymin = np.nanmin([np.nanmin(y), pymin])
    ymax = np.nanmax([np.nanmax(y), pymax])

    span = max(xmax - xmin, ymax - ymin)
    if not np.isfinite(span) or span <= 0:
        return bearing, tx_out, ty_out, np.nan

    if domain_margin_m is None:
        domain_margin_m = max(1000.0, 0.10 * span)

    xmin_d = xmin - domain_margin_m
    xmax_d = xmax + domain_margin_m
    ymin_d = ymin - domain_margin_m
    ymax_d = ymax + domain_margin_m

    res = float(grid_res_m)
    if res <= 0:
        raise ValueError("grid_res_m must be > 0")

    nx = int(np.ceil((xmax_d - xmin_d) / res)) + 1
    ny = int(np.ceil((ymax_d - ymin_d) / res)) + 1
    x_centers = xmin_d + np.arange(nx) * res
    y_centers = ymin_d + np.arange(ny) * res

    ray_max_m = float(np.hypot(xmax_d - xmin_d, ymax_d - ymin_d))

    mask_lake = polygon_mask_on_grid(lake_poly_m, x_centers, y_centers)

    try:
        from scipy.ndimage import distance_transform_edt, gaussian_filter
    except Exception as e:
        raise ImportError("smoothed-distance field requires scipy.ndimage") from e

    # EDT gives distance to nearest zero; we want zeros INSIDE lake => pass ~mask_lake
    # Use sampling for correctness even if res differs (here it doesn't).
    dist = distance_transform_edt(~mask_lake, sampling=(res, res)).astype(float)
    dist[mask_lake] = 0.0

    dist_s = gaussian_filter(dist, sigma=float(smooth_sigma_pix), mode="nearest")

    dd_dy, dd_dx = np.gradient(dist_s, res, res)
    tx = -dd_dx
    ty = -dd_dy
    norm = np.hypot(tx, ty)
    tx = np.where(norm > 0, tx / norm, np.nan)
    ty = np.where(norm > 0, ty / norm, np.nan)

    # Bilinear sample at points
    fx = (x - x_centers[0]) / res
    fy = (y - y_centers[0]) / res
    valid = (fx >= 0) & (fx <= (nx - 1)) & (fy >= 0) & (fy <= (ny - 1))
    if not np.any(valid):
        return bearing, tx_out, ty_out, ray_max_m

    fxv = fx[valid]
    fyv = fy[valid]
    ix0 = np.floor(fxv).astype(int)
    iy0 = np.floor(fyv).astype(int)
    ix1 = np.clip(ix0 + 1, 0, nx - 1)
    iy1 = np.clip(iy0 + 1, 0, ny - 1)
    wx = fxv - ix0
    wy = fyv - iy0

    def bilinear(F):
        f00 = F[iy0, ix0]
        f10 = F[iy0, ix1]
        f01 = F[iy1, ix0]
        f11 = F[iy1, ix1]
        return (
            (1 - wx) * (1 - wy) * f00 +
            (wx) * (1 - wy) * f10 +
            (1 - wx) * (wy) * f01 +
            (wx) * (wy) * f11
        )

    tx_s = bilinear(tx)
    ty_s = bilinear(ty)

    nn = np.hypot(tx_s, ty_s)
    ok = np.isfinite(nn) & (nn > 0)
    tx_s = np.where(ok, tx_s / nn, np.nan)
    ty_s = np.where(ok, ty_s / nn, np.nan)

    tx_out[valid] = tx_s
    ty_out[valid] = ty_s

    b = np.full(tx_s.size, np.nan, dtype=float)
    b[ok] = bearing_deg_from_xy(tx_s[ok], ty_s[ok])
    bearing[valid] = b

    return bearing, tx_out, ty_out, ray_max_m

# -------------------------
# Distance to FIRST lake intersection along a ray
# -------------------------
def distance_to_polygon_along_direction(
    x: np.ndarray,
    y: np.ndarray,
    tx: np.ndarray,
    ty: np.ndarray,
    lake_poly_m,
    *,
    ray_max_m: float,
):
    """
    For each point (x,y) and unit direction (tx,ty), compute the distance (meters)
    to the FIRST intersection of the ray with the lake polygon.

    Returns (N,) distances in meters (NaN if no intersection).
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    tx = np.asarray(tx, float)
    ty = np.asarray(ty, float)

    n = x.size
    out = np.full(n, np.nan, dtype=float)

    if n == 0 or not np.isfinite(ray_max_m) or ray_max_m <= 0:
        return out

    lake_prep = prep(lake_poly_m)

    def _collect_candidate_points(g):
        pts = []
        if g.is_empty:
            return pts
        gt = g.geom_type
        if gt == "Point":
            pts.append(g)
        elif gt == "MultiPoint":
            pts.extend(list(g.geoms))
        elif gt == "LineString":
            coords = list(g.coords)
            if coords:
                pts.append(Point(coords[0]))
                pts.append(Point(coords[-1]))
        elif gt == "MultiLineString":
            for part in g.geoms:
                coords = list(part.coords)
                if coords:
                    pts.append(Point(coords[0]))
                    pts.append(Point(coords[-1]))
        elif gt == "GeometryCollection":
            for part in g.geoms:
                pts.extend(_collect_candidate_points(part))
        else:
            for part in getattr(g, "geoms", []):
                pts.extend(_collect_candidate_points(part))
        return pts

    for i in range(n):
        if not (np.isfinite(x[i]) and np.isfinite(y[i]) and np.isfinite(tx[i]) and np.isfinite(ty[i])):
            continue
        if (tx[i] * tx[i] + ty[i] * ty[i]) < 1e-12:
            continue

        p0 = Point(float(x[i]), float(y[i]))

        # If already in/touching lake, set 0
        if lake_prep.contains(p0) or lake_poly_m.touches(p0):
            out[i] = 0.0
            continue

        p1 = (float(x[i] + tx[i] * ray_max_m), float(y[i] + ty[i] * ray_max_m))
        ray = LineString([(p0.x, p0.y), p1])

        inter = ray.intersection(lake_poly_m)
        if inter.is_empty:
            continue

        candidates = _collect_candidate_points(inter)
        if not candidates:
            continue

        best_s = np.inf
        for pt in candidates:
            s = (pt.x - p0.x) * tx[i] + (pt.y - p0.y) * ty[i]
            if s >= -1e-6 and s < best_s:
                best_s = s

        if np.isfinite(best_s):
            out[i] = float(max(0.0, best_s))

    return out

# ============================================================
# MAIN function you call
# ============================================================
def add_dist_and_angle_to_destination(
    rediscretized_ll_w_meta_df: pd.DataFrame,
    *,
    LAKE_GEOJSON = bats_lake_outline_path(),
    METRIC_EPSG: int = 32636,
    ROOST_LAT: float = 33.112373,
    ROOST_LON: float = 35.584135,
    # numerical knobs (not model params)
    lake_grid_res_m: float = 20.0,
    lake_domain_margin_m: float | None = None,
    lake_jacobi_tol: float = 1e-5,
    # NEW: choose target field for lake
    lake_field: str = "harmonic",          # {"harmonic","smoothed_distance"}
    lake_smooth_sigma_pix: float = 3.0,    # only used when lake_field="smoothed_distance"
    tag_col: str = "tag",
    time_col: str = "time",
    lat_col: str = "lat",
    lon_col: str = "lon",
):
    """
    Adds distances + bearings + angle error.

    Lake ("Agamon"):
      - dist_to_lake_m: Euclidean distance to shoreline
      - dist_to_lake_grad_m: distance along chosen field ray to first lake intersection
      - bearing_to_dest_deg: chosen field bearing

    Home:
      - distance + bearing to roost point
    """
    df = rediscretized_ll_w_meta_df.copy()

    need = [tag_col, time_col, lat_col, lon_col, "flyHome", "flyAgamon"]
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns: {missing}")

    df["flyHome"] = pd.to_numeric(df["flyHome"], errors="coerce").fillna(0).astype(int)
    df["flyAgamon"] = pd.to_numeric(df["flyAgamon"], errors="coerce").fillna(0).astype(int)
    df["dest"] = make_dest_label(df["flyHome"], df["flyAgamon"])

    # load lake polygon
    lake_4326 = gpd.read_file(LAKE_GEOJSON).to_crs(4326)
    lake_poly_4326 = lake_4326.unary_union
    lake_poly_m = gpd.GeoSeries([lake_poly_4326], crs=4326).to_crs(METRIC_EPSG).iloc[0]
    lake_boundary_m = lake_poly_m.boundary

    # roost point in metric CRS
    roost_pt_m = gpd.GeoSeries([Point(ROOST_LON, ROOST_LAT)], crs=4326).to_crs(METRIC_EPSG).iloc[0]
    ROOST_X, ROOST_Y = roost_pt_m.x, roost_pt_m.y

    # project track points to metric CRS
    df = df.dropna(subset=[lat_col, lon_col]).copy()
    pts_m = gpd.GeoSeries(
        gpd.points_from_xy(df[lon_col].astype(float), df[lat_col].astype(float)),
        crs=4326
    ).to_crs(METRIC_EPSG)

    df["x_m"] = pts_m.x.to_numpy()
    df["y_m"] = pts_m.y.to_numpy()

    # distance to roost
    dx_r = df["x_m"].to_numpy() - ROOST_X
    dy_r = df["y_m"].to_numpy() - ROOST_Y
    df["dist_to_roost_m"] = np.hypot(dx_r, dy_r)

    # true Euclidean distance to shoreline
    df["dist_to_lake_m"] = pts_m.distance(lake_boundary_m).to_numpy()

    # ---- choose lake target field ----
    field = str(lake_field).strip().lower()
    if field in ("harmonic", "laplace", "potential"):
        bearing_to_lake, lake_tx, lake_ty, ray_max_m = harmonic_target_vector_to_polygon(
            df["x_m"].to_numpy(),
            df["y_m"].to_numpy(),
            lake_poly_m,
            grid_res_m=lake_grid_res_m,
            domain_margin_m=lake_domain_margin_m,
            jacobi_tol=lake_jacobi_tol,
        )
    elif field in ("smoothed_distance", "smoothed", "edt", "distance"):
        bearing_to_lake, lake_tx, lake_ty, ray_max_m = smoothed_distance_target_vector_to_polygon(
            df["x_m"].to_numpy(),
            df["y_m"].to_numpy(),
            lake_poly_m,
            grid_res_m=lake_grid_res_m,
            domain_margin_m=lake_domain_margin_m,
            smooth_sigma_pix=lake_smooth_sigma_pix,
        )
    else:
        raise ValueError("lake_field must be one of {'harmonic','smoothed_distance'}")

    # distance to lake along chosen field ray
    df["dist_to_lake_grad_m"] = distance_to_polygon_along_direction(
        df["x_m"].to_numpy(),
        df["y_m"].to_numpy(),
        lake_tx,
        lake_ty,
        lake_poly_m,
        ray_max_m=ray_max_m,
    )

    # dist_to_dest: keep original Euclidean distance definition
    df["dist_to_dest_m"] = np.where(
        df["dest"].eq("Home"),
        df["dist_to_roost_m"],
        np.where(df["dest"].eq("Agamon"), df["dist_to_lake_m"], np.nan)
    )

    # alt distance consistent with chosen ray
    df["dist_to_dest_grad_m"] = np.where(
        df["dest"].eq("Home"),
        df["dist_to_roost_m"],
        np.where(df["dest"].eq("Agamon"), df["dist_to_lake_grad_m"], np.nan)
    )

    # step vectors and bearing (within tag)
    df = df.sort_values([tag_col, time_col]).reset_index(drop=True)
    df["dx_step_m"] = df.groupby(tag_col)["x_m"].diff().to_numpy()
    df["dy_step_m"] = df.groupby(tag_col)["y_m"].diff().to_numpy()
    df["step_len_m"] = np.hypot(df["dx_step_m"], df["dy_step_m"])

    df["step_bearing_deg"] = bearing_deg_from_xy(
        df["dx_step_m"].to_numpy(),
        df["dy_step_m"].to_numpy()
    )
    df.loc[~np.isfinite(df["step_len_m"]) | (df["step_len_m"] <= 0), "step_bearing_deg"] = np.nan

    # bearing to roost (point)
    to_roost_dx = ROOST_X - df["x_m"].to_numpy()
    to_roost_dy = ROOST_Y - df["y_m"].to_numpy()
    bearing_to_roost = bearing_deg_from_xy(to_roost_dx, to_roost_dy)

    # bearing to dest per row
    df["bearing_to_dest_deg"] = np.where(
        df["dest"].eq("Home"),
        bearing_to_roost,
        np.where(df["dest"].eq("Agamon"), bearing_to_lake, np.nan)
    )

    # angle error (circular, centered at 0)
    df["flightAngDiffDest_deg"] = wrap_deg(df["step_bearing_deg"] - df["bearing_to_dest_deg"])
    df["flightAngDiffDest_rad"] = np.deg2rad(df["flightAngDiffDest_deg"])

    return df
