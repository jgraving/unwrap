# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy", "pandas", "pyarrow", "rasterio", "pillow"]
# ///
"""Reproduce the zebra map-figure data files for the jgraving/unwrap-data repo.

Converts the original HerdHover ``circular_stats_data`` sources into the
git-friendly layout used by ``unwrap.datasets.load_zebra_mosaic`` and
``unwrap.datasets.load_zebra_postures`` (see the ``zebra_bearings/`` section
of the unwrap-data README):

- ``postures/observation088-{track}.parquet``
    The (18, T, 10, 2) posture-keypoint array of
    ``observation088-postures.npy`` split per track into wide float64
    parquets (frame + kp0_x..kp9_y) so every file stays within git's
    per-file size limit. Conversion is verified lossless: every parquet is
    re-read and compared element-for-element against the npy slice
    (``np.array_equal(..., equal_nan=True)``).

- ``mosaic/manifest.json`` + ``mosaic/tile_{i:03d}.png``
    A native-resolution lossless crop of the aerial orthomosaic
    (RGBA, 0.02707 m/px, UTM zone 37N) covering the window the escape-map
    figures use, padded for arrows, split into horizontal row-strips saved
    as PNG (lossless) with a JSON manifest carrying the extent, CRS and
    per-tile row offsets. Each PNG is checked to decode to exactly the
    pixels that were written.

The small per-frame *bearings* the map figures overlay are already in
``zebra_bearings/raw/*.parquet`` (columns ``compass_bearing`` and
``pre_scare``), so ``ob088_bearings.csv`` is not duplicated here.

Usage:

    uv run tools/prepare_zebra_map_data.py \\
        --postures-npy .../observation088-postures.npy \\
        --mosaic-tif .../riverside4_5_11_merge_transparent_mosaic_group1.tif \\
        --out-dir .../unwrap-data/zebra_bearings

The crop window is the one shared by every published escape-map panel
(x 267525-267680, y 37290-37770 m UTM 37N) padded by --pad-m.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from PIL import Image
from rasterio.windows import from_bounds

# UTM 37N meters; union of every window used by the published panels.
CROP_XMIN, CROP_XMAX = 267525.0, 267680.0
CROP_YMIN, CROP_YMAX = 37290.0, 37770.0

GIT_FILE_LIMIT = 100 * 1024 * 1024 - 1  # github hard per-file limit


def npy_to_parquets(npy_path: Path, out_dir: Path) -> None:
    """Split the posture npy into per-track float64 parquets, verify lossless."""
    arr = np.load(npy_path, mmap_mode="r")
    n_tracks, n_frames, n_kps, n_xy = arr.shape
    assert arr.dtype == np.float64 and n_xy == 2, \
        f"expected float64 (tracks, frames, keypoints, 2), got {arr.dtype} {arr.shape}"
    obs = npy_path.stem.split("-")[0]
    out_dir.mkdir(parents=True, exist_ok=True)
    kp_cols = [f"kp{k}_{c}" for k in range(n_kps) for c in ("x", "y")]
    for t in range(n_tracks):
        track = f"{obs}-{t}"
        pos = np.asarray(arr[t], dtype=np.float64)  # (frames, kps, 2)
        frame = np.arange(n_frames, dtype=np.int64)
        flat = pos.reshape(n_frames, n_kps * 2)
        df = pd.DataFrame(flat, columns=kp_cols)
        df.insert(0, "frame", frame)
        df.insert(1, "track_number", track)
        out = out_dir / f"{track}.parquet"
        df.to_parquet(out, index=False)
        got = pd.read_parquet(out)
        back = got[kp_cols].to_numpy(dtype=np.float64)
        assert np.array_equal(back, flat, equal_nan=True), \
            f"round-trip mismatch for {track}"
        assert list(got["frame"].to_numpy()) == list(frame)
        print(f"  {out.name}: {n_frames} frames, {out.stat().st_size/1e6:.1f} MB (verified)")


def crop_mosaic_to_tiles(tif_path: Path, out_dir: Path, *, pad_m: float,
                         tile_rows: int, tile_bytes_budget: int = 95 * 1024 * 1024) -> None:
    """Crop the mosaic at native resolution; write RGBA PNG row-strips + manifest."""
    with rasterio.open(tif_path) as src:
        assert src.count == 4, f"expected RGBA, got {src.count} bands"
        tr = src.transform
        win = from_bounds(CROP_XMIN - pad_m, CROP_YMIN - pad_m,
                          CROP_XMAX + pad_m, CROP_YMAX + pad_m, transform=tr)
        win = win.intersection(rasterio.windows.from_bounds(
            *src.bounds, transform=tr))
        img = src.read(window=win)  # (4, H, W) uint8
        wt = rasterio.windows.transform(win, tr)
        left, top = wt.c, wt.f
        n_cols = img.shape[2]
        n_rows = img.shape[1]
        extent = {
            "left": left,
            "right": left + wt.a * n_cols,
            "bottom": top + wt.e * n_rows,
            "top": top,
        }
        crs_epsg = src.crs.to_epsg()
        crs_wkt = src.crs.to_wkt()
        pixel_size = abs(tr.a)

    img = np.ascontiguousarray(np.transpose(img, (1, 2, 0)))  # (H, W, 4)
    print(f"  crop {n_cols}x{n_rows} px, extent "
          f"({extent['left']:.2f}, {extent['right']:.2f}) x "
          f"({extent['bottom']:.2f}, {extent['top']:.2f})")

    out_dir.mkdir(parents=True, exist_ok=True)
    tiles = []
    for i, row0 in enumerate(range(0, n_rows, tile_rows)):
        nrows = min(tile_rows, n_rows - row0)
        strip = Image.fromarray(img[row0:row0 + nrows], mode="RGBA")
        path = out_dir / f"tile_{i:03d}.png"
        strip.save(path, format="PNG", compress_level=9)
        back = np.asarray(Image.open(path).convert("RGBA"))
        assert np.array_equal(back, img[row0:row0 + nrows]), \
            f"PNG round-trip mismatch for {path.name}"
        size = path.stat().st_size
        assert size <= tile_bytes_budget, \
            f"{path.name} is {size/1e6:.1f} MB > {tile_bytes_budget/1e6:.0f} MB; raise --tile-rows"
        assert size <= GIT_FILE_LIMIT, f"{path.name} exceeds the git per-file limit"
        tiles.append({"file": path.name, "row0": row0, "nrows": nrows,
                      "bytes": size})
        print(f"  {path.name}: rows {row0}-{row0+nrows}, {size/1e6:.1f} MB (verified)")

    manifest = {
        "format": "zebra-mosaic-tiles-v1",
        "source": tif_path.name,
        "crs": {"epsg": crs_epsg, "wkt": crs_wkt},
        "pixel_size_m": pixel_size,
        "extent": extent,
        "n_rows": n_rows,
        "n_cols": n_cols,
        "dtype": "uint8",
        "bands": ["R", "G", "B", "A"],
        "tiles": tiles,
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
    print(f"  manifest.json written ({len(tiles)} tiles)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--postures-npy", type=Path, required=True)
    ap.add_argument("--mosaic-tif", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True,
                    help="target zebra_bearings/ directory in an unwrap-data checkout")
    ap.add_argument("--pad-m", type=float, default=15.0)
    ap.add_argument("--tile-rows", type=int, default=2400)
    args = ap.parse_args()

    print(f"postures: {args.postures_npy} -> {args.out_dir / 'postures'}")
    npy_to_parquets(args.postures_npy, args.out_dir / "postures")
    print(f"mosaic: {args.mosaic_tif} -> {args.out_dir / 'mosaic'}")
    crop_mosaic_to_tiles(args.mosaic_tif, args.out_dir / "mosaic",
                         pad_m=args.pad_m, tile_rows=args.tile_rows)


if __name__ == "__main__":
    main()
