"""Tests for the zebra escape-map figure path.

Everything runs on synthetic fixtures: a fake ``zebra_bearings`` data
directory (tiny mosaic tiles + posture parquets) and hand-built
trajectories, so no ``unwrap-data`` checkout is needed. Covers the two
map loaders in ``unwrap.datasets`` (tile stitching, window slicing,
downsampling; posture split loading), the distance-mode arrow sampler in
``zebra_plots``, and a smoke render of the map panel.
"""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from PIL import Image

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "notebooks"))

from unwrap.datasets import load_zebra_mosaic, load_zebra_postures  # noqa: E402
import zebra_plots as plots  # noqa: E402


# ---------------------------------------------------------------------------
# synthetic zebra_bearings tree
# ---------------------------------------------------------------------------

# mosaic: 10 cols x 4 rows px, pixel size 1 m, extent (0,10) x (0,4);
# pixel value encodes its (row, col) so slicing can be verified exactly.
N_ROWS, N_COLS, PIXEL_SIZE = 4, 10, 1.0
EXTENT = {"left": 0.0, "right": 10.0, "bottom": 0.0, "top": 4.0}
IMG = np.zeros((N_ROWS, N_COLS, 4), dtype=np.uint8)


def _fill():
    for r in range(N_ROWS):
        for c in range(N_COLS):
            IMG[r, c] = (r, c, (r * 10 + c) % 256, 255)


_fill()


@pytest.fixture()
def data_tree(tmp_path):
    zebra = tmp_path / "zebra_bearings"
    mos = zebra / "mosaic"
    mos.mkdir(parents=True)
    tiles = []
    for i, (r0, nrows) in enumerate([(0, 2), (2, 2)]):
        name = f"tile_{i:03d}.png"
        Image.fromarray(IMG[r0:r0 + nrows], mode="RGBA").save(mos / name)
        tiles.append({"file": name, "row0": r0, "nrows": nrows, "bytes": 0})
    (mos / "manifest.json").write_text(json.dumps({
        "format": "zebra-mosaic-tiles-v1",
        "crs": {"epsg": 32637, "wkt": ""},
        "pixel_size_m": PIXEL_SIZE,
        "extent": EXTENT,
        "n_rows": N_ROWS,
        "n_cols": N_COLS,
        "dtype": "uint8",
        "bands": ["R", "G", "B", "A"],
        "tiles": tiles,
    }))

    pos_dir = zebra / "postures"
    pos_dir.mkdir()
    frames = np.arange(10)
    for t in range(3):
        df = pd.DataFrame({
            "frame": frames,
            "track_number": f"observationT-{t}",
            "kp0_x": 1.0 * t + frames,
            "kp0_y": np.full(10, 2.0),
            "kp1_x": np.full(10, np.nan),
            "kp1_y": np.full(10, np.nan),
        })
        df.to_parquet(pos_dir / f"observationT-{t}.parquet", index=False)
    return zebra


# ---------------------------------------------------------------------------
# load_zebra_mosaic
# ---------------------------------------------------------------------------

def test_mosaic_full_stitch(data_tree):
    img, extent = load_zebra_mosaic(data_tree)
    assert img.shape == (N_ROWS, N_COLS, 4)
    assert img.dtype == np.uint8
    assert extent == (EXTENT["left"], EXTENT["right"], EXTENT["bottom"], EXTENT["top"])
    assert np.array_equal(img, IMG)


def test_mosaic_window_slices_exact_rows_and_cols(data_tree):
    # columns 2..5, rows 1..2 (y 1..3 with top=4, pixel size 1)
    img, extent = load_zebra_mosaic(data_tree, xlim=(2.0, 6.0), ylim=(1.0, 3.0))
    assert img.shape == (2, 4, 4)
    assert np.array_equal(img, IMG[1:3, 2:6])
    assert extent == (2.0, 6.0, 1.0, 3.0)


def test_mosaic_window_never_reads_outside_tiles(data_tree):
    # y 0..1.5 is pixel-aligned to rows 2..3 (row 2 covers y 1..2); extent
    # reports the pixel-aligned window, which may slightly exceed the request
    img, extent = load_zebra_mosaic(data_tree, ylim=(0.0, 1.5))
    assert np.array_equal(img, IMG[2:4])
    assert extent == (0.0, 10.0, 0.0, 2.0)


def test_mosaic_downsample_keeps_extent(data_tree):
    img, extent = load_zebra_mosaic(data_tree, max_side_px=5)
    assert max(img.shape[:2]) <= 5
    assert extent == (EXTENT["left"], EXTENT["right"], EXTENT["bottom"], EXTENT["top"])


def test_mosaic_window_outside_crop_errors(data_tree):
    with pytest.raises(ValueError, match="outside"):
        load_zebra_mosaic(data_tree, xlim=(200.0, 300.0))


# ---------------------------------------------------------------------------
# load_zebra_postures
# ---------------------------------------------------------------------------

def test_postures_tracks_filter_and_nan_preserved(data_tree):
    df = load_zebra_postures(data_tree, observation="observationT", tracks=[0, 2])
    assert set(df["track_number"].unique()) == {"observationT-0", "observationT-2"}
    assert len(df) == 20
    assert df["kp1_x"].isna().all()
    got = df[df["track_number"] == "observationT-2"]
    assert np.allclose(got["kp0_x"], 2.0 + np.arange(10))


def test_postures_all_tracks(data_tree):
    df = load_zebra_postures(data_tree, observation="observationT")
    assert df["track_number"].nunique() == 3


def test_postures_missing_observation_errors(data_tree):
    with pytest.raises(FileNotFoundError):
        load_zebra_postures(data_tree, observation="observationZZZ")


# ---------------------------------------------------------------------------
# frames_every_distance_greedy
# ---------------------------------------------------------------------------

def test_greedy_sampler_fixed_spacing():
    # start_offset 0 places the first arrow at the frame right after the
    # anchor, then every spacing_m from the last anchor (published behavior)
    x = np.arange(20, dtype=float)
    y = np.zeros_like(x)
    idx = plots.frames_every_distance_greedy(x, y, spacing_m=3.0)
    assert list(idx) == [1, 4, 7, 10, 13, 16, 19]


def test_greedy_sampler_start_offset():
    x = np.arange(20, dtype=float)
    y = np.zeros_like(x)
    idx = plots.frames_every_distance_greedy(x, y, spacing_m=3.0, start_offset_m=10.0)
    assert list(idx) == [10, 13, 16, 19]


def test_greedy_sampler_gap_reanchors():
    # 1 m apart, NaN gap on 5..7; after the gap the walk re-anchors at 8
    # and the next placement lands at 11 (8 + 3 m)
    x = np.arange(12, dtype=float)
    y = np.zeros_like(x)
    x[5:8] = np.nan
    y[5:8] = np.nan
    idx = plots.frames_every_distance_greedy(x, y, spacing_m=3.0)
    assert list(idx) == [1, 4, 11]


def test_greedy_sampler_all_nan_empty():
    x = np.full(10, np.nan)
    y = np.full(10, np.nan)
    assert len(plots.frames_every_distance_greedy(x, y, spacing_m=1.0)) == 0


def test_greedy_sampler_bad_spacing():
    with pytest.raises(ValueError):
        plots.frames_every_distance_greedy([1.0, 2.0], [0.0, 0.0], spacing_m=0.0)


# ---------------------------------------------------------------------------
def test_decimate_points_endpoints_and_spacing():
    rng = np.random.default_rng(7)
    x = np.cumsum(rng.normal(0.05, 0.1, 2000))
    y = np.cumsum(rng.normal(0.0, 0.1, 2000))
    sel = plots._decimate_points(x, y, spacing_m=1.0)
    assert sel[0] == 0 and sel[-1] == 1999
    assert np.all(np.diff(sel) > 0)
    pts = np.column_stack([x[sel], y[sel]])
    d = np.hypot(np.diff(pts[:, 0]), np.diff(pts[:, 1]))
    # every greedy-placed vertex is at least spacing_m from the previous
    # kept vertex (the final endpoint is appended unconditionally)
    assert np.all(d[:-1] >= 1.0 - 1e-9)


def test_decimate_points_nan_runs_select_nothing():
    x = np.arange(100, dtype=float)
    y = np.zeros(100)
    x[40:60] = np.nan
    sel = plots._decimate_points(x, y, 10.0)
    assert not np.isin(sel, np.arange(40, 60)).any()
    assert sel[0] == 0 and sel[-1] == 99


def test_decimate_points_two_valid_points():
    x = np.array([1.0, np.nan, np.nan, 4.0])
    y = np.array([0.0, np.nan, np.nan, 0.0])
    sel = plots._decimate_points(x, y, 10.0)
    assert list(sel) == [0, 3]


# plot_escape_map smoke
# ---------------------------------------------------------------------------

def _map_inputs():
    img, extent = IMG, (EXTENT["left"], EXTENT["right"], EXTENT["bottom"], EXTENT["top"])
    pos = pd.DataFrame({
        "frame": np.tile(np.arange(10), 3),
        "track_number": np.repeat([f"observationT-{t}" for t in range(3)], 10),
        "kp0_x": np.concatenate([t + np.arange(10) for t in range(3)]).astype(float),
        "kp0_y": np.full(30, 2.0),
    })
    bearings = pd.DataFrame({
        "track_number": np.repeat("observationT-0", 10),
        "frame": np.arange(10),
        "compass_bearing": np.linspace(0, 90, 10),
        "pre_scare": np.arange(10) < 5,
    })
    return img, extent, pos, bearings


def test_plot_escape_map_distance_mode_smoke():
    img, extent, pos, bearings = _map_inputs()
    fig = plots.plot_escape_map(
        img, extent, pos, bearings, focal_track="observationT-0",
        arrow_mode="distance", arrows_every_m=3.0, arrow_length=0.5,
        save=False,
    )
    ax = fig.axes[0]
    from matplotlib.collections import LineCollection
    from matplotlib.patches import FancyArrowPatch
    lc = [c for c in ax.collections if isinstance(c, LineCollection)]
    arrows = [p for p in ax.patches if isinstance(p, FancyArrowPatch)]
    # focal + non-focal LineCollections; arrows at frames 1,4,7 (sampler
    # candidates [1,4,7,10,...] intersected with the bearing frames 0..9)
    assert len(lc) >= 2
    assert len([a for a in arrows if a.get_zorder() == 4.0]) == 3
    plt.close(fig)


def test_draw_escape_map_distance_uses_frame_ids():
    img, extent, pos, bearings = _map_inputs()
    pos["frame"] = 100 + 2 * pos["frame"]
    bearings["frame"] = 100 + 2 * bearings["frame"]
    fig, ax = plt.subplots()
    try:
        plots.draw_escape_map(
            ax, img, extent, pos, bearings, focal_track="observationT-0",
            arrows_every_m=3.0, frame_max=108,
        )
        assert len(ax.patches) == 2  # focal rows 1 and 4, frames 102 and 108
    finally:
        plt.close(fig)


def test_plot_escape_map_time_mode_smoke():
    img, extent, pos, bearings = _map_inputs()
    fig = plots.plot_escape_map(
        img, extent, pos, bearings, focal_track="observationT-0",
        arrow_mode="time", fps=1.0, arrow_every_s=1.0, arrow_length=0.5,
        save=False,
    )
    from matplotlib.patches import FancyArrowPatch
    arrows = [p for p in fig.axes[0].patches
              if isinstance(p, FancyArrowPatch) and p.get_zorder() == 4.0]
    assert len(arrows) == 10  # every frame at 1 fps / 1 s spacing
    plt.close(fig)


def test_plot_escape_map_bad_mode():
    img, extent, pos, bearings = _map_inputs()
    with pytest.raises(ValueError, match="arrow_mode"):
        plots.plot_escape_map(img, extent, pos, bearings,
                              focal_track="observationT-0",
                              arrow_mode="sideways", save=False)
