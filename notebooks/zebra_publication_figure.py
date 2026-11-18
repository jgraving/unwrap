"""Reproduce the zebra escape publication figure from the saved fit.

This rebuilds the circstats paper's zebra composite figure (900 x 1080
points, committed as ``assets/reference.pdf``) from stored data only: the
aerial mosaic, the posture trajectories, the per-frame bearings and the
cached hierarchical von Mises fit, all pinned to a ``jgraving/unwrap-data``
revision with SHA-256 checksums. Nothing is sampled, refit, or downloaded
outside that revision.

Layout: the published composite's own page-point geometry (measured from
the reference PDF; see ``ZEBRA_FIGURE.md`` for the measurement notes, the
inset anchor frames, the asset provenance, and the fitted-number
baseline):

* a  the escape map (mosaic, trajectories, distance-placed orientation
     arrows) with the two circular callout insets and their body-heading
     arrows, drawn by the ``draw_escape_map`` stack from ``zebra_plots``
     (zebra map PR #154);
* b  the species photograph (provided raster asset);
* c  the equally weighted mean direction of the focal group (polar KDE);
* d  group-level alignment (resultant length of the 18 fitted individual
     mean directions);
* e  the 18 individual mean directions (per-track polar KDEs plus the
     observed circular-mean headings);
* f  individual-level concentration (paired means with 95% HDIs).

c-f are evaluated by :mod:`zebra_posthoc.summary_stats` from the cached
posterior; the pre/post markers and labels follow the house convention
(orange square = before the approach, purple circle = after).

Exports: SVG (fonts outlined), PDF (the original JPEG sources restored to
their image objects), a PNG preview rendered from the PDF, self-contained
HTML with the caption, and a provenance JSON recording fonts, package
versions, input hashes and the numerical summaries.

Run from the repository root; see ``ZEBRA_FIGURE.md`` for provenance.
"""
from __future__ import annotations

import argparse
import base64
import hashlib
import importlib.metadata
import io
import json
from pathlib import Path
import platform
import re
import subprocess
import urllib.request
import xml.etree.ElementTree as ET

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import (ArrowStyle, Circle, FancyArrowPatch,
                                FancyBboxPatch, Rectangle)
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ASSETS = Path(__file__).parent / "assets" / "zebra_publication"

#: jgraving/unwrap-data revision carrying the map data (PR #11) and the
#: zebra fit cache. The zebra fit blobs in this revision are byte-identical
#: to those of revision 31eb4bdb4a0fe9a0c0bd20df7f101dbd884a8dcc (the
#: beetle/locust figures' data pin).
DATA_REVISION = "86c41b5f43eb61191405b9bc85baf730a492a06f"

INPUT_HASHES = {
    'fit_cache/zebra_escape_fit.meta.json': '324d49f261c5deebdcbe0c6dfb3607a4bc34a68183dbf8856929b8390652e23a',
    'fit_cache/zebra_escape_fit.nc': '0c01dd8463651a817f90776e1089928e84b05a594bf1f86ad7e08bbdfe7ec71b',
    'zebra_bearings/mosaic/manifest.json': 'b705d1ea3df720ca212f510da9d00333b13c125c0831a5e74c9295d7e94e0bb3',
    'zebra_bearings/mosaic/tile_000.png': 'bebd6a8ce262377a130c90dc4f853e0109275ee6e07fea4ffdadc3ce3ef09bb8',
    'zebra_bearings/mosaic/tile_001.png': 'd49774d1977608374b97ee123573705155c8963e4a1dc338937d1046ba55705e',
    'zebra_bearings/mosaic/tile_002.png': '1798bae732618d9b8d3ef856ad3148cae88581d674914ae9b6a4dd1a33cfc0a2',
    'zebra_bearings/mosaic/tile_003.png': 'eb4729f4dbc851fd592a2e42e24f2112b307d0f23cde738794a21097dc3a3327',
    'zebra_bearings/mosaic/tile_004.png': '8d30110bad0939dc1938052477f632bcd78ccd94a23a7f79640f37f3e39a0c45',
    'zebra_bearings/mosaic/tile_005.png': 'ae3f1fc72c0e0019c5dc424f583386b86f7d9e48d97eae71b58ce8e9959100dc',
    'zebra_bearings/mosaic/tile_006.png': 'c471cc957edc439ba17a3288007e7bf10bcb4f91c92f6d6915b5e4fb03ccb507',
    'zebra_bearings/mosaic/tile_007.png': '9099290a363c3c379fb95a5a39b75302bf19a5a0b9bdb65a3f7c95c56f79ceeb',
    'zebra_bearings/postures/observation088-0.parquet': '2c57928ab2d74c8b044ffcf3098488cd2a05e3188878c18175b9e79b7b28881d',
    'zebra_bearings/postures/observation088-1.parquet': 'c118f5f35d7a71800c4bfc0bc29d86561e9f45c27b7f40ffc9fce10fc004af49',
    'zebra_bearings/postures/observation088-10.parquet': '2a28a61d208e83a5d37081763510d9bffe660d437c6e29fa9d162dbaded46cd8',
    'zebra_bearings/postures/observation088-11.parquet': 'd1a11113260e2302fcb35f549dfd2462d16c989d67478601f167bb28ff9b9ba3',
    'zebra_bearings/postures/observation088-12.parquet': 'e49500ef46b085de061d2a59e3e79502866de6bf9b80f67bec6dcc473b9bafc5',
    'zebra_bearings/postures/observation088-13.parquet': '6c9d9560e76ddfab2dd2bb7a2b18b9dee37c7157323a260a1c79e1b62af11b04',
    'zebra_bearings/postures/observation088-14.parquet': '8a3b8a7f02f6284f9621c892e99292ede5a167da39b984ee6b6d0fc255e32a58',
    'zebra_bearings/postures/observation088-15.parquet': '2bd1946e4e535d0f1c4e534a8c7be3d6b38c1359d97e3a4ac5011bb9757113a0',
    'zebra_bearings/postures/observation088-16.parquet': '1f4827307b4ba14733ac5eba7af424d426e62a06f4b1a8c8ad7a87d56e30760b',
    'zebra_bearings/postures/observation088-17.parquet': '774c47eb96726d0b57190c362476e69b745a4642e987be9a28e9388a2ca8d484',
    'zebra_bearings/postures/observation088-2.parquet': 'ef0a78d89288abcfccbf1073ea72f766556adcbb30f9d657d21e7a9b1f48f5a3',
    'zebra_bearings/postures/observation088-3.parquet': 'da6abd0a0b6c69e473901509e185b48fb44bb38df5b2b4a89588902a19cb2591',
    'zebra_bearings/postures/observation088-4.parquet': 'b4b1533242ab8c6ade1a6e649910799106e5ae693f6323e9d38cf593645f7c0e',
    'zebra_bearings/postures/observation088-5.parquet': '3655f9116b01d83262f7f29cdd83e112c0a104d5409a4ea44d9b08745e26487f',
    'zebra_bearings/postures/observation088-6.parquet': '6f7d6653b6bd5824c9eee19db77e98c2c0e706f1937ea9503b838f02da851aac',
    'zebra_bearings/postures/observation088-7.parquet': 'e213126655a389e9910ad9c7c5878dedff961cafeeb82477f562d237f2a075fe',
    'zebra_bearings/postures/observation088-8.parquet': '098bb98ad180678c00b23f9643831c96a5e42985e8ec848663cccf79bfc693f9',
    'zebra_bearings/postures/observation088-9.parquet': '8b3d2588ab80e2139fd47ef78f419428f35e9bbe20cd348fd36ee85abc36739b',
    'zebra_bearings/raw/observation088-0.parquet': '852ccd30371ae568e5f25d9688b946c912ed35e62d0f1b3058c3b16d223e017e',
    'zebra_bearings/raw/observation088-1.parquet': '7dbcafada56e5a115f7faaed3fc7bac4df4571fb6f4ee0475b46168858f7b7cf',
    'zebra_bearings/raw/observation088-10.parquet': 'a38fa10ae36f52087cc633e5d448fd925e4958a7c5610e8cfdcb93758d15952b',
    'zebra_bearings/raw/observation088-11.parquet': '5d4a20700eb3bff4a79ca54133cd2d2fd26532e16b64b074d6dfd06997a8d9a4',
    'zebra_bearings/raw/observation088-12.parquet': '8fdbff49f1abf9ee6219795d2ef30303a8cb718a40a7d1c317ca7558c8e49899',
    'zebra_bearings/raw/observation088-13.parquet': 'f0e3cb808f6c7921cd95dcbfd4fd655896a06dfe56783e081c6e2ace3a3f53eb',
    'zebra_bearings/raw/observation088-14.parquet': '6405d60e6037c2c66460586df41be54a563e30c3a5ceb240f2894fa1ac74a55c',
    'zebra_bearings/raw/observation088-15.parquet': '25bbadfd189c77e8dc0b1b8ef8b468a115eb3f9e86171858651a8b6f52860ca4',
    'zebra_bearings/raw/observation088-16.parquet': '5e54e45feace1d3b090605136a8c1782b0a520a1b4cf44e140db78daf2059538',
    'zebra_bearings/raw/observation088-17.parquet': '8543418cd2ca19e11aa1080b0f123f2620623c02190741230dde35cc8fc798b2',
    'zebra_bearings/raw/observation088-2.parquet': 'c4206d9b977aeb761d9c99dbf193cd6691ff46ec23f56ffb4b11a83060a19271',
    'zebra_bearings/raw/observation088-3.parquet': 'f1c22717c6168f11d925d2bb7a012b79ad02bdb3bf70fa78ef5f8e6bf42ba631',
    'zebra_bearings/raw/observation088-4.parquet': '989d098018db81079c3f4680416ea63a77cd0e471ca6ef072db64793dfba7af8',
    'zebra_bearings/raw/observation088-5.parquet': 'cf7c327ac064f8c331fa7002ada51ecdcebfaafb83a75fe5b1f43e0bb09cebe6',
    'zebra_bearings/raw/observation088-6.parquet': '58b07d05619eee249371d255aba6016057dd78afd66d85302c5587b92ba78f7b',
    'zebra_bearings/raw/observation088-7.parquet': 'f22b91eb524142b6b2d54c8d84b2b713f4edb19ce46db282957429444198eefa',
    'zebra_bearings/raw/observation088-8.parquet': '6a0175fd6e7151830c888eeb38cbf9414df728be4368dcb6c09f36771136fbd2',
    'zebra_bearings/raw/observation088-9.parquet': '9458cac455479de27905d650b7da76e0383c6dcce11e518d30358e890006c006',
    'zebra_bearings/zebra_bearings_blocks.parquet': '51beb1e2b92fd2d71e4620d7c4e2cbdf2c1b62ec5ecafbac86b8e205847a7e6d',
}

# ---------------------------------------------------------------------------
# Published map window (map section of notebooks/zebra_escape.ipynb)
# ---------------------------------------------------------------------------
OBS = "observation088"
FOCAL_TRACK = f"{OBS}-5"
FPS = 30.0
FRAME_MAX = 90698                 # last usable tracking frame
PANEL_XLIM = (267525.0, 267680.0)    # UTM 37N meters
PANEL_YLIM = (37300.0, 37770.0)
MOSAIC_YLIM = (37290.0, 37770.0)     # mosaic read window (10 m clip below)
MOSAIC_MAX_SIDE_PX = 2600       # downsampled mosaic side (house value)
#: The published composite embeds the mosaic as a 640 x 1940 px raster.
# Limit the input raster before clipping it to the displayed map window.
MOSAIC_EMBED_SIDE = 1940         # long side of the embedded mosaic, px
FOCAL_LW = 2.5                 # focal trajectory width (measured on the paper)
#: Trajectory vertex spacing (meters) for the composite's map panel:
# one meter is 2.26 points at the panel's 0.443 m/pt, so the
# distance-decimated polylines stay subpixel-faithful at panel
# scale. The 18 tracks cover ~48.4 km of travel; at this spacing the
# vector export keeps ~48k line segments, the size at which the SVG
# backend (one element per collection segment) stays in the house
# artwork range. The published composite flattens its own
# trajectories into a ~640-px-wide raster with no finer detail.
MAP_VERTEX_SPACING_M = 1.0

# ---------------------------------------------------------------------------
# Layout of the published composite, in page points, origin top-left
# (measured from the reference PDF; see ZEBRA_FIGURE.md).
# ---------------------------------------------------------------------------
PAGE_W, PAGE_H = 900.0, 1080.0
MAP_RECT = (15.9, 12.3, 365.8, 1073.0)
MAP_FRAME_RECT = (13.7, 10.1, 367.9, 1075.2)
PHOTO_RECT = (393.4, 12.3, 877.9, 335.3)
GRID_BACKING_RECT = (376.0, 370.7, 872.2, 1074.4)

POLAR_C = (503.6, 544.45)     # panel c polar axis centre
POLAR_E = (503.6, 902.3)      # panel e polar axis centre
POLAR_R_EDGE_PT = 82.7        # axis edge (r = 2.5 in KDE units)
POLAR_RMAX = 2.5
POLAR_RTICKS = (1.0, 2.0)     # grid circles at KDE density 0 and 1
POLAR_GRID_COLOR = "0.51"
POLAR_GRID_LW = 1.375
R_DOT = 2.25                  # radius carrying the mean-direction dots

BAR_D = (702.0, 439.2, 867.4, 649.7)
BAR_F = (702.0, 797.1, 867.4, 1007.6)

# degree labels: (text, centre x, centre y) measured from the PDF
DEGREE_LABELS_C = [
    ("0°", 503.6, 436.8), ("45°", 579.8, 466.8),
    ("90°", 611.3, 542.9), ("135°", 579.8, 619.1),
    ("180°", 503.6, 650.6), ("225°", 427.5, 619.1),
    ("270°", 395.9, 542.9), ("315°", 427.5, 466.8)]
DEGREE_LABELS_E = [(t, x, y + (POLAR_E[1] - POLAR_C[1]))
                   for t, x, y in DEGREE_LABELS_C]
COMPASS_LABELS = [("N", 0, -96), ("E", 96, 0),
                  ("S", 0, 96), ("W", -96, 0)]
RADIAL_LABEL_C = [("0", 520.5, 505.8), ("1.0", 539.4, 475.2)]
RADIAL_LABEL_E = [(t, x, y + (POLAR_E[1] - POLAR_C[1]))
                  for t, x, y in RADIAL_LABEL_C]

# titles: (text, centre x, centre y), 17.5 pt regular, centred
TITLES = [
    ("Group-level\nmean direction", 503.6, 388.4),
    ("Group-level alignment", 784.7, 383.1),
    ("Individual-level\nmean direction", 503.5, 746.3),
    ("Individual-level\nconcentration", 784.7, 744.3)]

# panel letters: (char, left x, centre y, color), 32 pt bold
LETTERS = [
    ("a", 31.8, 38.0, "white"),
    ("b", 408.8, 38.0, "white"),
    ("c", 390.6, 380.2, "black"),
    ("d", 649.1, 380.2, "black"),
    ("e", 390.6, 743.2, "black"),
    ("f", 653.6, 743.2, "black")]

GRID_COLOR = "black"
TICK_LW = 1.5

# Page coordinates (pt, y down). The full photographs extend beyond
# their circular masks. Rotations and arrow endpoints are from the artwork.
INSET_PRE = dict(center=(115.04, 331.45), r=65.25, side=208.7775, rotation=-90,
                 arrow=((81.36, 316.29), (129.57, 342.03)), anchor_frame=17623,
                 asset="drone_inset_pre.jpg")
INSET_POST = dict(center=(115.04, 743.02), r=67.16, side=214.9035, rotation=-45,
                  arrow=((149.63, 713.68), (108.58, 750.25)), anchor_frame=80780,
                  asset="drone_inset_post.jpg")
RING_LW = 5.0
CONE_LW = 3.0
ANCHOR_WINDOW_FRAMES = 15       # +/- 0.5 s at 30 fps

# Map legend centre x=103.35 pt in the 349.9-pt map panel.
MAP_LEGEND_ANCHOR = (0.25, 0.5)
# The reference scales the original map's point-based arrow styling by
# 1.519. The heading arrows' nominal length is still 12 UTM meters.
MAP_ARROW_SCALE = 1.5188

# scale bar and north arrow, measured
SCALE_RECTS = ((296.5, 1014.2, 22.5, 8.7), (319.0, 1014.2, 22.6, 8.7))
SCALE_RECT_COLORS = ("white", "black")
SCALE_LABEL = (319.05, 1039.35)
NORTH_X = 319.05
NORTH_Y = 911.4                 # "N" text centre
NORTH_TIP_Y = 930.48
NORTH_BASE_Y = 988.19
NORTH_HEAD_LEN = 21.27
NORTH_HEAD_WIDTH = 18.23
NORTH_TAIL_WIDTH = 3.80

# bottom Pre/Post legend (box, lines, markers, texts), measured
LEGEND_BOX = (542.8, 1042.6, 731.1, 1072.6)
LEGEND_EDGE = (0.8, 0.8, 0.8)
LEGEND_PRE_LINE = (549.8, 1056.0, 584.8, 1056.0)
LEGEND_POST_LINE = (643.6, 1056.0, 678.6, 1056.0)
LEGEND_PRE_MARK = (562.9, 1051.6, 8.8, 8.8)      # square (x, y, w, h)
LEGEND_POST_MARK = (656.7, 1051.6, 8.7, 8.7)     # circle (x, y, d, d)
LEGEND_TEXTS = [("Pre", 608.95, 1056.1), ("Post", 706.6, 1056.1)]

# typography, measured from the PDF
FS = dict(letter=32.0, title=17.5, degree=17.5, radial=15.0,
          tick=17.5, axis=17.5, bottom_legend=17.5, scale=15.2,
          north=18.2, map_legend=15.2)
COMPANION_FS = dict(letter=22.0, title=15.0, degree=13.0, radial=12.0,
                    tick=14.0, axis=14.0, bottom_legend=14.0, scale=13.0,
                    north=15.0, map_legend=13.0)

# marker sizes (points) and line widths, measured
MS_GROUP_PRE = 9.9          # c/d pre square side, pt
MS_GROUP_POST = 8.8         # c/d post circle: s = pi/4 d^2, d = 9.9 pt
MS_EMP_PRE = 7.6            # e pre square side, pt
MS_EMP_POST = 6.6           # e post circle, d = 7.4 pt
MS_CONC_PRE = 9.1           # f pre square side, pt
MS_CONC_POST = 8.1          # f post circle, d = 9.1 pt
LW_GROUP_LOOP = 4.0
LW_IND_LOOP = 0.8
IND_LOOP_ALPHA = 0.5        # house value; the paper's baseline ring is
#                            # the salmon blend of 36 overlapped loops
LW_HDI = 2.0
LW_PAIRED = 1.4

JITTER_SEED = 123               # per-track paired-offset seed (house value)
PREVIEW_SCALE = 2               # 1800 x 2160 px at the published page size
JITTER_SD = 0.15

#: The figure caption (fresh wording; the older manuscript caption called
#: panel c a "population mean", which the finite focal-group estimand is
#: not).
MODEL_ANGLE_NOTE = ("Panels c and e show model angles clockwise from 0° "
                    "at the top. Panel a uses geographic compass bearings.")
COMPASS_ANGLE_NOTE = "Panels a, c and e use geographic compass bearings."
CAPTION = (
    "Escape response of a focal group of 18 Grevy's zebra (*Equus grevyi*) "
    "filmed from overhead drones at the Mpala Conservancy, Kenya. "
    "(a) Georeferenced trajectories of the 18 tracked zebra on the aerial "
    "orthomosaic of the grazing site: the focal animal of this panel in "
    "black, the rest of the group in white, with the focal animal's body "
    "orientation every 5 m of straight-line travel as arrows (orange before "
    "the approach stimulus, purple after). The circular insets (orange / "
    "purple) magnify the moments indicated by each callout and show the "
    "body-heading arrow drawn at that moment. (b) A Grevy's zebra. "
    "(c) Posterior distribution of the equally weighted mean direction of "
    "the 18 focal zebras before and after the approach: the loop traces the "
    "relative density of the group-mean-direction posterior draws (radial "
    "scale 0 to 1) and the marker the posterior mean direction. "
    "(d) Group-level alignment: the resultant length of the 18 fitted "
    "individual mean directions, with the posterior mean before and after "
    "the approach. (e) Individual-level mean directions: the posterior "
    "density of each fitted individual mean direction (loops) and each "
    "animal's observed circular-mean heading (markers). " + MODEL_ANGLE_NOTE + " "
    "(f) Individual-level "
    "concentration: each animal's posterior mean resultant length with 95% "
    "HDI, before and after the approach, paired within animal across "
    "conditions. Panels c-f are evaluated from a saved hierarchical von "
    "Mises fit of the block-averaged body orientations (4 chains x 1,000 "
    "draws, zero divergent transitions); the 'group' in c-d is the finite "
    "focal group, not a population parameter."
)
COMPANION_CAPTION = CAPTION.replace(MODEL_ANGLE_NOTE, COMPASS_ANGLE_NOTE)

SCARE_FRAME_NOTE = (
    "The per-frame pre_scare flag (used for the map's pre/post arrow "
    "splitting, as in PR #154) flips at frame 53902 for all 18 "
    "observation088 tracks. scare_frames.csv lists observation088 at frame "
    "72703; the figure follows the per-frame flag and the discrepancy is "
    "recorded as published."
)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare_inputs(root, *, download=False):
    """Check the pinned unwrap-data inputs against :data:`INPUT_HASHES`.

    Missing files are fetched from ``jgraving/unwrap-data`` at
    :data:`DATA_REVISION` when ``download`` is set; existing files are
    verified in place and never overwritten.
    """
    root = Path(root)
    for relative, expected in INPUT_HASHES.items():
        path = root / relative
        if not path.exists():
            if not download:
                raise FileNotFoundError(
                    f"Missing {path}; supply --download or --data-root")
            url = (f"https://raw.githubusercontent.com/jgraving/"
                   f"unwrap-data/{DATA_REVISION}/{relative}")
            with urllib.request.urlopen(url, timeout=60) as response:
                data = response.read()
            if hashlib.sha256(data).hexdigest() != expected:
                raise ValueError(f"Downloaded checksum mismatch: {relative}")
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
        if sha256(path) != expected:
            raise ValueError(
                f"Checksum mismatch: {path}; existing files are never "
                f"overwritten")
    return root


def resolve_fonts(directory=None, family="Arial"):
    """Register ``directory`` fonts and record the resolved family faces.

    A missing face raises instead of silently substituting; the resolved
    font files and their hashes land in the provenance.
    """
    if directory:
        for path in font_manager.findSystemFonts(fontpaths=[str(directory)]):
            font_manager.fontManager.addfont(path)
    fonts = {}
    for weight in ("normal", "bold"):
        try:
            path = font_manager.findfont(
                font_manager.FontProperties(family=family, weight=weight),
                fallback_to_default=False)
        except ValueError as exc:
            raise ValueError(
                f"{family} unavailable; supply licensed fonts with "
                f"--font-dir") from exc
        fonts[weight] = {"file": Path(path).name, "sha256": sha256(path)}
    return fonts


def check_assets():
    """Verify the committed raster assets against their manifest."""
    manifest = json.loads((ASSETS / "images.json").read_text())
    for name, item in manifest["images"].items():
        path = ASSETS / name
        if not path.is_file():
            raise FileNotFoundError(f"Missing asset {path}")
        if sha256(path) != item["sha256"]:
            raise ValueError(f"Asset checksum mismatch: {path}")
    if sha256(ASSETS / "reference.pdf") != manifest["source_sha256"]:
        raise ValueError("Reference PDF checksum mismatch")
    return manifest


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
def load_fit(root):
    """Load the cached fit, reject divergent caches, read the fit metadata."""
    import arviz as az

    fit_path = root / "fit_cache" / "zebra_escape_fit.nc"
    idata = az.from_netcdf(fit_path)
    if "diverging" in idata.sample_stats:
        divergences = int(idata.sample_stats.diverging.values.sum())
    else:
        divergences = 0
    if divergences:
        raise ValueError(
            f"Cached fit has {divergences} divergences; review before plotting")
    metadata = json.loads((root / "fit_cache"
                           / "zebra_escape_fit.meta.json").read_text())
    return idata, metadata, divergences


def fit_diagnostics(idata):
    """Independent ArviZ reads of the cached fit's convergence statistics.

    These are checks of the existing fit (recorded for provenance), not a
    new sample. Plain split R-hat, rank-normalized R-hat, and bulk/tail
    ESS are each reported over every posterior variable and level.
    """
    import arviz as az

    posterior = idata.posterior

    def _to_dataset(obj):
        # arviz 1.x returns diagnostics as an xarray DataTree rooted at the
        # posterior (one level per variable and random-effect level)
        if hasattr(obj, "ds"):
            return obj.ds
        return obj

    return {
        "max_rhat_split": float(_to_dataset(az.rhat(posterior))
                                .to_array().max()),
        "max_rhat_rank": float(_to_dataset(az.rhat(posterior, method="rank"))
                               .to_array().max()),
        "min_ess_bulk": float(_to_dataset(az.ess(posterior, method="bulk"))
                              .to_array().min()),
        "min_ess_tail": float(_to_dataset(az.ess(posterior, method="tail"))
                              .to_array().min()),
        "arviz": importlib.metadata.version("arviz"),
    }


def summary_numbers(idata, df_blocks):
    """The numerical baseline behind panels c-f, for the provenance.

    Uses :func:`zebra_posthoc.summary_stats` (the same path the notebook
    figure uses), plus the observed-heading counts and the paired
    per-animal turn (circular mean of the per-track post-minus-pre
    observed-mean directions). Directions are reported in degrees
    on the figure's polar axis (zero at north, clockwise) -- the
    scale the loops and markers are drawn on -- so the JSON and
    the picture agree.
    """
    import zebra_posthoc as posthog
    from unwrap.circular import circ_mean as _circ_mean

    prep = posthog.summary_stats(idata, df_blocks)

    def axis_deg(angle_rad):
        # degrees on the figure's polar axis (zero at north,
        # clockwise -- the data angle as plotted on that scale)
        return float(np.rad2deg(angle_rad) % 360.0)

    diffs = []
    for m_pre, m_post in zip(prep["emp_mu_pre"], prep["emp_mu_post"]):
        if np.isnan(m_pre) or np.isnan(m_post):
            continue
        # wrap the per-track turn onto the circle before averaging
        diffs.append(float(np.angle(np.exp(1j * (m_post - m_pre)))))
    paired_turn = (axis_deg(_circ_mean(np.asarray(diffs)))
                   if diffs else None)

    return {
        "stats": prep,
        "group_mean_direction_deg": {
            "pre": axis_deg(prep["mu_group_pre_mean"]),
            "post": axis_deg(prep["mu_group_post_mean"]),
        },
        "group_alignment": {
            key: float(prep[key]) for key in
            ("R_group_pre_mean", "R_group_pre_lo", "R_group_pre_hi",
             "R_group_post_mean", "R_group_post_lo", "R_group_post_hi")},
        "individual_concentration": {
            "R_mean_pre": [float(v) for v in prep["R_ind_pre_mean"]],
            "R_mean_post": [float(v) for v in prep["R_ind_post_mean"]],
            "R_hdi_pre": [[float(a), float(b)]
                          for a, b in np.asarray(prep["R_ind_pre_hdi"])],
            "R_hdi_post": [[float(a), float(b)]
                           for a, b in np.asarray(prep["R_ind_post_hdi"])],
        },
        "empirical_mean_direction_deg": {
            "pre": [None if np.isnan(v) else axis_deg(v)
                    for v in np.asarray(prep["emp_mu_pre"], dtype=float)],
            "post": [None if np.isnan(v) else axis_deg(v)
                     for v in np.asarray(prep["emp_mu_post"], dtype=float)],
        },
        "observed_headings_n": {
            "pre": int((df_blocks["PostStimulus"] == 0).sum()),
            "post": int((df_blocks["PostStimulus"] == 1).sum()),
        },
        "paired_permanent_turn_deg": paired_turn,
        "tracks": [str(t) for t in np.sort(df_blocks["Track"].unique())],
        "nT": int(prep["nT"]),
    }


def build_context(data_root, *, mosaic_max_side_px=MOSAIC_MAX_SIDE_PX):
    """Load every input and assemble the drawing context.

    Returns a dict with the mosaic array + extent, the posture/bearing
    tables, the summary numbers (including the posterior-draw arrays
    needed to draw c-f), the inset anchor positions and arrow bearings,
    and the raster asset paths.
    """
    import sys

    sys.path.insert(0, str(Path(__file__).parent))
    from unwrap.datasets import (load_zebra_bearings, load_zebra_bearings_raw,
                                 load_zebra_mosaic, load_zebra_postures)

    root = Path(data_root)
    idata, fit_meta, divergences = load_fit(root)
    df_blocks = load_zebra_bearings(data_root=root).copy().reset_index(drop=True)
    numbers = summary_numbers(idata, df_blocks)

    img, extent = load_zebra_mosaic(
        data_root=root, xlim=PANEL_XLIM, ylim=MOSAIC_YLIM,
        max_side_px=mosaic_max_side_px)
    img = _downsample_side(img, MOSAIC_EMBED_SIDE)
    inputs = sorted(Path(name) for name in INPUT_HASHES)
    posture_tracks = [int(path.stem.removeprefix(f"{OBS}-")) for path in inputs
                      if path.parent == Path("zebra_bearings/postures")]
    bearing_tracks = [path.stem for path in inputs
                      if path.parent == Path("zebra_bearings/raw")]
    df_pos = load_zebra_postures(data_root=root, observation=OBS,
                                 tracks=posture_tracks)
    df_bearings = load_zebra_bearings_raw(data_root=root, tracks=bearing_tracks)

    pos = df_pos.loc[df_pos["track_number"] == FOCAL_TRACK] \
                .sort_values("frame")
    bear = df_bearings.loc[df_bearings["track_number"] == FOCAL_TRACK] \
                     .sort_values("frame").set_index("frame")
    anchors = {}
    for key, spec in (("pre", INSET_PRE), ("post", INSET_POST)):
        f0 = int(spec["anchor_frame"])
        if f0 not in set(pos["frame"]):
            raise ValueError(f"Anchor frame {f0} missing from the focal track")
        row = pos.loc[pos["frame"] == f0]
        x, y = float(row["kp0_x"].iloc[0]), float(row["kp0_y"].iloc[0])
        window = bear.loc[
            (bear.index >= f0 - ANCHOR_WINDOW_FRAMES)
            & (bear.index <= f0 + ANCHOR_WINDOW_FRAMES),
            "compass_bearing"].dropna().to_numpy(float)
        if window.size < 10:
            raise ValueError(f"Too few bearings around anchor frame {f0}")
        math = np.deg2rad(90.0 - window)
        mean_bearing = float(
            (90.0 - np.rad2deg(np.arctan2(np.sin(math).mean(),
                                          np.cos(math).mean()))) % 360.0)
        anchors[key] = {"frame": f0, "pos": [x, y],
                        "arrow_bearing_deg": mean_bearing,
                        "n_window": int(window.size)}

    return {
        "idata": idata, "img": img, "extent": tuple(float(v) for v in extent),
        "df_pos": df_pos, "df_bearings": df_bearings,
        "numbers": numbers, "anchors": anchors,
        "photo": ASSETS / "zebra_photo.jpg",
        "inset_assets": {key: ASSETS / spec["asset"]
                         for key, spec in (("pre", INSET_PRE),
                                           ("post", INSET_POST))},
        "fit_meta": fit_meta, "divergences": divergences,
    }


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------
def _downsample_side(img, max_side):
    """Resize ``img`` so its longest side is at most ``max_side`` px.

    The export only needs to match the published composite's mosaic
    raster; everything else is vector. No-op when already small.
    """
    if int(max(img.shape[:2])) <= int(max_side):
        return img
    from PIL import Image

    height, width = img.shape[:2]
    scale = max_side / float(max(height, width))
    resized = Image.fromarray(img).resize(
        (max(1, int(round(width * scale))), max(1, int(round(height * scale)))),
        Image.LANCZOS)
    return np.asarray(resized)


def _figure_rect(rect):
    """PDF page rect (x0, y0, x1, y1, y down) -> matplotlib figure fraction."""
    x0, y0, x1, y1 = rect
    return [x0 / PAGE_W, (PAGE_H - y1) / PAGE_H,
            (x1 - x0) / PAGE_W, (y1 - y0) / PAGE_H]


def _tangent_endpoints(apex, center, radius):
    """The two tangent points from ``apex`` to the circle (same units).

    In the right triangle apex-tangent-center (right angle at the tangent
    point) the angle at ``center`` between the direction to ``apex`` and
    the radius to the tangent point is acos(r/d), so the tangent points
    are ``center + r * R(+/-acos(r/d)) * unit(apex - center)``.
    """
    p = np.asarray(apex, float)
    c = np.asarray(center, float)
    d = float(np.linalg.norm(c - p))
    if d <= radius:
        raise ValueError("Anchor inside the inset circle; cone undefined")
    cos_t = radius / d
    sin_t = float(np.sqrt(1.0 - cos_t * cos_t))
    unit = (p - c) / d
    ends = []
    for sign in (1.0, -1.0):
        rot = np.array([[cos_t, -sign * sin_t],
                        [sign * sin_t, cos_t]])
        ends.append(c + radius * (rot @ unit))
    return ends


def _draw_inset(map_ax, overlay, fig, key, spec, anchor, color, *,
                cone_lw, ring_lw):
    """One circular callout: tangent cone on the map, inset image, ring,
    and body-heading arrow on the overlay plane."""
    cx, cy = spec["center"]
    r = spec["r"]

    # tangent cone, in the map's UTM coordinates (under the inset disc)
    apex = np.asarray(anchor["pos"], float)
    to_data = map_ax.transData.inverted()

    def page_to_data(px, py):
        disp = fig.transFigure.transform((px / PAGE_W, (PAGE_H - py) / PAGE_H))
        return to_data.transform(disp)

    center_data = np.asarray(page_to_data(cx, cy), float)
    radius_data = float(page_to_data(cx + r, cy)[0] - center_data[0])
    for end in _tangent_endpoints(apex, center_data, radius_data):
        map_ax.plot([apex[0], end[0]], [apex[1], end[1]], color=color,
                    lw=cone_lw, solid_capstyle="butt", zorder=6)

    # The JPEG will be drawn beneath this ring, in page coordinates.
    overlay.add_patch(Circle((cx, cy), r, facecolor="none", edgecolor=color,
                             lw=ring_lw, zorder=2, gid=f"inset-{key}"))
    overlay.add_patch(FancyArrowPatch(
        *spec["arrow"], arrowstyle=ArrowStyle.Simple(
            head_length=18.0, head_width=20.0, tail_width=5.0),
        mutation_scale=1.0, shrinkA=0, shrinkB=0,
        facecolor=color, edgecolor="none", zorder=3))


def _style_polar(ax, *, grid_color, grid_lw, compass):
    # The model uses mathematical angles (east = 0, increasing CCW).
    # Plot in those coordinates to display geographic directions.
    ax.set_theta_zero_location("E" if compass else "N")
    ax.set_theta_direction(1 if compass else -1)
    ax.set_rlim(0, POLAR_RMAX)
    ax.set_rticks(POLAR_RTICKS)
    angles = np.deg2rad(np.arange(0.0, 360.0, 90.0 if compass else 45.0))
    ax.set_xticks(angles)
    ax.set_xticklabels([""] * len(angles))  # spokes only; labels are manual
    ax.set_yticklabels([])                # radial labels are manual
    ax.grid(True, color=grid_color, linewidth=grid_lw)
    ax.spines["polar"].set_visible(compass)
    if compass:
        ax.spines["polar"].set_edgecolor(grid_color)
        ax.spines["polar"].set_linewidth(grid_lw)
    ax.tick_params(axis="x", which="major", length=0, pad=0)
    ax.tick_params(axis="y", which="major", length=0, pad=0)


def _style_bar(ax, *, tick_color, tick_lw):
    ax.set_facecolor("white")
    for side, spine in ax.spines.items():
        spine.set_visible(side in ("left", "bottom"))
        spine.set_linewidth(tick_lw)
        spine.set_color(tick_color)
    ax.tick_params(axis="both", which="major", length=5, width=tick_lw,
                   colors=tick_color)


SVG_NS = "http://www.w3.org/2000/svg"
XLINK = "{http://www.w3.org/1999/xlink}href"


def _embed_artwork(svg, ctx, *, caption):
    """Place JPEGs directly in SVG page space, shared by all four exports."""
    from PIL import Image

    ET.register_namespace("", SVG_NS)
    ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")
    ns = "{" + SVG_NS + "}"
    root = ET.fromstring(svg)
    # Older Matplotlib versions hash object addresses into polar clip IDs.
    # Preserve every clip and reference, naming them in document order.
    clip_urls = {}
    for index, clip in enumerate(root.iter(ns + "clipPath")):
        name = f"mpl-clip-{index}"
        clip_urls[f"url(#{clip.get('id')})"] = f"url(#{name})"
        clip.set("id", name)
    for el in root.iter():
        if el.get("clip-path") in clip_urls:
            el.set("clip-path", clip_urls[el.get("clip-path")])
    root.set("role", "img")
    root.set("aria-labelledby", "figure-title figure-description")
    for tag, name, value in (("title", "figure-title", "Zebra escape response"),
                             ("desc", "figure-description", caption)):
        el = ET.Element(ns + tag, id=name)
        el.text = value
        root.insert(0, el)

    def jpeg(parent, path, **geometry):
        el = ET.SubElement(parent, ns + "image", geometry)
        el.set("preserveAspectRatio", "none")
        el.set(XLINK, "data:image/jpeg;base64," +
               base64.b64encode(path.read_bytes()).decode("ascii"))

    # Keep the map renderer's georeferencing and clip. Its backdrop is
    # white, including any no-data pixels in the RGBA mosaic.
    mosaic = root.find(f".//{ns}image")
    raster = Image.open(io.BytesIO(base64.b64decode(mosaic.get(XLINK).split(",")[1])))
    rgba = raster.convert("RGBA")
    rgb = Image.new("RGB", rgba.size, "white")
    rgb.paste(rgba, mask=rgba.getchannel("A"))
    buf = io.BytesIO()
    rgb.save(buf, "JPEG", quality=92, subsampling=0)
    mosaic.set(XLINK, "data:image/jpeg;base64," +
               base64.b64encode(buf.getvalue()).decode("ascii"))
    mosaic.attrib.pop("style", None)

    x0, y0, x1, y1 = PHOTO_RECT
    jpeg(root.find(".//*[@id='page-background']"), ctx["photo"],
         x=str(x0), y=str(y0), width=str(x1 - x0), height=str(y1 - y0))
    defs = root.find(ns + "defs")
    for key, spec in (("pre", INSET_PRE), ("post", INSET_POST)):
        cx, cy = spec["center"]
        side = spec["side"]
        clip_id = f"crop-{key}"
        clip = ET.SubElement(defs, ns + "clipPath", id=clip_id)
        ET.SubElement(clip, ns + "circle", cx=str(cx), cy=str(cy), r=str(spec["r"]))
        group = ET.Element(ns + "g", {"clip-path": f"url(#{clip_id})"})
        root.find(f".//*[@id='inset-{key}']").insert(0, group)
        jpeg(group, ctx["inset_assets"][key], x=str(-side / 2), y=str(-side / 2),
             width=str(side), height=str(side),
             transform=f"translate({cx} {cy}) rotate({spec['rotation']})")
    return ET.tostring(root, encoding="unicode") + "\n"


def draw_figure(ctx, *, family="Arial", style="published"):
    """Render the composite to an SVG string (fonts outlined, no text)."""
    import zebra_plots as plots

    if style not in ("published", "companion"):
        raise ValueError(f"Unknown zebra figure style: {style!r}")
    companion = style == "companion"
    sizes = COMPANION_FS if companion else FS
    polar_grid = ({"grid_color": "0.81", "grid_lw": 0.85} if companion
                  else {"grid_color": POLAR_GRID_COLOR, "grid_lw": POLAR_GRID_LW})
    bar_axes = ({"tick_color": "0.22", "tick_lw": 1.0} if companion
                else {"tick_color": GRID_COLOR, "tick_lw": TICK_LW})
    interval_lw = 1.6 if companion else LW_HDI
    numbers = ctx["numbers"]
    prep = numbers["stats"]
    nT = int(prep["nT"])
    theta = prep["theta"]
    pre, post = plots.PRE_COLOR, plots.POST_COLOR

    rc = {"font.family": family,
          "font.size": sizes["tick"],
          "legend.fontsize": sizes["map_legend"],
          "svg.fonttype": "path",
          "svg.hashsalt": "unwrap-zebra-publication-v1"}
    with mpl.rc_context(rc):
        fig = plt.figure(figsize=(PAGE_W / 72.0, PAGE_H / 72.0), dpi=72)
        fig.patch.set_facecolor("white")

        def page_axes(zorder):
            ax = fig.add_axes([0.0, 0.0, 1.0, 1.0], zorder=zorder)
            ax.set_xlim(0.0, PAGE_W)
            ax.set_ylim(PAGE_H, 0.0)
            ax.axis("off")
            ax.set_facecolor("none")
            return ax

        # -- layer 0: page backing ----------------------------------------
        base = page_axes(0)
        base.set_gid("page-background")
        base.add_patch(Rectangle((0, 0), PAGE_W, PAGE_H,
                                 facecolor="white", edgecolor="none",
                                 zorder=1))
        for rect in (MAP_FRAME_RECT, GRID_BACKING_RECT):
            x0, y0, x1, y1 = rect
            base.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0,
                                     facecolor="white", edgecolor="none",
                                     zorder=2))

        # -- panel a: the escape map ----------------------------------------
        map_ax = fig.add_axes(_figure_rect(MAP_RECT), aspect="equal",
                              zorder=4)
        plots.draw_escape_map(
            map_ax, ctx["img"], ctx["extent"], ctx["df_pos"],
            ctx["df_bearings"], FOCAL_TRACK,
            arrow_mode="distance", arrows_every_m=5.0, arrow_length=12.0,
            tail_width_pt=1.2 * MAP_ARROW_SCALE,
            head_width_pt=6.0 * MAP_ARROW_SCALE,
            head_length_pt=7.5 * MAP_ARROW_SCALE,
            arrow_shrink_pt=2.0 * MAP_ARROW_SCALE,
            arrow_edge_width=0.7 * MAP_ARROW_SCALE, arrow_alpha=0.87,
            other_stride=3, focal_stride=2, frame_max=FRAME_MAX,
            decimate_m=MAP_VERTEX_SPACING_M)
        map_ax.images[0].set_interpolation("none")
        plots._finish_escape_map(
            map_ax, ctx["extent"], PANEL_XLIM, PANEL_YLIM, None,
            "black", FOCAL_LW, "white", pre, post,
            tight_layout=False, legend_loc="center",
            legend_anchor=MAP_LEGEND_ANCHOR, fancybox=True,
            legend_arrow_size=9 if companion else 8 * MAP_ARROW_SCALE,
            legend_arrow_edge_width=MAP_ARROW_SCALE,
            framealpha=0.85, scale_bar=False, compass=False)

        # -- panel c: group-level mean direction (polar) -------------------
        ax_c = fig.add_axes(_figure_rect(
            (POLAR_C[0] - POLAR_R_EDGE_PT, POLAR_C[1] - POLAR_R_EDGE_PT,
             POLAR_C[0] + POLAR_R_EDGE_PT, POLAR_C[1] + POLAR_R_EDGE_PT)),
            projection="polar", zorder=5)
        _style_polar(ax_c, **polar_grid, compass=companion)
        group_lw = 3.0 if companion else LW_GROUP_LOOP
        ax_c.plot(theta, prep["r_group_pre"], color=pre, lw=group_lw)
        ax_c.plot(theta, prep["r_group_post"], color=post, lw=group_lw)
        ax_c.scatter([prep["mu_group_pre_mean"]], [R_DOT],
                     s=MS_GROUP_PRE ** 2, facecolor=pre, edgecolor="none",
                     marker=plots.PRE_MARK, zorder=5)
        ax_c.scatter([prep["mu_group_post_mean"]], [R_DOT],
                     s=MS_GROUP_POST ** 2, facecolor=post, edgecolor="none",
                     marker=plots.POST_MARK, zorder=5)

        # -- panel d: group-level alignment ---------------------------------
        ax_d = fig.add_axes(_figure_rect(BAR_D), zorder=5)
        _style_bar(ax_d, **bar_axes)
        ax_d.set_xlim(-0.5, 1.5)
        ax_d.set_ylim(0.0, 1.0)
        ax_d.set_xticks([0.0, 1.0])
        ax_d.set_xticklabels(["Pre", "Post"])
        ax_d.set_yticks([0.0, 0.5, 1.0])
        ax_d.set_ylabel("Mean resultant length", labelpad=10, fontsize=sizes["axis"])
        ax_d.vlines(0.0, prep["R_group_pre_lo"], prep["R_group_pre_hi"],
                    color=pre, lw=interval_lw, alpha=0.95, zorder=2)
        ax_d.vlines(1.0, prep["R_group_post_lo"], prep["R_group_post_hi"],
                    color=post, lw=interval_lw, alpha=0.95, zorder=2)
        ax_d.scatter([0.0], [prep["R_group_pre_mean"]], s=MS_GROUP_PRE ** 2,
                     facecolor=pre, edgecolor="none", marker=plots.PRE_MARK,
                     zorder=3)
        ax_d.scatter([1.0], [prep["R_group_post_mean"]], s=MS_GROUP_POST ** 2,
                     facecolor=post, edgecolor="none",
                     marker=plots.POST_MARK, zorder=3)

        # -- panel e: individual-level mean directions (polar) --------------
        ax_e = fig.add_axes(_figure_rect(
            (POLAR_E[0] - POLAR_R_EDGE_PT, POLAR_E[1] - POLAR_R_EDGE_PT,
             POLAR_E[0] + POLAR_R_EDGE_PT, POLAR_E[1] + POLAR_R_EDGE_PT)),
            projection="polar", zorder=5)
        _style_polar(ax_e, **polar_grid, compass=companion)
        for j in range(nT):
            ax_e.plot(theta, prep["r_ind_pre"][j], color=pre, lw=LW_IND_LOOP,
                      alpha=IND_LOOP_ALPHA)
            ax_e.plot(theta, prep["r_ind_post"][j], color=post,
                      lw=LW_IND_LOOP, alpha=IND_LOOP_ALPHA)
        for m_pre, m_post in zip(prep["emp_mu_pre"], prep["emp_mu_post"]):
            if not np.isnan(m_pre):
                ax_e.scatter([m_pre], [R_DOT], s=MS_EMP_PRE ** 2,
                             facecolor=pre, edgecolor="none",
                             marker=plots.PRE_MARK, zorder=5)
            if not np.isnan(m_post):
                ax_e.scatter([m_post], [R_DOT], s=MS_EMP_POST ** 2,
                             facecolor=post, edgecolor="none",
                             marker=plots.POST_MARK, zorder=5)

        # -- panel f: individual-level concentration -------------------------
        ax_f = fig.add_axes(_figure_rect(BAR_F), zorder=5)
        _style_bar(ax_f, **bar_axes)
        ax_f.set_xlim(-0.5, 1.5)
        ax_f.set_ylim(0.0, 1.0)
        ax_f.set_xticks([0.0, 1.0])
        ax_f.set_xticklabels(["Pre", "Post"])
        ax_f.set_yticks([0.0, 0.5, 1.0])
        ax_f.set_ylabel("Mean resultant length", labelpad=10, fontsize=sizes["axis"])
        rng = np.random.default_rng(JITTER_SEED)
        delta = rng.normal(0.0, JITTER_SD, size=nT)
        for j in range(nT):
            ax_f.plot([delta[j], 1.0 + delta[j]],
                      [prep["R_ind_pre_mean"][j],
                       prep["R_ind_post_mean"][j]],
                      color="0.60" if companion else "black",
                      lw=0.8 if companion else LW_PAIRED, zorder=0)
        for j in range(nT):
            lo, hi = prep["R_ind_pre_hdi"][j]
            ax_f.vlines(delta[j], lo, hi, color=pre, lw=interval_lw, zorder=1)
            ax_f.scatter([delta[j]], [prep["R_ind_pre_mean"][j]],
                         s=MS_CONC_PRE ** 2, facecolor=pre, edgecolor="none",
                         marker=plots.PRE_MARK, zorder=2)
        for j in range(nT):
            lo, hi = prep["R_ind_post_hdi"][j]
            ax_f.vlines(1.0 + delta[j], lo, hi, color=post, lw=interval_lw,
                        zorder=1)
            ax_f.scatter([1.0 + delta[j]], [prep["R_ind_post_mean"][j]],
                         s=MS_CONC_POST ** 2, facecolor=post, edgecolor="none",
                         marker=plots.POST_MARK, zorder=2)

        # -- overlay: insets, scale bar, north arrow, legend, text ----------
        overlay = page_axes(10)
        # resolve the aspect-adjusted map box before mapping page points
        fig.canvas.draw()
        for key in ("pre", "post"):
            spec = INSET_PRE if key == "pre" else INSET_POST
            _draw_inset(map_ax, overlay, fig, key, spec,
                        ctx["anchors"][key],
                        pre if key == "pre" else post,
                        cone_lw=1.75 if companion else CONE_LW,
                        ring_lw=3.0 if companion else RING_LW)

        for (x, y, w, h), col in zip(SCALE_RECTS, SCALE_RECT_COLORS):
            overlay.add_patch(Rectangle((x, y), w, h, facecolor=col,
                                        edgecolor="black", linewidth=1.2,
                                        zorder=8))
        overlay.text(SCALE_LABEL[0], SCALE_LABEL[1], "20 m",
                     ha="center", va="center", fontsize=sizes["scale"],
                     color="black", zorder=8)
        overlay.text(NORTH_X, NORTH_Y, "N", ha="center", va="center",
                     fontsize=sizes["north"], color="black", zorder=8)
        overlay.add_patch(FancyArrowPatch(
            (NORTH_X, NORTH_BASE_Y), (NORTH_X, NORTH_TIP_Y),
            arrowstyle=ArrowStyle.Simple(head_length=NORTH_HEAD_LEN,
                                         head_width=NORTH_HEAD_WIDTH,
                                         tail_width=NORTH_TAIL_WIDTH),
            mutation_scale=1.0, shrinkA=0, shrinkB=0,
            facecolor="black", edgecolor="black", linewidth=MAP_ARROW_SCALE,
            joinstyle="miter", zorder=8))

        if companion:
            handles = [
                Line2D([], [], color=col, lw=2.7, marker=marker,
                       markersize=7, markeredgewidth=0, label=label)
                for col, marker, label in (
                    (pre, plots.PRE_MARK, "Pre"), (post, plots.POST_MARK, "Post"))
            ]
            overlay.legend(handles=handles, loc="center",
                           bbox_to_anchor=(0.71, 0.022), ncols=2,
                           frameon=False, fontsize=sizes["bottom_legend"])
        else:
            bx0, by0, bx1, by1 = LEGEND_BOX
            overlay.add_patch(FancyBboxPatch(
                (bx0, by0), bx1 - bx0, by1 - by0,
                boxstyle="round,pad=0,rounding_size=6.0", mutation_aspect=1.0,
                facecolor="white", edgecolor=LEGEND_EDGE, linewidth=2.5,
                zorder=8))
            for (x0, ya, x1, yb), col in zip((LEGEND_PRE_LINE,
                                              LEGEND_POST_LINE), (pre, post)):
                overlay.plot([x0, x1], [ya, yb], color=col, lw=4.0,
                             solid_capstyle="butt", zorder=9)
            x, y, w, h = LEGEND_PRE_MARK
            overlay.add_patch(Rectangle((x, y), w, h, facecolor=pre,
                                        edgecolor="none", zorder=9))
            x, y, d, _ = LEGEND_POST_MARK
            overlay.add_patch(Circle((x + d / 2, y + d / 2), d / 2,
                                     facecolor=post, edgecolor="none", zorder=9))
            for text, tx, ty in LEGEND_TEXTS:
                overlay.text(tx, ty, text, ha="center", va="center",
                             fontsize=sizes["bottom_legend"], color="black",
                             zorder=9)

        for index, (text, tx, ty) in enumerate(TITLES):
            if companion:
                ty = (371.0, 727.0)[index // 2]
            overlay.text(tx, ty, text, ha="center",
                         va="top" if companion else "center",
                         fontsize=sizes["title"], color="black",
                         weight="bold" if companion else "normal",
                         linespacing=1.05, zorder=9)
        for character, lx, ly, color in LETTERS:
            if companion and character in "cdef":
                ly = 380.0 if character in "cd" else 736.0
                if character in "df":
                    lx = 649.1
            overlay.text(lx, ly, character + ".", ha="left", va="center",
                         fontsize=sizes["letter"], weight="bold", color=color,
                         zorder=9)
        angle_labels = ([(text, cx + dx, cy + dy) for cx, cy in (POLAR_C, POLAR_E)
                         for text, dx, dy in COMPASS_LABELS]
                        if companion else DEGREE_LABELS_C + DEGREE_LABELS_E)
        for text, tx, ty in angle_labels:
            overlay.text(tx, ty, text, ha="center", va="center",
                         fontsize=sizes["degree"], color="black", zorder=9)
        for text, tx, ty in RADIAL_LABEL_C + RADIAL_LABEL_E:
            overlay.text(tx, ty, text, ha="center", va="center",
                         fontsize=sizes["radial"], color="black", zorder=9)

        buffer = io.StringIO()
        fig.savefig(buffer, format="svg",
                    metadata={"Date": None,
                              "Creator":
                                  "unwrap zebra_publication_figure.py"})
        plt.close(fig)
    caption = COMPANION_CAPTION if companion else CAPTION
    return _embed_artwork(buffer.getvalue(), ctx, caption=caption)


# ---------------------------------------------------------------------------
# PDF export helpers
# ---------------------------------------------------------------------------
_PDF_ID_RE = re.compile(rb"/ID\[<([0-9A-Fa-f]{32})>"
                        rb"<([0-9A-Fa-f]{32})>\]")


def _stabilize_pdf_id(data):
    """Pin the trailer /ID array to a content-derived hash.

    MuPDF derives the /ID from write-time details, the last
    non-deterministic bytes of the export. Both slots are set to the
    16-byte (slot-width) truncation of the SHA-256 of the document
    with the ID region zeroed, so the file is a pure function of its
    content.
    """
    matches = _PDF_ID_RE.findall(data)
    if len(matches) != 1:
        raise ValueError(f"expected one trailer /ID, found {len(matches)}")
    m = _PDF_ID_RE.search(data)
    region = m.end(2) - m.start(1)
    base = data[:m.start(1)] + b"0" * region + data[m.end(2):]
    digest = hashlib.sha256(base).hexdigest()[:32].upper().encode("ascii")
    out = data[:m.start(1)] + digest + b"><" + digest + data[m.end(2):]
    if len(out) != len(data):
        raise ValueError("ID stabilisation changed the byte length")
    return out


def restore_pdf_source_images(pdf_path, svg):
    """Restore each SVG JPEG only to a PDF image with identical decoded pixels.

    CairoSVG stores decoded RGB pixels in the PDF. Matching those pixels
    before copying the JPEG stream preserves every placement and clip,
    including rotated crops and the orientation of the mosaic.
    """
    import pymupdf as fitz
    from PIL import Image

    with fitz.open(pdf_path) as doc:
        page = doc[0]
        rasters = {}
        for row in page.get_images(full=True):
            pix = fitz.Pixmap(doc, row[0])
            key = (pix.width, pix.height, hashlib.sha256(pix.samples).digest())
            rasters.setdefault(key, set()).add(row[0])

        scratch = doc.new_page()
        used = set()
        for el in ET.fromstring(svg).iter("{" + SVG_NS + "}image"):
            mime, encoded = el.get(XLINK).split(",", 1)
            if mime != "data:image/jpeg;base64":
                raise ValueError("Expected an embedded JPEG in the zebra figure")
            data = base64.b64decode(encoded)
            with Image.open(io.BytesIO(data)) as im:
                key = (*im.size, hashlib.sha256(im.convert("RGB").tobytes()).digest())
            targets = rasters.get(key)
            if not targets:
                raise ValueError("No PDF raster matches the decoded SVG JPEG")
            src = scratch.insert_image(fitz.Rect(0, 0, 1, 1), stream=data)
            for dst in targets:
                doc.xref_copy(src, dst)
                doc.xref_set_key(dst, "Interpolate", "true")
            used.update(targets)
        if used != {xref for targets in rasters.values() for xref in targets}:
            raise ValueError("PDF contains an unexpected raster")
        doc.delete_page(-1)
        doc.set_metadata(dict(doc.metadata, creationDate="", modDate=""))
        data = doc.tobytes(garbage=4, deflate=True)
    pdf_path.write_bytes(_stabilize_pdf_id(data))


def _html_page(svg, *, caption):
    return f'''<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Zebra escape response</title>
<style>body{{margin:0;background:#f5f5f4;color:#262626;font:15px/1.65 Arial,sans-serif}}
main{{max-width:1180px;margin:40px auto;padding:0 24px}}figure{{margin:0}}.figure{{background:white;padding:20px;border:1px solid #e6e6e6}}svg{{display:block;width:100%;height:auto}}figcaption{{max-width:96ch;margin:24px auto}}a{{color:#665017}}</style>
<main><nav><a href="zebra_publication.svg">SVG</a> · <a href="zebra_publication.pdf">PDF</a></nav>
<figure><div class="figure">{svg}</div><figcaption>{caption}</figcaption></figure></main></html>'''


def write_outputs(output_dir, ctx, svg, metadata, *, family):
    """Write the SVG/PDF/PNG/HTML exports and the provenance JSON.

    Returns the output stem (the exports live at ``zebra_publication.{svg,
    pdf, png, html, json}`` next to it).
    """
    import cairosvg
    import pymupdf as fitz

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = output_dir / "zebra_publication"
    stem.with_suffix(".svg").write_text(svg)
    cairosvg.svg2pdf(bytestring=svg.encode(),
                     write_to=str(stem.with_suffix(".pdf")))
    restore_pdf_source_images(stem.with_suffix(".pdf"), svg)
    with fitz.open(stem.with_suffix(".pdf")) as document:
        document[0].get_pixmap(matrix=fitz.Matrix(PREVIEW_SCALE, PREVIEW_SCALE)).save(
            str(stem.with_suffix(".png")))
    caption = metadata.get("caption", CAPTION)
    stem.with_suffix(".html").write_text(_html_page(svg, caption=caption))
    metadata["outputs_sha256"] = {s: sha256(stem.with_suffix(s))
                                  for s in (".svg", ".pdf", ".png",
                                            ".html")}
    stem.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n")
    return stem


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Rebuild the zebra escape publication figure from the "
                    "saved fit and pinned data (no sampling).")
    parser.add_argument("--data-root", type=Path,
                        default=ROOT / ".cache/zebra-figure" / DATA_REVISION)
    parser.add_argument("--download", action="store_true",
                        help="fetch any missing pinned unwrap-data files")
    parser.add_argument("--output-dir", type=Path,
                        help="export directory (default: assets/zebra-publication "
                             "or assets/zebra-companion, depending on style)")
    parser.add_argument("--style", choices=("published", "companion"),
                        default="published",
                        help="published reference or companion figure styling")
    parser.add_argument("--font-dir", type=Path,
                        help="directory of licensed font files")
    parser.add_argument("--font-family", default="Arial")
    args = parser.parse_args(argv)

    fonts = resolve_fonts(args.font_dir, args.font_family)
    root = prepare_inputs(args.data_root, download=args.download)
    artwork = check_assets()
    ctx = build_context(root)

    # the notebook's convergence gate on the cached fit (same read as
    # zebra_escape.ipynb): no estimate may be flagged on bulk ESS or R-hat
    from unwrap.circular import circular_summary
    summary_df = circular_summary(
        ctx["idata"],
        var_names=["Intercept", "PostStimulus", "1|Track_kappa", "1|Track",
                   "PostStimulus|Track_kappa", "PostStimulus|Track",
                   "kappa_Intercept", "kappa_PostStimulus",
                   "kappa_1|Track_sigma", "kappa_1|Track",
                   "kappa_PostStimulus|Track_sigma",
                   "kappa_PostStimulus|Track"],
        circ_var_names=["Intercept", "PostStimulus", "1|Track"],
    )
    flagged = len(summary_df[(summary_df["ess_bulk"] < 400)
                             | (summary_df["r_hat"] > 1.01)])
    if flagged:
        raise ValueError(
            f"{flagged} cached-fit estimates fail the notebook convergence "
            f"gate (bulk ESS < 400 or R-hat > 1.01)")

    svg = draw_figure(ctx, family=args.font_family, style=args.style)
    post = ctx["idata"].posterior
    numbers = ctx["numbers"]
    metadata = {
        "figure": "zebra_circ_stats_figure",
        "style": args.style,
        "layout": ("Published 900 x 1080-point composite geometry, measured "
                   "from the reference PDF (see ZEBRA_FIGURE.md)"),
        "data_repository": "jgraving/unwrap-data",
        "data_revision": DATA_REVISION,
        "inputs_sha256": INPUT_HASHES,
        "code_revision": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "python": platform.python_version(),
        "font_family": args.font_family,
        "fonts": fonts,
        "packages": {name: importlib.metadata.version(name) for name in
                     ("numpy", "scipy", "pandas", "matplotlib", "xarray",
                      "h5netcdf", "h5py", "arviz", "cairosvg", "pymupdf")},
        "posterior_chains": int(post.sizes["chain"]),
        "posterior_draws_per_chain": int(post.sizes["draw"]),
        "fit_meta": ctx["fit_meta"],
        "divergences": ctx["divergences"],
        "convergence": fit_diagnostics(ctx["idata"]),
        "notebook_gate_flagged_rows": 0,
        "map_window": {"observation": OBS, "focal_track": FOCAL_TRACK,
                       "frame_max": FRAME_MAX,
                       "xlim": list(PANEL_XLIM), "ylim": list(PANEL_YLIM),
                       "mosaic_qlim": list(MOSAIC_YLIM),
                       "mosaic_max_side_px": MOSAIC_MAX_SIDE_PX,
                       "mosaic_embed_side_px": MOSAIC_EMBED_SIDE,
                       "trajectory_vertex_spacing_m": MAP_VERTEX_SPACING_M},
        "artwork": artwork,
        "inset_layout": {"pre": INSET_PRE, "post": INSET_POST},
        "inset_anchors": ctx["anchors"],
        "inset_anchor_note": (
            "Anchor frames back-solved from the published callout cone "
            "apexes: frame 17623 sits ~7.8 m from the apex fitted to the "
            "paper's pre-cone, frame 80780 ~5 m from the post-cone apex "
            "at the path bend. The paper's hand-drawn cone lines deviate "
            "2.5-4 degrees from the exact tangents from the anchor "
            "positions. Measured compass bearings are recorded in inset_anchors. "
            "The two inset arrows follow the animals in the rotated photographs "
            "as drawn in the reference."),
        "scare_frame_note": SCARE_FRAME_NOTE,
        "summary_numbers": {k: v for k, v in numbers.items()
                            if k != "stats"},
        "caption": COMPANION_CAPTION if args.style == "companion" else CAPTION,
    }
    if args.style == "companion":
        group = numbers["group_mean_direction_deg"]
        individual = numbers["empirical_mean_direction_deg"]
        metadata["plotted_compass_bearings_deg"] = {
            "group_mean": {phase: (90.0 - angle) % 360.0
                           for phase, angle in group.items()},
            "individual_observed_mean": {
                phase: [None if angle is None else (90.0 - angle) % 360.0
                        for angle in angles]
                for phase, angles in individual.items()},
        }
    output_dir = args.output_dir
    if output_dir is None:
        output_dir = ROOT / "assets" / (
            "zebra-companion" if args.style == "companion" else "zebra-publication")
    stem = write_outputs(output_dir, ctx, svg, metadata, family=args.font_family)
    print(f"Wrote SVG, PDF, PNG, HTML and provenance JSON: {stem}")


if __name__ == "__main__":
    main()
