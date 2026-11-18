"""Reproduce the Goldschtein et al. bat publication figure (583 x 626 points).

This rebuilds the circstats paper's composite figure (published in
``notebooks/assets/goldshtein_bats_publication/reference.pdf``, polished
in Inkscape) the way the zebra figure does: from the notebooks' own
plotting code, the pinned model fits, and the pinned data, so notebook
changes flow into the paper figure.

* panels a, c, d, e, f, g and h are each rendered by the
  ``goldshtein_bats_plots.py`` function the ``goldshtein_bats.ipynb``
  cell calls, with that cell's arguments (panel a:
  ``plot_destination_panels_compact``; c: ``plot_total_elpd``;
  d: ``plot_delta_total_and_by_treatment``; e:
  ``plot_mean_hdi_horizontal``; f: ``plot_mixture_polar_pdfs``;
  g: ``plot_population_curve_smooth_compact``; h:
  ``plot_tracks_pgoal_small_bigfont``);
* panel b (the tagged-bat photograph, the treatment-sense icon grid,
  the "Treatment groups" label, the "Placeholder" watermark and the
  panel letters a-h) has no notebook counterpart and is pinned
  verbatim from the reference (`artwork.json` + the photo JPEG).

Data and fits are pinned to a ``jgraving/unwrap-data`` revision with
SHA-256 checksums (:data:`DATA_REVISION`, :data:`INPUT_HASHES`): the two
tracking workbooks, the lake outline, and the two cached fits
(``bats_vm_v1`` baseline von Mises, ``bats_vmu_v1`` hierarchical
mixture, both 4 chains x 2,000 draws, seed 123). The slim fit
snapshots strip the per-observation log-likelihood (vMU also strips the
per-observation deterministics), exactly as the notebook cells
recompute them on load: the baseline's log-likelihood from the rebuilt
bambi graph (cell 41) and the mixture's deterministics and
log-likelihood from :func:`goldshtein_bats_posthoc.restore_per_obs_
deterministics` plus the rebuilt raw-PyMC graph (cell 53).

Layout: the Inkscape composite embeds each notebook panel figure as a
scaled vector page. The per-panel scales (measured from the
reference's line widths, see ``GOLDSHTEIN_BATS_FIGURE.md``: a 1.001x,
c/d/e 1.151x, f 1.352x, g 1.385x, h 0.837x) and the measured anchors
(the panel axes boxes, in page points) are recorded below. This
generator re-renders each panel at its current natural size and embeds
it at the reference scale anchored on the panel's own rendered axes
box, so each panel keeps the published composite's position and scale
while its content (numbers, labels, styles) follows the current
notebook code. The two flight-map panels take their CartoDB basemaps
from the published figure's own tile rasters (``image_45.jpeg`` and
``image_23.jpeg``) through the ``basemap_img`` parameter of the map
functions, so no network access is needed.

Expected, documented deviations from the reference (it was assembled
from an older notebook state; see ``GOLDSHTEIN_BATS_FIGURE.md``): the
fitted numbers move with the pinned fits (panels c/d ELPD and deltas,
e's density HDIs, f's loops, g's curves); panel e is labelled
"Predictive density (mixture)" and panel h's colorbar "Posterior
directness" with evi_post_mean-coloured tracks (the current cell
calls); panel h's scale bar and north arrow sit in the rightmost
subplot (the current code's placement); panels a/g/h basemaps and
tracks are vector where the Inkscape reference rasterises some of
them.

Exports (``assets/goldshtein-bats-publication/`` by default): SVG
(fonts as text), PDF (the original basemap/photo JPEG sources restored
to their image objects at their measured placements, dates cleared and
/ID pinned to the content hash), a 3x PNG preview rasterised from the
PDF, self-contained HTML with the caption, and a provenance JSON
recording fonts, package versions, input hashes, fit metadata and the
numerical summaries.

Run from any directory::

    python notebooks/goldshtein_bats_publication_figure.py \
        [--data-root DIR] [--download] [--output-dir DIR] \
        [--font-dir DIR] [--font-family FAMILY]

See ``GOLDSHTEIN_BATS_FIGURE.md`` for provenance and the measured
layout notes.
"""
from __future__ import annotations

import argparse
import base64
import contextlib
import hashlib
import importlib.metadata
import io
import json
import math
import os
import platform
import re
import subprocess
import sys
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.font_manager import FontProperties, findfont  # noqa: E402
from matplotlib.path import Path as MplPath  # noqa: E402
from matplotlib.patches import PathPatch  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
ASSET_DIR = Path(__file__).resolve().parent / "assets" / \
    "goldshtein_bats_publication"
REFERENCE_PDF = ASSET_DIR / "reference.pdf"
ARTWORK_PATH = ASSET_DIR / "artwork.json"

PAGE_W = 583.0          # points
PAGE_H = 626.0          # points
FIGURE_NAME = "goldshtein_bats_publication"
PREVIEW_SCALE = 3      # 1749 x 1878 px, the reference raster scale

#: jgraving/unwrap-data revision carrying the bats data and the two
#: cached fits (the same revision the zebra figure pins).
DATA_REVISION = "86c41b5f43eb61191405b9bc85baf730a492a06f"

#: SHA-256 pins for every data/fit input (see the asset README and
#: GOLDSHTEIN_BATS_FIGURE.md).
INPUT_HASHES = {
    'fit_cache/bats_vm_v1.meta.json': 'ceff73f3751593763d76f2dfcf27f12a231e3e758d8875583aa4df4f83b8ca5f',
    'fit_cache/bats_vm_v1.nc': '430996b9aa3bc744dda1a7234d6be32515c349f9eec1e75a655b928bdd9cd30a',
    'fit_cache/bats_vmu_v1.meta.json': '38ca16ba9f3e98c11ddda4c0d03cd193e29e393b462e9d2577086c3a8f81d48f',
    'fit_cache/bats_vmu_v1.nc': '4b09c61855bb42f005e415813ebecc225c2a1526f6427823670fa54661e5d1e4',
    'goldshtein_bats/agamon_hula_lake_outline.geojson': 'dfbb1ebe1bb862250275b9784db0250e48b29068ea3cbfb3810826644ec7bbc2',
    'goldshtein_bats/all_bats_info.xls': '9196cc953fdc7f32fa235b51264c1993c38c3be9166e91169b098971af4788ad',
    'goldshtein_bats/dataAll_trackDest.xls': '1d6194f9ca52f4ef9522f073f7b45ee4af16448a67e6bf6353527bdeb331c036',
}

# ---------------------------------------------------------------------------
# Committed assets (the reference artwork; see the asset README)
# ---------------------------------------------------------------------------
#: SHA-256 pins for the committed assets. The Inkscape reference's other
#: rasters (panels a/h track-strip overlays and their SMask channels) are
#: superseded by the code-drawn vector tracks and are not committed.
ASSET_SHA256 = {
    "artwork.json":
        "9da0eb2639fea22d6af64820e1072524de6b159705b42cfe8d12c05e3d851730",
    "image_23.jpeg":
        "f5f042faee6ad39359f3cf757db286072139733e5241f2b458050fdfb7f01d1c",
    "image_279.jpeg":
        "6ffabea4e91fc7e53c7a97a97cb651642c8f8c392b1e4add3ae2fc7fad269734",
    "image_45.jpeg":
        "95c074d54e570b973f1b0c969d6e54d0f9ae56636a2f568664723919bf49975b",
    "reference.pdf":
        "849d8f791e9c53ffbcf559f4fb5a839fb770d44709b63991245a29cf35ae264f",
}

#: artwork.json contents the generator replays (schema 1 of the
#: reference extraction): 459 vector drawings, 8 image objects (12
#: placements) and 64 text spans.
N_EXPECTED_DRAWINGS = 459
N_EXPECTED_IMAGES = 8
N_EXPECTED_PLACEMENTS = 12
N_EXPECTED_SPANS = 64

#: Panel b artwork (the pinned chrome of this figure): the reference
#: drawings whose bounding-box centre sits in the icon-grid band between
#: the photo (bottom 92.57 pt, y down) and panel g's top edge, under the
#: photo's x extent. Measured count: 69.
N_EXPECTED_PANEL_B_DRAWINGS = 69
PANEL_B_X_MIN = 204.74 - 2.0        # photo left, minus slack (y down)
PANEL_B_X_MAX = 305.06 + 2.0        # photo right, plus slack
PANEL_B_Y_MAX = 252.0               # panel g top edge, generous cap

#: Panel letters are pinned typography (they belong to the Inkscape
#: chrome, not to the notebook figures).
PANEL_LETTERS = tuple(f"{c}." for c in "abcdefgh")

CAPTION = (
    "Hierarchical von Mises-Uniform mixture modeling of sensory-deprivation "
    "homing trajectories in bats. Free-flying bats were translocated within "
    "their home range and tracked during homing flights to a known goal "
    "(either the roost \"Home\" or the Agamon Hula Lake region). Step "
    "headings were converted to a circular bearing error by subtracting the "
    "instantaneous distance-minimizing target bearing so that y=0 indicates "
    "perfect alignment toward the goal and y=+/-pi indicates movement "
    "directly away from the goal. We fit a hierarchical von Mises-Uniform "
    "mixture (vMU) in which each step is generated either by a directed von "
    "Mises component or an undirected uniform component, allowing "
    "treatments to affect the directed-component mean mu, "
    "directed-component concentration kappa, and the directed-step mixing "
    "weight omega. (a) Representative trajectories for roost and lake goals "
    "with translocation locations and signal dropouts indicated (scale bar, "
    "1000 m). (b) A photo of a tagged bat, with the four treatment-sense "
    "combinations indicated. (c) Predictive model comparison using "
    "PSIS-LOO expected log predictive density (ELPD), contrasting a "
    "unimodal von Mises likelihood with the vMU; higher ELPD indicates "
    "better expected out-of-sample predictive performance. "
    "(d) ELPD difference (vMU minus von Mises), shown overall and within "
    "each treatment group. (e) Posterior predictive density at the goal "
    "direction (y*=0) under the fitted vMU, reported by treatment; the "
    "dashed vertical line indicates the uniform baseline density 1/(2pi). "
    "(f) Treatment-level directed-component von Mises densities on the "
    "circle (posterior mean and interval band). "
    "(g) Directed-component probability summarized along flights as a "
    "function of cumulative path distance from the endpoint (computed post "
    "hoc by smoothing posterior step-wise responsibilities); shaded bands "
    "indicate posterior uncertainty. (h) Example trajectories (one per "
    "treatment) colored by posterior mean directed-component probability "
    "(step-wise mixture responsibility). Additional per-flight posterior "
    "predictive checks, sensitivity analyses, and supplementary trajectory "
    "summaries are provided in the Supplementary Figures."
)

# ---------------------------------------------------------------------------
# Measured layout of the published composite (page points, measured from
# the reference PDF; derivation in GOLDSHTEIN_BATS_FIGURE.md).
#
# Boxes are (x0, y0, x1, y1) with y DOWN from the top-left page corner,
# matching the artwork.json extraction.
# ---------------------------------------------------------------------------
#: The two flight-map axes boxes of panel a (their basemap image
#: placements, image_45 x2).
PANEL_A_AXES = (
    (17.15, 16.48, 105.12, 238.91),
    (106.81, 16.48, 194.78, 238.91),
)
#: Panel c axis box (left + bottom spines of the reference).
PANEL_C_BOX = (404.10, 25.87, 565.07, 90.61)
#: Panel d axis box.
PANEL_D_BOX = (404.10, 135.22, 565.07, 215.50)
#: Panel e axis box.
PANEL_E_BOX = (34.26, 268.96, 195.23, 338.88)
#: The four polar axes boxes of panel f: the outer gray ring of each
#: quadrant sits at data radius EDGE_RADIUS + AMP = 2.0 of RMAX = 2.03,
#: so each polar box is the ring centre at half of
#: ring_diameter * RMAX / 2.0.
PANEL_F_AXES = (
    (20.735, 385.740, 118.795, 483.800),
    (125.420, 385.740, 223.480, 483.800),
    (20.735, 493.605, 118.795, 591.665),
    (125.420, 493.605, 223.480, 591.665),
)
#: Panel g axis box.
PANEL_G_BOX = (281.68, 259.09, 575.00, 340.46)
#: The four flight-map axes boxes of panel h (their basemap image
#: placements, image_23 x4).
PANEL_H_AXES = (
    (249.72, 378.30, 330.70, 583.06),
    (333.04, 378.30, 414.02, 583.06),
    (416.35, 378.30, 497.34, 583.06),
    (499.67, 378.30, 580.65, 583.06),
)
#: Panel h colorbar gradient bar (the reference's colorbar image_36
#: placement); the regenerated colorbar's bar is snapped to it.
PANEL_H_CBAR_BAR = (364.46, 586.24, 465.74, 591.26)
#: The tagged-bat photograph (image_279 placement; 0.4684 page alpha).
PANEL_B_PHOTO = (204.74, 6.62, 305.06, 92.57)
PANEL_B_PHOTO_ALPHA = 0.4684

#: Inkscape embed scale per panel: the reference line width of the
#: panel's own chrome divided by the notebook code's line width of the
#: same stroke (panel a border 1.0 pt; spines 0.8 pt; panel f outer
#: ring 0.8 pt; panel h border 1.5 pt). Exact artwork.json values:
#: a 1.0014435..., c/d/e 0.9209632..., f 1.0818636...,
#: g 1.1081346..., h 1.2558923....
SCALE_A = 1.0014435052871704      # 1.0014435052871704 / 1.0
SCALE_C = 0.9209632873535156 / 0.8   # = 1.1512041091918945
SCALE_D = 0.9209632873535156 / 0.8
SCALE_E = 0.9209632873535156 / 0.8
SCALE_F = 1.0818636417388916 / 0.8   # = 1.3523295521736145
SCALE_G = 1.1081346273422241 / 0.8   # = 1.3851682841777802
SCALE_H = 1.2558923959732056 / 1.5   # = 0.8372615973154704

# Notebook constants replicated by the generator (the notebook cells'
# own values).
SEED = 123
ROOST_LAT, ROOST_LON = 33.112373, 35.584135
METRIC_EPSG = 32636
STEP_LENGTH_M = 60.0
LAST_N_STEPS = 110
RECENT_WINDOW = 10

#: The vM/vMU snapshot metadata exactly as the notebook cells recorded
#: it (validated against the pinned ``.meta.json`` files on load).
VM_SNAPSHOT_META = {
    "model": "vm_baseline",
    "seed": int(SEED),
    "tune": 2000,
    "draws": 2000,
    "target_accept": 0.9,
    "stripped": ["log_likelihood"],
}
VMU_SNAPSHOT_META = {
    "model": "vmu_mixture_full",
    "seed": SEED,
    "tune": 2000,
    "draws": 2000,
    "target_accept": 0.9,
    "log_likelihood": True,
    "stripped": ["log_likelihood", "mu", "kappa", "p_goal"],
}


# ---------------------------------------------------------------------------
# Input verification
# ---------------------------------------------------------------------------
def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare_inputs(root, *, download: bool = False) -> Path:
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


def check_assets() -> dict:
    """Verify the committed inputs against :data:`ASSET_SHA256`."""
    manifest = {}
    for name, expected in ASSET_SHA256.items():
        path = ASSET_DIR / name
        if not path.is_file():
            raise FileNotFoundError(f"Missing asset {path}")
        actual = sha256(path)
        if actual != expected:
            raise ValueError(f"Asset checksum mismatch: {path}")
        manifest[name] = {"sha256": actual,
                          "bytes": path.stat().st_size}
    return manifest


def load_artwork() -> dict:
    """Load the curated reference artwork and check its counts."""
    asset = json.loads(ARTWORK_PATH.read_text())
    if asset.get("schema") != 1:
        raise ValueError(f"unexpected artwork schema in {ARTWORK_PATH}")
    if tuple(asset["source"]["page_pt"]) != (PAGE_W, PAGE_H):
        raise ValueError("asset page size does not match the reference")
    drawings = asset["drawings"]
    if len(drawings) != N_EXPECTED_DRAWINGS:
        raise ValueError(f"expected {N_EXPECTED_DRAWINGS} pinned drawings, "
                         f"found {len(drawings)}")
    images = asset["images"]
    placements = sum(len(im["placements"]) for im in images.values())
    if len(images) != N_EXPECTED_IMAGES or \
            placements != N_EXPECTED_PLACEMENTS:
        raise ValueError(f"expected {N_EXPECTED_IMAGES} images / "
                         f"{N_EXPECTED_PLACEMENTS} placements, found "
                         f"{len(images)} / {placements}")
    if len(asset["text_spans"]) != N_EXPECTED_SPANS:
        raise ValueError(f"expected {N_EXPECTED_SPANS} text spans, found "
                         f"{len(asset['text_spans'])}")
    return asset


def _drawing_bbox(drawing) -> tuple:
    xs = [float(v) for it in drawing["items"] for v in it[1::2]]
    ys = [float(v) for it in drawing["items"] for v in it[2::2]]
    return min(xs), min(ys), max(xs), max(ys)     # y-down page points


def panel_b_drawings(asset) -> dict:
    """The pinned panel b artwork (icon grid drawings) from the asset.

    A drawing belongs to the panel b artwork when its bounding-box
    centre sits in the icon-grid band (below the photo, above panel
    g, under the photo's x extent). The count is checked against the
    measured value so a stale asset fails loudly.
    """
    chosen = {}
    for name, d in asset["drawings"].items():
        x0, y0, x1, y1 = _drawing_bbox(d)
        cx = 0.5 * (x0 + x1)
        cy = 0.5 * (y0 + y1)
        if (PANEL_B_X_MIN <= cx <= PANEL_B_X_MAX
                and PANEL_B_PHOTO[3] - 2.0 <= cy <= PANEL_B_Y_MAX):
            chosen[int(name)] = d
    if len(chosen) != N_EXPECTED_PANEL_B_DRAWINGS:
        raise ValueError(
            f"expected {N_EXPECTED_PANEL_B_DRAWINGS} panel b artwork "
            f"drawings, found {len(chosen)}")
    return chosen


def chrome_spans(asset) -> list:
    """The pinned chrome text spans: the panel letters, 'Treatment
    groups' and the 'Placeholder' photo watermark (y-down origins)."""
    wanted = set(PANEL_LETTERS) | {"Treatment groups", "Placeholder"}
    spans = [s for s in asset["text_spans"] if s["text"] in wanted]
    if len(spans) != len(wanted):
        missing = wanted - {s["text"] for s in spans}
        raise ValueError(f"artwork is missing chrome spans: {missing}")
    return spans


def resolve_fonts(directory=None, family="Arial"):
    """Register ``directory`` fonts and record the resolved family faces.

    A missing face raises instead of silently substituting; the resolved
    font files and their hashes land in the provenance.
    """
    from matplotlib import font_manager

    if directory:
        for path in font_manager.findSystemFonts(
                fontpaths=[str(directory)]):
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


# ---------------------------------------------------------------------------
# Data and fits
# ---------------------------------------------------------------------------
def _import_notebook_modules(data_root):
    """Import the notebook's sibling modules (like the notebook cell 2).

    Satisfies the plots module's import-time call of
    ``bats_lake_outline_path()`` (which clones the default unwrap-data
    cache when missing) by leaving it an existing checkout marker when
    the real checkout is absent: an empty ``.git`` directory makes
    ``ensure_clone`` trust the location in place, and every real read
    below uses the pinned data root, which the lake path is re-pointed
    at here.
    """
    sys.path.insert(0, str(Path(__file__).parent))
    from unwrap import datasets

    default_root = datasets.DEFAULT_DATA_ROOT / "unwrap-data"
    if not default_root.exists():
        (default_root / ".git").mkdir(parents=True)
    import goldshtein_bats_preprocess as pre      # noqa: E402
    import goldshtein_bats_posthoc as posthoc     # noqa: E402
    import goldshtein_bats_plots as plots         # noqa: E402
    plots.LAKE_GEOJSON = str(datasets.bats_lake_outline_path(
        data_root=Path(data_root)))
    return pre, posthoc, plots


def _load_snapshot(root, name, expect):
    """Load and validate one pinned fit snapshot (issue-#130 format).

    Reuses the notebooks' own ``fit_cache.load_snapshot`` (sha256 of the
    netCDF file, the caller-supplied fit metadata, and the posterior
    variable layout) pointed at the pinned data root's ``fit_cache/``.
    """
    import fit_cache

    env = dict(os.environ)
    saved_enable = env.get(fit_cache.ENV_ENABLE)
    saved_dir = env.get(fit_cache.ENV_DIR)
    env[fit_cache.ENV_ENABLE] = "1"
    env[fit_cache.ENV_DIR] = str(Path(root) / "fit_cache")
    try:
        with _env_patch(env):
            return fit_cache.load_snapshot(name, expect=expect)
    finally:
        pass


@contextlib.contextmanager
def _env_patch(env):
    saved = {}
    for key, value in env.items():
        saved[key] = os.environ.get(key)
        os.environ[key] = value
    try:
        yield
    finally:
        for key, value in saved.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def build_context(data_root):
    """Load every input and assemble the drawing context.

    Returns a dict with the modelling table, the restored fits (both
    with their per-observation log-likelihoods, as the notebook cells
    restore them), every panel's posthoc tables, the two basemap
    rasters, the fit metadata / diagnostics, and the curated panel b
    artwork.
    """
    import pandas as pd
    import pymc as pm

    from unwrap.datasets import (bats_lake_outline_path, load_bats_info,
                                 load_bats_tracks)

    root = Path(data_root)
    pre, posthoc, plots = _import_notebook_modules(root)

    # ---- the notebook's preprocessing chain (cells 4-15, 39) ----------
    meta_df = load_bats_info(data_root=root)
    raw_track_df = load_bats_tracks(data_root=root)
    rediscretized = pre.rediscretize_latlon_parallel(
        raw_track_df,
        spacing_m=STEP_LENGTH_M,
        n_workers=16,
        include_last=False,
        gap_factor=10.0,
        gap_time_s=120.0,
        time_unit="days",
    )
    merged = pre.merge_tag_metadata(rediscretized, meta_df)
    track_w_treat = pre.add_dist_and_angle_to_destination(
        merged,
        LAKE_GEOJSON=bats_lake_outline_path(data_root=root),
        METRIC_EPSG=METRIC_EPSG,
        ROOST_LAT=ROOST_LAT, ROOST_LON=ROOST_LON,
        lake_grid_res_m=20.0,
        lake_domain_margin_m=None,
        lake_jacobi_tol=1e-5,
        lake_field="harmonic",
    ).copy()
    df_track = posthoc.build_modeling_table(track_w_treat)
    treatment_order = (
        df_track["treatment"].cat.categories.tolist()
        if hasattr(df_track["treatment"], "cat")
        else sorted(df_track["treatment"].dropna().unique())
    )

    # ---- the pinned fits, restored as the notebook cells do ----------
    # cell 53's data prep (shared by both restorations below)
    tag_codes, tag_idx = np.unique(df_track["tag"].to_numpy(),
                                   return_inverse=True)
    treat_idx = pd.Categorical(df_track["treatment"]).codes

    # baseline von Mises (cell 41): the slim snapshot stores the draws
    # without the per-observation log-likelihood, which the notebook
    # recomputes from the rebuilt bambi graph against the stored draws.
    idata_treatment = _load_snapshot(root, "bats_vm_v1", VM_SNAPSHOT_META)
    idata_treatment = _restore_vm_log_likelihood(
        idata_treatment, df_track, plotmod=plots, posthoc=posthoc)

    # vMU mixture (cell 53): restore the stripped per-observation
    # deterministics from the latent draws, then recompute the
    # per-observation log-likelihood with the rebuilt model graph.
    idata_mix = _load_snapshot(root, "bats_vmu_v1", VMU_SNAPSHOT_META)
    idata_mix = posthoc.restore_per_obs_deterministics(
        idata_mix, tag_idx, treat_idx)
    if getattr(idata_mix, "log_likelihood", None) is None:
        model_full = _vmum_model_graph(df_track, sample=False)
        idata_mix = pm.compute_log_likelihood(
            idata_mix, model=model_full, progressbar=False)

    for label, idata in (("baseline von Mises", idata_treatment),
                         ("vMU mixture", idata_mix)):
        if int(idata.sample_stats.diverging.sum()):
            raise ValueError(
                f"{label} fit has divergent transitions; review before "
                f"plotting")

    fit_meta = {}
    for name, expect in (("vm_baseline", VM_SNAPSHOT_META),
                         ("vmu_mixture_full", VMU_SNAPSHOT_META)):
        meta_path = root / "fit_cache" / f"{_SNAPSHOT_NC[name]}.meta.json"
        fit_meta[name] = json.loads(meta_path.read_text())

    # ---- panel tables (the notebook figure cells' statistics) ---------
    # c (cell 57/58): PSIS-LOO ELPD of both models, best first
    res_total, best_model, _loo_cache = \
        posthoc.compute_psisloo_elpd_summary(
            models={"von Mises-Uniform": idata_mix,
                    "von Mises": idata_treatment},
            df_fit=df_track,
            var_name=None,
            pareto_thresh=0.7,
            order_best_first=True,
        )

    # d (cell 60/61): delta ELPD, total and within treatment
    delta_df = posthoc.elpd_diff_total_and_by_treatment(
        idata_vmu=idata_mix, idata_vm=idata_treatment,
        df_fit=df_track,
        treat_col="treatment", treat_name_col="treatmentName",
        var_name=None,
    )

    # e (cell 70/73): population vMU summaries, mixture predictive
    # density at the goal direction by treatment
    res_full = posthoc.summarize_pop_p_kappa_mu_density0_mass45_alignment0_by_treatment(
        idata_mix_full=idata_mix, df_fit=df_track,
        treat_col="treatment", treat_name_col="treatmentName",
        hdi_prob=0.95, max_draws=2000,
        kappa_uniform=posthoc.KAPPA_UNIFORM,
        mass_half_width_deg=45.0,
    )
    dens_xmax = float(max(res_full["dens0_mix_hi"].max(),
                          res_full["dens0_dir_hi"].max()) * 1.05)
    dens_xlim = (0.0, dens_xmax)

    # h (cells 94/98): per-step p(goal) and posterior-directedness
    # columns on the modelling table
    df_track["_step"] = df_track.groupby("tag", observed=False).cumcount()
    df_track, _ = posthoc.append_pgoal_post_to_df(
        df_track, idata_mix, y_col="HeadingRadians",
        max_draws=2000, hdi_prob=0.95)
    df_track, _ = posthoc.append_evidence_to_df(
        df_track, idata_mix, y_col="HeadingRadians",
        max_draws=2000, ci_prob=0.90)

    # g (cells 100/104): population directedness vs cumulative path
    # distance from the end (recency-weighted smoothing, W = 10)
    df_track_with_evi, evi_draws = posthoc.append_evidence_to_df(
        df_track, idata_mix, y_col="HeadingRadians",
        max_draws=2000, ci_prob=0.90)
    _treat_levels, out_steps_k = posthoc.population_curve_kernel_per_tag(
        df_track_with_evi, evi_draws,
        x_kind="steps_from_end",
        last_n_steps=LAST_N_STEPS,
        steps_x_units="m",
        step_length_m=STEP_LENGTH_M,
        window=RECENT_WINDOW,
        gap_col="gap_before",
        hdi_prob=0.90,
        min_tags_at_x=2,
    )

    # ---- basemap rasters (the published figure's own tile mosaics) ----
    from PIL import Image
    basemap_a = np.asarray(
        Image.open(ASSET_DIR / "image_45.jpeg").convert("RGB"))
    basemap_h = np.asarray(
        Image.open(ASSET_DIR / "image_23.jpeg").convert("RGB"))

    # ---- curated panel b artwork --------------------------------------
    asset = load_artwork()

    return {
        "track_w_treat": track_w_treat,
        "meta_df": meta_df,
        "df_track": df_track,
        "idata_vm": idata_treatment,
        "idata_vmu": idata_mix,
        "treatment_order": list(treatment_order),
        "res_total": res_total,
        "best_model": best_model,
        "delta_df": delta_df,
        "res_full": res_full,
        "dens_xlim": dens_xlim,
        "out_steps_k": out_steps_k,
        "basemap_a": basemap_a,
        "basemap_h": basemap_h,
        "photo": ASSET_DIR / "image_279.jpeg",
        "panel_b_drawings": panel_b_drawings(asset),
        "chrome_spans": chrome_spans(asset),
        "fit_meta": fit_meta,
        "n_obs": int(len(df_track)),
        "n_tags": int(df_track["tag"].nunique()),
    }


_SNAPSHOT_NC = {"vm_baseline": "bats_vm_v1",
                "vmu_mixture_full": "bats_vmu_v1"}


def _restore_vm_log_likelihood(idata, df_track, *, plotmod, posthoc):
    """The notebook cell 41 cache-hit path: rebuild the bambi baseline
    model on the modelling table and evaluate its graph's
    per-observation log-likelihood against the stored draws."""
    import bambi as bmb
    import numpy as _np
    import pymc as pm

    from unwrap.link.bambi import LINKS
    from unwrap.link.numpy import inverse_softplus

    formula = bmb.Formula(
        "HeadingRadians ~ 1 + treatment + (1 | tag)",   # mu submodel
        "kappa ~ 1 + treatment + (1 | tag)",            # kappa submodel
    )
    priors = {
        "Intercept": bmb.Prior("VonMises", mu=0, kappa=1.0),
        "treatment": bmb.Prior("VonMises", mu=0, kappa=2.0),
        "1|tag": bmb.Prior(
            "Normal", mu=0.0,
            sigma=bmb.Prior("HalfNormal", sigma=1.0 / _np.sqrt(30.0)),
        ),
        "kappa": {
            "Intercept": bmb.Prior("Normal", mu=inverse_softplus(3),
                                   sigma=1.0),
            "treatment": bmb.Prior("Normal", mu=0, sigma=1.0),
            "1|tag": bmb.Prior(
                "Normal", mu=0.0,
                sigma=bmb.Prior("HalfNormal", sigma=1.0)),
        },
    }
    model = bmb.Model(
        formula,
        df_track,
        family="vonmises",
        priors=priors,
        link={"mu": LINKS["circ_mod"], "kappa": LINKS["softplus"]},
        center_predictors=False,
        categorical=["tag", "treatment"],
        dropna=True,
    )
    model.build()
    return pm.compute_log_likelihood(
        posthoc.graph_free_rvs(idata, model),
        model=model.backend.model,
        progressbar=False,
    )


def _vmum_model_graph(df_track, *, sample: bool = False):
    """The notebook cell 53 ``fit_vmum_model_full`` model graph.

    A verbatim copy of the cell's model (raw PyMC) so a cache hit can
    rebuild the graph without sampling, exactly the cell's own cache-hit
    path. ``sample=True`` would run the fit (unused by the generator).
    """
    import warnings

    import pandas as pd
    import pymc as pm
    import pytensor.tensor as pt

    from unwrap.link.pymc import circular_modulo
    from scipy.special import logit

    LOG_UNIF = -np.log(2 * np.pi)
    y = df_track["HeadingRadians"].to_numpy()
    _tag_codes, tag_idx = np.unique(df_track["tag"].to_numpy(),
                                    return_inverse=True)
    n_tag = len(_tag_codes)
    treat_cat = pd.Categorical(df_track["treatment"])
    treat_idx = treat_cat.codes
    n_treat = len(treat_cat.categories)

    def inverse_softplus(x):
        return np.log(np.exp(x) - 1.0)

    with pm.Model() as vmum_full:
        y_obs = pm.Data("y_obs", y)
        tag = pm.Data("tag", tag_idx)
        treat = pm.Data("treat", treat_idx)

        mu0 = pm.VonMises("mu0", mu=0.0, kappa=1.0)
        mu_tag_kappa = pm.LogNormal("mu_tag_kappa", mu=np.log(30.0),
                                    sigma=0.1)
        mu_tag = pm.VonMises("mu_tag", mu=0.0, kappa=mu_tag_kappa,
                             shape=n_tag)
        mu_treat = pm.VonMises("mu_treat", mu=0.0, kappa=2.0,
                               shape=n_treat - 1)
        mu_treat_full = pt.concatenate([pt.zeros((1,)), mu_treat])
        mu = pm.Deterministic(
            "mu", circular_modulo(mu0 + mu_treat_full[treat] + mu_tag[tag]))

        eta_kappa0 = pm.Normal("eta_kappa0",
                               mu=inverse_softplus(5.0), sigma=1.0)
        sigma_kappa_tag = pm.HalfNormal("sigma_kappa_tag", sigma=1.0)
        eta_kappa_tag = pm.Normal("eta_kappa_tag", mu=0.0, sigma=1.0,
                                  shape=n_tag)
        eta_kappa_base = eta_kappa0 + eta_kappa_tag[tag] * sigma_kappa_tag
        kappa_treat = pm.Normal("kappa_treat", mu=0.0, sigma=1.0,
                                shape=n_treat - 1)
        kappa_treat_full = pt.concatenate([pt.zeros((1,)), kappa_treat])
        eta_kappa = eta_kappa_base + kappa_treat_full[treat]
        kappa = pm.Deterministic("kappa", pt.softplus(eta_kappa))

        logit_p0 = pm.Normal("logit_p0", mu=logit(0.8), sigma=1.0)
        sigma_p_tag = pm.HalfNormal("sigma_p_tag", sigma=1.0)
        logit_p_tag = pm.Normal("logit_p_tag", mu=0.0, sigma=1.0,
                                shape=n_tag)
        logit_p_base = logit_p0 + logit_p_tag[tag] * sigma_p_tag
        beta_treat = pm.Normal("beta_treat", mu=0.0, sigma=1.0,
                               shape=n_treat - 1)
        beta_treat_full = pt.concatenate([pt.zeros((1,)), beta_treat])
        logit_p = logit_p_base + beta_treat_full[treat]
        w = pm.Deterministic(
            "p_goal", pm.math.sigmoid(pt.stack([logit_p, -logit_p], axis=1)))

        comp_goal = pm.VonMises.dist(mu=mu, kappa=kappa)
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", message="DensityDist has been deprecated.*",
                category=FutureWarning)
            comp_unif = pm.DensityDist.dist(
                logp=lambda v: pt.full_like(v, LOG_UNIF))
        pm.Mixture("Heading_obs", w=w, comp_dists=[comp_goal, comp_unif],
                   observed=y_obs)

        if sample:
            pm.sample(tune=2000, draws=2000, target_accept=0.9,
                      random_seed=SEED, return_inferencedata=True,
                      idata_kwargs={"log_likelihood": True})

    return vmum_full


def fit_diagnostics(idata) -> dict:
    """Independent ArviZ reads of a cached fit's convergence statistics
    (recorded for provenance, not a new sample)."""
    import arviz as az

    posterior = idata.posterior

    def _to_dataset(obj):
        if hasattr(obj, "ds"):
            return obj.ds
        return obj

    try:
        return {
            "max_rhat_split": float(_to_dataset(az.rhat(posterior))
                                    .to_array().max()),
            "max_rhat_rank": float(_to_dataset(az.rhat(posterior,
                                                       method="rank"))
                                   .to_array().max()),
            "min_ess_bulk": float(_to_dataset(az.ess(posterior,
                                                     method="bulk"))
                                  .to_array().min()),
            "min_ess_tail": float(_to_dataset(az.ess(posterior,
                                                     method="tail"))
                                  .to_array().min()),
        }
    except Exception as exc:       # diagnostics must not gate exports
        return {"error": f"{type(exc).__name__}: {exc}"}


# ---------------------------------------------------------------------------
# Panel drawing (the notebook cells' calls, verbatim arguments)
# ---------------------------------------------------------------------------
def _panel_rc(family):
    """The house figure style (unwrap.plots DEFAULT_RC) with the resolved
    family substituted for Arial; the panels render at display dpi."""
    from unwrap import plots as _plots

    rc = dict(_plots._default_rc())
    rc["font.family"] = family
    rc["font.sans-serif"] = [family]
    return rc


def draw_panel_a(ctx, family):
    """Panel a (notebook cell 23), offline basemap from the asset."""
    plots = _PLOTS
    with mpl.rc_context(_panel_rc(family)):
        fig, _axes = plots.plot_destination_panels_compact(
            ctx["track_w_treat"], ctx["meta_df"],
            fig_w=2.5, fig_h=3.8, title_fs=9, legend_fs=6,
            basemap_img=ctx["basemap_a"],
            save=False, save_path=None)
        return fig


def draw_panel_c(ctx, family):
    """Panel c (notebook cell 58)."""
    plots = _PLOTS
    with mpl.rc_context(_panel_rc(family)):
        fig, _ax = plots.plot_total_elpd(
            ctx["res_total"], title=None, xlabel="ELPD",
            vline_at_best=False, save=False, save_path=None)
        return fig


def draw_panel_d(ctx, family):
    """Panel d (notebook cell 61)."""
    plots = _PLOTS
    with mpl.rc_context(_panel_rc(family)):
        fig, _ax = plots.plot_delta_total_and_by_treatment(
            ctx["delta_df"],
            treatment_order=ctx["treatment_order"],
            title=None, xlabel="ELPD difference", vline_at_zero=True,
            save=False, save_path=None)
        return fig


def draw_panel_e(ctx, family):
    """Panel e (notebook cell 73)."""
    plots = _PLOTS
    with mpl.rc_context(_panel_rc(family)):
        fig, _ax = plots.plot_mean_hdi_horizontal(
            ctx["res_full"], "dens0_mix_mean", "dens0_mix_lo",
            "dens0_mix_hi", "Predictive density (mixture)",
            xlim=ctx["dens_xlim"], vline_x=1.0 / (2.0 * np.pi),
            show_treatment_labels=False, show_treatment_ticks=False,
            save=False, save_path=None)
        return fig


def draw_panel_f(ctx, family):
    """Panel f (notebook cell 68)."""
    plots = _PLOTS
    with mpl.rc_context(_panel_rc(family)):
        fig = plots.plot_mixture_polar_pdfs(
            ctx["idata_vmu"], ctx["df_track"],
            save=False, save_path=None)
        return fig


def draw_panel_g(ctx, family):
    """Panel g (notebook cell 104)."""
    plots = _PLOTS
    with mpl.rc_context(_panel_rc(family)):
        fig = plots.plot_population_curve_smooth_compact(
            ctx["out_steps_k"], ctx["treatment_order"],
            xlabel="Cumulative path distance from end (m)",
            ylabel="Posterior\ndirectedness",
            invert_x=False,
            xlim=(-LAST_N_STEPS * STEP_LENGTH_M, 0),
            x_tick_abs=True,
            save=False, save_path=None)
        return fig


def _draw_panel_h(ctx, family, cbar_box=None):
    """Panel h (notebook cell 98), offline basemap from the asset.

    ``cbar_box`` (figure-fraction axes box) snaps the colorbar's
    gradient bar onto the published composite's measured bar when one
    is given.
    """
    plots = _PLOTS
    with mpl.rc_context(_panel_rc(family)):
        fig = plots.plot_tracks_pgoal_small_bigfont(
            ctx["track_w_treat"],
            p_col="evi_post_mean",
            df_prob=ctx["df_track"],
            prob_tag_col="tag",
            prob_step_col="_step",
            prob_p_col="evi_post_mean",
            gap_col="gap_before",
            gap_linestyle="--",
            gap_alpha=0.75,
            seg_lw=0.8,
            wspace=0.03,
            hspace=0.02,
            cbar_label=plots.EVIDENCE_LABEL,
            basemap_img=ctx["basemap_h"],
            cbar_box=cbar_box,
            save=False, save_path=None)
        return fig


def draw_panel_h(ctx, family):
    """Panel h, with the colorbar snapped to the measured bar box.

    The colorbar axis box is a figure fraction, so it is solved in two
    deterministic passes: an initial render at the cell's own fractions
    fixes the panel's translate (anchored on the first axes box), and
    the second render, when needed, passes the fraction that lands the
    gradient bar on :data:`PANEL_H_CBAR_BAR` exactly.
    """
    fig = _draw_panel_h(ctx, family)
    fig.canvas.draw()
    ax0 = fig.axes[0]
    tx, ty, _w, _h = _translate(fig, ax0, PANEL_H_AXES[0], SCALE_H)
    cax0 = fig.axes[-1]
    box = _figure_box_pt(fig, cax0)          # figure pt, y up
    w_pt, h_pt = (v * 72.0 for v in fig.get_size_inches())
    need = _bar_box_error(box, tx, ty, SCALE_H, h_pt, PANEL_H_CBAR_BAR)
    if need > 0.25:
        # the measured bar box (page pt, y down) as a figure fraction
        # (y up from the figure bottom edge): the part's top edge sits
        # at ty (y down), so a figure point p (y up) maps to
        # ty + SCALE_H * (h_pt - p)
        x0d, y0d, x1d, y1d = PANEL_H_CBAR_BAR
        frac = [
            (x0d - tx) / (SCALE_H * w_pt),
            1.0 - (y1d - ty) / (SCALE_H * h_pt),
            (x1d - x0d) / (SCALE_H * w_pt),
            (y1d - y0d) / (SCALE_H * h_pt),
        ]
        plt.close(fig)
        fig = _draw_panel_h(ctx, family, cbar_box=frac)
    return fig


PANELS = {
    "a": (draw_panel_a, SCALE_A, PANEL_A_AXES[0]),
    "c": (draw_panel_c, SCALE_C, PANEL_C_BOX),
    "d": (draw_panel_d, SCALE_D, PANEL_D_BOX),
    "e": (draw_panel_e, SCALE_E, PANEL_E_BOX),
    "f": (draw_panel_f, SCALE_F, PANEL_F_AXES[0]),
    "g": (draw_panel_g, SCALE_G, PANEL_G_BOX),
    "h": (draw_panel_h, SCALE_H, PANEL_H_AXES[0]),
}
#: Panel letters that have no code-drawn panel (b) sit between a and c.
PANEL_ORDER = ("a", "b", "c", "d", "e", "f", "g", "h")

# the plots module is imported once (build_context); the panel drawers
# read it through this indirection so importing this module stays cheap
# and side-effect free.
_PLOTS = None


# ---------------------------------------------------------------------------
# Chrome (panel b artwork + letters), the verbatim reference renderer
# ---------------------------------------------------------------------------
def path_from_items(items) -> MplPath:
    """Convert a PyMuPDF item list into a matplotlib path.

    Items come in page points with y down from the top-left corner; the
    flip to matplotlib y-up happens here (``PAGE_H - y``).  Kinds: ``l``
    line, ``c`` cubic Bezier, ``re`` rectangle, ``qu`` quad (the four
    corners in ul/ur/ll/lr order, replayed as a closed polygon).
    Subpaths that close on themselves are terminated with CLOSEPOLY so
    strokes join cleanly.
    """
    def flip(x, y):
        return (x, PAGE_H - y)

    verts: list[tuple[float, float]] = []
    codes: list[int] = []
    cur = None          # end of the current subpath (flipped coords)
    starts: list = []   # open subpath starts

    def begin(x: float, y: float):
        nonlocal cur
        p = flip(x, y)
        if cur is None or abs(cur[0] - p[0]) > 1e-9 \
                or abs(cur[1] - p[1]) > 1e-9:
            verts.append(p)
            codes.append(MplPath.MOVETO)
            starts.append([p])
            cur = p

    for it in items:
        kind = it[0]
        if kind == "l":
            x1, y1, x2, y2 = (float(v) for v in it[1:5])
            begin(x1, y1)
            verts.append(flip(x2, y2))
            codes.append(MplPath.LINETO)
            cur = verts[-1]
        elif kind == "c":
            (x1, y1), (x2, y2), (x3, y3), (x4, y4) = (
                (float(it[k]), float(it[k + 1])) for k in (1, 3, 5, 7))
            begin(x1, y1)
            verts.append(flip(x2, y2))
            codes.append(MplPath.CURVE4)
            verts.append(flip(x3, y3))
            codes.append(MplPath.CURVE4)
            verts.append(flip(x4, y4))
            codes.append(MplPath.CURVE4)
            cur = verts[-1]
        elif kind == "re":
            x0, y0, x1, y1 = (float(v) for v in it[1:5])
            if cur is not None:
                codes.append(MplPath.CLOSEPOLY)
                verts.append(cur)
                starts.pop()
            cur = None
            a = flip(x0, y0)
            b = flip(x1, y0)
            c = flip(x1, y1)
            dd = flip(x0, y1)
            verts.append(a)
            codes.append(MplPath.MOVETO)
            starts.append([a])
            for p in (b, c, dd):
                verts.append(p)
                codes.append(MplPath.LINETO)
            codes.append(MplPath.CLOSEPOLY)
            verts.append(a)
            starts.pop()
        elif kind == "qu":
            ul = (float(it[1]), float(it[2]))
            ur = (float(it[3]), float(it[4]))
            ll = (float(it[5]), float(it[6]))
            lr = (float(it[7]), float(it[8]))
            if cur is not None:
                codes.append(MplPath.CLOSEPOLY)
                verts.append(cur)
                starts.pop()
            cur = None
            pts = [flip(*p) for p in (ul, ur, lr, ll)]
            verts.append(pts[0])
            codes.append(MplPath.MOVETO)
            for p in pts[1:]:
                verts.append(p)
                codes.append(MplPath.LINETO)
            codes.append(MplPath.CLOSEPOLY)
            verts.append(pts[0])
            starts.pop()
        else:
            raise ValueError(f"unknown drawing item kind {kind!r}")

    for start in starts:
        if cur is not None and abs(cur[0] - start[0][0]) <= 1e-9 \
                and abs(cur[1] - start[0][1]) <= 1e-9:
            codes.append(MplPath.CLOSEPOLY)
            verts.append(cur)
    return MplPath(verts, codes)


def _color(c, opacity=None):
    if c is None:
        return "none"
    rgb = tuple(float(v) for v in c[:3])
    if opacity is None or opacity in (1.0, 1):
        return rgb
    return rgb + (float(opacity),)


def draw_drawing(ax, dr, zorder):
    """Render one pinned drawing (fill and/or stroke)."""
    path = path_from_items(dr["items"])
    t = dr["type"]
    face = _color(dr["fill"], dr.get("fill_opacity")) if "f" in t \
        else "none"
    if "s" in t:
        edge = _color(dr["color"], dr.get("stroke_opacity"))
        lw = dr.get("width") or 0.0
        ls = "solid"
        dash = dr.get("dashes")
        if dash and lw:
            ls = (0, (dash["on"] / lw, dash["off"] / lw))
        patch = PathPatch(path, facecolor=face, edgecolor=edge,
                          linewidth=lw, linestyle=ls,
                          capstyle="butt", joinstyle="miter",
                          zorder=zorder)
    else:
        patch = PathPatch(path, facecolor=face, edgecolor="none",
                          zorder=zorder)
    ax.add_patch(patch)
    return patch


def draw_text(ax, spans, *, family_override=None, zorder=1000):
    """Typeset every span at its measured baseline origin.

    Each span carries the exact baseline origin measured from the
    embedded face (page points, y down) and the line direction, so the
    substitute face only needs the same family metrics: the anchor sits
    at the origin (``rotation_mode="anchor"`` keeps ``ha``/``va`` in the
    unrotated frame) and the line of text follows the measured
    direction angle.
    """
    drawn = 0
    skipped = 0
    for s in spans:
        text = s["text"]
        if not text.strip(" \xa0"):
            skipped += 1
            continue
        bold = bool(s["flags"] & 16)
        stack = list(FONT_STACK)
        if family_override:
            stack = [family_override] + [
                f for f in stack if f != family_override]
        fp = FontProperties(family=stack, size=s["size"],
                            weight="bold" if bold else "normal")
        color = s["color"]
        if color.startswith("#00") and len(color) == 9:
            color = "#" + color[3:]  # drop the leading zero byte
        ox, oy = s["origin"]          # baseline start, y down
        dx, dy = s["dir"]
        rotation = math.degrees(math.atan2(-dy, dx))
        ax.text(ox, PAGE_H - oy, text, fontproperties=fp, color=color,
                ha="left", va="baseline", rotation=rotation,
                rotation_mode="anchor", alpha=s.get("fill_alpha", 1.0),
                zorder=zorder)
        drawn += 1
    return drawn, skipped


#: The reference embeds one subset-TrueType face (the paper's Arial
#: stack, embedded by Inkscape under a stripped name).  Substitute
#: metric equivalents: the resolved family, then Liberation Sans in CI,
#: then DejaVu Sans elsewhere.
FONT_STACK = ["Liberation Sans", "DejaVu Sans"]


def draw_chrome(ctx, *, family=None):
    """Render the page-sized chrome layer: the panel b artwork (icon
    grid), the panel letters and the chrome labels at their pinned
    positions. The photograph is drawn into the final exports directly
    from the original stream (its 0.4684 page alpha has no stored
    SMask in the reference)."""
    fig = plt.figure(figsize=(PAGE_W / 72.0, PAGE_H / 72.0), dpi=72)
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_xlim(0, PAGE_W)
    ax.set_ylim(0.0, PAGE_H)
    ax.axis("off")
    ax.set_facecolor("none")
    # transparent figure: the page background is drawn by the
    # composition itself, and a white full-page patch would hide the
    # photograph (and every part embedded after it) in the PDF export
    fig.patch.set_facecolor("none")

    drawings = ctx["panel_b_drawings"]
    for z, i in enumerate(sorted(drawings)):
        draw_drawing(ax, drawings[i], zorder=10 + z)
    drawn, skipped = draw_text(ax, ctx["chrome_spans"],
                               family_override=family)
    stats = {"drawings": len(drawings),
             "text_drawn": drawn,
             "text_skipped": skipped}
    return fig, stats


# ---------------------------------------------------------------------------
# Placement of a natural panel figure on the page
# ---------------------------------------------------------------------------
#: matplotlib figures render at 100 dpi for the panels; the SVG unit is
#: the 100-dpi pixel, the PDF unit the point. 1 page pt = FIG_PX_TO_PT^-1
#: SVG units... i.e. SVG units convert to points by this factor.
FIG_PX_TO_PT = 72.0 / 100.0


def _figure_box_pt(fig, ax):
    w_pt, h_pt = (v * 72.0 for v in fig.get_size_inches())
    p = ax.get_position()
    return (p.x0 * w_pt, p.y0 * h_pt,
            (p.x0 + p.width) * w_pt, (p.y0 + p.height) * h_pt)


#: This pymupdf's page-geometry APIs (insert_image, show_pdf_page,
#: draw_*) take page coordinates with the origin at the top-left
#: corner (y down; the PDF content streams themselves stay y up), so
#: the measured layout boxes (all y down) and the translates below are
#: passed to them as-is.
def _translate(fig, ax, measured_ydown, scale):
    """The (tx, ty) page-point translate (y DOWN: ty is the part's top
    edge on the page) that places the panel's rendered ``ax`` box on the
    measured reference box at ``scale``."""
    w_pt, h_pt = (v * 72.0 for v in fig.get_size_inches())
    p = ax.get_position()
    nx0 = p.x0 * w_pt
    ny0_down = h_pt - (p.y0 + p.height) * h_pt   # anchor top, y down
    mx0, my0, _mw, _mh = measured_ydown
    return (mx0 - scale * nx0, my0 - scale * ny0_down,
            scale * w_pt, scale * h_pt)


def _bar_box_error(box_pt, tx, ty_down, scale, h_pt, measured_ydown):
    """How far (points) the rendered colorbar bar is from the measured
    bar box, after the panel's (tx, ty_down, scale) embed. ``box_pt`` is
    the bar box in figure points, y up from the figure's bottom edge."""
    bx0, by0, bx1, by1 = box_pt
    placed = (tx + scale * bx0,
              ty_down + scale * (h_pt - by1),
              tx + scale * bx1,
              ty_down + scale * (h_pt - by0))
    mx0, my0, mx1, my1 = measured_ydown
    return max(abs(placed[0] - mx0), abs(placed[1] - my0),
               abs(placed[2] - mx1), abs(placed[3] - my1))


# ---------------------------------------------------------------------------
# Export composition
# ---------------------------------------------------------------------------
SVG_NS = "http://www.w3.org/2000/svg"
XLINK_NS = "http://www.w3.org/1999/xlink"
_URL_REF = re.compile(r"url\(#([^)]+)\)")


def _svg_bytes(fig) -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format="svg", metadata={"Date": None})
    plt.close(fig)
    return buf.getvalue()


def _pdf_bytes(fig) -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format="pdf", metadata={"Date": None})
    plt.close(fig)
    return buf.getvalue()


def _prefix_svg_ids(root, prefix: str):
    """Prefix every ``id`` and ``url(#id)`` reference in an SVG part so
    several parts can share one document."""
    for el in root.iter():
        id_ = el.get("id")
        if id_ is not None:
            el.set("id", prefix + id_)
        for key, value in el.attrib.items():
            if value.startswith("url(#"):
                m = _URL_REF.search(value)
                if m:
                    el.set(key, value.replace(
                        f"#{m.group(1)}", f"#{prefix}{m.group(1)}"))


def _split_svg(svg_bytes):
    """Split a part SVG into (defs, content); content is every root
    child except ``defs``/``metadata``/``desc``/``title``."""
    root = ET.fromstring(svg_bytes)
    defs = None
    content = []
    for child in root:
        tag = child.tag.split("}", 1)[-1]
        if tag == "defs":
            defs = child
        elif tag in ("metadata", "desc", "title"):
            continue
        else:
            content.append(child)
    return (defs, content)


def _swap_panel_rasters(content, expected_count, jpeg_path):
    """Re-emit the part's basemap <image> elements from the original
    JPEG (the matplotlib SVG backend resamples every raster).

    The basemaps are the only *large* rasters in a part SVG (the panel
    h colorbar gradient is a narrow strip and is left as-is), so every
    data-image taller than the threshold is replaced in place (geometry
    untouched: the tracks are drawn on the same axes and must stay
    aligned); a wrong count fails loudly.
    """
    ns = "{%s}" % SVG_NS
    b64 = base64.b64encode(jpeg_path.read_bytes()).decode("ascii")
    swapped = 0
    for el in content:
        for img in el.iter(ns + "image"):
            href = img.get("{%s}href" % XLINK_NS) or img.get("href")
            if not href or not href.startswith("data:image"):
                continue
            if float(img.get("height")) < 50.0:
                continue     # narrow strip: the colorbar gradient
            geom = {k: img.get(k) for k in ("x", "y", "width", "height")}
            img.attrib.clear()
            for key, value in geom.items():
                img.set(key, value)
            img.set("preserveAspectRatio", "none")
            img.set("{%s}href" % XLINK_NS, "data:image/jpeg;base64," + b64)
            swapped += 1
    if swapped != expected_count:
        raise ValueError(
            f"expected {expected_count} basemap rasters in the part "
            f"SVG, found {swapped}")


def compose_svg(ctx, parts, family):
    """Compose the page SVG: white backing, the photograph, the chrome
    layer, and each panel part wrapped in its measured
    translate/scale group.

    ``parts`` maps panel letter -> (svg_bytes, placement dict); the
    panel b slot (the chrome, page units at scale 1) is in the same
    table.
    """
    ET.register_namespace("", SVG_NS)
    ET.register_namespace("xlink", XLINK_NS)
    ns = "{%s}" % SVG_NS

    root = ET.Element(ns + "svg", {
        "width": f"{PAGE_W:g}pt", "height": f"{PAGE_H:g}pt",
        "viewBox": f"0 0 {PAGE_W:g} {PAGE_H:g}"})
    root.set("role", "img")
    root.set("aria-labelledby", "figure-title figure-description")
    for tag, name, value in (("title", "figure-title",
                              "Goldschtein bat figure"),
                             ("desc", "figure-description", CAPTION)):
        el = ET.SubElement(root, ns + tag, {"id": name})
        el.text = value

    defs = ET.SubElement(root, ns + "defs")
    backing = ET.SubElement(root, ns + "rect", {
        "x": "0", "y": "0", "width": f"{PAGE_W:g}",
        "height": f"{PAGE_H:g}", "fill": "white"})
    backing.set("id", "page-background")

    photo = ET.SubElement(root, ns + "image")
    x0, y0, x1, y1 = PANEL_B_PHOTO
    photo.set("id", "panel-b-photo")
    photo.set("x", "%.6f" % x0)
    photo.set("y", "%.6f" % y0)
    photo.set("width", "%.6f" % (x1 - x0))
    photo.set("height", "%.6f" % (y1 - y0))
    photo.set("preserveAspectRatio", "none")
    photo.set("opacity", f"{PANEL_B_PHOTO_ALPHA:.4f}")
    photo.set("{%s}href" % XLINK_NS,
              "data:image/jpeg;base64," +
              base64.b64encode(ctx["photo"].read_bytes()).decode("ascii"))

    for letter in PANEL_ORDER:
        svg_bytes, placement = parts[letter]
        part_defs, content = _split_svg(svg_bytes)
        prefix = f"panel-{letter}--"
        if part_defs is not None:
            _prefix_svg_ids(part_defs, prefix)
            for child in part_defs:
                defs.append(child)
        group = ET.SubElement(root, ns + "g", {"id": f"panel-{letter}"})
        if letter != "b":
            # the part SVG is the panel figure in 100-dpi px, y down
            # from its top-left corner; map part units (x, y) to page pt
            # (tx + k*x, ty + k*y) with ty the part's top edge on the
            # page (page pt, y down)
            tx, ty, k = placement["tx"], placement["ty"], \
                placement["scale"] * FIG_PX_TO_PT
            group.set("transform",
                      f"matrix({k:.9f} 0 0 {k:.9f} "
                      f"{tx:.6f} {ty:.6f})")
        # letter "b" is the page-sized chrome layer, embedded 1:1
        if letter in ("a", "h"):
            n_rasters = 2 if letter == "a" else 4
            img_path = ASSET_DIR / ("image_45.jpeg" if letter == "a"
                                    else "image_23.jpeg")
            _swap_panel_rasters(content, n_rasters, img_path)
        _prefix_svg_ids(group, prefix)
        for child in content:
            group.append(child)
    return ET.tostring(root, encoding="unicode") + "\n"


def compose_pdf(ctx, pdfs, placements) -> bytes:
    """Compose the page PDF: the photograph (original stream, uniform
    0.4684 SMask), the chrome layer at 1:1, and every panel at its
    measured scale; the original basemap JPEG streams are then restored
    to their image objects and the trailer is stabilised."""
    import pymupdf as fitz
    from PIL import Image

    doc = fitz.open()
    page = doc.new_page(width=PAGE_W, height=PAGE_H)

    # photograph: original stream + uniform soft mask at the reference's
    # ExtGState alpha (the reference stores no SMask for it); the rect
    # is y-down, matching the measured placement and pymupdf 2.x's
    # page-geometry convention
    x0, y0, x1, y1 = PANEL_B_PHOTO
    photo_rect = fitz.Rect(x0, y0, x1, y1)
    photo_xref = page.insert_image(photo_rect, stream=ctx["photo"].read_bytes())
    scratch = doc.new_page(width=1, height=1)
    mask = Image.new("L", (16, 16),
                     int(round(PANEL_B_PHOTO_ALPHA * 255)))
    buf = io.BytesIO()
    mask.save(buf, "PNG")
    mask_xref = scratch.insert_image(fitz.Rect(0, 0, 1, 1),
                                     stream=buf.getvalue())
    doc.xref_set_key(photo_xref, "SMask", f"{mask_xref} 0 R")
    doc.delete_page(len(doc) - 1)
    # the scratch page invalidated the page proxy (pymupdf 2.x):
    # re-resolve it before touching the page again
    page = doc[0]

    for letter in ("b", "a", "c", "d", "e", "f", "g", "h"):
        pl = placements[letter]
        if letter == "b":
            rect = fitz.Rect(0, 0, PAGE_W, PAGE_H)
        else:
            tx, ty, w, h = pl["rect"]
            rect = fitz.Rect(tx, ty, tx + w, ty + h)
        src = fitz.open(stream=pdfs[letter], filetype="pdf")
        # the rect's aspect is the part page's by construction, so the
        # stretch-to-fill is a uniform scale: no contain/centre slack;
        # rect is y-down (see _translate)
        page.show_pdf_page(rect, src, 0, keep_proportion=False)
        src.close()

    # insert_image preserves the image's own aspect ratio and centres
    # it in the rect; the reference stretches the photograph over the
    # full measured box, so rewrite the photo Do's CTM to the exact
    # placement (page content is y up: bottom = PAGE_H - y1)
    photo_name = next(n for n, r in _xobject_names(doc, page.xref).items()
                      if r == photo_xref)
    exact_do = ("q %.6f 0 0 %.6f %.6f %.6f cm /%s Do" % (
        x1 - x0, y1 - y0, x0, PAGE_H - y1, photo_name)).encode()
    for cxx in page.get_contents():
        s = doc.xref_stream(cxx)
        if b"/" + photo_name.encode() + b" Do" not in s:
            continue
        doc.update_stream(
            cxx, re.sub(
                rb"q\s+[\d.]+ 0 0 [\d.]+ [\d.]+ [\d.]+ cm\s+/"
                + photo_name.encode() + rb"\s+Do",
                exact_do, s, count=1))

    _restore_rasters(doc, page, ctx)
    meta = doc.metadata
    meta["creationDate"] = ""
    meta["modDate"] = ""
    doc.set_metadata(meta)
    data = doc.tobytes(garbage=4, deflate=True)
    doc.close()
    return _stabilize_pdf(data)


#: (page box, y up, left-bottom-right-top tuple (x0, yb, x1, yt)), asset
#: stream file) for every raster the composition must carry verbatim.
def _raster_targets(ctx):
    targets = []
    # basemaps: the panel axes boxes (the images are drawn over exactly
    # those data extents)
    for box in PANEL_A_AXES:
        x0, y0, x1, y1 = box
        targets.append(((x0, PAGE_H - y1, x1, PAGE_H - y0),
                        ASSET_DIR / "image_45.jpeg", None))
    for box in PANEL_H_AXES:
        x0, y0, x1, y1 = box
        targets.append(((x0, PAGE_H - y1, x1, PAGE_H - y0),
                        ASSET_DIR / "image_23.jpeg", None))
    # the photograph (with its uniform 0.4684 mask)
    x0, y0, x1, y1 = PANEL_B_PHOTO
    targets.append(((x0, PAGE_H - y1, x1, PAGE_H - y0),
                    ASSET_DIR / "image_279.jpeg", "uniform"))
    return targets


def _xobject_names(doc, xref) -> dict:
    """name -> resource xref of an object's /Resources/XObject."""
    kind, val = doc.xref_get_key(xref, "Resources/XObject")
    if kind == "xref":
        val = doc.xref_object(int(val.split()[0]), compressed=True)
    return {k.decode(): int(v) for k, v in
            re.findall(rb"/([A-Za-z0-9_]+) (\d+) 0 R", val.encode())}


def _form_matrix(doc, xref):
    kind, val = doc.xref_get_key(xref, "Matrix")
    if kind != "array":
        return (1.0, 0.0, 0.0, 1.0, 0.0, 0.0)
    return tuple(float(v) for v in val.strip("[]").split())


def _affine_compose(m, acc):
    """Compose affine (a b c d e f) maps: the result applies ``acc``
    first, then ``m``."""
    a, b, c, d, e, f = m
    A, B, C, D, E, F = acc
    return (a * A + c * B, b * A + d * B,
            a * C + c * D, b * C + d * D,
            a * E + c * F + e, b * E + d * F + f)


def _affine_box(t, box):
    a, b, c, d, e, f = t
    x0, y0, x1, y1 = box
    xs = [a * x + c * y + e for x in (x0, x1) for y in (y0, y1)]
    ys = [b * x + d * y + f for x in (x0, x1) for y in (y0, y1)]
    return min(xs), min(ys), max(xs), max(ys)     # y up


#: content streams may space the image Do's on newlines (q / cm / Do
#: each on their own line), so every separator is \s+
_IMAGE_RE = re.compile(
    rb"q\s+([\d.\-]+)\s+0\s+0\s+([\d.\-]+)\s+([\d.\-]+)\s"
    rb"+([\d.\-]+)\s+cm\s+/([A-Za-z0-9_]+)\s+Do\s+Q")
#: child form invocations carry either an inline cm matrix (q 1 2 3 4
#: 5 6 cm /child Do) or none (/fullpage Do); the whole "numbers + cm"
#: prefix is one optional group (a repeat capture would keep only its
#: last number), followed by the name
_CHILD_RE = re.compile(
    rb"((?:-?[\d.]+ ){5}-?[\d.]+ cm\s)?/([A-Za-z0-9_]+)\s+Do\b")
_IDENTITY = (1.0, 0.0, 0.0, 1.0, 0.0, 0.0)


def _restore_rasters(doc, page, ctx):
    """Copy each committed original JPEG stream onto the matching image
    object of the composed page (the PDF backend resamples every raster
    to the figure dpi), re-attaching the photograph's SMask.

    Panel parts are embedded through ``show_pdf_page``: each part page
    lives in a wrapper form XObject (its placement matrix) containing
    the part's content form (identity matrix), where the basemap
    image Do's sit, in part page coordinates. The walk composes the
    chain of matrices to bring every image placement into page
    coordinates (y up) and matches them against the measured
    placements.
    """
    import pymupdf as fitz

    targets = _raster_targets(ctx)    # page boxes, y up (see there)
    used = set()
    to_swap = []            # (image xref, target index)

    def match(box_yup_local, t):
        box_page = _affine_box(t, box_yup_local)   # y up
        for idx, (tbox, _file, _mask) in enumerate(targets):
            if idx in used:
                continue
            if all(abs(a - b) < 1.5 for a, b in zip(box_page, tbox)):
                used.add(idx)
                return idx
        return None

    def handle(content, name2xref, t):
        for m in _IMAGE_RE.finditer(content):
            w, h, x, y = (float(m.groups()[i]) for i in range(4))
            name = m.groups()[4].decode()
            xref = name2xref.get(name)
            idx = match((x, y, x + w, y + h), t)
            if idx is not None:
                to_swap.append((xref, idx))
            elif h < 12.0:
                continue     # narrow strip: the panel h colorbar
                             # gradient survives at its own placement
            elif xref is None:
                continue
            else:
                raise ValueError(
                    f"unexpected raster placement in the composed page: "
                    f"/{name} local box=({x:.2f}, {y:.2f}, {x + w:.2f},"
                    f" {y + h:.2f})")
        for m in _CHILD_RE.finditer(content):
            name = m.groups()[1].decode()
            xref = name2xref.get(name)
            if xref is None:
                continue
            kind, sub = doc.xref_get_key(xref, "Subtype")
            if sub != "/Form":
                continue
            t_child = t
            if m.group(1):
                # the inline cm maps into this content space first
                inline = tuple(float(v) for v in m.group(1).split()[:-1])
                t_child = _affine_compose(t_child, inline)
            # the form's own /Matrix maps its BBox space into that
            # (page content stays y up throughout the walk)
            t_child = _affine_compose(t_child, _form_matrix(doc, xref))
            handle(doc.xref_stream(xref), _xobject_names(doc, xref),
                   t_child)

    handle(b"\n".join(doc.xref_stream(cx) for cx in page.get_contents()),
           _xobject_names(doc, page.xref), _IDENTITY)

    scratch = doc.new_page(width=1, height=1)
    for xref, idx in to_swap:
        _file, mask_kind = targets[idx][1], targets[idx][2]
        src = scratch.insert_image(
            fitz.Rect(0, 0, 1, 1), stream=_file.read_bytes())
        if src != xref:
            # insert_image dedupes identical byte streams, so the
            # re-insert of a stream that already sits on the target
            # (the photo, placed from the same file in compose_pdf)
            # hands the target's own xref back; a self-copy would
            # destroy the object
            doc.xref_copy(src, xref)
        doc.xref_set_key(xref, "Interpolate", "true")
        if mask_kind == "uniform":
            from PIL import Image
            mask = Image.new("L", (16, 16),
                             int(round(PANEL_B_PHOTO_ALPHA * 255)))
            buf = io.BytesIO()
            mask.save(buf, "PNG")
            msrc = scratch.insert_image(fitz.Rect(0, 0, 1, 1),
                                        stream=buf.getvalue())
            doc.xref_set_key(xref, "SMask", f"{msrc} 0 R")

    doc.delete_page(len(doc) - 1)


def _stabilize_pdf(data: bytes) -> bytes:
    """Pin the trailer dates and /ID so the file is a content function.

    MuPDF's ``tobytes`` derives the /ID from write-time details, the
    last non-deterministic bytes of the export.  The ID region is first
    normalised to a constant, the SHA-256 of the result is computed,
    and the digest is written into both ID slots.
    """
    data = re.sub(rb"/CreationDate \(.*?\)",
                  b"/CreationDate (D:20000101000000Z)", data)
    data = re.sub(rb"/ModDate \(.*?\)", b"/ModDate (D:20000101000000Z)",
                  data)
    _ID_RE = rb"/ID\s*\[[^\]]*\]"
    data = re.sub(_ID_RE, b"/ID [<0> <0>]", data)
    digest = hashlib.sha256(data).hexdigest().encode("ascii")
    data = re.sub(_ID_RE, b"/ID [<" + digest + b"> <" + digest + b">]",
                  data)
    return data


def _html_page(svg: str) -> str:
    return (
        "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n"
        "<meta charset=\"utf-8\">\n"
        "<title>Goldschtein bat figure</title>\n"
        "<style>\n"
        "body{margin:0;background:#fff;display:flex;flex-direction:"
        "column;align-items:center;justify-content:center;min-height:"
        "100vh}\n"
        "figure{margin:2rem auto}\n"
        "img{max-width:95vw;height:auto}\n"
        "figcaption{max-width:900px;margin-top:1rem;color:#231f20;"
        "font:14px/1.5 Helvetica,Arial,sans-serif}\n"
        "</style>\n</head>\n<body>\n<figure>\n"
        "<img src=\"data:image/svg+xml;base64,"
        + base64.b64encode(svg.encode()).decode("ascii")
        + "\" alt=\"Hierarchical von Mises-Uniform mixture modeling of "
           "homing trajectories in bats.\">\n"
        f"<figcaption>{CAPTION}</figcaption>\n"
        "</figure>\n</body>\n</html>\n")


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------
def render(ctx, output_dir, *, family="Arial") -> dict:
    """Render all exports into ``output_dir``; returns provenance."""
    global _PLOTS

    import pymupdf

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    _PLOTS = sys.modules["goldshtein_bats_plots"]

    mpl.rcParams["svg.hashsalt"] = "unwrap-goldshtein-bats-publication-v1"

    # ---- panels (natural figures) + placement on the page ------------
    parts = {}
    pdfs = {}
    placements = {}
    stats = {}

    chrome_fig, chrome_stats = draw_chrome(ctx, family=family)
    parts["b"] = (_svg_bytes(chrome_fig), {"scale": 1.0})
    # the SVG save consumed the figure; render the chrome once more for
    # the PDF part
    chrome_fig2, _ = draw_chrome(ctx, family=family)
    pdfs["b"] = _pdf_bytes(chrome_fig2)
    placements["b"] = {"rect": (0.0, 0.0, PAGE_W, PAGE_H), "scale": 1.0}
    stats["chrome"] = chrome_stats

    for letter in ("a", "c", "d", "e", "f", "g", "h"):
        fn, scale, measured = PANELS[letter]
        fig = fn(ctx, family)
        # the white figure patch would paint over the page chrome
        # (the panel letters sit inside the parts' figure pages) and
        # the photograph: the authored reference is transparent outside
        # the axes, whose faces stay white
        fig.patch.set_facecolor("none")
        fig.canvas.draw()
        anchor = fig.axes[0]
        tx, ty, w, h = _translate(fig, anchor, measured, scale)
        placements[letter] = {
            "tx": tx, "ty": ty, "scale": scale,
            "rect": (tx, ty, w, h),
            "measured_box_ydown": measured,
            "natural_anchor_pt": _figure_box_pt(fig, anchor),
        }
        parts[letter] = (None, placements[letter])
        svg_buf = io.BytesIO()
        fig.savefig(svg_buf, format="svg", metadata={"Date": None})
        pdf_buf = io.BytesIO()
        fig.savefig(pdf_buf, format="pdf", metadata={"Date": None})
        plt.close(fig)
        parts[letter] = (svg_buf.getvalue(), placements[letter])
        pdfs[letter] = pdf_buf.getvalue()
        stats[letter] = {
            "figsize_pt": list(placements[letter]["rect"][2:]),
            "anchor_ydown": measured,
        }

    svg = compose_svg(ctx, parts, family)
    svg_path = output_dir / f"{FIGURE_NAME}.svg"
    svg_path.write_text(svg)

    pdf_bytes = compose_pdf(ctx, pdfs, placements)
    pdf_path = output_dir / f"{FIGURE_NAME}.pdf"
    pdf_path.write_bytes(pdf_bytes)

    png_path = output_dir / f"{FIGURE_NAME}.png"
    with pymupdf.open(stream=pdf_bytes, filetype="pdf") as doc:
        pix = doc.load_page(0).get_pixmap(
            matrix=pymupdf.Matrix(PREVIEW_SCALE, PREVIEW_SCALE))
        pix.save(str(png_path))

    html_path = output_dir / f"{FIGURE_NAME}.html"
    html_path.write_text(_html_page(svg))

    provenance = {
        "figure": FIGURE_NAME,
        "caption": CAPTION,
        "page_pt": [PAGE_W, PAGE_H],
        "layout": {
            "scales": {"a": SCALE_A, "c": SCALE_C, "d": SCALE_D,
                       "e": SCALE_E, "f": SCALE_F, "g": SCALE_G,
                       "h": SCALE_H},
            "panels": {letter: {
                "scale": placements[letter]["scale"],
                "measured_box_pt_ydown":
                    list(placements[letter].get(
                        "measured_box_ydown", PANEL_B_PHOTO)),
                "placement_pt": {
                    "tx": placements[letter]["tx"] if "tx"
                    in placements[letter] else 0.0,
                    "ty": placements[letter]["ty"] if "ty"
                    in placements[letter] else 0.0,
                    "width": placements[letter]["rect"][2],
                    "height": placements[letter]["rect"][3],
                },
                "natural_anchor_pt":
                    list(placements[letter].get("natural_anchor_pt",
                                                (0, 0, 0, 0))),
            } for letter in ("a", "c", "d", "e", "f", "g", "h")},
            "panel_b_photo": {"box_pt_ydown": list(PANEL_B_PHOTO),
                              "alpha": PANEL_B_PHOTO_ALPHA},
            "panel_h_colorbar_bar_pt_ydown": list(PANEL_H_CBAR_BAR),
        },
        "stats": stats,
        "fonts": {"family": family,
                  "chrome_stack": FONT_STACK},
        "env": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "matplotlib": importlib.metadata.version("matplotlib"),
            "numpy": importlib.metadata.version("numpy"),
        },
        "exports": {
            "svg_sha256": sha256(svg_path),
            "pdf_sha256": sha256(pdf_path),
            "png_sha256": sha256(png_path),
            "html_sha256": sha256(html_path),
        },
    }
    (output_dir / f"{FIGURE_NAME}.json").write_text(
        json.dumps(provenance, indent=2) + "\n")
    return provenance


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Rebuild the Goldschtein bat publication figure "
                    "from the saved fits and pinned data (no sampling).")
    parser.add_argument("--data-root", type=Path,
                        default=ROOT / ".cache/goldshtein-bats-figure"
                        / DATA_REVISION)
    parser.add_argument("--download", action="store_true",
                        help="fetch any missing pinned unwrap-data files")
    parser.add_argument("--output-dir", type=Path,
                        default=ROOT / "assets"
                        / "goldshtein-bats-publication")
    parser.add_argument("--font-dir", type=Path,
                        help="directory of licensed font files")
    parser.add_argument("--font-family", default="Arial")
    args = parser.parse_args(argv)

    fonts = resolve_fonts(args.font_dir, args.font_family)
    root = prepare_inputs(args.data_root, download=args.download)
    assets = check_assets()
    ctx = build_context(root)

    ctx["fit_meta"]["env"] = {
        "data_repository": "jgraving/unwrap-data",
        "data_revision": DATA_REVISION,
        "inputs_sha256": INPUT_HASHES,
    }
    ctx["fit_diagnostics"] = {
        "vm_baseline": fit_diagnostics(ctx["idata_vm"]),
        "vmu_mixture_full": fit_diagnostics(ctx["idata_vmu"]),
    }

    prov = render(ctx, args.output_dir, family=args.font_family)
    prov["inputs"] = {
        "assets_sha256": {k: v["sha256"] for k, v in assets.items()},
        "data_repository": "jgraving/unwrap-data",
        "data_revision": DATA_REVISION,
        "data_inputs_sha256": INPUT_HASHES,
    }
    prov["fonts"] = {**prov["fonts"], "resolved": fonts}
    prov["fit_meta"] = ctx["fit_meta"]
    prov["fit_diagnostics"] = ctx["fit_diagnostics"]
    prov["n_obs"] = ctx["n_obs"]
    prov["n_tags"] = ctx["n_tags"]
    out_json = args.output_dir / f"{FIGURE_NAME}.json"
    out_json.write_text(json.dumps(prov, indent=2) + "\n")

    for name, h in prov["exports"].items():
        print(f"{name}  {h[:16]}...")
    print(f"Wrote {FIGURE_NAME}.{{svg,pdf,png,html,json}} to "
          f"{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
