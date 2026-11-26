"""Shared pinned inputs and drawing helpers for the bat homing figure.

The publication composite (``goldshtein_bats_layout_review.py``) builds on
this module for its pinned data/fit loading (:func:`build_context`,
:func:`prepare_inputs`), the curated Inkscape artwork of the treatment-group
panel (:func:`load_artwork`, :func:`panel_b_drawings`, :func:`draw_drawing`),
font resolution and the pinned input hashes recorded in the render
provenance.  Everything the composite needs comes from the pinned unwrap-data
snapshot (``DATA_REVISION``/``INPUT_HASHES``) and the notebook's own
preprocess/posthoc/plots modules."""

from __future__ import annotations

import contextlib
import hashlib
import importlib.metadata
import json
import os
import re
import sys
import urllib.request
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
from matplotlib.path import Path as MplPath  # noqa: E402
from matplotlib.patches import PathPatch  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
ASSET_DIR = Path(__file__).resolve().parent / "assets" / \
    "goldshtein_bats_figure"
REFERENCE_PDF = ASSET_DIR / "reference.pdf"
ARTWORK_PATH = ASSET_DIR / "artwork.json"

#: The reference original's page size (points); panel-b artwork boxes are
#: measured against it (see the layout note below).
PAGE_W = 583.0          # points
PAGE_H = 626.0          # points

#: jgraving/unwrap-data revision carrying the bats data and the two
#: cached fits (the same revision the zebra figure pins).
DATA_REVISION = "86c41b5f43eb61191405b9bc85baf730a492a06f"

#: SHA-256 pins for every data/fit input (see the asset README).
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
    "basemap_osm.geojson.gz":
        "f8eecc4a968875ada1dace186651d49239aeb94955a832b453bd63f16c4392d7",
    "image_279.jpeg":
        "6ffabea4e91fc7e53c7a97a97cb651642c8f8c392b1e4add3ae2fc7fad269734",
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
# Measured layout of the reference original (page points, measured from
# notebooks/assets/goldshtein_bats_figure/reference.pdf by the Inkscape pass).
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
        "panel_b_drawings": panel_b_drawings(asset),
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


