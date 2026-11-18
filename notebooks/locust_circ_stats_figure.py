"""Programmatic reconstruction of the locust circular-statistics composite
figure (panels a-d).

Panels:
  a. Virtual-swarm schematic, rebuilt from the published artwork geometry
     (``notebooks/assets/locust_circ_stats/``: the virtual-swarm sprite PNG
     plus the image placements, vector paths and text spans extracted from
     the published PDF, replayed through matplotlib).
  b. LOO ELPD comparison - the house ``plot_loo_elpd`` panel, drawn onto the
     composite page by the house function itself (see ``_house_panel``).
  c. Population-level predictive check - the house
     ``locust_posthoc.population_validation_stats`` estimates drawn in the
     publication style (Model vs Data only, plus the optic-flow and
     vectorial-encoding arrows).
  d. Hypothesis testing - the house ``plot_hypothesis_density`` panel,
     drawn onto the composite page.

Data and fits: trajectories from ``unwrap-data`` through
``unwrap.datasets.load_directional`` + ``rediscretize_directional``; the
bimodal and unimodal fits load from the published fit cache in the same
repository (``fit_cache/locust_direction_{bimodal,unimodal}.nc``). The
models are rebuilt only to attach the per-observation log-likelihood the
LOO/ELPD post-hoc reads (PyMC 6 replacement for ``idata_kwargs=
{"log_likelihood": True}``), mirroring notebook cells 10 and 14.

The page is 900 x 600 pt (12.5 x 8.33 in), drawn in PDF points through a
full-page stage axes; positions/colors/sizes reproduce the published
layout. Saves to ``figures/locust/locust_circ_stats_figure.pdf``.
"""
from __future__ import annotations

import json
import math
import os
import pathlib
import re
import sys

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.patches import FancyArrowPatch, PathPatch
from matplotlib.path import Path

import numpy as np
import pandas as pd
import pymc as pm
import pytensor.tensor as pt
import arviz as az
import PIL.Image
import PIL.ImageChops

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import unwrap  # noqa: E402
from unwrap.datasets import load_directional  # noqa: E402
from unwrap.link.pymc import circular_modulo  # noqa: E402
from unwrap.preprocess import rediscretize_directional  # noqa: E402

import fit_cache  # noqa: E402
import locust_posthoc as posthoc  # noqa: E402
import locust_direction_plots as plots  # noqa: E402

ASSETS = HERE / "assets" / "locust_circ_stats"
PAGE_W, PAGE_H = 900.0, 600.0  # PDF points, origin top-left

# ---------------------------------------------------------------------------
# layout and type (measured from the published PDF, in pt)
# ---------------------------------------------------------------------------
PANEL_B = (146.2, 401.8, 284.8, 98.0)   # axes rectangle (x, y_top, w, h)
PANEL_D = (577.2, 401.8, 296.0, 97.9)
POLAR_CENTER = (615.1, 206.3)           # polar center
POLAR_R_OUT = 118.7                     # outer grid circle radius (r = 1.0)

#: The published page carries heavier type and lines than a standalone house
#: panel: tick caps measure 12.2 pt tall (Arial cap height 0.716 em -> ~17 pt),
#: spines/ticks/error bars 2.09 pt, ticks 6.3 pt long. Scaling the house
#: rcParams (8 pt / 1.0 pt) by those measured numbers is what makes the
#: house panels sit on the same visual scale as panel a.
TYPE_SZ = 17.0
LEGEND_SZ = 15.0
LINE_W = 2.09
TICK_LEN = 6.3
MARK_SCALE = 2.09   # marker/cap/elinewidth factor in the published dots

#: axis ticks as printed in the published figure. The pads reproduce the
#: published margins from the plotted statistics: panel b's whisker extents
#: sit 10 ELPD inside the axis (published left/right edges -3549.6/-3355.7
#: against whiskers at -3540.8/-3365.2), panel d's 0.0115 density units
#: inside (0.0462/0.3152 against 0.0576/0.3034).
PANEL_B_XTICKS = [-3500.0, -3450.0, -3400.0]
PANEL_B_PAD = 10.0
PANEL_D_XTICKS = [0.1, 0.2, 0.3]
PANEL_D_PAD = 0.0115

#: dash rhythm of panel d's dashed "Uniform" reference line, measured off the
#: published PDF (7.74 pt on, 3.35 pt off at a 2.09 pt stroke)
UNIFORM_DASH = (7.74302, 3.3483332)


# ---------------------------------------------------------------------------
# data + cached fits
# ---------------------------------------------------------------------------
def load_data():
    raw = load_directional()
    df = rediscretize_directional(raw)
    df["animal_id"] = pd.factorize(df["animal_id"], sort=True)[0]
    return df


def build_models(df):
    """Rebuild the two model graphs (mirrors notebook cells 10 and 14) so
    the cached draws can be scored for log-likelihood."""
    directions = df["direction"].values
    individuals = df["animal_id"].values.astype(np.int32)
    n_ind = df["animal_id"].nunique()

    def inverse_softplus(x):
        return pt.log(pt.exp(x) - 1)

    def sort(x, sort_idx):
        return pt.take(x, sort_idx, axis=-1)

    with pm.Model() as model_bimodal:
        mu_pop_unsorted = pm.VonMises("mu_pop_unsorted", mu=0, kappa=1e-6, shape=(2,))
        sort_idx = pt.argsort(circular_modulo(mu_pop_unsorted))
        mu_pop = pm.Deterministic("mu_pop", sort(mu_pop_unsorted, sort_idx))
        mu_kappa_unsorted = pm.LogNormal("mu_kappa_unsorted", mu=np.log(30.0), sigma=0.1, shape=(2,))
        mu_kappa = pm.Deterministic("mu_kappa", sort(mu_kappa_unsorted, sort_idx))
        mu_unsorted = (pm.VonMises("mu_unsorted", mu=0, kappa=mu_kappa_unsorted,
                                   shape=(n_ind, 2)) + mu_pop_unsorted)
        mu = pm.Deterministic("mu", circular_modulo(sort(mu_unsorted, sort_idx)))
        kappa_log_pop_unsorted = pm.Normal("kappa_log_pop_unsorted", mu=inverse_softplus(2.0),
                                           sigma=0.5, shape=(2,))
        kappa_pop = pm.Deterministic("kappa_pop", sort(pt.softplus(kappa_log_pop_unsorted), sort_idx))
        kappa_sigma_unsorted = pm.HalfNormal("kappa_sigma_unsorted", sigma=1.5, shape=(2,))
        kappa_sigma = pm.Deterministic("kappa_sigma", sort(kappa_sigma_unsorted, sort_idx))
        kappa_log_unsorted = (pm.Normal("kappa_log_unsorted", mu=0, sigma=1.0, shape=(n_ind, 2))
                              * kappa_sigma_unsorted + kappa_log_pop_unsorted)
        kappa = pm.Deterministic("kappa", sort(pt.softplus(kappa_log_unsorted), sort_idx))
        weights_logit_pop = pm.Normal("weights_logit_pop", mu=0, sigma=0.5)
        weights_pop_unsorted = pm.math.sigmoid(pt.stack([weights_logit_pop, -weights_logit_pop]))
        weights_pop = pm.Deterministic("weights_pop", sort(weights_pop_unsorted, sort_idx))
        logit_sigma = pm.HalfNormal("logit_sigma", sigma=3.0)
        weights_logit = (pm.Normal("weights_logit_unsorted", mu=0, sigma=1.0, shape=(n_ind,))
                         * logit_sigma + weights_logit_pop)
        weights_unsorted = pm.math.sigmoid(pt.stack([weights_logit, -weights_logit], axis=-1))
        weights = pm.Deterministic("weights", sort(weights_unsorted, sort_idx))
        components = pm.VonMises.dist(mu=mu[individuals], kappa=kappa[individuals])
        pm.Mixture("directions_obs", w=weights[individuals], comp_dists=components,
                   observed=directions)

    with pm.Model() as model_unimodal:
        mu_pop = pm.VonMises("mu_pop", mu=0.0, kappa=1e-6)
        mu_kappa = pm.LogNormal("mu_kappa", mu=np.log(30.0), sigma=0.1)
        mu_raw = pm.VonMises("mu_raw", mu=0.0, kappa=mu_kappa, shape=n_ind)
        mu = pm.Deterministic("mu", circular_modulo(mu_raw + mu_pop))
        kappa_log_pop = pm.Normal("kappa_log_pop", mu=inverse_softplus(2.0), sigma=0.5)
        kappa_sigma = pm.HalfNormal("kappa_sigma", sigma=1.5)
        kappa_log = (pm.Normal("kappa_log", mu=0.0, sigma=1.0, shape=n_ind)
                     * kappa_sigma + kappa_log_pop)
        kappa = pm.Deterministic("kappa", pt.softplus(kappa_log))
        pm.VonMises("directions_obs", mu=mu[individuals], kappa=kappa[individuals],
                    observed=directions)

    return model_bimodal, model_unimodal


def load_fits(model_bimodal, model_unimodal):
    """Load the published bimodal/unimodal fit snapshots from unwrap-data
    and attach the per-observation log-likelihood."""
    data_root = pathlib.Path(os.environ.get("UNWRAP_DATA_ROOT",
                                            str(pathlib.Path.home() / ".cache" / "unwrap-data")))
    os.environ["UNWRAP_PY_FIT_CACHE"] = "1"
    os.environ["UNWRAP_PY_FIT_CACHE_DIR"] = str(data_root / "unwrap-data" / "fit_cache")

    out = []
    for name, model in (("locust_direction_bimodal", model_bimodal),
                        ("locust_direction_unimodal", model_unimodal)):
        snap = fit_cache.load_snapshot(name)
        idata = az.from_netcdf(snap) if not isinstance(snap, az.InferenceData) else snap
        with model:
            idata = pm.compute_log_likelihood(idata, progressbar=False)
        out.append(idata)
    return out


# ---------------------------------------------------------------------------
# panel a: replay the published artwork geometry
# ---------------------------------------------------------------------------
def _path_from_items(items):
    """Turn extracted PDF path items into a matplotlib Path.

    PDF path construction is continuous: every segment starts where the
    previous one ended and only a discontinuity (or a ``re``/``qu`` box)
    opens a new subpath. Emitting a fresh ``MOVETO`` per segment would
    break closed outlines (the black locust, the arrow heads) into
    non-fillable slivers.
    """
    verts, codes = [], []
    cur = None

    def move_to(x, y):
        if cur is None or abs(x - cur[0]) > 0.01 or abs(y - cur[1]) > 0.01:
            verts.append((x, y))
            codes.append(Path.MOVETO)

    for item in items:
        kind = item[0]
        if kind == "l":
            (x0, y0), (x1, y1) = item[1], item[2]
            move_to(x0, y0)
            verts.append((x1, y1))
            codes.append(Path.LINETO)
            cur = (x1, y1)
        elif kind == "c":
            (x0, y0), (x1, y1), (x2, y2), (x3, y3) = item[1:5]
            move_to(x0, y0)
            verts += [(x1, y1), (x2, y2), (x3, y3)]
            codes += [Path.CURVE4, Path.CURVE4, Path.CURVE4]
            cur = (x3, y3)
        elif kind == "re":
            x0, y0, x1, y1 = item[1]
            verts += [(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)]
            codes += [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO,
                      Path.CLOSEPOLY]
            cur = None
        elif kind == "qu":
            ul, ur, ll, lr = item[1:5]
            verts += [tuple(ul), tuple(ur), tuple(lr), tuple(ll), tuple(ul)]
            codes += [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO,
                      Path.CLOSEPOLY]
            cur = None
    return Path(verts, codes) if verts else None


#: resolution of the baked photo layer, as a multiple of the sprite's own
#: pixel grid (2 -> ~450 dpi at print size). The windows are stitched in
#: memory so panel a is one image with no clip rectangles.
PHOTO_OVERSAMPLE = 2.0


#: the photo background is this bright in every channel; the published panel
#: clips each locust to its silhouette, so these pixels are made transparent
#: here (the arrows and rules behind the photos show through them).
PHOTO_BG_FLOOR = 250


def _photo_alpha(sprite, floor=PHOTO_BG_FLOOR):
    """Alpha for the sprite: 0 over the photographic background.

    The published figure clips every placement to the locust's outline, so
    the grey flight arrows read through the gaps around the insect. Tracing
    that outline from the sprite itself is what the white background makes
    easy: the background is one bright connected region reaching the sprite
    border, so the region touching the border becomes transparent and
    anything enclosed (the pale markings on the body) stays opaque.
    """
    from scipy import ndimage

    a = np.asarray(sprite.convert("L"))
    near_white = a >= floor
    labels, _ = ndimage.label(near_white)
    border = set(labels[0, :]) | set(labels[-1, :]) | set(labels[:, 0]) | set(labels[:, -1])
    border.discard(0)
    background = np.isin(labels, list(border))
    return np.where(background, 0, 255).astype("uint8")


def _photo_layer(spec, oversample=PHOTO_OVERSAMPLE):
    """Stitch the locust photos into one RGBA layer.

    The published panel places twelve clipped windows of a single
    virtual-swarm sprite, each window further clipped to the locust inside
    it. Compositing them here (sub-pixel positioned, transparent outside the
    insect) gives one image to draw, which reproduces the print exactly:
    matplotlib's own per-image clipping leaves a 1 px seam along every
    window edge when the PDF is rendered.

    Returns ``(layer, rect)`` with rect the layer's page rectangle
    ``(left, top, right, bottom)`` in points, top-left origin.
    """
    sprite = PIL.Image.open(ASSETS / spec["sprite"]).convert("RGB")
    alpha = PIL.Image.fromarray(_photo_alpha(sprite), mode="L")
    w_px, h_px = sprite.size

    # page-space bounds of all the windows (top-left origin)
    left = min(pl["clip"][0] for pl in spec["placements"])
    right = max(pl["clip"][0] + pl["clip"][2] for pl in spec["placements"])
    top = min(PAGE_H - (pl["clip"][1] + pl["clip"][3]) for pl in spec["placements"])
    bottom = max(PAGE_H - pl["clip"][1] for pl in spec["placements"])
    scale = (w_px / spec["placements"][0]["cm"][0]) * oversample   # px per pt
    size = (math.ceil((right - left) * scale), math.ceil((bottom - top) * scale))

    layer = PIL.Image.new("RGBA", size, (0, 0, 0, 0))
    for pl in spec["placements"]:
        a, _, _, d, tx, ty = pl["cm"]           # PDF user space, origin bottom-left
        cx, cy, cw, ch = pl["clip"]             # clip rectangle, same space
        # canvas pixel -> sprite pixel (the affine maps output onto input)
        sx = w_px / (a * scale)
        bx = (left - tx) * w_px / a
        sy = h_px / (d * scale)
        by = h_px * (1.0 - (PAGE_H - ty - top) / d)
        coeffs = (sx, 0, bx, 0, sy, by)
        src = sprite.transform(size, PIL.Image.AFFINE, coeffs, resample=PIL.Image.BICUBIC)
        window = PIL.Image.new("L", size, 0)
        window.paste(255, (round((cx - left) * scale), round((PAGE_H - cy - ch - top) * scale),
                           round((cx + cw - left) * scale), round((PAGE_H - cy - top) * scale)))
        silhouette = alpha.transform(size, PIL.Image.AFFINE, coeffs,
                                     resample=PIL.Image.NEAREST)
        layer.alpha_composite(PIL.Image.merge(
            "RGBA", (*src.split(), PIL.ImageChops.darker(window, silhouette))))
    return layer, (left, top, right, bottom)


def draw_panel_a(ax, spec):
    """Replay the published panel-a artwork: photo layer, vector paths, text."""
    layer, (left, top, right, bottom) = _photo_layer(spec)
    ax.imshow(np.asarray(layer), extent=(left, right, top, bottom),
              origin="lower", zorder=50, interpolation="lanczos")

    for i, d in enumerate(spec["paths"]):
        p = _path_from_items(d["items"])
        if p is None:
            continue
        kw = dict(zorder=10 + i)
        face = "none"
        if d["fill"] is not None:
            face = tuple(d["fill"])
            if d.get("fill_opacity") is not None:
                face += (d["fill_opacity"],)
        edge = "none"
        lw = None
        if d["stroke"] is not None:
            edge = tuple(d["stroke"])
            if d.get("stroke_opacity") is not None:
                edge += (d["stroke_opacity"],)
            lw = d.get("width") or 1.0
        patch = PathPatch(p, facecolor=face, edgecolor=edge, linewidth=lw, **kw)
        cap, join = d.get("lineCap"), d.get("lineJoin")
    for i, d in enumerate(spec["paths"]):
        p = _path_from_items(d["items"])
        if p is None:
            continue
        kw = dict(zorder=10 + i)
        face = "none"
        if d["fill"] is not None:
            face = tuple(d["fill"])
            if d.get("fill_opacity") is not None:
                face += (d["fill_opacity"],)
        edge = "none"
        lw = None
        if d["stroke"] is not None:
            edge = tuple(d["stroke"])
            if d.get("stroke_opacity") is not None:
                edge += (d["stroke_opacity"],)
            lw = d.get("width") or 1.0
        patch = PathPatch(p, facecolor=face, edgecolor=edge, linewidth=lw, **kw)
        cap, join = d.get("lineCap"), d.get("lineJoin")
        cap = cap[0] if isinstance(cap, list) else cap
        join = join[0] if isinstance(join, list) else join
        if cap is not None:
            patch.set_capstyle(["butt", "round", "projecting"][int(cap)])
        if join is not None:
            patch.set_joinstyle(["miter", "round", "bevel"][int(join)])
        if d.get("dashes"):
            seq = [float(v) for v in re.findall(r"[-\d.]+", d["dashes"])]
            if any(seq):
                patch.set_linestyle((0, seq))
        ax.add_patch(patch)

    for s in spec["spans"]:
        color = f"#{s['color']:06x}"
        text = s["text"].replace("\u2028", "\n")
        ax.text(s["bbox"][0], s["bbox"][1], text,
                fontsize=s["size"], color=color, va="top", ha="left",
                family="Arial",
                fontweight="bold" if "Bold" in s["font"] else "normal",
                zorder=100)


# ---------------------------------------------------------------------------
# panel c: population predictive check in the publication style
# ---------------------------------------------------------------------------
# strokes measured from the published PDF (pt): model line 2.9, data line 2.9
# dashed 10.74/4.64, bands at alpha 0.25, grid circles 0.69 grey at 1.55, the
# r = 0.5 reference circle black at 0.97 and the optic-flow / vectorial-
# encoding arrows 4.0 long reaching r = 0.30.
GRID_GREY = "0.69"
GRID_LW = 1.55
MODEL_LW = 2.9
DATA_LW = 2.9
DATA_DASHES = (10.736325 / DATA_LW, 4.6427357 / DATA_LW)   # matplotlib scales by lw
BAND_ALPHA = 0.25
EDGE_CIRCLE_LW = 0.97
ARROW_LW = 4.0
ARROW_R = 0.469         # FancyArrowPatch trims the shaft by the head, so the
ARROW_HEAD = 42.0       # path end sits beyond the published tip at r = 0.423;
                        # 42 gives the measured 16.8 pt head at lw 4.0

#: legend box measured from the published PDF: 97.9 x 48.5 pt, the handle band
#: 30.95 pt across (the house handler draws the line 10 % longer), rows 21.9 pt
#: apart. matplotlib's legend bbox starts half a font size inside the anchor.
LEGEND_BOX = (741.44, 63.42, 839.35, 111.89)
LEGEND_HANDLE_LEN = 30.95 / LEGEND_SZ
LEGEND_TEXT_PAD = 0.902          # 13.5 pt gap between handle and label
LEGEND_BORDERPAD = 0.31          # 4.65 pt inside the frame
LEGEND_LABELSPACING = 0.52
LEGEND_ANCHOR_PAD = 0.5 * LEGEND_SZ
LEGEND_FRAME_LW = 1.93
LEGEND_FRAME_COLOUR = "0.8"


def draw_population_polar(fig, trace_bimodal, df):
    """Polar population predictive check: the same estimates as the house
    ``plot_directional_validation_polar`` (``locust_posthoc
    .population_validation_stats``) drawn in the publication style - Model
    and Data bands only, the two hypothesis arrows, and the page's grid."""
    stats = posthoc.population_validation_stats(trace_bimodal, df)
    angles = stats["angles"]
    density_median, density_lower, density_upper = (
        stats["density_median"], stats["density_lower"], stats["density_upper"])
    kde_angles, kde_median = stats["kde_angles"], stats["kde_median"]
    kde_lower_bounds, kde_upper_bounds = stats["kde_lower_bounds"], stats["kde_upper_bounds"]

    edge_radius = 0.5
    ylim_max = edge_radius + np.max(density_upper) * 1.2

    # Square box centred on the published polar centre; the size is then
    # calibrated against the axes transform so the r = 1.0 grid circle lands
    # exactly on POLAR_R_OUT.
    ax = fig.add_axes([0, 0, 1, 1], projection="polar")
    ax.set_ylim(0, ylim_max)
    ax.set_xlim(-np.pi, np.pi)
    size = 2.0 * POLAR_R_OUT * max(ylim_max, 1.0)
    _place(ax, (POLAR_CENTER[0] - size / 2, POLAR_CENTER[1] - size / 2, size, size))
    centre = np.asarray(ax.transData.transform((0.0, 0.0)), dtype=float)
    rim = np.asarray(ax.transData.transform((0.0, 1.0)), dtype=float)
    r_pt = np.hypot(*(rim - centre)) / fig.dpi * 72.0
    size *= POLAR_R_OUT / r_pt
    _place(ax, (POLAR_CENTER[0] - size / 2, POLAR_CENTER[1] - size / 2, size, size))

    line_model, = ax.plot(angles, edge_radius + density_median, color=plots.MODEL_RED,
                          linewidth=MODEL_LW, label="Model", zorder=4)
    ax.fill_between(angles, edge_radius + density_lower, edge_radius + density_upper,
                    color=plots.MODEL_RED, alpha=BAND_ALPHA, edgecolor="none", zorder=2)
    line_data, = ax.plot(kde_angles, edge_radius + kde_median, color="black",
                         linestyle="--", linewidth=DATA_LW, dashes=DATA_DASHES,
                         label="Data", zorder=5)
    ax.fill_between(kde_angles, edge_radius + kde_lower_bounds,
                    edge_radius + kde_upper_bounds, color="black", alpha=BAND_ALPHA,
                    edgecolor="none", zorder=3)

    ax.grid(color=GRID_GREY, linewidth=GRID_LW)
    ax.add_artist(PathPatch(Path.circle((0, 0), edge_radius),
                            transform=ax.transData._b, fill=False,
                            edgecolor="black", linewidth=EDGE_CIRCLE_LW))
    ax.spines["polar"].set_visible(False)
    ax.set_yticks([edge_radius, edge_radius + 0.5])
    ax.set_yticklabels(["", ""])
    ax.set_xticks([-np.pi / 2, 0, np.pi / 2, np.pi])
    ax.set_xticklabels(["−90°", "0°", "90°", "±180°"], fontsize=TYPE_SZ)
    ax.tick_params(axis="x", pad=-4)

    # optic-flow (forward, +x) and vectorial-encoding (sideways, ±y) arrows
    for theta in (0.0, np.pi / 2, -np.pi / 2):
        color = plots.UNIMODAL if theta == 0.0 else plots.BIMODAL
        ax.add_patch(FancyArrowPatch((theta, 0.0), (theta, ARROW_R),
                                     transform=ax.transData,
                                     arrowstyle="-|>", mutation_scale=ARROW_HEAD,
                                     color=color, linewidth=ARROW_LW,
                                     zorder=6, shrinkA=0, shrinkB=0))

    leg = ax.legend(
        [line_model, line_data], ["Model", "Data"],
        handler_map={
            line_model: plots._density_legend_handler(plots.MODEL_RED, "-", BAND_ALPHA, MODEL_LW),
            line_data: plots._density_legend_handler("black", "--", BAND_ALPHA, DATA_LW),
        },
        loc="upper left", frameon=True, fontsize=LEGEND_SZ,
        handlelength=LEGEND_HANDLE_LEN, borderpad=LEGEND_BORDERPAD,
        handletextpad=LEGEND_TEXT_PAD, labelspacing=LEGEND_LABELSPACING,
        bbox_to_anchor=((LEGEND_BOX[0] - LEGEND_ANCHOR_PAD) / PAGE_W,
                        (PAGE_H - LEGEND_BOX[1] + LEGEND_ANCHOR_PAD) / PAGE_H),
        bbox_transform=fig.transFigure,
    )
    leg.get_frame().set_edgecolor(LEGEND_FRAME_COLOUR)
    leg.get_frame().set_linewidth(LEGEND_FRAME_LW)
    return ax


# ---------------------------------------------------------------------------
# composite
# ---------------------------------------------------------------------------
def _place(ax, rect):
    """Position an axes with a PDF-point rectangle (x, y_top, w, h)."""
    x, y, w, h = rect
    ax.set_position([x / PAGE_W, (PAGE_H - y - h) / PAGE_H, w / PAGE_W, h / PAGE_H])


class _house_panel:
    """Run a house plotting function so it draws onto ``fig`` at ``rect``.

    The house panels make a figure of their own (``plt.subplots(figsize=…,
    dpi=…)``), then ``plt.tight_layout()``, and set their type and stroke
    sizes for a standalone 3-inch panel. Inside this context those calls are
    redirected so the same drawing code lands on the composite page:

    - ``plt.subplots`` adds one axes to the composite at ``rect`` (PDF
      points, top-left origin) and ``plt.tight_layout`` is a no-op, so the
      page fixes the layout;
    - the module type size ``locust_direction_plots.FONT_SZ`` (used for the
      in-axes "Uniform" label) becomes the page's 17 pt;
    - dot / error-bar / reference-line strokes are scaled by ``mark_scale``
      to the published weights.

    ``save=False`` in the call keeps the house function from writing its own
    file.
    """

    def __init__(self, fig, rect, mark_scale=MARK_SCALE, dash=None):
        self.fig, self.rect, self.mark_scale, self.dash = fig, rect, mark_scale, dash
        self.ax = None

    def __enter__(self):
        fig, rect, scale = self.fig, self.rect, self.mark_scale

        def subplots(*args, **kwargs):
            nrows = args[0] if args else 1
            ncols = args[1] if len(args) > 1 else 1
            if nrows * ncols != 1:
                raise NotImplementedError(
                    f"_house_panel places one axes; {nrows}x{ncols} requested")
            projection = (kwargs.get("subplot_kw") or {}).get("projection")
            ax = fig.add_axes([0, 0, 1, 1], projection=projection)
            _place(ax, rect)
            ax.errorbar = _scaled(ax.errorbar, scale,
                                  ("markersize", "capsize", "elinewidth",
                                   "markeredgewidth"))
            ax.axvline = _styled_vline(ax.axvline, scale, self.dash)
            self.ax = ax
            return fig, ax

        self._saved = plt.subplots, plt.tight_layout, plots.FONT_SZ
        plt.subplots = subplots
        plt.tight_layout = lambda *a, **k: None
        plots.FONT_SZ = TYPE_SZ
        return self

    def __exit__(self, *exc):
        plt.subplots, plt.tight_layout, plots.FONT_SZ = self._saved


def _scaled(method, factor, keys):
    """Wrap a drawing call so the named style kwargs are scaled up."""
    def wrapper(*args, **kwargs):
        for key in keys:
            if kwargs.get(key) is not None:
                kwargs[key] = kwargs[key] * factor
        return method(*args, **kwargs)
    return wrapper


def _styled_vline(method, factor, dash):
    """Wrap ``axvline``: scale the line width, and use the published dash
    rhythm (``dash`` is a period in points) for any dashed reference line."""
    def wrapper(*args, **kwargs):
        if kwargs.get("lw") is not None:
            kwargs["lw"] *= factor
        dashed = kwargs.pop("ls", None) or kwargs.pop("linestyle", None)
        if dashed and dash is not None:
            kwargs["linestyle"] = (0, dash)
        return method(*args, **kwargs)
    return wrapper


def main():
    plots.set_style()
    # the page is far larger than a standalone house panel, so the house
    # 8 pt / 1.0 pt style moves to the published scale measured above
    mpl.rcParams.update({
        "font.size": TYPE_SZ,
        "axes.labelsize": TYPE_SZ,
        "xtick.labelsize": TYPE_SZ,
        "ytick.labelsize": TYPE_SZ,
        "legend.fontsize": LEGEND_SZ,
        "axes.linewidth": LINE_W,
        "xtick.major.width": LINE_W,
        "ytick.major.width": LINE_W,
        "xtick.major.size": TICK_LEN,
        "ytick.major.size": TICK_LEN,
        "lines.linewidth": LINE_W,
    })

    df = load_data()
    model_bimodal, model_unimodal = build_models(df)
    trace_bimodal, trace_unimodal = load_fits(model_bimodal, model_unimodal)

    fig = Figure(figsize=(PAGE_W / 72.0, PAGE_H / 72.0), dpi=72)
    stage = fig.add_axes([0, 0, 1, 1], zorder=0)
    stage.set_xlim(0, PAGE_W)
    stage.set_ylim(PAGE_H, 0)  # top-left origin
    stage.set_axis_off()

    # panel a from the published artwork geometry
    spec = json.loads((ASSETS / "panel_a.json").read_text())
    draw_panel_a(stage, spec)

    # panels b and d: the house functions, drawn straight onto the page
    with _house_panel(fig, PANEL_B) as panel_b:
        models = {"Bimodal": trace_bimodal, "Unimodal": trace_unimodal}
        plots.plot_loo_elpd(models, save=False)
    rows = posthoc.loo_elpd_rows(models)
    lo = min(r[1] - r[2] for r in rows)
    hi = max(r[1] + r[2] for r in rows)
    panel_b.ax.set_xlim(lo - PANEL_B_PAD, hi + PANEL_B_PAD)
    panel_b.ax.set_xticks(PANEL_B_XTICKS)

    with _house_panel(fig, PANEL_D, dash=UNIFORM_DASH) as panel_d:
        plots.plot_hypothesis_density(trace_bimodal, save=False)
    stats = posthoc.hypothesis_density_stats(trace_bimodal)
    lo_d = float(np.min(stats["centers"] - stats["xerr"][:, 0]))
    hi_d = float(np.max(stats["centers"] + stats["xerr"][:, 1]))
    panel_d.ax.set_xlim(lo_d - PANEL_D_PAD, hi_d + PANEL_D_PAD)
    panel_d.ax.set_xticks(PANEL_D_XTICKS)

    # panel c
    ax_c = draw_population_polar(fig, trace_bimodal, df)
    ax_c.set_zorder(5)
    stage.set_zorder(0)

    out_dir = plots.FIG_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "locust_circ_stats_figure.pdf")
    fig.savefig(out_dir / "locust_circ_stats_figure.png", dpi=300)
    print(f"saved {out_dir / 'locust_circ_stats_figure.pdf'}")
    return fig


if __name__ == "__main__":
    main()
