"""Programmatic reconstruction of the locust circular-statistics composite
figure, extended with a third row and with the panel letters following the
figure's own flow of information: setup (a), model selection (b) and
predictive check (c), then the hypothesis test at three levels (d, e, f).
The reference page drew its four panels as a, b, c (bottom left, empty
here), d; the panels it drew as c and d are d and e here, with their
measured geometry kept exactly.

Panels:
  a. Virtual-swarm schematic, rebuilt from the reference artwork geometry
     (``notebooks/assets/locust_circ_stats/``: the virtual-swarm sprite PNG
     plus the image placements, vector paths and text spans extracted from
     the reference PDF, replayed through matplotlib).
  b. LOO ELPD comparison - the house ``plot_loo_elpd`` panel, drawn onto the
     composite page by the house function itself (see ``_house_panel``).
  c. Predictive check as a pair of polar plots, one per model: each model's
     replicated posterior-predictive direction density against the observed
     bootstrap KDE (drawn solid; dashed is the uniform reference on this
     page), on d's polar axes and radius convention, KDE estimate only. The
     house twin is ``plot_predictive_models_polar``.
  d. Hypothesis testing, predictive density - the bimodal fit's
     posterior-predictive direction density (``locust_posthoc
     .population_validation_stats``, faint, for context) with the two
     hypothesis densities of panel e marked at their directions against the
     uniform reference (the dashed circle), plus the optic-flow and
     vectorial-encoding arrows.
  e. Hypothesis testing, population level - the house ``plot_hypothesis_density``
     panel, drawn onto the composite page.
  f. The same hypothesis test as e, per animal: each animal's ELPD at the
     two hypotheses as a swarm, one dot per animal. The house twin is
     ``plot_hypothesis_elpd_swarm``.

Data and fits: trajectories from ``unwrap-data`` through
``unwrap.datasets.load_directional`` + ``rediscretize_directional``; the
bimodal and unimodal fits load from the pinned fit cache in the same
repository (``fit_cache/locust_direction_{bimodal,unimodal}.nc``). The
models are rebuilt only to attach the per-observation log-likelihood the
LOO/ELPD post-hoc reads (PyMC 6 replacement for ``idata_kwargs=
{"log_likelihood": True}``), mirroring notebook cells 10 and 14.

The page is 900 x 900 pt, drawn in PDF points through a full-page stage axes;
positions/colors/sizes reproduce the reference layout for the panels the
reference page has (a, b and the two hypothesis panels, here d and e), which
fill its 900 x 600 pt, and the predictive-check pair (c) and the individual
hypothesis panel (f) fill a third row of 300 pt below it. The composite being
matched is the manuscript's figure, so it is called the reference figure -
and its PDF the reference PDF - throughout. Saves to
``figures/locust/locust_circ_stats_figure.pdf``.

``--snapshot-dir`` (default ``assets/locust-circ-stats``) additionally writes
the review snapshots - PDF, PNG, SVG, a self-contained HTML review page and a
provenance JSON - the same set the beetle publication figure ships, so the
rendered figure can be reviewed straight from the pull request.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import os
import pathlib
import platform
import re
import subprocess
import sys

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.patches import PathPatch, Polygon
from matplotlib.path import Path

import numpy as np
import pandas as pd
import pymc as pm
import pytensor.tensor as pt
import arviz as az
import PIL.Image
import PIL.ImageChops

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(HERE))

import unwrap  # noqa: E402
from unwrap.datasets import load_directional  # noqa: E402
from unwrap.link.pymc import circular_modulo  # noqa: E402
from unwrap.preprocess import rediscretize_directional  # noqa: E402

import fit_cache  # noqa: E402
import locust_posthoc as posthoc  # noqa: E402
import locust_direction_plots as plots  # noqa: E402

ASSETS = HERE / "assets" / "locust_circ_stats"
PAGE_W, PAGE_H = 900.0, 900.0  # PDF points, origin top-left (the reference's
# 900 x 600 plus a 300 pt third row; see ROW_3_TOP below)
#: The reference PDF's own height. Panel a's photo placements are recorded in
#: that page's bottom-left space, so they are flipped against it rather than
#: against this (taller) page.
REF_PAGE_H = 600.0

#: review snapshots, in the same shape as ``assets/beetles-publication``
SNAPSHOT_DIR = REPO / "assets" / "locust-circ-stats"
STEM = "locust_circ_stats_figure"
CAPTION = (
    "Locust circular statistics. (a) Virtual-swarm schematic of the free-flight "
    "arena. (b) Leave-one-out cross-validated predictive density of the bimodal "
    "and unimodal course models (+/- one standard error). (c) Predictive check, "
    "one polar plot per model: each model's replicated posterior-predictive "
    "direction density (von Mises KDE, median and 2.5/97.5 band) against the "
    "observed bootstrap KDE (solid line). (d) The hypothesis test on the "
    "predictive-density scale: the bimodal model's posterior-predictive "
    "direction density (faint, for context) and the optic-flow (forward) and "
    "vectorial-encoding (sideways) hypothesis densities marked at their "
    "directions against the uniform reference, as the dashed circle. (e) The "
    "same test in numbers: population posterior-predictive density at the "
    "hypothesized directions against the uniform reference. (f) The test of e "
    "taken down to the individual: each animal's ELPD at the hypothesized "
    "directions, one dot per animal, against the uniform reference (open "
    "markers: the pooled values)."
)

# ---------------------------------------------------------------------------
# layout and type (measured from the reference PDF, in pt)
# ---------------------------------------------------------------------------
PANEL_B = (146.2, 401.8, 284.8, 98.0)   # axes rectangle (x, y_top, w, h)
PANEL_D = (577.2, 401.8, 296.0, 97.9)
POLAR_CENTER = (615.1, 206.3)           # polar center
POLAR_R_OUT = 118.7                     # outer grid circle radius (r = 1.0)

#: The reference page is 900 x 600 pt and carries the reference figure's four
#: panels. This figure tells a longer story, so the page grows by a third row
#: and the letters follow the flow of information rather than the reference
#: page's order:
#:
#:   left column   a setup, b model selection, c predictive check (a pair of
#:                 polar plots, one per model)
#:   right column  d, e, f: the same hypothesis test at three levels -
#:                 predictive density at the hypotheses (polar), the two
#:                 hypotheses compared at population level, then at individual
#:                 level
#:
#: So the panels the reference page drew as c and d are letters d and e here,
#: and keep their measured geometry exactly. Panels a and b keep theirs. Only c
#: (a pair of polars, where the reference page had nothing) and f (panel e's
#: inset, taller, because the swarm needs the height for its stacks) are new.
ROW_3_TOP = 600.0
ROW_3_H = PAGE_H - ROW_3_TOP
#: where the reference page's second row starts (the y of its b. and d.
#: letters); panel e sits 43.1 pt below it, and panel f the same below row three
ROW_2_TOP = 358.743
PANEL_F_INSET_Y = PANEL_D[1] - ROW_2_TOP

#: The swarm's dots. ``SWARM_DOT_PT`` is the drawn diameter in points and
#: ``SWARM_PACK_GAP`` the packing clearance as a multiple of it: at gap 1.0 the
#: dots only just touch, and the densely-populated vectorial-encoding row - all
#: 52 animals within 1.2 ELPD - merges into a solid block of squares.
#: PANEL_F's height is set from the packing itself at draw time (see build):
#: the house function works in points, so any other height draws a packing it
#: did not compute.
SWARM_DOT_PT = 4.5
SWARM_PACK_GAP = 1.25
PANEL_F = (577.2, ROW_3_TOP + PANEL_F_INSET_Y, 296.0, 1.0)

#: The predictive-check pair: one polar per model, each showing that model's
#: replicated posterior-predictive heading density against the observed KDE, so
#: the reading is "does each model reproduce the data's shape", one model at a
#: time.  They are much smaller than panel d's polar (68 pt against 118.7)
#: because the two have to share the page's left half with their horizontal
#: tick labels: the left "±180°" label (43 pt wide at this type size) has to
#: start at the same 30-odd pt margin as the panel letters, the right "0°"
#: label has to stop short of panel f's row labels (x ~ 472), and the left
#: polar's "0°" has to clear the right polar's "±180°" in the middle.
POLAR_PAIR_R_OUT = 68.0
POLAR_PAIR_CENTRES = ((152.0, 759.0), (373.0, 759.0))
POLAR_PAIR_DATA_KEY = (262.5, 873.0)      # "Data" key, centred under the pair

#: Titles and panel letters.  The reference page's letters and titles come out
#: of ``panel_a.json`` in Helvetica Neue at the reference's weight; they are
#: filtered at draw time and redrawn here, because the flow (and so three of
#: the six titles) is ours.  Style as the reference: 32 pt bold letters and
#: 17 pt titles at 13.74 and 21.465 pt from the row top - but bold, matching
#: the bat composite's panel titles.
PANEL_TITLES = (
    ("Model selection", 292.8, 366.465),
    ("Predictive check", 0.5 * (POLAR_PAIR_CENTRES[0][0] + POLAR_PAIR_CENTRES[1][0]),
     ROW_3_TOP + 21.465),
    ("Hypothesis testing: predictive density", 657.9, 21.465),
    ("Hypothesis testing: population level", 716.3, 366.465),
    ("Hypothesis testing: individual level", 716.3, ROW_3_TOP + 21.465),
)
PANEL_LETTERS = (
    ("a.", 31.823, 13.743),
    ("b.", 30.948, ROW_2_TOP + 13.74),
    ("c.", 31.823, ROW_3_TOP + 13.74),
    ("d.", 459.922, 12.639),
    ("e.", 459.047, ROW_2_TOP + 13.74),
    ("f.", 459.047, ROW_3_TOP + 13.74),
)
#: the reference page's own letters and titles, filtered out of the panel-a
#: spans so only the labels inside its artwork ("Virtual swarm", "Optic flow",
#: ...) are replayed as drawn
SPAN_DROP_TEXTS = {"a.", "b.", "c.", "d.", "Model selection",
                   "Population-level predictive check", "Hypothesis testing"}
TITLE_SZ, LETTER_SZ = 17.0, 32.0
EDGE_RADIUS = 0.5                       # density zero, on the reference circle

#: The value axis the swarm is packed against.  The packing clearance is
#: SWARM_DOT_PT * SWARM_PACK_GAP (above), so the dots stay legible however the
#: panel is laid out; the axis limits only fix the scale the packing is
#: measured on, so they are passed to the house function rather than set
#: afterwards.
PANEL_F_XLIM = (-6.5, -0.5)
PANEL_F_XTICKS = [-6.0, -5.0, -4.0, -3.0, -2.0, -1.0]

#: The reference page carries heavier type and lines than a standalone house
#: panel: tick caps measure 12.2 pt tall (Arial cap height 0.716 em -> ~17 pt),
#: spines/ticks/error bars 2.09 pt, ticks 6.3 pt long. Scaling the house
#: rcParams (8 pt / 1.0 pt) by those measured numbers is what makes the
#: house panels sit on the same visual scale as panel a.
TYPE_SZ = 17.0
LEGEND_SZ = 15.0
KEY_SZ = 13.0   # the predictive-check pair's centre labels and "Data" key
LINE_W = 2.09
TICK_LEN = 6.3
MARK_SCALE = 2.09   # marker/cap/elinewidth factor in the reference figure's dots

#: axis ticks as printed in the reference figure. The pads reproduce that
#: figure's margins from the plotted statistics: panel b's whisker extents
#: sit 10 ELPD inside the axis (reference left/right edges -3549.6/-3355.7
#: against whiskers at -3540.8/-3365.2), panel d's 0.0115 density units
#: inside (0.0462/0.3152 against 0.0576/0.3034).
#: The reference page prints three ELPD ticks because its axis starts at
#: -3549.6; ours starts at -3581.7 (our unimodal whisker reaches -3571.7, see
#: the PR note), so a fourth tick keeps the spacing even across the axis.
PANEL_B_XTICKS = [-3550.0, -3500.0, -3450.0, -3400.0]
PANEL_B_PAD = 10.0
PANEL_D_XTICKS = [0.1, 0.2, 0.3]
PANEL_D_PAD = 0.0115

#: dash rhythm of panel d's dashed "Uniform" reference line, measured off the
#: reference PDF (7.74 pt on, 3.35 pt off at a 2.09 pt stroke)
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
    """Load the pinned bimodal/unimodal fit snapshots from unwrap-data
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
# panel a: replay the reference artwork geometry
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


#: the photo background is this bright in every channel; the reference panel
#: clips each locust to its silhouette, so these pixels are made transparent
#: here (the arrows and rules behind the photos show through them).
PHOTO_BG_FLOOR = 250


def _photo_alpha(sprite, floor=PHOTO_BG_FLOOR):
    """Alpha for the sprite: 0 over the photographic background.

    The reference figure clips every placement to the locust's outline, so
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

    The reference panel places twelve clipped windows of a single
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

    # page-space bounds of all the windows (top-left origin).  The placements
    # are the only part of panel_a.json still in the PDF's bottom-left space,
    # so they convert against the page the artwork was read from, not against
    # this page (which is taller, panels e and f being ours).
    left = min(pl["clip"][0] for pl in spec["placements"])
    right = max(pl["clip"][0] + pl["clip"][2] for pl in spec["placements"])
    top = min(REF_PAGE_H - (pl["clip"][1] + pl["clip"][3]) for pl in spec["placements"])
    bottom = max(REF_PAGE_H - pl["clip"][1] for pl in spec["placements"])
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
        by = h_px * (1.0 - (REF_PAGE_H - ty - top) / d)
        coeffs = (sx, 0, bx, 0, sy, by)
        src = sprite.transform(size, PIL.Image.AFFINE, coeffs, resample=PIL.Image.BICUBIC)
        window = PIL.Image.new("L", size, 0)
        window.paste(255, (round((cx - left) * scale),
                           round((REF_PAGE_H - cy - ch - top) * scale),
                           round((cx + cw - left) * scale),
                           round((REF_PAGE_H - cy - top) * scale)))
        silhouette = alpha.transform(size, PIL.Image.AFFINE, coeffs,
                                     resample=PIL.Image.NEAREST)
        layer.alpha_composite(PIL.Image.merge(
            "RGBA", (*src.split(), PIL.ImageChops.darker(window, silhouette))))
    return layer, (left, top, right, bottom)


def draw_panel_a(ax, spec):
    """Replay the reference panel-a artwork: photo layer, vector paths, text."""
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
        if s["text"] in SPAN_DROP_TEXTS:      # letters and titles are redrawn
            continue                          # by draw_panel_titles, in our flow
        text = s["text"].replace("\u2028", "\n")
        ax.text(s["bbox"][0], s["bbox"][1], text,
                fontsize=s["size"], color=color, va="top", ha="left",
                family="Arial",
                fontweight="bold" if "Bold" in s["font"] else "normal",
                zorder=100)


def draw_panel_titles(stage):
    """Panel letters and titles for the six-panel page.

    The reference page draws its own letters and titles, and its geometry is
    replayed for panels a, b and the two hypothesis panels; the texts are
    filtered out of that replay and redrawn here, because the flow of the
    figure (and so the wording of three of the titles) is ours.  Sizes and
    offsets are the reference page's - 32 pt bold letters, 17 pt titles at
    13.74 / 21.465 pt from the row top - but the titles are bold, as the bat
    composite's are."""
    for text, x, y in PANEL_TITLES:
        stage.text(x, y, text, fontsize=TITLE_SZ, family="Arial", ha="center",
                   va="top", fontweight="bold", zorder=100)
    for text, x, y in PANEL_LETTERS:
        stage.text(x, y, text, fontsize=LETTER_SZ, family="Arial",
                   fontweight="bold", ha="left", va="top", zorder=100)


# ---------------------------------------------------------------------------
# panels c and e: the polar panels, in the publication style
# ---------------------------------------------------------------------------
# strokes measured from the reference PDF (pt): model line 2.9, data line 2.9
# dashed 10.74/4.64, bands at alpha 0.25, grid circles 0.69 grey at 1.55, the
# r = 0.5 reference circle black at 0.97 and the optic-flow / vectorial-
# encoding arrows 4.0 long reaching r = 0.30.
GRID_GREY = "0.69"
GRID_LW = 1.55
MODEL_LW = 2.9
DATA_LW = 2.9
DATA_DASHES = (10.736325 / DATA_LW, 4.6427357 / DATA_LW)   # matplotlib scales by lw
BAND_ALPHA = 0.25
DATA_BAND_ALPHA = 0.12   # panels c and d: the observed band lies under its line
#: panel d's model band is context behind the hypothesis markers, so it sits
#: lighter than panel c's, where the bands are the subject.  Its dashed
#: uniform circle carries panel e's line weight (house 1.0 pt at MARK_SCALE).
MODEL_BAND_ALPHA = 0.15
MODEL_LINE_LW = 1.6
MODEL_LINE_ALPHA = 0.45
UNIFORM_LW = MARK_SCALE
#: panel d's hypothesis markers, which draw panel e's two estimates.  Panel e
#: draws them 12.5 pt across, but there one density unit is 118.7 pt against
#: panel e's 1096: at 12.5 pt the optic-flow marker's edge would come within
#: 3.5 pt of the dashed uniform circle (centre 9.8 pt away), so it would read as
#: sitting on the reference.  10 pt keeps a 4.8 pt gap there and 4.0 pt from
#: the r = 0.5 circle.
HYP_MARK_PT = 10.0
EDGE_CIRCLE_LW = 0.97
ARROW_LW = 4.0
#: The page's arrow style, measured off the reference PDF: a 4.0 pt shaft
#: closed by a triangular head whose base and length are both 16.8 pt (panel
#: a's swarm arrows are the same shape at 3.0 pt shaft / 13.2 pt head).
#: Panel d's tips sit 47.96 pt (forward) and 50.22 / 49.97 pt (the two
#: sideways heads) out from the polar centre.  These arrows are drawn as
#: explicit paths so panel d uses panel a's arrow rather than matplotlib's
#: "-|>" head, which is blunter and shorter at the same nominal size.
ARROW_HEAD = 16.8
ARROW_TIP_R = {0.0: 47.96, np.pi / 2: 50.22, -np.pi / 2: 49.97}

#: Polar tick labels, as offsets from the polar centre measured at the
#: reference radius POLAR_R_OUT (their reference page positions are
#: (617.00, 73.56), (739.28, 204.72), (615.02, 335.00) and (492.56, 204.72)
#: about (615.1, 206.3); the nearest label edges clear the r = 1.0 grid
#: circle by 3.8-8.2 pt).  A panel with a different radius scales the offsets
#: with it.  matplotlib places polar tick labels on the grid circle, which
#: would set "±180°" on top of the 180° grid line; the house polar gets the
#: clearance by padding its label strings with spaces, too blunt a tool at
#: this type size.
POLAR_TICK_LABELS = (
    (1.90, -132.74, "center", "90°"),
    (124.18, -1.58, "left", "0°"),
    (-0.08, 128.70, "center", "−90°"),
    (-122.54, -1.58, "right", "±180°"),
)


def _page_polar_axes(fig, centre, r_out, ylim_max, tick_labels=None):
    """A polar axes on the page's grid, calibrated so that the r = 1.0 grid
    circle is a circle of ``r_out`` points about ``centre``.

    A square box is placed on the target and its size then corrected against
    the axes transform: the polar projection fits its radius to the box and
    the y limit together, so the first placement is only an estimate. The
    reference's black r = 0.5 circle (density zero), grid and tick labels come
    with it; the labels are figure text, so they clear the grid the way they do
    on the reference page.  ``tick_labels`` is a (dx, dy, ha, text) sequence in
    points from the centre; it defaults to ``POLAR_TICK_LABELS`` (measured at
    the reference radius) scaled with the radius, which suits a polar of a
    similar size - a smaller panel sharing its column should pass labels
    hugging its own circle instead."""
    ax = fig.add_axes([0, 0, 1, 1], projection="polar")
    ax.set_ylim(0, ylim_max)
    ax.set_xlim(-np.pi, np.pi)
    size = 2.0 * r_out * max(ylim_max, 1.0)
    _place(ax, (centre[0] - size / 2, centre[1] - size / 2, size, size))
    centre_pt = np.asarray(ax.transData.transform((0.0, 0.0)), dtype=float)
    rim = np.asarray(ax.transData.transform((0.0, 1.0)), dtype=float)
    r_pt = np.hypot(*(rim - centre_pt)) / fig.dpi * 72.0
    size *= r_out / r_pt
    _place(ax, (centre[0] - size / 2, centre[1] - size / 2, size, size))

    ax.grid(color=GRID_GREY, linewidth=GRID_LW)
    ax.add_artist(PathPatch(Path.circle((0, 0), EDGE_RADIUS),
                            transform=ax.transData._b, fill=False,
                            edgecolor="black", linewidth=EDGE_CIRCLE_LW))
    ax.spines["polar"].set_visible(False)
    ax.set_yticks([EDGE_RADIUS, EDGE_RADIUS + 0.5])
    ax.set_yticklabels(["", ""])
    ax.set_xticks([-np.pi / 2, 0, np.pi / 2, np.pi])
    ax.set_xticklabels(["", "", "", ""])
    ax.tick_params(axis="x", pad=-4)

    if tick_labels is None:
        scale = r_out / POLAR_R_OUT
        tick_labels = [(dx * scale, dy * scale, ha, text)
                       for dx, dy, ha, text in POLAR_TICK_LABELS]
    for dx, dy, ha, label in tick_labels:
        fig.text((centre[0] + dx) / PAGE_W,
                 (PAGE_H - (centre[1] + dy)) / PAGE_H, label,
                 ha=ha, va="center", fontsize=TYPE_SZ, zorder=8)
    return ax


def draw_hypothesis_polar(fig, trace_bimodal, df):
    """Panel d: the hypothesis test on the polar axes - the bimodal fit's
    posterior-predictive direction density (``locust_posthoc
    .population_validation_stats``, faint, for context) with the two
    hypothesis densities of panel e marked at the directions they describe -
    optic flow (forward) as a blue circle at 0 degrees, vectorial encoding
    (sideways) as orange squares at both +/-90 degrees - against the uniform
    reference, which here is the dashed circle at 1/(2*pi). The markers,
    colours and shapes are panel e's, so d reads as e's polar twin; the
    intervals stay in panel e, where the scale can show them."""
    stats = posthoc.population_validation_stats(trace_bimodal, df)
    hyp = posthoc.hypothesis_density_stats(trace_bimodal)
    angles = stats["angles"]
    density_median, density_lower, density_upper = (
        stats["density_median"], stats["density_lower"], stats["density_upper"])

    ylim_max = EDGE_RADIUS + np.max(density_upper) * 1.2
    ax = _page_polar_axes(fig, POLAR_CENTER, POLAR_R_OUT, ylim_max)

    # the model's predictive density behind the hypotheses. The observed-data
    # band this panel used to carry is gone: panel c checks each model's
    # predictive density against the data, and this panel's job is where the
    # two hypotheses sit.
    ax.fill_between(angles, EDGE_RADIUS + density_lower,
                    EDGE_RADIUS + density_upper, color=plots.BIMODAL,
                    alpha=MODEL_BAND_ALPHA, edgecolor="none", zorder=2)
    ax.plot(angles, EDGE_RADIUS + density_median, color=plots.BIMODAL,
            linewidth=MODEL_LINE_LW, alpha=MODEL_LINE_ALPHA, zorder=3)

    # the uniform reference: panel e's dashed line, as the polar twin
    ax.add_artist(PathPatch(
        Path.circle((0, 0), EDGE_RADIUS + float(hyp["p_uniform"])),
        transform=ax.transData._b, fill=False, edgecolor="black",
        linewidth=UNIFORM_LW, linestyle=(0, UNIFORM_DASH), zorder=4))

    # the hypothesis densities of panel e at their directions; the
    # perpendicular hypothesis is symmetric, so it marks both sideways axes.
    # Point estimates only: one density unit is 118.7 pt here against panel
    # e's 1096, so the +/-0.03 HDIs would be 3.6 pt long - shorter than the
    # marker - and panel e is where the intervals are read.
    centers = np.asarray(hyp["centers"], dtype=float)
    for thetas, i, color, marker in (((0.0,), 0, plots.UNIMODAL, "o"),
                                     ((-np.pi / 2, np.pi / 2), 1, plots.BIMODAL, "s")):
        ax.plot(thetas, np.full(len(thetas), EDGE_RADIUS + centers[i]), marker=marker,
                linestyle="", color=color, markersize=HYP_MARK_PT,
                markeredgecolor=color, markeredgewidth=MARK_SCALE, zorder=6)

    # optic-flow (forward, +x) and vectorial-encoding (sideways, ±y) arrows,
    # drawn in the page's own arrow style (see ARROW_HEAD) so panel d reads
    # like panel a.  They go on a page-point overlay rather than on the polar
    # axes: the measured triangle heads come out exact in page coordinates,
    # and nothing needs to survive a polar transform.
    overlay = fig.add_axes([0, 0, 1, 1], zorder=7)
    overlay.set_xlim(0, PAGE_W)
    overlay.set_ylim(PAGE_H, 0)  # top-left origin, like the panel-a stage
    overlay.set_axis_off()
    overlay.patch.set_visible(False)
    cx, cy = POLAR_CENTER
    for theta, tip_r in ARROW_TIP_R.items():
        color = plots.UNIMODAL if theta == 0.0 else plots.BIMODAL
        ux, uy = math.cos(theta), -math.sin(theta)   # page y grows downward
        nx, ny = -uy, ux                             # along the head's base
        base_r = tip_r - ARROW_HEAD
        overlay.plot([cx, cx + ux * base_r], [cy, cy + uy * base_r], color=color,
                     linewidth=ARROW_LW, solid_capstyle="butt", zorder=6)
        head = [(cx + ux * tip_r, cy + uy * tip_r)]
        for side in (1.0, -1.0):
            head.append((cx + ux * base_r + side * nx * ARROW_HEAD / 2,
                         cy + uy * base_r + side * ny * ARROW_HEAD / 2))
        overlay.add_patch(Polygon(head, closed=True, facecolor=color,
                                  edgecolor="none", zorder=6))

    # no legend: the marker colours and shapes are panel e's row labels
    # ("Optic flow" / "Vec. encoding"), directly below, and the dashed circle
    # is the same reference panel e labels.
    plotted = {
        "model_density_max": float(np.max(density_median)),
        "model_density_max_angle_deg": float(
            np.degrees(angles[int(np.argmax(density_median))])),
        "hypothesis_densities": [
            {"hypothesis": label, "at_deg": at, "density": float(center),
             "hdi_low": float(center - err[0]), "hdi_high": float(center + err[1])}
            for label, at, center, err in (
                ("optic flow (forward)", 0.0, np.asarray(hyp["centers"])[0],
                 np.asarray(hyp["xerr"])[0]),
                ("vectorial encoding (perpendicular)", 90.0,
                 np.asarray(hyp["centers"])[1], np.asarray(hyp["xerr"])[1]))
        ],
        "uniform_density_reference": float(hyp["p_uniform"]),
        "note": "observed-data band omitted; panel c draws it for each model",
    }
    return ax, plotted


def _pair_tick_labels(r_out):
    """The predictive-check pair's polar tick labels, a fixed distance outside
    their own r = 1.0 grid circle.  The reference page's offsets
    (``POLAR_TICK_LABELS``) are measured at its radius; scaled down by the
    pair's smaller radius they would run the left "±180°" label off the page
    and slide the right "0°" label under panel f's row labels.  The clearances
    here are the reference page's own: its labels sit 3.8-5.5 pt beyond the
    circle horizontally and 9.5-14 pt beyond it vertically."""
    return ((0.0, -(r_out + 14.0), "center", "90°"),
            (r_out + 5.5, 0.0, "left", "0°"),
            (0.0, r_out + 9.5, "center", "−90°"),
            (-(r_out + 3.8), 0.0, "right", "±180°"))


def draw_predictive_pair(fig, stage, trace_bimodal, trace_unimodal, df):
    """Panel c: the predictive check as a pair of polar plots, one per model.

    Each polar carries one model's replicated posterior-predictive direction
    density (von Mises KDE median + 2.5/97.5 band across the shared
    replicates, ``locust_posthoc.both_models_ppc_stats``) against the observed
    bootstrap KDE, on panel d's radius convention - density zero on the
    r = 0.5 circle.  Splitting the models is the point: on a single axis the
    two bands lie so close together that the only reading is which is nearer
    the data, and the shape each model reproduces is lost.

    The observed KDE is drawn as a *solid* line, where the reference page and
    the house panels draw data dashed.  On this page dashed is the uniform
    reference - panel d's circle, panels e and f's vertical lines - and a
    dashed data line reads as a second uniform sitting next to the one it is
    being compared to.

    Drawn here rather than by the house ``plot_predictive_models_polar`` for
    the same reason as panel d: the page's grid weights, band alphas and
    tick-label positions are measured, and the house polar puts its labels on
    the grid circle."""
    stats = posthoc.both_models_ppc_stats(trace_bimodal, trace_unimodal, df)
    obs = posthoc.bootstrap_von_mises_kde_stats(df)
    angles = stats["population_angles"]

    ylim_max = EDGE_RADIUS + 1.2 * max(np.max(stats["Bimodal"]["kde_upper"]),
                                       np.max(stats["Unimodal"]["kde_upper"]),
                                       np.max(obs["kde_upper_bounds"]))

    def _peak(median, grid):
        return {"max": float(np.max(median)),
                "max_angle_deg": float(np.degrees(grid[int(np.argmax(median))]))}

    plotted = {}
    for (name, colour), centre in zip((("Bimodal", plots.BIMODAL),
                                       ("Unimodal", plots.UNIMODAL)),
                                      POLAR_PAIR_CENTRES):
        model = stats[name]
        ax = _page_polar_axes(fig, centre, POLAR_PAIR_R_OUT, ylim_max,
                              tick_labels=_pair_tick_labels(POLAR_PAIR_R_OUT))
        ax.set_zorder(5)
        ax.fill_between(angles, EDGE_RADIUS + model["kde_lower"],
                        EDGE_RADIUS + model["kde_upper"], color=colour,
                        alpha=BAND_ALPHA, edgecolor="none", zorder=2)
        ax.plot(angles, EDGE_RADIUS + model["kde_median"], color=colour,
                linewidth=MODEL_LW, zorder=5)
        # the data over the model band, solid: see the docstring
        ax.fill_between(obs["kde_angles"], EDGE_RADIUS + obs["kde_lower_bounds"],
                        EDGE_RADIUS + obs["kde_upper_bounds"], color="black",
                        alpha=DATA_BAND_ALPHA, edgecolor="none", zorder=3)
        ax.plot(obs["kde_angles"], EDGE_RADIUS + obs["kde_median"],
                color="black", linewidth=DATA_LW, zorder=6)
        # r = 0 sits below the axis limit (density zero is r = 0.5), so the
        # label is explicitly unclipped
        ax.text(0.0, 0.0, name, color=colour, fontsize=KEY_SZ, ha="center",
                va="center", zorder=7, clip_on=False,
                bbox=dict(facecolor="white", edgecolor="none", pad=1.2))
        plotted[name.lower()] = _peak(model["kde_median"], angles)

    # one key for the pair: what the solid black line is
    key_x, key_y = POLAR_PAIR_DATA_KEY
    stage.plot([key_x - 15.5, key_x + 15.5], [key_y, key_y], color="black",
               linewidth=DATA_LW, solid_capstyle="butt", zorder=6,
               clip_on=False)
    stage.text(key_x + 21.0, key_y, "Data", fontsize=KEY_SZ, family="Arial",
               ha="left", va="center", zorder=6)

    plotted["data"] = _peak(obs["kde_median"], obs["kde_angles"])
    plotted["estimate"] = "von Mises KDE, histogram estimate omitted"
    plotted["data_line"] = "solid (dashed is the uniform reference on this page)"
    return plotted


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
      to the reference weights.

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
    """Wrap ``axvline``: scale the line width, and use the reference dash
    rhythm (``dash`` is a period in points) for any dashed reference line."""
    def wrapper(*args, **kwargs):
        if kwargs.get("lw") is not None:
            kwargs["lw"] *= factor
        dashed = kwargs.pop("ls", None) or kwargs.pop("linestyle", None)
        if dashed and dash is not None:
            kwargs["linestyle"] = (0, dash)
        return method(*args, **kwargs)
    return wrapper


def build_figure():
    """Compose the page. Returns ``(fig, plotted)`` where ``plotted`` records
    the statistics the figure draws, for the provenance snapshot."""
    plots.set_style()
    # the page is far larger than a standalone house panel, so the house
    # 8 pt / 1.0 pt style moves to the reference figure's scale measured above
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

    # panel a from the reference artwork geometry
    spec = json.loads((ASSETS / "panel_a.json").read_text())
    draw_panel_a(stage, spec)
    draw_panel_titles(stage)

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

    # panel d (the reference page's c): the hypothesis test on the
    # predictive-density scale, polar, with the two hypothesis densities
    # marked at their angles
    ax_d, polar_plotted = draw_hypothesis_polar(fig, trace_bimodal, df)
    ax_d.set_zorder(5)

    # panel c: the predictive check, one polar per model against the data;
    # sets its own axes zorder
    predictive_plotted = draw_predictive_pair(fig, stage, trace_bimodal,
                                              trace_unimodal, df)
    stage.set_zorder(0)

    # panel f: the hypothesis test of panel e, one dot per animal.  The
    # uniform reference is drawn but not labelled here - panel e labels the
    # same line (on the density scale) directly above.  The axes starts 1 pt
    # tall: the house swarm sets its y limits in points from its own packing,
    # and a polar/axes in points only draws 1 unit per point when the box
    # height equals the limit span, so the panel height is read back from the
    # limits and the axes re-placed.  (The height cannot be predicted from
    # the dot size alone: the greedy packer's stack reach is not monotonic in
    # the packing distance.)
    with _house_panel(fig, PANEL_F, dash=UNIFORM_DASH) as panel_f:
        plots.plot_hypothesis_elpd_swarm(
            trace_bimodal, save=False, dot_pt=SWARM_DOT_PT,
            pack_gap=SWARM_PACK_GAP, xlim=PANEL_F_XLIM, xticks=PANEL_F_XTICKS,
            uniform_label=False)
    f_x, f_y, f_w, _ = PANEL_F
    _place(panel_f.ax, (f_x, f_y, f_w, abs(panel_f.ax.get_ylim()[0]) * 2))

    hip = posthoc.hypothesis_elpd_individual_stats(trace_bimodal)
    n_ind = int(hip["n_ind"])
    p_prefers_ve = np.asarray(hip["p_gt_all"])[:n_ind]
    individual_hypothesis = {
        "n_individuals": n_ind,
        "reference": "ELPD (log predictive density) per animal at each "
                     "hypothesised direction; dots are animals, open markers "
                     "the pooled values",
        "log_uniform_reference": float(hip["logp_uniform"]),
        "by_hypothesis": [
            {"hypothesis": label,
             "pooled_elpd": float(pooled), "pooled_se": float(se),
             "min": float(np.min(vals)), "median": float(np.median(vals)),
             "max": float(np.max(vals))}
            for label, vals, pooled, se in (
                ("optic flow (forward)", np.asarray(hip["elpd_f_all"])[:n_ind],
                 hip["elpd_f_all"][n_ind], hip["se_f_all"][n_ind]),
                ("vectorial encoding (perpendicular)",
                 np.asarray(hip["elpd_p_all"])[:n_ind], hip["elpd_p_all"][n_ind],
                 hip["se_p_all"][n_ind]))
        ],
        "p_VE_gt_OF_median": float(np.median(p_prefers_ve)),
        "frac_individuals_P_VE_gt_OF_above_half":
            float(np.mean(p_prefers_ve > 0.5)),
    }

    plotted = {
        "loo_elpd": [
            {"model": name, "elpd": value, "se": se} for name, value, se in rows
        ],
        "hypothesis_density": [
            {"hypothesis": label, "density": float(center),
             "hdi_low": float(center - err[0]), "hdi_high": float(center + err[1])}
            for label, center, err in zip(("forward", "perpendicular"),
                                          np.asarray(stats["centers"]),
                                          np.asarray(stats["xerr"]))
        ],
        "uniform_density_reference": float(stats["p_uniform"]),
        "hypothesis_polar": polar_plotted,
        "predictive_pair": predictive_plotted,
        "hypothesis_elpd_individual": individual_hypothesis,
        "observations": int(len(df)),
        "individuals": int(pd.Series(df["animal_id"]).nunique()),
    }
    return fig, plotted


# ---------------------------------------------------------------------------
# output: generated figures + committed review snapshots
# ---------------------------------------------------------------------------
def _sha256(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def _git_revision(*where):
    for root in where:
        try:
            return subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=str(root), text=True,
                stderr=subprocess.DEVNULL).strip() or None
        except (OSError, subprocess.CalledProcessError):
            continue
    return None


def _fit_spec(name):
    """The sampling settings recorded alongside the cached snapshot."""
    meta = fit_cache.cache_dir() / f"{name}.meta.json"
    if not meta.exists():
        return {"snapshot": None}
    info = json.loads(meta.read_text())
    return {"snapshot": str(fit_cache.cache_dir() / f"{name}.nc"),
            "sha256": info.get("sha256"), **info.get("extra", {})}


def provenance(plotted):
    data_root = pathlib.Path(os.environ.get("UNWRAP_DATA_ROOT",
                                            str(pathlib.Path.home() / ".cache" / "unwrap-data")))
    packages = ("numpy", "scipy", "pandas", "matplotlib", "arviz", "pymc",
                "pytensor", "xarray", "pillow")
    return {
        "figure": "locust circular-statistics composite (the reference page's "
                  "panels with the letters following this figure's information "
                  "flow: a setup, b model selection, c predictive check, d-f "
                  "hypothesis testing)",
        "page_points": [PAGE_W, PAGE_H],
        "code_revision": _git_revision(REPO),
        "data_repository": "jgraving/unwrap-data",
        "data_revision": _git_revision(data_root / "unwrap-data"),
        "font_family": mpl.rcParams["font.family"][0],
        "python": platform.python_version(),
        "packages": {name: importlib.metadata.version(name) for name in packages},
        "panel_a_source": {
            "kind": "extracted from the reference figure PDF",
            "assets": ["notebooks/assets/locust_circ_stats/panel_a.json",
                       "notebooks/assets/locust_circ_stats/virtual_swarm.png"],
            "note": "image placements, vector paths and text spans read from the "
                    "reference PDF and replayed; the sprite is the virtual-swarm "
                    "photograph inside it",
        },
        "fits": {name: _fit_spec(name) for name in
                 ("locust_direction_bimodal", "locust_direction_unimodal")},
        "plotted_values": plotted,
        "corrections": [
            "panel a records 34 path primitives: the first extraction pass kept "
            "paths above y = 300 pt only and missed the two grey flight arrows of "
            "the bottom swarm's second row (y = 301-321 pt).",
            "the panel letters and titles are not replayed from panel_a.json: "
            "the reference page's four letters and three titles are filtered out "
            "of the span replay and redrawn, because the flow (and so the "
            "wording of three of the titles) is ours.  Sizes and offsets are the "
            "reference's - 32 pt bold letters and 17 pt titles at 13.74 and "
            "21.465 pt from the row top - but the titles are bold, matching the "
            "bat composite's panel titles.",
            "panel d's arrows are drawn as a shaft plus a 16.8 pt triangular head "
            "in page points, i.e. panel a's arrow style; matplotlib's \"-|>\" head "
            "is blunter and shorter at the same nominal size.",
            "panel d's polar tick labels are placed at their reference page "
            "positions: matplotlib puts polar tick labels on the grid circle, which "
            "lands \"±180°\" on top of the 180° grid line.",
            "panel b prints four ELPD ticks where the reference page prints three, "
            "because our axis spans to -3581.7 rather than -3549.6.",
            "panel d (the reference page's c) is not the reference page's panel "
            "any more: the observed-data band is gone (that comparison is panel "
            "c's, where each model's predictive density meets the data), the "
            "bimodal predictive density stays as faint context, and the two "
            "hypothesis densities of panel e are marked at their directions - "
            "optic flow at 0 degrees, vectorial encoding at both +/-90 - against "
            "the uniform reference, drawn as a dashed circle.  The markers are "
            "point estimates only: their +/-HDIs would be 3.6 pt long at the "
            "polar's scale (shorter than a marker), so panel e remains the place "
            "the intervals are read.",
            "panel c is a pair of polar plots, one per model, each against the "
            "observed bootstrap KDE - the reference page had nothing there.  The "
            "observed KDE is drawn solid rather than the dashed line of the "
            "reference and the house panels: on this page dashed is the uniform "
            "reference (panel d's circle, panels e and f's vertical lines), and "
            "a dashed data line reads as a second uniform.  Splitting the models "
            "is the point: on a shared axis the two bands lie so close together "
            "that the only reading is which is nearer the data.",
            "panel c's polars are 68 pt radius against panel d's 118.7, so the "
            "two can share the page's left half with their horizontal tick "
            "labels; their labels sit at the reference page's own clearances "
            "(3.8-5.5 pt beyond the circle horizontally, 9.5-14 pt vertically) "
            "rather than its offsets scaled down, which would run the left "
            "\"±180°\" label off the page.",
            "panel c's polars are 85 pt radius against panel d's 118.7, so the "
            "two can share the page's left half with their horizontal tick "
            "labels; their labels sit at the reference page's own clearances "
            "(3.8-5.5 pt beyond the circle horizontally, 9.5-14 pt vertically) "
            "rather than its offsets scaled down, which would run the left "
            "\"±180°\" label off the page.",
            "panel f's axes height is read from the house swarm's own y limits "
            "after it draws (4 * the stack reach + 4 * the packing clearance = "
            "202.5 pt at dot 4.5 pt and gap 1.25) and the axes re-placed to it: "
            "the swarm sets its limits in points, so the axes height has to "
            "equal that number or the packing the function computes is not the "
            "packing that gets drawn - and the height is not predictable from "
            "the dot size alone, because the greedy packer's stack reach is not "
            "monotonic in the packing distance.  Dots are packed 1.25 times "
            "their own diameter apart: packed at their own size they only touch, "
            "and the vectorial-encoding row, where all 52 animals fall within "
            "1.2 ELPD, then reads as one solid block.",
        ],
        "known_differences": [
            "The page is 900 x 900 pt against the reference's 900 x 600: panels "
            "c and f are new in a third row of 300 pt, and the letters follow "
            "this figure's flow of information (setup, model selection, "
            "predictive check, then hypothesis testing at three levels), so the "
            "reference page's c and d are d and e here.  Their measured geometry "
            "is kept exactly, as are panels a and b; all six sit within the "
            "reference's 900 x 600 plus the new row.",
            "Text is set in Arial; the reference page uses Helvetica Neue, whose "
            "glyph widths differ by 1-4% at the same size and position.",
            "Panel b's cached unimodal fit gives ELPD -3533.40 +/- 38.34 where "
            "the reference page prints -3502.0 +/- 38.8; the bimodal fit matches "
            "(-3404.51 +/- 39.27 against -3403.8 +/- 39.2).",
        ],
        "verification": {
            "reference_region_against_the_previous_build": {
                "method": "text spans and drawing primitives within the reference "
                          "page's own 900 x 600, compared between the previously "
                          "committed PDF and this one; the panel letters and "
                          "titles are excluded, they are redrawn to this figure's "
                          "flow",
                "text_spans": "26 of 26 identical to 0.000 pt",
                "drawing_primitives": "88 of 89 identical to 0.05 pt; the one that "
                                      "moves is the page background rectangle "
                                      "(945 -> 900 pt)",
                "png_note": "over the same region the committed PNG does differ: "
                            "matplotlib places every axes from the bottom of the "
                            "page, so changing the page height re-lays the raster "
                            "phase of the whole figure.  900 pt is exactly 3750 px "
                            "at 300 dpi and 945 pt is 3937.5, and no taller page "
                            "reproduces the old grid, so this one is not "
                            "recoverable.",
            },
            "mean_abs_pixel_difference_against_reference_pdf": {
                "panel_a": 5.2, "panel_b": 12.7, "panel_d": 12.9, "panel_e": 12.3,
                "regions": "panel a: x 28-448, y 5-350; panel b: 28-448, 360-520; "
                           "panel d: 455-880, 5-350; panel e: 455-880, 360-520 "
                           "(pt, top-left origin); the regions include the redrawn "
                           "bold titles, whose wording differs from the reference "
                           "page's in d and e",
                "method": "300 dpi, with the best shift within +/-4 px taken up.  "
                          "Panels c and f are new and have no reference panel.",
            },
        },
    }


def write_html(stem):
    return f"""<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Locust circular-statistics figure</title>
<style>body{{margin:0;background:#f5f5f4;color:#262626;font:15px/1.65 Arial,sans-serif}}
main{{max-width:1100px;margin:40px auto;padding:0 24px}}figure{{margin:0}}
.figure{{background:white;padding:20px;border:1px solid #e6e6e6}}img{{display:block;width:100%;height:auto}}
figcaption{{max-width:94ch;margin:24px auto}}a{{color:#665017}}nav{{margin-bottom:18px}}</style>
<main><nav><a href="{stem}.svg">SVG</a> · <a href="{stem}.pdf">PDF</a> · <a href="{stem}.json">provenance</a></nav>
<figure><div class="figure"><img src="{stem}.png" alt="Locust circular-statistics figure"></div>
<figcaption>{CAPTION}</figcaption></figure></main>
"""


def write_outputs(fig, plotted, output_dir, snapshot_dir=None):
    """Save PDF/PNG/SVG to ``output_dir``; when ``snapshot_dir`` is given, also
    write the committed review set (formats + provenance JSON) there."""
    outputs = {}
    targets = [pathlib.Path(output_dir)]
    if snapshot_dir is not None:
        targets.append(pathlib.Path(snapshot_dir))
    for target in targets:
        target.mkdir(parents=True, exist_ok=True)
        stem = target / STEM
        fig.savefig(stem.with_suffix(".pdf"),
                    metadata={"Creator": "unwrap locust_circ_stats_figure.py",
                              "CreationDate": None})
        fig.savefig(stem.with_suffix(".png"), dpi=300)
        with mpl.rc_context({"svg.fonttype": "path",
                             "svg.hashsalt": "unwrap-locust-circ-stats-v1"}):
            fig.savefig(stem.with_suffix(".svg"), format="svg",
                        metadata={"Date": None,
                                  "Creator": "unwrap locust_circ_stats_figure.py"})
        outputs[target] = {suffix: _sha256(stem.with_suffix(suffix))
                           for suffix in (".pdf", ".png", ".svg")}
    if snapshot_dir is not None:
        target = pathlib.Path(snapshot_dir)
        info = provenance(plotted)
        info["outputs_sha256"] = outputs[target]
        (target / f"{STEM}.json").write_text(json.dumps(info, indent=2) + "\n")
        (target / f"{STEM}.html").write_text(write_html(STEM), encoding="utf-8")
    return outputs


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output-dir", type=pathlib.Path, default=plots.FIG_DIR,
                        help="where to write the generated figure (default: figures/locust)")
    parser.add_argument("--snapshot-dir", type=pathlib.Path, default=SNAPSHOT_DIR,
                        help="review snapshot directory, or '' to skip it "
                             "(default: assets/locust-circ-stats)")
    args = parser.parse_args(argv)

    fig, plotted = build_figure()
    outputs = write_outputs(fig, plotted, args.output_dir,
                            args.snapshot_dir or None)
    for target, hashes in outputs.items():
        print(f"saved {target / (STEM + '.pdf')} "
              f"(png {hashes['.png'][:8]}, pdf {hashes['.pdf'][:8]})")
    return fig


if __name__ == "__main__":
    main()
