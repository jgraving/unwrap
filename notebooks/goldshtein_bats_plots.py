"""Figure helpers for the echolocating-bat homing notebook.

The notebook keeps the statistics (the model fits, the PSIS-LOO
comparisons, the population/tag-level summaries, the posterior
predictive checks and the per-step goal-component responsibilities, in
``goldshtein_bats_posthoc.py``); each ``plot_*`` function takes the fit
results and the modeling table, computes any remaining curve
quantities through ``goldshtein_bats_posthoc`` where needed, and is a
single call from the notebook, followed by ``plt.show()``.

Style and saving follow the shared house convention from
``unwrap.plots``: figures render at display resolution (100 dpi) in the
notebook cell and are saved as 300 dpi vector files under
``<repo>/figures/bats/``. Every public figure function accepts a
``save`` flag (on by default) and an explicit repo-relative
``save_path`` (the notebook always passes it); the written files are
gitignored.
"""
from pathlib import Path

import math
import os

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import Point, LineString, Polygon, MultiPolygon
from shapely.prepared import prep
from matplotlib import colors as mcolors
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.transforms import offset_copy
from matplotlib.patches import (
    ArrowStyle,
    Circle,
    FancyArrowPatch,
    Patch,
    Rectangle,
)
from matplotlib.ticker import FuncFormatter
import contextily as cx
from pyproj import Transformer
import seaborn as sns
import arviz as az

from unwrap import plots as _plots
from unwrap.datasets import bats_lake_outline_path

import goldshtein_bats_posthoc as posthoc
from goldshtein_bats_posthoc import (
    KAPPA_UNIFORM,
    circular_modulo,
    sigmoid,
    softplus,
    _stack_samples_first,
    hdi_1d,
    hdi_along_samples,
    _format_treatment_name,
    log_i0_stable,
    log_vonmises_pdf,
    build_res_treat_from_vonmises_idata,
    summarize_band,
    format_polar,
    format_treatment_title,
)

#: Publication figures for this study live in
#: <repo>/figures/bats/, anchored on this file's location so the output
#: does not depend on the kernel cwd (the notebook may run from notebooks/
#: or from the repository root).
FIG_DIR = Path(__file__).resolve().parent.parent / "figures" / "bats"

_EXTRA_RC = {}


def set_style():
    """One figure style for every figure in this notebook (the shared
    house convention; figures are saved under ``FIG_DIR`` while
    rendering 100 dpi inline)."""
    _plots.set_style(FIG_DIR, _EXTRA_RC)


def _save(fig, save, save_path, transparent=False, dpi=300):
    """Write ``fig`` to ``save_path`` (a repo-relative path such as
    ``figures/bats/x.pdf``) through the house two-res saver. The file
    format follows the path's suffix."""
    if not save:
        return
    p = Path(save_path)
    if p.suffix.lower() == ".svg":
        _plots.save_figure(fig, FIG_DIR, p.stem, pdf=False, png=False, svg=True,
                           transparent=transparent, dpi=dpi)
    else:
        _plots.save_figure(fig, FIG_DIR, p.stem, pdf=True, png=False, svg=False,
                           transparent=transparent, dpi=dpi)


# ------------------------------------------------------------------ anchors
#: map anchors and the Web Mercator EPSG for every basemap figure
LAKE_GEOJSON = bats_lake_outline_path()
WEBMERC_EPSG = 3857
ROOST_LAT, ROOST_LON = 33.112373, 35.584135
TRANS_POINTS = [(33.135136, 35.592712), (33.081932, 35.604462)]

# radial scaling of the polar figure families (the annotated notebook
# used slightly different constants per cell; the values below are the
# cell-local ones of the figure they belong to, set in that function)
DPI = 300
LINE_LW = 0.8


N_COLS = 8  # grid columns of the individual PPC polar grid

#: radial constants of the 2x2 mixture-polar family (plot_mixture_polar_pdfs and its _polar_format_mix helper)
EDGE_RADIUS = 1.0
AMP = 1.0
OUTER_LW = 0.8

def hide_axes(ax):
    """Remove axis labels, ticks, tick labels, and spines."""
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

def _nice_scale_length_m(width_m: float) -> float:
    # choose a nice 1/2/5 * 10^k length around ~20% of map width
    target = max(width_m * 0.2, 1.0)
    exp = np.floor(np.log10(target))
    base = target / (10 ** exp)
    if base < 1.5:
        nice = 1.0
    elif base < 3.5:
        nice = 2.0
    elif base < 7.5:
        nice = 5.0
    else:
        nice = 10.0
    return float(nice * (10 ** exp))

def _pick_order_cols(df):
    """Pick the first reasonable ordering column that exists."""
    for c in ["frame", "step", "t", "time", "timestamp"]:
        if c in df.columns:
            return c
    return None

def _plot_gdf_capture_artists(gdf, ax, **kwargs):
    """
    Plot a GeoDataFrame and return the matplotlib artists added by this call.
    Useful for selectively rasterizing only the heavy trajectory layers.
    """
    n_col = len(ax.collections)
    n_line = len(ax.lines)
    n_patch = len(ax.patches)

    gdf.plot(ax=ax, **kwargs)

    new_artists = []
    new_artists.extend(list(ax.collections[n_col:]))
    new_artists.extend(list(ax.lines[n_line:]))
    new_artists.extend(list(ax.patches[n_patch:]))
    return new_artists

def _set_rasterized(artists, rasterized=True):
    for a in artists:
        try:
            a.set_rasterized(rasterized)
        except Exception:
            pass

def _style_trajectory_artists(artists):
    """
    Make rasterized trajectories look smoother in vector outputs (PDF/SVG).
    """
    for a in artists:
        for method_name, val in [
            ("set_antialiased", True),
            ("set_antialiaseds", [True]),
            ("set_capstyle", "round"),
            ("set_joinstyle", "round"),
            ("set_snap", False),
        ]:
            if hasattr(a, method_name):
                try:
                    getattr(a, method_name)(val)
                except Exception:
                    pass

def add_scale_bar_dual(
    ax,
    total_length_m=20.0,
    n_segments=2,
    segment_colors=("#f2f2f2", "#222222"),
    segment_height_pt=8,
    edgecolor="#f2f2f2",
    edgewidth=0.9,
    label=True,
    label_above=False,
    label_kwargs=None,
    ax_xy=(0.95, 0.08),
    offset_pt=(0, -41),
    align="right",
    direction="left",
    pad_text_pt=4,
    return_info=False
):
    """
    Two-tone horizontal scale bar (meters, data units).
    Placement uses axes coords + point offset. Drawn as rectangles in data coords.
    """
    if label_kwargs is None:
        label_kwargs = dict(
            color="#f2f2f2",
            fontsize=9.5,
            ha="center",
            va="top",
            fontname="Arial",
            fontstretch="normal",
        )

    fig = ax.figure
    if fig.canvas is not None:
        fig.canvas.draw()

    base = offset_copy(ax.transAxes, fig=fig, x=offset_pt[0], y=offset_pt[1], units="points")
    x_disp, y_disp = base.transform(ax_xy)
    x_anchor, y_anchor = ax.transData.inverted().transform((x_disp, y_disp))

    seg_len = total_length_m / float(max(n_segments, 1))
    if direction == "left":
        x0 = (x_anchor if align == "left" else x_anchor - total_length_m)
        xdir = +1
    elif direction == "right":
        x0 = (x_anchor if align == "left" else x_anchor + total_length_m)
        if align == "right":
            x0 = x_anchor - total_length_m
        xdir = +1
    else:
        raise ValueError("direction must be 'left' or 'right'")

    # Convert bar height (points) to data units using axes pixel height
    bbox = ax.get_window_extent()
    axpix_h = bbox.height
    y_min, y_max = ax.get_ylim()
    pt_to_data_y = (y_max - y_min) / axpix_h
    half_h = 0.5 * segment_height_pt * pt_to_data_y
    yb, yt = y_anchor - half_h, y_anchor + half_h

    colors = (segment_colors * ((n_segments + len(segment_colors) - 1) // len(segment_colors)))[:n_segments]

    patches = []
    for i in range(n_segments):
        xa = x0 + xdir * (i * seg_len)
        xb = x0 + xdir * ((i + 1) * seg_len)
        p = plt.Rectangle(
            (min(xa, xb), yb),
            abs(xb - xa),
            yt - yb,
            facecolor=colors[i],
            edgecolor=edgecolor,
            linewidth=edgewidth,
            zorder=5
        )
        ax.add_patch(p)
        patches.append(p)

    x_center = x0 + xdir * (total_length_m / 2.0)
    txt = None
    if label:
        txt = ax.text(x_center, y_anchor, f"{int(total_length_m)} m", **label_kwargs)
        dy = pad_text_pt if label_above else -pad_text_pt
        txt.set_transform(offset_copy(ax.transData, fig=fig, x=0, y=dy, units="points"))
        txt.set_zorder(6)

    if return_info:
        return {
            "center_data": (x_center, y_anchor),
            "y_anchor_data": y_anchor,
            "segment_height_pt": segment_height_pt,
            "offset_pt": offset_pt,
            "ax_xy": ax_xy,
            "total_length_m": total_length_m,
            "patches": patches,
            "text": txt,
        }
    return None

def add_compass(
    ax,
    ax_xy=(0.92, 0.18),
    offset_pt=(0, -95),
    arrow_len_pt=36,
    head_width_pt=10,
    head_length_pt=12,
    tail_width_pt=2.1,
    color="#f2f2f2",
    edgecolor="#f2f2f2",
    edgewidth=1.0,
    alpha=1.0,
    draw_N=True,
    N_fontsize=10,
    N_offset_pt=(0, 3)
):
    """
    North arrow drawn in AXES coordinates so it stays upright and consistent size.
    """
    fig = ax.figure
    if fig.canvas is not None:
        fig.canvas.draw()

    base = offset_copy(ax.transAxes, fig=fig, x=offset_pt[0], y=offset_pt[1], units="points")

    bbox = ax.get_window_extent()
    axpix_h = bbox.height
    dy_axes = (arrow_len_pt * fig.dpi / 72.0) / axpix_h

    x_axes, y_axes = ax_xy
    start = (x_axes, y_axes)
    end = (x_axes, y_axes + dy_axes)

    style = ArrowStyle.Simple(head_length=head_length_pt, head_width=head_width_pt, tail_width=tail_width_pt)
    patch = FancyArrowPatch(
        start, end,
        transform=ax.transAxes,
        arrowstyle=style,
        mutation_scale=1.0,
        linewidth=edgewidth,
        edgecolor=edgecolor,
        facecolor=color,
        alpha=alpha,
        joinstyle="miter",
        capstyle="butt",
        zorder=6
    )
    patch.set_transform(base)
    ax.add_patch(patch)

    if draw_N:
        tx = ax.text(
            x_axes, y_axes + dy_axes, "N",
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=N_fontsize,
            color=color,
            zorder=6,
            fontname="Arial"
        )
        tx.set_transform(offset_copy(
            ax.transAxes, fig=fig,
            x=N_offset_pt[0] + offset_pt[0],
            y=N_offset_pt[1] + offset_pt[1],
            units="points"
        ))

def _add_scale_and_compass(ax):
    """
    Bottom-right scalebar + compass centered above it.
    Styled for a LIGHT basemap (e.g., Positron): dark annotation colors.
    """
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    corner_pad_px = 16  # tighter for compact layout
    pad_pt = corner_pad_px * 72.0 / fig.dpi
    seg_h_pt = 7
    scale_text_breathe_pt = 7

    # scale length from axis width (meters; EPSG:3857)
    x0, x1 = ax.get_xlim()
    scale_len_m = _nice_scale_length_m(abs(x1 - x0))

    # DARK (black) text for Positron/light basemap
    label_kwargs = dict(
        color="#000000",
        fontsize=9.0,
        ha="center",
        va="top",
        fontname="Arial",
        fontstretch="normal"
    )

    def _draw_scale(offset_y_pt):
        return add_scale_bar_dual(
            ax,
            total_length_m=scale_len_m,
            n_segments=2,
            # black/white scalebar for light background
            segment_colors=("#ffffff", "#000000"),
            segment_height_pt=seg_h_pt,
            edgecolor="#000000",
            edgewidth=0.9,
            ax_xy=(1.0, 0.0),
            offset_pt=(-pad_pt, offset_y_pt),
            align="right",
            direction="left",
            label=True,
            label_above=False,
            pad_text_pt=scale_text_breathe_pt,
            label_kwargs=label_kwargs,
            return_info=True
        )

    # initial guess then adjust to exact pixel padding
    offset_y_pt = pad_pt + (seg_h_pt / 2.0) + scale_text_breathe_pt + 4.0
    sb = _draw_scale(offset_y_pt)
    fig.canvas.draw()

    txt = sb["text"]
    if txt is not None:
        ax_bbox = ax.get_window_extent(renderer=renderer)
        txt_bbox = txt.get_window_extent(renderer=renderer)
        current_bottom_pad_px = txt_bbox.y0 - ax_bbox.y0
        dy_px = corner_pad_px - current_bottom_pad_px
        dy_pt = dy_px * 72.0 / fig.dpi

        if abs(dy_pt) > 0.25:
            for p in sb["patches"]:
                p.remove()
            txt.remove()
            sb = _draw_scale(offset_y_pt + dy_pt)
            fig.canvas.draw()

    # compass centered above scale bar
    (cx_data, y_anchor_data) = sb["center_data"]
    cx_disp, cy_disp = ax.transData.transform((cx_data, y_anchor_data))
    compass_gap_pt = (seg_h_pt / 2.0) + 8.0
    cy_disp_up = cy_disp + compass_gap_pt * fig.dpi / 72.0
    cx_axes, cy_axes = ax.transAxes.inverted().transform((cx_disp, cy_disp_up))

    add_compass(
        ax,
        ax_xy=(cx_axes, cy_axes),
        offset_pt=(0, 0),
        arrow_len_pt=30,
        head_width_pt=9,
        head_length_pt=10,
        tail_width_pt=1.9,
        # DARK compass for light basemap
        color="#000000",
        edgecolor="#000000",
        edgewidth=0.9,
        draw_N=True,
        N_fontsize=9,
        N_offset_pt=(0, 2),
    )

def _ensure_dest_cols(track_df, meta_df, tag_col="tag"):
    """Ensure flyHome/flyAgamon exist on track_df without creating _x/_y duplicates."""
    out = track_df.copy()

    need = ["flyHome", "flyAgamon"]
    missing = [c for c in need if c not in out.columns]
    if missing:
        meta_dest = (
            meta_df[[tag_col, "flyHome", "flyAgamon"]]
            .drop_duplicates(subset=[tag_col], keep="first")
        )
        out = out.merge(meta_dest, on=tag_col, how="left", validate="many_to_one")

    out["flyHome"] = pd.to_numeric(out["flyHome"], errors="coerce").fillna(0).astype(int)
    out["flyAgamon"] = pd.to_numeric(out["flyAgamon"], errors="coerce").fillna(0).astype(int)
    return out

def _dest_label_from_flags(flyHome, flyAgamon):
    fh = int(flyHome == 1)
    fa = int(flyAgamon == 1)
    if fh and not fa:
        return "Home"
    if fa and not fh:
        return "Agamon"
    if fh and fa:
        return "Ambiguous"
    return "Unknown"

def _make_segment_and_gap_lines_from_df(
    df,
    *,
    tag_col="tag", lat_col="lat", lon_col="lon",
    treat_col="treatment", treat_name_col="treatmentName",
    segment_col="segment_id",
    gap_col="gap_before",
):
    """
    Build:
      (1) One LineString per contiguous segment (solid)
      (2) One LineString per detected gap (dashed connector)
    """
    d = df.dropna(subset=[tag_col, lat_col, lon_col]).copy()
    d[lat_col] = pd.to_numeric(d[lat_col], errors="coerce")
    d[lon_col] = pd.to_numeric(d[lon_col], errors="coerce")
    d = d[np.isfinite(d[lat_col]) & np.isfinite(d[lon_col])].copy()

    order_col = _pick_order_cols(d)

    seg_records = []
    gap_records = []

    has_segments = segment_col in d.columns
    has_gaps = gap_col in d.columns

    for tag, g_tag in d.groupby(tag_col, sort=False):
        g_tag = g_tag.sort_values(order_col) if order_col else g_tag.sort_index()

        flyHome = int(g_tag["flyHome"].iloc[0]) if "flyHome" in g_tag.columns else 0
        flyAgamon = int(g_tag["flyAgamon"].iloc[0]) if "flyAgamon" in g_tag.columns else 0
        dest = _dest_label_from_flags(flyHome, flyAgamon)

        trt = g_tag[treat_col].iloc[0] if treat_col in g_tag.columns else np.nan
        trt_name = g_tag[treat_name_col].iloc[0] if treat_name_col in g_tag.columns else "Unknown"

        # solid segment lines
        if has_segments:
            for seg_id, g_seg in g_tag.groupby(segment_col, sort=True):
                coords = list(zip(g_seg[lon_col].to_numpy(), g_seg[lat_col].to_numpy()))
                if len(coords) < 2:
                    continue
                seg_records.append({
                    "tag": tag,
                    "treatment": trt,
                    "treatmentName": trt_name,
                    "dest": dest,
                    "segment_id": int(seg_id),
                    "geometry": LineString(coords),
                })
        else:
            coords = list(zip(g_tag[lon_col].to_numpy(), g_tag[lat_col].to_numpy()))
            if len(coords) >= 2:
                seg_records.append({
                    "tag": tag,
                    "treatment": trt,
                    "treatmentName": trt_name,
                    "dest": dest,
                    "segment_id": 0,
                    "geometry": LineString(coords),
                })

        # dashed gap connectors
        if has_gaps:
            gb = pd.to_numeric(g_tag[gap_col], errors="coerce").fillna(0).astype(int).to_numpy()
            gap_starts = np.where(gb == 1)[0]  # this row is first point after a gap
            for j in gap_starts:
                if j <= 0:
                    continue
                x0, y0 = float(g_tag.iloc[j - 1][lon_col]), float(g_tag.iloc[j - 1][lat_col])
                x1, y1 = float(g_tag.iloc[j][lon_col]), float(g_tag.iloc[j][lat_col])
                if not (np.isfinite(x0) and np.isfinite(y0) and np.isfinite(x1) and np.isfinite(y1)):
                    continue
                gap_records.append({
                    "tag": tag,
                    "treatment": trt,
                    "treatmentName": trt_name,
                    "dest": dest,
                    "geometry": LineString([(x0, y0), (x1, y1)]),
                })

    if len(seg_records):
        seg_gdf = gpd.GeoDataFrame(seg_records, geometry="geometry", crs="EPSG:4326")
    else:
        seg_gdf = gpd.GeoDataFrame(
            {"tag": [], "treatment": [], "treatmentName": [], "dest": [], "segment_id": [], "geometry": []},
            geometry="geometry", crs="EPSG:4326"
        )

    if len(gap_records):
        gap_gdf = gpd.GeoDataFrame(gap_records, geometry="geometry", crs="EPSG:4326")
    else:
        gap_gdf = gpd.GeoDataFrame(
            {"tag": [], "treatment": [], "treatmentName": [], "dest": [], "geometry": []},
            geometry="geometry", crs="EPSG:4326"
        )

    seg_gdf = seg_gdf.to_crs(epsg=WEBMERC_EPSG)
    gap_gdf = gap_gdf.to_crs(epsg=WEBMERC_EPSG)
    return seg_gdf, gap_gdf

def plot_destination_panels(
    resampled_df,
    meta_df,
    *,
    tag_col="tag",
    lat_col="lat",
    lon_col="lon",
    treat_col="treatment",
    treat_name_col="treatmentName",

    # colorblind-ish cyan replacement for tab:cyan
    lake_facecolor="#56B4E9",
    lake_alpha=0.35,
    # style knobs
    res_lw: float = 1.6,
    res_alpha: float = 0.95,
    # basemap knobs
    basemap_alpha: float = 1.0,
    basemap_zoom: int | None = None,         # higher zoom => more detailed tiles (but more tiles)
    basemap_retina: bool = True,             # use @2x tiles when provider supports {r}
    basemap_interpolation: str = "bilinear", # smoothing of raster display
    # gap style
    gap_lw: float | None = None,
    gap_alpha: float = 0.85,
    gap_linestyle: str = "--",
    # legend
    show_treatment_legend: bool = True,
    save: bool = True,
    save_path=None,
):
    # Apply your preferred export/style settings ONLY inside this plot,
    # so other figures in the notebook keep "normal" defaults.
    style_rc = {
        "font.family": "Arial",
        "font.sans-serif": ["Arial"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",  # keep text as text in SVG
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    }

    with mpl.rc_context(rc=style_rc):
        # 1) Ensure dest flags exist on rediscretized df
        res_df = _ensure_dest_cols(resampled_df, meta_df, tag_col=tag_col)

        # 2) Build GeoDataFrames: solid segments + dashed gaps
        res_segs_gdf, res_gaps_gdf = _make_segment_and_gap_lines_from_df(
            res_df,
            tag_col=tag_col, lat_col=lat_col, lon_col=lon_col,
            treat_col=treat_col, treat_name_col=treat_name_col,
            segment_col="segment_id",
            gap_col="gap_before",
        )

        # (optional) debug prints
        # print("rediscretized segment-lines:", len(res_segs_gdf))
        # print("gap connectors:", len(res_gaps_gdf))

        if len(res_segs_gdf) == 0:
            raise ValueError("No lines to plot. Check that lat/lon are present and not all NaN.")

        if gap_lw is None:
            gap_lw = max(0.9, res_lw * 0.9)

        # 3) Points layer (reference + translocations)
        pts = [{"site": "Reference", "geometry": Point(ROOST_LON, ROOST_LAT)}]
        pts += [{"site": "Translocation", "geometry": Point(lon, lat)} for (lat, lon) in TRANS_POINTS]
        pts_gdf = gpd.GeoDataFrame(pts, crs="EPSG:4326").to_crs(epsg=WEBMERC_EPSG)

        # 4) Lake polygon
        lake_gdf_3857 = gpd.read_file(LAKE_GEOJSON).to_crs(epsg=WEBMERC_EPSG)

        # 5) Treatment colors (seaborn colorblind)
        treatments = sorted(pd.unique(res_segs_gdf[treat_col].dropna())) if treat_col in res_segs_gdf.columns else []
        if len(treatments):
            treat_colors = sns.color_palette("colorblind", n_colors=len(treatments))
            treat_to_color = {t: mcolors.to_hex(treat_colors[i]) for i, t in enumerate(treatments)}
        else:
            treat_to_color = {}

        # labels (use treatmentName if present)
        if treat_name_col in res_segs_gdf.columns and len(treatments):
            name_map = res_segs_gdf.groupby(treat_col)[treat_name_col].first().to_dict()
            treat_to_label = {t: _format_treatment_name(name_map.get(t, "Unknown")) for t in treatments}
        else:
            treat_to_label = {t: str(t) for t in treatments}

        # basemap provider, optionally "retina" (@2x) if supported
        basemap_source = cx.providers.CartoDB.DarkMatter
        if basemap_retina and hasattr(basemap_source, "build_url"):
            try:
                # returns a URL string with {r} filled (provider supports retina via {r})
                basemap_source = basemap_source.build_url(r="@2x")
            except Exception:
                pass

        # 6) Figure panels
        fig, axes = plt.subplots(1, 2, figsize=(7, 8), constrained_layout=False)
        panels = [("Home", axes[0], "Home"), ("Agamon", axes[1], "Agamon Hula Lake")]

        # unified extent across everything
        layers_for_bounds = [g for g in [res_segs_gdf, res_gaps_gdf, pts_gdf, lake_gdf_3857] if len(g)]
        minx = min(g.total_bounds[0] for g in layers_for_bounds)
        miny = min(g.total_bounds[1] for g in layers_for_bounds)
        maxx = max(g.total_bounds[2] for g in layers_for_bounds)
        maxy = max(g.total_bounds[3] for g in layers_for_bounds)

        pad_x = (maxx - minx) * 0.06 if maxx > minx else 1000
        pad_y = (maxy - miny) * 0.06 if maxy > miny else 1000
        xlim = (minx - pad_x, maxx + pad_x)
        ylim = (miny - pad_y, maxy + pad_y)

        for dest, ax, title in panels:
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)

            # basemap (higher zoom and/or retina improves apparent raster detail)
            try:
                cx.add_basemap(
                    ax,
                    source=basemap_source,
                    attribution=False,
                    alpha=basemap_alpha,
                    zoom=basemap_zoom,              # None = auto; set e.g. 14/15/16 for sharper tiles
                    interpolation=basemap_interpolation,
                )
            except Exception as e:
                print("Basemap failed (continuing without tiles):", e)

            # lake polygon ONLY on Agamon panel
            if dest == "Agamon" and len(lake_gdf_3857):
                lake_gdf_3857.plot(
                    ax=ax,
                    facecolor=lake_facecolor,
                    edgecolor=lake_facecolor,
                    alpha=lake_alpha,
                    linewidth=1.2,
                    zorder=12
                )
                lake_gdf_3857.boundary.plot(
                    ax=ax,
                    color=lake_facecolor,
                    linewidth=1.4,
                    alpha=1.0,
                    zorder=13
                )

            # --- Dashed gap connectors (draw FIRST, then solid segments on top) ---
            if len(res_gaps_gdf):
                sub_gap = res_gaps_gdf.loc[res_gaps_gdf["dest"] == dest]
                for t in treatments:
                    sub_gap_t = sub_gap.loc[sub_gap[treat_col] == t]
                    if len(sub_gap_t):
                        sub_gap_t.plot(
                            ax=ax,
                            color=treat_to_color[t],
                            linewidth=gap_lw,
                            alpha=gap_alpha,
                            linestyle=gap_linestyle,
                            zorder=26
                        )

            # --- Solid rediscretized segments (colored) ---
            sub_seg = res_segs_gdf.loc[res_segs_gdf["dest"] == dest]
            for t in treatments:
                sub_seg_t = sub_seg.loc[sub_seg[treat_col] == t]
                if len(sub_seg_t):
                    sub_seg_t.plot(
                        ax=ax,
                        color=treat_to_color[t],
                        linewidth=res_lw,
                        alpha=res_alpha,
                        zorder=28
                    )

            # reference + translocation markers
            ref = pts_gdf.loc[pts_gdf["site"] == "Reference"]
            ax.scatter(ref.geometry.x, ref.geometry.y,
                       marker="*", s=220, facecolors="white", edgecolors="black",
                       linewidths=0.9, zorder=50)

            trans = pts_gdf.loc[pts_gdf["site"] == "Translocation"]
            ax.scatter(trans.geometry.x, trans.geometry.y,
                       marker="X", s=140, facecolors="white", edgecolors="black",
                       linewidths=0.9, zorder=50)

            ax.set_aspect("equal", adjustable="box")
            ax.set_title(title, fontsize=16)
            hide_axes(ax)

        # --- Legend ---
        treat_handles = (
            [Patch(color=treat_to_color[t], label=treat_to_label[t]) for t in treatments]
            if show_treatment_legend and len(treatments) else
            []
        )

        other_handles = [
            Line2D([0], [0], marker="*", color="none", markerfacecolor="white",
                   markeredgecolor="black", markersize=12, label="Roost"),
            Line2D([0], [0], marker="X", color="none", markerfacecolor="white",
                   markeredgecolor="black", markersize=9, label="Translocation"),
        ]
        if len(res_gaps_gdf):
            other_handles.append(Line2D([0, 1], [0, 0], color="black", lw=gap_lw,
                                        linestyle=gap_linestyle, alpha=gap_alpha, label="Signal Dropout"))

        fig.legend(
            handles=treat_handles + other_handles,
            loc="center",
            bbox_to_anchor=(1.02, 0.5),
            frameon=True,
            framealpha=0.9,
        )

        plt.tight_layout(rect=[0, 0, 0.86, 1])

        # scale + compass on RIGHT panel (Agamon)
        _add_scale_and_compass(axes[1])

        _save(fig, save, save_path)
        return fig

def plot_destination_panels_compact(
    resampled_df,
    meta_df,
    *,
    tag_col="tag",
    lat_col="lat",
    lon_col="lon",
    treat_col="treatment",
    treat_name_col="treatmentName",
    lake_facecolor="tab:cyan",
    lake_alpha=0.35,
    # line style
    res_lw=1.35,
    res_alpha=0.95,
    basemap_alpha=1.0,
    # gap style
    gap_lw=None,
    gap_alpha=0.85,
    gap_linestyle="--",
    # compact figure + larger text
    fig_w=5.6,
    fig_h=3.8,
    title_fs=12.5,
    legend_fs=7.2,
    pdf_raster_dpi=300,                 # dpi used for rasterized layers inside PDF
    svg_raster_dpi=600,  # raster dpi of the single saved PDF export
    save: bool = True,
    save_path=None,
):
    # 1) Ensure destination flags exist
    res_df = _ensure_dest_cols(resampled_df, meta_df, tag_col=tag_col)

    # 2) Build segment/gap geodataframes
    res_segs_gdf, res_gaps_gdf = _make_segment_and_gap_lines_from_df(
        res_df,
        tag_col=tag_col, lat_col=lat_col, lon_col=lon_col,
        treat_col=treat_col, treat_name_col=treat_name_col,
        segment_col="segment_id",
        gap_col="gap_before",
    )

    print("rediscretized segment-lines:", len(res_segs_gdf))
    print("gap connectors:", len(res_gaps_gdf))

    if len(res_segs_gdf) == 0:
        raise ValueError("No lines to plot. Check that lat/lon are present and not all NaN.")

    if gap_lw is None:
        gap_lw = max(0.8, res_lw * 0.85)

    # 3) Reference + translocation points
    pts = [{"site": "Reference", "geometry": Point(ROOST_LON, ROOST_LAT)}]
    pts += [{"site": "Translocation", "geometry": Point(lon, lat)} for (lat, lon) in TRANS_POINTS]
    pts_gdf = gpd.GeoDataFrame(pts, crs="EPSG:4326").to_crs(epsg=WEBMERC_EPSG)

    # 4) Lake polygon
    lake_gdf_3857 = gpd.read_file(LAKE_GEOJSON).to_crs(epsg=WEBMERC_EPSG)

    # 5) Treatment colors (seaborn colorblind)
    treatments = sorted(pd.unique(res_segs_gdf[treat_col].dropna())) if treat_col in res_segs_gdf.columns else []
    if len(treatments):
        treat_colors = sns.color_palette("colorblind", n_colors=len(treatments))
        treat_to_color = {t: mcolors.to_hex(treat_colors[i]) for i, t in enumerate(treatments)}
    else:
        treat_to_color = {}

    if treat_name_col in res_segs_gdf.columns and len(treatments):
        name_map = res_segs_gdf.groupby(treat_col)[treat_name_col].first().to_dict()
        treat_to_label = {t: _format_treatment_name(name_map.get(t, "Unknown")) for t in treatments}
    else:
        treat_to_label = {t: str(t) for t in treatments}

    # 6) Figure: compact + larger fonts relative to size
    fig, axes = plt.subplots(1, 2, figsize=(fig_w, fig_h), dpi=DPI, constrained_layout=False)
    panels = [("Home", axes[0], "Home"), ("Agamon", axes[1], "Agamon Hula Lake")]

    # Track only trajectory artists (solid segments + gap connectors)
    trajectory_artists = []

    # Unified extent
    layers_for_bounds = [g for g in [res_segs_gdf, res_gaps_gdf, pts_gdf, lake_gdf_3857] if len(g)]
    minx = min(g.total_bounds[0] for g in layers_for_bounds)
    miny = min(g.total_bounds[1] for g in layers_for_bounds)
    maxx = max(g.total_bounds[2] for g in layers_for_bounds)
    maxy = max(g.total_bounds[3] for g in layers_for_bounds)

    pad_x = (maxx - minx) * 0.05 if maxx > minx else 1000
    pad_y = (maxy - miny) * 0.05 if maxy > miny else 1000
    xlim = (minx - pad_x, maxx + pad_x)
    ylim = (miny - pad_y, maxy + pad_y)

    for dest, ax, title in panels:
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)

        # basemap (already raster by nature)
        try:
            cx.add_basemap(ax, source=cx.providers.CartoDB.PositronNoLabels, attribution=False, alpha=basemap_alpha)
        except Exception as e:
            print("Basemap failed (continuing without tiles):", e)

        # lake ONLY on Agamon panel (keep vector)
        if dest == "Agamon" and len(lake_gdf_3857):
            lake_gdf_3857.plot(
                ax=ax,
                facecolor=lake_facecolor,
                edgecolor=lake_facecolor,
                alpha=lake_alpha,
                linewidth=1.0,
                zorder=12,
            )
            lake_gdf_3857.boundary.plot(
                ax=ax,
                color=lake_facecolor,
                linewidth=1.1,
                alpha=1.0,
                zorder=13,
            )

        # dashed gap connectors (heavy layer -> candidate for rasterization)
        if len(res_gaps_gdf):
            sub_gap = res_gaps_gdf.loc[res_gaps_gdf["dest"] == dest]
            for t in treatments:
                sub_gap_t = sub_gap.loc[sub_gap[treat_col] == t]
                if len(sub_gap_t):
                    artists = _plot_gdf_capture_artists(
                        sub_gap_t,
                        ax=ax,
                        color=treat_to_color[t],
                        linewidth=gap_lw,
                        alpha=gap_alpha,
                        linestyle=gap_linestyle,
                        zorder=26,
                    )
                    trajectory_artists.extend(artists)

        # solid rediscretized segments (heavy layer -> candidate for rasterization)
        sub_seg = res_segs_gdf.loc[res_segs_gdf["dest"] == dest]
        for t in treatments:
            sub_seg_t = sub_seg.loc[sub_seg[treat_col] == t]
            if len(sub_seg_t):
                artists = _plot_gdf_capture_artists(
                    sub_seg_t,
                    ax=ax,
                    color=treat_to_color[t],
                    linewidth=res_lw,
                    alpha=res_alpha,
                    zorder=28,
                )
                trajectory_artists.extend(artists)

        # Improve trajectory appearance whether vector or rasterized
        _style_trajectory_artists(trajectory_artists)

        # reference + translocation markers (keep vector)
        ref = pts_gdf.loc[pts_gdf["site"] == "Reference"]
        ax.scatter(
            ref.geometry.x, ref.geometry.y,
            marker="*", s=150, facecolors="white", edgecolors="black",
            linewidths=0.8, zorder=50
        )

        trans = pts_gdf.loc[pts_gdf["site"] == "Translocation"]
        ax.scatter(
            trans.geometry.x, trans.geometry.y,
            marker="X", s=90, facecolors="white", edgecolors="black",
            linewidths=0.8, zorder=50
        )

        ax.set_aspect("equal", adjustable="box")
        ax.set_title(title, fontsize=title_fs, pad=2.0)
        hide_axes(ax)
        
        # black outline around the map panel
        rect = plt.Rectangle(
            (0, 0), 1, 1,
            transform=ax.transAxes,
            fill=False,
            edgecolor="black",
            linewidth=1.0,
            zorder=100,
            clip_on=False
        )
        ax.add_patch(rect)

    # --- Legend: top-right INSIDE right plot, compact and properly fitting ---
    other_handles = [
        Line2D([0], [0], marker="*", color="none",
               markerfacecolor="white", markeredgecolor="black",
               markeredgewidth=0.8, markersize=7.2, label="Roost"),
        Line2D([0], [0], marker="X", color="none",
               markerfacecolor="white", markeredgecolor="black",
               markeredgewidth=0.8, markersize=6.2, label="Translocation"),
    ]
    if len(res_gaps_gdf):
        other_handles.append(
            Line2D([0, 1], [0, 0], color="black", lw=1.0,
                   linestyle=gap_linestyle, alpha=0.95, label="Signal dropout")
        )

    axes[1].legend(
        handles=other_handles,
        loc="upper right",
        bbox_to_anchor=(0.985, 0.985),
        borderaxespad=0.2,
        frameon=True,
        framealpha=0.92,
        facecolor="white",
        edgecolor="0.8",
        fontsize=legend_fs,
        handlelength=1.4,
        handletextpad=0.45,
        labelspacing=0.35,
        borderpad=0.35,
        columnspacing=0.7,
        markerscale=0.9,
    )

    # Tight/compact layout
    fig.subplots_adjust(left=0.01, right=0.995, bottom=0.02, top=0.90, wspace=0.02)

    # scale + compass on RIGHT panel (now styled dark for Positron)
    _add_scale_and_compass(axes[0])

    _rasterize = str(save_path).endswith(".pdf")
    _set_rasterized(trajectory_artists, _rasterize)
    _save(fig, save, save_path, dpi=pdf_raster_dpi if _rasterize else svg_raster_dpi)
    return fig, axes

def _nice_scale_length_m_flow(width_m: float, target_frac: float = 0.22) -> float:
    width_m = float(abs(width_m))
    if not np.isfinite(width_m) or width_m <= 0:
        return 1000.0
    target = max(width_m * float(target_frac), 1.0)
    exp = np.floor(np.log10(target))
    base = target / (10 ** exp)
    if base < 1.5:
        nice = 1.0
    elif base < 3.5:
        nice = 2.0
    elif base < 7.5:
        nice = 5.0
    else:
        nice = 10.0
    return float(nice * (10 ** exp))

def _format_scale_label(total_m: float) -> str:
    total_m = float(total_m)
    if np.isfinite(total_m) and total_m >= 1000 and abs(total_m % 1000) < 1e-6:
        return f"{int(round(total_m/1000))} km"
    return f"{int(round(total_m))} m"

def _fallback_add_scale_and_compass(ax, *, target_frac=0.22):
    x0, x1 = ax.get_xlim()
    width_m = abs(x1 - x0)
    L = _nice_scale_length_m_flow(width_m, target_frac=target_frac)

    # scale bar (axes coords)
    xA0, yA0 = 0.05, 0.06
    frac = (L / width_m) if width_m > 0 else 0.2
    bar = min(0.40, max(0.14, frac))

    ax.plot([xA0, xA0 + bar], [yA0, yA0], transform=ax.transAxes, color="white",
            lw=6, solid_capstyle="butt", zorder=300)
    ax.plot([xA0, xA0 + bar / 2], [yA0, yA0], transform=ax.transAxes, color="black",
            lw=6, solid_capstyle="butt", zorder=301)
    ax.text(xA0 + bar / 2, yA0 - 0.02, _format_scale_label(L),
            transform=ax.transAxes, ha="center", va="top", color="white", fontsize=9, zorder=302)

    ax.annotate("N", xy=(0.92, 0.20), xytext=(0.92, 0.08),
                xycoords="axes fraction", textcoords="axes fraction",
                ha="center", va="bottom", color="white", fontsize=10,
                arrowprops=dict(arrowstyle="-|>", color="white", lw=1.4),
                zorder=303)

def _safe_normalize(dx, dy):
    nrm = np.hypot(dx, dy)
    outx = np.full_like(dx, np.nan, dtype=float)
    outy = np.full_like(dy, np.nan, dtype=float)
    m = np.isfinite(nrm) & (nrm > 0)
    outx[m] = dx[m] / nrm[m]
    outy[m] = dy[m] / nrm[m]
    return outx, outy

def _make_lines_from_df(df, *, tag_col="tag", lat_col="lat", lon_col="lon"):
    d = df.dropna(subset=[tag_col, lat_col, lon_col]).copy()
    d[lat_col] = pd.to_numeric(d[lat_col], errors="coerce")
    d[lon_col] = pd.to_numeric(d[lon_col], errors="coerce")
    d = d[np.isfinite(d[lat_col]) & np.isfinite(d[lon_col])].copy()

    order_col = _pick_order_cols(d)
    recs = []
    for tag, g_tag in d.groupby(tag_col, sort=False):
        g_tag = g_tag.sort_values(order_col) if order_col else g_tag.sort_index()
        flyHome = int(g_tag["flyHome"].iloc[0])
        flyAgamon = int(g_tag["flyAgamon"].iloc[0])
        dest = _dest_label_from_flags(flyHome, flyAgamon)
        coords = list(zip(g_tag[lon_col].to_numpy(), g_tag[lat_col].to_numpy()))
        if len(coords) < 2:
            continue
        recs.append({"tag": tag, "dest": dest, "geometry": LineString(coords)})

    if recs:
        return gpd.GeoDataFrame(recs, geometry="geometry", crs="EPSG:4326").to_crs(epsg=WEBMERC_EPSG)

    return gpd.GeoDataFrame({"tag": [], "dest": [], "geometry": []}, geometry="geometry", crs="EPSG:4326").to_crs(
        epsg=WEBMERC_EPSG
    )

def lake_outer_polygon_largest(geom):
    """
    Return a Polygon using the EXTERIOR ring of the largest polygon in `geom`.
    This is used ONLY for mask/field calculations (stable outer boundary).
    """
    if geom is None:
        return None

    if isinstance(geom, Polygon):
        main = geom
    elif isinstance(geom, MultiPolygon):
        main = max(list(geom.geoms), key=lambda p: float(p.area))
    else:
        try:
            return lake_outer_polygon_largest(geom.buffer(0))
        except Exception:
            return None

    return Polygon(main.exterior.coords)  # drop interiors/holes

def polygon_mask_on_grid(poly, x_centers, y_centers, chunk=200_000):
    X, Y = np.meshgrid(x_centers, y_centers)
    xx, yy = X.ravel(), Y.ravel()

    try:
        import shapely
        pts = shapely.points(xx, yy)
        inside = shapely.contains(poly, pts)
        return np.asarray(inside, dtype=bool).reshape(Y.shape)
    except Exception:
        pass

    poly_p = prep(poly)
    inside = np.zeros(xx.size, dtype=bool)
    for a in range(0, xx.size, int(chunk)):
        b = min(xx.size, a + int(chunk))
        inside[a:b] = [poly_p.contains(Point(x, y)) for x, y in zip(xx[a:b], yy[a:b])]
    return inside.reshape(Y.shape)

def _adaptive_grid_res(xlim, ylim, base_res_m, max_n=420):
    x0, x1 = xlim
    y0, y1 = ylim
    w = float(abs(x1 - x0))
    h = float(abs(y1 - y0))
    res = float(base_res_m)
    if not np.isfinite(res) or res <= 0:
        res = max(w, h) / max_n if max(w, h) > 0 else 200.0
    res = max(res, w / max_n if w > 0 else res, h / max_n if h > 0 else res)
    return res

def solve_laplace_dirichlet(mask_lake, *, tol=1e-5, max_iter=20000):
    ny, nx = mask_lake.shape
    u = np.ones((ny, nx), float)
    u[mask_lake] = 0.0
    u[0, :] = u[-1, :] = u[:, 0] = u[:, -1] = 1.0

    free = (~mask_lake).copy()
    free[0, :] = free[-1, :] = free[:, 0] = free[:, -1] = False
    free_flat = np.flatnonzero(free.ravel())

    try:
        import scipy.sparse as sp
        import scipy.sparse.linalg as spla

        idx = -np.ones((ny, nx), int)
        idx.ravel()[free_flat] = np.arange(free_flat.size)
        n = free_flat.size
        if n == 0:
            return u

        A = sp.lil_matrix((n, n), float)
        b = np.zeros(n, float)

        def add(eq, r, c):
            j = idx[r, c]
            if j >= 0:
                A[eq, j] -= 1.0
            else:
                b[eq] += u[r, c]

        rr, cc = np.where(free)
        for r, c in zip(rr, cc):
            eq = idx[r, c]
            A[eq, eq] = 4.0
            add(eq, r - 1, c)
            add(eq, r + 1, c)
            add(eq, r, c - 1)
            add(eq, r, c + 1)

        sol, _ = spla.cg(A.tocsr(), b, atol=0.0, rtol=1e-8, maxiter=3000)
        u2 = u.ravel()
        u2[free_flat] = sol
        return u2.reshape((ny, nx))
    except Exception:
        u_old = u.copy()
        for _ in range(int(max_iter)):
            u_new = u_old.copy()
            u_new[1:-1, 1:-1] = 0.25 * (
                u_old[:-2, 1:-1] + u_old[2:, 1:-1] + u_old[1:-1, :-2] + u_old[1:-1, 2:]
            )
            u_new[mask_lake] = 0.0
            u_new[0, :] = u_new[-1, :] = u_new[:, 0] = u_new[:, -1] = 1.0
            if np.nanmax(np.abs(u_new - u_old)) < tol:
                return u_new
            u_old = u_new
        return u_old

def harmonic_field_on_extent(xlim, ylim, lake_poly_calc_3857, *, grid_res_m=160.0, tol=1e-5, max_grid_n=420):
    grid_res_m = _adaptive_grid_res(xlim, ylim, grid_res_m, max_n=max_grid_n)
    x0, x1 = xlim
    y0, y1 = ylim
    x_centers = np.arange(min(x0, x1), max(x0, x1) + grid_res_m, grid_res_m)
    y_centers = np.arange(min(y0, y1), max(y0, y1) + grid_res_m, grid_res_m)

    mask = polygon_mask_on_grid(lake_poly_calc_3857, x_centers, y_centers)
    u = solve_laplace_dirichlet(mask, tol=tol)

    du_dy, du_dx = np.gradient(u, grid_res_m, grid_res_m)
    tx, ty = _safe_normalize(-du_dx, -du_dy)

    X, Y = np.meshgrid(x_centers, y_centers)
    return X, Y, tx, ty, mask, grid_res_m

def smoothed_distance_field_on_extent(
    xlim, ylim, lake_poly_calc_3857, *, grid_res_m=160.0, sigma_pix=3.0, max_grid_n=420
):
    grid_res_m = _adaptive_grid_res(xlim, ylim, grid_res_m, max_n=max_grid_n)
    x0, x1 = xlim
    y0, y1 = ylim
    x_centers = np.arange(min(x0, x1), max(x0, x1) + grid_res_m, grid_res_m)
    y_centers = np.arange(min(y0, y1), max(y0, y1) + grid_res_m, grid_res_m)

    mask = polygon_mask_on_grid(lake_poly_calc_3857, x_centers, y_centers)

    from scipy.ndimage import distance_transform_edt, gaussian_filter

    dist = distance_transform_edt(~mask, sampling=(grid_res_m, grid_res_m)).astype(float)
    dist[mask] = 0.0
    dist_s = gaussian_filter(dist, sigma=float(sigma_pix), mode="nearest")

    du_dy, du_dx = np.gradient(dist_s, grid_res_m, grid_res_m)
    tx, ty = _safe_normalize(-du_dx, -du_dy)

    X, Y = np.meshgrid(x_centers, y_centers)
    return X, Y, tx, ty, mask, grid_res_m

def roost_field_on_extent(xlim, ylim, roost_xy_3857, *, grid_res_m=220.0, max_grid_n=420):
    grid_res_m = _adaptive_grid_res(xlim, ylim, grid_res_m, max_n=max_grid_n)
    x0, x1 = xlim
    y0, y1 = ylim
    x_centers = np.arange(min(x0, x1), max(x0, x1) + grid_res_m, grid_res_m)
    y_centers = np.arange(min(y0, y1), max(y0, y1) + grid_res_m, grid_res_m)

    X, Y = np.meshgrid(x_centers, y_centers)
    rx, ry = float(roost_xy_3857[0]), float(roost_xy_3857[1])
    dx = rx - X
    dy = ry - Y
    tx, ty = _safe_normalize(dx, dy)
    return X, Y, tx, ty, grid_res_m

def plot_destination_flow_fields(
    resampled_df,
    meta_df,
    *,
    tag_col="tag",
    lat_col="lat",
    lon_col="lon",
    figsize=(8, 9),                     # EXACT match to your trajectory panels
    lake_facecolor="tab:cyan",
    lake_alpha=0.35,
    lake_edge_lw=1.4,
    basemap_alpha: float = 1.0,
    basemap_source=cx.providers.CartoDB.DarkMatter,
    show_basemap: bool = True,
    show_ui_scale_compass: bool = True, # NEW: turn off scale+compass
    show_legend: bool = True,

    # ---- choose lake field ----
    lake_field: str = "harmonic",       # {"harmonic","smoothed_distance"}

    # ---- field resolution knobs (separate for lake vs roost) ----
    lake_grid_res_m: float = 160.0,
    lake_max_grid_n: int = 520,
    lake_jacobi_tol: float = 1e-5,
    lake_smooth_sigma_pix: float = 3.0,

    roost_grid_res_m: float = 220.0,
    roost_max_grid_n: int = 520,

    # ---- quiver display knobs (bigger arrows -> lower scale, larger width/head sizes) ----
    quiver_stride_lake: int = 1,
    quiver_stride_roost: int = 1,

    quiver_scale_lake: float = 22.0,    # LOWER = BIGGER arrows
    quiver_scale_roost: float = 22.0,   # LOWER = BIGGER arrows

    quiver_width_lake: float = 0.0042,
    quiver_width_roost: float = 0.0042,

    quiver_alpha_lake: float = 0.95,
    quiver_alpha_roost: float = 0.95,

    quiver_color_lake: str = "white",
    quiver_color_roost: str = "white",

    quiver_headwidth: float = 5.5,
    quiver_headlength: float = 7.0,
    quiver_headaxislength: float = 6.2,
    quiver_minshaft: float = 2.2,
    save: bool = True,
    save_path=None,
):
    # 1) Ensure dest flags exist (extent matching)
    res_df = _ensure_dest_cols(resampled_df, meta_df, tag_col=tag_col)

    # 2) Build GeoDataFrame lines ONLY to replicate extent computation (same as trajectory plot)
    res_lines_gdf = _make_lines_from_df(res_df, tag_col=tag_col, lat_col=lat_col, lon_col=lon_col)
    if len(res_lines_gdf) == 0:
        raise ValueError("No rediscretized lines found for extent. Check lat/lon.")

    # 3) Points layer (reference + translocations)
    pts = [{"site": "Reference", "geometry": Point(ROOST_LON, ROOST_LAT)}]
    pts += [{"site": "Translocation", "geometry": Point(lon, lat)} for (lat, lon) in TRANS_POINTS]
    pts_gdf = gpd.GeoDataFrame(pts, crs="EPSG:4326").to_crs(epsg=WEBMERC_EPSG)

    # 4) Lake polygon (PLOT FULL geometry; COMPUTE field on largest outer polygon)
    lake_gdf_3857 = gpd.read_file(LAKE_GEOJSON).to_crs(epsg=WEBMERC_EPSG)
    if len(lake_gdf_3857) == 0:
        raise ValueError("Lake GeoJSON is empty or could not be read.")

    lake_full_union = (
        lake_gdf_3857.geometry.union_all()
        if hasattr(lake_gdf_3857.geometry, "union_all")
        else lake_gdf_3857.unary_union
    )
    lake_poly_calc_3857 = lake_outer_polygon_largest(lake_full_union)
    if lake_poly_calc_3857 is None:
        raise ValueError("Could not derive a largest outer polygon for lake calculations.")

    # 5) unified extent (SAME logic as trajectory plot)
    layers_for_bounds = [g for g in [res_lines_gdf, pts_gdf, lake_gdf_3857] if len(g)]
    minx = min(g.total_bounds[0] for g in layers_for_bounds)
    miny = min(g.total_bounds[1] for g in layers_for_bounds)
    maxx = max(g.total_bounds[2] for g in layers_for_bounds)
    maxy = max(g.total_bounds[3] for g in layers_for_bounds)

    pad_x = (maxx - minx) * 0.06 if maxx > minx else 1000
    pad_y = (maxy - miny) * 0.06 if maxy > miny else 1000
    xlim = (minx - pad_x, maxx + pad_x)
    ylim = (miny - pad_y, maxy + pad_y)

    # roost in 3857
    roost_3857 = gpd.GeoSeries([Point(ROOST_LON, ROOST_LAT)], crs="EPSG:4326").to_crs(epsg=WEBMERC_EPSG).iloc[0]
    roost_xy = (roost_3857.x, roost_3857.y)

    # 6) compute flow fields on SAME extent
    XR, YR, txR, tyR, _ = roost_field_on_extent(
        xlim, ylim, roost_xy, grid_res_m=float(roost_grid_res_m), max_grid_n=int(roost_max_grid_n)
    )

    field = str(lake_field).strip().lower()
    if field in ("harmonic", "laplace", "potential"):
        XL, YL, txL, tyL, maskL, _ = harmonic_field_on_extent(
            xlim, ylim, lake_poly_calc_3857,
            grid_res_m=float(lake_grid_res_m),
            tol=float(lake_jacobi_tol),
            max_grid_n=int(lake_max_grid_n),
        )
        lake_title = "Agamon Hula Lake"
    elif field in ("smoothed_distance", "smoothed", "edt", "distance"):
        XL, YL, txL, tyL, maskL, _ = smoothed_distance_field_on_extent(
            xlim, ylim, lake_poly_calc_3857,
            grid_res_m=float(lake_grid_res_m),
            sigma_pix=float(lake_smooth_sigma_pix),
            max_grid_n=int(lake_max_grid_n),
        )
        lake_title = f"Agamon Hula Lake"
    else:
        raise ValueError("lake_field must be one of {'harmonic','smoothed_distance'}")

    # 7) plot (same layout/size)
    fig, axes = plt.subplots(1, 2, figsize=figsize, constrained_layout=False)
    panels = [("Home", axes[0], "Home (roost flow)"), ("Agamon", axes[1], lake_title)]

    for dest, ax, title in panels:
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)

        if show_basemap:
            try:
                cx.add_basemap(ax, source=basemap_source, attribution=False, alpha=basemap_alpha)
            except Exception as e:
                print("Basemap failed (continuing without tiles):", e)

        # Plot FULL lake polygon only on Agamon panel (as before)
        if dest == "Agamon":
            lake_gdf_3857.plot(
                ax=ax,
                facecolor=lake_facecolor,
                edgecolor=lake_facecolor,
                alpha=lake_alpha,
                linewidth=1.2,
                zorder=12,
            )
            lake_gdf_3857.boundary.plot(
                ax=ax,
                color=lake_facecolor,
                linewidth=lake_edge_lw,
                alpha=1.0,
                zorder=13,
            )

        # --- flow fields ---
        if dest == "Home":
            s = max(1, int(quiver_stride_roost))
            ax.quiver(
                XR[::s, ::s], YR[::s, ::s],
                txR[::s, ::s], tyR[::s, ::s],
                scale=float(quiver_scale_roost),
                color=quiver_color_roost,
                width=float(quiver_width_roost),
                alpha=float(quiver_alpha_roost),
                headwidth=float(quiver_headwidth),
                headlength=float(quiver_headlength),
                headaxislength=float(quiver_headaxislength),
                minshaft=float(quiver_minshaft),
                zorder=25,
            )
        else:
            s = max(1, int(quiver_stride_lake))
            # IMPORTANT: do NOT draw inside lake (mask from the *largest outer polygon*)
            m = maskL[::s, ::s]
            tx_plot = np.ma.array(txL[::s, ::s], mask=m)
            ty_plot = np.ma.array(tyL[::s, ::s], mask=m)
            ax.quiver(
                XL[::s, ::s], YL[::s, ::s],
                tx_plot, ty_plot,
                scale=float(quiver_scale_lake),
                color=quiver_color_lake,
                width=float(quiver_width_lake),
                alpha=float(quiver_alpha_lake),
                headwidth=float(quiver_headwidth),
                headlength=float(quiver_headlength),
                headaxislength=float(quiver_headaxislength),
                minshaft=float(quiver_minshaft),
                zorder=25,
            )

        # reference + translocation markers
        ref = pts_gdf.loc[pts_gdf["site"] == "Reference"]
        ax.scatter(ref.geometry.x, ref.geometry.y,
                   marker="*", s=220, facecolors="white", edgecolors="black",
                   linewidths=0.9, zorder=50)

        trans = pts_gdf.loc[pts_gdf["site"] == "Translocation"]
        ax.scatter(trans.geometry.x, trans.geometry.y,
                   marker="X", s=140, facecolors="white", edgecolors="black",
                   linewidths=0.9, zorder=50)

        ax.set_aspect("equal", adjustable="box")
        ax.set_title(title)
        hide_axes(ax)

    # Legend (markers only)
    if show_legend:
        other_handles = [
            Line2D([0], [0], marker="*", color="none", markerfacecolor="white",
                   markeredgecolor="black", markersize=12, label="Roost"),
            Line2D([0], [0], marker="X", color="none", markerfacecolor="white",
                   markeredgecolor="black", markersize=9, label="Translocation"),
        ]
        fig.legend(
            handles=other_handles,
            loc="center",
            bbox_to_anchor=(1.0, 0.5),
            frameon=True,
            framealpha=0.9,
        )

    plt.tight_layout(rect=[0, 0, 0.86, 1])

    # Optional scale + compass on RIGHT panel
    if show_ui_scale_compass:
        _add_scale_and_compass(axes[1])

    _save(fig, save, save_path)
    return fig

def plot_per_tag_polar_kde(
    track_w_treat,
    angle_col="flightAngDiffDest_deg",
    tag_col="tag",
    treat_col="treatment",
    treat_name_col="treatmentName",
    kappa=30,
    grid_n=360,
    r0=1.0,              # density=0 -> r=1.0
    r_scale=1.0,         # density=1 -> r=2.0
    lw=2.0,
    ring_color="0.60",
    ring_lw=0.9,
    ring_alpha=0.65,
    spoke_color="0.75",
    spoke_lw=0.8,
    spoke_alpha=0.7,
    zero_seg_color="0.65",
    zero_seg_lw=1.6,
    zero_seg_alpha=0.95,
    title_fs=9,
    save: bool = True,
    save_path=None,
):
    """
    - theta shown as (-180, 180]
    - ONLY two density rings: r=1 (dens=0) and r=2 (dens=1)
    - keep ONLY the 4 main direction spokes + tick labels
    - add a gray segment at 0° from r=1 to r=2
    - remove polar spine + default grid + radial ticks
    """
    df = track_w_treat[[tag_col, treat_col, treat_name_col, angle_col]].copy()
    df = df.dropna(subset=[tag_col, treat_col, angle_col])

    # degrees -> radians, wrapped to (-pi, pi]
    ang = df[angle_col].astype(float).to_numpy()
    ang_wrapped = ((ang + 180.0) % 360.0) - 180.0
    df["_theta"] = np.deg2rad(ang_wrapped)

    tag_meta = (
        df.groupby(tag_col, as_index=False)
          .agg(**{
              treat_col: (treat_col, "first"),
              treat_name_col: (treat_name_col, "first"),
          })
          .sort_values([treat_col, tag_col])
    )
    tags_ordered = tag_meta[tag_col].tolist()

    treatments = sorted(tag_meta[treat_col].unique())
    cmap = mpl.colormaps.get_cmap("tab10")
    treat_to_color = {t: cmap(i % 10) for i, t in enumerate(treatments)}

    n = len(tags_ordered)
    ncols = math.ceil(math.sqrt(n))
    nrows = math.ceil(n / ncols)

    fig, axes = plt.subplots(
        nrows, ncols,
        subplot_kw={"projection": "polar"},
        figsize=(3.2 * ncols, 3.2 * nrows),
    )
    axes = np.array(axes).reshape(-1)

    # KDE grid in (-pi, pi]
    theta_grid = np.linspace(-np.pi, np.pi, grid_n, endpoint=False)
    norm = 2 * np.pi * np.i0(kappa)

    def vonmises_kde(sample_theta: np.ndarray) -> np.ndarray:
        if sample_theta.size == 0:
            return np.zeros_like(theta_grid)
        return np.exp(kappa * np.cos(theta_grid[:, None] - sample_theta[None, :])).mean(axis=1) / norm

    dens_by_tag = {}
    max_d_global = 0.0
    for tag in tags_ordered:
        th = df.loc[df[tag_col] == tag, "_theta"].to_numpy()
        dens = vonmises_kde(th)
        dens_by_tag[tag] = dens
        md = float(np.nanmax(dens)) if dens.size else 0.0
        if np.isfinite(md):
            max_d_global = max(max_d_global, md)

    # only two rings
    ring_radii = [r0, r0 + r_scale * 1.0]
    ring_theta = np.linspace(-np.pi, np.pi, 500)

    top_density = max(1.0, float(max_d_global))
    rmax = r0 + r_scale * top_density + 0.20

    # 4 main directions in (-180, 180]: -90, 0, 90, 180
    xticks_deg = [-90, 0, 90, 180]
    xticks = np.deg2rad(xticks_deg)
    xtick_labels = [r"$-90^\circ$", r"$0^\circ$", r"$90^\circ$", r"$\pm180^\circ$"]

    for i, tag in enumerate(tags_ordered):
        ax = axes[i]
        trt = tag_meta.loc[tag_meta[tag_col] == tag, treat_col].iloc[0]
        color = treat_to_color[trt]

        dens = dens_by_tag[tag]
        r = r0 + r_scale * dens

        # KDE curve
        ax.plot(theta_grid, r, color=color, linewidth=lw, alpha=0.95)

        # density rings at r=1 and r=2 only
        for rr in ring_radii:
            ax.plot(
                ring_theta, np.full_like(ring_theta, rr),
                color=ring_color, linewidth=ring_lw, alpha=ring_alpha, zorder=0
            )

        # remove default decorations
        ax.grid(False)
        ax.spines["polar"].set_visible(False)
        ax.set_rlim(0, rmax)
        ax.set_rticks([])
        ax.set_yticklabels([])

        # show -180..180
        ax.set_thetamin(-180)
        ax.set_thetamax(180)

        # tick labels
        ax.set_xticks(xticks)
        ax.set_xticklabels(xtick_labels)

        # draw the 4 angle spokes (minimal)
        for th in xticks:
            ax.plot([th, th], [r0, rmax], color=spoke_color, linewidth=spoke_lw, alpha=spoke_alpha, zorder=0)

        # gray segment at 0° from r=1 to r=2
        ax.plot([0.0, 0.0], [r0, r0 + r_scale * 1.0],
                color=zero_seg_color, linewidth=zero_seg_lw, alpha=zero_seg_alpha, zorder=1)

        ax.set_title(f"tag {tag}", fontsize=title_fs, va="bottom")

    for j in range(n, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    _save(fig, save, save_path)
    return fig

def plot_treatment_coefficients_mu_kappa(idata_treatment, df_track, *, save=True, save_path=None):
    """Population mu/kappa coefficients of the Bambi von Mises model (intercept + treatment effects) with circular/linear mean + 95% HDI."""
    post = idata_treatment.posterior

    def _has(varname: str) -> bool:
        return varname in post.data_vars

    def _stack_1d(varname: str) -> np.ndarray:
        if not _has(varname):
            raise KeyError(f"Missing '{varname}' in posterior. Available: {list(post.data_vars)[:40]} ...")
        return post[varname].stack(sample=("chain", "draw")).values.astype(float)

    def circular_modulo(x):
        return (x + np.pi) % (2 * np.pi) - np.pi

    def circular_mean(theta):
        z = np.exp(1j * np.asarray(theta, float))
        return float(circular_modulo(np.angle(z.mean())))

    def circular_hdi(theta, hdi_prob=0.95):
        theta = np.asarray(theta, float)
        mu_c = circular_mean(theta)
        centered = circular_modulo(theta - mu_c)
        h = az.hdi(centered, prob=hdi_prob, skipna=True)
        low_c, high_c = float(h[0]), float(h[1])
        low = float(circular_modulo(low_c + mu_c))
        high = float(circular_modulo(high_c + mu_c))
        return low, high, mu_c

    def linear_hdi_mean(x, hdi_prob=0.95):
        x = np.asarray(x, float)
        h = az.hdi(x, prob=hdi_prob, skipna=True)
        return float(np.nanmean(x)), float(h[0]), float(h[1])

    def circ_dist(a, b):
        return float(circular_modulo(b - a))

    def _cat_dim(varname: str):
        da = post[varname]
        dims = [d for d in da.dims if d not in ("chain", "draw")]
        if len(dims) != 1:
            raise ValueError(f"Expected 1 non-(chain/draw) dim for '{varname}', got {da.dims}")
        return dims[0]

    def _sel_level(varname: str, lvl, nonref_levels):
        dim = _cat_dim(varname)
        coords = list(post[varname].coords[dim].values)
        if lvl in coords:
            return {dim: lvl}
        idx = nonref_levels.index(lvl)
        return {dim: coords[idx]}

    def _format_treatment_name(name: str) -> str:
        s = str(name).strip().lower()
        parts = [p.strip() for p in s.split(",") if p.strip()]
        if "control" in parts and len(parts) > 1:
            parts = [p for p in parts if p != "control"]
        if not parts:
            parts = ["control"]
        parts = sorted(parts)
        parts = [p[:1].upper() + p[1:] for p in parts]
        return ", ".join(parts)

    # ---- treatment levels in fit order (for consistent colors/labels) ----
    treat_levels = (
        df_track["treatment"].cat.categories.tolist()
        if hasattr(df_track["treatment"], "cat")
        else sorted(df_track["treatment"].unique())
    )
    ref = treat_levels[0]
    nonref = treat_levels[1:]

    # ---- nice labels ----
    nice = {t: str(t) for t in treat_levels}
    if "treatmentName" in df_track.columns:
        name_map = df_track.groupby("treatment")["treatmentName"].first().to_dict()
        for t in treat_levels:
            nice[t] = _format_treatment_name(name_map.get(t, t))
    else:
        for t in treat_levels:
            nice[t] = _format_treatment_name(t)

    # ---- detect which coefficient blocks exist ----
    mu_has_intercept = _has("Intercept")
    mu_has_treat     = _has("treatment")

    k_has_intercept  = _has("kappa_Intercept")
    k_has_treat      = _has("kappa_treatment")

    have_mu_any = mu_has_intercept or mu_has_treat
    have_k_any  = k_has_intercept or k_has_treat
    if not (have_mu_any or have_k_any):
        raise ValueError(
            "Neither μ nor κ coefficients found in posterior. "
            "Expected at least one of: 'Intercept', 'treatment', 'kappa_Intercept', 'kappa_treatment'."
        )

    # ---- rows: include treatments if either μ or κ has treatment effects ----
    rows_have_treat = (mu_has_treat or k_has_treat) and (len(nonref) > 0)
    if rows_have_treat:
        terms_order = ["Control (Intercept)"] + [nice[t] for t in nonref]
        row_treatment = [ref] + nonref
    else:
        terms_order = ["Control (Intercept)"]
        row_treatment = [ref]

    y_base = np.arange(len(terms_order), dtype=float)

    # ---- colors by treatment row ----
    cmap = mpl.colormaps.get_cmap("tab10")
    treat_to_color = {t: cmap(i % 10) for i, t in enumerate(treat_levels)}

    # ---- allocate summaries (NaN = not present / not plotted) ----
    n = len(terms_order)
    mu_mean = np.full(n, np.nan)
    mu_lo   = np.full(n, np.nan)
    mu_hi   = np.full(n, np.nan)

    k_mean  = np.full(n, np.nan)
    k_lo    = np.full(n, np.nan)
    k_hi    = np.full(n, np.nan)

    # ---- fill intercept row (if present) ----
    if mu_has_intercept:
        mu_draws = circular_modulo(_stack_1d("Intercept"))
        mu_lo[0], mu_hi[0], mu_mean[0] = circular_hdi(mu_draws, hdi_prob=0.95)

    if k_has_intercept:
        k_draws = _stack_1d("kappa_Intercept")
        k_mean[0], k_lo[0], k_hi[0] = linear_hdi_mean(k_draws, hdi_prob=0.95)

    # ---- fill treatment rows (only for the parameters that exist) ----
    if rows_have_treat:
        for i, lvl in enumerate(nonref, start=1):
            if mu_has_treat:
                mu_draws = (
                    post["treatment"]
                    .sel(**_sel_level("treatment", lvl, nonref))
                    .stack(sample=("chain", "draw"))
                    .values.astype(float)
                )
                mu_draws = circular_modulo(mu_draws)
                mu_lo[i], mu_hi[i], mu_mean[i] = circular_hdi(mu_draws, hdi_prob=0.95)

            if k_has_treat:
                k_draws = (
                    post["kappa_treatment"]
                    .sel(**_sel_level("kappa_treatment", lvl, nonref))
                    .stack(sample=("chain", "draw"))
                    .values.astype(float)
                )
                k_mean[i], k_lo[i], k_hi[i] = linear_hdi_mean(k_draws, hdi_prob=0.95)

    # ---- xerr (NaN-safe; we’ll skip plotting where missing) ----
    mu_xerr = np.full((2, n), np.nan)
    for i in range(n):
        if np.isfinite(mu_mean[i]) and np.isfinite(mu_lo[i]) and np.isfinite(mu_hi[i]):
            mu_xerr[:, i] = [
                abs(circ_dist(mu_lo[i], mu_mean[i])),
                abs(circ_dist(mu_mean[i], mu_hi[i])),
            ]

    k_xerr = np.full((2, n), np.nan)
    for i in range(n):
        if np.isfinite(k_mean[i]) and np.isfinite(k_lo[i]) and np.isfinite(k_hi[i]):
            k_xerr[:, i] = [k_mean[i] - k_lo[i], k_hi[i] - k_mean[i]]

    # ---- y offsets: only separate μ/κ if both are present anywhere ----
    if have_mu_any and have_k_any:
        dy = 0.16
        y_mu = y_base - dy
        y_k  = y_base + dy
    else:
        dy = 0.0
        y_mu = y_base
        y_k  = y_base

    # ---- plot ----
    fig = plt.figure(figsize=(1.8, 0.40 * len(terms_order)), dpi=300)
    ax = plt.gca()

    for i in range(n):
        c = treat_to_color[row_treatment[i]]

        # μ
        if np.isfinite(mu_mean[i]) and np.all(np.isfinite(mu_xerr[:, i])):
            ax.errorbar(
                mu_mean[i], y_mu[i], xerr=mu_xerr[:, [i]],
                fmt="D", color=c, ecolor=c,
                capsize=3, elinewidth=1, capthick=1,
                markersize=5,
                label=None,
            )

        # κ
        if np.isfinite(k_mean[i]) and np.all(np.isfinite(k_xerr[:, i])):
            ax.errorbar(
                k_mean[i], y_k[i], xerr=k_xerr[:, [i]],
                fmt="^", color=c, ecolor=c,
                capsize=3, elinewidth=1, capthick=1,
                markersize=5,
                label=None,
            )

    ax.axvline(0.0, linestyle="--", linewidth=1, color="k", zorder=0)
    ax.set_yticks(y_base)
    ax.set_yticklabels(terms_order)
    ax.invert_yaxis()
    ax.set_xlabel("Coefficient value")

    # ---- legend: only include what we actually plotted ----
    leg_c = "0.2"
    handles = []
    if have_mu_any:
        handles.append(Line2D([0], [0], marker="D", linestyle="None", color=leg_c, markersize=5, label="μ"))
    if have_k_any:
        handles.append(Line2D([0], [0], marker="^", linestyle="None", color=leg_c, markersize=5, label="κ"))
    if handles:
        ax.legend(handles=handles, loc="lower right", frameon=True)

    _save(fig, save, save_path)
    return fig

def plot_treatment_pred_density_at0_mu_kappa(idata_treatment, df_track, *, save=True, save_path=None):
    """Predictive density at 0 deg of the Bambi von Mises model for the population intercept and each treatment shift (mean + 95% HDI, uniform reference line)."""
    post = idata_treatment.posterior

    # ----------------------
    # user knobs
    # ----------------------
    theta0_deg = 0.0
    theta0 = np.deg2rad(theta0_deg)  # radians

    # hide ONLY y ticks + labels (keep y-axis spine)
    HIDE_TREATMENT_TICKS_AND_LABELS = True
    show_treatment_labels = not HIDE_TREATMENT_TICKS_AND_LABELS
    show_treatment_ticks  = not HIDE_TREATMENT_TICKS_AND_LABELS

    # figure sizing (auto-narrow if both labels + ticks hidden)
    fig_width = None            # set a float to force width (e.g., 2.1)
    fig_width_with_labels = 3.3
    fig_width_no_labels   = 2.3

    # marker/interval style to match your other plots
    CAP = 0.18
    LW  = 1.2
    MS  = 5
    HDI_PROB = 0.95
    PALETTE_NAME = "colorblind"

    # ----------------------
    # helpers
    # ----------------------
    def _has(varname: str) -> bool:
        return varname in post.data_vars

    def circular_modulo(x):
        return (x + np.pi) % (2 * np.pi) - np.pi

    def softplus(x):
        x = np.asarray(x, float)
        return np.where(x > 30, x, np.log1p(np.exp(x)))

    def vonmises_pdf(theta, mu, kappa):
        # pdf(theta | mu, kappa) = exp(kappa cos(theta-mu)) / (2π I0(kappa))
        return np.exp(kappa * np.cos(theta - mu)) / (2.0 * np.pi * np.i0(kappa))

    def _stack_samples(da):
        return da.stack(sample=("chain", "draw")).values.astype(float)

    def _format_treatment_name(name: str) -> str:
        s = str(name).strip().lower()
        parts = [p.strip() for p in s.split(",") if p.strip()]
        if "control" in parts and len(parts) > 1:
            parts = [p for p in parts if p != "control"]
        if not parts:
            parts = ["control"]
        parts = sorted(parts)
        parts = [p[:1].upper() + p[1:] for p in parts]
        return ", ".join(parts)

    def _cat_dim(varname: str):
        da = post[varname]
        dims = [d for d in da.dims if d not in ("chain", "draw")]
        if len(dims) != 1:
            raise ValueError(f"Expected 1 non-(chain/draw) dim for '{varname}', got {da.dims}")
        return dims[0]

    def _sel_level(varname: str, lvl, nonref_levels):
        dim = _cat_dim(varname)
        coords = list(post[varname].coords[dim].values)
        if lvl in coords:
            return {dim: lvl}
        idx = nonref_levels.index(lvl)
        return {dim: coords[idx]}

    # ----------------------
    # require intercepts for μ and κ (density needs both)
    # ----------------------
    if not _has("Intercept"):
        raise KeyError(f"Missing 'Intercept' (needed for μ). Available: {list(post.data_vars)[:40]} ...")
    if not _has("kappa_Intercept"):
        raise KeyError(f"Missing 'kappa_Intercept' (needed for κ). Available: {list(post.data_vars)[:40]} ...")

    mu_has_treat = _has("treatment")
    k_has_treat  = _has("kappa_treatment")

    # ----------------------
    # treatment levels (fit order)
    # ----------------------
    treat_levels = (
        df_track["treatment"].cat.categories.tolist()
        if hasattr(df_track["treatment"], "cat")
        else sorted(df_track["treatment"].unique())
    )
    ref = treat_levels[0]
    nonref = treat_levels[1:]

    # ----------------------
    # nice labels
    # ----------------------
    nice = {t: str(t) for t in treat_levels}
    if "treatmentName" in df_track.columns:
        # observed=False avoids dropping unused categories if treatment is categorical
        name_map = df_track.groupby("treatment", observed=False)["treatmentName"].first().to_dict()
        for t in treat_levels:
            nice[t] = _format_treatment_name(name_map.get(t, t))
    else:
        for t in treat_levels:
            nice[t] = _format_treatment_name(t)

    # ----------------------
    # rows to plot
    # include treatment rows if either μ or κ has treatment effects
    # ----------------------
    rows_have_treat = (mu_has_treat or k_has_treat) and (len(nonref) > 0)
    if rows_have_treat:
        labels = ["Control"] + [nice[t] for t in nonref]
        row_treatment = [ref] + nonref
    else:
        labels = ["Control"]
        row_treatment = [ref]

    # ----------------------
    # pull population-level intercept draws
    # ----------------------
    b0_mu = _stack_samples(post["Intercept"])        # (S,)
    b0_k  = _stack_samples(post["kappa_Intercept"])  # (S,)

    # ----------------------
    # compute density at theta0 for each row
    # ----------------------
    dens_draws_by_row = []
    for t in row_treatment:
        eta_mu = b0_mu.copy()
        eta_k  = b0_k.copy()

        if t != ref:
            if mu_has_treat:
                b_trt_mu = _stack_samples(post["treatment"].sel(**_sel_level("treatment", t, nonref)))
                eta_mu = eta_mu + b_trt_mu
            if k_has_treat:
                b_trt_k = _stack_samples(post["kappa_treatment"].sel(**_sel_level("kappa_treatment", t, nonref)))
                eta_k = eta_k + b_trt_k

        mu    = circular_modulo(eta_mu)
        kappa = softplus(eta_k)

        dens0 = vonmises_pdf(theta0, mu, kappa)  # (S,)
        dens_draws_by_row.append(dens0)

    # ----------------------
    # summaries: mean + HDI
    # ----------------------
    mean = np.array([d.mean() for d in dens_draws_by_row], float)
    hdi  = np.array([az.hdi(d, prob=HDI_PROB) for d in dens_draws_by_row], float)  # (K, 2)
    lo, hi = hdi[:, 0], hdi[:, 1]

    # ----------------------
    # seaborn colorblind palette (stable order by treat_levels)
    # ----------------------
    palette = sns.color_palette(PALETTE_NAME, n_colors=len(treat_levels))
    treat_to_color = {t: palette[i] for i, t in enumerate(treat_levels)}
    colors = [treat_to_color[t] for t in row_treatment]

    # ----------------------
    # figure width (auto-narrow if hiding both y labels + y ticks)
    # ----------------------
    if fig_width is None:
        if (not show_treatment_labels) and (not show_treatment_ticks):
            fig_width = fig_width_no_labels
        else:
            fig_width = fig_width_with_labels

    fig_height = 0.35 * len(labels)

    # ----------------------
    # plot (match horizontal mean ± HDI style)
    # ----------------------
    y = np.arange(len(labels), dtype=float)

    fig = plt.figure(figsize=(fig_width, fig_height), dpi=DPI)
    ax = plt.gca()

    for i in range(len(labels)):
        c = colors[i]
        ax.hlines(y[i], lo[i], hi[i], color=c, linewidth=LW, zorder=2)
        ax.vlines([lo[i], hi[i]], y[i] - CAP, y[i] + CAP, color=c, linewidth=LW, zorder=2)
        ax.plot(mean[i], y[i], "o", color=c, markersize=MS, markeredgewidth=0.0, zorder=3)

    # uniform reference (constant density 1/(2π))
    unif = 1.0 / (2.0 * np.pi)
    ax.axvline(unif, color="black", linestyle="--", linewidth=1)

    # y-axis: keep spine, optionally hide ONLY ticks/labels
    ax.set_yticks(y)
    if show_treatment_labels:
        ax.set_yticklabels(labels, fontsize=8)
    else:
        ax.set_yticklabels([""] * len(y))

    ax.tick_params(
        axis="y",
        which="both",
        left=show_treatment_ticks,
        length=3 if show_treatment_ticks else 0,
        labelleft=show_treatment_labels,
    )

    ax.invert_yaxis()
    ax.margins(y=0.12)
    ax.set_xlabel("Predictive density")

    # x-lims (include uniform line)
    xmin = float(min(np.min(lo), unif))
    xmax = float(max(np.max(hi), unif))
    pad = 0.10 * (xmax - xmin if xmax > xmin else 1.0)
    ax.set_xlim(xmin - pad, xmax + pad)

    # style clean-up to match other plot
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)

    plt.tight_layout()
    _save(fig, save, save_path)
    return fig

def plot_total_elpd(
    summary_df,      # output from compute_psisloo_elpd_summary()
    *,
    save=True,
    save_path=None,
    title=None,
    xlabel="PSIS-LOO ELPD",
    dpi=300,
    vline_at_best=True,
):
    out = summary_df.copy()

    # ensure a stable display order (assumes your summary already best-first)
    out = out.reset_index(drop=True)

    y = np.arange(len(out), dtype=float)

    fig = plt.figure(figsize=(3.0, 1.25), dpi=dpi)
    ax = plt.gca()

    if vline_at_best and "is_best" in out.columns and out["is_best"].any():
        best_elpd = float(out.loc[out["is_best"], "elpd"].iloc[0])
        ax.axvline(best_elpd, linestyle="--", linewidth=1.0, color="k", zorder=0)

    # plot ALL models with error bars (same style, no colors)
    for yy, (_, r) in zip(y, out.iterrows()):
        x = float(r["elpd"])
        se = float(r["se_elpd"]) if np.isfinite(r.get("se_elpd", np.nan)) else np.nan

        if np.isfinite(se):
            ax.errorbar(
                x, yy, xerr=se,
                fmt="o",
                color="black", ecolor="black",
                capsize=3, elinewidth=1, capthick=1,
                markersize=6, markeredgewidth=0.8,
                markeredgecolor="black",
                zorder=3
            )
        else:
            ax.plot(
                x, yy, "o",
                color="black",
                markersize=6,
                markeredgecolor="black",
                markeredgewidth=0.8,
                zorder=3
            )

    ax.set_yticks(y)
    ax.set_yticklabels(out["model"].tolist(), fontsize=8)
    ax.margins(y=0.4)
    ax.invert_yaxis()

    ax.set_xlabel(xlabel)
    if title:
        ax.set_title(title)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)

    plt.tight_layout()
    _save(fig, save, save_path)
    return fig, ax

def plot_delta_total_and_by_treatment(
    delta_df,
    *,
    treatment_order=None,
    save=True,
    save_path=None,
    title=None,
    xlabel="ELPD difference (vMU - vM)",
    dpi=300,
    vline_at_zero=True,
    y_step=1.35,
    y_margins=0.10,
):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    out = delta_df.copy().reset_index(drop=True)

    treat_rows = out[out["row_type"] == "treatment"].copy()
    total_row = out[out["row_type"] == "total"].copy()

    if treatment_order is None:
        treatment_order = (
            treat_rows["treatment"].cat.categories.tolist()
            if len(treat_rows) and hasattr(treat_rows["treatment"], "cat")
            else sorted(treat_rows["treatment"].dropna().unique())
        )

    treat_rows["treatment"] = pd.Categorical(
        treat_rows["treatment"], categories=list(treatment_order), ordered=True
    )
    treat_rows = treat_rows.sort_values("treatment")

    out = pd.concat([total_row, treat_rows], ignore_index=True)

    cmap = mpl.colormaps.get_cmap("tab10")
    treat_to_color = {t: cmap(i % 10) for i, t in enumerate(treatment_order)}

    y = np.arange(len(out), dtype=float) * float(y_step)

    fig_h = max(1.35, (0.20 * y_step) * len(out) + 0.20)
    fig = plt.figure(figsize=(3.0, fig_h), dpi=dpi)
    ax = plt.gca()

    if vline_at_zero:
        ax.axvline(0.0, linestyle="--", linewidth=1.0, color="k", zorder=0)

    for yy, (_, r) in zip(y, out.iterrows()):
        x = float(r["delta_elpd"])
        se = float(r["se_delta_elpd"]) if np.isfinite(r.get("se_delta_elpd", np.nan)) else np.nan

        if r["row_type"] == "total":
            color = "black"
            ecolor = "black"
        else:
            t = r["treatment"]
            color = treat_to_color.get(t, "black")
            ecolor = color

        if np.isfinite(se):
            # error bars + caps (NO marker here)
            ax.errorbar(
                x, yy, xerr=se,
                fmt="none",
                ecolor=ecolor,
                capsize=3, capthick=1,   # <-- caps WILL show
                elinewidth=1,
                zorder=2
            )

        # marker (no edge)
        ax.plot(
            x, yy, "o",
            color=color,
            markersize=6,
            markeredgecolor="none",     # <-- no outline
            markeredgewidth=0,
            zorder=3
        )

    ax.set_yticks(y)
    ax.set_yticklabels(out["label"].tolist(), fontsize=8)
    ax.margins(y=float(y_margins))
    ax.invert_yaxis()

    ax.set_xlabel(xlabel)
    if title:
        ax.set_title(title)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)

    plt.tight_layout()
    _save(fig, save, save_path)
    return fig, ax

def _ax_total_elpd_black(
    summary_df, ax, *, xlabel="ELPD", vline_at_best=False, show_model_names=True
):
    out = summary_df.copy().reset_index(drop=True)
    y = np.arange(len(out), dtype=float)

    if vline_at_best and "is_best" in out.columns and out["is_best"].any():
        best_elpd = float(out.loc[out["is_best"], "elpd"].iloc[0])
        ax.axvline(best_elpd, linestyle="--", linewidth=1.0, color="k", zorder=0)

    for yy, (_, r) in zip(y, out.iterrows()):
        x = float(r["elpd"])
        se = float(r["se_elpd"]) if np.isfinite(r.get("se_elpd", np.nan)) else np.nan

        if np.isfinite(se):
            ax.errorbar(
                x, yy, xerr=se,
                fmt="o",
                color="black", ecolor="black",
                capsize=3, elinewidth=1, capthick=1,
                markersize=6, markeredgewidth=0.8,
                markeredgecolor="black",
                zorder=3
            )
        else:
            ax.plot(x, yy, "o", color="black", markersize=6, zorder=3)

    ax.set_yticks(y)
    if show_model_names:
        ax.set_yticklabels(out["model"].tolist(), fontsize=8)
    else:
        ax.set_yticklabels([])
        ax.tick_params(axis="y", left=False, length=0)

    ax.margins(y=0.4)
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)

def _ax_delta_bracket(
    delta_df,
    ax,
    *,
    treatment_order=None,
    xlabel="ELPD difference",
    vline_at_zero=True,
    y_step=1.35,
    y_margins=0.10,
    bracket_label="Within-treatment",
    # NEW DEFAULTS:
    bracket_x=0.0,         # touch the y-axis spine line
    bracket_cap=0.020,     # cap length INTO the plot
    bracket_tick=0.030,    # tick length OUT to the label (left)
    label_pad=0.006,       # additional left padding for label (axes fraction)
):
    out = delta_df.copy().reset_index(drop=True)

    treat_rows = out[out["row_type"] == "treatment"].copy()
    total_row = out[out["row_type"] == "total"].copy()

    if treatment_order is None:
        treatment_order = (
            treat_rows["treatment"].cat.categories.tolist()
            if len(treat_rows) and hasattr(treat_rows["treatment"], "cat")
            else sorted(treat_rows["treatment"].dropna().unique())
        )

    treat_rows["treatment"] = pd.Categorical(
        treat_rows["treatment"], categories=list(treatment_order), ordered=True
    )
    treat_rows = treat_rows.sort_values("treatment")
    out = pd.concat([total_row, treat_rows], ignore_index=True)

    n_treat = max(1, len(treatment_order))
    treat_colors = sns.color_palette("colorblind", n_colors=n_treat)
    treat_to_color = {t: treat_colors[i % len(treat_colors)] for i, t in enumerate(treatment_order)}

    y = np.arange(len(out), dtype=float) * float(y_step)

    if vline_at_zero:
        ax.axvline(0.0, linestyle="--", linewidth=1.0, color="k", zorder=0)

    for yy, (_, r) in zip(y, out.iterrows()):
        x = float(r["delta_elpd"])
        se = float(r["se_delta_elpd"]) if np.isfinite(r.get("se_delta_elpd", np.nan)) else np.nan

        color = "black" if r["row_type"] == "total" else treat_to_color.get(r["treatment"], "black")

        if np.isfinite(se):
            ax.errorbar(
                x, yy, xerr=se,
                fmt="none",
                ecolor=color,
                capsize=3, capthick=1,
                elinewidth=1,
                zorder=2
            )
        ax.plot(
            x, yy, "o",
            color=color,
            markersize=6,
            markeredgecolor="none",
            markeredgewidth=0,
            zorder=3
        )

    # Total tick+label (only tick we keep)
    if len(y):
        ax.set_yticks([float(y[0])])
        ax.set_yticklabels(["Total"], fontsize=8)
    else:
        ax.set_yticks([])
        ax.set_yticklabels([])

    # bracket over treatment rows, with endpoints exactly at the treatment row centers (y[1]..y[-1])
    if len(y) >= 2:
        y0 = float(y[1])     # first treatment row center
        y1 = float(y[-1])    # last treatment row center
        ymid = 0.5 * (y0 + y1)

        # x in axes-fraction, y in data coords:
        trans = ax.get_yaxis_transform()

        x = float(bracket_x)       # 0.0 touches y-axis line
        cap = float(bracket_cap)   # into plot (right)
        tick = float(bracket_tick) # out to label (left)

        # square bracket "[" sitting on the y-axis spine line
        ax.plot([x, x], [y0, y1], transform=trans, color="black", lw=1.0, clip_on=False, zorder=10)
        ax.plot([x, x + cap], [y0, y0], transform=trans, color="black", lw=1.0, clip_on=False, zorder=10)
        ax.plot([x, x + cap], [y1, y1], transform=trans, color="black", lw=1.0, clip_on=False, zorder=10)

        # tick from bracket to the label (left)
        ax.plot([x, x - tick], [ymid, ymid], transform=trans, color="black", lw=1.0, clip_on=False, zorder=10)

        # label at the end of the tick
        ax.text(
            x - tick - label_pad, ymid, bracket_label,
            transform=trans,
            ha="right", va="center",
            fontsize=8, color="black",
            clip_on=False, zorder=11
        )

    ax.tick_params(axis="y", length=3)  # keep Total tick visible
    ax.margins(y=float(y_margins))
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)

def plot_psisloo_combined(res_total, delta_df, treatment_order, *, save=True, save_path=None):
    """Two stacked panels: total PSIS-LOO ELPD of both models, and the vMU - vM
    ELPD difference total + within-treatment bracket."""
    # ----------------------
    # PLOT (edit freely)
    # ----------------------
    SHOW_MODEL_NAMES = True

    top_h = 1.25
    y_step = 1.35
    bottom_h = max(1.35, 0.20 * y_step * len(delta_df) + 0.20)
    fig_h = top_h + bottom_h

    fig, (ax1, ax2) = plt.subplots(
        nrows=2,
        figsize=(3.2, fig_h),
        dpi=300,
        gridspec_kw={"height_ratios": [top_h, bottom_h]},
    )

    _ax_total_elpd_black(
        res_total, ax1,
        xlabel="ELPD",
        vline_at_best=False,
        show_model_names=SHOW_MODEL_NAMES,
    )

    _ax_delta_bracket(
        delta_df, ax2,
        treatment_order=treatment_order,
        xlabel="ELPD difference",
        vline_at_zero=True,
        y_step=y_step,
        bracket_label="Within-treatment",
        bracket_x=-0.02,       # touches y-axis line
        bracket_cap=0.020,   # caps into plot
        bracket_tick=0.030,  # tick out to label
        label_pad=0.006,
    )

    plt.tight_layout()

    _save(fig, save, save_path)
    return fig

def plot_delta_subplots_by_treatment(res_ind, treatment_order, *, save=True, save_path=None):
    """
    One row of subplots, one per treatment.
    Within each subplot, y-axis = individuals (sorted by delta_elpd).
    """
    res = res_ind.copy()
    res["treatment"] = pd.Categorical(res["treatment"], categories=list(treatment_order), ordered=True)

    # seaborn colorblind palette (stable order)
    palette = sns.color_palette("colorblind", n_colors=len(treatment_order))
    treat_to_color = {t: palette[i] for i, t in enumerate(treatment_order)}
    k = len(treatment_order)

    # global x-limits incl. errorbars
    global_min, global_max = np.inf, -np.inf
    for t in treatment_order:
        sub = res[res["treatment"] == t].sort_values("delta_elpd", ascending=False)
        if sub.empty:
            continue
        x = sub["delta_elpd"].to_numpy()
        xe = sub["se_delta_elpd"].to_numpy()
        xe = np.where(np.isfinite(xe), xe, 0.0)
        global_min = min(global_min, float(np.min(x - xe)))
        global_max = max(global_max, float(np.max(x + xe)))

    if not np.isfinite(global_min):
        global_min, global_max = -1.0, 1.0
    pad = 0.06 * (global_max - global_min if global_max > global_min else 1.0)
    xlim = (global_min - pad, global_max + pad)

    # compact sizing
    max_n = int(res.groupby("treatment").size().max()) if len(res) else 1
    fig_w = max(5.5, 1.6 * k)
    fig_h = max(2.0, 0.18 * max_n + 0.8)

    fig, axes = plt.subplots(1, k, figsize=(fig_w, fig_h), dpi=300, sharex=True)
    if k == 1:
        axes = [axes]

    marker_size = 3
    cap_size = 2
    err_lw = 0.8
    title_fs = 9
    ytick_fs = 7
    xtick_fs = 8

    for i, t in enumerate(treatment_order):
        ax = axes[i]
        sub = (res[res["treatment"] == t]
               .sort_values("delta_elpd", ascending=False)
               .reset_index(drop=True))

        label = sub["treatment_label"].iloc[0] if len(sub) else str(t)
        ax.set_title(label, fontsize=title_fs)
        ax.axvline(0.0, linestyle="--", linewidth=0.9, color="k", alpha=0.8, zorder=0)

        if len(sub):
            y = np.arange(len(sub))
            x = sub["delta_elpd"].to_numpy()
            xe = sub["se_delta_elpd"].to_numpy()
            color = treat_to_color.get(t, "black")

            for yy, xx, xerr in zip(y, x, xe):
                if np.isfinite(xerr):
                    ax.errorbar(
                        xx, yy, xerr=xerr,
                        fmt="o", capsize=cap_size, elinewidth=err_lw,
                        markersize=marker_size,
                        color=color, markeredgecolor=color,
                        markeredgewidth=0.5, zorder=3
                    )
                else:
                    ax.plot(
                        xx, yy, "o",
                        markersize=marker_size,
                        color=color, markeredgecolor=color,
                        markeredgewidth=0.5, zorder=3
                    )

            ax.set_yticks(y)
            ax.set_yticklabels(sub["individual"].astype(str).tolist(), fontsize=ytick_fs)
            ax.invert_yaxis()
        else:
            ax.set_yticks([])

        ax.set_xlim(*xlim)
        ax.tick_params(axis="x", labelsize=xtick_fs)
        ax.grid(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        if i != 0:
            ax.tick_params(axis="y", length=0)
            ax.set_ylabel("")

    fig.supxlabel("ELPD difference (vMU - vM)", fontsize=9)
    plt.tight_layout()

    _save(fig, save, save_path)
    return fig

def plot_mean_hdi_horizontal(
    res, value_col, lo_col, hi_col, xlabel, *,
    xlim=None, save=True, save_path=None,
    cap=0.18,
    lw=1.2,
    ms=5,
    palette_name="colorblind",
    show_treatment_labels=True,  # show/hide y tick labels (treatment names)
    show_treatment_ticks=True,   # show/hide y tick marks
    fig_width=None,              # override width (optional)
    fig_width_with_labels=3.3,   # default width when labels/ticks shown
    fig_width_no_labels=2.3,     # default width when labels+ticks hidden
    vline_x=None,                # optional dashed vertical reference
    vline_color="black",
    vline_ls="--",
    vline_lw=1.0,
):
    colors = sns.color_palette(palette_name, n_colors=len(res))
    y = np.arange(len(res), dtype=float)

    # --- auto width: narrower when BOTH labels and ticks are hidden ---
    if fig_width is None:
        if (not show_treatment_labels) and (not show_treatment_ticks):
            fig_width = fig_width_no_labels
        else:
            fig_width = fig_width_with_labels

    fig_height = 0.35 * len(res)
    fig = plt.figure(figsize=(fig_width, fig_height), dpi=DPI)
    ax = plt.gca()

    for yy, (_, r), c in zip(y, res.iterrows(), colors):
        x  = float(r[value_col])
        lo = float(r[lo_col])
        hi = float(r[hi_col])

        ax.hlines(yy, lo, hi, color=c, linewidth=lw, zorder=2)
        ax.vlines([lo, hi], yy - cap, yy + cap, color=c, linewidth=lw, zorder=2)
        ax.plot(x, yy, "o", color=c, markersize=ms, markeredgewidth=0.0, zorder=3)

    if vline_x is not None:
        ax.axvline(float(vline_x), color=vline_color, linestyle=vline_ls, linewidth=vline_lw, zorder=1)

    # keep y-axis spine; hide only ticks/labels if requested
    ax.set_yticks(y)

    if show_treatment_labels:
        ax.set_yticklabels(res["label"].tolist(), fontsize=8)
    else:
        ax.set_yticklabels([""] * len(y))

    ax.tick_params(
        axis="y",
        which="both",
        left=show_treatment_ticks,
        length=3 if show_treatment_ticks else 0,
        labelleft=show_treatment_labels,
    )

    ax.invert_yaxis()
    ax.margins(y=0.12)

    ax.set_xlabel(xlabel)
    if xlim is not None:
        ax.set_xlim(*xlim)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)

    plt.tight_layout()
    _save(fig, save, save_path)
    return fig, ax

def plot_mean_hdi_horizontal_on_ax(
    ax,
    res,
    *,
    value_col,
    lo_col,
    hi_col,
    xlabel,
    xlim=None,
    vline_x=None,
    palette_name="colorblind",
    show_y=False,          # NO treatment labels/ticks
    cap=0.18,
    lw=1.2,
    ms=5,
):
    colors = sns.color_palette(palette_name, n_colors=len(res))
    y = np.arange(len(res), dtype=float)

    for yy, (_, r), c in zip(y, res.iterrows(), colors):
        x  = float(r[value_col])
        lo = float(r[lo_col])
        hi = float(r[hi_col])

        ax.hlines(yy, lo, hi, color=c, linewidth=lw, zorder=2)
        ax.vlines([lo, hi], yy - cap, yy + cap, color=c, linewidth=lw, zorder=2)
        ax.plot(x, yy, "o", color=c, markersize=ms, markeredgewidth=0.0, zorder=3)

    if vline_x is not None:
        ax.axvline(float(vline_x), color="black", linestyle="--", linewidth=1.0, zorder=1)

    ax.set_yticks(y)
    if show_y:
        ax.set_yticklabels(res["label"].tolist(), fontsize=8)
        ax.tick_params(axis="y", left=True, labelleft=True, length=3)
    else:
        ax.set_yticklabels([""] * len(y))
        ax.tick_params(axis="y", left=False, labelleft=False, length=0)

    ax.invert_yaxis()
    ax.margins(y=0.12)

    ax.set_xlabel(xlabel)
    if xlim is not None:
        ax.set_xlim(*xlim)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)

def plot_psisloo_plus_params_2x2(
    delta_df, treatment_order, idata_mix_full, df_track,
    *, hdi_prob=0.95, max_draws=2000, save=True, save_path=None,
):
    """ELPD-difference panel on top of the full-model population summaries
    (mean angle, mixture weight, MRL, directed-component density at 0) in a 2x2 grid."""

    # --- 1) Summarize population-level params by treatment for 2x2 ---
    res_params = posthoc.summarize_pop_params_for_2x2_by_treatment(
        idata_mix_full=idata_mix_full,
        df_fit=df_track,
        treat_col="treatment",
        treat_name_col="treatmentName",
        treatment_order=treatment_order,
        hdi_prob=0.95,
        max_draws=2000,
    )

    # --- 2) Choose x-limits / reference lines ---
    # mu: auto-pad from HDIs
    mu_min = float(np.nanmin(res_params["mu_deg_lo"]))
    mu_max = float(np.nanmax(res_params["mu_deg_hi"]))
    mu_pad = 0.08 * (mu_max - mu_min + 1e-9)
    mu_xlim = (mu_min - mu_pad, mu_max + mu_pad)

    # mixture weight, MRL
    p_xlim = (0.0, 1.0)
    mrl_xlim = (0.0, 1.0)

    # directed-component predictive density at 0
    dens_xmax = float(np.nanmax(res_params["dens0_dir_hi"]) * 1.08)
    dens_xlim = (0.0, max(0.5, dens_xmax))
    uniform_density = 1.0 / (2.0 * np.pi)

    # --- 3) Figure geometry (keeps left margin wide enough for bracket label) ---
    n_treat = len(res_params)
    y_step = 1.35

    top_h = 1.45
    row_h = 0.35 * n_treat + 0.45     # height per row of the 2x2 block
    fig_h = top_h + 2.0 * row_h

    fig_w = 5  # adjust to taste

    fig = plt.figure(figsize=(fig_w, fig_h), dpi=DPI)

    # Grid: 3 rows x 2 cols (top spans both cols; bottom = 2x2)
    gs = fig.add_gridspec(
        nrows=3, ncols=2,
        height_ratios=[top_h, row_h, row_h],
        hspace=0.55, wspace=0.32
    )

    ax_top = fig.add_subplot(gs[0, :])
    ax_mu   = fig.add_subplot(gs[1, 0])
    ax_p    = fig.add_subplot(gs[1, 1], sharey=ax_mu)
    ax_mrl  = fig.add_subplot(gs[2, 0], sharey=ax_mu)
    ax_den  = fig.add_subplot(gs[2, 1], sharey=ax_mu)

    # IMPORTANT: reserve left margin for the top bracket label so bottom plots
    # "extend to match" the full width of the top axis + bracket label space.
    fig.subplots_adjust(left=0.30, right=0.985, top=0.98, bottom=0.08)

    # --- 4) Top: ELPD difference (was ax2; now ax1/top) ---
    _ax_delta_bracket(
        delta_df, ax_top,
        treatment_order=treatment_order,
        xlabel="ELPD difference",
        vline_at_zero=True,
        y_step=y_step,
        bracket_label="Within-treatment",
        bracket_x=0.0,       # bracket on the y-axis spine
        bracket_cap=0.020,
        bracket_tick=0.055,  # longer tick to label
        label_pad=0.012,
    )

    # --- 5) Bottom: 2x2 parameter summaries (NO treatment labels) ---
    plot_mean_hdi_horizontal_on_ax(
        ax_mu, res_params,
        value_col="mu_deg_mean", lo_col="mu_deg_lo", hi_col="mu_deg_hi",
        xlabel="Mean angle (deg)",
        xlim=mu_xlim,
        vline_x=0.0,          # goal direction reference
        show_y=False,
    )

    plot_mean_hdi_horizontal_on_ax(
        ax_p, res_params,
        value_col="p_mean", lo_col="p_lo", hi_col="p_hi",
        xlabel="Mixture weight",
        xlim=p_xlim,
        vline_x=0.5,
        show_y=False,
    )

    plot_mean_hdi_horizontal_on_ax(
        ax_mrl, res_params,
        value_col="mrl_mean", lo_col="mrl_lo", hi_col="mrl_hi",
        xlabel="MRL (directed component)",
        xlim=mrl_xlim,
        vline_x=0.0,
        show_y=False,
    )

    plot_mean_hdi_horizontal_on_ax(
        ax_den, res_params,
        value_col="dens0_dir_mean", lo_col="dens0_dir_lo", hi_col="dens0_dir_hi",
        xlabel="Predictive density",
        xlim=dens_xlim,
        vline_x=uniform_density,   # uniform baseline density
        show_y=False,
    )

    # Optional: slightly reduce clutter on shared y across bottom
    for ax in (ax_p, ax_mrl, ax_den):
        ax.spines["left"].set_visible(True)  # keep spine for alignment

    _save(fig, save, save_path)
    return fig

def _ax_pred_density(
    res,
    ax,
    *,
    xlabel="Predictive density",
    xlim=None,
    vline_x=None,
    title=None,
    title_fontsize=8,
    show_treatment_labels=True,
    show_treatment_ticks=True,
    y_step=1.35,
    y_margins=0.10,
    cap=0.18,
    lw=1.0,
    ms=6,
    palette_name="colorblind",
):
    out = res.copy().reset_index(drop=True)
    y = np.arange(len(out), dtype=float) * float(y_step)

    colors = sns.color_palette(palette_name, n_colors=len(out))

    if vline_x is not None:
        ax.axvline(float(vline_x), linestyle="--", linewidth=1.0, color="k", zorder=0)

    for yy, (_, r), c in zip(y, out.iterrows(), colors):
        x  = float(r["dens0_mean"])
        lo = float(r["dens0_lo"])
        hi = float(r["dens0_hi"])

        ax.hlines(yy, lo, hi, color=c, linewidth=lw, zorder=2)
        ax.vlines([lo, hi], yy - cap, yy + cap, color=c, linewidth=lw, zorder=2)
        ax.plot(
            x, yy, "o",
            color=c,
            markersize=ms,
            markeredgecolor="none",
            markeredgewidth=0,
            zorder=3
        )

    ax.set_yticks(y)
    if show_treatment_labels:
        ax.set_yticklabels(out["label"].tolist(), fontsize=8)
    else:
        ax.set_yticklabels([""] * len(y))

    ax.tick_params(
        axis="y",
        which="both",
        left=show_treatment_ticks,
        length=3 if show_treatment_ticks else 0,
        labelleft=show_treatment_labels,
    )

    ax.margins(y=float(y_margins))
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    if xlim is not None:
        ax.set_xlim(*xlim)

    if title is not None:
        ax.set_title(title, fontsize=title_fontsize, pad=3)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)

def plot_predictive_density_vmu_directed_vs_vm(res_full, df_track, idata_vmn, *, save=True, save_path=None):
    """Matched two-panel comparison at the population level: the vMU model's
   Directed-component-only predictive density at 0 vs the von Mises model's."""
    required_dir_cols = ["label", "dens0_dir_mean", "dens0_dir_lo", "dens0_dir_hi"]
    missing_dir_cols = [c for c in required_dir_cols if c not in res_full.columns]
    if missing_dir_cols:
        raise KeyError(
            "res_full is missing directed-density columns needed for matched comparison: "
            f"{missing_dir_cols}. Make sure your full-model summary code produced dens0_dir_*."
        )

    res_full_matched = res_full.copy()
    res_full_matched["dens0_mean"] = res_full_matched["dens0_dir_mean"]
    res_full_matched["dens0_lo"]   = res_full_matched["dens0_dir_lo"]
    res_full_matched["dens0_hi"]   = res_full_matched["dens0_dir_hi"]

    # Optional sanity check message (can remove)
    print("Using vMU directed-component-only predictive density columns for matched comparison:")
    print("  dens0_mean <- dens0_dir_mean")
    print("  dens0_lo   <- dens0_dir_lo")
    print("  dens0_hi   <- dens0_dir_hi")

    # Build vM summaries
    res_treat = posthoc.build_res_treat_from_vonmises_idata(
        idata_vmn,
        df_track=df_track,
        theta0_deg=0.0,
        hdi_prob=0.95,
    )

    # ============================================================
    # TWO-PANEL FIGURE (matched comparison: vMU directed component vs vM)
    # ============================================================

    # knobs (match reference)
    HIDE_TREATMENT_TICKS_AND_LABELS = True
    SHOW_TREATMENT_LABELS = not HIDE_TREATMENT_TICKS_AND_LABELS
    SHOW_TREATMENT_TICKS  = not HIDE_TREATMENT_TICKS_AND_LABELS

    y_step = 1.35
    FIG_W = 3.2

    # uniform reference density line
    uniform_density_at_0 = 1.0 / (2.0 * np.pi)

    # shared x-limits across panels (include uniform line)
    xmin = min(
        float(np.nanmin(res_full_matched["dens0_lo"])),
        float(np.nanmin(res_treat["dens0_lo"])),
        uniform_density_at_0,
    )
    xmax = max(
        float(np.nanmax(res_full_matched["dens0_hi"])),
        float(np.nanmax(res_treat["dens0_hi"])),
        uniform_density_at_0,
    )
    pad = 0.10 * (xmax - xmin if xmax > xmin else 1.0)
    shared_xlim = (xmin - pad, xmax + pad)

    # panel heights
    top_h = max(1.35, 0.20 * y_step * len(res_full_matched) + 0.20)
    bottom_h = max(1.35, 0.20 * y_step * len(res_treat) + 0.20)
    fig_h = top_h + bottom_h

    fig, (ax1, ax2) = plt.subplots(
        nrows=2,
        figsize=(FIG_W, fig_h),
        dpi=DPI,
        sharex=True,
        gridspec_kw={"height_ratios": [top_h, bottom_h]},
    )

    # Top panel: vMU model (DIRECTED COMPONENT ONLY, matched comparison)
    _ax_pred_density(
        res_full_matched,
        ax1,
        xlabel="Predictive density",
        xlim=shared_xlim,
        vline_x=uniform_density_at_0,
        title="von Mises-Uniform (directed component only)",
        title_fontsize=8,
        show_treatment_labels=SHOW_TREATMENT_LABELS,
        show_treatment_ticks=SHOW_TREATMENT_TICKS,
        y_step=y_step,
        y_margins=0.10,
        cap=0.18,
        lw=1.0,
        ms=6,
        palette_name="colorblind",
    )

    # Bottom panel: vM model (von Mises)
    _ax_pred_density(
        res_treat,
        ax2,
        xlabel="Predictive density",
        xlim=shared_xlim,
        vline_x=uniform_density_at_0,
        title="von Mises",
        title_fontsize=8,
        show_treatment_labels=SHOW_TREATMENT_LABELS,
        show_treatment_ticks=SHOW_TREATMENT_TICKS,
        y_step=y_step,
        y_margins=0.10,
        cap=0.18,
        lw=1.0,
        ms=6,
        palette_name="colorblind",
    )

    plt.tight_layout()
    _save(fig, save, save_path)
    return fig

def plot_psisloo_plus_predictive_density_3panel(res_total, delta_df, treatment_order, res_full, *, save=True, save_path=None):
    """Three stacked panels: total ELPD, vMU - vM ELPD difference, and the vMU
    mixture predictive density at 0 by treatment."""
    SHOW_MODEL_NAMES = True

    # predictive panel labels/ticks
    HIDE_PRED_TREATMENT_TICKS_AND_LABELS = True
    SHOW_PRED_LABELS = not HIDE_PRED_TREATMENT_TICKS_AND_LABELS
    SHOW_PRED_TICKS  = not HIDE_PRED_TREATMENT_TICKS_AND_LABELS

    FIG_W = 3.2
    DPI = 300
    y_step = 1.35

    # NEW: larger x-axis label/tick fonts (apply to all three panels)
    X_LABEL_FS = 10
    X_TICK_FS  = 9

    # uniform reference density (for mixture predictive density panel)
    uniform_density_at_0 = 1.0 / (2.0 * np.pi)

    # x-limits for predictive panel (include uniform line)
    xmin_pred = min(float(np.nanmin(res_full["dens0_lo"])), uniform_density_at_0)
    xmax_pred = max(float(np.nanmax(res_full["dens0_hi"])), uniform_density_at_0)
    pad_pred = 0.10 * (xmax_pred - xmin_pred if xmax_pred > xmin_pred else 1.0)
    xlim_pred = (xmin_pred - pad_pred, xmax_pred + pad_pred)

    # ----------------------
    # panel heights (same style/formula as your ELPD figure)
    # ----------------------
    h1 = 1.25  # total ELPD panel
    h2 = max(1.35, 0.20 * y_step * len(delta_df) + 0.20)   # ELPD delta panel
    h3 = max(1.35, 0.20 * y_step * len(res_full) + 0.20)   # predictive density panel
    fig_h = h1 + h2 + h3

    fig, (ax1, ax2, ax3) = plt.subplots(
        nrows=3,
        figsize=(FIG_W, fig_h),
        dpi=DPI,
        gridspec_kw={"height_ratios": [h1, h2, h3]},
    )

    # ----------------------
    # Panel 1: Total ELPD
    # ----------------------
    _ax_total_elpd_black(
        res_total, ax1,
        xlabel="ELPD",
        vline_at_best=False,
        show_model_names=SHOW_MODEL_NAMES,
    )

    # ----------------------
    # Panel 2: ELPD difference (total + within-treatment)
    # ----------------------
    _ax_delta_bracket(
        delta_df, ax2,
        treatment_order=treatment_order,
        xlabel="ELPD difference",
        vline_at_zero=True,
        y_step=y_step,
        bracket_label="Within-treatment",
        bracket_x=-0.02,
        bracket_cap=0.020,
        bracket_tick=0.030,
        label_pad=0.006,
    )

    # ----------------------
    # Panel 3: Predictive density (mixture model only; no title)
    # ----------------------
    _ax_pred_density(
        res_full,
        ax3,
        xlabel="Predictive density",
        xlim=xlim_pred,
        vline_x=uniform_density_at_0,
        title=None,  # <- no subplot title
        show_treatment_labels=SHOW_PRED_LABELS,
        show_treatment_ticks=SHOW_PRED_TICKS,
        y_step=y_step,
        y_margins=0.10,
        cap=0.18,
        lw=1.0,
        ms=6,
        palette_name="colorblind",
    )

    # ----------------------
    # NEW: increase x-axis label + tick font sizes across all panels
    # ----------------------
    for ax in (ax1, ax2, ax3):
        ax.xaxis.label.set_size(X_LABEL_FS)
        ax.tick_params(axis="x", which="major", labelsize=X_TICK_FS)
        ax.tick_params(axis="x", which="minor", labelsize=X_TICK_FS)

    plt.tight_layout()

    _save(fig, save, save_path)
    return fig

def plot_individual_subplots_by_treatment(
    res_ind,
    treatment_order,
    *,
    value_col,
    lo_col,
    hi_col,
    xlabel,
    pop_df=None,
    pop_mean_col=None,
    pop_lo_col=None,
    pop_hi_col=None,
    xlim=None,
    save=True,
    save_path=None,

    # optional dashed reference line (uniform / 0.5)
    ref_vline_x=None,
    ref_vline_color="black",
    ref_vline_ls="--",
    ref_vline_lw=1.0,

    # population mean line (solid black)
    pop_mean_line_color="black",
    pop_mean_line_ls="-",
    pop_mean_line_lw=1.0,

    # style
    palette_name="colorblind",
    point_ms=3.0,
    hdi_lw=0.9,
    capsize_pts=2.5,   # cap size in points (consistent across panels)
    title_fs=9,
    ytick_fs=7,
    xtick_fs=8,

):
    """
    One row of subplots, one per treatment.
    Within each subplot: individuals on y-axis, horizontal mean ± HDI.
    Optional population HDI band + solid population mean line.
    Optional dashed reference line (e.g., 0.5 for mixture weight, 1/(2π) for density).

    Notes:
    - Error-bar caps are drawn in points (constant visual size).
    - Individual intervals/points are drawn with clip_on=False so caps/points
      are not visually cut off at fixed x-limits.
    """
    res = res_ind.copy()

    # robust treatment categorical ordering
    res["treatment"] = pd.Categorical(
        res["treatment"],
        categories=list(treatment_order),
        ordered=True
    )

    # seaborn colorblind palette (stable treatment order)
    colors = sns.color_palette(palette_name, n_colors=len(treatment_order))
    treat_to_color = {t: colors[i] for i, t in enumerate(treatment_order)}
    k = len(treatment_order)

    # infer label columns
    treat_label_col = "label" if "label" in res.columns else ("treatment_label" if "treatment_label" in res.columns else None)
    indiv_col = "individual" if "individual" in res.columns else ("tag" if "tag" in res.columns else None)

    # global x-limits (include all individual intervals + population band + ref line)
    gmin, gmax = np.inf, -np.inf
    for t in treatment_order:
        sub = res.loc[res["treatment"] == t].copy()
        if len(sub):
            sub = sub.sort_values(value_col, ascending=False)
            x = pd.to_numeric(sub[value_col], errors="coerce").to_numpy(float)
            lo = pd.to_numeric(sub[lo_col], errors="coerce").to_numpy(float)
            hi = pd.to_numeric(sub[hi_col], errors="coerce").to_numpy(float)
            finite = np.isfinite(x) & np.isfinite(lo) & np.isfinite(hi)
            if np.any(finite):
                gmin = min(gmin, float(np.min(lo[finite])))
                gmax = max(gmax, float(np.max(hi[finite])))

        if pop_df is not None and pop_mean_col and pop_lo_col and pop_hi_col:
            prow = pop_df.loc[pop_df["treatment"] == t]
            if len(prow):
                plo = float(prow.iloc[0][pop_lo_col])
                phi = float(prow.iloc[0][pop_hi_col])
                gmin = min(gmin, plo)
                gmax = max(gmax, phi)

    if ref_vline_x is not None:
        gmin = min(gmin, float(ref_vline_x))
        gmax = max(gmax, float(ref_vline_x))

    if not np.isfinite(gmin) or not np.isfinite(gmax):
        gmin, gmax = -1.0, 1.0

    if xlim is None:
        pad = 0.06 * (gmax - gmin if gmax > gmin else 1.0)
        xlim_use = (gmin - pad, gmax + pad)
    else:
        xlim_use = xlim

    # compact sizing
    max_n = 1
    if len(res):
        counts = res.groupby("treatment", observed=False).size()
        if len(counts):
            max_n = int(counts.max())
    fig_w = max(5.5, 1.6 * k)
    fig_h = max(2.0, 0.18 * max_n + 0.8)

    fig, axes = plt.subplots(1, k, figsize=(fig_w, fig_h), dpi=DPI, sharex=True)
    if k == 1:
        axes = [axes]

    for i, t in enumerate(treatment_order):
        ax = axes[i]
        sub = res.loc[res["treatment"] == t].copy().sort_values(value_col, ascending=False).reset_index(drop=True)

        # title label
        if len(sub) and (treat_label_col is not None):
            title = str(sub[treat_label_col].iloc[0])
        else:
            title = str(t)
        ax.set_title(title, fontsize=title_fs)

        # dashed reference line (uniform / 0.5)
        if ref_vline_x is not None:
            ax.axvline(
                float(ref_vline_x),
                color=ref_vline_color,
                linestyle=ref_vline_ls,
                linewidth=ref_vline_lw,
                zorder=0,
            )

        # population band + solid mean line
        if pop_df is not None and pop_mean_col and pop_lo_col and pop_hi_col:
            prow = pop_df.loc[pop_df["treatment"] == t]
            if len(prow):
                pm  = float(prow.iloc[0][pop_mean_col])
                plo = float(prow.iloc[0][pop_lo_col])
                phi = float(prow.iloc[0][pop_hi_col])

                ax.axvspan(plo, phi, color="0.85", alpha=1.0, zorder=0.5)
                ax.axvline(
                    pm,
                    color=pop_mean_line_color,
                    linestyle=pop_mean_line_ls,
                    linewidth=pop_mean_line_lw,
                    zorder=1,
                )

        # individual intervals + points
        if len(sub):
            y = np.arange(len(sub), dtype=float)
            x = pd.to_numeric(sub[value_col], errors="coerce").to_numpy(float)
            lo = pd.to_numeric(sub[lo_col], errors="coerce").to_numpy(float)
            hi = pd.to_numeric(sub[hi_col], errors="coerce").to_numpy(float)
            c = treat_to_color.get(t, "black")

            for yy, xx, xlo, xhi in zip(y, x, lo, hi):
                x_ok  = np.isfinite(xx)
                lo_ok = np.isfinite(xlo)
                hi_ok = np.isfinite(xhi)

                if x_ok and lo_ok and hi_ok:
                    # asymmetric x-errors from mean to HDI bounds
                    xerr_left  = max(0.0, xx - xlo)
                    xerr_right = max(0.0, xhi - xx)

                    eb = ax.errorbar(
                        xx, yy,
                        xerr=np.array([[xerr_left], [xerr_right]]),
                        fmt="o",
                        color=c,
                        ecolor=c,
                        markersize=point_ms,
                        markeredgecolor=c,
                        markeredgewidth=0.4,
                        elinewidth=hdi_lw,
                        capsize=capsize_pts,   # points => consistent cap size
                        capthick=hdi_lw,
                        zorder=3,
                    )

                    # Prevent visible clipping at fixed x-lims (keep limits unchanged)
                    try:
                        # marker line
                        if eb.lines and eb.lines[0] is not None:
                            eb.lines[0].set_clip_on(False)
                        # cap lines
                        if len(eb.lines) > 1:
                            for capline in eb.lines[1]:
                                capline.set_clip_on(False)
                        # errorbar line collections
                        if len(eb.lines) > 2:
                            for barcol in eb.lines[2]:
                                barcol.set_clip_on(False)
                    except Exception:
                        pass

                elif x_ok:
                    pt = ax.plot(
                        xx, yy, "o",
                        markersize=point_ms,
                        color=c, markeredgecolor=c,
                        markeredgewidth=0.4, zorder=3
                    )[0]
                    pt.set_clip_on(False)

            ax.set_yticks(y)
            if indiv_col is not None:
                ax.set_yticklabels(sub[indiv_col].astype(str).tolist(), fontsize=ytick_fs)
            else:
                ax.set_yticklabels([str(j) for j in range(len(sub))], fontsize=ytick_fs)
            ax.invert_yaxis()
        else:
            ax.set_yticks([])

        ax.set_xlim(*xlim_use)  # keep exact user-requested limits
        ax.tick_params(axis="x", labelsize=xtick_fs)
        ax.grid(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        if i != 0:
            ax.tick_params(axis="y", length=0)
            ax.set_ylabel("")

    fig.supxlabel(xlabel, fontsize=9)
    plt.tight_layout()

    _save(fig, save, save_path)
    return fig

def plot_ppc_by_treatment(
    treat_levels, angles, bin_centers, obs_boot, ppc, treatment_label_map=None,
    *, save=True, save_path=None,
):
    """Polar population-level PPC per treatment: bootstrap KDE bands of the
    observed goal-relative headings vs the posterior-predictive KDE bands."""
    EDGE_RADIUS = 0.55
    BAND_ALPHA = 0.25
    HDI_PROB = 0.95
    k = len(treat_levels)
    fig, axes = plt.subplots(
        1, k,
        figsize=(0.95 * k, 1.25),
        dpi=DPI,
        subplot_kw={"projection": "polar"},
        squeeze=True
    )
    if k == 1:
        axes = [axes]

    # --- swapped tab10 -> seaborn "colorblind" (keep cycling behavior) ---
    base_palette = sns.color_palette("colorblind", 10)
    colors = [base_palette[i % 10] for i in range(k)]

    # ---- convert theta to [0, 2pi) and sort for plotting ----
    theta_plot = (angles + 2 * np.pi) % (2 * np.pi)
    order = np.argsort(theta_plot)
    theta_plot = theta_plot[order]

    bins_plot = (bin_centers + 2 * np.pi) % (2 * np.pi)
    b_order = np.argsort(bins_plot)
    bins_plot = bins_plot[b_order]

    # ---- close the seam at 0° to avoid a thin renderer gap ----
    def _close_circle(theta, y):
        theta_c = np.concatenate([theta, [2.0 * np.pi]])
        y_c     = np.concatenate([y,     [y[0]]])
        return theta_c, y_c

    # Determine global density max (for consistent r scaling)
    global_hi = 0.0
    for t in treat_levels:
        _, o_lo, o_hi = summarize_band(obs_boot[t]["kde"], HDI_PROB)
        _, p_lo, p_hi = summarize_band(ppc[t]["kde"], HDI_PROB)
        global_hi = max(global_hi, float(np.nanmax(o_hi)), float(np.nanmax(p_hi)))

        _, _, o_hhi = summarize_band(obs_boot[t]["hist"], HDI_PROB)
        _, _, p_hhi = summarize_band(ppc[t]["hist"], HDI_PROB)
        global_hi = max(global_hi, float(np.nanmax(o_hhi)), float(np.nanmax(p_hhi)))

    # outer ring corresponds roughly to the upper density envelope (leave headroom)
    outer_ring = EDGE_RADIUS + 0.5
    rmax = EDGE_RADIUS + 0.5 * 1.15

    for ax, t, col in zip(axes, treat_levels, colors):
        format_polar(ax, edge_radius=EDGE_RADIUS, outer_radius=outer_ring)
        ax.set_ylim(0, rmax)

        # ---- Observed bootstrap KDE (black dashed + band) ----
        o_med, o_lo, o_hi = summarize_band(obs_boot[t]["kde"], HDI_PROB)
        o_med, o_lo, o_hi = o_med[order], o_lo[order], o_hi[order]

        theta_c, o_med_c = _close_circle(theta_plot, o_med)
        _,       o_lo_c  = _close_circle(theta_plot, o_lo)
        _,       o_hi_c  = _close_circle(theta_plot, o_hi)

        ax.fill_between(
            theta_c, EDGE_RADIUS + o_lo_c, EDGE_RADIUS + o_hi_c,
            color="black", alpha=BAND_ALPHA, edgecolor="none", zorder=1
        )
        ax.plot(
            theta_c, EDGE_RADIUS + o_med_c,
            color="black", ls="--", lw=LINE_LW, zorder=3
        )

        # ---- PPC KDE (color solid + band) ----
        p_med, p_lo, p_hi = summarize_band(ppc[t]["kde"], HDI_PROB)
        p_med, p_lo, p_hi = p_med[order], p_lo[order], p_hi[order]

        theta_c, p_med_c = _close_circle(theta_plot, p_med)
        _,       p_lo_c  = _close_circle(theta_plot, p_lo)
        _,       p_hi_c  = _close_circle(theta_plot, p_hi)

        ax.fill_between(
            theta_c, EDGE_RADIUS + p_lo_c, EDGE_RADIUS + p_hi_c,
            color=col, alpha=BAND_ALPHA, edgecolor="none", zorder=2
        )
        ax.plot(
            theta_c, EDGE_RADIUS + p_med_c,
            color=col, lw=LINE_LW, zorder=4
        )

        # # ---- Observed bootstrap histogram (points + error bars) ----
        # h_med, h_lo, h_hi = summarize_band(obs_boot[t]["hist"], HDI_PROB)
        # h_med, h_lo, h_hi = h_med[b_order], h_lo[b_order], h_hi[b_order]
        # yerr = np.vstack([h_med - h_lo, h_hi - h_med])
        #
        # # Optional seam-close for histogram too:
        # bins_c, h_med_c = _close_circle(bins_plot, h_med)
        # _,      h_lo_c  = _close_circle(bins_plot, h_lo)
        # _,      h_hi_c  = _close_circle(bins_plot, h_hi)
        # yerr_c = np.vstack([h_med_c - h_lo_c, h_hi_c - h_med_c])
        #
        # ax.errorbar(
        #     bins_c,
        #     EDGE_RADIUS + h_med_c,
        #     yerr=yerr_c,
        #     fmt="o",
        #     ecolor="black",
        #     markeredgecolor="black",
        #     markerfacecolor="black",
        #     lw=0.8,
        #     capsize=2,
        #     markersize=2,
        #     zorder=5
        # )

        title = str(t)
        if treatment_label_map is not None:
            title = format_treatment_title(t, treatment_label_map)

        ax.set_title(title, pad=6)

    # Shared legend
    handles = [
        mpl.lines.Line2D([0], [0], color="black", ls="--", lw=1.2, label="Data"),
        mpl.lines.Line2D([0], [0], color="black", lw=1.2, label="Model"),
    ]

    fig.legend(handles=handles, loc="center", bbox_to_anchor=(0.5, 0.02),
               ncol=3, frameon=True, handlelength=1.6, columnspacing=1.2)

    fig.subplots_adjust(left=0.03, right=0.99, top=0.93, bottom=0.18, wspace=0.22)

    _save(fig, save, save_path)
    return fig

def plot_individual_ppc_grid(ppc, *, save=True, save_path=None):
    """Polar grid (one panel per tag) of each animal's bootstrap KDE bands vs
    its posterior-predictive KDE bands; ``ppc`` is the dict from
    ``posthoc.individual_ppc``."""
    results = ppc["results"]
    angles = ppc["angles"]
    to_r = ppc["to_r"]
    RMAX = ppc["RMAX"]
    EDGE_RADIUS = ppc["EDGE_RADIUS"]
    AMP = ppc["AMP"]
    treat_levels = ppc["treat_levels"]
    treat_to_color = ppc["treat_to_color"]
    nice_treat_label = ppc["nice_treat_label"]
    treat_to_color = {
        t: c for t, c in zip(
            treat_levels,
            sns.color_palette("colorblind", n_colors=len(treat_levels))
        )
    }

    def _polar_format_grid(ax, rmax):
        ax.spines["polar"].set_visible(False)
        ax.grid(False)
        ax.set_yticks([])
        ax.set_yticklabels([])

        # Only show 0° label, with NO tick mark/line
        ax.set_xticks([0.0])
        ax.set_xticklabels(["0°"])
        ax.tick_params(axis="x", pad=-3, length=0, width=0)

        ax.set_ylim(0, rmax)
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)

        # Rings (outermost removed)
        ax.add_artist(Circle((0, 0), EDGE_RADIUS, transform=ax.transData._b,
                             color="0.80", fill=False, linewidth=0.6))
        ax.add_artist(Circle((0, 0), EDGE_RADIUS + 0.5 * AMP, transform=ax.transData._b,
                             color="0.86", fill=False, linewidth=0.6))

        # Gray radial reference line at 0°
        ax.plot([0, 0], [EDGE_RADIUS, rmax], color="0.85", lw=0.8, zorder=0)

    # ============================================================
    # Plot grid
    # ============================================================
    n = len(results)
    ncols = N_COLS
    nrows = int(np.ceil(n / ncols))

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(0.86 * ncols, 0.86 * nrows),
        dpi=DPI,
        subplot_kw={"projection": "polar"},
        squeeze=False
    )

    # Delete unused axes
    for k in range(n, nrows * ncols):
        fig.delaxes(axes[k // ncols, k % ncols])

    for i, r in enumerate(results):
        ax = axes[i // ncols, i % ncols]
        _polar_format_grid(ax, rmax=RMAX)

        treat = r["treat"]
        col = treat_to_color.get(treat, "tab:blue")

        # Observed bootstrap KDE
        ax.fill_between(
            angles, to_r(r["obs_lo"]), to_r(r["obs_hi"]),
            facecolor="0.2", edgecolor="none", alpha=0.18, zorder=1
        )
        ax.plot(angles, to_r(r["obs_med"]), color="0.05", lw=0.9, ls="--", zorder=3)

        # PPC KDE
        ax.fill_between(
            angles, to_r(r["ppc_lo"]), to_r(r["ppc_hi"]),
            facecolor=col, edgecolor="none", alpha=0.28, zorder=2
        )
        ax.plot(angles, to_r(r["ppc_med"]), color=col, lw=1.0, ls="-", zorder=4)

        # Title (tight)
        ax.text(
            0.5, 1.2, f"{r['tag']}",
            transform=ax.transAxes, ha="center", va="bottom",
            fontsize=6, fontweight="bold",
        )

    # Legend
    handles_treat = [Line2D([0], [0], color=treat_to_color[t], lw=2.0) for t in treat_levels]
    labels_treat  = [nice_treat_label[t] for t in treat_levels]

    handles = [
        Line2D([0], [0], color="black", lw=1.6, ls="--"),
        Line2D([0], [0], color="black", lw=1.6, ls="-"),
    ] + handles_treat
    labels = ["Data", "Model"] + labels_treat

    fig.legend(
        handles, labels,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.012),
        ncol=min(6, len(labels)),
        frameon=True,
        handlelength=1.6,
        columnspacing=0.9
    )

    fig.subplots_adjust(
        left=0.015, right=0.985, top=0.985, bottom=0.085,
        wspace=0.25, hspace=0.45
    )

    _save(fig, save, save_path)
    return fig

def _polar_format_mix(ax, *, rmax):
    ax.spines["polar"].set_visible(False)
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.grid(False)
    ax.set_yticks([])
    ax.set_yticklabels([])

    ax.plot([0, 0], [EDGE_RADIUS, rmax], color="0.80", lw=0.9, ls="-", zorder=0)

    ax.add_artist(Circle((0, 0), EDGE_RADIUS,       transform=ax.transData._b,
                         color="0.75", fill=False, linewidth=OUTER_LW))
    ax.add_artist(Circle((0, 0), EDGE_RADIUS + AMP, transform=ax.transData._b,
                         color="0.85", fill=False, linewidth=OUTER_LW))

def _vonmises_pdf_mix(theta, mu, kappa):
    return np.exp(kappa * np.cos(theta - mu)) / (2 * np.pi * np.i0(kappa))

def plot_mixture_polar_pdfs(idata_mix_full, df_track, *, save=True, save_path=None):
    """2x2 polar panels of the full vMU mixture's predictive density on the
    goal-relative circle, population level, one panel per treatment band."""
    HDI_PROB = 0.95
    MAX_DRAWS = 4000
    EDGE_RADIUS = 1.0
    AMP = 1.0
    BAND_ALPHA = 0.35
    KAPPA_UNIFORM = 1e-6
    THETA_N = 900
    SHOW_LEGEND = False
    post = idata_mix_full.posterior

    # required base vars
    required = ["mu0", "eta_kappa0", "logit_p0"]
    missing = [v for v in required if v not in post.data_vars]
    if missing:
        raise KeyError(f"Missing {missing}. Available: {list(post.data_vars)[:120]} ...")

    # optional treatment vars (present in full model)
    has_mu_treat    = "mu_treat" in post.data_vars
    has_kappa_treat = "kappa_treat" in post.data_vars
    has_beta_treat  = "beta_treat" in post.data_vars

    # treatment levels/order (reference = first level)
    treat_levels = (
        df_track["treatment"].cat.categories.tolist()
        if hasattr(df_track["treatment"], "cat")
        else sorted(df_track["treatment"].unique())
    )
    ref = treat_levels[0]
    nonref = treat_levels[1:]

    # labels
    nice = {t: _format_treatment_name(t) for t in treat_levels}
    if "treatmentName" in df_track.columns:
        name_map = df_track.groupby("treatment", observed=False)["treatmentName"].first().to_dict()
        for t in treat_levels:
            nice[t] = _format_treatment_name(name_map.get(t, t))

    # panels: include all treatments if any treatment effect exists; else just control
    has_any_treat = (has_mu_treat or has_kappa_treat or has_beta_treat) and (len(nonref) > 0)
    if has_any_treat:
        panel_levels = [ref] + nonref
        panel_labels = ["Control"] + [nice[t] for t in nonref]
    else:
        panel_levels = [ref]
        panel_labels = ["Control"]

    # fixed 2x2 layout supports up to 4 panels
    if len(panel_levels) > 4:
        raise ValueError(
            f"2x2 layout supports up to 4 panels, but got {len(panel_levels)}. "
            "Filter panel_levels or switch to a larger grid."
        )

    # colors: seaborn colorblind (stable ordering by treat_levels)
    palette = sns.color_palette("colorblind", n_colors=max(3, len(treat_levels)))
    treat_to_color = {t: palette[i % len(palette)] for i, t in enumerate(treat_levels)}
    panel_colors = [treat_to_color[t] for t in panel_levels]

    # posterior draws (population-level only)
    rng = np.random.default_rng(123)
    theta = np.linspace(-np.pi, np.pi, THETA_N, endpoint=False)

    mu0 = circular_modulo(_stack_samples_first(post["mu0"]))          # (S,)
    ek0 = _stack_samples_first(post["eta_kappa0"]).reshape(-1)        # (S,)
    lp0 = _stack_samples_first(post["logit_p0"]).reshape(-1)          # (S,)

    if has_mu_treat:
        mu_tr = _stack_samples_first(post["mu_treat"])                # (S, Kmu)
        if mu_tr.ndim == 1:
            mu_tr = mu_tr[:, None]

    if has_kappa_treat:
        kt = _stack_samples_first(post["kappa_treat"])                # (S, Kk)
        if kt.ndim == 1:
            kt = kt[:, None]

    if has_beta_treat:
        bt = _stack_samples_first(post["beta_treat"])                 # (S, Kp)
        if bt.ndim == 1:
            bt = bt[:, None]

    # subsample once
    S = mu0.shape[0]
    if MAX_DRAWS is not None and S > MAX_DRAWS:
        idx = rng.choice(S, size=MAX_DRAWS, replace=False)
        mu0 = mu0[idx]
        ek0 = ek0[idx]
        lp0 = lp0[idx]
        if has_mu_treat:
            mu_tr = mu_tr[idx, :]
        if has_kappa_treat:
            kt = kt[idx, :]
        if has_beta_treat:
            bt = bt[idx, :]
        S = MAX_DRAWS

    # uniform component
    pdf_unif = _vonmises_pdf_mix(theta[None, :], 0.0, KAPPA_UNIFORM)  # (1, T)

    summ = {}
    peak_max = 0.0

    for t, lab in zip(panel_levels, panel_labels):
        is_ref = (t == ref)

        # mu(t): circular add on linear scale then wrap
        mu_t = mu0.copy()
        if (not is_ref) and has_mu_treat:
            j = nonref.index(t)
            mu_t = circular_modulo(mu_t + mu_tr[:, j])

        # kappa(t): eta_kappa0 + kappa_treat then softplus
        ek_t = ek0.copy()
        if (not is_ref) and has_kappa_treat:
            j = nonref.index(t)
            ek_t = ek_t + kt[:, j]
        k_t = softplus(ek_t)

        # p(t): logit_p0 + beta_treat then sigmoid
        lp_t = lp0.copy()
        if (not is_ref) and has_beta_treat:
            j = nonref.index(t)
            lp_t = lp_t + bt[:, j]
        p_goal = sigmoid(lp_t)

        # mixture pdf
        pdf_goal = _vonmises_pdf_mix(theta[None, :], mu_t[:, None], k_t[:, None])  # (S, T)
        pdf_mix  = p_goal[:, None] * pdf_goal + (1.0 - p_goal[:, None]) * pdf_unif  # (S, T)

        med = np.median(pdf_mix, axis=0)
        lo, hi = hdi_along_samples(pdf_mix, prob=HDI_PROB)

        summ[t] = dict(med=med, lo=lo, hi=hi, label=lab)
        peak_max = max(peak_max, float(np.nanmax(hi)))

    def to_r(y):
        return EDGE_RADIUS + AMP * (y / (peak_max + 1e-12))

    RMAX = EDGE_RADIUS + AMP * 1.03

    # ----------------------------
    # plot (fixed 2x2, compact)
    # ----------------------------
    k = len(panel_levels)
    nrows, ncols = 2, 2

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(2.55, 2.35),  # compact overall figure
        dpi=DPI,
        subplot_kw={"projection": "polar"},
        squeeze=False
    )
    axes = axes.ravel()

    # hide unused panels (if k < 4)
    for ax in axes[k:]:
        ax.set_visible(False)

    for ax, t, col, lab in zip(axes[:k], panel_levels, panel_colors, panel_labels):
        _polar_format_mix(ax, rmax=RMAX)
        ax.set_ylim(0, RMAX)

        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)

        res = summ[t]
        ax.fill_between(
            theta, to_r(res["lo"]), to_r(res["hi"]),
            facecolor=col, edgecolor="none", alpha=BAND_ALPHA, zorder=1
        )
        ax.plot(theta, to_r(res["med"]), color=col, lw=LINE_LW, zorder=3)

        # Panel title helps in a 2x2 layout
        # ax.set_title(lab, pad=2, fontsize=FONT_SIZE)

    # legend (optional)
    if SHOW_LEGEND:
        handles = [Line2D([0], [0], color=treat_to_color[t], lw=1.6) for t in panel_levels]
        labels = panel_labels
        fig.legend(
            handles, labels,
            loc="lower center", bbox_to_anchor=(0.5, 0.01),
            ncol=min(2, k), frameon=True, handlelength=1.2, columnspacing=0.8
        )
        bottom = 0.14
    else:
        bottom = 0.05

    fig.subplots_adjust(
        left=0.04, right=0.99, top=0.95, bottom=bottom,
        wspace=0.06, hspace=0.10
    )

    _save(fig, save, save_path)
    return fig

def _pick_order_col(df):
    for c in ["frame", "step", "t", "time", "timestamp"]:
        if c in df.columns:
            return c
    return None

def _ensure_step_index(df, tag_col="tag"):
    d = df.copy()
    order_col = _pick_order_col(d)
    if order_col is None:
        d["_step"] = d.groupby(tag_col, observed=False).cumcount()
    else:
        d = d.sort_values([tag_col, order_col]).copy()
        d["_step"] = d.groupby(tag_col, observed=False).cumcount()
    return d

def _make_points_gdf():
    pts = [{"site": "Roost", "geometry": Point(ROOST_LON, ROOST_LAT)}]
    pts += [{"site": "Translocation", "geometry": Point(lon, lat)} for (lat, lon) in TRANS_POINTS]
    return gpd.GeoDataFrame(pts, crs="EPSG:4326").to_crs(epsg=WEBMERC_EPSG)

def hide_axes_pgoal(ax):
    ax.set_xlabel(""); ax.set_ylabel("")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

def add_scale_bar_dual_pgoal(
    ax,
    total_length_m=20.0,
    n_segments=2,
    segment_colors=("#ffffff", "#000000"),  # BLACK/WHITE
    segment_height_pt=8,
    edgecolor="#000000",                    # BLACK
    edgewidth=0.9,
    label=True,
    label_above=False,
    label_kwargs=None,
    # precise placement:
    ax_xy=(1.0, 0.0),       # anchor in axes coords
    offset_pt=(-22, 40),    # (dx, dy) in points from anchor
    align="right",
    direction="left",
    pad_text_pt=4,
    return_info=False
):
    """
    Two-tone horizontal scale bar (meters, data units).
    Placement uses axes coords + point offset. Rectangles drawn in data coords.
    """
    if label_kwargs is None:
        label_kwargs = dict(
            color="#000000", fontsize=10, ha="center", va="top",
            fontname="Arial", fontstretch="normal"
        )

    fig = ax.figure
    if fig.canvas is not None:
        fig.canvas.draw()

    base = offset_copy(ax.transAxes, fig=fig, x=offset_pt[0], y=offset_pt[1], units="points")

    x_disp, y_disp = base.transform(ax_xy)
    x_anchor, y_anchor = ax.transData.inverted().transform((x_disp, y_disp))

    seg_len = total_length_m / float(max(n_segments, 1))
    if direction == "left":
        x0 = (x_anchor if align == "left" else x_anchor - total_length_m)
        xdir = +1
    elif direction == "right":
        x0 = (x_anchor if align == "left" else x_anchor + total_length_m)
        if align == "right":
            x0 = x_anchor - total_length_m
        xdir = +1
    else:
        raise ValueError("direction must be 'left' or 'right'")

    # Convert bar height (points) to data units using axes pixel height
    bbox = ax.get_window_extent()
    axpix_h = bbox.height
    y_min, y_max = ax.get_ylim()
    pt_to_data_y = (y_max - y_min) / axpix_h
    half_h = 0.5 * segment_height_pt * pt_to_data_y
    yb, yt = y_anchor - half_h, y_anchor + half_h

    colors = (segment_colors * ((n_segments + len(segment_colors) - 1) // len(segment_colors)))[:n_segments]

    patches = []
    for i in range(n_segments):
        xa = x0 + xdir * (i * seg_len)
        xb = x0 + xdir * ((i + 1) * seg_len)
        p = plt.Rectangle(
            (min(xa, xb), yb),
            abs(xb - xa),
            yt - yb,
            facecolor=colors[i],
            edgecolor=edgecolor,
            linewidth=edgewidth,
            zorder=90
        )
        ax.add_patch(p)
        patches.append(p)

    x_center = x0 + xdir * (total_length_m / 2.0)
    txt = None
    if label:
        txt = ax.text(x_center, y_anchor, f"{int(total_length_m)} m", **label_kwargs)
        dy = pad_text_pt if label_above else -pad_text_pt
        txt.set_transform(offset_copy(ax.transData, fig=fig, x=0, y=dy, units="points"))
        txt.set_zorder(91)

    if return_info:
        return {
            "center_data": (x_center, y_anchor),
            "y_anchor_data": y_anchor,
            "segment_height_pt": segment_height_pt,
            "offset_pt": offset_pt,
            "ax_xy": ax_xy,
            "total_length_m": total_length_m,
            "patches": patches,
            "text": txt
        }
    return None

def add_compass_pgoal(
    ax,
    ax_xy=(0.95, 0.12),
    offset_pt=(0, 0),
    arrow_len_pt=30,        # default smaller now
    head_width_pt=9,
    head_length_pt=10,
    tail_width_pt=2.0,
    color="#000000",
    edgecolor="#000000",
    edgewidth=1.0,
    alpha=1.0,
    draw_N=True,
    N_fontsize=10,
    N_offset_pt=(0, 3)
):
    """
    North arrow drawn in AXES coordinates (constant size, upright).
    """
    fig = ax.figure
    if fig.canvas is not None:
        fig.canvas.draw()

    base = offset_copy(ax.transAxes, fig=fig, x=offset_pt[0], y=offset_pt[1], units="points")

    bbox = ax.get_window_extent()
    axpix_h = bbox.height
    dy_axes = (arrow_len_pt * fig.dpi / 72.0) / axpix_h

    x_axes, y_axes = ax_xy
    start = (x_axes, y_axes)
    end   = (x_axes, y_axes + dy_axes)

    style = ArrowStyle.Simple(head_length=head_length_pt, head_width=head_width_pt, tail_width=tail_width_pt)
    patch = FancyArrowPatch(
        start, end,
        transform=ax.transAxes,
        arrowstyle=style,
        mutation_scale=1.0,
        linewidth=edgewidth,
        edgecolor=edgecolor,
        facecolor=color,
        alpha=alpha,
        joinstyle="miter", capstyle="butt",
        zorder=92
    )
    patch.set_transform(base)
    ax.add_patch(patch)

    if draw_N:
        tx = ax.text(
            x_axes, y_axes + dy_axes, "N",
            transform=ax.transAxes, ha="center", va="bottom",
            fontsize=N_fontsize, color="#000000", zorder=93,
            fontname="Arial"
        )
        tx.set_transform(offset_copy(
            ax.transAxes, fig=fig,
            x=N_offset_pt[0] + offset_pt[0],
            y=N_offset_pt[1] + offset_pt[1],
            units="points"
        ))

def _nice_scale_length_m_pgoal(width_m: float) -> float:
    target = max(width_m * 0.2, 1.0)
    exp = np.floor(np.log10(target))
    base = target / (10 ** exp)
    if base < 1.5:
        nice = 1.0
    elif base < 3.5:
        nice = 2.0
    elif base < 7.5:
        nice = 5.0
    else:
        nice = 10.0
    return float(nice * (10 ** exp))

def _add_scale_and_compass_bottom_right_equal_pad(ax, corner_pad_px=22):
    """
    Adds a black/white scalebar + compass to the *bottom-right* of `ax`
    with equal padding to right & bottom borders (in pixels).
    Compass is centered above the scale bar.
    """
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    pad_pt = corner_pad_px * 72.0 / fig.dpi
    seg_h_pt = 8

    # ↓↓↓ make scale text closer to the bar (was 10)
    scale_text_breathe_pt = 6

    # scale length from axis width (meters; EPSG:3857)
    x0, x1 = ax.get_xlim()
    scale_len_m = _nice_scale_length_m_pgoal(abs(x1 - x0))

    label_kwargs = dict(
        color="#000000", fontsize=10, ha="center", va="top",
        fontname="Arial", fontstretch="normal"
    )

    def _draw_scale(offset_y_pt):
        return add_scale_bar_dual_pgoal(
            ax,
            total_length_m=scale_len_m,
            n_segments=2,
            segment_colors=("#ffffff", "#000000"),
            segment_height_pt=seg_h_pt,
            edgecolor="#000000",
            edgewidth=0.9,
            ax_xy=(1.0, 0.0),
            offset_pt=(-pad_pt, offset_y_pt),
            align="right",
            direction="left",
            label=True,
            label_above=False,
            # ↓↓↓ also tighten label offset (was 10 via breathe=10)
            pad_text_pt=scale_text_breathe_pt,
            label_kwargs=label_kwargs,
            return_info=True
        )

    # initial guess; then corrected so label bottom padding matches corner_pad_px
    offset_y_pt = pad_pt + (seg_h_pt / 2.0) + scale_text_breathe_pt + 6.0
    sb = _draw_scale(offset_y_pt)
    fig.canvas.draw()

    txt = sb["text"]
    if txt is not None:
        ax_bbox = ax.get_window_extent(renderer=renderer)
        txt_bbox = txt.get_window_extent(renderer=renderer)
        current_bottom_pad_px = txt_bbox.y0 - ax_bbox.y0
        dy_px = corner_pad_px - current_bottom_pad_px
        dy_pt = dy_px * 72.0 / fig.dpi

        if abs(dy_pt) > 0.25:
            for p in sb["patches"]:
                p.remove()
            txt.remove()
            sb = _draw_scale(offset_y_pt + dy_pt)
            fig.canvas.draw()

    # compass centered above scale bar
    (cx_data, y_anchor_data) = sb["center_data"]
    cx_disp, cy_disp = ax.transData.transform((cx_data, y_anchor_data))

    # ↓↓↓ bring compass slightly closer to the scale bar (was (seg_h/2)+14)
    compass_gap_pt = (seg_h_pt / 2.0)

    cy_disp_up = cy_disp + compass_gap_pt * fig.dpi / 72.0
    cx_axes, cy_axes = ax.transAxes.inverted().transform((cx_disp, cy_disp_up))

    add_compass_pgoal(
        ax,
        ax_xy=(cx_axes, cy_axes),
        offset_pt=(0, 0),
        # ↓↓↓ make arrow smaller (was 42 / 12 / 14 / 2.5)
        arrow_len_pt=20,
        head_width_pt=9,
        head_length_pt=10,
        tail_width_pt=2.0,
        color="#000000",
        edgecolor="#000000",
        edgewidth=1.0,
        draw_N=True,
        # ↓↓↓ slightly smaller N to match
        N_fontsize=10,
        # ↓↓↓ keep N close to arrow head
        N_offset_pt=(0, 3),
    )

def plot_tracks_pgoal_viridis(
    tracks_df,
    *,
    tag_col="tag",
    lat_col="lat",
    lon_col="lon",
    treat_col="treatment",
    treat_name_col="treatmentName",
    p_col="p_goal_post_mean",
    gap_col="gap_before",
    df_prob=None,
    prob_tag_col="tag",
    prob_step_col="_step",
    prob_p_col="p_goal_post_mean",
    # style
    basemap_source=cx.providers.CartoDB.Positron,
    basemap_alpha=1.0,
    lake_facecolor="tab:cyan",
    lake_alpha=0.25,
    lake_edge_lw=1.0,
    seg_lw=0.9,
    seg_alpha=0.92,
    gap_lw=None,
    gap_alpha=0.80,
    gap_linestyle="--",
    # subsampling
    max_traj_per_panel=None,
    seed=123,
    # layout tweaks
    fig_w=8.2,
    fig_h=7.8,
    wspace=0.01,
    hspace=0.03,
    # colorbar tweaks
    cbar_width=0.012,
    cbar_pad=0.010,
    cbar_ticksize=8,
    cbar_labelsize=9,
    save: bool = True,
    save_path=None,
):
    d = tracks_df.copy()
    for c in [tag_col, lat_col, lon_col, treat_col]:
        if c not in d.columns:
            raise KeyError(f"tracks_df missing required column '{c}'")

    d[lat_col] = pd.to_numeric(d[lat_col], errors="coerce")
    d[lon_col] = pd.to_numeric(d[lon_col], errors="coerce")
    d = d[np.isfinite(d[lat_col]) & np.isfinite(d[lon_col])].copy()
    d = d.dropna(subset=[tag_col, treat_col]).copy()
    d = _ensure_step_index(d, tag_col=tag_col)

    # attach p if needed
    if p_col not in d.columns:
        if df_prob is None:
            raise KeyError(
                f"tracks_df missing '{p_col}'. Either compute it on tracks_df "
                f"or pass df_prob (with columns {prob_tag_col}, {prob_step_col}, {prob_p_col})."
            )
        if prob_step_col not in df_prob.columns:
            raise KeyError(f"df_prob missing '{prob_step_col}'. Make sure it has the per-tag step index.")
        tmp = df_prob[[prob_tag_col, prob_step_col, prob_p_col]].copy()
        tmp = tmp.rename(columns={prob_tag_col: tag_col, prob_step_col: "_step", prob_p_col: p_col})
        d = d.merge(tmp, on=[tag_col, "_step"], how="left", validate="many_to_one")

    d[p_col] = pd.to_numeric(d[p_col], errors="coerce")
    d = d[np.isfinite(d[p_col])].copy()
    d[p_col] = d[p_col].clip(0.0, 1.0)

    has_gaps = gap_col in d.columns
    if has_gaps:
        d[gap_col] = pd.to_numeric(d[gap_col], errors="coerce").fillna(0).astype(int)

    # treatment order + labels
    if pd.api.types.is_categorical_dtype(d[treat_col]):
        treat_levels = list(d[treat_col].cat.categories)
    else:
        treat_levels = sorted(d[treat_col].dropna().unique().tolist())
    if len(treat_levels) == 0:
        raise ValueError("No treatment levels found.")

    nice = {t: _format_treatment_name(t) for t in treat_levels}
    if treat_name_col in d.columns:
        name_map = d.groupby(treat_col, observed=False)[treat_name_col].first().to_dict()
        for t in treat_levels:
            nice[t] = _format_treatment_name(name_map.get(t, t))

    # NEW: subplot border colors (seaborn colorblind, same order as treat_levels)
    treat_border_colors = sns.color_palette("colorblind", n_colors=len(treat_levels))

    # project coords to WebMercator
    transformer = Transformer.from_crs("EPSG:4326", f"EPSG:{WEBMERC_EPSG}", always_xy=True)
    xs, ys = transformer.transform(d[lon_col].to_numpy(), d[lat_col].to_numpy())
    d["_x"] = xs
    d["_y"] = ys

    # lake + points
    lake_gdf = gpd.read_file(LAKE_GEOJSON).to_crs(epsg=WEBMERC_EPSG)
    pts_gdf = _make_points_gdf()

    # common bounds
    layers_bounds = []
    if len(lake_gdf):
        layers_bounds.append(lake_gdf.total_bounds)
    if len(pts_gdf):
        layers_bounds.append(pts_gdf.total_bounds)
    if len(d):
        layers_bounds.append(np.array([d["_x"].min(), d["_y"].min(), d["_x"].max(), d["_y"].max()]))

    minx = min(b[0] for b in layers_bounds)
    miny = min(b[1] for b in layers_bounds)
    maxx = max(b[2] for b in layers_bounds)
    maxy = max(b[3] for b in layers_bounds)

    pad_x = (maxx - minx) * 0.06 if maxx > minx else 1000
    pad_y = (maxy - miny) * 0.06 if maxy > miny else 1000
    xlim = (minx - pad_x, maxx + pad_x)
    ylim = (miny - pad_y, maxy + pad_y)

    # layout: 1 row, 4 cols (as in your code)
    ncols = 4
    nrows = 1

    mpl.rcParams.update({
        "font.family": "Arial",
        "font.sans-serif": ["Arial"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    })

    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_w, fig_h), dpi=300)
    axes = np.atleast_1d(axes).ravel()

    cmap = mpl.colormaps.get_cmap("viridis")
    norm = mcolors.Normalize(vmin=0.0, vmax=1.0)

    if gap_lw is None:
        gap_lw = max(0.7, seg_lw * 0.95)

    rng = np.random.default_rng(seed)

    for i, t in enumerate(treat_levels):
        ax = axes[i]
        ax.set_xlim(*xlim); ax.set_ylim(*ylim)

        # basemap
        try:
            cx.add_basemap(ax, source=basemap_source, attribution=False, alpha=basemap_alpha)
        except Exception as e:
            print("Basemap failed (continuing without tiles):", e)

        # lake
        if len(lake_gdf):
            lake_gdf.plot(ax=ax, facecolor=lake_facecolor, edgecolor=lake_facecolor,
                          alpha=lake_alpha, linewidth=lake_edge_lw, zorder=12)
            lake_gdf.boundary.plot(ax=ax, color=lake_facecolor, linewidth=lake_edge_lw + 0.2,
                                   alpha=1.0, zorder=13)

        # roost + trans points
        ref = pts_gdf.loc[pts_gdf["site"] == "Roost"]
        ax.scatter(ref.geometry.x, ref.geometry.y, marker="*", s=200,
                   facecolors="white", edgecolors="black", linewidths=0.8, zorder=50)

        trans = pts_gdf.loc[pts_gdf["site"] == "Translocation"]
        ax.scatter(trans.geometry.x, trans.geometry.y, marker="X", s=120,
                   facecolors="white", edgecolors="black", linewidths=0.8, zorder=50)

        # tracks for this treatment
        dd = d.loc[d[treat_col] == t].copy()
        tags = sorted(dd[tag_col].unique().tolist())

        if (max_traj_per_panel is not None) and (len(tags) > max_traj_per_panel):
            tags = rng.choice(tags, size=max_traj_per_panel, replace=False).tolist()
            dd = dd.loc[dd[tag_col].isin(tags)].copy()

        all_segments, all_values = [], []
        gap_segments, gap_values = [], []

        order_col = _pick_order_col(dd)
        for tg, g in dd.groupby(tag_col, sort=False):
            g = g.sort_values(order_col) if order_col else g.sort_values("_step")

            x = g["_x"].to_numpy(float)
            y = g["_y"].to_numpy(float)
            p = g[p_col].to_numpy(float)

            if len(x) < 2:
                continue

            pts_xy = np.column_stack([x, y])
            segs = np.stack([pts_xy[:-1], pts_xy[1:]], axis=1)
            pseg = 0.5 * (p[:-1] + p[1:])

            if has_gaps and gap_col in g.columns:
                gb = g[gap_col].to_numpy(int)
                cross_gap = (gb[1:] == 1)

                if np.any(~cross_gap):
                    all_segments.append(segs[~cross_gap])
                    all_values.append(pseg[~cross_gap])

                if np.any(cross_gap):
                    gap_segments.append(segs[cross_gap])
                    gap_values.append(pseg[cross_gap])
            else:
                all_segments.append(segs)
                all_values.append(pseg)

        # dashed gap connectors (draw first)
        if len(gap_segments):
            gsegs = np.concatenate(gap_segments, axis=0)
            gvals = np.concatenate(gap_values, axis=0)
            glc = LineCollection(
                gsegs,
                cmap=cmap,
                norm=norm,
                linewidths=gap_lw,
                alpha=gap_alpha,
                linestyles=gap_linestyle,
                zorder=28,
            )
            glc.set_array(gvals)
            ax.add_collection(glc)

        # solid segments (colored by p)
        if len(all_segments):
            segs = np.concatenate(all_segments, axis=0)
            vals = np.concatenate(all_values, axis=0)
            lc = LineCollection(
                segs,
                cmap=cmap,
                norm=norm,
                linewidths=seg_lw,
                alpha=seg_alpha,
                zorder=30
            )
            lc.set_array(vals)
            ax.add_collection(lc)

        ax.set_aspect("equal", adjustable="box")

        # REMOVED: treatment names
        # ax.set_title(nice[t], loc="left", fontsize=11)

        hide_axes_pgoal(ax)

        # NEW: colored box around each subplot (seaborn colorblind, same order)
        border_color = treat_border_colors[i]
        rect = plt.Rectangle(
            (0, 0), 1, 1,
            transform=ax.transAxes,
            fill=False,
            edgecolor=border_color,
            linewidth=2.2,
            zorder=99,
            clip_on=False
        )
        ax.add_patch(rect)

    # off unused axes
    for j in range(len(treat_levels), len(axes)):
        axes[j].axis("off")

    # --- add compass + scale bar to the RIGHTMOST SUBPLOT (bottom-right corner, equal pad) ---
    rightmost_idx = min(ncols - 1, max(len(treat_levels) - 1, 0))
    if axes[rightmost_idx].axison:
        _add_scale_and_compass_bottom_right_equal_pad(axes[rightmost_idx], corner_pad_px=22)

    # --- thin colorbar on the far right ---
    fig.subplots_adjust(
        left=0.01,
        right=0.91,
        top=0.99,
        bottom=0.06,
        wspace=wspace,
        hspace=hspace
    )

    cax_left = 0.91 + cbar_pad
    cax = fig.add_axes([cax_left, 0.20, cbar_width, 0.62])

    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label("Directed-component probability", fontsize=cbar_labelsize)
    cbar.ax.tick_params(labelsize=cbar_ticksize, width=0.6, length=2.0)
    cbar.outline.set_linewidth(0.6)

    # legend (+ optional gap entry)
    legend_handles = [
        Patch(facecolor=lake_facecolor, edgecolor=lake_facecolor, alpha=lake_alpha, label="Lake"),
        Line2D([0], [0], marker="*", color="none", markerfacecolor="white",
               markeredgecolor="black", markersize=11, label="Roost"),
        Line2D([0], [0], marker="X", color="none", markerfacecolor="white",
               markeredgecolor="black", markersize=8, label="Translocation"),
    ]
    if has_gaps and d[gap_col].sum() > 0:
        legend_handles.append(Line2D([0, 1], [0, 0], color="black", lw=gap_lw,
                                     linestyle=gap_linestyle, alpha=gap_alpha, label="Signal Dropout"))

    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=len(legend_handles),
        frameon=True,
        framealpha=0.9,
        bbox_to_anchor=(0.46, 0.10),
        borderpad=1.2,
        labelspacing=0.9,
        columnspacing=1.6,
        handletextpad=0.8,
        handlelength=2.0,
    )

    _save(fig, save, save_path)
    return fig

def _add_scale_and_compass_small_pad(ax, corner_pad_px=14):
    """
    Adds a black/white scale bar + compass to the bottom-right of `ax`
    with equal padding to right & bottom borders (in pixels).
    Compass is centered above the scale bar.
    """
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    pad_pt = corner_pad_px * 72.0 / fig.dpi
    seg_h_pt = 8
    scale_text_breathe_pt = 5  # keep label close to bar

    # scale length from axis width (meters; EPSG:3857)
    x0, x1 = ax.get_xlim()
    scale_len_m = _nice_scale_length_m_pgoal(abs(x1 - x0))

    # BIGGER than before
    label_kwargs = dict(
        color="#000000", fontsize=10, ha="center", va="top",
        fontname="Arial", fontstretch="normal"
    )

    def _draw_scale(offset_y_pt):
        return add_scale_bar_dual_pgoal(
            ax,
            total_length_m=scale_len_m,
            n_segments=2,
            segment_colors=("#ffffff", "#000000"),
            segment_height_pt=seg_h_pt,
            edgecolor="#000000",
            edgewidth=0.8,
            ax_xy=(1.0, 0.0),
            offset_pt=(-pad_pt, offset_y_pt),
            align="right",
            direction="left",
            label=True,
            label_above=False,
            pad_text_pt=scale_text_breathe_pt,
            label_kwargs=label_kwargs,
            return_info=True
        )

    # initial guess; then correct so label bottom padding matches corner_pad_px
    offset_y_pt = pad_pt + (seg_h_pt / 2.0) + scale_text_breathe_pt + 4.0
    sb = _draw_scale(offset_y_pt)
    fig.canvas.draw()

    txt = sb["text"]
    if txt is not None:
        ax_bbox = ax.get_window_extent(renderer=renderer)
        txt_bbox = txt.get_window_extent(renderer=renderer)
        current_bottom_pad_px = txt_bbox.y0 - ax_bbox.y0
        dy_px = corner_pad_px - current_bottom_pad_px
        dy_pt = dy_px * 72.0 / fig.dpi

        if abs(dy_pt) > 0.25:
            for p in sb["patches"]:
                p.remove()
            txt.remove()
            sb = _draw_scale(offset_y_pt + dy_pt)
            fig.canvas.draw()

    # compass centered above scale bar
    (cx_data, y_anchor_data) = sb["center_data"]
    cx_disp, cy_disp = ax.transData.transform((cx_data, y_anchor_data))

    compass_gap_pt = (seg_h_pt / 2.0)
    cy_disp_up = cy_disp + compass_gap_pt * fig.dpi / 72.0
    cx_axes, cy_axes = ax.transAxes.inverted().transform((cx_disp, cy_disp_up))

    add_compass_pgoal(
        ax,
        ax_xy=(cx_axes, cy_axes),
        offset_pt=(0, 0),
        arrow_len_pt=18,
        head_width_pt=8,
        head_length_pt=9,
        tail_width_pt=1.8,
        color="#000000",
        edgecolor="#000000",
        edgewidth=0.9,
        draw_N=True,
        N_fontsize=10,   # BIGGER
        N_offset_pt=(0, 2),
    )

def plot_tracks_pgoal_small_bigfont(
    tracks_df,
    *,
    tag_col="tag",
    lat_col="lat",
    lon_col="lon",
    treat_col="treatment",
    treat_name_col="treatmentName",
    p_col="p_goal_post_mean",
    gap_col="gap_before",
    df_prob=None,
    prob_tag_col="tag",
    prob_step_col="_step",
    prob_p_col="p_goal_post_mean",
    # style
    basemap_source=cx.providers.CartoDB.Positron,
    basemap_alpha=1.0,
    lake_facecolor="tab:cyan",
    lake_alpha=0.25,
    lake_edge_lw=0.8,
    seg_lw=0.8,
    seg_alpha=0.90,
    gap_lw=None,
    gap_alpha=0.75,
    gap_linestyle="--",
    # subsampling
    max_traj_per_panel=None,
    seed=123,
    # layout tweaks (SMALLER than before, still ~4:3)
    fig_w=5.6,
    fig_h=4.2,
    wspace=0.01,
    hspace=0.02,
    # bottom colorbar (horizontal) -- bigger text/thicker bar
    cbar_height=0.034,
    cbar_bottom=0.060,
    cbar_left=0.20,
    cbar_width=0.60,
    cbar_ticksize=8,      # BIGGER
    cbar_labelsize=9,     # BIGGER
    pdf_raster_dpi=300,            # dpi for rasterized track layers in PDF
    svg_raster_dpi=600,  # raster dpi of a saved SVG export
    save: bool = True,
    save_path=None,
):
    d = tracks_df.copy()
    for c in [tag_col, lat_col, lon_col, treat_col]:
        if c not in d.columns:
            raise KeyError(f"tracks_df missing required column '{c}'")

    d[lat_col] = pd.to_numeric(d[lat_col], errors="coerce")
    d[lon_col] = pd.to_numeric(d[lon_col], errors="coerce")
    d = d[np.isfinite(d[lat_col]) & np.isfinite(d[lon_col])].copy()
    d = d.dropna(subset=[tag_col, treat_col]).copy()
    d = _ensure_step_index(d, tag_col=tag_col)

    # attach p if needed
    if p_col not in d.columns:
        if df_prob is None:
            raise KeyError(
                f"tracks_df missing '{p_col}'. Either compute it on tracks_df "
                f"or pass df_prob (with columns {prob_tag_col}, {prob_step_col}, {prob_p_col})."
            )
        if prob_step_col not in df_prob.columns:
            raise KeyError(f"df_prob missing '{prob_step_col}'. Make sure it has the per-tag step index.")
        tmp = df_prob[[prob_tag_col, prob_step_col, prob_p_col]].copy()
        tmp = tmp.rename(columns={prob_tag_col: tag_col, prob_step_col: "_step", prob_p_col: p_col})
        d = d.merge(tmp, on=[tag_col, "_step"], how="left", validate="many_to_one")

    d[p_col] = pd.to_numeric(d[p_col], errors="coerce")
    d = d[np.isfinite(d[p_col])].copy()
    d[p_col] = d[p_col].clip(0.0, 1.0)

    has_gaps = gap_col in d.columns
    if has_gaps:
        d[gap_col] = pd.to_numeric(d[gap_col], errors="coerce").fillna(0).astype(int)

    # treatment order + labels
    if pd.api.types.is_categorical_dtype(d[treat_col]):
        treat_levels = list(d[treat_col].cat.categories)
    else:
        treat_levels = sorted(d[treat_col].dropna().unique().tolist())
    if len(treat_levels) == 0:
        raise ValueError("No treatment levels found.")

    if treat_name_col in d.columns:
        name_map = d.groupby(treat_col, observed=False)[treat_name_col].first().to_dict()
        _ = {t: _format_treatment_name(name_map.get(t, t)) for t in treat_levels}

    # colored subplot borders (seaborn colorblind, stable order)
    treat_border_colors = sns.color_palette("colorblind", n_colors=len(treat_levels))

    # project coords to WebMercator
    transformer = Transformer.from_crs("EPSG:4326", f"EPSG:{WEBMERC_EPSG}", always_xy=True)
    xs, ys = transformer.transform(d[lon_col].to_numpy(), d[lat_col].to_numpy())
    d["_x"] = xs
    d["_y"] = ys

    # lake + points
    lake_gdf = gpd.read_file(LAKE_GEOJSON).to_crs(epsg=WEBMERC_EPSG)
    pts_gdf = _make_points_gdf()

    # common bounds
    layers_bounds = []
    if len(lake_gdf):
        layers_bounds.append(lake_gdf.total_bounds)
    if len(pts_gdf):
        layers_bounds.append(pts_gdf.total_bounds)
    if len(d):
        layers_bounds.append(np.array([d["_x"].min(), d["_y"].min(), d["_x"].max(), d["_y"].max()]))

    minx = min(b[0] for b in layers_bounds)
    miny = min(b[1] for b in layers_bounds)
    maxx = max(b[2] for b in layers_bounds)
    maxy = max(b[3] for b in layers_bounds)

    pad_x = (maxx - minx) * 0.05 if maxx > minx else 1000
    pad_y = (maxy - miny) * 0.05 if maxy > miny else 1000
    xlim = (minx - pad_x, maxx + pad_x)
    ylim = (miny - pad_y, maxy + pad_y)

    # layout: 1 row, 4 cols
    ncols = 4
    nrows = 1

    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_w, fig_h), dpi=DPI)
    axes = np.atleast_1d(axes).ravel()

    cmap = mpl.colormaps.get_cmap("viridis")
    norm = mcolors.Normalize(vmin=0.0, vmax=1.0)

    if gap_lw is None:
        gap_lw = max(0.55, seg_lw * 0.9)

    rng = np.random.default_rng(seed)

    # collect only heavy track artists (solid + dashed collections)
    track_artists = []

    for i, t in enumerate(treat_levels):
        ax = axes[i]
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)

        # basemap (already raster by nature)
        try:
            cx.add_basemap(ax, source=basemap_source, attribution=False, alpha=basemap_alpha)
        except Exception as e:
            print("Basemap failed (continuing without tiles):", e)

        # lake (keep vector)
        if len(lake_gdf):
            lake_gdf.plot(
                ax=ax,
                facecolor=lake_facecolor,
                edgecolor=lake_facecolor,
                alpha=lake_alpha,
                linewidth=lake_edge_lw,
                zorder=12
            )
            lake_gdf.boundary.plot(
                ax=ax,
                color=lake_facecolor,
                linewidth=lake_edge_lw + 0.15,
                alpha=1.0,
                zorder=13
            )

        # roost + trans points (keep vector)
        ref = pts_gdf.loc[pts_gdf["site"] == "Roost"]
        ax.scatter(
            ref.geometry.x, ref.geometry.y,
            marker="*", s=105,
            facecolors="white", edgecolors="black",
            linewidths=0.7, zorder=50
        )

        trans = pts_gdf.loc[pts_gdf["site"] == "Translocation"]
        ax.scatter(
            trans.geometry.x, trans.geometry.y,
            marker="X", s=62,
            facecolors="white", edgecolors="black",
            linewidths=0.7, zorder=50
        )

        # tracks for this treatment
        dd = d.loc[d[treat_col] == t].copy()
        tags = sorted(dd[tag_col].unique().tolist())

        if (max_traj_per_panel is not None) and (len(tags) > max_traj_per_panel):
            tags = rng.choice(tags, size=max_traj_per_panel, replace=False).tolist()
            dd = dd.loc[dd[tag_col].isin(tags)].copy()

        all_segments, all_values = [], []
        gap_segments, gap_values = [], []

        order_col = _pick_order_col(dd)
        for tg, g in dd.groupby(tag_col, sort=False):
            g = g.sort_values(order_col) if order_col else g.sort_values("_step")

            x = g["_x"].to_numpy(float)
            y = g["_y"].to_numpy(float)
            p = g[p_col].to_numpy(float)

            if len(x) < 2:
                continue

            pts_xy = np.column_stack([x, y])
            segs = np.stack([pts_xy[:-1], pts_xy[1:]], axis=1)
            pseg = 0.5 * (p[:-1] + p[1:])

            if has_gaps and gap_col in g.columns:
                gb = g[gap_col].to_numpy(int)
                cross_gap = (gb[1:] == 1)

                if np.any(~cross_gap):
                    all_segments.append(segs[~cross_gap])
                    all_values.append(pseg[~cross_gap])

                if np.any(cross_gap):
                    gap_segments.append(segs[cross_gap])
                    gap_values.append(pseg[cross_gap])
            else:
                all_segments.append(segs)
                all_values.append(pseg)

        # dashed gap connectors (draw first) -> heavy layer candidate for rasterization
        if len(gap_segments):
            gsegs = np.concatenate(gap_segments, axis=0)
            gvals = np.concatenate(gap_values, axis=0)
            glc = LineCollection(
                gsegs,
                cmap=cmap,
                norm=norm,
                linewidths=gap_lw,
                alpha=gap_alpha,
                linestyles=gap_linestyle,
                zorder=28,
            )
            glc.set_array(gvals)
            _style_trajectory_artists([glc])
            ax.add_collection(glc)
            track_artists.append(glc)

        # solid segments (colored by p) -> heavy layer candidate for rasterization
        if len(all_segments):
            segs = np.concatenate(all_segments, axis=0)
            vals = np.concatenate(all_values, axis=0)
            lc = LineCollection(
                segs,
                cmap=cmap,
                norm=norm,
                linewidths=seg_lw,
                alpha=seg_alpha,
                zorder=30
            )
            lc.set_array(vals)
            _style_trajectory_artists([lc])
            ax.add_collection(lc)
            track_artists.append(lc)

        ax.set_aspect("equal", adjustable="box")
        hide_axes_pgoal(ax)

        # colored border per treatment (keep vector)
        rect = plt.Rectangle(
            (0, 0), 1, 1,
            transform=ax.transAxes,
            fill=False,
            edgecolor=treat_border_colors[i],
            linewidth=1.5,
            zorder=99,
            clip_on=False
        )
        ax.add_patch(rect)

    # off unused axes
    for j in range(len(treat_levels), len(axes)):
        axes[j].axis("off")

    # scale bar + compass in rightmost active subplot (keep vector)
    rightmost_idx = min(ncols - 1, max(len(treat_levels) - 1, 0))
    if axes[rightmost_idx].axison:
        _add_scale_and_compass_small_pad(axes[rightmost_idx], corner_pad_px=14)

    # compact layout + horizontal colorbar on bottom
    fig.subplots_adjust(
        left=0.01,
        right=0.99,
        top=0.995,
        bottom=max(0.13, cbar_bottom + cbar_height + 0.045),  # reserve room for larger text
        wspace=wspace,
        hspace=hspace
    )

    cax = fig.add_axes([cbar_left, cbar_bottom, cbar_width, cbar_height])

    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cbar.set_label("Directed-component probability", fontsize=cbar_labelsize, labelpad=2)
    cbar.ax.tick_params(labelsize=cbar_ticksize, width=0.6, length=2.2, pad=1)
    cbar.outline.set_linewidth(0.6)

    # no legend (paper-ready compact version)

    _rasterize = str(save_path).endswith(".pdf")
    _set_rasterized(track_artists, _rasterize)
    _save(fig, save, save_path, dpi=pdf_raster_dpi if _rasterize else svg_raster_dpi)
    return fig

def plot_population_curve_smooth(
    out_dict,
    treat_levels,
    *,
    xlabel,
    ylabel="Population mean p(goal | step)",
    band_alpha=0.20,
    line_lw=1.4,
    fig_w=4.4,
    fig_h=2.6,
    legend_fontsize=7,
    invert_x=False,
    xlim=None,
    x_tick_abs=False,
    save=True,
    save_path=None,
    palette_name="colorblind",
):
    # seaborn colorblind palette (stable treatment order)
    colors = sns.color_palette(palette_name, n_colors=len(treat_levels))
    treat_to_color = {t: colors[i] for i, t in enumerate(treat_levels)}

    fig = plt.figure(figsize=(fig_w, fig_h), dpi=DPI)
    ax = plt.gca()
    ax.axhline(0.5, color="0.90", lw=0.9, ls="--", zorder=0)

    for t in treat_levels:
        if t not in out_dict:
            continue
        s = out_dict[t]
        x = np.asarray(s["x_mid"], float)
        m = np.asarray(s["mean"], float)
        lo = np.asarray(s["lo"], float)
        hi = np.asarray(s["hi"], float)

        order = np.argsort(x)
        x, m, lo, hi = x[order], m[order], lo[order], hi[order]

        c = treat_to_color[t]
        ax.fill_between(x, lo, hi, color=c, alpha=band_alpha, linewidth=0)
        ax.plot(x, m, color=c, lw=line_lw, label=s.get("label", str(t)))

    ax.set_ylim(0, 1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if xlim is not None:
        ax.set_xlim(*xlim)
    if invert_x:
        ax.invert_xaxis()

    if x_tick_abs:
        def _abs_tick_formatter(val, pos):
            aval = abs(val)
            if np.isfinite(aval) and np.isclose(aval, np.round(aval), atol=1e-9):
                return f"{int(np.round(aval)):,}"
            if np.isfinite(aval):
                return f"{aval:.2f}"
            return ""
        ax.xaxis.set_major_formatter(FuncFormatter(_abs_tick_formatter))

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)

    ax.legend(
        frameon=True,
        fontsize=legend_fontsize,
        ncol=1,
        handlelength=1.4,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        borderaxespad=0.0,
    )
    fig.subplots_adjust(right=0.78)

    _save(fig, save, save_path)
    return fig

def plot_population_curve_smooth_compact(
    out_dict,
    treat_levels,
    *,
    xlabel,
    ylabel="Population mean p(goal | step)",
    band_alpha=0.18,
    line_lw=1.2,
    fig_w=3.40,
    fig_h=1.20,
    invert_x=False,
    xlim=None,
    x_tick_abs=False,
    save=True,
    save_path=None,
    palette_name="colorblind",

    # tiny-height friendly text
    axis_label_fs=7,
    tick_label_fs=6,
    y_labelpad=1,
    x_labelpad=1,
):
    colors = sns.color_palette(palette_name, n_colors=len(treat_levels))
    treat_to_color = {t: colors[i] for i, t in enumerate(treat_levels)}

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=DPI)

    ax.axhline(0.5, color="0.88", lw=0.8, ls="--", zorder=0)

    for t in treat_levels:
        if t not in out_dict:
            continue
        s = out_dict[t]
        x = np.asarray(s["x_mid"], float)
        m = np.asarray(s["mean"], float)
        lo = np.asarray(s["lo"], float)
        hi = np.asarray(s["hi"], float)

        order = np.argsort(x)
        x, m, lo, hi = x[order], m[order], lo[order], hi[order]

        c = treat_to_color[t]
        ax.fill_between(x, lo, hi, color=c, alpha=band_alpha, linewidth=0)
        ax.plot(x, m, color=c, lw=line_lw)

    ax.set_ylim(0, 1)
    ax.set_xlabel(xlabel, fontsize=axis_label_fs, labelpad=x_labelpad)
    ax.set_ylabel(ylabel, fontsize=axis_label_fs, labelpad=y_labelpad)

    if xlim is not None:
        ax.set_xlim(*xlim)
    if invert_x:
        ax.invert_xaxis()

    if x_tick_abs:
        def _abs_tick_formatter(val, pos):
            aval = abs(val)
            if np.isfinite(aval) and np.isclose(aval, np.round(aval), atol=1e-9):
                return f"{int(np.round(aval)):,}"
            if np.isfinite(aval):
                return f"{aval:.2f}"
            return ""
        ax.xaxis.set_major_formatter(FuncFormatter(_abs_tick_formatter))

    ax.tick_params(axis="both", which="major", labelsize=tick_label_fs, length=2.5, pad=1)
    ax.tick_params(axis="both", which="minor", labelsize=tick_label_fs, length=2.0, pad=1)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)

    # No legend; use manual margins because tight_layout can struggle at tiny height
    fig.subplots_adjust(left=0.13, right=0.995, bottom=0.30, top=0.98)

    _save(fig, save, save_path)
    return fig

def plot_tag_curves_faceted_4panel(
    out_by_treat,
    treat_levels,
    *,
    xlabel,
    ylabel="p(goal | step)",
    invert_x=False,
    xlim=None,
    x_tick_abs=False,
    band_alpha=0.10,
    line_alpha=0.75,
    line_lw=0.9,
    max_tags_per_panel=None,  # None = all tags; else randomly sample
    seed=123,
    save=True,
    save_path=None,
    fig_w=10.5,
    fig_h=7.6,
    palette_name="colorblind",
):
    rng = np.random.default_rng(seed)

    # seaborn colorblind palette (stable treatment order)
    colors = sns.color_palette(palette_name, n_colors=len(treat_levels))
    treat_to_color = {t: colors[i] for i, t in enumerate(treat_levels)}

    tl = [t for t in treat_levels if t in out_by_treat]
    pages = [tl[i:i+4] for i in range(0, len(tl), 4)]
    figs = []

    for p_i, page in enumerate(pages, start=1):
        fig, axes = plt.subplots(2, 2, figsize=(fig_w, fig_h), dpi=DPI, sharey=True)
        axes = axes.ravel()

        for ax_i, ax in enumerate(axes):
            if ax_i >= len(page):
                ax.axis("off")
                continue

            t = page[ax_i]
            block = out_by_treat[t]
            x = np.asarray(block["x_mid"], float)
            order = np.argsort(x)
            x = x[order]

            ax.axhline(0.5, color="0.90", lw=0.9, ls="--", zorder=0)

            c = treat_to_color[t]
            tags = list(block["tags"].items())
            tags.sort(key=lambda kv: kv[1].get("n_obs", 0), reverse=True)

            if (max_tags_per_panel is not None) and (len(tags) > int(max_tags_per_panel)):
                idx = rng.choice(len(tags), size=int(max_tags_per_panel), replace=False)
                idx = np.sort(idx)
                tags = [tags[i] for i in idx]

            for tg, s in tags:
                m = np.asarray(s["mean"], float)[order]
                lo = np.asarray(s["lo"], float)[order]
                hi = np.asarray(s["hi"], float)[order]

                # NaNs (from truncation) naturally break lines and bands
                ax.fill_between(x, lo, hi, color=c, alpha=band_alpha, linewidth=0)
                ax.plot(x, m, color=c, alpha=line_alpha, lw=line_lw)

            ax.set_ylim(0, 1)

            if xlim is not None:
                ax.set_xlim(*xlim)
            if invert_x:
                ax.invert_xaxis()

            if x_tick_abs:
                def _abs_tick_formatter(val, pos):
                    aval = abs(val)
                    if np.isfinite(aval) and np.isclose(aval, np.round(aval), atol=1e-9):
                        return f"{int(np.round(aval)):,}"
                    if np.isfinite(aval):
                        return f"{aval:,.2f}"
                    return ""
                ax.xaxis.set_major_formatter(FuncFormatter(_abs_tick_formatter))

            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.grid(False)

            # --- font sizes ---
            AXIS_LABEL_FS = 8
            TICK_LABEL_FS = 6
            TITLE_FS      = 11

            ax.set_title(block.get("label", str(t)), fontsize=TITLE_FS)
            ax.set_xlabel(xlabel, fontsize=AXIS_LABEL_FS)
            if ax_i % 2 == 0:
                ax.set_ylabel(ylabel, fontsize=AXIS_LABEL_FS)

            ax.tick_params(axis="both", which="major", labelsize=TICK_LABEL_FS)
            ax.tick_params(axis="both", which="minor", labelsize=TICK_LABEL_FS)

        fig.tight_layout()

        if save:
            root, ext = os.path.splitext(save_path)
            save_page = f"{root}_p{p_i}{ext}" if len(pages) > 1 else save_path
            _save(fig, save, save_page)

        figs.append(fig)

    return figs

def _treatment_plot_df(dfp, value_col):
    """Clean ``dfp`` to finite ``value_col`` values and derive a stable
    treatment order, tab10 colors, and pretty labels."""
    dfp = dfp.copy()
    dfp = dfp.dropna(subset=[value_col, "treatment"])
    dfp[value_col] = pd.to_numeric(dfp[value_col], errors="coerce")
    dfp = dfp[np.isfinite(dfp[value_col])]

    if pd.api.types.is_categorical_dtype(dfp["treatment"]):
        treatments = list(dfp["treatment"].cat.categories)
    else:
        treatments = sorted(dfp["treatment"].unique())

    cmap = mpl.colormaps.get_cmap("tab10")
    treat_to_color = {t: cmap(i % 10) for i, t in enumerate(treatments)}

    if "treatmentName" in dfp.columns:
        name_map = dfp.groupby("treatment")["treatmentName"].first().to_dict()
        treat_label = {t: _format_treatment_name(name_map.get(t, t)) for t in treatments}
    else:
        treat_label = {t: str(t) for t in treatments}

    return dfp, treatments, treat_to_color, treat_label


def _dots_median_interval(ax, dfp, value_col, treatments, treat_to_color,
                          treat_label, lo_q=25, hi_q=75):
    """Dots (individuals) + median and percentile interval per treatment."""
    med_q = 50
    xpos = np.arange(len(treatments))
    jitter = 0.10
    rng = np.random.default_rng(0)

    for i, t in enumerate(treatments):
        sub = dfp[dfp["treatment"] == t]
        y = sub[value_col].to_numpy(dtype=float)
        y = y[np.isfinite(y)]
        if y.size == 0:
            continue

        # dots (individuals)
        x = np.full(y.size, xpos[i], dtype=float) + rng.uniform(-jitter, jitter, size=y.size)
        ax.plot(
            x, y, "o",
            markersize=3.5,
            color=treat_to_color[t],
            markeredgecolor="none",
            markeredgewidth=0.4,
            alpha=0.9,
            zorder=3
        )

        # median + percentile interval
        qlo, qmed, qhi = np.percentile(y, [lo_q, med_q, hi_q])

        ax.errorbar(
            xpos[i], qmed,
            yerr=[[qmed - qlo], [qhi - qmed]],
            fmt="o",
            markersize=5.0,
            color=treat_to_color[t],
            ecolor="black",
            elinewidth=1.0,
            capsize=3,
            markeredgecolor="black",
            markeredgewidth=0.6,
            zorder=4
        )

    ax.set_xticks(xpos)
    ax.set_xticklabels([treat_label[t] for t in treatments], rotation=25, ha="right")


def plot_kappa_hat(kappa_df, *, save=True, save_path=None):
    r"""Estimated population concentration $\hat{\kappa}$ by treatment,
    dots + median + 25--75 percentile band."""
    value_col = "kappa_hat"  # or "Rbar"
    dfp, treatments, treat_to_color, treat_label = _treatment_plot_df(kappa_df, value_col)
    lo_q, hi_q = 25, 75   # change to (5, 95) if you prefer

    fig, ax = plt.subplots(figsize=(2.6, 2.2), dpi=300)
    _dots_median_interval(ax, dfp, value_col, treatments, treat_to_color,
                          treat_label, lo_q, hi_q)
    ax.set_ylabel(r"Estimated $\hat{\kappa}$")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, axis="y", alpha=0.25)
    # plt.ylim(0, 25)
    # plt.yscale("log")
    plt.ylim(0, 10)
    plt.tight_layout()
    _save(fig, save, save_path)
    return fig


def plot_rbar(kappa_df, *, save=True, save_path=None):
    """Mean resultant length $\vec{R}$ by treatment, dots + median +
    25--75 percentile band."""
    value_col = "Rbar"
    dfp, treatments, treat_to_color, treat_label = _treatment_plot_df(kappa_df, value_col)
    lo_q, hi_q = 25, 75   # change to (5, 95) if you prefer

    fig, ax = plt.subplots(figsize=(2.6, 2.2), dpi=300)
    _dots_median_interval(ax, dfp, value_col, treatments, treat_to_color,
                          treat_label, lo_q, hi_q)
    ax.set_ylabel("Mean\nresultant length")
    ax.set_ylim(0, 1.02)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.ylim(-0.05, 1.05)
    plt.tight_layout()
    _save(fig, save, save_path)
    return fig


def plot_alignment(align_df, *, save=True, save_path=None):
    r"""Alignment $\mathbb{E}[\cos(y)]$ of individual mean bearings to
    $0^\circ$ by treatment, dots + median + 25--75 percentile band."""
    value_col = "alignment0"
    dfp, treatments, treat_to_color, treat_label = _treatment_plot_df(align_df, value_col)
    lo_q, hi_q = 25, 75   # change to (5, 95) if you prefer

    fig, ax = plt.subplots(figsize=(2.6, 2.2), dpi=300)
    _dots_median_interval(ax, dfp, value_col, treatments, treat_to_color,
                          treat_label, lo_q, hi_q)
    ax.set_ylabel("Alignment\n$\\mathbb{E}[\\cos(y)]$")

    # alignment is in [-1, 1]
    ax.set_ylim(0, 1.02)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    _save(fig, save, save_path)
    return fig
