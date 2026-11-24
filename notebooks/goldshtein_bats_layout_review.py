"""Render the bat composite from the pinned tracks and model fits.

The flight maps (panels a and d) sit on a vector basemap drawn from the
pinned OpenStreetMap extract in the figure materials (water, roads,
buildings as true vectors — sharp at any zoom; attribution and
provenance in the asset README and the render JSON; see
:func:`draw_track_map` for the live CartoDB raster and plain-face
fallbacks); every other input — tracks, fits, the lake outline, the
artwork — is the pinned material :mod:`goldshtein_bats_layout_data`
validates."""
from __future__ import annotations

import argparse
from functools import lru_cache
import gzip
import json
from pathlib import Path

import goldshtein_bats_figure_core as source
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Patch, PathPatch
from matplotlib.path import Path as MplPath
from matplotlib.textpath import TextToPath
from matplotlib.ticker import MultipleLocator, StrMethodFormatter
from matplotlib.transforms import Affine2D, ScaledTranslation
import numpy as np
from PIL import Image
import pymupdf
from scipy.special import expit
from shapely.geometry.polygon import orient

from goldshtein_bats_layout_data import (
    PANEL_C_HDI_PROB, _population_parameters, _read_context, assemble_layout_data,
)


INK = "#222222"
MAP_FACE = "#f5f6f5"
HEATMAP_NORM = mpl.colors.TwoSlopeNorm(vmin=0, vcenter=.5, vmax=1)
PALETTE = ("#0072B2", "#BC7300", "#008664", "#C95118")
LABELS = ("Control", "Eye cover", "+ Magnet", "+ Anosmia")
PAGE_W, PAGE_H = 532, 890
#: Tile zoom for the live CartoDB raster fallback (used only when the
#: vector extract below is absent): pinned rather than left to
#: contextily's auto-zoom, which matches screen-pixel sizes and on a
#: print page undershoots into visibly blocky tiles.
TILE_ZOOM = 16

#: Pinned OpenStreetMap extract (water/roads/buildings clipped to the map
#: extents) drawn as the vector basemap; © OpenStreetMap contributors
#: (ODbL), Geofabrik israel-and-palestine extract, clipped by
#: notebooks/build_bats_vector_basemap.py.
OSM_BASEMAP_PATH = source.ASSET_DIR / "basemap_osm.geojson.gz"
OSM_ATTRIBUTION = ("Basemap features © OpenStreetMap contributors (ODbL), "
                   "Geofabrik israel-and-palestine extract")
OSM_WATER_FILL, OSM_WATER_EDGE = "#d9e9ed", "#a3bcc3"
OSM_BUILDING_FILL = "#dcdcdc"
OSM_ROAD_COLOR = "#cfcdca"
OSM_ROAD_WIDTHS = {"primary": .6, "secondary": .55, "secondary_link": .5,
                   "tertiary": .5, "tertiary_link": .45,
                   "unclassified": .32, "residential": .32, "road": .32}
OUTPUT_DIR = source.ROOT / "assets" / "goldshtein-bats-layout-review"
STEM = "bat_composite"
ONSET_INK = "#262626"
DIRECTED_FILL = "#87a1af"
UNDIRECTED_FILL = "#d4d5d5"
DIRECTEDNESS_CMAP = mpl.colors.ListedColormap(
    np.loadtxt(source.ASSET_DIR / "vanimo.txt"), name="vanimo"
).reversed().with_extremes(bad="white")
MAP_DROPOUT_COLOR = "#959b9c"
MAP_DROPOUT_STYLE = (0, (2, 2))
ONSET_MARKERS = {
    "bout_onset": (ONSET_INK, "white", 2.9, "Directed bout onset and 90% interval"),
    "flight_start": ("white", ONSET_INK, 2.6, "Directed at first usable step of flight"),
}


def first_directed_bout(scores, eligible, *, run=5):
    """Find the first qualifying run's starting index in each joint draw."""
    starts = np.full(scores.shape[1], -1, dtype=int)
    lengths = np.zeros_like(starts)
    for i, usable in enumerate(eligible):
        lengths = np.where(usable & (scores[i] > .5), lengths + 1, 0)
        first = (starts < 0) & (lengths == run)
        starts[first] = i - run + 1
    return starts


def onset_kind(start, usable, gaps):
    """Flag bouts detected at the flight's first usable step."""
    if start == np.argmax(usable) and not np.any(gaps[:start]):
        return "flight_start"
    return "bout_onset"


def control_heading_evidence(posterior, levels, posthoc, bearing_edges):
    """Average baseline-adjusted evidence over the control population draws."""
    mu, kappa, _ = _population_parameters(posterior, levels, posthoc)
    centers = ((bearing_edges[:-1] + bearing_edges[1:]) / 2)[:, None]
    return expit(posthoc.evidence_log_ratio(centers, mu[:, 0], kappa[:, 0])).mean(axis=1)


def flight_directedness(ctx, posthoc, levels):
    """Summarize the same recorded steps in routes, heatmaps, and treatment curves."""
    observed = ctx["df_track"]
    tracks = ctx["track_w_treat"]
    posterior = ctx["idata_vmu"].posterior
    # The observation-level deterministics are constant within each tag.
    tags, first = np.unique(observed["tag"], return_index=True)
    tag_columns = {tag: i for i, tag in enumerate(tags)}
    mu = posthoc._stack_samples_first(posterior["mu"].isel(mu_dim_0=first))
    kappa = posthoc._stack_samples_first(posterior["kappa"].isel(kappa_dim_0=first))
    headings = observed["HeadingRadians"].reindex(tracks.index)
    route_scores = headings.rename("directedness")
    lengths = observed.groupby("tag", observed=True).size()
    heatmaps, onsets, population = {}, {}, {}
    n_steps = int(lengths.max())
    step_km = source.STEP_LENGTH_M / 1000
    x_mid = np.arange(1 - n_steps, 1) * source.STEP_LENGTH_M
    for level in levels:
        flights = tracks.loc[tracks["treatment"].eq(level)].groupby("tag", observed=True)
        rows = sorted(flights, key=lambda item: (-lengths.loc[item[0]], str(item[0])))
        heatmap = np.full((len(rows), n_steps), np.nan)
        heatmaps[level] = heatmap
        onsets[level] = []
        sums = np.zeros((n_steps, mu.shape[0]))
        counts = np.zeros(n_steps, dtype=int)
        for row, (tag, flight) in enumerate(rows):
            theta = headings.loc[flight.index].to_numpy()
            valid = np.isfinite(theta)
            gaps = flight["gap_before"].to_numpy() == 1
            usable = valid & ~gaps
            i = tag_columns[tag]
            scores = expit(posthoc.evidence_log_ratio(theta[:, None], mu[:, i], kappa[:, i]))
            smoothed = posthoc.smooth_directedness(scores, valid, gaps, source.RECENT_WINDOW)
            means = np.full(len(theta), np.nan)
            means[usable] = smoothed[usable].mean(axis=1)
            route_scores.loc[flight.index] = means
            positions = np.flatnonzero(valid)
            offset = n_steps - len(positions)
            heatmap[row, offset:] = means[positions]
            bins = offset + np.flatnonzero(usable[positions])
            sums[bins] += smoothed[usable]
            counts[bins] += 1

            # Match the model card's first sustained goalward bout rule.
            starts = first_directed_bout(smoothed, usable & (np.abs(theta) < np.pi / 2))
            qualifying = starts[starts >= 0]
            if not len(qualifying):
                continue
            indices = np.quantile(qualifying, (.05, .5, .95), method="inverted_cdf")
            early, median, late = (len(positions) - 1 - np.searchsorted(positions, indices)) * step_km
            kind = onset_kind(int(indices[1]), usable, gaps)
            onsets[level].append((row, early, median, late, kind))
        average = np.full_like(sums, np.nan)
        supported = counts >= 2
        average[supported] = sums[supported] / counts[supported, None]
        lo, hi = posthoc.hdi_along_samples(average.T, prob=.90)
        population[level] = {"x_mid": x_mid, "mean": average.mean(axis=1), "lo": lo, "hi": hi,
                             "n_tags_by_x": counts}
    return route_scores, heatmaps, onsets, population


def _route_shades(color, n_routes):
    """Alternate dark and light tones within one treatment's color."""
    rgb = np.asarray(mpl.colors.to_rgb(color))
    ramp = mpl.colors.LinearSegmentedColormap.from_list(
        'route_shades', (rgb * .48, rgb, rgb + .35 * (1 - rgb))
    )
    order = np.empty(n_routes, dtype=int)
    order[::2] = np.arange((n_routes + 1) // 2)
    order[1::2] = np.arange(n_routes - 1, (n_routes - 1) // 2, -1)
    return ramp(np.linspace(0, 1, n_routes)[order])


def draw_flight_heatmap(ax, scores):
    """Plot one row per individual with fine white row separators."""
    n_bats, n_steps = scores.shape
    step_km = source.STEP_LENGTH_M / 1000
    far, end = (n_steps - .5) * step_km, -.5 * step_km
    ax.imshow(scores, extent=(far, end, n_bats, 0), aspect='auto',
              interpolation='nearest', cmap=DIRECTEDNESS_CMAP, norm=HEATMAP_NORM)
    ax.hlines(np.arange(1, n_bats), end, far, colors='white', linewidth=.42)
    ax.set_xlim(far, end)


def draw_flight_colorbar(fig, cax):
    """Label the direction of evidence on either side of the neutral midpoint."""
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=HEATMAP_NORM, cmap=DIRECTEDNESS_CMAP),
                      cax=cax, orientation='horizontal')
    cb.solids.set_rasterized(False)
    cb.solids.set_edgecolor('face')
    cb.outline.set_linewidth(.5)
    cb.set_ticks((0, .5, 1), labels=('0\nUndirected', '0.5\nEqual evidence', '1\nDirected'))
    cb.ax.tick_params(labelsize=7.3, length=2.25, pad=2)
    cb.ax.set_title('Posterior directedness', fontsize=9, pad=3)


def _draw_treatment_icons(ax, row_centers, colors):
    """Draw the four senses in each treatment row, fading the manipulated senses."""
    # Centers and disk radius in source figure points, matching the original grid.
    column_centers = (8.858, 28.272, 47.342, 66.656)
    disk_radius = 8.88

    def path(start, segments, *, fill="none", width=.9, tilt=0):
        vertices = [start]
        codes = [MplPath.MOVETO]
        for segment in segments:
            if len(segment) == 2:
                vertices.append(segment)
                codes.append(MplPath.LINETO)
            else:
                cx, cy, ex, ey = segment
                vertices.extend(((cx, cy), (ex, ey)))
                codes.extend((MplPath.CURVE3, MplPath.CURVE3))
        patch = PathPatch(MplPath(vertices, codes), facecolor=fill, edgecolor=INK,
                          linewidth=width, joinstyle="round", capstyle="round",
                          transform=Affine2D().rotate_deg(tilt) + transform, alpha=alpha)
        ax.add_patch(patch)

    for row, (y, color) in enumerate(zip(row_centers, colors, strict=True)):
        for column, x in enumerate(column_centers):
            muted = column < row
            alpha = .30 if muted else 1
            ax.add_patch(Circle((x, y), disk_radius, facecolor=color, edgecolor="none",
                                alpha=.15 if muted else .76, clip_on=False))
            transform = Affine2D().translate(x, y) + ax.transData
            if column == 0:
                # Almond-shaped eye with a round pupil.
                path((-6.7, 0), [(0, -7.6, 6.7, 0), (0, 7.6, -6.7, 0)], width=1.15)
                ax.add_patch(Circle((x, y), 2.4, facecolor=INK, edgecolor="none", alpha=alpha))
            elif column == 1:
                # The two squared pole tips keep the horseshoe recognizable at print size.
                path((-4.7, -5.3), [(-4.7, 1.2), (-4.7, 6, 0, 6),
                     (4.7, 6, 4.7, 1.2), (4.7, -5.3), (1.9, -5.3), (1.9, 1.2),
                     (1.9, 3.25, 0, 3.25), (-1.9, 3.25, -1.9, 1.2), (-1.9, -5.3),
                     (-4.7, -5.3)], fill=INK, width=.75, tilt=23)
                for left in (-4.7, 1.9):
                    path((left, -5.3), [(left + 2.8, -5.3), (left + 2.8, -3.0), (left, -3.0),
                         (left, -5.3)], fill="white", width=.75, tilt=23)
            elif column == 2:
                # Paired nostrils and nasal bridge.
                path((-1.8, -5.7), [(-2.0, -2.4, -3.5, -.5), (-6.2, 1.5, -4.8, 3.6),
                     (-4.0, 4.6, -2.6, 4.0)], width=1.05)
                path((1.8, -5.7), [(2.0, -2.4, 3.5, -.5), (6.2, 1.5, 4.8, 3.6),
                     (4.0, 4.6, 2.6, 4.0)], width=1.05)
                path((-3.2, 2.5), [(-1.7, 1.8, -.9, 3.6), (0, 5.2, .9, 3.6),
                     (1.7, 1.8, 3.2, 2.5)], width=.95)
            else:
                # Echo wavefronts arriving at the outer ear and its tragus.
                path((-4.6, -3.0), [(-1.8, -.9, -2.8, 2.2)], width=1.0)
                path((-6.7, -3.8), [(-3.3, -.8, -4.7, 3.3)], width=1.0)
                path((-1.25, .25), [(-1.7, -3.3, 1.1, -5.0), (2.5, -5.8, 4.0, -6.1),
                     (6.0, -3.7, 5.3, -.9), (5.1, 1.6, 3.35, 3.0),
                     (2.35, 3.8, 2.05, 5.0), (1.55, 6.7, .15, 6.0),
                     (-.9, 5.65, -.95, 4.6)], width=1.0)
                path((3.55, -3.7), [(1.20, -3.6, .55, -1.5), (.05, -.05, 1.45, .15),
                     (3.0, .3, 2.35, 1.85), (1.95, 2.6, 1.2, 2.95)], width=.95)


def draw_polar_density(ax, bearings, profile, color, radius_max, half_width_deg, *,
                       zero_radius=1.0, line_scale=1.0):
    """Draw the full mixture around a common zero-density reference ring."""
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set(ylim=(0, radius_max), xticks=[], yticks=[])
    ax.spines['polar'].set_visible(False)
    sector = np.abs(bearings) <= np.deg2rad(half_width_deg)
    ax.fill_between(bearings, 0, zero_radius + profile['mean'], where=sector,
                    color=color, alpha=.12, linewidth=0)
    ax.plot(bearings, np.full_like(bearings, zero_radius), color='#aaaaaa', lw=.55 * line_scale)
    upper_at_goal = profile['hi'][np.argmin(np.abs(bearings))]
    ax.plot([0, 0], [zero_radius, zero_radius + 1.08 * upper_at_goal], color='#9a9a9a', lw=.6 * line_scale)
    ax.fill_between(bearings, zero_radius + profile['lo'], zero_radius + profile['hi'],
                    color=color, alpha=.26, linewidth=0)
    ax.plot(bearings, zero_radius + profile['mean'], color=color, lw=1.05 * line_scale)


def draw_predictive_panel(data, colors, axes_at, page_text):
    """Compare goal-sector probabilities between groups and within each bat’s assigned group."""
    levels = data["levels"]
    bearings = data["density_bearings"]
    profiles = data["density_profiles"]
    examples = data["individual_density_profiles"]
    half_width = data["goal_window_half_width_degrees"]
    left, width = 120, 396
    column_width = width / len(levels)
    population_radius = 1 + 1.08 * max(p["hi"].max() for p in profiles.values())
    for i, (level, label) in enumerate(zip(levels, LABELS, strict=True)):
        center = left + (i + .5) * column_width
        page_text(center, 283, label, ha="center", color=colors[level], fontsize=10, weight="bold")
        polar = axes_at(center - 32, 292, 64, 64, projection="polar", facecolor="none")
        draw_polar_density(polar, bearings, profiles[level], colors[level], population_radius, half_width)

    page_text(18, 322, "Population density", fontsize=9.5, weight="bold")
    page_text(18, 363, "Treatment − control", fontsize=9.5, weight="bold")
    page_text(18, 376, "Sector probability", fontsize=9)
    contrast = axes_at(left, 344, width, 54)
    contrast.axhline(0, color="#777777", lw=.45, ls=(0, (3, 2)), zorder=0)
    contrast.plot(0, 0, "o", color=colors[levels[0]], ms=3.7)
    for i, level in enumerate(levels[1:], start=1):
        mean, lo, hi = data["goal_window_contrasts"][level]
        contrast.errorbar(i, mean, yerr=[[mean - lo], [hi - mean]], fmt="o", color=colors[level],
                          ms=4, elinewidth=1.05, capsize=2.3, capthick=.65)
    contrast.set(xlim=(-.5, len(levels) - .5), ylim=(-.206, .008), xticks=[], yticks=(-.2, -.1, 0))
    contrast.spines[["top", "right", "bottom"]].set_visible(False)
    contrast.yaxis.set_major_formatter(StrMethodFormatter("{x:g}"))
    contrast.tick_params(axis="y", labelsize=8.3, length=2.5, pad=3)

    # Rank bats only within their observed treatment. Horizontal position encodes
    # rank, while the vertical axis is the fitted probability for each bat.
    individual = axes_at(left, 408, width, 44)
    rows = data["goal_window_by_bat"]
    uniform = 2 * half_width / 360
    individual.set(xlim=contrast.get_xlim(), ylim=(0, 1), xticks=[], yticks=(0, .5, 1))
    individual.spines[["top", "right", "bottom"]].set_visible(False)
    individual.yaxis.set_major_formatter(StrMethodFormatter("{x:g}"))
    individual.tick_params(axis="y", labelsize=8.3, length=2.5, pad=3)
    individual.axhline(uniform, color="#494949", lw=.9, ls=(0, (3, 2)), zorder=1)
    page_text(18, 421, "Individual bats", fontsize=9.5, weight="bold")
    page_text(18, 434, "Sector probability", fontsize=9)
    for i, level in enumerate(levels):
        bats = sorted((r for r in rows if r["treatment"] == level), key=lambda r: (r["mean"], r["tag"]))
        count = len(bats)
        step = min(.055, .83 / max(count - 1, 1))
        x = i + (np.arange(count) - (count - 1) / 2) * step
        individual.vlines(x, [r["lo"] for r in bats], [r["hi"] for r in bats],
                          color=colors[level], alpha=.78, lw=.78, zorder=2)
        individual.scatter(x, [r["mean"] for r in bats], color=colors[level], s=7.4,
                           edgecolors="none", zorder=3)

    # The thumbnails share a radial scale. A smaller zero-density ring leaves
    # room for the individual lobes at this print size. Their centers fall
    # under the low, middle and high ranks with a gutter between treatments.
    page_text(18, 472, "Individual density", fontsize=9.5, weight="bold")
    zero_radius = .60
    largest_band = max(p["hi"].max() for group in examples.values() for p in group.values())
    radius_max = zero_radius + 1.08 * largest_band
    size = 2 * 15 * radius_max  # 15 print points per density unit in every thumbnail
    for i, level in enumerate(levels):
        group = examples[level]
        offsets = np.linspace(-.30, .30, len(group))
        for offset, profile in zip(offsets, group.values(), strict=True):
            center = left + (i + .5 + offset) * column_width
            polar = axes_at(center - size / 2, 477 - size / 2, size, size,
                            projection="polar", facecolor="none")
            draw_polar_density(polar, bearings, profile, colors[level], radius_max, half_width,
                               zero_radius=zero_radius, line_scale=2 / 3)


def draw_likelihood_diagrams(ax_vm, ax_vmu):
    """Illustrate two normalized likelihoods with similar goalward shapes."""
    bearings = np.linspace(-np.pi, np.pi, 361)
    cosine = np.cos(bearings)
    vm = np.exp(1.1 * cosine) / (2 * np.pi * np.i0(1.1))
    vmu_directed = .55 * np.exp(2.0 * cosine) / (2 * np.pi * np.i0(2.0))
    vmu_uniform = .45 / (2 * np.pi)
    reference = .24
    radius_max = reference + 1.06 * max(vm.max(), (vmu_directed + vmu_uniform).max())
    for ax, directed, flat in ((ax_vm, vm, 0), (ax_vmu, vmu_directed, vmu_uniform)):
        baseline = reference + flat
        top = baseline + directed
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set(ylim=(0, radius_max), xticks=[], yticks=[])
        ax.spines["polar"].set_visible(False)
        ax.fill_between(bearings, reference, baseline, color=UNDIRECTED_FILL, linewidth=0)
        ax.fill_between(bearings, baseline, top, color=DIRECTED_FILL, linewidth=0)
        ax.plot(bearings, top, color=INK, lw=.75)
        ax.plot(bearings, np.full_like(bearings, reference), color="#a0a0a0", lw=.45)


def center_zero_digit(label):
    """Align the zero in a centered degree label to its bearing tick."""
    measure = TextToPath().get_text_width_height_descent
    font = label.get_fontproperties()
    extra_width = measure("0°", font, ismath=False)[0] - measure("0", font, ismath=False)[0]
    shift = ScaledTranslation(extra_width / (2 * 72), 0, label.figure.dpi_scale_trans)
    label.set_transform(label.get_transform() + shift)


def draw_heading_evidence(ax, bearings, scores):
    """Use a fixed-width colored ring to map evidence around the circle."""
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_ylim(0, 1)
    ax.set_axis_off()
    ring = ax.pcolormesh(bearings, (.70, .96), scores[None, :], shading="flat",
                         cmap=DIRECTEDNESS_CMAP, norm=HEATMAP_NORM,
                         antialiased=False, linewidth=.18)
    ring.set_edgecolor("face")
    ticks = ((0, .48, "0°"), (90, .37, "90°"), (180, .43, "±180°"), (-90, .37, "−90°"))
    for angle, radius, label in ticks:
        tick = ax.text(np.deg2rad(angle), radius, label, ha="center", va="center", fontsize=6)
        if angle == 0:
            center_zero_digit(tick)


def draw_elpd_difference(plots, delta, ax, colors):
    """Use the notebook's paired ELPD differences inside the composite."""
    labels = dict(zip(colors, LABELS, strict=True))
    delta_labels = delta["treatment"].map(labels).fillna("Total")
    plots.plot_delta_total_and_by_treatment(
        delta.assign(label=delta_labels), treatment_order=list(colors),
        xlabel="ELPD difference (vMU − vM)",
        save=False, ax=ax, treatment_colors=colors,
    )
    for tick, color in zip(ax.get_yticklabels()[1:], colors.values(), strict=True):
        tick.set_color(color)
    ax.xaxis.set_major_locator(MultipleLocator(50))
    ax.xaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=5)
    ax.tick_params(axis="x", length=3, width=.65, labelsize=8.5, pad=2)


def draw_population_directedness(plots, ax, population, colors, distance_limits, distance_ticks):
    """Plot the full-range treatment means on the heatmaps' distance scale."""
    levels = list(colors)
    curves = {}
    for level, label in zip(levels, LABELS, strict=True):
        series = population[level]
        curves[level] = {**series, "x_mid": -np.asarray(series["x_mid"]) / 1000, "label": label}
    plots.plot_population_curve_smooth_compact(
        curves, levels, xlabel="Path distance remaining (km)", ylabel="Mean posterior\ndirectedness",
        xlim=distance_limits, axis_label_fs=8.7, tick_label_fs=7.8, y_labelpad=3,
        x_labelpad=3, save=False, ax=ax, treatment_colors=colors,
    )
    ax.set(xticks=distance_ticks, yticks=(0, .5, 1), yticklabels=("0", "0.5", "1"))
    ax.legend(loc="lower right", bbox_to_anchor=(.99, .04), frameon=False,
              ncols=2, fontsize=7.5, labelcolor="linecolor", handlelength=1.55, handletextpad=.42,
              columnspacing=1.0, borderaxespad=0, labelspacing=.25)


@lru_cache(maxsize=1)
def _vector_basemap():
    """The pinned OSM extract, projected (cached across maps)."""
    import geopandas as gpd
    features = json.load(gzip.open(OSM_BASEMAP_PATH, "rt"))["features"]
    return gpd.GeoDataFrame.from_features(
        features, crs="EPSG:4326").to_crs(f"EPSG:{source.METRIC_EPSG}")


def _ring_path(ring):
    xy = np.asarray(ring.xy).T
    codes = [MplPath.MOVETO] + [MplPath.LINETO] * (len(xy) - 2) + [MplPath.CLOSEPOLY]
    return MplPath(xy, codes)


def _polygon_path(polygon):
    """Compound path over a polygon's rings (holes oriented to punch out),
    flattened over MultiPolygon parts."""
    geoms = polygon.geoms if polygon.geom_type == "MultiPolygon" else [polygon]
    rings = []
    for geom in geoms:
        geom = orient(geom)
        rings.extend([geom.exterior, *geom.interiors])
    return MplPath.make_compound_path(*(_ring_path(ring) for ring in rings))


def draw_vector_basemap(ax, xlim, ylim):
    """Draw the pinned OpenStreetMap extract (water, roads, buildings) as a
    vector basemap — sharp at any zoom, no network.  Returns False when the
    extract is missing, leaving the caller to try raster tiles."""
    if not OSM_BASEMAP_PATH.exists():
        return False
    from shapely.geometry import box

    osm = _vector_basemap()
    extent = box(*xlim, *ylim)
    sel = osm[osm.intersects(extent)].clip(extent)

    for row in sel.itertuples():
        geom = row.geometry
        if geom is None or geom.is_empty:
            continue
        if row.kind == "water":
            ax.add_patch(PathPatch(_polygon_path(geom), facecolor=OSM_WATER_FILL,
                                   edgecolor=OSM_WATER_EDGE, lw=.35, zorder=1.0))

    road_segments = {}
    for row in sel.itertuples():
        if row.kind != "road":
            continue
        geom = row.geometry
        if geom is None or geom.is_empty:
            continue
        width = OSM_ROAD_WIDTHS.get(getattr(row, "fclass", None), .3)
        parts = geom.geoms if geom.geom_type == "MultiLineString" else [geom]
        for part in parts:
            if part.geom_type == "LineString" and len(part.coords) > 1:
                road_segments.setdefault(width, []).append(np.asarray(part.xy).T)
    for width, segments in road_segments.items():
        ax.add_collection(LineCollection(segments, colors=OSM_ROAD_COLOR,
                                         linewidths=width, zorder=1.4))

    rings = []
    for row in sel.itertuples():
        if row.kind != "building" or row.geometry is None or row.geometry.is_empty:
            continue
        geoms = (row.geometry.geoms
                 if row.geometry.geom_type == "MultiPolygon" else [row.geometry])
        for geom in geoms:
            geom = orient(geom)
            rings.extend([geom.exterior, *geom.interiors])
    if rings:
        ax.add_patch(PathPatch(MplPath.make_compound_path(*(_ring_path(r) for r in rings)),
                               facecolor=OSM_BUILDING_FILL, edgecolor="none", zorder=1.6))
    return True


def draw_track_map(ax, data, tracks, color, xlim, ylim, *, scores=None, show_scale=False,
                   basemap=True, zoom="auto"):
    """Show routes in individual treatment shades or on the directedness scale.

    ``basemap`` draws the landscape underneath.  The default source is the
    pinned OpenStreetMap extract (water/roads/buildings as true vectors,
    :func:`draw_vector_basemap` — sharp at any zoom, works offline).  When
    that extract is absent it falls back to live CartoDB "Positron" raster
    tiles (tiles fetched and reprojected onto the map's metric CRS; a
    ``CARTO_API_KEY`` in the environment, see
    :func:`goldshtein_bats_plots._carto_auth`, gets clean tiles rather than
    the watermarked anonymous pool); without either, the plain ``MAP_FACE``
    remains and the renderer says so.  ``zoom`` selects the fallback tile
    zoom (default: contextily's pick from the axis size, which overshoots
    for small panels and floods the tile service with requests)."""
    ax.set(xlim=xlim, ylim=ylim, xticks=[], yticks=[], facecolor=MAP_FACE, aspect="equal")
    for spine in ax.spines.values():
        spine.set_color(color)
        spine.set_linewidth(.65)
    if basemap and not draw_vector_basemap(ax, xlim, ylim):
        try:
            import contextily as cx
            from goldshtein_bats_plots import _carto_auth, _cartodb_quiet
            with _cartodb_quiet():
                cx.add_basemap(
                    ax, crs=f"EPSG:{source.METRIC_EPSG}", zorder=1, zoom=zoom,
                    source=_carto_auth(cx.providers.CartoDB.PositronNoLabels),
                    attribution=False)
        except Exception:  # offline renderer: keep the flat face
            # Never echo the exception text: tile failures quote the request
            # URL, which carries the CARTO key when one is configured.
            print("basemap tiles unavailable; drawing without tiles")
    lake = data["lake"]
    from goldshtein_bats_plots import LAKE_EDGE, LAKE_FILL  # data-root-safe local import
    polygons = [lake] if lake.geom_type == "Polygon" else lake.geoms
    for polygon in polygons:
        # Interior rings (holes/islands) are compound subpaths; orienting
        # the rings (exterior CCW, interiors CW) makes the fill rule punch
        # them out instead of painting over them.
        polygon = orient(polygon)
        rings = [polygon.exterior, *polygon.interiors]
        path = MplPath.make_compound_path(
            *(MplPath(np.asarray(ring.xy).T) for ring in rings))
        ax.add_patch(PathPatch(path, facecolor=LAKE_FILL, edgecolor=LAKE_EDGE,
                               lw=.45, zorder=2))

    missing, gaps, route, values, route_colors = [], [], [], [], []
    flights = tracks.groupby("tag", sort=False, observed=True)
    shades = (_route_shades(color, flights.ngroups) if scores is None
              else [None] * flights.ngroups)
    for (_, flight), shade in zip(flights, shades, strict=True):
        xy = flight[["x_m", "y_m"]].to_numpy()
        segments = np.stack((xy[:-1], xy[1:]), axis=1)
        crossed = flight["gap_before"].to_numpy()[1:] == 1
        gaps.extend(segments[crossed])
        if shade is not None:
            connected = segments[~crossed]
            route.extend(connected)
            route_colors.extend([shade] * len(connected))
            continue
        # The heading at each sample describes the segment that ends there.
        means = scores.loc[flight.index].to_numpy()[1:]
        has_score = np.isfinite(means) & ~crossed
        missing.extend(segments[~has_score & ~crossed])
        route.extend(segments[has_score])
        values.extend(means[has_score])

    ax.add_collection(LineCollection(missing, colors="#a7aeae", lw=.47, zorder=3))
    ax.add_collection(LineCollection(gaps, colors=MAP_DROPOUT_COLOR, lw=.45,
                                     linestyles=MAP_DROPOUT_STYLE, zorder=3))
    if scores is None:
        colored = LineCollection(route, colors=route_colors, lw=.67, capstyle="round", zorder=4)
        marker_style = {"mec": "black", "mfc": "white"}
    else:
        colored = LineCollection(route, cmap=DIRECTEDNESS_CMAP, norm=HEATMAP_NORM,
                                 lw=.88, capstyle="round", zorder=4)
        colored.set_array(np.asarray(values))
        marker_style = {"mec": "white", "mfc": INK}
    ax.add_collection(colored)
    for site in data["sites"]:
        ax.plot(*site, marker="X", ms=3.8, mew=.65, zorder=5, **marker_style)
    ax.plot(*data["roost"], marker="*", ms=5.2, mew=.65, zorder=5, **marker_style)
    if show_scale:
        x, y = xlim[1] - 460, ylim[0] + 590
        ax.plot([x - 1000, x], [y, y], color=INK, lw=1.1, zorder=6)
        ax.annotate("1 km", (x - 500, y), textcoords="offset points", xytext=(0, 3),
                    ha="center", fontsize=7)
        north_tip = (.90, .90)
        ax.annotate("", north_tip, (.90, .83), xycoords="axes fraction", zorder=6,
                    arrowprops={"arrowstyle": "-|>", "lw": .8, "color": INK,
                                "shrinkA": 0, "shrinkB": 0})
        ax.annotate("N", north_tip, xycoords="axes fraction", textcoords="offset points",
                    xytext=(0, 2), ha="center", fontsize=7, zorder=6)


def draw_onsets(ax, onsets):
    """Mark each bout start's median and conditional 90% interval."""
    for row, early, median, late, kind in onsets:
        y = row + .5
        ax.hlines(y, late, early, colors="white", linewidth=1.15, zorder=4)
        ax.hlines(y, late, early, colors=ONSET_INK, linewidth=.58, zorder=5)
        fill, edge, size, _ = ONSET_MARKERS[kind]
        ax.plot(median, y, marker="o", ms=size, mfc=fill, mec=edge, mew=.5, zorder=6)


def draw_composite(data, delta, flights, plots, *, photo_path=None):
    """Place the model panels above two matching grids of flight views."""
    route_scores, heatmaps, onsets, population = flights
    levels = data["levels"]
    colors = dict(zip(levels, PALETTE, strict=True))
    tracks = data["tracks"]
    treatments = {level: tracks.loc[tracks["treatment"].eq(level)] for level in levels}
    features = np.vstack([tracks[["x_m", "y_m"]].to_numpy(), data["sites"], data["roost"],
                          np.reshape(data["lake"].bounds, (2, 2))])
    xlim = features[:, 0].min() - 225, features[:, 0].max() + 225
    ylim = features[:, 1].min() - 300, features[:, 1].max() + 300
    map_ratio = np.ptp(xlim) / np.ptp(ylim)
    fig = plt.figure(figsize=(PAGE_W / 72, PAGE_H / 72), facecolor="white")

    def axes_at(x, y, width, height, **kwargs):
        return fig.add_axes((x / PAGE_W, 1 - (y + height) / PAGE_H,
                             width / PAGE_W, height / PAGE_H), **kwargs)

    def page_text(x, y, text, **kwargs):
        return fig.text(x / PAGE_W, 1 - y / PAGE_H, text, va="center", **kwargs)

    def letter(x, y, value, **kwargs):
        page_text(x, y, value + ".", fontsize=16, weight="bold", **kwargs)

    route_top, route_h = 31, 109
    route_center = route_top + route_h / 2
    photo_h = 102
    photo_top = route_center - photo_h / 2
    photo = axes_at(18, photo_top, 108, photo_h)
    with Image.open(photo_path or source.ASSET_DIR / "image_279.jpeg") as image:
        photo.imshow(image.convert("RGB"))
    photo.axis("off")
    letter(18, route_top - 9, "a")
    page_text(39, route_top - 9, "Homing experiment", weight="bold", fontsize=11)
    page_text(139, photo_top + 4, "Treatment groups", fontsize=10.7)
    icon_top = photo_top + 14
    # The original row centers and 0.74 artwork scale expressed in figure points.
    row_centers = (8.747, 32.538, 56.344, 80.135)
    icons = axes_at(139, icon_top, 75.4, 88.9)
    icons.set(xlim=(0, 75.4), ylim=(88.9, 0), aspect="equal")
    icons.axis("off")
    _draw_treatment_icons(icons, row_centers, colors.values())
    for level, label, y in zip(levels, LABELS, row_centers, strict=True):
        if label.startswith("+ "):
            page_text(223, icon_top + y, "+", ha="right", color=colors[level], fontsize=9.2)
        page_text(226, icon_top + y, label.removeprefix("+ "), color=colors[level], fontsize=9.2)

    route_w = route_h * map_ratio
    for center, level, label in zip((296, 358, 420, 482), levels, LABELS, strict=True):
        color = colors[level]
        page_text(center, route_top - 9, label, ha="center", color=color, fontsize=8.7)
        map_ax = axes_at(center - route_w / 2, route_top, route_w, route_h)
        draw_track_map(map_ax, data, treatments[level], color, xlim, ylim, show_scale=level == levels[-1],
                       zoom=TILE_ZOOM)
    map_key = plots._compact_destination_legend(
        show_dropout=tracks["gap_before"].eq(1).any(), gap_linestyle=MAP_DROPOUT_STYLE,
        gap_color=MAP_DROPOUT_COLOR,
    )
    fig.legend(handles=map_key, frameon=False, loc="center", ncols=len(map_key),
               bbox_to_anchor=(389 / PAGE_W, 1 - 150 / PAGE_H), fontsize=8,
               handlelength=1.45, handletextpad=.45, columnspacing=1.25, markerscale=.85)

    body_top = 170
    letter(18, body_top, "b")
    page_text(39, body_top, "Model selection", weight="bold", fontsize=11)
    vm_center, vmu_center = 48, 134
    model_top, model_size = 197, 35
    page_text(vm_center, 190, "von Mises (vM)", fontsize=8.7, ha="center")
    page_text(vmu_center, 190, "vM–Uniform mixture", fontsize=8.7, ha="center")
    draw_likelihood_diagrams(
        axes_at(vm_center - model_size / 2, model_top, model_size, model_size, projection="polar"),
        axes_at(vmu_center - model_size / 2, model_top, model_size, model_size, projection="polar"),
    )
    fig.legend(handles=[Patch(facecolor=DIRECTED_FILL, label="Directed"),
                        Patch(facecolor=UNDIRECTED_FILL, label="Undirected")],
               frameon=False, loc="center left", bbox_to_anchor=(177 / PAGE_W, 1 - 214 / PAGE_H),
               fontsize=8.6, handlelength=1.3, handletextpad=.45, borderaxespad=0)
    delta_ax = axes_at(320, 187, 193, 49)
    draw_elpd_difference(plots, delta, delta_ax, colors)

    letter(18, 268, "c")
    page_text(39, 268, "Goal-relative headings", weight="bold", fontsize=11)
    draw_predictive_panel(data, colors, axes_at, page_text)
    bearings = data["density_bearings"]

    flight_top = 504
    letter(18, flight_top, "d")
    page_text(39, flight_top, "Directedness along homing routes", weight="bold", fontsize=11)
    map_h = 160
    map_w = map_h * map_ratio
    map_x = 18
    map_lefts = (map_x, map_x + map_w + 14)
    grid_top = flight_top + 27
    map_tops = (grid_top, grid_top + map_h + 16)
    row_height = 4
    upper_rows = max(data["counts"][level] for level in levels[:2])
    lower_rows = max(data["counts"][level] for level in levels[2:])
    heatmap_top = grid_top + 111
    heatmap_tops = (heatmap_top, heatmap_top + upper_rows * row_height + 40)
    tile_x, tile_w, tile_gap = 187, 159, 11
    tile_lefts = (tile_x, tile_x + tile_w + tile_gap)
    right_width = 2 * tile_w + tile_gap
    step_km = source.STEP_LENGTH_M / 1000
    # Keep the half-step padding so the cells at both ends are fully visible.
    distance_limits = ((source.LAST_N_STEPS + .5) * step_km, -.5 * step_km)
    distance_ticks = (6, 3, 0)
    population_x = 221
    population_ax = axes_at(population_x, grid_top, tile_x + right_width - population_x, 69)
    draw_population_directedness(plots, population_ax, population, colors, distance_limits, distance_ticks)

    for i, (level, label) in enumerate(zip(levels, LABELS, strict=True)):
        row, col = divmod(i, 2)
        map_top = map_tops[row]
        left = map_lefts[col]
        tile_left = tile_lefts[col]
        color = colors[level]
        n_bats = data["counts"][level]
        page_text(left + map_w / 2, map_top - 8, label, ha="center", color=color,
                  fontsize=9.1, weight="bold")
        map_ax = axes_at(left, map_top, map_w, map_h)
        draw_track_map(map_ax, data, treatments[level], color, xlim, ylim,
                       scores=route_scores, show_scale=level == levels[-1],
                       zoom=TILE_ZOOM)

        tile_top = heatmap_tops[row]
        page_text(tile_left + tile_w / 2, tile_top - 8, label, ha="center", color=color,
                  fontsize=9.3, weight="bold")
        tile = axes_at(tile_left, tile_top, tile_w, n_bats * row_height)
        draw_flight_heatmap(tile, heatmaps[level])
        draw_onsets(tile, onsets[level])
        tile.set(xlim=distance_limits, ylim=(n_bats, 0), yticks=[], xticks=distance_ticks)
        tile.spines[["left", "top", "right"]].set_visible(False)
        tile.spines["bottom"].set_visible(True)
        if row == 1:
            tile.spines["bottom"].set_position(("outward", (lower_rows - n_bats) * row_height))
        tile.tick_params(axis="x", labelsize=8, width=.6, length=2, pad=2,
                         bottom=True, labelbottom=True)

    # Use the empty parts of the heatmap grid for its keys.
    evidence_polar = axes_at(tile_lefts[0] + 10, heatmap_top + 60, 45, 45, projection="polar")
    draw_heading_evidence(evidence_polar, bearings, data["heading_evidence"])
    key_center, key_w = tile_lefts[1] + tile_w / 2, 110
    key_top = heatmap_top + data["counts"][levels[1]] * row_height + 26
    draw_flight_colorbar(fig, axes_at(key_center - key_w / 2, key_top, key_w, 5.5))

    onset_key = [Line2D([], [], color=ONSET_INK, linewidth=1.1,
                        linestyle="-" if fill == ONSET_INK else "none",
                        marker="o", markersize=size + .8, markerfacecolor=fill,
                        markeredgecolor=edge, markeredgewidth=.5, label=label)
                 for fill, edge, size, label in ONSET_MARKERS.values()]
    fig.legend(handles=onset_key, frameon=False, loc="upper center",
               bbox_to_anchor=(key_center / PAGE_W, 1 - (key_top + 32) / PAGE_H),
               fontsize=7.6, handletextpad=.5, handlelength=1.65, borderaxespad=0)

    heatmap_bottom = heatmap_tops[1] + lower_rows * row_height
    tile_center = tile_x + right_width / 2
    page_text(tile_center, heatmap_bottom + 18, "Path distance remaining (km)", ha="center", fontsize=8.7)
    return fig


def save_figure(fig, output, stem):
    """Save a vector PDF and rasterize it for the review preview."""
    pdf = output / f"{stem}.pdf"
    fig.savefig(pdf, metadata={"Title": stem.replace("_", " "),
                               "CreationDate": None, "ModDate": None})
    with pymupdf.open(pdf) as document:
        document[0].get_pixmap(matrix=pymupdf.Matrix(2.7, 2.7), alpha=False).save(
            str(output / f"{stem}.png"))
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path,
                        default=source.ROOT / ".cache/goldshtein-bats-figure" / source.DATA_REVISION)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--font-dir", type=Path)
    parser.add_argument("--font-family", default="Arial")
    parser.add_argument("--photo", type=Path)
    args = parser.parse_args(argv)
    root = source.prepare_inputs(args.data_root, download=args.download)
    fonts = source.resolve_fonts(args.font_dir, args.font_family)
    ctx, posthoc, plots, lake_path = _read_context(root)
    data = assemble_layout_data(ctx, posthoc, plots, lake_path)
    data["heading_evidence"] = control_heading_evidence(
        ctx["idata_vmu"].posterior, data["levels"], posthoc, data["density_bearings"]
    )
    flights = flight_directedness(ctx, posthoc, data["levels"])
    args.output_dir.mkdir(parents=True, exist_ok=True)
    style = {"font.family": args.font_family, "font.size": 9, "pdf.fonttype": 42,
             "text.color": INK, "axes.labelcolor": INK,
             "axes.edgecolor": INK, "axes.linewidth": .65}
    with mpl.rc_context(style):
        fig = draw_composite(data, ctx["delta_df"], flights, plots, photo_path=args.photo)
        save_figure(fig, args.output_dir, STEM)
    summary = {
        "inputs_sha256": source.INPUT_HASHES,
        "data_revision": source.DATA_REVISION,
        "basemap": {
            "vector_extract": "notebooks/assets/goldshtein_bats_figure/basemap_osm.geojson.gz",
            "vector_extract_sha256": (source.sha256(OSM_BASEMAP_PATH)
                                      if OSM_BASEMAP_PATH.exists() else None),
            "attribution": OSM_ATTRIBUTION,
            "raster_fallback": ("live CartoDB Positron tiles via contextily when the "
                                "vector extract is absent (CARTO_API_KEY for clean tiles)"),
            "raster_tile_zoom": TILE_ZOOM,
        },
        "fonts": fonts,
        "n_bats_by_treatment": data["counts"],
        "n_model_steps": len(data["evidence"]),
        "elpd": data["scores"],
        "elpd_difference_vmu_minus_vm": {
            ("total" if row.row_type == "total" else str(int(row.treatment))): {
                "delta_elpd": float(row.delta_elpd), "se": float(row.se_delta_elpd)
            }
            for row in ctx["delta_df"].itertuples()
        },
        "panel_c_hdi_probability": PANEL_C_HDI_PROB,
        "goal_window_by_bat": data["goal_window_by_bat"],
        "individual_density_tags_by_treatment": {
            level: list(profiles) for level, profiles in data["individual_density_profiles"].items()
        },
        "density_at_zero_by_treatment": data["density_at_zero_by_treatment"],
        "density_at_zero_treatment_contrasts": data["full_contrasts"],
        "goal_window_half_width_degrees": data["goal_window_half_width_degrees"],
        "goal_window_rates": data["goal_window_groups"],
        "goal_window_treatment_contrasts": data["goal_window_contrasts"],
        "goal_window_pooled_contrast": data["pooled_goal_window_contrast"],
    }
    (args.output_dir / "bat_layout_review.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
