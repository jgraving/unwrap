"""Render a compact bat composite with paired maps and flight heatmaps."""
from __future__ import annotations

import argparse
from pathlib import Path

import goldshtein_bats_layout_review as composite
import goldshtein_bats_publication_figure as source
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.textpath import TextToPath
from matplotlib.ticker import MultipleLocator, StrMethodFormatter
from matplotlib.transforms import ScaledTranslation
import numpy as np
from PIL import Image
from scipy.special import expit

from goldshtein_bats_layout_data import _population_parameters, _read_context, assemble_layout_data
from locust_density_figure import add_compact_unit_label


PAGE_W, PAGE_H = 532, 736
OUTPUT_DIR = source.ROOT / "assets" / "goldshtein-bats-figure-ideas"
STEM = "bat_composite_directedness"
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
        ax.plot(bearings, top, color=composite.INK, lw=.75)
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
                         cmap=DIRECTEDNESS_CMAP, norm=composite.HEATMAP_NORM,
                         antialiased=False, linewidth=.18)
    ring.set_edgecolor("face")
    ticks = ((0, .48, "0°"), (90, .37, "90°"), (180, .43, "±180°"), (-90, .37, "−90°"))
    for angle, radius, label in ticks:
        tick = ax.text(np.deg2rad(angle), radius, label, ha="center", va="center", fontsize=6)
        if angle == 0:
            center_zero_digit(tick)


def draw_elpd_difference(plots, delta, ax, colors):
    """Use the notebook's paired ELPD differences inside the composite."""
    labels = dict(zip(colors, composite.LABELS, strict=True))
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
    for level, label in zip(levels, composite.LABELS, strict=True):
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


def draw_track_map(ax, data, tracks, color, xlim, ylim, *, scores=None, show_scale=False):
    """Show routes in individual treatment shades or on the directedness scale."""
    ax.set(xlim=xlim, ylim=ylim, xticks=[], yticks=[], facecolor=composite.MAP_FACE, aspect="equal")
    for spine in ax.spines.values():
        spine.set_color(color)
        spine.set_linewidth(.65)
    lake = data["lake"]
    polygons = [lake] if lake.geom_type == "Polygon" else lake.geoms
    for polygon in polygons:
        ax.fill(*polygon.exterior.xy, fc="#d9e9ed", ec="#a3bcc3", lw=.45, zorder=2)

    missing, gaps, route, values, route_colors = [], [], [], [], []
    flights = tracks.groupby("tag", sort=False, observed=True)
    shades = (composite._route_shades(color, flights.ngroups) if scores is None
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
        colored = LineCollection(route, cmap=DIRECTEDNESS_CMAP, norm=composite.HEATMAP_NORM,
                                 lw=.88, capstyle="round", zorder=4)
        colored.set_array(np.asarray(values))
        marker_style = {"mec": "white", "mfc": composite.INK}
    ax.add_collection(colored)
    for site in data["sites"]:
        ax.plot(*site, marker="X", ms=3.8, mew=.65, zorder=5, **marker_style)
    ax.plot(*data["roost"], marker="*", ms=5.2, mew=.65, zorder=5, **marker_style)
    if show_scale:
        x, y = xlim[1] - 460, ylim[0] + 590
        ax.plot([x - 1000, x], [y, y], color=composite.INK, lw=1.1, zorder=6)
        ax.annotate("1 km", (x - 500, y), textcoords="offset points", xytext=(0, 3),
                    ha="center", fontsize=7)
        north_tip = (.90, .90)
        ax.annotate("", north_tip, (.90, .83), xycoords="axes fraction", zorder=6,
                    arrowprops={"arrowstyle": "-|>", "lw": .8, "color": composite.INK,
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
    colors = dict(zip(levels, composite.PALETTE, strict=True))
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
    if photo_path is None:
        photo.text(.5, .075, "Photo placeholder", transform=photo.transAxes,
                   ha="center", va="center", color="white", fontsize=8,
                   bbox={"facecolor": "#272727", "edgecolor": "none", "alpha": .58, "pad": 2.3})
    letter(18, route_top - 9, "a")
    page_text(39, route_top - 9, "Homing experiment", weight="bold", fontsize=11)
    page_text(139, photo_top + 4, "Treatment groups", fontsize=10.7)
    ix0, iy0, ix1, iy1 = 202.9, 114.76, 304.79, 234.88
    icon_top, icon_scale = photo_top + 14, .74
    icons = axes_at(139, icon_top, (ix1 - ix0) * icon_scale, (iy1 - iy0) * icon_scale)
    icons.set(xlim=(ix0, ix1), ylim=(source.PAGE_H - iy1, source.PAGE_H - iy0))
    icons.axis("off")
    icon_rows = (126.58, 158.73, 190.90, 223.05)
    drawings = source.panel_b_drawings(source.load_artwork())
    composite._draw_treatment_icons(icons, drawings, icon_rows, list(colors.values()), active_alpha=.76)
    for level, label, y in zip(levels, composite.LABELS, icon_rows, strict=True):
        label_y = icon_top + (y - iy0) * icon_scale
        if label.startswith("+ "):
            page_text(223, label_y, "+", ha="right", color=colors[level], fontsize=9.2)
        page_text(226, label_y, label.removeprefix("+ "), color=colors[level], fontsize=9.2)

    route_w = route_h * map_ratio
    for center, level, label in zip((296, 358, 420, 482), levels, composite.LABELS, strict=True):
        color = colors[level]
        page_text(center, route_top - 9, label, ha="center", color=color, fontsize=8.7)
        map_ax = axes_at(center - route_w / 2, route_top, route_w, route_h)
        draw_track_map(map_ax, data, treatments[level], color, xlim, ylim, show_scale=level == levels[-1])
    map_key = plots._compact_destination_legend(
        show_dropout=tracks["gap_before"].eq(1).any(), gap_linestyle=MAP_DROPOUT_STYLE,
        gap_color=MAP_DROPOUT_COLOR,
    )
    fig.legend(handles=map_key, frameon=False, loc="center", ncols=len(map_key),
               bbox_to_anchor=(389 / PAGE_W, 1 - 150 / PAGE_H), fontsize=8,
               handlelength=1.45, handletextpad=.45, columnspacing=1.25, markerscale=.85)

    body_top = 168
    letter(18, body_top, "b")
    page_text(39, body_top, "Model selection", weight="bold", fontsize=11)
    page_text(63, body_top + 16, "von Mises (vM)", fontsize=8.5, ha="center")
    page_text(151, body_top + 16, "von Mises–Uniform\nMixture (vMU)", fontsize=8.3, ha="center")
    draw_likelihood_diagrams(axes_at(41.5, body_top + 28, 43, 43, projection="polar"),
                             axes_at(129.5, body_top + 28, 43, 43, projection="polar"))
    fig.legend(handles=[Patch(facecolor=DIRECTED_FILL, label="Directed"),
                        Patch(facecolor=UNDIRECTED_FILL, label="Undirected")],
               frameon=False, loc="center", bbox_to_anchor=(107 / PAGE_W, 1 - (body_top + 69) / PAGE_H),
               ncols=2, fontsize=7.8, handlelength=1.3, handletextpad=.45, columnspacing=1.2)
    delta_ax = axes_at(69, body_top + 80, 126, 51)
    draw_elpd_difference(plots, delta, delta_ax, colors)

    letter(211, body_top, "c")
    page_text(232, body_top, "Predictive heading density", weight="bold", fontsize=11)
    bearings = data["density_bearings"]
    profiles = data["density_profiles"]
    radius_max = 1 + 1.08 * max(profile["hi"].max() for profile in profiles.values())
    contrast_x, contrast_w = 250, 266
    column_w = contrast_w / len(levels)
    polar_top, polar_size = body_top + 18, 64
    for i, level in enumerate(levels):
        center = contrast_x + (i + .5) * column_w
        polar = axes_at(center - polar_size / 2, polar_top, polar_size, polar_size,
                        projection="polar")
        composite.draw_polar_density(polar, bearings, profiles[level], colors[level], radius_max)
        center_zero_digit(page_text(center, polar_top - 3, "0°", ha="center", fontsize=8))

    contrast_top, contrast_h = body_top + 80, 63
    contrast = axes_at(contrast_x, contrast_top, contrast_w, contrast_h)
    composite.draw_density_contrast(contrast, data, colors, vertical=True)
    # The open area is bounded by the zero line and the first treatment interval.
    lo, hi = contrast.get_ylim()
    zero_y = contrast_top + contrast_h * hi / (hi - lo)
    open_right = contrast_x + 1.5 * column_w - 3
    label_x = (contrast_x + open_right) / 2 - 5
    label_center_y = (zero_y + contrast_top + contrast_h) / 2
    label = page_text(label_x, label_center_y - 4, "Density difference", ha="center", fontsize=8.2)
    prop = label.get_fontproperties()
    add_compact_unit_label(fig, label_x, PAGE_H - (label_center_y + 7), "vs. control", "rad", -1,
                           fontsize=prop.get_size_in_points(), prop=prop)

    flight_top = 350
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
    heatmap_top = grid_top + 119
    heatmap_tops = (heatmap_top, heatmap_top + upper_rows * row_height + 32)
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

    for i, (level, label) in enumerate(zip(levels, composite.LABELS, strict=True)):
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
                       scores=route_scores, show_scale=level == levels[-1])

        tile_top = heatmap_tops[row]
        if row == 0:
            # Align the upper axes below the heatmaps, with room for the keys above Eye cover.
            tile_top += (upper_rows - n_bats) * row_height
        page_text(tile_left + tile_w / 2, tile_top - 8, label, ha="center", color=color,
                  fontsize=9.3, weight="bold")
        tile = axes_at(tile_left, tile_top, tile_w, n_bats * row_height)
        composite.draw_flight_heatmap(tile, heatmaps[level], cmap=DIRECTEDNESS_CMAP)
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
    key_top = heatmap_top - 13
    composite.draw_flight_colorbar(fig, axes_at(key_center - key_w / 2, key_top, key_w, 5.5),
                                   cmap=DIRECTEDNESS_CMAP, tick_fontsize=7.3)

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
    source.resolve_fonts(args.font_dir, args.font_family)
    ctx, posthoc, plots, lake_path = _read_context(root)
    data = assemble_layout_data(ctx, posthoc, plots, lake_path)
    data["heading_evidence"] = control_heading_evidence(
        ctx["idata_vmu"].posterior, data["levels"], posthoc, data["density_bearings"]
    )
    flights = flight_directedness(ctx, posthoc, data["levels"])
    args.output_dir.mkdir(parents=True, exist_ok=True)
    style = {"font.family": args.font_family, "font.size": 9, "pdf.fonttype": 42,
             "text.color": composite.INK, "axes.labelcolor": composite.INK,
             "axes.edgecolor": composite.INK, "axes.linewidth": .65}
    with mpl.rc_context(style):
        fig = draw_composite(data, ctx["delta_df"], flights, plots, photo_path=args.photo)
        composite.save_figure(fig, args.output_dir, STEM)


if __name__ == "__main__":
    main()
