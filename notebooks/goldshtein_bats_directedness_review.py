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


PAGE_W, PAGE_H = 532, 663
OUTPUT_DIR = source.ROOT / "assets" / "goldshtein-bats-figure-ideas"
STEM = "bat_composite_directedness"
ONSET_GOLD = "#FFC850"
ONSET_EDGE = "#664512"
DIRECTED_FILL = "#87a1af"
UNDIRECTED_FILL = "#d4d5d5"
DIRECTEDNESS_CMAP = composite._directedness_colormap(
    (335, 204), lightness=.82, chroma=.135, neutral_lightness=.29
)
ONSET_MARKERS = {
    "bout_onset": ("o", ONSET_GOLD, 2.5, "Directed bout onset and 90% interval"),
    "flight_start": ("o", "white", 2.6, "Directed at first usable step of flight"),
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
    """Average control's baseline-adjusted evidence over posterior draws."""
    mu, kappa, _ = _population_parameters(posterior, levels, posthoc)
    centers = (bearing_edges[:-1] + bearing_edges[1:]) / 2
    log_ratio = posthoc.evidence_log_ratio(centers[:, None], mu[:, 0], kappa[:, 0])
    return expit(log_ratio).mean(axis=1)


def flight_directedness(ctx, posthoc, levels):
    """Build matched route and heatmap scores, and draw-wise bout onsets."""
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
    heatmaps, onsets = {}, {}
    n_steps = source.LAST_N_STEPS + 1
    step_km = source.STEP_LENGTH_M / 1000
    for level in levels:
        flights = tracks.loc[tracks["treatment"].eq(level)].groupby("tag", observed=True)
        rows = sorted(flights, key=lambda item: (-lengths.loc[item[0]], str(item[0])))
        heatmap = np.full((len(rows), n_steps), np.nan)
        heatmaps[level] = heatmap
        onsets[level] = []
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
            visible = means[valid][-n_steps:]
            heatmap[row, -len(visible):] = visible

            # Match the model card's first sustained goalward bout rule.
            starts = first_directed_bout(smoothed, usable & (np.abs(theta) < np.pi / 2))
            qualifying = starts[starts >= 0]
            if not len(qualifying):
                continue
            indices = np.quantile(qualifying, (.05, .5, .95), method="inverted_cdf")
            positions = np.flatnonzero(valid)
            early, median, late = (len(positions) - 1 - np.searchsorted(positions, indices)) * step_km
            kind = onset_kind(int(indices[1]), usable, gaps)
            onsets[level].append((row, early, median, late, kind))
    return route_scores, heatmaps, onsets


def draw_likelihood_diagrams(ax_vm, ax_vmu):
    """Illustrate two normalized likelihoods with a shared directed component."""
    bearings = np.linspace(-np.pi, np.pi, 361)
    kappa = 2.0
    vm = np.exp(kappa * np.cos(bearings)) / (2 * np.pi * np.i0(kappa))
    reference = .24
    for ax, directed_weight in ((ax_vm, 1.0), (ax_vmu, .55)):
        flat = (1 - directed_weight) / (2 * np.pi)
        baseline = reference + flat
        top = baseline + directed_weight * vm
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set(ylim=(0, reference + 1.06 * vm.max()), xticks=[], yticks=[])
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
    ring = ax.pcolormesh(bearings, (.66, .96), scores[None, :], shading="flat",
                         cmap=DIRECTEDNESS_CMAP, norm=composite.HEATMAP_NORM,
                         antialiased=False, linewidth=.18)
    ring.set_edgecolor("face")
    for angle, label in ((0, "0°"), (90, "90°"), (180, "180°"), (-90, "−90°")):
        tick = ax.text(np.deg2rad(angle), .45, label, ha="center", va="center", fontsize=6.7)
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
    ax.add_collection(LineCollection(gaps, colors="#959b9c", lw=.45,
                                     linestyles=(0, (2, 2)), zorder=3))
    if scores is None:
        colored = LineCollection(route, colors=route_colors, lw=.67, capstyle="round", zorder=4)
    else:
        colored = LineCollection(route, cmap=DIRECTEDNESS_CMAP, norm=composite.HEATMAP_NORM,
                                 lw=.88, capstyle="round", zorder=4)
        colored.set_array(np.asarray(values))
    ax.add_collection(colored)
    for site in data["sites"]:
        ax.plot(*site, marker="X", ms=3.8, mec="white", mew=.65, mfc=composite.INK, zorder=5)
    ax.plot(*data["roost"], marker="*", ms=5.2, mec="white", mew=.65, mfc=composite.INK, zorder=5)
    if show_scale:
        x, y = xlim[1] - 460, ylim[0] + 590
        ax.plot([x - 1000, x], [y, y], color=composite.INK, lw=1.1, zorder=6)
        ax.annotate("1 km", (x - 500, y), textcoords="offset points", xytext=(0, 3),
                    ha="center", fontsize=7)
        arrow_top = ylim[1] - 1580
        ax.annotate("", (x, arrow_top), (x, arrow_top - 720), zorder=6,
                    arrowprops={"arrowstyle": "-|>", "lw": .8, "color": composite.INK,
                                "shrinkA": 0, "shrinkB": 0})
        ax.annotate("N", (x, arrow_top), textcoords="offset points", xytext=(0, 2),
                    ha="center", fontsize=7)


def draw_onsets(ax, onsets):
    """Mark each bout start's median and conditional 90% interval."""
    for row, early, median, late, kind in onsets:
        y = row + .5
        ax.hlines(y, late, early, colors=ONSET_EDGE, linewidth=1.05, zorder=4)
        ax.hlines(y, late, early, colors=ONSET_GOLD, linewidth=.58, zorder=5)
        marker, fill, size, _ = ONSET_MARKERS[kind]
        ax.plot(median, y, marker=marker, ms=size, mfc=fill, mec=ONSET_EDGE, mew=.45, zorder=6)


def draw_composite(data, delta, flights, plots, *, photo_path=None):
    """Place the model panels above two matching grids of flight views."""
    route_scores, heatmaps, onsets = flights
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

    photo = axes_at(18, 18, 108, 102)
    with Image.open(photo_path or source.ASSET_DIR / "image_279.jpeg") as image:
        photo.imshow(image.convert("RGB"))
    photo.axis("off")
    if photo_path is None:
        photo.text(.5, .075, "Photo placeholder", transform=photo.transAxes,
                   ha="center", va="center", color="white", fontsize=8,
                   bbox={"facecolor": "#272727", "edgecolor": "none", "alpha": .58, "pad": 2.3})
    letter(26, 32, "a", color="white")
    page_text(139, 22, "Treatment groups", fontsize=10.7, weight="bold")
    ix0, iy0, ix1, iy1 = 202.9, 114.76, 304.79, 234.88
    icon_top, icon_scale = 32, .74
    icons = axes_at(139, icon_top, (ix1 - ix0) * icon_scale, (iy1 - iy0) * icon_scale)
    icons.set(xlim=(ix0, ix1), ylim=(source.PAGE_H - iy1, source.PAGE_H - iy0))
    icons.axis("off")
    icon_rows = (126.58, 158.73, 190.90, 223.05)
    drawings = source.panel_b_drawings(source.load_artwork())
    composite._draw_treatment_icons(icons, drawings, icon_rows, list(colors.values()))
    for level, label, y in zip(levels, composite.LABELS, icon_rows, strict=True):
        page_text(220, icon_top + (y - iy0) * icon_scale, label, color=colors[level], fontsize=9.2)

    route_h = 109
    route_w = route_h * map_ratio
    for center, level, label in zip((296, 358, 420, 482), levels, composite.LABELS, strict=True):
        color = colors[level]
        page_text(center, 22, label, ha="center", color=color, fontsize=8.7)
        map_ax = axes_at(center - route_w / 2, 31, route_w, route_h)
        draw_track_map(map_ax, data, treatments[level], color, xlim, ylim)

    letter(18, 151, "b")
    page_text(39, 151, "Model comparison", weight="bold", fontsize=11)
    page_text(63, 166, "von Mises (vM)", fontsize=8.5, ha="center")
    page_text(151, 166, "von Mises–Uniform\nMixture (vMU)", fontsize=8.3, ha="center")
    draw_likelihood_diagrams(axes_at(36, 177, 54, 54, projection="polar"),
                             axes_at(124, 177, 54, 54, projection="polar"))
    fig.legend(handles=[Patch(facecolor=DIRECTED_FILL, label="Directed"),
                        Patch(facecolor=UNDIRECTED_FILL, label="Undirected")],
               frameon=False, loc="center", bbox_to_anchor=(107 / PAGE_W, 1 - 233 / PAGE_H),
               ncols=2, fontsize=7.8, handlelength=1.3, handletextpad=.45, columnspacing=1.2)
    delta_ax = axes_at(69, 242, 122, 54)
    draw_elpd_difference(plots, delta, delta_ax, colors)

    letter(211, 151, "c")
    page_text(232, 151, "Predictive heading density", weight="bold", fontsize=11)
    bearings = data["density_bearings"]
    profiles = data["density_profiles"]
    radius_max = 1 + 1.08 * max(profile["hi"].max() for profile in profiles.values())
    polar_size = 66
    slots = ((250, 172), (318, 172), (250, 246), (318, 246))
    for (center, y), level in zip(slots, levels, strict=True):
        polar = axes_at(center - polar_size / 2, y, polar_size, polar_size, projection="polar")
        composite.draw_polar_density(polar, bearings, profiles[level], colors[level], radius_max)
        center_zero_digit(page_text(center, y - 3, "0°", ha="center", fontsize=8))

    x, y, width, height = 390, 200, 102, 73
    contrast = axes_at(x, y, width, height)
    composite.draw_density_contrast(contrast, data, colors)
    prop = contrast.xaxis.label.get_fontproperties()
    add_compact_unit_label(fig, x + width / 2, PAGE_H - (y + height + 23),
                           "Density difference vs. control", "radian", -1,
                           fontsize=prop.get_size_in_points(), prop=prop)

    letter(18, 329, "d")
    page_text(39, 329, "Directedness by individual", weight="bold", fontsize=11)
    map_h = 132
    map_w = map_h * map_ratio
    map_centers = (52, 127)
    tile_lefts = (194, 363)
    tops = (357, 507)
    row_height = 4
    lower_rows = max(data["counts"][level] for level in levels[2:])
    tile_w = 151
    for i, (level, label) in enumerate(zip(levels, composite.LABELS, strict=True)):
        row, col = divmod(i, 2)
        top = tops[row]
        map_center = map_centers[col]
        tile_left = tile_lefts[col]
        color = colors[level]
        n_bats = data["counts"][level]
        page_text(map_center, top - 8, label, ha="center", color=color, fontsize=9.1, weight="bold")
        map_ax = axes_at(map_center - map_w / 2, top, map_w, map_h)
        draw_track_map(map_ax, data, treatments[level], color, xlim, ylim,
                       scores=route_scores, show_scale=level == levels[-1])

        page_text(tile_left + tile_w / 2, top - 8, label, ha="center", color=color,
                  fontsize=9.3, weight="bold")
        tile = axes_at(tile_left, top, tile_w, n_bats * row_height)
        composite.draw_flight_heatmap(tile, heatmaps[level], cmap=DIRECTEDNESS_CMAP)
        draw_onsets(tile, onsets[level])
        tile.set(ylim=(n_bats, 0), yticks=[], xticks=(6, 3, 0))
        tile.spines[["left", "top", "right"]].set_visible(False)
        tile.spines["bottom"].set_visible(row == 1)
        if row == 1:
            tile.spines["bottom"].set_position(("outward", (lower_rows - n_bats) * row_height))
        tile.tick_params(axis="x", labelsize=8, width=.6, length=2, pad=2,
                         bottom=row == 1, labelbottom=row == 1)

    evidence_center = tile_lefts[1] + tile_w / 2
    evidence_polar = axes_at(evidence_center - 29, 421, 58, 58, projection="polar")
    draw_heading_evidence(evidence_polar, bearings, data["control_heading_evidence"])
    page_text(180, 470, "Individual", rotation=90, ha="center", fontsize=9.2)
    page_text(354, 588, "Path distance remaining (km)", ha="center", fontsize=9.5)
    composite.draw_flight_colorbar(fig, axes_at(201, 613, 122, 5.5), cmap=DIRECTEDNESS_CMAP)
    onset_key = [Line2D([], [], color=ONSET_GOLD, linewidth=1.1,
                        linestyle="-" if fill == ONSET_GOLD else "none",
                        marker=marker, markersize=size + .8, markerfacecolor=fill,
                        markeredgecolor=ONSET_EDGE, markeredgewidth=.5, label=label)
                 for marker, fill, size, label in ONSET_MARKERS.values()]
    fig.legend(handles=onset_key, frameon=False,
               loc="center left", bbox_to_anchor=(344 / PAGE_W, 1 - 619 / PAGE_H),
               fontsize=8.3, handletextpad=.55, handlelength=1.8, borderaxespad=0)
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
    data["control_heading_evidence"] = control_heading_evidence(
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
