"""Render two candidate homing panels from the pinned bat model and tracks."""
from __future__ import annotations

import argparse
from pathlib import Path

import goldshtein_bats_layout_review as composite
import goldshtein_bats_publication_figure as source
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Arc, Wedge
from matplotlib.ticker import MultipleLocator, StrMethodFormatter
import numpy as np

from goldshtein_bats_layout_data import _goal_window_mass, _population_parameters, _read_context


PAGE_W, PAGE_H = composite.PAGE_W, 317
ABSENT = "#eeeeee"
OUTPUT_DIR = source.ROOT / "assets" / "goldshtein-bats-figure-ideas"


def goal_cone_curves(posterior, levels, posthoc):
    """Summarize heading mass across cone widths and paired control contrasts."""
    half_angle = np.arange(0.0, 181.0, 3.0)
    tolerance = np.deg2rad(half_angle)
    mu, kappa, weight = _population_parameters(posterior, levels, posthoc)
    draws = [
        _goal_window_mass(tolerance, mu[:, i, None], kappa[:, i, None],
                          weight[:, i, None], posthoc.KAPPA_UNIFORM)
        for i in range(len(levels))
    ]
    means = {level: group.mean(axis=0) for level, group in zip(levels, draws, strict=True)}
    contrasts = {}
    for level, group in zip(levels[1:], draws[1:], strict=True):
        difference = group - draws[0]
        lo, hi = posthoc.hdi_along_samples(difference, prob=0.95)
        contrasts[level] = (difference.mean(axis=0), lo, hi)
    return half_angle, means, contrasts


def flight_directedness(tracks, levels, posthoc):
    """Arrange smoothed step scores by treatment, aligned on the last step."""
    n_steps = tracks.groupby("tag", observed=True).size().max()
    heatmaps = {}
    for level in levels:
        flights = tracks[tracks["treatment"].eq(level)].groupby("tag", observed=True)
        rows = sorted(flights, key=lambda item: (-len(item[1]), str(item[0])))
        heatmap = np.full((len(rows), n_steps), np.nan)
        for row, (_, flight) in enumerate(rows):
            scores = flight["evi_post_mean"].to_numpy()
            gaps = flight["gap_before"].to_numpy() == 1
            # The smoother is linear and uses the same weights for every draw.
            heatmap[row, -len(scores):] = posthoc.smooth_directedness(
                scores[:, None], np.isfinite(scores), gaps, source.RECENT_WINDOW
            )[:, 0]
        heatmaps[level] = heatmap
    return heatmaps


def draw_ideas(cones, heatmaps, levels):
    """Draw the goal-cone and flight-by-flight views at composite width."""
    fig = plt.figure(figsize=(PAGE_W / 72, PAGE_H / 72), facecolor="white")
    colors = dict(zip(levels, composite.PALETTE, strict=True))

    def axes_at(x, y, w, h):
        return fig.add_axes([x / PAGE_W, 1 - (y + h) / PAGE_H, w / PAGE_W, h / PAGE_H])

    def page_text(x, y, text, **kwargs):
        return fig.text(x / PAGE_W, 1 - y / PAGE_H, text, va="center", **kwargs)

    page_text(23, 24, "A", fontsize=15, weight="bold")
    page_text(43, 24, "Heading within the goal cone", fontsize=12, weight="bold")
    page_text(336, 24, "B", fontsize=15, weight="bold")
    page_text(356, 24, "Homing, flight by flight", fontsize=12, weight="bold")

    cone = axes_at(51, 38, 87, 45)
    cone.set(xlim=(-0.07, 2.35), ylim=(-0.69, 0.71), aspect="equal")
    cone.axis("off")
    example_angle = 36
    cone.add_patch(Wedge((0, 0), 1, -example_angle, example_angle,
                         facecolor="#e3edf2", edgecolor="#7894a3", linewidth=.65))
    cone.annotate("", (1.09, 0), (0, 0), arrowprops={"arrowstyle": "->", "lw": .85})
    cone.add_patch(Arc((0, 0), 1.06, 1.06, theta1=0, theta2=example_angle,
                       color=composite.INK, linewidth=.9))
    angle_mid = np.deg2rad(example_angle / 2)
    cone.text(.72 * np.cos(angle_mid), .72 * np.sin(angle_mid), "α",
              ha="center", va="center", fontsize=10)
    cone.text(1.17, 0, "Goal bearing", va="center", fontsize=8)

    fig.legend([Line2D([], [], color=colors[level], lw=1.5) for level in levels],
               composite.LABELS, loc="center", bbox_to_anchor=(216 / PAGE_W, 1 - 57 / PAGE_H),
               ncols=2, frameon=False, columnspacing=1.3, handlelength=1.45, fontsize=8.5)

    angles, means, contrasts = cones
    share = axes_at(58, 91, 248, 85)
    delta = axes_at(58, 197, 248, 82)
    for ax in (share, delta):
        ax.spines[["top", "right"]].set_visible(False)
        ax.set(xlim=(0, 180), xticks=(0, 45, 90, 135, 180))
        ax.tick_params(length=3, width=.65)

    for level, curve in means.items():
        share.plot(angles, curve, color=colors[level], lw=1.35)
    share.set(ylim=(0, 1.02), yticks=(0, .25, .5, .75, 1), ylabel="Predictive probability")
    share.yaxis.set_major_formatter(StrMethodFormatter("{x:g}"))
    share.tick_params(axis="x", labelbottom=False)

    delta.axhline(0, color="#666666", lw=.65, ls=(0, (3, 2)), zorder=2)
    for level, (mean, lo, hi) in contrasts.items():
        delta.fill_between(angles, lo, hi, color=colors[level], alpha=.16, lw=0)
        delta.plot(angles, mean, color=colors[level], lw=1.35, zorder=3)
    delta.yaxis.set_major_locator(MultipleLocator(.05))
    delta.set(xlabel="Half-angle α from the goal bearing (°)",
              ylabel="Difference in probability\nvs. control")

    page_text(356, 43, "One row per bat. Gray marks distance beyond its modeled path.",
              fontsize=8.5, color="#555555")
    cmap = mpl.colormaps["viridis"].with_extremes(bad=ABSENT)
    norm = mpl.colors.Normalize(0, 1)
    step_km = source.STEP_LENGTH_M / 1000
    far, end = (heatmaps[levels[0]].shape[1] - .5) * step_km, -.5 * step_km
    top = 57
    for level, label in zip(levels, composite.LABELS, strict=True):
        heatmap = heatmaps[level]
        n_bats = len(heatmap)
        height = 2.77 * n_bats
        ax = axes_at(356, top, 202, height)
        ax.imshow(heatmap, extent=(far, end, n_bats, 0), aspect="auto",
                  interpolation="nearest", cmap=cmap, norm=norm)
        ax.set(yticks=[])
        ax.spines[:].set_visible(False)
        page_text(568, top + height / 2 - 4.4, label, color=colors[level],
                  fontsize=9.5, weight="bold")
        page_text(568, top + height / 2 + 6.4, f"{n_bats} bats", fontsize=8.5, color="#555555")
        if level == levels[-1]:
            ax.set(xticks=np.arange(0, far, 2), xlabel="Modeled path remaining (km)")
            ax.tick_params(axis="x", length=2.6, width=.6, pad=2)
            ax.spines["bottom"].set_visible(True)
        else:
            ax.set(xticks=[])
        top += height + 5.2

    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=axes_at(386, 295, 112, 5.5),
                      orientation="horizontal", ticks=(0, .5, 1), format="%g")
    cb.outline.set_linewidth(.5)
    cb.ax.tick_params(labelsize=8, length=2.25, pad=2)
    page_text(506, 298, "Posterior directedness", fontsize=8.5)
    return fig


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path,
                        default=source.ROOT / ".cache/goldshtein-bats-figure" / source.DATA_REVISION)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--font-dir", type=Path)
    parser.add_argument("--font-family", default="Arial")
    args = parser.parse_args(argv)

    source.resolve_fonts(args.font_dir, args.font_family)
    root = source.prepare_inputs(args.data_root, download=args.download)
    ctx, posthoc, _, _ = _read_context(root)
    levels = [int(level) for level in ctx["treatment_order"]]
    cones = goal_cone_curves(ctx["idata_vmu"].posterior, levels, posthoc)
    heatmaps = flight_directedness(ctx["df_track"], levels, posthoc)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    style = {"font.family": args.font_family, "font.size": 9, "pdf.fonttype": 42,
             "text.color": composite.INK, "axes.labelcolor": composite.INK,
             "axes.edgecolor": composite.INK, "axes.linewidth": .65}
    with mpl.rc_context(style):
        fig = draw_ideas(cones, heatmaps, levels)
        composite.save_figure(fig, args.output_dir, "bat_figure_ideas")


if __name__ == "__main__":
    main()
