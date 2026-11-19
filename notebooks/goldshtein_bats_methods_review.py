"""Render the compact-map composite and a fitted-parameter plot for review."""
from __future__ import annotations

import argparse
from pathlib import Path

import goldshtein_bats_figure_ideas as ideas
import goldshtein_bats_layout_review as composite
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import numpy as np
from scipy.stats import vonmises

from goldshtein_bats_layout_data import _read_context, assemble_layout_data


def heading_score_summary(ctx):
    """Bin the paired stepwise PSIS-LOO score difference by observed error."""
    scores = ctx["pointwise_elpd"]
    difference = scores["von Mises-Uniform"] - scores["von Mises"]
    error = np.rad2deg(np.abs(ctx["df_track"]["HeadingRadians"].to_numpy()))
    if difference.shape != error.shape or not np.isfinite(difference).all():
        raise ValueError("Expected one finite paired ELPD difference per modeled step")
    edges = np.arange(0, 181, 30)
    indices = np.minimum(np.searchsorted(edges, error, side="right") - 1, len(edges) - 2)
    n = len(edges) - 1
    counts = np.bincount(indices, minlength=n)
    sums = np.bincount(indices, weights=difference, minlength=n)
    if np.any(counts == 0):
        raise ValueError("Expected paired scores in all six heading-error bins")
    return {"centers": (edges[:-1] + edges[1:]) / 2, "means": sums / counts}


def mixture_parameter_summary(ctx, posthoc):
    """Summarize fitted component weight and conditional precision by bat."""
    rows = ctx["df_track"]
    _, first_row = np.unique(rows["tag"].to_numpy(), return_index=True)
    posterior = ctx["idata_vmu"].posterior
    kappa = posthoc._stack_samples_first(posterior["kappa"].isel(kappa_dim_0=first_row))
    weight = posthoc._stack_samples_first(
        posterior["p_goal"].isel(p_goal_dim_0=first_row, p_goal_dim_1=0)
    )
    # Center the tolerance on the directed component mean to isolate precision.
    precision = 2 * vonmises.cdf(np.pi / 6, kappa) - 1
    return {
        "treatment": rows["treatment"].to_numpy()[first_row].astype(int),
        "weight": weight.mean(axis=0),
        "precision": precision.mean(axis=0),
    }


def draw_mixture_parameters(parameters, levels):
    """Show two model quantities for each of the 67 bats."""
    width, height = 326, 228
    fig = plt.figure(figsize=(width / 72, height / 72), facecolor="white")
    fig.text(19 / width, 1 - 20 / height, "Fitted mixture parameters by bat",
             fontsize=12, weight="bold", va="center")
    ax = fig.add_axes([55 / width, (height - 191) / height, 252 / width, 131 / height])
    for level, label, color in zip(levels, composite.LABELS, composite.PALETTE, strict=True):
        selected = parameters["treatment"] == level
        ax.scatter(parameters["weight"][selected], parameters["precision"][selected],
                   color=color, edgecolors="white", linewidth=.45, s=23, alpha=.88,
                   label=label)
    fig.legend(*ax.get_legend_handles_labels(), loc="center",
               bbox_to_anchor=(.54, 1 - 41 / height), ncols=4, frameon=False,
               fontsize=8.5, handletextpad=.4, columnspacing=1.2,
               labelcolor=composite.PALETTE)
    ax.set(xlim=(.15, 1.02), ylim=(.2, .95), xticks=(.2, .4, .6, .8, 1),
           yticks=(.2, .4, .6, .8), xlabel="Directed-component weight",
           ylabel="Within ±30° of directed mean\n(probability)")
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_formatter(StrMethodFormatter("{x:g}"))
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=3, width=.65)
    return fig


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=composite.source.ROOT /
                        ".cache/goldshtein-bats-figure" / composite.source.DATA_REVISION)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=ideas.OUTPUT_DIR)
    parser.add_argument("--font-dir", type=Path)
    parser.add_argument("--font-family", default="Arial")
    args = parser.parse_args(argv)
    root = composite.source.prepare_inputs(args.data_root, download=args.download)
    ctx, posthoc, plots, lake_path = _read_context(root)
    layout = assemble_layout_data(ctx, posthoc, plots, lake_path)
    heatmaps = ideas.flight_directedness(ctx["df_track"], layout["levels"], posthoc)
    scores = heading_score_summary(ctx)
    parameters = mixture_parameter_summary(ctx, posthoc)
    composite.source.resolve_fonts(args.font_dir, args.font_family)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    style = {"font.family": args.font_family, "font.size": 9, "pdf.fonttype": 42,
             "text.color": composite.INK, "axes.labelcolor": composite.INK,
             "axes.edgecolor": composite.INK, "axes.linewidth": .65}
    with mpl.rc_context(style):
        fig = composite.draw_composite(layout, flight_heatmaps=heatmaps, heading_scores=scores)
        composite.save_figure(fig, args.output_dir, "bat_composite_methods")
        fig = draw_mixture_parameters(parameters, layout["levels"])
        composite.save_figure(fig, args.output_dir, "bat_mixture_parameters")


if __name__ == "__main__":
    main()
