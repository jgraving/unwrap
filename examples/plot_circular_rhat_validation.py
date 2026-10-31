"""Plot a compact circular-diagnostic summary from saved validation JSON.

Run from the repository root; no sampling or model dependencies are needed::

    python examples/plot_circular_rhat_validation.py

Use ``--preview /path/to/preview.png`` to also render a PNG for inspection.
"""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


STYLES = {
    "raw_rank_rhat": ("Raw rank", "#667085", "s"),
    "joint_circular": ("Joint circular", "#007F73", "o"),
    "rank_sin": ("Rank sine", "#D55E00", "^"),
    "rank_cos": ("Rank cosine", "#D55E00", "v"),
}
ANGULAR_METRICS = ("joint_circular", "rank_sin", "rank_cos")


def style_axis(axis):
    axis.spines[["top", "right", "left"]].set_visible(False)
    axis.spines["bottom"].set_color("#D0D5DD")
    axis.tick_params(axis="y", length=0, pad=8)
    axis.tick_params(axis="x", length=3, color="#D0D5DD", pad=6)
    axis.set_axisbelow(True)
    axis.grid(axis="x", color="#EAECF0", linewidth=0.7)


def points(axis, metrics, positions, values, text_offset):
    for metric, position in zip(metrics, positions):
        _, color, marker = STYLES[metric]
        value = values[metric]
        axis.plot(
            value, position, marker=marker, color=color, markersize=7, linestyle="none"
        )
        axis.text(
            value + text_offset, position, f"{value:.4f}", va="center", fontsize=9.5
        )


def make_figure(result):
    """Keep fitted-posterior values and synthetic medians explicitly separate."""
    fit = next(
        (
            row
            for row in result["mcmc"]
            if row["scenario"] == "period_separated_intercepts"
        ),
        None,
    )
    if fit is None:
        raise ValueError(
            "The summary figure requires saved period-separated PyMC fit results."
        )
    draws = max(row["draws_per_chain"] for row in result["summaries"])
    summaries = {
        row["scenario"]: row["metrics"]
        for row in result["summaries"]
        if row["draws_per_chain"] == draws
    }

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9.5,
            "text.color": "#243244",
            "axes.labelcolor": "#243244",
            "xtick.color": "#667085",
            "ytick.color": "#243244",
            "svg.fonttype": "none",
            "svg.hashsalt": "unwrap-circular-rhat-summary",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )
    figure, (fitted, synthetic) = plt.subplots(1, 2, figsize=(8.8, 5.0))
    figure.subplots_adjust(left=0.14, right=0.98, bottom=0.29, top=0.69, wspace=0.60)
    figure.text(
        0.04,
        0.95,
        "Full turns preserve direction; multimodal differences can be missed",
        fontsize=12.5,
        weight="bold",
    )
    figure.text(
        0.04,
        0.89,
        "Saved validation results · joint covariance and rank statistics remain distinct",
        fontsize=9.5,
        color="#667085",
    )

    fitted.set_title("A  Different full-turn modes", loc="left", fontsize=10.5, pad=38)
    fitted.text(
        0,
        1.10,
        f"One PyMC fit · {fit['chains']} chains × {fit['draws_per_chain']:,} retained draws",
        transform=fitted.transAxes,
        fontsize=9,
        color="#667085",
    )
    fitted_metrics = tuple(STYLES)
    fitted_positions = (3, 2, 1, 0)
    points(fitted, fitted_metrics, fitted_positions, fit, text_offset=0.07)
    fitted.set_yticks(
        fitted_positions, [STYLES[metric][0] for metric in fitted_metrics]
    )
    fitted.set_xlim(0.9, 3.42)
    fitted.set_ylim(-0.55, 3.55)
    fitted.set_xticks([1.0, 1.5, 2.0, 2.5, 3.0])
    fitted.set_xlabel("Diagnostic value", labelpad=9)
    style_axis(fitted)

    synthetic.set_title("B  Multimodal blind spots", loc="left", fontsize=10.5, pad=38)
    synthetic.text(
        0,
        1.10,
        f"Medians · {result['replicates']} replicates · {draws:,} draws/chain",
        transform=synthetic.transAxes,
        fontsize=9,
        color="#667085",
    )
    scenarios = (
        ("different_antipodal_pairs", "Different antipodal pairs", (5.6, 4.6, 3.6)),
        (
            "same_projection_marginals",
            "Matching projection marginals",
            (1.6, 0.6, -0.4),
        ),
    )
    for scenario, label, positions in scenarios:
        values = {
            metric: summaries[scenario][metric]["median"] for metric in ANGULAR_METRICS
        }
        points(synthetic, ANGULAR_METRICS, positions, values, text_offset=0.03)
        synthetic.text(
            0.0,
            positions[0] + 0.8,
            label,
            transform=synthetic.get_yaxis_transform(),
            fontsize=9,
            weight="bold",
        )
    positions = [position for _, _, group in scenarios for position in group]
    synthetic.set_yticks(
        positions, [STYLES[metric][0] for metric in ANGULAR_METRICS] * 2
    )
    synthetic.set_xlim(0.95, 1.96)
    synthetic.set_ylim(-0.95, 6.95)
    synthetic.set_xticks([1.0, 1.25, 1.5, 1.75])
    synthetic.set_xlabel("Diagnostic value", labelpad=9)
    style_axis(synthetic)

    figure.text(
        0.04,
        0.105,
        "Compare directions modulo a turn; retain raw\nchecks when the real coefficient matters.",
        fontsize=9.5,
    )
    figure.text(
        0.54,
        0.105,
        "Neither approach detects every angular pattern.\nInspect per-chain angular distributions.",
        fontsize=9.5,
    )
    figure.text(
        0.04,
        0.025,
        "No shared convergence cutoff is implied. Full results and definitions are in the validation report.",
        fontsize=8.5,
        color="#667085",
    )
    return figure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path, default=Path("docs/validation/circular_rhat.json")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("docs/validation/circular_rhat_summary.svg")
    )
    parser.add_argument(
        "--preview", type=Path, help="Optional PNG path for visual inspection"
    )
    args = parser.parse_args()
    result = json.loads(args.input.read_text())
    figure = make_figure(result)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, metadata={"Date": None})
    if args.output.suffix.lower() == ".svg":
        # Matplotlib emits trailing spaces in SVG styles; keep the artifact diff-clean.
        args.output.write_text(
            "\n".join(line.rstrip() for line in args.output.read_text().splitlines())
            + "\n"
        )
    if args.preview is not None:
        args.preview.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(args.preview, dpi=160)
    plt.close(figure)


if __name__ == "__main__":
    main()
