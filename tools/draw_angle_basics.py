"""Regenerate the circular geometry and von Mises teaching figures."""

from io import StringIO
from pathlib import Path
from xml.sax.saxutils import escape

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle
from scipy.special import i0e, i1e
from scipy.stats import vonmises

ASSETS = Path(__file__).resolve().parents[1] / "docs/assets"
INK = "#193131"
TEAL = "#08726d"
RUST = "#a54b2d"
RULE = "#adc3c0"


def save(fig, name, title, description):
    svg = StringIO()
    fig.savefig(svg, format="svg", metadata={"Title": title, "Date": None})
    content = svg.getvalue().replace(
        "<svg ", '<svg role="img" aria-labelledby="title desc" ', 1,
    ).replace(
        "<title>", f'<desc id="desc">{escape(description)}</desc>\n <title id="title">', 1,
    )
    (ASSETS / f"{name}.svg").write_text(
        "\n".join(line.rstrip() for line in content.splitlines()) + "\n", encoding="utf-8",
    )
    plt.close(fig)


def dial(ax):
    """Use Cartesian angles throughout: zero at the right, positive anticlockwise."""
    ax.set(aspect="equal", xlim=(-1.25, 1.25), ylim=(-1.2, 1.2))
    ax.axis("off")
    ax.axhline(0, color=RULE, lw=0.8)
    ax.axvline(0, color=RULE, lw=0.8)
    ax.add_patch(Circle((0, 0), 1, edgecolor=RULE, facecolor="none", lw=1.3))


def arrow(ax, endpoint, color, width):
    ax.annotate(
        "", xy=endpoint, xytext=(0, 0),
        arrowprops=dict(arrowstyle="-|>", color=color, lw=width,
                        mutation_scale=11, shrinkA=0, shrinkB=0),
    )


def angle_views():
    fig = plt.figure(figsize=(5.6, 2.55), layout="constrained")
    strips, circle = fig.subplots(1, 2, width_ratios=[1.55, 1])
    strips.set(xlim=(-0.08, 1.08), ylim=(-0.35, 1.4))
    strips.axis("off")
    for y, heading, ticks, values in (
        (1, "Cut at 0°", ("0°", "180°", "360°"), (10 / 360, 350 / 360)),
        (0, "Cut at −180°", ("−180°", "0°", "180°"), (170 / 360, 190 / 360)),
    ):
        strips.text(0, y + 0.26, heading, weight="bold")
        strips.plot([0, 1], [y, y], lw=1.5, color=RULE)
        for x, label in zip((0, 0.5, 1), ticks):
            strips.plot([x, x], [y - 0.045, y + 0.045], lw=1.2, color=INK)
            strips.text(x, y - 0.13, label, ha="center", va="top", fontsize=12.5)
        strips.scatter(values, [y, y], color=TEAL, s=42, zorder=3)
    dial(circle)
    circle.set_title("One circle", weight="bold", pad=6)
    for angle, label, dy in ((10, "10°", 0.21), (350, "350°", -0.27)):
        x, y = np.cos(np.deg2rad(angle)), np.sin(np.deg2rad(angle))
        circle.plot([0, x], [0, y], color=TEAL, lw=1.6)
        circle.scatter(x, y, color=TEAL, s=35)
        circle.text(x - 0.10, y + dy, label, ha="right", color=TEAL, weight="bold")
    circle.text(-0.06, 0.85, "90°", ha="right", va="top", fontsize=12)
    circle.text(-1.02, -0.20, "180°", ha="left", fontsize=12)
    save(fig, "angle-basics", "Two headings on a number line and on the circle",
         "The same 350° and 10° points appear near opposite ends of a line cut at zero, "
         "together near the middle of a line cut at minus 180°, and 20° apart on a circle.")


def resultant_lengths():
    fig, axes = plt.subplots(1, 3, figsize=(5.5, 2.1), layout="constrained")
    for ax, spread in zip(axes, (10, 60, 90)):
        dial(ax)
        angles = np.deg2rad([-spread, spread])
        vectors = np.column_stack((np.cos(angles), np.sin(angles)))
        for v in vectors:
            ax.plot([0, v[0]], [0, v[1]], lw=1.6, color=TEAL)
            ax.scatter(*v, s=32, color=TEAL, zorder=3)
        mean = vectors.mean(axis=0)
        length = np.linalg.norm(mean)
        if length > 1e-12:
            arrow(ax, mean, RUST, 2.4)
        else:
            ax.scatter(0, 0, s=32, color=RUST, zorder=3)
        ax.set_title(f"±{spread}°", pad=0, weight="bold")
        ax.text(0, -1.28, f"MRL = {length:.2f}", ha="center", va="top", color=RUST)
    save(fig, "mean-resultant", "Mean unit vectors for three pairs of headings",
         "Observed directions at plus and minus 10°, 60°, and 90° have mean vectors "
         "of length 0.98, 0.50, and zero. The last pair cancels exactly.")


def density_curves():
    fig, ax = plt.subplots(figsize=(5.6, 3.55), layout="constrained")
    degrees = np.linspace(-180, 180, 721)
    for mean, concentration, color, style in (
        (0, 1, TEAL, "--"),
        (0, 4, TEAL, "-"),
        (90, 4, RUST, "-."),
    ):
        density = vonmises.pdf(np.deg2rad(degrees - mean), concentration)
        ax.plot(degrees, density, color=color, ls=style, lw=2.4,
                label=f"μ = {mean}°, κ = {concentration}")
    ax.set(xlim=(-180, 180), ylim=(0, 0.82), yticks=[0, 0.2, 0.4, 0.6, 0.8],
           xlabel="Heading", ylabel="Density per radian")
    ax.set_xticks([-180, -90, 0, 90, 180], ["−180°", "−90°", "0°", "90°", "180°"])
    ax.legend(loc="upper left", frameon=False, fontsize=12.5, handlelength=2.1)
    ax.grid(axis="y", color=RULE, alpha=0.45)
    save(fig, "von-mises-density", "Rotation and concentration of a von Mises density",
         "Two curves peak at zero. Raising kappa from one to four makes the peak "
         "narrower and taller. Changing the mean to 90° shifts that same narrow peak.")


def concentration_curve():
    fig, ax = plt.subplots(figsize=(5.6, 3.1), layout="constrained")
    kappa = np.linspace(0, 10, 401)
    ax.plot(kappa, i1e(kappa) / i0e(kappa), color=TEAL, lw=2.5)
    for k in (1, 2, 4, 8):
        r = i1e(k) / i0e(k)
        ax.vlines(k, 0, r, color=RULE, ls=":", lw=1)
        ax.scatter(k, r, color=TEAL, s=35, zorder=3)
        ax.annotate(f"{r:.2f}", (k, r), xytext=(5, -11), textcoords="offset points",
                    va="top", color=TEAL, size=12.5)
    ax.set(xlim=(0, 10), ylim=(0, 1.02), xticks=[0, 1, 2, 4, 6, 8, 10],
           yticks=[0, 0.25, 0.5, 0.75, 1], xlabel="Concentration κ",
           ylabel="Population MRL")
    ax.grid(axis="y", color=RULE, alpha=0.45)
    save(fig, "von-mises-mrl", "Concentration mapped to the population mean resultant length",
         "The curve starts at zero and approaches one. Kappa values one, two, four, "
         "and eight correspond to mean resultant lengths 0.45, 0.70, 0.86, and 0.94.")


def main():
    with mpl.rc_context({
        "svg.hashsalt": "unwrap-angle-basics",
        "svg.fonttype": "none",
        "text.color": INK,
        "axes.labelcolor": INK,
        "axes.edgecolor": RULE,
        "xtick.color": INK,
        "ytick.color": INK,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "font.size": 13,
        "font.family": "DejaVu Sans",
        "axes.axisbelow": True,
    }):
        angle_views()
        resultant_lengths()
        density_curves()
        concentration_curve()


if __name__ == "__main__":
    main()
