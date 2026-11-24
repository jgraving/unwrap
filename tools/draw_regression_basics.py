"""Draw conditional heading densities and an exact circular partial-pooling example."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy.stats import vonmises


ASSETS = Path(__file__).resolve().parents[1] / "docs" / "assets"
TEAL = "#08665F"
RUST = "#A74A25"
PURPLE = "#513F7D"
INK = "#203331"
GRID = "#C2CECA"


def compass(ax, radius):
    """Bearings clockwise from the top, signed on the same compass for every panel."""
    ax.set_theta_direction(-1)
    ax.set_theta_zero_location("N")
    ax.set_xticks(np.deg2rad([0, 90, 180, 270]), ["0°", "90°", "±180°", "−90°"])
    ax.tick_params(axis="x", pad=0, labelsize=10, colors=INK)
    ax.set_ylim(0, radius)
    ax.grid(color=GRID, linewidth=0.7)
    ax.spines["polar"].set_color(GRID)
    ax.set_axisbelow(True)


def save(fig, basename, description):
    path = ASSETS / f"{basename}.svg"
    fig.savefig(
        path,
        metadata={"Date": None, "Title": basename.replace("-", " "), "Description": description},
        facecolor="white",
    )
    plt.close(fig)
    path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines()) + "\n")


def draw_cue_effects():
    theta = np.linspace(-np.pi, np.pi, 721)
    baseline = vonmises.pdf(theta, 2, loc=np.deg2rad(-45))
    cases = [("Turn", 45, 2), ("Tighten", -45, 6), ("Turn and tighten", 45, 6)]
    fig, axes = plt.subplots(1, 3, figsize=(8.2, 3.5), subplot_kw={"projection": "polar"})
    fig.subplots_adjust(left=0.055, right=0.945, bottom=0.24, top=0.80, wspace=0.49)
    for ax, (title, mu_deg, kappa) in zip(axes, cases, strict=True):
        compass(ax, 1)
        ax.set_yticks([0.5, 1.0], ["", ""])
        density = vonmises.pdf(theta, kappa, loc=np.deg2rad(mu_deg))
        ax.fill(theta, density, color=TEAL, alpha=0.10)
        ax.plot(theta, density, color=TEAL, linewidth=2.5, label="Cue on")
        ax.plot(theta, baseline, color=RUST, linestyle="--", linewidth=1.8, label="Cue off")
        ax.set_title(title, fontsize=13, color=INK, weight="semibold", pad=14)
        ax.text(0.5, -0.23, f"On: μ = {mu_deg}°, κ = {kappa}", transform=ax.transAxes,
                ha="center", va="center", fontsize=11, color=INK)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles[::-1], labels[::-1], loc="lower center", bbox_to_anchor=(0.5, 0.005),
               ncol=2, frameon=False, handlelength=2.8)
    save(fig, "regression-cue-effects", "Three circular density plots share a cue-off density. "
         "The cue can rotate the peak, sharpen it, or do both.")


def draw_pooling():
    population = np.deg2rad(165)
    kappa_between, kappa_within = 5, 2
    samples = [[-30, 30], [-45, -35, -25, -20, -15, -5, 5, 15, 20, 25, 35, 45]]
    fig, axes = plt.subplots(1, 2, figsize=(6.7, 4.3), subplot_kw={"projection": "polar"})
    fig.subplots_adjust(left=0.075, right=0.925, bottom=0.34, top=0.80, wspace=0.49)
    for ax, offsets in zip(axes, samples, strict=True):
        y = np.deg2rad(-135 + np.array(offsets))
        # Independent von Mises observations and a von Mises prior on the animal mean.
        data_vector = kappa_within * np.exp(1j * y).sum()
        posterior = kappa_between * np.exp(1j * population) + data_vector
        center = np.angle(posterior)
        half_width = vonmises.ppf(0.90, abs(posterior))
        compass(ax, 1.05)
        ax.set_yticks([])
        ax.plot([population, population], [0, 0.92], color=TEAL, linestyle="--", linewidth=1.8)
        ax.plot([center, center], [0, 0.76], color=PURPLE, linewidth=2.6)
        arc = np.linspace(center - half_width, center + half_width, 181)
        ax.plot(arc, np.full_like(arc, 0.76), color=PURPLE, linewidth=5, solid_capstyle="round")
        ax.scatter(y, np.full_like(y, 0.96), s=25, color=RUST, zorder=4)
        ax.set_title(f"{len(y)} headings", color=INK, fontsize=14, weight="semibold", pad=14)
        ax.text(0.5, -0.24, f"Pooled mean {np.rad2deg(center):.0f}°", transform=ax.transAxes,
                ha="center", va="center", color=PURPLE, fontsize=12)
    labels = [Line2D([], [], marker="o", color=RUST, linestyle="none", label="Observed heading"),
              Line2D([], [], color=TEAL, linestyle="--", label="Population mean"),
              Line2D([], [], color=PURPLE, linewidth=3, label="Animal mean and 80% credible arc")]
    fig.legend(handles=labels[:2], loc="lower center", bbox_to_anchor=(0.5, 0.075),
               frameon=False, ncol=2)
    fig.legend(handles=labels[2:], loc="lower center", bbox_to_anchor=(0.5, 0.005), frameon=False)
    save(fig, "regression-partial-pooling", "Two animals with the same observed mean borrow "
         "different amounts of information from their shared population direction.")


if __name__ == "__main__":
    # Keep labels selectable and regenerated paths stable in version control.
    with mpl.rc_context({"svg.fonttype": "none", "svg.hashsalt": "unwrap-regression-basics"}):
        draw_cue_effects()
        draw_pooling()
