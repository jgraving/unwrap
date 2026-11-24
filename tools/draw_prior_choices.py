"""Regenerate the figures for docs/priors.md using NumPy, SciPy and Matplotlib.

Source: jgraving/circstats, elife-template.tex at 9bef710, sections
"Priors for Circular GLMs" and "Mixture weights — priors and identifiability".
"""

from pathlib import Path

import matplotlib
import numpy as np
from scipy.special import expit
from scipy.stats import halfnorm, lognorm, norm, vonmises

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ASSETS = Path(__file__).resolve().parents[1] / "docs" / "assets"
TEAL, RUST, VIOLET = "#126a62", "#96543b", "#564b84"
STYLES = ((TEAL, "-"), (RUST, "--"), (VIOLET, "-."))
DEGREES = np.linspace(-180, 180, 721)
RADIANS = np.deg2rad(DEGREES)


def panels(left, right):
    fig, axes = plt.subplots(1, 2, figsize=(8.7, 3.2), layout="constrained")
    for ax, title in zip(axes, (left, right), strict=True):
        ax.set_title(title, loc="left", fontweight="bold", fontsize=12, pad=9)
        ax.grid(axis="y", alpha=0.17)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
    return fig, axes


def degree_axis(ax, label):
    ax.set(xlim=(-180, 180), ylim=(0, None), xlabel=label, ylabel="Density per radian")
    ticks = [-180, -90, 0, 90, 180]
    ax.set_xticks(ticks, [f"{t}°" for t in ticks])


def save(fig, name, title):
    destination = ASSETS / f"prior-{name}.svg"
    fig.savefig(destination, metadata={"Title": title, "Date": None})
    # Matplotlib inserts trailing spaces in SVG paths. Keep the assets diff-friendly.
    lines = destination.read_text(encoding="utf-8").splitlines()
    destination.write_text("\n".join(line.rstrip() for line in lines) + "\n", encoding="utf-8")
    plt.close(fig)


def directions():
    fig, (left, right) = panels("Which baseline direction?", "How many turns across the design?")
    for kappa, center, label, (color, style) in zip(
        (np.finfo(float).eps, 4.0), (0, np.pi / 4), ("Any direction", "Cue at 45°"), STYLES
    ):
        left.plot(DEGREES, vonmises.pdf(RADIANS, kappa, loc=center),
                  color=color, ls=style, label=label)
    degree_axis(left, "Candidate mean direction")
    left.legend(loc="upper left")

    x = np.linspace(0, 1, 401)
    for rate, (color, style) in zip((90, 450), STYLES):
        angle = (rate * x + 180) % 360 - 180
        # Break the line at the plotting seam, where the physical path is continuous.
        jumps = np.flatnonzero(np.abs(np.diff(angle)) > 180) + 1
        angle[jumps] = np.nan
        right.plot(x, angle, color=color, ls=style, label=f"Slope {rate}° per unit")
    right.set(xlabel="Predictor", ylabel="Wrapped mean direction", xlim=(0, 1),
              ylim=(-195, 195), xticks=(0, 0.25, 0.5, 0.75, 1))
    ticks = [-180, -90, 0, 90, 180]
    right.set_yticks(ticks, [f"{t}°" for t in ticks])
    right.legend(loc="upper right")
    save(fig, "angles", "Baseline directions and rotations over a predictor")


def concentrations():
    fig, (left, right) = panels("A prior before the link", "What spread does it predict?")
    # Inverse softplus of the desired median concentration, with Normal SD 1.
    center = np.log(np.expm1(2.0))
    eta = np.linspace(-3.2, 5.5, 300)
    left.plot(eta, np.logaddexp(0, eta), color="#3a4746")
    quantiles = norm.ppf([0.05, 0.5, 0.95], loc=center, scale=1)
    kappas = np.logaddexp(0, quantiles)
    left.axvspan(quantiles[0], quantiles[-1], color=TEAL, alpha=0.10)
    for q, k, pct, (color, style) in zip(quantiles, kappas, [5, 50, 95], STYLES, strict=True):
        left.plot(q, k, "o", color=color)
        left.annotate(f"{pct}%", (q, k), xytext=(-6, 12), textcoords="offset points",
                      ha="right", color=color)
        right.plot(DEGREES, vonmises.pdf(RADIANS, k), color=color, ls=style,
                   label=f"{pct}th percentile  ·  κ = {k:.1f}")
    left.set(xlabel="Concentration intercept before softplus",
             ylabel="Observation concentration κ", ylim=(0, None))
    degree_axis(right, "Heading relative to mean direction")
    right.legend(loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=1)
    save(fig, "concentration", "A pre-softplus Normal prior and its implied heading spreads")


def alignment():
    fig, (left, right) = panels("Plausible between-animal alignment", "Where could an animal mean point?")
    width = np.log(100 / 10) / (2 * norm.ppf(0.975))
    hyperpriors = (
        (halfnorm(scale=2), "Half-Normal (scale 2)"),
        (lognorm(s=width, scale=np.sqrt(10 * 100)), "LogNormal (95% from 10 to 100)"),
        (lognorm(s=0.1, scale=30), "LogNormal (log SD 0.1)"),
    )
    kappa = np.geomspace(0.1, 200)
    # Integrate the conditional offset density over each hyperprior by its quantiles.
    probabilities = (np.arange(1024) + 0.5) / 1024
    for (prior, label), (color, style) in zip(hyperpriors, STYLES, strict=True):
        left.plot(kappa, prior.cdf(kappa), color=color, ls=style, label=label)
        draws = prior.ppf(probabilities)
        offsets = vonmises.pdf(RADIANS[:, None], draws).mean(axis=1)
        right.plot(DEGREES, offsets, color=color, ls=style)
    left.set(xscale="log", xlim=(0.1, 200), ylim=(0, 1.03),
             xlabel="Between-animal concentration $\\kappa_\\mu$ (log scale)",
             ylabel="Prior cumulative probability")
    degree_axis(right, "Offset from population direction")
    fig.legend(*left.get_legend_handles_labels(), loc="outside lower center", ncol=2)
    save(fig, "alignment", "Hyperpriors for between-animal concentration and angular offsets")


def mixture_weights():
    fig, (left, right) = panels("How much of each mode?", "What could one animal generate?")
    scales = [0.5, 1, 2]
    for row, scale in enumerate(scales):
        low, high = expit(norm.ppf([0.025, 0.975], scale=scale))
        left.plot([low, high], [row, row], color=VIOLET, lw=3, marker="|", ms=12)
        left.plot(0.5, row, "o", color=VIOLET)
        left.text(0.5, row + 0.12, f"{low:.2f} to {high:.2f}", ha="center", fontsize=9)
    left.set(yticks=range(3), yticklabels=["0.5", "1", "2"], ylim=(-0.35, 2.42),
             xlim=(0, 1), xticks=(0, 0.25, 0.5, 0.75, 1), ylabel="SD on the logit scale",
             xlabel="Probability of choosing +90° mode")
    positive = vonmises.pdf(RADIANS, 6, loc=np.pi / 2)
    negative = vonmises.pdf(RADIANS, 6, loc=-np.pi / 2)
    for weight, (color, style) in zip((0.50, 0.98), STYLES):
        density = weight * positive + (1 - weight) * negative
        right.plot(DEGREES, density, color=color, ls=style, label=f"Weight {weight:.2f}")
    degree_axis(right, "Heading")
    right.legend(loc="upper left")
    save(fig, "mixture", "A prior on animal mixture weights and its consequences for headings")


if __name__ == "__main__":
    with plt.rc_context({
        "svg.fonttype": "none", "svg.hashsalt": "unwrap-priors", "font.size": 10,
        "text.color": "#192725", "axes.labelcolor": "#192725", "axes.edgecolor": "#66736e",
        "xtick.color": "#43534f", "ytick.color": "#43534f",
        "legend.frameon": False, "legend.fontsize": 8.5, "lines.linewidth": 2,
    }):
        directions()
        concentrations()
        alignment()
        mixture_weights()
