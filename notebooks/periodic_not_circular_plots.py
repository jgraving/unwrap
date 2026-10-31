"""Figure helpers for the periodic_not_circular example notebook.

The notebook keeps the statistical code (simulating the rhythm, fitting the
periodic model) and the raw posterior predictive draws; the summarised
quantity (median / HDI) and all matplotlib mechanics live here.
"""
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import arviz as az


def set_style():
    """Global figure style used by this notebook (Arial text and math)."""
    plt.rcParams.update({
        "font.family": "Arial",
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
    })


def plot_beta_fit_ppc(time, activity, y):
    """A 24-hour activity rhythm: the simulated observed activity with the
    model's posterior predictive median and 90% HDI overlaid.

    Parameters
    ----------
    time : array of shape (T,)
        Hours in the day.
    activity : array of shape (T,)
        The simulated observed activity.
    y : array of shape (S, T)
        Posterior predictive draws for the activity (one row per draw).
    """
    median = np.median(y, axis=0)
    hdi90 = az.hdi(y, hdi_prob=0.90)

    fig, ax = plt.subplots(figsize=(4.2, 2.4), dpi=300)

    # Transparent backgrounds
    fig.patch.set_facecolor("none")
    fig.patch.set_alpha(0)
    ax.set_facecolor("none")
    ax.patch.set_facecolor("none")
    ax.patch.set_alpha(0)

    # Observations: tiny points, light
    ax.plot(time, activity, ".", color="black", markersize=1.5, alpha=0.35, zorder=2)
    # 90% HDI band (warm so it does not clash with a coloured background)
    ax.fill_between(time, hdi90[:, 0], hdi90[:, 1], color="tab:orange",
                    alpha=0.18, linewidth=0, zorder=1)
    # Median line
    ax.plot(time, median, color="black", linewidth=1.3, zorder=3)

    # Minimal cosmetics
    ax.set_xlabel("Time (hours)", fontsize=8)
    ax.set_ylabel("Activity", fontsize=8)
    ax.tick_params(labelsize=7, width=0.8, length=3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(0, 24)
    ax.set_xticks([0, 6, 12, 18, 24])

    fig.tight_layout()
    out = Path(__file__).resolve().parent / "figures" / "periodic_beta_fit.pdf"
    out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out, bbox_inches="tight", transparent=True)
    plt.show()
