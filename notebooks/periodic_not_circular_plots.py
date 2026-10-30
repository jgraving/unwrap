"""Figures for notebooks/periodic_not_circular.ipynb.

The study's single figure — the median 24-hour activity curve with the
90% HDI band over the ``activity_obs`` posterior-predictive samples, plus
the faint observed points — is drawn with the shared house look from
``unwrap.plots`` (100 dpi inline display / 300 dpi publication file under
``figures/periodic_not_circular/``). The PDF is saved transparent so the
figure composites onto the warm appendix background it was designed for.
"""

from pathlib import Path

import arviz as az
import matplotlib.pyplot as plt
import numpy as np

from unwrap import plots as _plots

#: Publication figure for this study lives in <repo>/figures/periodic_not_circular/.
FIG_DIR = Path(__file__).resolve().parent.parent / "figures" / "periodic_not_circular"

# Study-specific additions on top of the house style: the mathtext fonts
# the original notebook set (the shared style already covers font.family).
_EXTRA_RC = {
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
}


def set_style():
    """Apply the shared figure style for this study (figures under ``FIG_DIR``)."""
    _plots.set_style(FIG_DIR, _EXTRA_RC)


def plot_periodic_beta_fit(ppc, time, activity, name="periodic_beta_fit"):
    """Median activity curve + 90% HDI over the posterior predictive.

    ``ppc`` is the ``activity_obs`` posterior-predictive draw array
    (``(chain, draw, time)`` or ``(draw, time)``); ``time`` (hours) and
    ``activity`` are the simulated data the fit was run on. The figure is
    saved as a transparent vector PDF under ``figures/periodic_not_circular/``
    and returned.
    """
    # posterior predictive samples: (chain, draw, time) OR sometimes (draw, time) depending on config
    y = np.asarray(ppc["activity_obs"])
    if y.ndim == 3:
        y = y.reshape(-1, y.shape[-1])  # (S, T)
    elif y.ndim == 2:
        pass
    else:
        raise ValueError(f"Unexpected posterior predictive shape: {y.shape}")

    # Use a warm/orange for bands so it doesn't clash with blue background in appendix
    hdi_color = "tab:orange"

    median = np.median(y, axis=0)
    hdi90 = az.hdi(y, hdi_prob=0.90)  # (T, 2)

    fig, ax = plt.subplots(figsize=(4.2, 2.4))

    # Transparent backgrounds
    fig.patch.set_facecolor("none")
    fig.patch.set_alpha(0)
    ax.set_facecolor("none")
    ax.patch.set_facecolor("none")
    ax.patch.set_alpha(0)

    # Observations: tiny points, light
    ax.plot(time, activity, ".", color="black", markersize=1.5, alpha=0.35, zorder=2)

    # 90% HDI band
    ax.fill_between(time, hdi90[:, 0], hdi90[:, 1], color=hdi_color, alpha=0.18, linewidth=0, zorder=1)

    # Median line
    ax.plot(time, median, color="black", linewidth=1.3, zorder=3)

    # Minimal cosmetics
    ax.set_xlabel("Time (hours)", fontsize=8)
    ax.set_ylabel("Activity", fontsize=8)
    ax.tick_params(labelsize=7, width=0.8, length=3)

    # Clean spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Nicer x ticks for 24h
    ax.set_xlim(0, 24)
    ax.set_xticks([0, 6, 12, 18, 24])

    fig.tight_layout()
    _plots.save_figure(fig, FIG_DIR, name, pdf=True, png=False, transparent=True)
    return fig
