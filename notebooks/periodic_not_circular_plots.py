"""Figure helpers for the periodic_not_circular example notebook.

The notebook keeps the statistical code (simulating the rhythm, building
and fitting the periodic Beta model); this module holds the matplotlib
"how" (axes, bands, cosmetics) and the export of the figure.
``plot_beta_fit_ppc`` additionally samples the posterior predictive draws
it needs and summarises them (median / HDI), so its cell is a single call
in the notebook.

Style and saving follow the shared house convention from ``unwrap.plots``:
figures render at display resolution (100 dpi) in the notebook cell and are
saved as 300 dpi vector files under ``<repo>/figures/periodic_not_circular/``.
Every public figure function accepts a ``save`` flag (on by default) and an
explicit repo-relative ``save_path`` (the notebook always passes it); the
written files are gitignored.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import arviz as az
import pymc as pm

from unwrap import plots as _plots

#: Publication figures for this study live in
#: <repo>/figures/periodic_not_circular/, anchored on this file's location
#: so the output does not depend on the kernel cwd (the notebook may run
#: from notebooks/ or from the repository root).
FIG_DIR = Path(__file__).resolve().parent.parent / "figures" / "periodic_not_circular"


def set_style():
    """Consistent figure style across the notebook (the shared house
    convention from ``unwrap.plots``; figures under ``FIG_DIR``)."""
    _plots.set_style(FIG_DIR)


def _save(fig, save, save_path):
    """Write ``fig`` to ``save_path`` (a repo-relative path such as
    ``figures/periodic_not_circular/x.pdf``) through the house two-res
    saver. The file format follows the path's suffix."""
    if not save:
        return
    p = Path(save_path)
    if p.suffix.lower() == ".svg":
        _plots.save_figure(fig, FIG_DIR, p.stem, pdf=False, png=False, svg=True)
    else:
        _plots.save_figure(fig, FIG_DIR, p.stem, pdf=True, png=False, svg=False)


def plot_beta_fit_ppc(time, activity, model, trace, *, save=True,
                      save_path="figures/periodic_not_circular/periodic_beta_fit.pdf"):
    """A 24-hour activity rhythm: the simulated observed activity with the
    model's posterior predictive median and 50%/90% HDI overlaid.

    The posterior predictive draws are sampled here (one row per draw) and
    summarised, so the notebook cell is a single call.

    Parameters
    ----------
    time : array of shape (T,)
        Hours in the day.
    activity : array of shape (T,)
        The simulated observed activity.
    model : pm.Model
        The fitted model (the posterior predictive is drawn from it).
    trace : InferenceData
        The fitted model's posterior.
    """
    ppc = pm.sample_posterior_predictive(
        trace, model, var_names=["activity_obs"],
        return_inferencedata=False, progressbar=False,
    )
    # posterior predictive samples: (chain, draw, time) OR (draw, time)
    y = np.asarray(ppc["activity_obs"])
    if y.ndim == 3:
        y = y.reshape(-1, y.shape[-1])  # (S, T)
    elif y.ndim == 2:
        pass
    else:
        raise ValueError(f"Unexpected posterior predictive shape: {y.shape}")

    median = np.median(y, axis=0)
    # arviz 1.3 reads a 2-D array as (n, samples): the draws axis must be
    # last, so the time axis is put in front (one row per time point, HDI
    # across posterior draws per point). arviz 0.23 interpreted the same
    # (S, T) array the other way round, with a FutureWarning that the
    # interpretation would change - hence the transpose.
    hdi50 = az.hdi(y.T, prob=0.50)  # arviz 1.x: hdi_prob -> prob
    hdi90 = az.hdi(y.T, prob=0.90)

    fig, ax = plt.subplots(figsize=(4.2, 2.4))

    # Observations: tiny points, light
    ax.plot(time, activity, ".", color="black", markersize=1.5, alpha=0.35, zorder=2)
    # 50%/90% HDI bands (outer 90% first, inner 50% on top)
    ax.fill_between(time, hdi90[:, 0], hdi90[:, 1], color="tab:orange",
                    alpha=0.18, linewidth=0, zorder=1)
    ax.fill_between(time, hdi50[:, 0], hdi50[:, 1], color="tab:orange",
                    alpha=0.35, linewidth=0, zorder=1)
    # Median line
    ax.plot(time, median, color="black", linewidth=1.3, zorder=3)

    # Minimal cosmetics (house style supplies the rest: font, facecolor, dpi)
    ax.set_xlabel("Time (hours)", fontsize=8)
    ax.set_ylabel("Activity", fontsize=8)
    ax.tick_params(labelsize=7, width=0.8, length=3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(0, 24)
    ax.set_xticks([0, 6, 12, 18, 24])

    fig.tight_layout()
    _save(fig, save, save_path)
    return fig
