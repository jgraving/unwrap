"""Figures for notebooks/circular_regression.ipynb (the tutorial notebook).

The three tutorial figures, on the shared house style from ``unwrap.plots``
(100 dpi inline display / 300 dpi vector PDFs under
``figures/circular_regression/``):

* :func:`plot_trace` — ArviZ traces of the fixed effects and the two kappa
  coefficients (compact layout, as in the original cell);
* :func:`plot_posterior_curves` — posterior circular mean direction and
  mean concentration over the predictor grid, against the generating
  curves;
* :func:`plot_replicated_mean_resultant` — the mean-resultant length of the
  replicated datasets with the observed value.

The posterior bookkeeping these figures display (prediction grids,
link-inverted draws, the replicated-data statistics) stays in the notebook
as post-hoc analysis; only the figure construction lives here.
"""

from pathlib import Path

import arviz as az
import matplotlib.pyplot as plt
import numpy as np

from unwrap import plots as _plots

#: Publication figures for this study live in <repo>/figures/circular_regression/.
FIG_DIR = Path(__file__).resolve().parent.parent / "figures" / "circular_regression"


def set_style():
    """Apply the shared figure style for this study (figures under ``FIG_DIR``).

    No study-specific rcParams: the original tutorial used the plain house
    look, which the shared style already is.
    """
    _plots.set_style(FIG_DIR)


def plot_trace(idata, var_names, name="trace"):
    """ArviZ traces for the fit, compact layout, saved under ``FIG_DIR``.

    ``var_names`` lists the parameters to trace (scalars: no per-index
    expansion), so a vector PDF save stays small.
    """
    az.plot_trace(idata, var_names=var_names, compact=True)
    fig = plt.gcf()
    fig.tight_layout()
    _plots.save_figure(fig, FIG_DIR, name)
    return fig


def plot_posterior_curves(data, grid, mean_direction, concentration_draws,
                          name="posterior_curves"):
    """Posterior direction / concentration curves against the generators.

    ``grid`` is the predictor grid, ``mean_direction`` the circular mean of
    the predicted directions at the grid, ``concentration_draws`` the
    posterior softplus concentration draws at the grid, and ``data`` the
    observed ``(x, y)`` frame.
    """
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].scatter(data.x, data.y, s=12, alpha=0.3, color="gray", label="Observed")
    axes[0].plot(grid, mean_direction, label="Posterior circular mean")
    axes[0].plot(grid, grid, "--", label="Generating mean")
    axes[0].set(xlabel="Predictor x", ylabel="Direction (radians)", ylim=(-np.pi, np.pi))
    axes[1].plot(grid, concentration_draws.mean(axis=0), label="Posterior mean")
    axes[1].plot(grid, np.logaddexp(0, 2 + 0.5 * grid), "--", label="Generating concentration")
    axes[1].set(xlabel="Predictor x", ylabel="Concentration kappa")
    for ax in axes:
        ax.legend(fontsize=8)
    fig.tight_layout()
    _plots.save_figure(fig, FIG_DIR, name)
    return fig


def plot_replicated_mean_resultant(replicated_R, observed_R,
                                   name="replicated_mean_resultant"):
    """Histogram of the replicated datasets' mean resultant length.

    ``replicated_R`` is the per-replicate mean resultant length
    (``(S,)``) and ``observed_R`` the scalar value for the observed data.
    """
    fig, ax = plt.subplots(figsize=(6, 3))
    ax.hist(replicated_R.ravel(), bins=30, alpha=0.7, label="Replicated datasets")
    ax.axvline(observed_R, color="black", label="Observed dataset")
    ax.set(xlabel="Mean resultant length", ylabel="Replicated dataset count")
    ax.legend()
    fig.tight_layout()
    _plots.save_figure(fig, FIG_DIR, name)
    return fig
