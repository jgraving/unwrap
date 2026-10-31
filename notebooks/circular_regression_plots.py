"""Figure helpers for the circular_regression example notebook.

The notebook keeps the statistical code (data generation, the ``bambi`` von
Mises fit, and the posterior/predictive computations); this module holds the
matplotlib "how" (axes, lines, legends, layout) and the export of the
figures.

Style and saving follow the shared house convention from ``unwrap.plots``:
figures render at display resolution (100 dpi) in the notebook cell and are
saved as 300 dpi vector files under ``<repo>/figures/circular_regression/``.
Every public figure function accepts a ``save`` flag (on by default) and an
explicit repo-relative ``save_path`` (the notebook always passes it); the
written files are gitignored.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import arviz as az

from unwrap import plots as _plots

#: Publication figures for this study live in
#: <repo>/figures/circular_regression/, anchored on this file's location so
#: the output does not depend on the kernel cwd (the notebook may run from
#: notebooks/ or from the repository root).
FIG_DIR = Path(__file__).resolve().parent.parent / "figures" / "circular_regression"

#: Study-specific rcParams on top of the house style (the house style already
#: sets Arial, the 100/300 dpi tiers, type-42 PDF outlines and hairline
#: axes/ticks).
_EXTRA_RC = {
    "svg.fonttype": "none",
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
}


def set_style():
    """Consistent figure style across the notebook (house style + the study's
    SVG/hairline additions; figures under ``FIG_DIR``)."""
    _plots.set_style(FIG_DIR, _EXTRA_RC)


def _save(fig, save, save_path):
    """Write ``fig`` to ``save_path`` (a repo-relative path such as
    ``figures/circular_regression/x.svg``) through the house two-res saver.
    The file format follows the path's suffix."""
    if not save:
        return
    p = Path(save_path)
    if p.suffix.lower() == ".svg":
        _plots.save_figure(fig, FIG_DIR, p.stem, pdf=False, png=False, svg=True)
    else:
        _plots.save_figure(fig, FIG_DIR, p.stem, pdf=True, png=False, svg=False)


def plot_trace(idata, var_names, *, save=True,
               save_path="figures/circular_regression/trace.svg"):
    """ArviZ trace plot of the posterior for the named variables.

    arviz 1.x: ``az.plot_trace`` returns a PlotCollection (the old grid-based
    layout and its ``compact`` argument were removed). The collection renders
    itself with ``pc.show()``, so the cell is the single call below (no
    separate ``plt.show()``).
    """
    pc = az.plot_trace(idata, var_names=var_names)
    if save:
        # Anchor on FIG_DIR exactly like _save (the kernel cwd is not
        # reliable): the written file is FIG_DIR/<name of save_path>.
        # savefig.dpi (300) comes from the house style set by set_style.
        p = FIG_DIR / Path(save_path).name
        p.parent.mkdir(parents=True, exist_ok=True)
        pc.savefig(p)
        print(f"Saved figure to: {_plots.repo_rel(p)}")
    pc.show()
    return pc


def plot_predictions(data, grid, mean_direction, concentration_draws, *,
                     save=True,
                     save_path="figures/circular_regression/predictions.svg"):
    """Observed points plus the posterior mean-direction and concentration
    curves across the predictor. The dashed reference lines use the known
    generating values (mean = x, kappa = softplus(2 + 0.5 x))."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].scatter(data.x, data.y, s=12, alpha=0.3, color="gray", label="Observed")
    axes[0].plot(grid, mean_direction, label="Posterior circular mean")
    axes[0].plot(grid, grid, "--", label="Generating mean")
    axes[0].set(xlabel="Predictor x", ylabel="Direction (radians)", ylim=(-np.pi, np.pi))
    axes[1].plot(grid, concentration_draws.mean(axis=0), label="Posterior mean")
    axes[1].plot(grid, np.logaddexp(0, 2 + 0.5 * grid), "--",
                 label="Generating concentration")
    axes[1].set(xlabel="Predictor x", ylabel="Concentration kappa")
    for ax in axes:
        ax.legend(fontsize=8)
    fig.tight_layout()
    _save(fig, save, save_path)
    return fig


def plot_ppc(replicated_R, observed_R, *, save=True,
             save_path="figures/circular_regression/ppc.svg"):
    """Histogram of replicated-dataset mean resultant lengths against the
    observed value (a circular-concentration-style PPC)."""
    fig, ax = plt.subplots(figsize=(6, 3))
    ax.hist(replicated_R.ravel(), bins=30, alpha=0.7, label="Replicated datasets")
    ax.axvline(observed_R, color="black", label="Observed dataset")
    ax.set(xlabel="Mean resultant length", ylabel="Replicated dataset count")
    ax.legend()
    fig.tight_layout()
    _save(fig, save, save_path)
    return fig
