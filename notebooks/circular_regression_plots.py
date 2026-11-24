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
from unwrap.circular import wrap_pi

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
    _plots.save_by_path(fig, FIG_DIR, save_path)


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


_FIT = "#086b72"
_TRUTH = "#a33b12"


def _bearing_degrees(theta):
    """Lift one turn to [-90, 270) degrees for this experiment's bearing axis."""
    return np.degrees(np.pi / 2 + wrap_pi(theta - np.pi / 2))


def _direction_summary(draws):
    """Circular center and central 90% arc, in display degrees at each x."""
    center = np.angle(np.exp(1j * draws).mean(axis=0))
    offsets = np.quantile(wrap_pi(draws - center), [0.05, 0.95], axis=0)
    center_deg = _bearing_degrees(center)
    lo, hi = center_deg + np.degrees(offsets)
    return center_deg, lo, hi


def _cue_axes(axes):
    for ax in axes:
        ax.set(xlabel="Cue strength (x)", xticks=[-1, -0.5, 0, 0.5, 1])
        ax.set_axisbelow(True)
        ax.grid(axis="y", color="0.90", linewidth=0.7)
        ax.spines[["right", "top"]].set_visible(False)
    axes[0].set(ylim=(-90, 270), yticks=[-90, 0, 90, 180, 270],
                ylabel="Bearing (degrees)")


def plot_predictions(data, grid, mu_draws, kappa_draws, true_mu, true_kappa, *,
                     save=True, save_path="figures/circular_regression/predictions.svg"):
    """Compare both fitted curves and pointwise credible bands with the simulation."""
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    _cue_axes(axes)
    axes[0].set_title("Where do animals point?", loc="left")
    axes[1].set_title("How tightly do they align?", loc="left")
    axes[1].set_ylabel("Concentration (κ)")
    points = axes[0].scatter(data.x, _bearing_degrees(data.y), s=8, linewidths=0,
                             alpha=0.25, color="0.38", label="Individual headings")
    mu_summary = _direction_summary(mu_draws)
    lo, med, hi = np.quantile(kappa_draws, [0.05, 0.5, 0.95], axis=0)
    for ax, (center, low, high), truth in zip(
        axes, (mu_summary, (med, lo, hi)), (_bearing_degrees(true_mu), true_kappa)
    ):
        band = ax.fill_between(grid, low, high, color=_FIT, alpha=0.24,
                               linewidth=0, label="90% credible band")
        estimate, = ax.plot(grid, center, color=_FIT, linewidth=2, label="Fitted")
        reference, = ax.plot(grid, truth, "--", color=_TRUTH, linewidth=1.5,
                             label="Generating curve")
    axes[1].set_ylim(bottom=0)
    fig.legend(handles=[points, estimate, band, reference], loc="lower center",
               ncols=4, frameon=False, fontsize=9)
    fig.tight_layout(rect=(0, 0.10, 1, 1))
    _save(fig, save, save_path)
    return fig


def plot_ppc(levels, observed, predicted, *, save=True,
             save_path="figures/circular_regression/ppc.svg"):
    """Observed vs replicated circular means and resultant lengths by cue strength."""
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    _cue_axes(axes)
    axes[0].set_title("Group mean bearing", loc="left")
    axes[1].set_title("Within-group alignment", loc="left")
    axes[1].set(ylabel="Mean resultant length", ylim=(0, 1))
    mu_summary = _direction_summary(np.angle(predicted))
    lo, med, hi = np.quantile(np.abs(predicted), [0.05, 0.5, 0.95], axis=0)
    for ax, (center, low, high), obs in zip(
        axes, (mu_summary, (med, lo, hi)), (_bearing_degrees(np.angle(observed)), np.abs(observed))
    ):
        interval = ax.vlines(levels, low, high, color=_FIT, linewidth=3,
                             label="90% predictive interval")
        dots = ax.scatter(levels, center, color=_FIT, s=18, zorder=3,
                          label="Predictive center")
        diamonds = ax.scatter(levels, obs, color=_TRUTH, edgecolors="white",
                              linewidths=0.6, marker="D", s=33, zorder=4, label="Observed group")
    fig.legend(handles=[interval, dots, diamonds], loc="lower center", ncols=3,
               frameon=False, fontsize=9)
    fig.tight_layout(rect=(0, 0.10, 1, 1))
    _save(fig, save, save_path)
    return fig
