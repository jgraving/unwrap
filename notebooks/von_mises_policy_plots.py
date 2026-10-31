"""Figure helpers for the von_mises_policy example notebook.

The notebook defines the von Mises steering policy (goal heading, the von Mises
density, and the stepwise Monte-Carlo sampling) and produces, for each
concentration, the radial density profile and the simulated trajectories. This
module only draws the two-row figure from those quantities.
"""
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.cm import get_cmap

DPI = 300
FONT_SIZE = 7
AX_LW = 0.8
LINE_LW = 1.0
TICK_LW = 0.8
TICK_LEN = 2.5


def set_style():
    """Global figure style: Arial for normal text, STIX for math."""
    mpl.rcParams.update({
        "text.usetex": False,
        "font.family": "Arial",
        "mathtext.fontset": "stix",
        "mathtext.default": "it",
        "font.size": FONT_SIZE,
        "axes.labelsize": FONT_SIZE,
        "axes.titlesize": FONT_SIZE,
        "legend.fontsize": FONT_SIZE,
        "xtick.labelsize": FONT_SIZE,
        "ytick.labelsize": FONT_SIZE,
        "axes.linewidth": AX_LW,
        "lines.linewidth": LINE_LW,
        "xtick.major.width": TICK_LW,
        "ytick.major.width": TICK_LW,
        "xtick.major.size": TICK_LEN,
        "ytick.major.size": TICK_LEN,
        "savefig.dpi": DPI,
        "figure.dpi": DPI,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "path",
    })


def _clean_axis(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def plot_navigation(x_start, target, kappas, theta_grid, radials, trajectories,
                    R_policy=1.75, amp=1.0):
    """Two-row figure: the exact von Mises density as a radial modulation
    (top) and 10 stepwise Monte-Carlo trajectories (bottom), one column per
    concentration in ``kappas``.

    Parameters
    ----------
    x_start, target : array of shape (2,)
        Policy origin and goal.
    kappas : array of shape (K,)
        Concentrations, one column each.
    theta_grid : array of shape (T,)
        Angles on which the radial profile was evaluated.
    radials : list of K arrays of shape (T,)
        Radii of the density modulation at each ``theta_grid`` value.
    trajectories : list of K arrays of shape (T_steps + 1, M, 2)
        Simulated trajectories for each concentration.
    """
    n = len(kappas)
    mu0 = np.arctan2(target[1] - x_start[1], target[0] - x_start[0])
    Rmax = R_policy + amp
    goal_dir = np.array([np.cos(mu0), np.sin(mu0)], dtype=np.float32)
    target_top = x_start + (Rmax + 0.22) * goal_dir  # marker just outside the curve

    fig, axes = plt.subplots(2, n, figsize=(1.10 * n, 2.20), dpi=DPI, squeeze=False)

    # --- Row 0: exact von Mises PDF as a radial modulation ---
    for i, kappa in enumerate(kappas):
        ax = axes[0, i]
        xs = x_start[0] + radials[i] * np.cos(theta_grid)
        ys = x_start[1] + radials[i] * np.sin(theta_grid)
        ax.plot(xs, ys, color="k", lw=1.1, zorder=3)
        ax.add_patch(Circle(xy=x_start, radius=R_policy, fill=False, lw=AX_LW,
                            ec="k", alpha=0.95))
        ax.plot(x_start[0], x_start[1], "o", ms=2.2, c="k", zorder=5)
        ax.plot(target_top[0], target_top[1], "^", ms=3.5, c="k", zorder=5)
        ax.set_title(rf"$\kappa={float(kappa):g}$", pad=2)
        ax.set_aspect("equal", adjustable="box")
        _clean_axis(ax)
        pad = 0.30
        ax.set_xlim(x_start[0] - Rmax - pad, x_start[0] + Rmax + pad)
        ax.set_ylim(x_start[1] - Rmax - pad, x_start[1] + Rmax + pad)

    # --- Row 1: trajectories by stepwise von Mises sampling ---
    label_idx = 1
    for i, kappa in enumerate(kappas):
        ax = axes[1, i]
        pos = trajectories[i]
        M = pos.shape[1]
        cmap = get_cmap("tab10", M) if M <= 10 else get_cmap("tab20", M)
        for j in range(M):
            ax.plot(pos[:, j, 0], pos[:, j, 1], lw=0.55, alpha=0.95,
                    color=cmap(j), zorder=2)
        ax.plot(x_start[0], x_start[1], "o", ms=3.4, c="k", zorder=5)
        ax.plot(target[0], target[1], "^", ms=4.0, c="k", zorder=5)
        if i == label_idx:
            ax.annotate(r"Start $\mathbf{x}$", xy=(x_start[0], x_start[1]),
                        xytext=(8, 0), textcoords="offset points", ha="left",
                        va="center")
            ax.annotate(r"Target $\mathbf{g}$", xy=(target[0], target[1]),
                        xytext=(8, 0), textcoords="offset points", ha="left",
                        va="center")
        ax.set_aspect("equal", adjustable="box")
        _clean_axis(ax)
        ax.set_xlim(-7.5, 7.5)
        ax.set_ylim(-3.5, 10.8)

    axes[0, 0].set_ylabel("Policy density", labelpad=4)
    axes[1, 0].set_ylabel("Trajectories", labelpad=4)
    fig.subplots_adjust(left=0.10, right=0.995, bottom=0.06, top=0.93,
                        wspace=0.22, hspace=0.10)

    out_dir = Path(__file__).resolve().parent / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    save_path = out_dir / "vonmises_navigation.pdf"
    save_path_svg = out_dir / "vonmises_navigation.svg"
    fig.savefig(save_path, bbox_inches="tight")
    fig.savefig(save_path_svg, bbox_inches="tight")
    plt.show()
    print(f"Saved figure to: {save_path}")
    print(f"Saved figure to: {save_path_svg}")
