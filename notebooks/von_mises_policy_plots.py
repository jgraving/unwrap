"""Figure for notebooks/von_mises_policy.ipynb.

A single 2x4 figure: the top row renders the exact von Mises policy
densities (kappa in 0, 0.5, 2, 5) about the goal direction from the start
point, as radial modulations of a fixed step radius; the bottom row shows
the stepwise stochastic trajectories drawn by sampling headings around the
goal direction at each step. There is no dataset and no model fit here, so
this study has no preprocess module — the seeded trajectory simulation that
draws the bottom row lives with the figure.

Drawn with the shared house look from ``unwrap.plots`` (100 dpi inline
display / 300 dpi publication files under ``figures/von_mises_policy/``);
a vector PDF and a vector SVG are saved, the SVG for the notebook's
Inkscape workflow.
"""

from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

from unwrap import plots as _plots

#: Publication figures for this study live in <repo>/figures/von_mises_policy/.
FIG_DIR = Path(__file__).resolve().parent.parent / "figures" / "von_mises_policy"

# -----------------------------
# Config
# -----------------------------
X_START = np.array([0.0, 0.0], dtype=np.float32)
TARGET = np.array([0.0, 10.0], dtype=np.float32)

KAPPAS = np.array([0, 0.5, 2, 5], dtype=np.float32)

T_STEPS = 300
STEP_LEN = 0.15
NUM_TRAJ = 10
RNG_SEED = 12345

R_POLICY = 1.75
AMP = 1.0

# Study-specific additions on top of the house style: the STIX mathtext
# fonts and the compact 7 pt sizing of the original notebook (Arial family,
# fonttype-42 PDFs and the two resolution tiers come from the shared style).
_EXTRA_RC = {
    "text.usetex": False,
    # Math text (kappa, x, g, etc.)
    "mathtext.fontset": "stix",
    "mathtext.default": "it",
    "font.size": 7,
    "axes.labelsize": 7,
    "axes.titlesize": 7,
    "legend.fontsize": 7,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.0,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    # If you import SVG into Inkscape, paths are most robust
    "svg.fonttype": "path",
}


def set_style():
    """Apply the shared figure style for this study (figures under ``FIG_DIR``)."""
    _plots.set_style(FIG_DIR, _EXTRA_RC)


def goal_heading(x, y):
    d = y - x
    return np.arctan2(d[1], d[0])


def vonmises_pdf(theta, mu, kappa):
    return np.exp(kappa * np.cos(theta - mu)) / (2.0 * np.pi * np.i0(kappa))


def _clean_axis(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def plot_vonmises_policy(name="vonmises_navigation"):
    """Policy densities (top row) + stochastic trajectories (bottom row).

    Saved as a 300 dpi vector PDF and an SVG under
    ``figures/von_mises_policy/`` and returned.
    """
    # Figure: top row PDFs, bottom row trajectories
    n = len(KAPPAS)
    fig, axes = plt.subplots(
        2, n,
        figsize=(1.10 * n, 2.20),   # more compact overall
        squeeze=False,
    )

    theta_grid = np.linspace(0, 2*np.pi, 1024, endpoint=False)
    mu0 = goal_heading(X_START, TARGET)

    # --- Make top row compact: bound the view by the policy radius plus the
    #     modulation amplitude, NOT the literal target at y=10 (which
    #     creates huge whitespace).
    Rmax = R_POLICY + AMP

    # -------- Row 0: exact von Mises PDF rendered as a radial modulation --------
    for i, kappa in enumerate(KAPPAS):
        ax = axes[0, i]

        pdf = vonmises_pdf(theta_grid, mu=mu0, kappa=float(kappa)).astype(np.float32)
        pdf_norm = (pdf - pdf.min()) / (pdf.max() - pdf.min() + 1e-12)

        rad = R_POLICY + AMP * pdf_norm
        xs = X_START[0] + rad * np.cos(theta_grid)
        ys = X_START[1] + rad * np.sin(theta_grid)

        ax.plot(xs, ys, color="k", lw=1.1, zorder=3)
        ax.add_patch(Circle(xy=X_START, radius=R_POLICY, fill=False, lw=0.8, ec="k", alpha=0.95))

        ax.plot(X_START[0], X_START[1], "o", ms=2.2, c="k", zorder=5)

        ax.set_title(rf"$\kappa={float(kappa):g}$", pad=2)

        ax.set_aspect("equal", adjustable="box")
        _clean_axis(ax)

        pad = 0.30
        ax.set_xlim(X_START[0] - Rmax - pad, X_START[0] + Rmax + pad)
        ax.set_ylim(X_START[1] - Rmax - pad, X_START[1] + Rmax + pad)

    # -------- Row 1: trajectories by stepwise von Mises sampling --------
    label_idx = 1

    for i, kappa in enumerate(KAPPAS):
        ax = axes[1, i]

        rng = np.random.default_rng(RNG_SEED + i)
        M = NUM_TRAJ
        pos = np.empty((T_STEPS + 1, M, 2), dtype=np.float32)
        pos[0, :, :] = X_START[None, :]

        for t in range(T_STEPS):
            x_t = pos[t]
            mu_t = np.arctan2(TARGET[None, 1] - x_t[:, 1], TARGET[None, 0] - x_t[:, 0])
            th_t = rng.vonmises(mu=mu_t, kappa=float(kappa), size=M).astype(np.float32)

            step = np.stack([np.cos(th_t), np.sin(th_t)], axis=-1) * STEP_LEN
            pos[t + 1] = x_t + step

        # mpl.colormaps is the current API (matplotlib.cm.get_cmap is
        # deprecated); same tab10/tab20 colormaps.
        cmap = mpl.colormaps["tab10"] if M <= 10 else mpl.colormaps["tab20"]
        for j in range(M):
            ax.plot(pos[:, j, 0], pos[:, j, 1], lw=0.55, alpha=0.95, color=cmap(j), zorder=2)

        ax.plot(X_START[0], X_START[1], "o", ms=3.4, c="k", zorder=5)
        ax.plot(TARGET[0],  TARGET[1],  "^", ms=4.0, c="k", zorder=5)

        if i == label_idx:
            ax.annotate(r"Start $\mathbf{x}$",
                        xy=(X_START[0], X_START[1]),
                        xytext=(8, 0), textcoords="offset points",
                        ha="left", va="center")
            ax.annotate(r"Target $\mathbf{g}$",
                        xy=(TARGET[0], TARGET[1]),
                        xytext=(8, 0), textcoords="offset points",
                        ha="left", va="center")

        ax.set_aspect("equal", adjustable="box")
        _clean_axis(ax)

        # slightly tighter than before
        ax.set_xlim(-7.5, 7.5)
        ax.set_ylim(-3.5, 10.8)

    # Row labels (Arial)
    axes[0, 0].set_ylabel("Policy density", labelpad=4)
    axes[1, 0].set_ylabel("Trajectories",  labelpad=4)

    # Compact spacing between panels/rows
    fig.subplots_adjust(left=0.10, right=0.995, bottom=0.06, top=0.93, wspace=0.22, hspace=0.10)

    _plots.save_figure(fig, FIG_DIR, name, pdf=True, svg=True)
    return fig
