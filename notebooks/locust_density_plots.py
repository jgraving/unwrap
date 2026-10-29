"""Plots for the locust density/order notebook (``locust_density_order.ipynb``).

The notebook keeps the statistical code (loading, model fitting, per-cell
summaries). This module holds the matplotlib drawing. Import it from the
notebook's header cell (which also inserts the notebook dir on ``sys.path``)::

    import locust_density_plots as plots
    plots.set_style()

and call the ``plot_*`` functions in place of the old inline drawing code.
Each function mirrors the drawing block of its source cell and is passed only
the pre-computed values it needs.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.legend_handler import HandlerBase

FIG_DIR = Path(__file__).resolve().parent.parent / "figures" / "locust"
DPI = 300
FONT_SZ = 8
AX_LW = 1.0
TICK_W = 1.0


def set_style():
    """One figure style for every figure in this notebook."""
    mpl.rcParams.update({
        "font.family": "Arial",
        "font.size": FONT_SZ,
        "figure.dpi": DPI,
        "savefig.dpi": DPI,
        "axes.linewidth": AX_LW,
        "xtick.major.width": TICK_W,
        "ytick.major.width": TICK_W,
        "lines.linewidth": 1.0,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    print(f"Figures will be saved to: {FIG_DIR.relative_to(FIG_DIR.parent.parent)}")


def _ensure_fig_dir():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    return FIG_DIR


def _save(fig, name, pdf=True, png=True, transparent=False):
    """Save to FIG_DIR as .pdf and/or .png. Returns the list of paths written."""
    out = _ensure_fig_dir()
    written = []
    if pdf:
        p = out / (name + ".pdf")
        p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(p, dpi=DPI, bbox_inches="tight", pad_inches=0.01,
                    transparent=transparent)
        written.append(str(p))
    if png:
        p = out / (name + ".png")
        fig.savefig(p, dpi=DPI)
        written.append(str(p))
    if written:
        print("Saved figure" + ("s" if len(written) > 1 else "") + " to: " +
              ("; ".join(written) if len(written) > 1 else written[0]))
    return written


# ---------------------------------------------------------------------------
# Shared helpers reused by several figure cells in this notebook.
# ---------------------------------------------------------------------------
def _sorted_unique(series):
    vals = series.dropna().unique().tolist()
    try:
        return sorted(vals, key=float)
    except Exception:
        return sorted(vals)


def _lerp(a, b, t):
    return (1.0 - t) * np.asarray(a, float) + t * np.asarray(b, float)


def cell_color(order_idx, density_idx, n_rows, n_cols,
               tmin=0.10, tmax=0.95,
               blue=(0.10, 0.25, 1.00),
               red=(1.00, 0.15, 0.15),
               gamma=1.15,
               mix="max"):
    """Deep-purple 2D colormap: black->red (order, top->bottom) + black->blue
    (density, left->right). Returns an RGB tuple."""
    to = order_idx / max(n_rows - 1, 1)
    td = density_idx / max(n_cols - 1, 1)
    to = np.power(to, gamma)
    td = np.power(td, gamma)
    to = tmin + (tmax - tmin) * to
    td = tmin + (tmax - tmin) * td
    black = np.array([0.0, 0.0, 0.0], float)
    c_order = _lerp(black, red, to)   # black -> red
    c_dens = _lerp(black, blue, td)   # black -> blue
    if mix == "add":
        c = c_order + c_dens
    else:  # "max": channelwise max keeps the field deep
        c = np.maximum(c_order, c_dens)
    return tuple(np.clip(c, 0, 1))


def nice_ceil(v):
    if v <= 0:
        return 1.0
    v = v * 1.10
    exp = np.floor(np.log10(v))
    frac = v / (10 ** exp)
    for nice in (1.0, 2.0, 2.5, 5.0, 10.0):
        if frac <= nice:
            return nice * (10 ** exp)
    return 10.0 * (10 ** exp)


def kappa_from_order(order, eps=1e-6):
    """Map 'order' (mean resultant length R) -> kappa."""
    from unwrap.bessel import inverse_bessel_ratio
    R = float(np.clip(order, 0.0, 1.0 - eps))
    return inverse_bessel_ratio(R)


def _color_key_axes(fig, n_rows, n_cols, x=0.88, size=0.12, key="cell_color"):
    """Draw the shared square 2D order/density color key off to the right of
    ``fig``, matching the ``cell_color`` panel palette (row 0 on top)."""
    key_ax = fig.add_axes([x, 0.5 - size / 2, size, size])
    grid = np.zeros((n_rows, n_cols, 3), float)
    for rr in range(n_rows):
        for cc in range(n_cols):
            grid[rr, cc, :] = cell_color(rr, cc, n_rows, n_cols)
    key_ax.imshow(grid, origin="upper", aspect="equal")
    key_ax.set_xlabel("Density", labelpad=2)
    key_ax.set_ylabel("Order", labelpad=2)
    key_ax.set_xticks([0, n_cols - 1] if n_cols > 1 else [0])
    key_ax.set_yticks([0, n_rows - 1] if n_rows > 1 else [0])
    key_ax.set_xticklabels(["low", "high"] if n_cols > 1 else [""])
    key_ax.set_yticklabels(["low", "high"] if n_rows > 1 else [""])
    for s in key_ax.spines.values():
        s.set_linewidth(1)
    return key_ax


# ---------------------------------------------------------------------------
# Cell 12: trajectories per animal, one panel per (order, density) condition,
# colored by the shared 2D purple map; 0° stimulus line in each panel.
# ---------------------------------------------------------------------------
def plot_trajectory_grid(merged_df, orders, densities):
    XLIM = (-8, 8)
    YLIM = (-8, 8)
    n_rows = len(orders)
    n_cols = len(densities)

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(n_cols * 1.4 + 1.0, n_rows * 1.4),
        dpi=DPI, squeeze=False,
    )

    for r, o in enumerate(orders):
        for c, d in enumerate(densities):
            ax = axes[r, c]
            df_sub = merged_df.loc[(merged_df["order"] == o) & (merged_df["density"] == d)]
            col = cell_color(r, c, n_rows, n_cols)

            ax.set_xlim(*XLIM)
            ax.set_ylim(*YLIM)
            ax.set_xticks([-8, 0, 8])
            ax.set_yticks([-8, 0, 8])
            ax.set_aspect("equal")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.set_title(f"Density = {d}\nOrder = {o}", fontsize=8, pad=2)

            if df_sub.empty:
                ax.text(0.5, 0.5, "no data", ha="center", va="center",
                        transform=ax.transAxes, alpha=0.6)
                ax.plot([0, XLIM[1]], [0, 0], color="black", linewidth=1.0, zorder=10)
                continue

            for _, grp in df_sub.groupby("fname", sort=False):
                ax.plot(grp["X"].values, grp["Y"].values,
                        color=col, linewidth=0.8, alpha=0.75, zorder=2)

            ax.plot([0, XLIM[1]], [0, 0], color="black", linewidth=1.0, zorder=10)

    fig.tight_layout(rect=[0.0, 0.0, 0.86, 1.0])
    _color_key_axes(fig, n_rows, n_cols, x=0.88, size=0.12)

    _save(fig, "trajectories_order_by_density_grid")
    plt.show()


# ---------------------------------------------------------------------------
# Cell 13: mean displacement vector per individual, one panel per
# (order, density) condition, in step-length units (r = 1 = one step).
# ---------------------------------------------------------------------------
def plot_mean_vector_grid(means_df, step_length):
    SCALE = 1 / step_length

    orders = _sorted_unique(means_df["order"])
    densities = _sorted_unique(means_df["density"])
    n_rows = len(orders)
    n_cols = len(densities)

    max_abs_component = float(np.max(np.abs(means_df[["dX", "dY"]].to_numpy())))
    L = nice_ceil(max_abs_component * SCALE)
    L = max(L, 1.15)
    XLIM = (-L, L)
    YLIM = (-L, L)

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(n_cols * 1.25 + 1.0, n_rows * 1.25),
        dpi=DPI, squeeze=False,
    )

    theta = np.linspace(0, 2 * np.pi, 256)
    circ_x = np.cos(theta)
    circ_y = np.sin(theta)

    for r, o in enumerate(orders):
        for c, d in enumerate(densities):
            ax = axes[r, c]
            panel = means_df[(means_df["order"] == o) & (means_df["density"] == d)]
            col = cell_color(r, c, n_rows, n_cols)

            ax.set_xlim(*XLIM)
            ax.set_ylim(*YLIM)
            ax.set_aspect("equal", "box")
            ax.set_xticks([])
            ax.set_yticks([])
            for s in ax.spines.values():
                s.set_visible(False)
            ax.set_title(f"Density = {d}\nOrder = {o}", fontsize=8, pad=2)
            ax.plot(circ_x, circ_y, color="black", linewidth=0.8, zorder=1)

            if panel.empty:
                ax.text(0.5, 0.5, "no data", ha="center", va="center",
                        transform=ax.transAxes, alpha=0.6)
                ax.plot([0, 1], [0, 0], color="black", linewidth=0.9, zorder=10)
                continue

            for row in panel.itertuples(index=False):
                U = float(row.dX) * SCALE
                V = float(row.dY) * SCALE
                ax.plot([0.0, U], [0.0, V],
                        linewidth=1.1, alpha=0.80, color=col, zorder=3)

            ax.plot(0, 0, ".", ms=2, color="k", zorder=4)
            ax.plot([0, 1], [0, 0], color="black", linewidth=0.9, zorder=10)
            ax.plot([0.96, 1.02], [0.0, 0.0], color="black", lw=1.0, zorder=10,
                    solid_capstyle="round")

    fig.tight_layout(rect=[0.0, 0.0, 0.86, 1.0])
    _color_key_axes(fig, n_rows, n_cols, x=0.88, size=0.12)

    _save(fig, "mean_vectors_order_by_density_grid")
    plt.show()


# ---------------------------------------------------------------------------
# Cell 14: synthetic schematic — per row, DENSITY dots (left, Reds) and
# ORDER trajectories (right, Purples). Self-contained illustration.
# ---------------------------------------------------------------------------
def plot_order_density_trajectories(
        order_density_pairs=None, n_traj=12, n_steps=100, step_L=0.04,
        cmap_min=0.45, cmap_max=0.95, seed=123,
        fname="order_density_trajectories"):
    mpl.rcParams.update({
        "font.size": 7,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "font.family": "Arial",
    })
    if order_density_pairs is None:
        # largest at the top
        orders = np.array([1.0, 0.8, 0.6, 0.4, 0.2, 0.0])
        densities = np.array([64, 32, 16, 8, 4, 1])
    else:
        orders, densities = order_density_pairs

    n_rows = len(orders)
    n_cols = 2

    rng = np.random.default_rng(seed)
    cmap_density = plt.cm.Reds
    cmap_order = plt.cm.Purples

    dens_levels = np.sort(np.unique(densities))
    order_levels = np.sort(np.unique(orders))
    dens_pos = np.linspace(cmap_min, cmap_max, len(dens_levels))
    order_pos = np.linspace(cmap_min, cmap_max, len(order_levels))
    dens_color_map = {d: cmap_density(dens_pos[i]) for i, d in enumerate(dens_levels)}
    order_color_map = {o: cmap_order(order_pos[i]) for i, o in enumerate(order_levels)}

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(1.0, 3.0),
        dpi=DPI, squeeze=False,
    )

    XLIM = (0.0, 1.0)
    YLIM = (0.0, 1.0)

    for row, (order, dens) in enumerate(zip(orders, densities)):
        # LEFT: density dots (Reds)
        ax_dens = axes[row, 0]
        n = int(dens)
        dot_col = dens_color_map[dens]
        if n == 1:
            xs = np.array([0.5]); ys = np.array([0.5])
        elif n == 4:
            xs, ys = np.meshgrid([0.25, 0.75], [0.25, 0.75])
            xs = xs.ravel(); ys = ys.ravel()
        else:
            n_side = int(np.ceil(np.sqrt(n)))
            xs, ys = np.meshgrid(np.linspace(0.1, 0.9, n_side),
                                 np.linspace(0.1, 0.9, n_side))
            xs = xs.ravel()[:n]; ys = ys.ravel()[:n]
        jitter_scale = 0.025
        xs = np.clip(xs + rng.normal(0, jitter_scale, size=xs.size), 0.075, 0.925)
        ys = np.clip(ys + rng.normal(0, jitter_scale, size=ys.size), 0.075, 0.925)

        ax_dens.scatter(xs, ys, s=1, color=dot_col)
        ax_dens.set_xlim(*XLIM)
        ax_dens.set_ylim(*YLIM)
        ax_dens.set_aspect("equal", "box")
        ax_dens.set_xticks([])
        ax_dens.set_yticks([])
        ax_dens.text(-0.08, 0.5, f"{dens:.0f}", ha="right", va="center",
                     transform=ax_dens.transAxes)

        # RIGHT: trajectories (Purples)
        ax_traj = axes[row, 1]
        kappa = kappa_from_order(order)
        traj_col = order_color_map[order]
        if order == 0.0:
            x0 = rng.uniform(0.1, 0.9, size=n_traj)
            y0 = rng.uniform(0.1, 0.9, size=n_traj)
        else:
            x0 = np.full(n_traj, 0.05)
            y0 = np.linspace(0.15, 0.85, n_traj)

        for i in range(n_traj):
            if order == 1.0:
                angles = np.zeros(n_steps)
            else:
                angles = rng.vonmises(mu=0.0, kappa=kappa, size=n_steps)
            dx = step_L * np.cos(angles)
            dy = step_L * np.sin(angles)
            xs = np.concatenate([[x0[i]], x0[i] + np.cumsum(dx)])
            ys = np.concatenate([[y0[i]], y0[i] + np.cumsum(dy)])
            ax_traj.plot(xs, ys, color=traj_col, linewidth=0.7)

        ax_traj.set_xlim(*XLIM)
        ax_traj.set_ylim(*YLIM)
        ax_traj.set_aspect("equal", "box")
        ax_traj.set_xticks([])
        ax_traj.set_yticks([])
        ax_traj.text(1.08, 0.5, f"{order:.1f}", ha="left", va="center",
                     transform=ax_traj.transAxes)

    _save(fig, fname)
    plt.show()


# ---------------------------------------------------------------------------
# Cell 15: synthetic 6x6 grid — n_traj == density trajectories per panel,
# colored by the shared 2D purple map; 2D key on the right.
# ---------------------------------------------------------------------------
def plot_order_density_pairwise_trajectories(
        orders=None, densities=None, n_steps=60, step_L=0.04,
        seed=123, fname="order_density_pairwise_trajectories_only"):
    mpl.rcParams.update({
        "font.size": 7,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "font.family": "Arial",
    })
    if orders is None:
        orders = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    if densities is None:
        densities = np.array([1, 4, 8, 16, 32, 64])

    n_rows = len(orders)
    n_cols = len(densities)

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(6.0 + 1.0, 6.0),
        dpi=DPI, squeeze=False,
    )
    XLIM = (0.0, 1.0)
    YLIM = (0.0, 1.0)

    for r, order in enumerate(orders):
        for c, dens in enumerate(densities):
            ax = axes[r, c]
            cell_rng = np.random.default_rng(seed + 1000 * r + 10 * c)
            n_traj = int(dens)
            kappa = kappa_from_order(order)
            col = cell_color(r, c, n_rows, n_cols)

            if order == 0.0:
                x0 = cell_rng.uniform(0.10, 0.90, size=n_traj)
                y0 = cell_rng.uniform(0.10, 0.90, size=n_traj)
            else:
                x0 = np.full(n_traj, 0.0)
                if n_traj == 1:
                    y0 = np.array([0.5])
                else:
                    y0 = np.linspace(0.05, 0.95, n_traj)

            lw = 0.30 if n_traj >= 16 else 0.55
            for i in range(n_traj):
                if order == 1.0:
                    angles = np.zeros(n_steps)
                else:
                    angles = cell_rng.vonmises(mu=0.0, kappa=kappa, size=n_steps)
                dx = step_L * np.cos(angles)
                dy = step_L * np.sin(angles)
                xs = np.concatenate([[x0[i]], x0[i] + np.cumsum(dx)])
                ys = np.concatenate([[y0[i]], y0[i] + np.cumsum(dy)])
                ax.plot(xs, ys, color=col, linewidth=lw, alpha=0.85)

            ax.set_xlim(*XLIM)
            ax.set_ylim(*YLIM)
            ax.set_aspect("equal", "box")
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_color("black")
            if c == 0:
                ax.text(-0.20, 0.5, f"{order:.1f}", ha="right", va="center",
                        transform=ax.transAxes)
            if r == 0:
                ax.set_title(f"{dens:.0f}", pad=3)

    fig.text(0.5, 0.95, "Density (m$^{-2}$)", ha="center", va="top")
    fig.text(0.07, 0.5, "Order", ha="left", va="center", rotation=90)
    _color_key_axes(fig, n_rows, n_cols, x=0.95, size=0.13)

    _save(fig, fname)
    plt.show()


# ---------------------------------------------------------------------------
# Cell 21: two-panel predictive alignment = A(kappa)*cos(mu).
#   Panel i: vs Order (colored by density),  Panel ii: vs Density (colored by order).
# The notebook computes pred_i / pred_ii ({level: {mean, low, high}}) first.
# ---------------------------------------------------------------------------
def plot_alignment_prediction(order_grid, density_grid, dens_levels, order_levels,
                              pred_i, pred_ii):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5, 2), dpi=DPI)

    cmap_density = plt.cm.Reds
    cmap_order = plt.cm.Purples

    colors_d = np.linspace(0.95, 0.3, len(dens_levels))
    h_i, l_i = [], []
    for i, d in enumerate(dens_levels):
        col = cmap_density(colors_d[i])
        y = pred_i[d]["mean"]; lo = pred_i[d]["low"]; hi = pred_i[d]["high"]
        line, = ax1.plot(order_grid, y, color=col, lw=1.5)
        ax1.fill_between(order_grid, lo, hi, color=col, alpha=0.25, linewidth=0)
        h_i.append(line); l_i.append(f"{d:g}")
    ax1.set_xlabel("Order")
    ax1.set_ylabel("Directional alignment")
    ax1.set_xlim(order_grid.min(), order_grid.max())
    ax1.set_ylim(-0.05, 1.05)
    ax1.set_yticks([0.0, 0.5, 1.0])
    ax1.legend(h_i, l_i, title="Density", loc="center left",
               bbox_to_anchor=(1.02, 0.5), frameon=False)

    colors_o = np.linspace(0.95, 0.3, len(order_levels))
    h_ii, l_ii = [], []
    for i, o in enumerate(order_levels):
        col = cmap_order(colors_o[i])
        y = pred_ii[o]["mean"]; lo = pred_ii[o]["low"]; hi = pred_ii[o]["high"]
        line, = ax2.plot(density_grid, y, color=col, lw=1.5)
        ax2.fill_between(density_grid, lo, hi, color=col, alpha=0.25, linewidth=0)
        h_ii.append(line); l_ii.append(f"{o:.1f}" if (o % 1) else f"{int(o)}")
    ax2.set_xscale("log", base=2)
    ax2.set_xlabel("Density")
    ax2.set_xlim(density_grid.min(), density_grid.max())
    ax2.set_xticks([1, 2, 4, 8, 16, 32, 64])
    ax2.set_xticklabels([1, 2, 4, 8, 16, 32, 64])
    ax2.set_ylim(-0.05, 1.05)
    ax2.set_yticks([0.0, 0.5, 1.0])
    ax2.legend(h_ii, l_ii, title="Order", loc="center left",
               bbox_to_anchor=(1.02, 0.5), frameon=False)

    fig.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# Cell 25: model selection -- delta ELPD for the matched model sweep, split
# into "order included" vs "order excluded" models.
# ---------------------------------------------------------------------------
def plot_model_selection(res_df, use_se):
    import textwrap

    def _wrap_label(s, width=40):
        return textwrap.fill(str(s), width=width)

    y = np.arange(len(res_df))
    y_labels = [_wrap_label(n, 40) for n in res_df["pretty_name"]]

    g_incl = res_df[res_df["has_order"]].copy()
    g_excl = res_df[~res_df["has_order"]].copy()
    g_incl["y"] = g_incl.index.to_numpy()
    g_excl["y"] = g_excl.index.to_numpy()

    fig = plt.figure(figsize=(1.8, 0.33 * len(res_df)), dpi=DPI)
    if use_se:
        fig.gca().errorbar(g_incl["delta_elpd"], g_incl["y"], xerr=g_incl["se_delta_elpd"],
                           fmt="o", capsize=3, elinewidth=1, markersize=5,
                           label="Order included")
        fig.gca().errorbar(g_excl["delta_elpd"], g_excl["y"], xerr=g_excl["se_delta_elpd"],
                           fmt="s", capsize=3, elinewidth=1, markersize=5,
                           label="Order excluded")
    else:
        fig.gca().plot(g_incl["delta_elpd"], g_incl["y"], "o", label="Order included")
        fig.gca().plot(g_excl["delta_elpd"], g_excl["y"], "s", label="Order excluded")

    fig.gca().axvline(0.0, linestyle="--", linewidth=1, color="k", zorder=0)
    fig.gca().set_yticks(y, y_labels)
    fig.gca().invert_yaxis()
    fig.gca().set_xlabel("ELPD difference")
    fig.legend(loc="upper left", frameon=True)

    _save(fig, "order_density_model_selection")
    plt.show()


def _style_axis(ax, ax_lw=0.8):
    for s in ax.spines.values():
        s.set_linewidth(ax_lw)


# ---------------------------------------------------------------------------
# Cell 43: two-panel predictive alignment (hierarchical, centered predictors).
#   Panel i: vs Order colored by Density; Panel ii: vs Density colored by Order.
# The notebook computes pred_i / pred_ii first (mean / low / high per level).
# ---------------------------------------------------------------------------
def plot_alignment_prediction_centered(order_grid, density_grid, dens_levels,
                                       order_levels, pred_i, pred_ii):
    mpl.rcParams.update({
        "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7,
        "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
        "axes.linewidth": 0.8, "lines.linewidth": 1.0,
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    })
    fig_w = 1.8 * 2 + 0.6
    fig_h = 1.55
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(fig_w, fig_h), dpi=DPI)

    cmap_density = plt.cm.Reds
    cmap_order = plt.cm.Purples
    fill_alpha = 0.35
    cmap_min, cmap_max = 0.45, 0.95

    colors_d = np.linspace(cmap_max, cmap_min, len(dens_levels))
    h_i, l_i = [], []
    for i, d in enumerate(dens_levels):
        col = cmap_density(colors_d[i])
        y = pred_i[d]["mean"]; lo = pred_i[d]["low"]; hi = pred_i[d]["high"]
        ax1.fill_between(order_grid, lo, hi, color=col, alpha=fill_alpha,
                         edgecolor="none", zorder=1)
        line, = ax1.plot(order_grid, y, color=col, lw=1.0, zorder=2)
        h_i.append(line); l_i.append(f"{d:g}")
    ax1.set_xlabel("Order")
    ax1.set_ylabel("Alignment")
    ax1.set_xlim(order_grid.min(), order_grid.max())
    ax1.set_ylim(-0.0, 1.0)
    ax1.set_yticks([0.0, 0.5, 1.0])
    ax1.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax1.set_xticklabels(["0", "0.2", "0.4", "0.6", "0.8", "1.0"])
    _style_axis(ax1)
    ax1.legend(h_i, l_i, title="Density", loc="center left",
               bbox_to_anchor=(1.02, 0.5), frameon=True, handlelength=1.2,
               borderpad=0.3, labelspacing=0.25, handletextpad=0.4)

    colors_o = np.linspace(cmap_max, cmap_min, len(order_levels))
    h_ii, l_ii = [], []
    for i, o in enumerate(order_levels):
        col = cmap_order(colors_o[i])
        y = pred_ii[o]["mean"]; lo = pred_ii[o]["low"]; hi = pred_ii[o]["high"]
        ax2.fill_between(density_grid, lo, hi, color=col, alpha=fill_alpha,
                         edgecolor="none", zorder=1)
        line, = ax2.plot(density_grid, y, color=col, lw=1.0, zorder=2)
        h_ii.append(line)
        l_ii.append(f"{o:.1f}" if (o % 1) else f"{int(o)}")
    ax2.set_xscale("log", base=2)
    ax2.set_xlabel("Density")
    ax2.set_xlim(density_grid.min(), density_grid.max())
    ax2.set_xticks([1, 2, 4, 8, 16, 32, 64])
    ax2.set_xticklabels([1, 2, 4, 8, 16, 32, 64])
    ax2.set_yticks([0.0, 0.5, 1.0])
    ax2.set_ylim(-0.0, 1.0)
    _style_axis(ax2)
    ax2.legend(h_ii, l_ii, title="Order", loc="center left",
               bbox_to_anchor=(1.02, 0.5), frameon=True, handlelength=1.2,
               borderpad=0.3, labelspacing=0.25, handletextpad=0.4)

    fig.tight_layout()
    _save(fig, "order_density_alignment")
    plt.show()


# ---------------------------------------------------------------------------
# Cell 44: two-panel posterior predictive density at y=0 (computed in log
# space in the notebook, exponentiated to density here). Uniform reference
# line at 1/(2*pi). Transparent axes.
# ---------------------------------------------------------------------------
def plot_elpd0_prediction_centered(order_grid, density_grid, dens_levels,
                                   order_levels, pred_i, pred_ii):
    mpl.rcParams.update({
        "font.family": "Arial",
        "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7,
        "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
        "axes.linewidth": 0.8, "lines.linewidth": 1.0,
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    })
    fig_w = 1.8 * 2 + 0.6
    fig_h = 1.55
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(fig_w, fig_h), dpi=DPI)

    fig.patch.set_facecolor("none")
    fig.patch.set_alpha(0)
    for ax in (ax1, ax2):
        ax.set_facecolor("none")
        ax.patch.set_facecolor("none")
        ax.patch.set_alpha(0)

    cmap_density = plt.cm.Reds
    cmap_order = plt.cm.Purples
    fill_alpha = 0.35
    cmap_min, cmap_max = 0.45, 0.95
    p_uniform = 1.0 / (2.0 * np.pi)

    colors_d = np.linspace(cmap_max, cmap_min, len(dens_levels))
    h_i, l_i = [], []
    for i, d in enumerate(dens_levels):
        col = cmap_density(colors_d[i])
        y = pred_i[d]["mean"]; lo = pred_i[d]["low"]; hi = pred_i[d]["high"]
        ax1.fill_between(order_grid, lo, hi, color=col, alpha=fill_alpha,
                         edgecolor="none", zorder=1)
        line, = ax1.plot(order_grid, y, color=col, lw=1.0, zorder=2)
        h_i.append(line); l_i.append(f"{d:g}")
    ax1.set_xlabel("Order")
    ax1.set_ylabel("Predictive density")
    ax1.set_xlim(order_grid.min(), order_grid.max())
    ax1.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax1.set_xticklabels(["0", "0.2", "0.4", "0.6", "0.8", "1.0"])
    ax1.axhline(p_uniform, ls="--", color="k", lw=0.9, zorder=0)
    _style_axis(ax1)
    ax1.legend(h_i, l_i, title="Density", loc="center left",
               bbox_to_anchor=(1.02, 0.5), frameon=True, handlelength=1.2,
               borderpad=0.3, labelspacing=0.25, handletextpad=0.4)

    colors_o = np.linspace(cmap_max, cmap_min, len(order_levels))
    h_ii, l_ii = [], []
    for i, o in enumerate(order_levels):
        col = cmap_order(colors_o[i])
        y = pred_ii[o]["mean"]; lo = pred_ii[o]["low"]; hi = pred_ii[o]["high"]
        ax2.fill_between(density_grid, lo, hi, color=col, alpha=fill_alpha,
                         edgecolor="none", zorder=1)
        line, = ax2.plot(density_grid, y, color=col, lw=1.0, zorder=2)
        h_ii.append(line)
        l_ii.append(f"{o:.1f}" if (o % 1) else f"{int(o)}")
    ax2.set_xscale("log", base=2)
    ax2.set_xlabel("Density")
    ax2.set_xlim(density_grid.min(), density_grid.max())
    ax2.set_xticks([1, 2, 4, 8, 16, 32, 64])
    ax2.set_xticklabels([1, 2, 4, 8, 16, 32, 64])
    ax2.axhline(p_uniform, ls="--", color="k", lw=0.9, zorder=0)
    _style_axis(ax2)
    ax2.legend(h_ii, l_ii, title="Order", loc="center left",
               bbox_to_anchor=(1.02, 0.5), frameon=True, handlelength=1.2,
               borderpad=0.3, labelspacing=0.25, handletextpad=0.4)

    ax1.set_ylim(0.0, 1.0)
    ax2.set_ylim(0.0, 1.0)

    fig.tight_layout()
    _save(fig, "order_density_elpd0", transparent=True)
    plt.show()


# ---------------------------------------------------------------------------
# Cell 45: two rows of polar small-multiples, multiple von Mises predictive
# PDFs per panel. Row 1: columns = order, curves = density (Reds);
# Row 2: columns = density, curves = order (Purples). The notebook computes
# ``summ`` {(o, d): {med, lo, hi, ...}} and ``peak_max`` first.
# ---------------------------------------------------------------------------
def plot_polar_pdfs_grid(theta, summ, order_levels, dens_levels, peak_max,
                         hdi_prob=0.90, fname="order_density_polar_pdfs_grid"):
    from matplotlib.patches import Circle
    from matplotlib.lines import Line2D

    mpl.rcParams.update({
        "font.size": 5, "axes.labelsize": 5, "axes.titlesize": 5,
        "legend.fontsize": 5, "xtick.labelsize": 5, "ytick.labelsize": 5,
        "font.family": "Arial", "svg.fonttype": "none",
        "axes.linewidth": 0.8, "lines.linewidth": 0.5,
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    })

    EDGE_RADIUS = 1.0
    AMP = 1.0
    OUTER_LW = 0.6
    cmap_density = plt.cm.Reds
    cmap_order = plt.cm.Purples
    cmap_min, cmap_max = 0.45, 0.95

    RMAX = EDGE_RADIUS + AMP * 1.03

    dens_pos = np.linspace(cmap_max, cmap_min, len(dens_levels))
    order_pos = np.linspace(cmap_max, cmap_min, len(order_levels))
    dens_colors = {d: cmap_density(dens_pos[i]) for i, d in enumerate(dens_levels)}
    order_colors = {o: cmap_order(order_pos[i]) for i, o in enumerate(order_levels)}

    def to_r(y):
        return EDGE_RADIUS + AMP * (y / (peak_max + 1e-12))

    def format_polar(ax):
        ax.set_theta_zero_location("E")
        ax.set_theta_direction(1)
        ax.spines["polar"].set_visible(False)
        ax.set_xticks([])
        ax.set_xticklabels([])
        ax.grid(False)
        ax.set_yticks([])
        ax.set_yticklabels([])
        ax.plot([0, 0], [EDGE_RADIUS, RMAX], color="0.80", lw=0.9, ls="-", zorder=0)
        ax.add_artist(Circle((0, 0), EDGE_RADIUS, transform=ax.transData._b,
                             color="0.75", fill=False, linewidth=OUTER_LW))
        ax.add_artist(Circle((0, 0), EDGE_RADIUS + AMP, transform=ax.transData._b,
                             color="0.85", fill=False, linewidth=OUTER_LW))

    n1 = len(order_levels)
    n2 = len(dens_levels)
    ncol = max(n1, n2)
    panel_w = 0.8
    fig_h = 2.0
    fig, axes = plt.subplots(
        2, ncol, figsize=(panel_w * ncol, fig_h), dpi=DPI,
        subplot_kw={"projection": "polar"}, squeeze=True,
    )
    axes = np.atleast_2d(axes)
    for j in range(ncol):
        if j >= n1:
            axes[0, j].axis("off")
        if j >= n2:
            axes[1, j].axis("off")

    # Row 1: panels = order; curves = density
    for j, o in enumerate(order_levels[::-1]):
        ax = axes[0, j]
        format_polar(ax)
        ax.set_ylim(0, RMAX)
        for i, d in enumerate(dens_levels):
            col = dens_colors[d]
            res = summ[(float(o), float(d))]
            ax.fill_between(theta, to_r(res["lo"]), to_r(res["hi"]),
                            facecolor=col, edgecolor="none", zorder=1, alpha=0.5)
        for i, d in enumerate(dens_levels):
            col = dens_colors[d]
            res = summ[(float(o), float(d))]
            ax.plot(theta, to_r(res["med"]), color=col, lw=0.5, zorder=3)
        ax.text(0.5, 0.5, f"Order\n{float(o):g}", ha="center", va="center",
                fontsize=5, color="k", transform=ax.transAxes, zorder=10)

    # Row 2: panels = density; curves = order
    for j, d in enumerate(dens_levels[::-1]):
        ax = axes[1, j]
        format_polar(ax)
        ax.set_ylim(0, RMAX)
        for i, o in enumerate(order_levels):
            col = order_colors[o]
            res = summ[(float(o), float(d))]
            ax.fill_between(theta, to_r(res["lo"]), to_r(res["hi"]),
                            facecolor=col, edgecolor="none", zorder=1, alpha=0.5)
        for i, o in enumerate(order_levels):
            col = order_colors[o]
            res = summ[(float(o), float(d))]
            ax.plot(theta, to_r(res["med"]), color=col, lw=0.5, zorder=3)
        ax.text(0.5, 0.5, f"Density\n{float(d):g}", ha="center", va="center",
                fontsize=5, color="k", transform=ax.transAxes, zorder=10)

    dens_handles = [Line2D([0], [0], color=dens_colors[d], lw=1.6) for d in dens_levels]
    dens_labels = [f"{d:g}" for d in dens_levels]
    order_handles = [Line2D([0], [0], color=order_colors[o], lw=1.6) for o in order_levels]
    order_labels = [f"{o:g}" for o in order_levels]

    fig.subplots_adjust(left=0.03, right=0.90, top=0.92, bottom=0.06,
                        wspace=0.22, hspace=-0.15)
    fig.legend(dens_handles, dens_labels, title="Density", loc="center",
               bbox_to_anchor=(0.95, 0.69), frameon=True, handlelength=1.2,
               borderpad=0.3, labelspacing=0.25, handletextpad=0.4)
    fig.legend(order_handles, order_labels, title="Order", loc="center",
               bbox_to_anchor=(0.95, 0.29), frameon=True, handlelength=1.2,
               borderpad=0.3, labelspacing=0.25, handletextpad=0.4)

    out = _ensure_fig_dir()
    (out / (fname + ".svg")).unlink(missing_ok=True)
    (out / (fname + ".pdf")).unlink(missing_ok=True)
    (out / (fname + ".png")).unlink(missing_ok=True)
    fig.savefig(out / (fname + ".svg"), format="svg")
    fig.savefig(out / (fname + ".pdf"), bbox_inches="tight")
    fig.savefig(out / (fname + ".png"), dpi=DPI)
    print("Saved figure to:", str(out / (fname + ".pdf")))
    plt.show()


# ---------------------------------------------------------------------------
# Cell 46: coefficient posteriors (fixed effects only), mu (circles) vs
# kappa (triangles) stacked per term, with 95% HDI error bars. The notebook
# computes the per-term means and asymmetric xerr arrays first.
# ---------------------------------------------------------------------------
def plot_coefficient_posteriors(terms_order, mu_mean, mu_xerr, k_mean, k_xerr,
                                fname="order_density_coefficients"):
    # This cell originally inherited the polar-cell style block (font 5).
    mpl.rcParams.update({
        "font.size": 5, "axes.labelsize": 5, "axes.titlesize": 5,
        "legend.fontsize": 5, "xtick.labelsize": 5, "ytick.labelsize": 5,
        "font.family": "Arial", "svg.fonttype": "none",
        "axes.linewidth": 0.8, "lines.linewidth": 0.5,
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    })
    tab10 = mpl.cm.get_cmap("tab10")
    mu_color = tab10(5)   # brown
    k_color = tab10(2)    # green

    y_base = np.arange(len(terms_order), dtype=float)
    dy = 0.16
    y_mu = y_base - dy
    y_k = y_base + dy

    fig = plt.figure(figsize=(1.8, 0.4 * len(terms_order)), dpi=DPI)
    ax = fig.gca()
    ax.errorbar(mu_mean, y_mu, xerr=mu_xerr,
                fmt="D", color=mu_color, ecolor=mu_color,
                capsize=3, elinewidth=1, capthick=1,
                markersize=5, markeredgewidth=0.8, label="μ")
    ax.errorbar(k_mean, y_k, xerr=k_xerr,
                fmt="^", color=k_color, ecolor=k_color,
                capsize=3, elinewidth=1, capthick=1,
                markersize=5, markeredgewidth=0.8, label="κ")
    ax.axvline(0.0, linestyle="--", linewidth=1, color="k", zorder=0)
    ax.set_yticks(y_base, terms_order)
    ax.invert_yaxis()
    ax.set_xlabel("Coefficient value")
    ax.legend(loc="upper right", frameon=True)

    _save(fig, fname)
    plt.show()


# ---------------------------------------------------------------------------
# Cell 47: two-panel figure. LEFT = matched-model sweep ΔELPD,
# RIGHT = coefficient posteriors (mu vs kappa). The notebook computes the
# data (res_df slicing/labels, per-term means + asymmetric xerr) first.
# ---------------------------------------------------------------------------
def plot_model_selection_and_coefficients(
        g_incl, g_excl, y_m, y_labels, use_se,
        terms_order, y_base, y_mu, y_k,
        mu_mean, mu_xerr, k_mean, k_xerr,
        fname="order_density_model_selection_and_coefficients"):
    FONT_SIZE = 7
    AX_LW     = 0.8
    TICK_LW   = 0.8
    TICK_LEN  = 2.5
    LINE_LW   = 1.0
    ERR_LW    = 1.0
    CAPSIZE   = 2.5
    CAPTHICK  = 1.0
    MS        = 4.5
    MEW       = 0.8

    mpl.rcParams.update({
        "font.size": FONT_SIZE, "axes.labelsize": FONT_SIZE,
        "axes.titlesize": FONT_SIZE, "legend.fontsize": FONT_SIZE,
        "xtick.labelsize": FONT_SIZE, "ytick.labelsize": FONT_SIZE,
        "axes.linewidth": AX_LW, "lines.linewidth": LINE_LW,
        "xtick.major.width": TICK_LW, "ytick.major.width": TICK_LW,
        "xtick.major.size": TICK_LEN, "ytick.major.size": TICK_LEN,
        "savefig.dpi": DPI, "figure.dpi": DPI,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })

    def _style_spines(ax):
        for s in ax.spines.values():
            s.set_linewidth(AX_LW)

    def _finish_axis(ax, xlabel, yticks=None, yticklabels=None, legend_loc=None, vline0=True):
        if vline0:
            ax.axvline(0.0, linestyle="--", linewidth=ERR_LW, color="k", zorder=0)
        if yticks is not None:
            ax.set_yticks(yticks)
        if yticklabels is not None:
            ax.set_yticklabels(yticklabels)
        ax.invert_yaxis()
        ax.set_xlabel(xlabel)
        _style_spines(ax)
        if legend_loc is not None:
            leg = ax.legend(loc=legend_loc, frameon=True,
                            handlelength=1.2, borderpad=0.3,
                            labelspacing=0.3, handletextpad=0.4)
            leg.get_frame().set_linewidth(AX_LW)

    tab10 = mpl.cm.get_cmap("tab10")
    ELPD_INCL_COLOR = tab10(0)
    ELPD_EXCL_COLOR = tab10(1)
    MU_COLOR        = tab10(5)
    K_COLOR         = tab10(2)

    ROW_H = 0.36
    # row count mirrors the notebook: max(len(res_df), len(terms_order)); len(y_m) == len(res_df)
    n_rows = max(len(y_m), len(terms_order))
    fig_h = ROW_H * n_rows

    fig, (ax_elpd, ax_coef) = plt.subplots(1, 2, figsize=(1.8 * 2 + 0.6, fig_h), dpi=DPI)

    if use_se:
        ax_elpd.errorbar(
            g_incl["delta_elpd"], g_incl["y"], xerr=g_incl["se_delta_elpd"],
            fmt="o", color=ELPD_INCL_COLOR, ecolor=ELPD_INCL_COLOR,
            ms=MS, mew=MEW, capsize=CAPSIZE, elinewidth=ERR_LW, capthick=CAPTHICK,
            label="Order included")
        ax_elpd.errorbar(
            g_excl["delta_elpd"], g_excl["y"], xerr=g_excl["se_delta_elpd"],
            fmt="s", color=ELPD_EXCL_COLOR, ecolor=ELPD_EXCL_COLOR,
            ms=MS, mew=MEW, capsize=CAPSIZE, elinewidth=ERR_LW, capthick=CAPTHICK,
            label="Order excluded")
    else:
        ax_elpd.plot(g_incl["delta_elpd"], g_incl["y"], "o",
                     color=ELPD_INCL_COLOR, ms=MS, mew=MEW, label="Order included")
        ax_elpd.plot(g_excl["delta_elpd"], g_excl["y"], "s",
                     color=ELPD_EXCL_COLOR, ms=MS, mew=MEW, label="Order excluded")

    _finish_axis(ax_elpd, xlabel="ELPD difference", yticks=y_m, yticklabels=y_labels,
                 legend_loc="upper left", vline0=True)

    ax_coef.errorbar(mu_mean, y_mu, xerr=mu_xerr,
                     fmt="D", color=MU_COLOR, ecolor=MU_COLOR,
                     ms=MS, mew=MEW, capsize=CAPSIZE, elinewidth=ERR_LW, capthick=CAPTHICK,
                     label="μ")
    ax_coef.errorbar(k_mean, y_k, xerr=k_xerr,
                     fmt="^", color=K_COLOR, ecolor=K_COLOR,
                     ms=MS, mew=MEW, capsize=CAPSIZE, elinewidth=ERR_LW, capthick=CAPTHICK,
                     label="κ")
    _finish_axis(ax_coef, xlabel="Coefficient value", yticks=y_base, yticklabels=terms_order,
                 legend_loc="upper right", vline0=True)

    fig.tight_layout()
    _save(fig, fname)
    plt.show()


class _HandlerPatchOverLine(HandlerBase):
    def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height,
                       fontsize, trans):
        band, line = orig_handle
        r = Rectangle(
            (xdescent, ydescent), width, height,
            facecolor=band.get_facecolor(), edgecolor="none",
            alpha=band.get_alpha(), transform=trans)
        l = Line2D(
            [xdescent, xdescent + width],
            [ydescent + 0.5 * height, ydescent + 0.5 * height],
            color=line.get_color(), lw=line.get_linewidth(), transform=trans)
        return [r, l]


# ---------------------------------------------------------------------------
# Cell 48: two-column small multiples. Left: alignment vs Order, one row per
# density. Right: alignment vs Density, one row per order. Model prediction
# curves + HDI band (pred_i / pred_ii) overlaid with the empirical bootstrap
# means/HDIs (emp) at the observed cells. The notebook computes pred_i,
# pred_ii and emp first.
# ---------------------------------------------------------------------------
def plot_alignment_ppc_bootstrap(dens_levels, order_levels, order_grid,
                                 density_grid, pred_i, pred_ii, emp,
                                 fname="order_density_alignment_ppc_bootstrap"):
    FONT_SIZE = 7
    AX_LW     = 0.8
    TICK_LW   = 0.8
    TICK_LEN  = 2.5

    mpl.rcParams.update({
        "font.size": FONT_SIZE, "axes.labelsize": FONT_SIZE,
        "axes.titlesize": FONT_SIZE, "legend.fontsize": FONT_SIZE,
        "xtick.labelsize": FONT_SIZE, "ytick.labelsize": FONT_SIZE,
        "axes.linewidth": AX_LW,
        "xtick.major.width": TICK_LW, "ytick.major.width": TICK_LW,
        "xtick.major.size": TICK_LEN, "ytick.major.size": TICK_LEN,
        "savefig.dpi": DPI, "figure.dpi": DPI,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })

    def emp_for_density(d):
        sub = emp.loc[np.isclose(emp["density"], float(d))].sort_values("order")
        return None if sub.empty else sub

    def emp_for_order(o):
        sub = emp.loc[np.isclose(emp["order"], float(o))].sort_values("density")
        return None if sub.empty else sub

    n_rows = max(len(dens_levels), len(order_levels))
    fig_h = max(2.2, 1.25 * n_rows)
    fig = plt.figure(figsize=(7.4, fig_h), dpi=DPI)
    gs = fig.add_gridspec(n_rows, 2, wspace=0.28, hspace=0.35, width_ratios=[1.0, 1.0])

    dens_cols = [plt.cm.Reds(v)    for v in np.linspace(0.95, 0.58, len(dens_levels))]
    ord_cols  = [plt.cm.Purples(v) for v in np.linspace(0.95, 0.58, len(order_levels))]

    Y_LIM = (-0.55, 1.05)
    Y_TICKS = [-0.5, 0.0, 0.5, 1.0]
    BAND_ALPHA = 0.28
    LINE_LW = 1.8
    DATA_LW = 0.9

    for r in range(n_rows):
        ax = fig.add_subplot(gs[r, 0])
        if r >= len(dens_levels):
            ax.axis("off"); continue
        d = float(dens_levels[r]); col = dens_cols[r]

        y = np.asarray(pred_i[d]["mean"], float)
        lo = np.asarray(pred_i[d]["low"], float)
        hi = np.asarray(pred_i[d]["high"], float)

        ax.plot(order_grid, y, color=col, lw=LINE_LW, zorder=3)
        ax.fill_between(order_grid, lo, hi, color=col, alpha=BAND_ALPHA, linewidth=0, zorder=2)

        sub = emp_for_density(d)
        if sub is not None:
            x = sub["order"].to_numpy(float)
            m = sub["emp_mean"].to_numpy(float)
            el = m - sub["emp_low"].to_numpy(float)
            eh = sub["emp_high"].to_numpy(float) - m
            nA = sub["n_animals"].to_numpy(int)
            ms = 2.6 + 0.12 * np.sqrt(nA)
            ax.errorbar(x, m, yerr=np.vstack([el, eh]),
                        fmt="o", ms=0, lw=DATA_LW, color=col, zorder=4, capsize=2)
            ax.scatter(x, m, s=(ms**2), facecolors="white", edgecolors=col,
                       linewidths=DATA_LW, zorder=5)

        ax.set_title(f"Density = {d:g}", loc="center", pad=2)
        ax.set_xlim(-0.05, 1.05); ax.set_ylim(*Y_LIM); ax.set_yticks(Y_TICKS)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        ax.axhline(0, linestyle="--", color="k", lw=0.8, zorder=1)
        if r < n_rows - 1:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel("Order"); ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax.set_ylabel("Alignment" if r == 0 else "")

    for r in range(n_rows):
        ax = fig.add_subplot(gs[r, 1])
        if r >= len(order_levels):
            ax.axis("off"); continue
        o = float(order_levels[r]); col = ord_cols[r]

        y = np.asarray(pred_ii[o]["mean"], float)
        lo = np.asarray(pred_ii[o]["low"], float)
        hi = np.asarray(pred_ii[o]["high"], float)

        ax.plot(density_grid, y, color=col, lw=LINE_LW, zorder=3)
        ax.fill_between(density_grid, lo, hi, color=col, alpha=BAND_ALPHA, linewidth=0, zorder=2)

        sub = emp_for_order(o)
        if sub is not None:
            x = sub["density"].to_numpy(float)
            m = sub["emp_mean"].to_numpy(float)
            el = m - sub["emp_low"].to_numpy(float)
            eh = sub["emp_high"].to_numpy(float) - m
            nA = sub["n_animals"].to_numpy(int)
            ms = 2.6 + 0.12 * np.sqrt(nA)
            ax.errorbar(x, m, yerr=np.vstack([el, eh]),
                        fmt="o", ms=0, lw=DATA_LW, color=col, zorder=4, capsize=2)
            ax.scatter(x, m, s=(ms**2), facecolors="white", edgecolors=col,
                       linewidths=DATA_LW, zorder=5)

        ax.set_title(f"Order = {o:g}", loc="center", pad=2)
        ax.set_xscale("log", base=2); ax.set_xlim(0.9, 70); ax.set_ylim(*Y_LIM)
        ax.set_yticks(Y_TICKS)
        ax.set_xticks([1, 2, 4, 8, 16, 32, 64])
        ax.set_xticklabels([1, 2, 4, 8, 16, 32, 64])
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        ax.axhline(0, linestyle="--", color="k", lw=0.8, zorder=1)
        ax.set_ylabel("")
        if r < n_rows - 1:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel("Density")
        if r == 5:
            model_handle = (Patch(facecolor="0.70", alpha=BAND_ALPHA, edgecolor="none"),
                            Line2D([0], [0], color="k", lw=LINE_LW))
            data_handle = Line2D([0], [0], marker="o", linestyle="none",
                                 markerfacecolor="white", markeredgecolor="k",
                                 color="k", markersize=5)
            ax.legend([model_handle, data_handle], ["Model", "Data"],
                      handler_map={tuple: _HandlerPatchOverLine()},
                      frameon=True, loc="upper right", borderaxespad=0.1,
                      handlelength=1.2)

    fig.tight_layout()
    _save(fig, fname, png=False)
    plt.show()


def _rc_params_46_49():
    mpl.rcParams.update({
        "font.size": 8, "ytick.major.width": 1, "xtick.major.width": 1,
        "axes.linewidth": 1, "font.family": "Arial",
    })


# ---------------------------------------------------------------------------
# Cell 49: small multiples of dots + errorbars, model (filled) vs data (open),
# at the observed (density, order) cells. Left: vs Order per density;
# Right: vs Density per order. ``emp`` = data bootstrap, ``mod`` = model PPC.
# ---------------------------------------------------------------------------
def plot_cell_ppc_dots(dens_levels, order_levels, emp, mod,
                       fname="order_density_cell_ppc_dots"):
    _rc_params_46_49()

    n_rows = max(len(dens_levels), len(order_levels))
    fig_h = max(2.2, 1.25 * n_rows)
    fig = plt.figure(figsize=(7.4, fig_h), dpi=300)
    gs = fig.add_gridspec(n_rows, 2, wspace=0.28, hspace=0.35, width_ratios=[1.0, 1.0])

    dens_cols = [plt.cm.Reds(v) for v in np.linspace(0.90, 0.45, len(dens_levels))]
    ord_cols  = [plt.cm.Purples(v) for v in np.linspace(0.90, 0.45, len(order_levels))]

    Y_LIM = (-0.55, 1.05)
    Y_TICKS = [-0.5, 0.0, 0.5, 1.0]
    DX_ORDER = 0.02
    DENS_FAC = 1.08

    for r in range(n_rows):
        ax = fig.add_subplot(gs[r, 0])
        if r >= len(dens_levels):
            ax.axis("off"); continue
        d = float(dens_levels[r]); col = dens_cols[r]

        subd = emp.loc[np.isclose(emp["density"], d)].sort_values("order")
        subm = mod.loc[np.isclose(mod["density"], d)].sort_values("order")
        if subd.empty and subm.empty:
            ax.axis("off"); continue
        sub = pd.merge(subd, subm, on=["density", "order"], how="outer",
                       suffixes=("", "_model")).sort_values("order")

        x = sub["order"].to_numpy(float)
        x_data = x - DX_ORDER; x_model = x + DX_ORDER
        nA = sub["n_animals"].to_numpy()
        if np.all(pd.isna(nA)):
            nA = sub["n_animals_model"].to_numpy()
        nA = np.where(pd.isna(nA), 1.0, nA).astype(float)
        ms = 2.6 + 0.12 * np.sqrt(nA)

        if "emp_mean" in sub.columns:
            m = sub["emp_mean"].to_numpy(float)
            lo = sub["emp_low"].to_numpy(float); hi = sub["emp_high"].to_numpy(float)
            ok = ~np.isnan(m) & ~np.isnan(lo) & ~np.isnan(hi)
            if np.any(ok):
                ax.errorbar(x_data[ok], m[ok], yerr=np.vstack([m[ok]-lo[ok], hi[ok]-m[ok]]),
                            fmt="o", ms=0, lw=0.9, color=col, zorder=3, capsize=2)
                ax.scatter(x_data[ok], m[ok], s=(ms[ok]**2), facecolors="white",
                           edgecolors=col, linewidths=0.9, zorder=4)
        if "model_mean" in sub.columns:
            m = sub["model_mean"].to_numpy(float)
            lo = sub["model_low"].to_numpy(float); hi = sub["model_high"].to_numpy(float)
            ok = ~np.isnan(m) & ~np.isnan(lo) & ~np.isnan(hi)
            if np.any(ok):
                ax.errorbar(x_model[ok], m[ok], yerr=np.vstack([m[ok]-lo[ok], hi[ok]-m[ok]]),
                            fmt="o", ms=0, lw=1.0, color=col, zorder=5, capsize=2)
                ax.scatter(x_model[ok], m[ok], s=(ms[ok]**2)*0.9, facecolors=col,
                           edgecolors=col, linewidths=0.8, zorder=6)

        ax.set_title(f"Density = {d:g}", loc="center", pad=2)
        ax.set_xlim(-0.05, 1.05); ax.set_ylim(*Y_LIM); ax.set_yticks(Y_TICKS)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        if r < n_rows - 1:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel("Order"); ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax.set_ylabel("Alignment" if r == 0 else "")
        ax.axhline(0, linestyle="--", color="k", zorder=0)

    for r in range(n_rows):
        ax = fig.add_subplot(gs[r, 1])
        if r >= len(order_levels):
            ax.axis("off"); continue
        o = float(order_levels[r]); col = ord_cols[r]

        subd = emp.loc[np.isclose(emp["order"], o)].sort_values("density")
        subm = mod.loc[np.isclose(mod["order"], o)].sort_values("density")
        if subd.empty and subm.empty:
            ax.axis("off"); continue
        sub = pd.merge(subd, subm, on=["density", "order"], how="outer",
                       suffixes=("", "_model")).sort_values("density")

        x = sub["density"].to_numpy(float)
        x_data = x / DENS_FAC; x_model = x * DENS_FAC
        nA = sub["n_animals"].to_numpy()
        if np.all(pd.isna(nA)):
            nA = sub["n_animals_model"].to_numpy()
        nA = np.where(pd.isna(nA), 1.0, nA).astype(float)
        ms = 2.6 + 0.12 * np.sqrt(nA)

        if "emp_mean" in sub.columns:
            m = sub["emp_mean"].to_numpy(float)
            lo = sub["emp_low"].to_numpy(float); hi = sub["emp_high"].to_numpy(float)
            ok = ~np.isnan(m) & ~np.isnan(lo) & ~np.isnan(hi)
            if np.any(ok):
                ax.errorbar(x_data[ok], m[ok], yerr=np.vstack([m[ok]-lo[ok], hi[ok]-m[ok]]),
                            fmt="o", ms=0, lw=0.9, color=col, zorder=3, capsize=2)
                ax.scatter(x_data[ok], m[ok], s=(ms[ok]**2), facecolors="white",
                           edgecolors=col, linewidths=0.9, zorder=4)
        if "model_mean" in sub.columns:
            m = sub["model_mean"].to_numpy(float)
            lo = sub["model_low"].to_numpy(float); hi = sub["model_high"].to_numpy(float)
            ok = ~np.isnan(m) & ~np.isnan(lo) & ~np.isnan(hi)
            if np.any(ok):
                ax.errorbar(x_model[ok], m[ok], yerr=np.vstack([m[ok]-lo[ok], hi[ok]-m[ok]]),
                            fmt="o", ms=0, lw=1.0, color=col, zorder=5, capsize=2)
                ax.scatter(x_model[ok], m[ok], s=(ms[ok]**2)*0.9, facecolors=col,
                           edgecolors=col, linewidths=0.8, zorder=6)

        ax.set_title(f"Order = {o:g}", loc="center", pad=2)
        if r == 5:
            data_proxy = Line2D([0], [0], marker="o", linestyle="none",
                                markerfacecolor="white", markeredgecolor="black",
                                markeredgewidth=1.0, color="black", markersize=5)
            model_proxy = Line2D([0], [0], marker="o", linestyle="none",
                                 markerfacecolor="black", markeredgecolor="black",
                                 color="black", markersize=5)
            ax.legend([model_proxy, data_proxy], ["Model", "Data"],
                      loc="upper right", frameon=True, borderaxespad=0.1,
                      handlelength=0.8, fontsize=7)
        ax.set_xscale("log", base=2); ax.set_xlim(0.8, 75); ax.set_ylim(*Y_LIM)
        ax.set_yticks(Y_TICKS)
        ax.set_xticks([1, 2, 4, 8, 16, 32, 64])
        ax.set_xticklabels([1, 2, 4, 8, 16, 32, 64])
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        if r < n_rows - 1:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel("Density")
        ax.set_ylabel("")
        ax.axhline(0, linestyle="--", color="k", zorder=0)

    fig.tight_layout()
    _save(fig, fname)
    plt.show()


# ---------------------------------------------------------------------------
# Cell 50: grid of subplots by (density, order); x = individual,
# y = alignment. Data = open black dots + bootstrap CI; Model = filled red
# dots + PPC CI. The notebook computes ``ppc_df`` first.
# ---------------------------------------------------------------------------
def plot_individual_cell_ppc(ppc_df,
                             fname="order_density_individual_cell_ppc"):
    _rc_params_46_49()

    dens_levels = np.sort(ppc_df["density"].unique())
    order_levels = np.sort(ppc_df["order"].unique())
    nrows = len(dens_levels); ncols = len(order_levels)

    fig_w = max(6.8, 2.2 * ncols)
    fig_h = max(3.2, 1.8 * nrows)
    fig, axs = plt.subplots(nrows, ncols, figsize=(fig_w, fig_h), dpi=300, sharey=True)
    axs = np.atleast_2d(axs).reshape(nrows, ncols)

    Y_LIM = (-1.05, 1.05)
    Y_TICKS = [-1.0, -0.5, 0.0, 0.5, 1.0]

    for ri, d in enumerate(dens_levels):
        for ci, o in enumerate(order_levels):
            ax = axs[ri, ci]
            sub = ppc_df.loc[(np.isclose(ppc_df["density"], d)) &
                             (np.isclose(ppc_df["order"], o))].copy()
            if sub.empty:
                ax.axis("off"); continue
            sub = sub.sort_values("animal_id").reset_index(drop=True)
            x = np.arange(len(sub), dtype=float)
            x_data = x - 0.14; x_model = x + 0.14

            ax.vlines(x_data, sub["data_low"], sub["data_high"], colors="black", lw=0.9, zorder=2)
            ax.scatter(x_data, sub["data_mean"], s=18, facecolors="white", edgecolors="black",
                       linewidths=0.9, zorder=3)
            ax.vlines(x_model, sub["model_low"], sub["model_high"], colors="tab:red", lw=1.0, zorder=1)
            ax.scatter(x_model, sub["model_mean"], s=16, facecolors="tab:red", edgecolors="tab:red",
                       linewidths=0.8, zorder=2)

            ax.set_ylim(*Y_LIM); ax.set_yticks(Y_TICKS)
            ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
            ax.set_xticks(x)
            ax.set_xticklabels(sub["animal_id"].tolist(), rotation=90, fontsize=6)
            if ri == 0:
                ax.set_title(f"Order={o:g}", pad=2)
            if ci == 0:
                ax.set_ylabel(f"Density={d:g}\nAlignment")
            else:
                ax.set_ylabel("")
            ax.axhline(0, linestyle="--", color="k", zorder=0)

    data_proxy = Line2D([0], [0], marker="o", linestyle="none",
                        markerfacecolor="white", markeredgecolor="black",
                        markeredgewidth=1.0, color="black", markersize=5)
    model_proxy = Line2D([0], [0], marker="o", linestyle="none",
                         markerfacecolor="tab:red", markeredgecolor="tab:red",
                         color="tab:red", markersize=5)
    axs[-1, -1].legend([model_proxy, data_proxy], ["Model", "Data"],
                       loc="lower right", frameon=True, borderaxespad=0.2,
                       handletextpad=0.4)

    fig.tight_layout(rect=[0, 0, 0.97, 1])
    _save(fig, fname)
    plt.show()
