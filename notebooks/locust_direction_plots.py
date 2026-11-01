"""Figure helpers for the locust directional mixture-model notebook.

The notebook keeps the statistical work (loading the data, fitting the
hierarchical bimodal and unimodal von Mises mixture models); the pointwise-LOO
aggregation, the population density estimates (posterior predictive density,
bootstrap histogram, von Mises KDE), the per-individual PPC arrays and the
hypothesis tests live in ``locust_posthoc``. Each ``plot_*`` function takes
the fit result and the modelling table, calls into ``locust_posthoc`` for
the quantities it draws, and is a single call from the notebook.

Style and saving follow the shared house convention from ``unwrap.plots``:
figures render at display resolution (100 dpi) in the notebook cell and are
saved as 300 dpi vector files under ``<repo>/figures/locust/``. Every public
figure function accepts a ``save`` flag (on by default) and an explicit
repo-relative ``save_path`` (the notebook always passes it); the written
files are gitignored.

Palette: Bimodal / "Vec. encoding" is orange, Unimodal / "Optic flow" is
blue, the posterior predictive model curves are red (colorblind-safe,
matching the paper figures).
"""
from pathlib import Path

import matplotlib.lines as mlines
import matplotlib.legend_handler as ml_handler
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
from scipy.special import i0, iv
from scipy.stats import vonmises

from unwrap import plots as _plots

import locust_posthoc as posthoc

#: Publication figures for this study live in
#: <repo>/figures/locust/, anchored on this file's location so the output
#: does not depend on the kernel cwd (the notebook may run from notebooks/
#: or from the repository root).
FIG_DIR = Path(__file__).resolve().parent.parent / "figures" / "locust"

FONT_SZ = 8
AX_LW = 1.0
TICK_W = 1.0
LINE_W = 1.0

BIMODAL = "tab:orange"   # Vec. encoding
UNIMODAL = "tab:blue"    # Optic flow
MODEL_RED = "tab:red"    # posterior predictive model curves

_EXTRA_RC = {
    "font.size": FONT_SZ,
    "axes.linewidth": AX_LW,
    "xtick.major.width": TICK_W,
    "ytick.major.width": TICK_W,
    "lines.linewidth": LINE_W,
}

# Per-individual two-panel figures: ROW_GAP / Y_OFF keep the rows spaced and
# the two markers offset per row.
ROW_GAP = 1.2
Y_OFF = 0.23


def set_style():
    """One figure style for every figure in this notebook (the shared house
    convention from ``unwrap.plots`` plus the study's type size; figures are
    saved 300 dpi under ``FIG_DIR`` while rendering 100 dpi inline)."""
    _plots.set_style(FIG_DIR, _EXTRA_RC)


def _save(fig, save, save_path):
    """Write ``fig`` to ``save_path`` (a repo-relative path such as
    ``figures/locust/x.pdf``) through the house two-res saver. The file
    format follows the path's suffix."""
    if not save:
        return
    p = Path(save_path)
    if p.suffix.lower() == ".svg":
        _plots.save_figure(fig, FIG_DIR, p.stem, pdf=False, png=False, svg=True)
    else:
        _plots.save_figure(fig, FIG_DIR, p.stem, pdf=True, png=False, svg=False)


def _despine(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for sp in ax.spines.values():
        sp.set_linewidth(AX_LW)
    ax.tick_params(axis="both", which="both",
                   bottom=True, top=False, left=True, right=False,
                   width=TICK_W)


def _asym_xerr(center, hdi):
    return np.vstack([center - hdi[:, 0], hdi[:, 1] - center])


def _indiv_figure(n_rows, row_gap=ROW_GAP):
    fig_h = max(2.8, 0.22 * n_rows)
    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=(8.2, fig_h), dpi=_plots.DEFAULT_DISPLAY_DPI,
        sharey=True, gridspec_kw={"width_ratios": [1.35, 1.0]},
    )
    y_base = np.arange(n_rows) * row_gap
    y_sep = y_base[:-1] + row_gap / 2.0
    for ax in (axA, axB):
        for ys in y_sep:
            ax.axhline(ys, color="0.90", lw=0.8, zorder=0)
    return fig, axA, axB, y_base


def _panel_a_frame(axA, y_base, ytick_labels, xlabel, y_of, y_ve):
    axA.set_yticks(y_base)
    axA.set_yticklabels(ytick_labels)
    axA.invert_yaxis()
    axA.set_ylabel("Individual")
    axA.set_xlabel(xlabel)
    axA.spines["top"].set_visible(False)
    axA.spines["right"].set_visible(False)
    ymin = y_of.min() - 0.35
    ymax = y_ve.max() + 0.35
    axA.set_ylim(ymax, ymin)
    axA.text(-0.12, 1.02, "a.", transform=axA.transAxes, fontweight="bold",
             fontfamily="Arial", fontsize=16, va="bottom", ha="left")


def _panel_b_p_gt(axB, p_gt_all, y_base, xlabel):
    axB.plot(p_gt_all, y_base, "o", color="k", markersize=3.5, zorder=3)
    axB.axvline(0.5, ls="--", c="k", lw=1, zorder=1)
    axB.axvline(1.0, ls="-", c="k", lw=1, zorder=1)
    axB.set_xlim(0, 1.02)
    axB.set_xlabel(xlabel)
    axB.tick_params(axis="y", left=False, labelleft=False)
    axB.spines["top"].set_visible(False)
    axB.spines["right"].set_visible(False)
    axB.text(-0.12, 1.02, "b.", transform=axB.transAxes, fontweight="bold",
             fontfamily="Arial", fontsize=16, va="bottom", ha="left")


def _horizontal_dots(vals, xerr, xlabel, uniform=None, invert=False,
                     xlim=None, uniform_label="Uniform", uniform_text_off=0.01):
    """Small horizontal-dot figure (3.0 x 1.25 in) shared by the population
    hypothesis / alignment figures:
      row 0 = "Optic flow" (forward, blue, circle)
      row 1 = "Vec. encoding" (perpendicular, orange, square)
    ``uniform`` is an optional vertical reference line."""
    labels = ["Optic flow", "Vec. encoding"]
    colors = [UNIMODAL, BIMODAL]
    markers = ["o", "s"]
    y = np.arange(2)
    vals = np.asarray(vals, dtype=float)
    xerr = np.asarray(xerr, dtype=float)
    symmetric = xerr.ndim == 1  # (2,) mean+-, or (2,2) [left, right]

    fig, ax = plt.subplots(figsize=(3.0, 1.25), dpi=_plots.DEFAULT_DISPLAY_DPI)
    for yi, m, xe, col, mk in zip(y, vals, xerr, colors, markers):
        xerr_arg = xe if symmetric else np.array([[xe[0]], [xe[1]]])
        ax.errorbar(
            m, yi, xerr=xerr_arg, fmt=mk,
            color=col, ecolor=col,
            elinewidth=1, capsize=3, capthick=1,
            markersize=6, markeredgecolor=col, markeredgewidth=1, zorder=3,
        )

    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    for lab, col in zip(ax.get_yticklabels(), colors):
        lab.set_color(col)

    if invert:
        ax.invert_yaxis()
    ax.margins(y=0.35)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("")
    ax.set_title("")
    _despine(ax)

    if uniform is not None:
        ax.axvline(uniform, ls="--", c="k", lw=1, zorder=1)
        if uniform_label:
            ax.text(
                uniform + uniform_text_off, 0.8, uniform_label,
                rotation=90, va="top", ha="center", fontsize=FONT_SZ, color="k",
                transform=ax.get_xaxis_transform(), clip_on=False,
            )
    if xlim is not None:
        ax.set_xlim(*xlim)

    plt.tight_layout()
    return fig, ax


def _density_legend_handler(color, linestyle=None, alpha=0.25, linewidth=1.5):
    """Legend handler: filled patch (alpha) + centre line, for density bands."""
    class HandlerLineInPatch(ml_handler.HandlerBase):
        def create_artists(self, legend, orig_handle, xdescent, ydescent,
                           width, height, fontsize, trans):
            patch = mpatches.Rectangle(
                xy=(-xdescent, -ydescent), width=width, height=height,
                facecolor=color, edgecolor="none", alpha=alpha, transform=trans,
            )
            line = mlines.Line2D(
                [-0.05 * width, 1.05 * width], [0.5 * height, 0.5 * height],
                color=color, linestyle=linestyle if linestyle else "-",
                linewidth=linewidth, transform=trans,
            )
            return [patch, line]
    return HandlerLineInPatch()


# ---------------------------------------------------------------------------
# Trajectories: two XY scatter/line panels of every locust.
# ---------------------------------------------------------------------------
def plot_trajectories(df, *, save=True,
                      save_path="figures/locust/locust_directional_trajectories.pdf"):
    """Two side-by-side XY panels of every locust's path — panel 0 zoomed
    out (-30..30 m, `','` markers), panel 1 zoomed in (-5..5 m, `'o'`
    markers); one line per animal, colours cycled through the tab20
    colormap."""
    animal_ids = list(df["animal_id"].unique())
    colors = plt.cm.tab20.colors

    fig, axs = plt.subplots(1, 2, figsize=(12, 10), constrained_layout=True)

    for idx, animal_id in enumerate(animal_ids):
        individual_df = df[df["animal_id"] == animal_id]

        axs[0].plot(individual_df["X"], individual_df["Y"], marker=",",
                    linestyle="-", color=colors[idx % 20])
        axs[0].set_xlim(-30, 30)
        axs[0].set_ylim(-30, 30)
        axs[0].set_xlabel("X Coordinate (meters)")
        axs[0].set_ylabel("Y Coordinate (meters)")
        axs[0].set_aspect(1)

        axs[1].plot(individual_df["X"], individual_df["Y"], marker="o",
                    linestyle="-", color=colors[idx % 20])
        axs[1].set_xlim(-5, 5)
        axs[1].set_ylim(-5, 5)
        axs[1].set_xlabel("X Coordinate (meters)")
        axs[1].set_ylabel("Y Coordinate (meters)")
        axs[1].set_aspect(1)

    _save(fig, save, save_path)


# ---------------------------------------------------------------------------
# LOO ELPD for bimodal vs unimodal (dot + errorbar).
# ---------------------------------------------------------------------------
def plot_loo_elpd(models, *, save=True,
                  save_path="figures/locust/locust_bimodal_unimodal_elpd.pdf"):
    """LOO expected log predictive density for the two models. ``models`` is
    the name->fit dictionary in [Bimodal, Unimodal] order (the ELPD + SE are
    computed in ``locust_posthoc``); the labels, colours (orange = Bimodal,
    blue = Unimodal) and markers (square, circle) are fixed. Horizontal
    dot+errorbar, 3.0 x 1.25 in."""
    rows = posthoc.loo_elpd_rows(models)

    # Fixed order: Bimodal (orange, square) then Unimodal (blue, circle).
    labels = ["Bimodal\n(Vec. encoding)", "Unimodal\n(Optic flow)"]
    colors = [BIMODAL, UNIMODAL]
    markers = ["s", "o"]

    x = np.array([r[1] for r in rows], dtype=float)
    xerr = np.array([r[2] for r in rows], dtype=float)
    y = np.arange(len(rows))

    fig, ax = plt.subplots(figsize=(3.0, 1.25), dpi=_plots.DEFAULT_DISPLAY_DPI)
    for yi, xi, sei, col, mk in zip(y, x, xerr, colors, markers):
        ax.errorbar(
            xi, yi, xerr=sei,
            fmt=mk,
            color=col, ecolor=col,
            elinewidth=1, capsize=3, capthick=1,
            markersize=6,
            markeredgecolor=col,
            markeredgewidth=1,
        )

    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    for lab, col in zip(ax.get_yticklabels(), colors):
        lab.set_color(col)

    ax.margins(y=0.35)
    ax.set_xlabel("ELPD")
    ax.set_ylabel("")
    ax.set_title("")
    _despine(ax)

    plt.tight_layout()
    _save(fig, save, save_path)


# ---------------------------------------------------------------------------
# Per-individual stratified ELPD + dELPD.
# ---------------------------------------------------------------------------
def plot_elpd_by_individual(models, df, *, save=True,
                            save_path="figures/locust/locust_elpd_by_individual_ab.pdf"):
    """Panel a: per-individual stratified ELPD (orange square = Bimodal,
    blue circle = Unimodal, offset per individual). Panel b: dELPD
    (Bimodal - Unimodal) per individual, black dots, zero dashed
    reference. The pointwise LOO arithmetic lives in
    ``locust_posthoc.elpd_by_individual``."""
    stats = posthoc.elpd_by_individual(models, df)
    elpd_bi, elpd_uni = stats["elpd_bi"], stats["elpd_uni"]
    se_bi, se_uni = stats["se_bi"], stats["se_uni"]
    delta, se_delta = stats["delta"], stats["se_delta"]

    n_ind = len(elpd_bi)
    y = np.arange(n_ind)
    fig_h = max(2.8, 0.18 * n_ind)
    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=(8.2, fig_h), dpi=_plots.DEFAULT_DISPLAY_DPI,
        sharey=True, gridspec_kw={"width_ratios": [1.35, 1.0]},
    )
    for ax in (axA, axB):
        for ys in y[:-1] + 0.5:
            ax.axhline(ys, color="0.90", lw=0.8, zorder=0)

    axA.errorbar(elpd_bi, y + 0.15, xerr=se_bi, fmt="s", color=BIMODAL,
                 ecolor=BIMODAL, elinewidth=1, capsize=2, capthick=1,
                 markersize=4.5, markeredgewidth=1, zorder=3,
                 label="Bimodal (Vec. encoding)")
    axA.errorbar(elpd_uni, y - 0.15, xerr=se_uni, fmt="o", color=UNIMODAL,
                 ecolor=UNIMODAL, elinewidth=1, capsize=2, capthick=1,
                 markersize=4.5, markeredgewidth=1, zorder=3,
                 label="Unimodal (Optic flow)")
    axA.set_yticks(y)
    axA.set_yticklabels([str(i + 1) for i in range(n_ind)])
    axA.invert_yaxis()
    axA.set_xlabel("ELPD contribution (summed within individual)")
    axA.set_ylabel("Individual")
    axA.spines["top"].set_visible(False)
    axA.spines["right"].set_visible(False)
    leg = axA.legend(frameon=True, loc="lower right")
    leg.get_frame().set_facecolor("white")
    leg.get_frame().set_alpha(1.0)
    axA.text(-0.12, 1.02, "a.", transform=axA.transAxes, fontweight="bold",
             fontfamily="Arial", fontsize=16, va="bottom", ha="left")

    axB.errorbar(delta, y, xerr=se_delta, fmt="o", color="k", ecolor="k",
                 elinewidth=1, capsize=2, capthick=1, markersize=4.5, zorder=3)
    axB.axvline(0.0, ls="--", c="k", lw=1, zorder=1)
    axB.set_xlabel(r"ELPD difference (Bimodal $-$ Unimodal)")
    axB.tick_params(axis="y", left=False, labelleft=False)
    axB.spines["top"].set_visible(False)
    axB.spines["right"].set_visible(False)
    axB.text(-0.12, 1.02, "b.", transform=axA.transAxes, fontweight="bold",
             fontfamily="Arial", fontsize=16, va="bottom", ha="left")

    plt.tight_layout()
    _save(fig, save, save_path)


# ---------------------------------------------------------------------------
# Population-level directional validation figures (cartesian + polar).
#
# Both plot the same three estimates of the *population* directional density
# (posterior predictive model, histogram with bootstrap CI, Von-Mises KDE)
# computed in ``locust_posthoc.population_validation_stats``: the cartesian
# two-panel figure and the polar figure differ only in projection.
# ---------------------------------------------------------------------------
def plot_directional_validation_cartesian(trace_bimodal, df, *, save=True,
                                          save_path="figures/locust/population_directional_validation_cartesian.pdf"):
    """Two stacked panels: (top) model density band + histogram bootstrap
    error bars; (bottom) adds the Von-Mises KDE. X-axis in radians with
    pi-fraction ticks."""
    stats = posthoc.population_validation_stats(trace_bimodal, df)
    angles = stats["angles"]
    density_median, density_lower, density_upper = (
        stats["density_median"], stats["density_lower"], stats["density_upper"])
    kde_angles, kde_median = stats["kde_angles"], stats["kde_median"]
    kde_lower_bounds, kde_upper_bounds = (
        stats["kde_lower_bounds"], stats["kde_upper_bounds"])
    hist_bin_centers, hist_median = stats["hist_bin_centers"], stats["hist_median"]
    hist_errors = stats["hist_errors"]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(1.1 * 2, 1.1 * 2),
                                   sharex=True,
                                   dpi=_plots.DEFAULT_DISPLAY_DPI)

    for ax in (ax1, ax2):
        line_plot, = ax.plot(angles, density_median, MODEL_RED,
                             label="Posterior Density")
        ax.fill_between(angles, density_lower, density_upper, color=MODEL_RED,
                        alpha=0.25)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    line_plot_kde, = ax2.plot(kde_angles, kde_median, "black", linestyle="--",
                              label="Median KDE")
    ax2.fill_between(kde_angles, kde_lower_bounds, kde_upper_bounds,
                     color="black", alpha=0.25)

    errorbar_plot = ax1.errorbar(
        hist_bin_centers, hist_median, yerr=hist_errors, fmt="o", ecolor="black",
        markeredgecolor="black", markerfacecolor="black", linewidth=1, capsize=2,
        markersize=2, label="Bootstrap Density\n(Median ± 95% C.I.)")

    ax2.set_xlabel("Direction (radians)")
    ax2.set_xticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi],
                   [r"$-\pi$", r"$-\frac{\pi}{2}$", "0", r"$\frac{\pi}{2}$", r"$\pi$"])
    plt.xlim(-np.pi, np.pi)
    fig.text(-0.1, 0.5, "Density", va="center", rotation="vertical")

    plt.legend(
        [line_plot, errorbar_plot, line_plot_kde],
        ["Model", "Histogram", "Kernel"],
        handler_map={line_plot_kde: _density_legend_handler("black", "--"),
                     line_plot: _density_legend_handler(MODEL_RED)},
        bbox_to_anchor=(1, 1), frameon=False, loc="center left",
    )

    _save(fig, save, save_path)


def plot_directional_validation_polar(trace_bimodal, df, *, save=True,
                                      save_path="figures/locust/population_directional_validation.pdf"):
    """Polar view of the population directional density: model band, KDE
    band and histogram bootstrap dots stacked on a radius offset, with an
    outer reference circle and a custom patch+line legend for the two
    density bands."""
    stats = posthoc.population_validation_stats(trace_bimodal, df)
    angles = stats["angles"]
    density_median, density_lower, density_upper = (
        stats["density_median"], stats["density_lower"], stats["density_upper"])
    kde_angles, kde_median = stats["kde_angles"], stats["kde_median"]
    kde_lower_bounds, kde_upper_bounds = (
        stats["kde_lower_bounds"], stats["kde_upper_bounds"])
    hist_bin_centers, hist_median = (stats["hist_bin_centers"], stats["hist_median"])
    hist_lower_bounds, hist_upper_bounds = (
        stats["hist_lower_bounds"], stats["hist_upper_bounds"])

    edge_radius = 0.5

    fig, ax = plt.subplots(subplot_kw={"projection": "polar"},
                           figsize=(3, 3), dpi=_plots.DEFAULT_DISPLAY_DPI)

    line_plot, = ax.plot(angles, edge_radius + density_median, color=MODEL_RED,
                         linewidth=1.5, label="Model")
    ax.fill_between(angles, edge_radius + density_lower,
                    edge_radius + density_upper,
                    color=MODEL_RED, alpha=0.25, edgecolor="none")

    line_plot_kde, = ax.plot(kde_angles, edge_radius + kde_median, color="black",
                             linestyle="--", linewidth=1.5, label="Kernel")
    ax.fill_between(kde_angles, edge_radius + kde_lower_bounds,
                    edge_radius + kde_upper_bounds, color="black", alpha=0.25,
                    edgecolor="none")

    errorbar_plot = ax.errorbar(
        hist_bin_centers, edge_radius + hist_median,
        yerr=[(hist_median + edge_radius) - (hist_lower_bounds + edge_radius),
              (hist_upper_bounds + edge_radius) - (hist_median + edge_radius)],
        fmt="o", ecolor="black", markeredgecolor="black", markerfacecolor="black",
        linewidth=1, capsize=2, markersize=2, label="Histogram")

    outer_circle = Circle((0, 0), edge_radius, transform=ax.transData._b,
                          color="black", fill=False, linewidth=0.5)
    ax.add_artist(outer_circle)

    ax.spines["polar"].set_visible(False)
    ax.set_ylim(0, edge_radius + np.max(density_upper) * 1.2)
    ax.set_xlim(-np.pi, np.pi)
    ax.set_yticks([edge_radius, edge_radius + 0.5])
    ax.set_yticklabels(["", ""])
    ax.set_xticks([-np.pi / 2, 0, np.pi / 2, np.pi])
    ax.set_xticklabels(["−90°", "0°", " 90°", "±180°     "])
    ax.tick_params(axis="x", pad=-4)

    ax.legend(
        [line_plot, errorbar_plot, line_plot_kde],
        ["Model", "Histogram", "Kernel"],
        handler_map={
            line_plot: _density_legend_handler(MODEL_RED, "-", 0.25, 1.5),
            line_plot_kde: _density_legend_handler("black", "--", 0.25, 1.5),
        },
        bbox_to_anchor=(1, 1), loc="center left", frameon=True,
    )

    plt.tight_layout()
    _save(fig, save, save_path)


# ---------------------------------------------------------------------------
# 52-panel individual-level PPC polar grid.
#
# Draws, per individual, the model posterior 95% density band + mean (red),
# the bootstrap Von-Mises KDE band + median (gray), observed directions as
# outward ticks, "residual" inward ticks (inverse mean predictive density)
# and the two mixture-component mean vectors as red arrows, on a polar
# projection with density=0 / density=1 reference circles. The array
# computations stay in locust_posthoc.individual_ppc_stats.
# ---------------------------------------------------------------------------
def plot_individual_ppc_grid(trace_bimodal, df, *, save=True,
                             save_path="figures/locust/locust_directional_individual_ppc.pdf"):
    stats = posthoc.individual_ppc_stats(trace_bimodal, df)

    num_individuals = stats["n_individuals"]
    per_ind_data = stats["per_ind_data"]
    angles = stats["angles"]
    pdf_lower, pdf_upper, pdf_mean_vals = (
        stats["pdf_lower"], stats["pdf_upper"], stats["pdf_mean_vals"])
    kde_med, kde_lo, kde_hi = stats["kde_med"], stats["kde_lo"], stats["kde_hi"]
    mu_mean, kappa_mean, weights_mean = (
        stats["mu_mean"], stats["kappa_mean"], stats["weights_mean"])
    max_pdf_global = stats["max_pdf_global"]

    ncols = 7
    nrows = int(np.ceil(num_individuals / ncols))
    edge_radius = 1.0

    ref_line_color = "gray"
    kde_band_color = "0.5"
    kde_line_color = "0.1"
    model_band_color = "tab:red"
    model_line_color = "tab:red"
    data_tick_color = "black"
    residual_color = "tab:cyan"

    EPS_local = 1e-12

    fig, axs = plt.subplots(
        nrows=nrows, ncols=ncols,
        subplot_kw={"projection": "polar"},
        figsize=(2.0 * ncols, 2.0 * nrows), dpi=_plots.DEFAULT_DISPLAY_DPI,
    )
    axs = np.array(axs).reshape(-1)

    rmax = edge_radius + max_pdf_global * 1.15

    for idx in range(num_individuals):
        ax = axs[idx]

        for a in (0, np.pi, np.pi / 2, -np.pi / 2):
            ax.plot([a, a], [0, rmax], color=ref_line_color,
                    linestyle="--", linewidth=0.8, alpha=0.8)

        individual_data = per_ind_data[idx]

        if np.isfinite(kde_med[idx]).any():
            ax.fill_between(angles, edge_radius + kde_lo[idx],
                            edge_radius + kde_hi[idx], color=kde_band_color,
                            alpha=0.25, linewidth=0, zorder=1)
            ax.plot(angles, edge_radius + kde_med[idx], color=kde_line_color,
                    linestyle="--", linewidth=1.2,
                    label="Kernel (bootstrap 95%)" if idx == 0 else None, zorder=2)

        bottom = pdf_lower[idx]
        top = pdf_upper[idx]
        mean_pdf = pdf_mean_vals[idx]
        ax.fill_between(angles, edge_radius + bottom, edge_radius + top,
                        color=model_band_color, alpha=0.35, linewidth=0, zorder=3)
        ax.plot(angles, edge_radius + mean_pdf, color=model_line_color,
                linewidth=1.6, label="Model (posterior 95%)" if idx == 0 else None,
                zorder=4)

        line_length = 0.10
        for ang in individual_data:
            ax.plot([ang, ang], [edge_radius, edge_radius + line_length],
                    color=data_tick_color, alpha=0.9, linewidth=1.0, zorder=5)

        if individual_data.size > 0:
            mu_ind = mu_mean[idx]
            kappa_ind = kappa_mean[idx]
            weights_ind = np.clip(weights_mean[idx], EPS_local, 1.0)
            weights_ind = weights_ind / weights_ind.sum()

            log_c1_mean = np.log(weights_ind[0]) + vonmises.logpdf(angles, kappa_ind[0], loc=mu_ind[0])
            log_c2_mean = np.log(weights_ind[1]) + vonmises.logpdf(angles, kappa_ind[1], loc=mu_ind[1])
            log_pdf_mean_grid = np.logaddexp(log_c1_mean, log_c2_mean)
            inv_pdf_grid = np.exp(-log_pdf_mean_grid)
            norm = (np.trapezoid if hasattr(np, "trapezoid") else np.trapz)(inv_pdf_grid, angles)

            for ang in individual_data:
                log_pdf_c1 = np.log(weights_ind[0]) + vonmises.logpdf(ang, kappa_ind[0], loc=mu_ind[0])
                log_pdf_c2 = np.log(weights_ind[1]) + vonmises.logpdf(ang, kappa_ind[1], loc=mu_ind[1])
                log_pdf_val = np.logaddexp(log_pdf_c1, log_pdf_c2)
                inv_val = np.exp(-log_pdf_val)
                residual = -(inv_val / (norm + EPS_local))
                ax.plot([ang, ang], [edge_radius, edge_radius + residual],
                        color=residual_color, alpha=0.9, linewidth=1.0, zorder=6)

        mu_ind = mu_mean[idx]
        kappa_ind = kappa_mean[idx]
        weights_ind = np.clip(weights_mean[idx], EPS_local, 1.0)
        weights_ind = weights_ind / weights_ind.sum()
        R1 = iv(1, kappa_ind[0]) / i0(kappa_ind[0])
        R2 = iv(1, kappa_ind[1]) / i0(kappa_ind[1])
        arrow_length_1 = edge_radius * float(R1) * float(weights_ind[0])
        arrow_length_2 = edge_radius * float(R2) * float(weights_ind[1])

        ax.annotate("", xy=(mu_ind[0], arrow_length_1), xytext=(0, 0),
                    arrowprops=dict(facecolor=model_line_color, edgecolor=model_line_color,
                                    arrowstyle="->", lw=1.2, shrinkA=0, shrinkB=0),
                    zorder=7)
        ax.annotate("", xy=(mu_ind[1], arrow_length_2), xytext=(0, 0),
                    arrowprops=dict(facecolor=model_line_color, edgecolor=model_line_color,
                                    arrowstyle="->", lw=1.2, shrinkA=0, shrinkB=0),
                    zorder=7)

        ax.add_artist(Circle((0, 0), edge_radius, transform=ax.transData._b,
                             color="black", fill=False, linewidth=0.5, zorder=10))
        ax.add_artist(Circle((0, 0), edge_radius + 1.0, transform=ax.transData._b,
                             color="black", fill=False, linewidth=0.5, zorder=10))

        ax.set_ylim(0, rmax)
        ax.set_title(f"Individual {idx + 1}", fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines["polar"].set_visible(False)

    for ax in axs[num_individuals:]:
        ax.set_visible(False)

    handles, labels = axs[0].get_legend_handles_labels()
    if len(handles) > 0:
        fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(1.02, 1.02), frameon=True)

    plt.tight_layout()
    _save(fig, save, save_path)


# ---------------------------------------------------------------------------
# Population hypothesis / alignment dot figures.
# ---------------------------------------------------------------------------
def plot_hypothesis_density(trace_bimodal, *, save=True,
                            save_path="figures/locust/locust_directional_hypothesis_horizontal.pdf"):
    """Posterior predictive density at the hypothesized directions (forward
    vs perpendicular, computed in
    ``locust_posthoc.hypothesis_density_stats``): densities with asymmetric
    HDI error bars and the uniform 1/(2*pi) reference (dashed)."""
    stats = posthoc.hypothesis_density_stats(trace_bimodal)
    fig, _ = _horizontal_dots(stats["centers"], stats["xerr"],
                              "Predictive density",
                              uniform=stats["p_uniform"], invert=True,
                              uniform_text_off=0.01)
    _save(fig, save, save_path)


def plot_hypothesis_elpd(trace_bimodal, *, save=True,
                         save_path="figures/locust/locust_directional_hypothesis_horizontal_elpd.pdf"):
    """ELPD at the hypothesized directions (mean +/- 1.96*SE, computed in
    ``locust_posthoc.hypothesis_elpd_stats``) with the log-uniform
    reference."""
    stats = posthoc.hypothesis_elpd_stats(trace_bimodal)
    fig, _ = _horizontal_dots(stats["elpds"], stats["xerr"], "ELPD",
                              uniform=stats["logp_uniform"], invert=True,
                              uniform_text_off=0.05)
    _save(fig, save, save_path)


def plot_alignment(trace_bimodal, *, save=True,
                   save_path="figures/locust/locust_directional_alignment_horizontal.pdf"):
    """Posterior mean alignment (computed in
    ``locust_posthoc.population_alignment_stats``): forward cosine
    similarity (OF) vs axial alignment (VE), asymmetric HDI error bars,
    zero baseline; scale kept at -0.05..1.0."""
    stats = posthoc.population_alignment_stats(trace_bimodal)
    fig, _ = _horizontal_dots(stats["means"], stats["xerr"],
                              "Posterior mean alignment",
                              uniform=0.0, uniform_label="",
                              xlim=(-0.05, 1.0))
    _save(fig, save, save_path)


# ---------------------------------------------------------------------------
# Per-individual hypothesis / alignment two-panel (a/b) figures.
# ---------------------------------------------------------------------------
def plot_hypothesis_individual(trace_bimodal, *, save=True,
                               save_path="figures/locust/locust_directional_hypothesis_individual_ab.pdf"):
    """Panel a: per-individual HDI of posterior predictive density (blue
    circle = Optic flow, orange square = Vec. encoding), uniform reference,
    plus a Pooled row. Panel b: P(VE > OF) per individual + pooled. Arrays
    computed in ``locust_posthoc.hypothesis_individual_stats``."""
    stats = posthoc.hypothesis_individual_stats(trace_bimodal)
    n_ind = stats["n_ind"]
    center_f_all, center_p_all = stats["center_f_all"], stats["center_p_all"]
    hdi_f_all, hdi_p_all = stats["hdi_f_all"], stats["hdi_p_all"]
    p_gt_all = stats["p_VE_gt_OF_all"]

    n_rows = n_ind + 1
    fig, axA, axB, y_base = _indiv_figure(n_rows)
    y_of = y_base - Y_OFF
    y_ve = y_base + Y_OFF

    axA.errorbar(center_f_all, y_of, xerr=_asym_xerr(center_f_all, hdi_f_all),
                 fmt="o", color=UNIMODAL, ecolor=UNIMODAL, elinewidth=1,
                 capsize=2, markersize=3.5, alpha=0.95, label="Optic flow", zorder=3)
    axA.errorbar(center_p_all, y_ve, xerr=_asym_xerr(center_p_all, hdi_p_all),
                 fmt="s", color=BIMODAL, ecolor=BIMODAL, elinewidth=1,
                 capsize=2, markersize=3.5, alpha=0.95, label="Vec. encoding", zorder=3)
    axA.axvline(stats["p_uniform"], ls="--", c="k", lw=1, zorder=1)
    _panel_a_frame(axA, y_base, [str(i + 1) for i in range(n_ind)] + ["Pooled"],
                   "Predictive density", y_of, y_ve)
    leg = axA.legend(frameon=True, loc="upper right")
    leg.get_frame().set_facecolor("white")
    leg.get_frame().set_alpha(1.0)

    _panel_b_p_gt(axB, p_gt_all, y_base, "P(Vec. encoding > Optic flow)")

    fig.tight_layout()
    _save(fig, save, save_path)


def plot_hypothesis_elpd_individual(trace_bimodal, *, save=True,
                                    save_path="figures/locust/locust_directional_hypothesis_individual_elpd_ab.pdf"):
    """Panel a: per-individual ELPD at the two hypotheses, mean +/- 1.96*SE
    (blue circle = OF, orange square = VE), log-uniform reference, plus a
    Pooled row. Panel b: P(VE > OF) per individual + pooled. Arrays
    computed in ``locust_posthoc.hypothesis_elpd_individual_stats``."""
    stats = posthoc.hypothesis_elpd_individual_stats(trace_bimodal)
    n_ind = stats["n_ind"]
    elpd_f_all, elpd_p_all = stats["elpd_f_all"], stats["elpd_p_all"]
    se_f_all, se_p_all = stats["se_f_all"], stats["se_p_all"]
    p_gt_all = stats["p_gt_all"]

    n_rows = n_ind + 1
    fig, axA, axB, y_base = _indiv_figure(n_rows)
    y_of = y_base - Y_OFF
    y_ve = y_base + Y_OFF

    axA.errorbar(elpd_f_all, y_of, xerr=posthoc.Z * se_f_all, fmt="o",
                 color=UNIMODAL, ecolor=UNIMODAL, elinewidth=1, capsize=2,
                 markersize=3.5, alpha=0.95, label="Optic flow", zorder=3)
    axA.errorbar(elpd_p_all, y_ve, xerr=posthoc.Z * se_p_all, fmt="s",
                 color=BIMODAL, ecolor=BIMODAL, elinewidth=1, capsize=2,
                 markersize=3.5, alpha=0.95, label="Vec. encoding", zorder=3)
    axA.axvline(stats["logp_uniform"], ls="--", c="k", lw=1, zorder=1)
    _panel_a_frame(axA, y_base, [str(i + 1) for i in range(n_ind)] + ["Pooled"],
                   "ELPD", y_of, y_ve)
    leg = axA.legend(frameon=True, loc="upper left")
    leg.get_frame().set_facecolor("white")
    leg.get_frame().set_alpha(1.0)

    _panel_b_p_gt(axB, p_gt_all, y_base, "P(Vec. encoding > Optic flow)")

    fig.tight_layout()
    _save(fig, save, save_path)


def plot_alignment_individual(trace_bimodal, *, save=True,
                              save_path="figures/locust/locust_directional_alignment_individual_ab.pdf"):
    """Panel a: per-individual posterior mean alignment HDI (blue circle =
    OF, orange square = VE axial), zero baseline, scale -1..1. Panel b:
    P(VE > OF) per individual. Arrays computed in
    ``locust_posthoc.individual_alignment_stats``."""
    stats = posthoc.individual_alignment_stats(trace_bimodal)
    n_ind = stats["n_ind"]
    mean_f, mean_p = stats["mean_f"], stats["mean_p"]
    hdi_f, hdi_p = stats["hdi_f"], stats["hdi_p"]
    p_gt = stats["p_VE_gt_OF"]

    fig, axA, axB, y_base = _indiv_figure(n_ind)
    y_of = y_base - Y_OFF
    y_ve = y_base + Y_OFF

    axA.errorbar(mean_f, y_of, xerr=_asym_xerr(mean_f, hdi_f),
                 fmt="o", color=UNIMODAL, ecolor=UNIMODAL, elinewidth=1,
                 capsize=2, markersize=3.5, alpha=0.95, label="Optic flow", zorder=3)
    axA.errorbar(mean_p, y_ve, xerr=_asym_xerr(mean_p, hdi_p),
                 fmt="s", color=BIMODAL, ecolor=BIMODAL, elinewidth=1,
                 capsize=2, markersize=3.5, alpha=0.95, label="Vec. encoding", zorder=3)
    axA.axvline(0.0, ls="--", c="k", lw=1, zorder=1)
    _panel_a_frame(axA, y_base, [str(i + 1) for i in range(n_ind)],
                   "Posterior mean alignment", y_of, y_ve)
    axA.set_xlim(-1, 1)
    leg = axA.legend(frameon=True, loc="upper left")
    leg.get_frame().set_facecolor("white")
    leg.get_frame().set_alpha(1.0)

    _panel_b_p_gt(axB, p_gt, y_base, "P(Vec. encoding > Optic flow)")

    fig.tight_layout()
    _save(fig, save, save_path)
