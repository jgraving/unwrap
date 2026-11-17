"""Figure helpers for the zebra_escape example notebook.

The notebook keeps the statistical work (loading the block-averaged data and
fitting the hierarchical von Mises GLM); the post-hoc computations live in
``zebra_posthoc`` and the matplotlib "how" (palette, marker scheme, gridspecs,
panels, legends, export) lives here. Each ``plot_*`` function takes the fit
result (and, for the PPC figure, the model) and the modelling table, calls
into ``zebra_posthoc`` for the quantities it draws, and is a single call from
the notebook.

Style and saving follow the shared house convention from ``unwrap.plots``:
figures render at display resolution (100 dpi) in the notebook cell and are
saved as 300 dpi vector files under ``<repo>/figures/zebra/``. Every public
figure function accepts a ``save`` flag (on by default) and an explicit
repo-relative ``save_path`` (the notebook always passes it); the written
files are gitignored.

Palette: pre-stimulus is orange, post-stimulus is purple (colorblind-safe).
"""
from pathlib import Path

import warnings

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import ArrowStyle, FancyArrowPatch
from matplotlib.transforms import offset_copy

from unwrap import plots as _plots
from unwrap.circular import circ_xerr_from_interval as _circ_xerr

from zebra_posthoc import ppc_bootstrap_stats, summary_stats

#: Publication figures for this study live in
#: <repo>/figures/zebra/, anchored on this file's location so the output does
#: not depend on the kernel cwd (the notebook may run from notebooks/ or from
#: the repository root).
FIG_DIR = Path(__file__).resolve().parent.parent / "figures" / "zebra"

TICK_LEN = 3
TICK_W = 0.6
AX_LW = 0.6
FONT_SZ = 7

_EXTRA_RC = {
    "font.size": FONT_SZ,
    "axes.labelsize": FONT_SZ,
    "axes.titlesize": FONT_SZ,
    "legend.fontsize": FONT_SZ,
    "xtick.labelsize": FONT_SZ,
    "ytick.labelsize": FONT_SZ,
    "axes.linewidth": AX_LW,
    "xtick.major.width": TICK_W,
    "ytick.major.width": TICK_W,
    "xtick.major.size": TICK_LEN,
    "ytick.major.size": TICK_LEN,
}


def set_style():
    """Consistent figure style across the notebook (the shared house
    convention from ``unwrap.plots`` plus the study's small type additions;
    figures under ``FIG_DIR``)."""
    _plots.set_style(FIG_DIR, _EXTRA_RC)


_CB = sns.color_palette("colorblind")
PRE_COLOR = _CB[1]  # orange
POST_COLOR = _CB[4]  # purple

# Marker scheme.
PRE_MARK = "s"  # square
POST_MARK = "o"  # circle

MARKER_SHRINK = 0.70
SCALE_S = MARKER_SHRINK**2

S_GROUP_POLAR = 18 * SCALE_S
S_GROUP_ALIGN = 22 * SCALE_S
S_EMP_POLAR = 8 * SCALE_S
S_INDIV_POINTS = 8 * SCALE_S

LEGEND_MS = 5 * MARKER_SHRINK


def _save(fig, save, save_path):
    """Write ``fig`` to ``save_path`` (a repo-relative path such as
    ``figures/zebra/x.pdf``) through the house two-res saver. The file
    format follows the path's suffix."""
    if not save:
        return
    _plots.save_by_path(fig, FIG_DIR, save_path)


def _tight_layout(fig):
    """``fig.tight_layout()`` without the benign 'Axes are not compatible'
    UserWarning: this grid carries inset/twin axes that tight_layout skips
    around, so the layout is correct and the warning is noise the published
    notebook should not surface."""
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="This figure includes Axes that are not compatible with "
                    "tight_layout",
            category=UserWarning)
        fig.tight_layout()



def _restore_ticks(ax):
    _plots.tick_style(ax, width=TICK_W, length=TICK_LEN)


def _despine(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(AX_LW)
    _restore_ticks(ax)


# ---------------------------------------------------------------------------
# Figure 1: 2x2 posterior summary (group/individual, direction/concentration)
# ---------------------------------------------------------------------------
def plot_summary(idata, df_blocks, *, save=True,
                 save_path="figures/zebra/zebra_main_summary.pdf"):
    """Two-by-two posterior summary: group-level mean direction (polar KDE),
    group-level alignment (mean R with 95% HDI), individual-level mean
    directions (polar KDE per track plus empirical dots), and individual-level
    concentration (paired means with 95% HDIs). Pre is orange, post is
    purple.

    The per-draw per-track quantities, KDE radii, alignment, concentration
    means and HDIs are all computed here from the fit result via
    :func:`zebra_posthoc.summary_stats`, so the notebook cell is a single
    call.
    """
    prep = summary_stats(idata, df_blocks)
    nT = prep["nT"]
    theta = prep["theta"]
    r_dot = 2.25
    rmax = 2.5
    FIG_W, FIG_H = 3.2, 4.1
    fig = plt.figure(figsize=(FIG_W, FIG_H))
    gs = fig.add_gridspec(2, 2, wspace=0.7, hspace=0.70)

    ax_mu_group = fig.add_subplot(gs[0, 0], projection="polar")
    ax_R_group = fig.add_subplot(gs[0, 1])
    ax_mu_ind = fig.add_subplot(gs[1, 0], projection="polar")
    ax_R_ind = fig.add_subplot(gs[1, 1])

    # --- polar formatting ---
    def style_polar(ax):
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set_rlim(0, rmax)
        ax.set_rticks([1.0, 2.0])
        ax.set_yticklabels(["0", "1.0"], color="black")
        ax.grid(True, color="0.25", linewidth=0.55, alpha=0.65)
        ax.spines["polar"].set_visible(False)
        ax.tick_params(axis="x", which="major", labelsize=7, pad=0, colors="black")
        ax.tick_params(axis="y", which="major", labelsize=6, pad=0, colors="black")

    style_polar(ax_mu_group)
    style_polar(ax_mu_ind)

    # (1) Group-level mean direction (polar)
    ax_mu_group.plot(theta, prep["r_group_pre"], color=PRE_COLOR, lw=1.6)
    ax_mu_group.plot(theta, prep["r_group_post"], color=POST_COLOR, lw=1.6)
    ax_mu_group.scatter(prep["mu_group_pre_mean"], r_dot, s=S_GROUP_POLAR,
                        color=PRE_COLOR, marker=PRE_MARK, zorder=5)
    ax_mu_group.scatter(prep["mu_group_post_mean"], r_dot, s=S_GROUP_POLAR,
                        color=POST_COLOR, marker=POST_MARK, zorder=5)
    ax_mu_group.set_title("Group-level\nmean direction", pad=10)

    # (2) Group alignment (vertical)
    x_pre, x_post = 0, 1
    m_pre = prep["R_group_pre_mean"]
    lo_pre, hi_pre = prep["R_group_pre_lo"], prep["R_group_pre_hi"]
    m_post = prep["R_group_post_mean"]
    lo_post, hi_post = prep["R_group_post_lo"], prep["R_group_post_hi"]

    ax_R_group.vlines(x_pre, lo_pre, hi_pre, color=PRE_COLOR, lw=1.4, alpha=0.95, zorder=2)
    ax_R_group.vlines(x_post, lo_post, hi_post, color=POST_COLOR, lw=1.4, alpha=0.95, zorder=2)
    ax_R_group.scatter([x_pre], [m_pre], s=S_GROUP_ALIGN, color=PRE_COLOR, marker=PRE_MARK, zorder=3)
    ax_R_group.scatter([x_post], [m_post], s=S_GROUP_ALIGN, color=POST_COLOR, marker=POST_MARK, zorder=3)

    ax_R_group.set_xlim(-0.5, 1.5)
    ax_R_group.set_xticks([x_pre, x_post])
    ax_R_group.set_xticklabels(["Pre", "Post"])
    ax_R_group.set_ylim(0, 1.0)
    ax_R_group.set_yticks([0.0, 0.5, 1.0])
    ax_R_group.set_ylabel("Mean resultant length")
    ax_R_group.set_title("Group-level alignment", pad=20)
    _despine(ax_R_group)
    _restore_ticks(ax_R_group)

    # (3) Individual-level mean directions (polar)
    for j in range(nT):
        ax_mu_ind.plot(theta, prep["r_ind_pre"][j], color=PRE_COLOR, lw=0.8, alpha=0.5)
        ax_mu_ind.plot(theta, prep["r_ind_post"][j], color=POST_COLOR, lw=0.8, alpha=0.5)
    for m_pre, m_post in zip(prep["emp_mu_pre"], prep["emp_mu_post"]):
        if not np.isnan(m_pre):
            ax_mu_ind.scatter(m_pre, r_dot, s=S_EMP_POLAR, color=PRE_COLOR,
                              marker=PRE_MARK, alpha=1.0, zorder=5)
        if not np.isnan(m_post):
            ax_mu_ind.scatter(m_post, r_dot, s=S_EMP_POLAR, color=POST_COLOR,
                              marker=POST_MARK, alpha=1.0, zorder=5)
    ax_mu_ind.set_title("Individual-level\nmean direction", pad=10)

    # (4) Individual concentration (paired means + 95% HDIs)
    rng = np.random.default_rng(123)
    delta = rng.normal(0.0, 0.15, size=nT)
    x_pre = 0.0 + delta
    x_post = 1.0 + delta
    for j in range(nT):
        ax_R_ind.plot([x_pre[j], x_post[j]], [prep["R_ind_pre_mean"][j], prep["R_ind_post_mean"][j]],
                      color="k", lw=0.8, alpha=1.0, zorder=0)
    for j in range(nT):
        lo, hi = prep["R_ind_pre_hdi"][j]
        ax_R_ind.vlines(x_pre[j], lo, hi, color=PRE_COLOR, lw=0.8, alpha=1.0, zorder=1)
        ax_R_ind.scatter(x_pre[j], prep["R_ind_pre_mean"][j], s=S_INDIV_POINTS,
                         color=PRE_COLOR, marker=PRE_MARK, zorder=2)
    for j in range(nT):
        lo, hi = prep["R_ind_post_hdi"][j]
        ax_R_ind.vlines(x_post[j], lo, hi, color=POST_COLOR, lw=0.8, alpha=1.0, zorder=1)
        ax_R_ind.scatter(x_post[j], prep["R_ind_post_mean"][j], s=S_INDIV_POINTS,
                         color=POST_COLOR, marker=POST_MARK, zorder=2)

    ax_R_ind.set_xlim(-0.5, 1.5)
    ax_R_ind.set_xticks([0, 1])
    ax_R_ind.set_xticklabels(["Pre", "Post"])
    ax_R_ind.set_ylim(0, 1.0)
    ax_R_ind.set_yticks([0.0, 0.5, 1.0])
    ax_R_ind.set_ylabel("Mean resultant length")
    ax_R_ind.set_title("Individual-level\nconcentration", pad=15)
    _despine(ax_R_ind)

    _restore_ticks(ax_R_group)
    _restore_ticks(ax_R_ind)

    # Legend (compact, no overlap)
    handles = [
        Line2D([0], [0], color=PRE_COLOR, lw=1.6, marker=PRE_MARK,
               markersize=LEGEND_MS, label="Pre"),
        Line2D([0], [0], color=POST_COLOR, lw=1.6, marker=POST_MARK,
               markersize=LEGEND_MS, label="Post"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=True,
               bbox_to_anchor=(0.5, 0.01), columnspacing=1.2, handletextpad=0.6)

    _tight_layout(fig)
    _save(fig, save, save_path)
    return fig


# ---------------------------------------------------------------------------
# Figure 2: supplemental PPC - bootstrap (data) vs posterior (model).
# ---------------------------------------------------------------------------


def plot_ppc(idata, model, df_blocks, *, save=True,
             save_path="figures/zebra/zebra_ppc_bootstrap_vs_posterior.pdf"):
    """Supplemental PPC figure comparing a nonparametric bootstrap of the data
    against the model posterior, for group mean direction, group alignment,
    and per-track mean direction and concentration, pre vs post.

    The bootstrap resampling, the per-observation response parameters (read
    from the PyMC graph, see :func:`zebra_posthoc.ppc_bootstrap_stats`) and
    every 95% HDI summary are computed here, so the notebook cell is a
    single call.
    """
    prep = ppc_bootstrap_stats(idata, model, df_blocks)
    tracks = prep["tracks"]
    nT = prep["nT"]
    mu_tbl = prep["mu_tbl"]
    R_tbl = prep["R_tbl"]

    fig_w = 8.8
    fig_h = max(3.4, 0.28 * nT + 2.2)
    bottom_ratio = max(1.2, 0.55 + 0.05 * nT)

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = fig.add_gridspec(2, 4, height_ratios=[0.55, bottom_ratio],
                          wspace=0.45, hspace=0.38)

    ax_g_mu = fig.add_subplot(gs[0, 0:2])
    ax_g_R = fig.add_subplot(gs[0, 2:4])
    ax_mu_pre = fig.add_subplot(gs[1, 0])
    ax_mu_post = fig.add_subplot(gs[1, 1])
    ax_R_pre = fig.add_subplot(gs[1, 2])
    ax_R_post = fig.add_subplot(gs[1, 3])

    # panel labels (a-f): bold Arial, readable size
    PANEL_FS = 11

    def _panel_label(ax, lab, xoff=0.05):
        ax.text(xoff, 1.0, lab, transform=ax.transAxes, ha="left", va="top",
                fontsize=PANEL_FS, fontweight="bold", fontfamily="Arial",
                clip_on=False, zorder=10)

    _panel_label(ax_g_mu, "a.", 0.02)
    _panel_label(ax_g_R, "b.", 0.02)
    _panel_label(ax_mu_pre, "c.")
    _panel_label(ax_mu_post, "d.")
    _panel_label(ax_R_pre, "e.")
    _panel_label(ax_R_post, "f.")

    y_track = np.arange(nT, dtype=float)
    dy_src = 0.14

    def add_row_boxes(ax):
        ax.set_axisbelow(True)
        for y in np.arange(-0.5, nT, 1.0):
            ax.axhline(y, color="0.90", lw=0.5, zorder=0)

    for ax in (ax_mu_pre, ax_mu_post, ax_R_pre, ax_R_post):
        add_row_boxes(ax)

    # (1) group mean direction (deg)
    y_base = np.array([0.0, 1.0])  # Pre, Post
    y_boot = y_base - 0.12
    y_post = y_base + 0.12

    def _plot_group_mu_point(ax, m, lo, hi, y, color, marker, filled):
        wlo, whi = _circ_xerr(m, lo, hi)
        ax.errorbar(np.degrees(m), y, xerr=np.array([[wlo], [whi]]), fmt=marker,
                    color=color, ecolor=color, capsize=2, ms=3.2,
                    mfc=(color if filled else "white"), mec=color, mew=0.9,
                    lw=(1.2 if filled else 1.0), zorder=3)

    _plot_group_mu_point(ax_g_mu, prep["muG_pre_m_b"], prep["muG_pre_lo_b"], prep["muG_pre_hi_b"], y_boot[0], PRE_COLOR, PRE_MARK, filled=False)
    _plot_group_mu_point(ax_g_mu, prep["muG_pre_m_p"], prep["muG_pre_lo_p"], prep["muG_pre_hi_p"], y_post[0], PRE_COLOR, PRE_MARK, filled=True)
    _plot_group_mu_point(ax_g_mu, prep["muG_post_m_b"], prep["muG_post_lo_b"], prep["muG_post_hi_b"], y_boot[1], POST_COLOR, POST_MARK, filled=False)
    _plot_group_mu_point(ax_g_mu, prep["muG_post_m_p"], prep["muG_post_lo_p"], prep["muG_post_hi_p"], y_post[1], POST_COLOR, POST_MARK, filled=True)

    ax_g_mu.set_yticks(y_base)
    ax_g_mu.set_yticklabels(["Pre", "Post"])
    ax_g_mu.set_xlim(-180, 180)
    ax_g_mu.set_xticks([-180, -90, 0, 90, 180])
    ax_g_mu.set_xlabel("Group mean direction (deg)")
    ax_g_mu.set_title("Group mean direction")
    ax_g_mu.spines["top"].set_visible(False)
    ax_g_mu.spines["right"].set_visible(False)
    _restore_ticks(ax_g_mu)

    # (2) group alignment (R): Pre/Post
    rows = [
        ("Pre", PRE_COLOR, PRE_MARK, prep["R_pre_m_b"], prep["R_pre_lo_b"], prep["R_pre_hi_b"], prep["R_pre_m_p"], prep["R_pre_lo_p"], prep["R_pre_hi_p"]),
        ("Post", POST_COLOR, POST_MARK, prep["R_post_m_b"], prep["R_post_lo_b"], prep["R_post_hi_b"], prep["R_post_m_p"], prep["R_post_lo_p"], prep["R_post_hi_p"]),
    ]
    y_base2 = np.arange(len(rows), dtype=float)
    y_boot2 = y_base2 - 0.12
    y_post2 = y_base2 + 0.12

    for i, (lab, col, mark, mb, lob, hib, mp, lop, hip) in enumerate(rows):
        if np.isfinite(mb):
            ax_g_R.errorbar(mb, y_boot2[i], xerr=np.array([[mb - lob], [hib - mb]]),
                            fmt=mark, color=col, ecolor=col, capsize=2, ms=3.2,
                            mfc="white", mec=col, mew=0.9, lw=1.0, zorder=3)
        if np.isfinite(mp):
            ax_g_R.errorbar(mp, y_post2[i], xerr=np.array([[mp - lop], [hip - mp]]),
                            fmt=mark, color=col, ecolor=col, capsize=2, ms=3.2,
                            mfc=col, mec=col, mew=0.9, lw=1.2, zorder=3)

    ax_g_R.set_yticks(y_base2)
    ax_g_R.set_yticklabels([r[0] for r in rows])
    ax_g_R.set_xlim(0, 1.0)
    ax_g_R.set_xticks([0.0, 0.5, 1.0])
    ax_g_R.set_xlabel("Mean resultant length")
    ax_g_R.set_title("Group alignment")
    ax_g_R.spines["top"].set_visible(False)
    ax_g_R.spines["right"].set_visible(False)
    _restore_ticks(ax_g_R)

    # bottom panels: one condition per axis
    def _plot_mu_condition(ax, cond, cond_color, cond_marker):
        sub_b = mu_tbl[(mu_tbl["Condition"] == cond) & (mu_tbl["Source"] == "bootstrap")].set_index("Track")
        sub_p = mu_tbl[(mu_tbl["Condition"] == cond) & (mu_tbl["Source"] == "posterior")].set_index("Track")
        for j, t in enumerate(tracks):
            if t in sub_b.index:
                m = float(sub_b.loc[t, "mean"])
                lo = float(sub_b.loc[t, "low"])
                hi = float(sub_b.loc[t, "high"])
                wlo, whi = _circ_xerr(m, lo, hi)
                ax.errorbar(np.degrees(m), y_track[j] - dy_src, xerr=np.array([[wlo], [whi]]),
                            fmt=cond_marker, color=cond_color, ecolor=cond_color,
                            capsize=2, ms=2.7, mfc="white", mec=cond_color, mew=0.8, lw=0.9, zorder=3)
            if t in sub_p.index:
                m = float(sub_p.loc[t, "mean"])
                lo = float(sub_p.loc[t, "low"])
                hi = float(sub_p.loc[t, "high"])
                wlo, whi = _circ_xerr(m, lo, hi)
                ax.errorbar(np.degrees(m), y_track[j] + dy_src, xerr=np.array([[wlo], [whi]]),
                            fmt=cond_marker, color=cond_color, ecolor=cond_color,
                            capsize=2, ms=2.7, mfc=cond_color, mec=cond_color, mew=0.8, lw=1.0, zorder=3)

    def _plot_R_condition(ax, cond, cond_color, cond_marker):
        sub_b = R_tbl[(R_tbl["Condition"] == cond) & (R_tbl["Source"] == "bootstrap")].set_index("Track")
        sub_p = R_tbl[(R_tbl["Condition"] == cond) & (R_tbl["Source"] == "posterior")].set_index("Track")
        for j, t in enumerate(tracks):
            if t in sub_b.index:
                m = float(sub_b.loc[t, "mean"])
                lo = float(sub_b.loc[t, "low"])
                hi = float(sub_b.loc[t, "high"])
                ax.errorbar(m, y_track[j] - dy_src, xerr=np.array([[m - lo], [hi - m]]),
                            fmt=cond_marker, color=cond_color, ecolor=cond_color,
                            capsize=2, ms=2.7, mfc="white", mec=cond_color, mew=0.8, lw=0.9, zorder=3)
            if t in sub_p.index:
                m = float(sub_p.loc[t, "mean"])
                lo = float(sub_p.loc[t, "low"])
                hi = float(sub_p.loc[t, "high"])
                ax.errorbar(m, y_track[j] + dy_src, xerr=np.array([[m - lo], [hi - m]]),
                            fmt=cond_marker, color=cond_color, ecolor=cond_color,
                            capsize=2, ms=2.7, mfc=cond_color, mec=cond_color, mew=0.8, lw=1.0, zorder=3)

    _plot_mu_condition(ax_mu_pre, "Pre", PRE_COLOR, PRE_MARK)
    _plot_mu_condition(ax_mu_post, "Post", POST_COLOR, POST_MARK)
    _plot_R_condition(ax_R_pre, "Pre", PRE_COLOR, PRE_MARK)
    _plot_R_condition(ax_R_post, "Post", POST_COLOR, POST_MARK)

    for ax in (ax_mu_pre, ax_mu_post, ax_R_pre, ax_R_post):
        ax.set_yticks(y_track)
        ax.invert_yaxis()
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        _restore_ticks(ax)

    ax_mu_pre.set_yticklabels(tracks)
    ax_mu_post.set_yticklabels([""] * nT)
    ax_R_pre.set_yticklabels([""] * nT)
    ax_R_post.set_yticklabels([""] * nT)

    ax_mu_pre.set_xlim(-220, 220)
    ax_mu_post.set_xlim(-220, 220)
    ax_mu_pre.set_xticks([-180, -90, 0, 90, 180])
    ax_mu_post.set_xticks([-180, -90, 0, 90, 180])

    ax_R_pre.set_xlim(0, 1.0)
    ax_R_post.set_xlim(0, 1.0)
    ax_R_pre.set_xticks([0.0, 0.5, 1.0])
    ax_R_post.set_xticks([0.0, 0.5, 1.0])

    def _pair_xlabel(fig, ax_left, ax_right, text, y_pad=0.045):
        pL = ax_left.get_position()
        pR = ax_right.get_position()
        x = 0.5 * (pL.x0 + pR.x1)
        y = min(pL.y0, pR.y0) - y_pad
        fig.text(x, y, text, ha="center", va="top")

    ax_mu_pre.set_xlabel("")
    ax_mu_post.set_xlabel("")
    ax_R_pre.set_xlabel("")
    ax_R_post.set_xlabel("")

    _pair_xlabel(fig, ax_mu_pre, ax_mu_post, "Mean direction (deg)", y_pad=0.045)
    _pair_xlabel(fig, ax_R_pre, ax_R_post, "Mean resultant length", y_pad=0.045)

    def _pair_label(ax_left, ax_right, text, y_pad=0.006):
        pL = ax_left.get_position()
        pR = ax_right.get_position()
        x = 0.5 * (pL.x0 + pR.x1)
        y = max(pL.y1, pR.y1) + y_pad
        fig.text(x, y, text, ha="center", va="bottom")

    _pair_label(ax_mu_pre, ax_mu_post, "Individual mean direction")
    _pair_label(ax_R_pre, ax_R_post, "Individual concentration")

    legend_handles = [
        Line2D([0], [0], marker="o", linestyle="none", color="k",
               markerfacecolor="white", markeredgecolor="k", markeredgewidth=0.9,
               markersize=4.2, label="Bootstrap (data)"),
        Line2D([0], [0], marker="o", linestyle="none", color="k",
               markerfacecolor="k", markeredgecolor="k", markeredgewidth=0.9,
               markersize=4.2, label="Posterior (model)"),
        Line2D([0], [0], color=PRE_COLOR, lw=1.5, marker=PRE_MARK, markersize=5, label="Pre"),
        Line2D([0], [0], color=POST_COLOR, lw=1.5, marker=POST_MARK, markersize=5, label="Post"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=4, frameon=True,
               bbox_to_anchor=(0.5, -0.01), columnspacing=1.2, handletextpad=0.6)

    _save(fig, save, save_path)
    return fig


# ---------------------------------------------------------------------------
# Escape-response map figure (the map section of zebra_escape.ipynb).
#
# Group trajectories plus the focal animal's body-orientation arrows drawn
# on the aerial orthomosaic crop. Data comes from the loaders in
# ``unwrap.datasets`` (``load_zebra_mosaic`` for the backdrop,
# ``load_zebra_postures`` for positions, ``load_zebra_bearings_raw`` for
# the per-frame bearings); everything is plotted as recorded -- no
# smoothing anywhere. Arrow fills follow the house pre/post palette
# (orange pre-stimulus, purple post-stimulus).
# ---------------------------------------------------------------------------


def _save_map(fig, save, save_path):
    """Write ``fig`` under ``figures/zebra/`` (house two-res saver)."""
    if not save:
        return
    _plots.save_by_path(fig, FIG_DIR, save_path)


def frames_every_distance_greedy(x, y, spacing_m, *, start_offset_m=0.0):
    """Pick frame indices at ~fixed straight-line spacing (greedy walk).

    Start at the first valid (non-NaN) sample; walk forward and place a
    point whenever the straight-line distance from the last anchor exceeds
    the threshold -- ``start_offset_m`` for the first placement,
    ``spacing_m`` thereafter. A NaN gap re-anchors at the next valid frame.
    Returns indices into ``x``/``y``.
    """
    x = np.asarray(x)
    y = np.asarray(y)
    assert x.ndim == y.ndim == 1 and x.shape == y.shape
    spacing_m = float(spacing_m)
    start_offset_m = float(start_offset_m)
    if spacing_m <= 0:
        raise ValueError("spacing_m must be > 0")

    n = x.size
    valid = np.isfinite(x) & np.isfinite(y)
    out = []

    i = 0
    while i < n and not valid[i]:
        i += 1
    if i >= n:
        return np.array([], dtype=int)

    anchor = i
    threshold = start_offset_m
    placed_any = False

    i += 1
    while i < n:
        if not valid[i]:
            i += 1
            while i < n and not valid[i]:
                i += 1
            if i >= n:
                break
            anchor = i
            threshold = start_offset_m if not placed_any else spacing_m
            i += 1
            continue

        dx = x[i] - x[anchor]
        dy = y[i] - y[anchor]
        dist = (dx * dx + dy * dy) ** 0.5
        if dist >= threshold:
            out.append(i)
            anchor = i
            placed_any = True
            threshold = spacing_m
            i += 1
        else:
            i += 1

    return np.array(out, dtype=int)


def add_scale_bar_dual(
    ax,
    total_length_m=20.0,
    n_segments=2,
    segment_colors=("white", "black"),
    segment_height_pt=8,
    edgecolor="black",
    edgewidth=0.8,
    label=True,
    label_above=False,
    label_kwargs=None,
    ax_xy=(0.95, 0.08),
    offset_pt=(0, -41),
    align="right",
    direction="left",
    pad_text_pt=4,
):
    """Two-tone horizontal scale bar with a length in METERS (data units).

    Placement uses axes coords plus a point offset; the bar is drawn as
    rectangles in data coordinates.
    """
    if label_kwargs is None:
        label_kwargs = dict(color="black", fontsize=10, ha="center", va="top")

    base = offset_copy(ax.transAxes, fig=ax.figure, x=offset_pt[0],
                       y=offset_pt[1], units="points")
    x_disp, y_disp = base.transform(ax_xy)
    x_anchor, y_anchor = ax.transData.inverted().transform((x_disp, y_disp))

    seg_len = total_length_m / float(max(n_segments, 1))
    if direction == "left":
        x0 = x_anchor if align == "left" else x_anchor - total_length_m
    elif direction == "right":
        if align == "right":
            x0 = x_anchor - total_length_m
        else:
            x0 = x_anchor
    else:
        raise ValueError("direction must be 'left' or 'right'")
    xdir = +1

    # height in points -> data units
    fig_dpi = ax.figure.dpi
    _, axpix_h = ax.figure.get_size_inches() * fig_dpi
    y_min, y_max = ax.get_ylim()
    pt_to_data_y = (y_max - y_min) / axpix_h
    half_h = 0.5 * segment_height_pt * pt_to_data_y
    yb, yt = y_anchor - half_h, y_anchor + half_h

    colors = (segment_colors *
              ((n_segments + len(segment_colors) - 1) // len(segment_colors))
              )[:n_segments]
    for i in range(n_segments):
        xa = x0 + xdir * (i * seg_len)
        xb = x0 + xdir * ((i + 1) * seg_len)
        ax.add_patch(plt.Rectangle(
            (min(xa, xb), yb), abs(xb - xa), yt - yb,
            facecolor=colors[i], edgecolor=edgecolor, linewidth=edgewidth,
            zorder=5))

    if label:
        x_center = x0 + xdir * (total_length_m / 2.0)
        txt = ax.text(x_center, y_anchor, f"{int(total_length_m)} m",
                      **label_kwargs)
        dy = pad_text_pt if label_above else -pad_text_pt
        txt.set_transform(offset_copy(ax.transData, fig=ax.figure, x=0,
                                      y=dy, units="points"))


def add_compass(
    ax,
    ax_xy=(0.92, 0.18),
    offset_pt=(0, -95),
    arrow_len_pt=42,
    head_width_pt=12,
    head_length_pt=14,
    tail_width_pt=2.5,
    color="black",
    edgecolor="black",
    edgewidth=1.0,
    alpha=1.0,
    draw_N=True,
    N_fontsize=12,
    N_offset_pt=(0, 4),
):
    """North arrow drawn in AXES coordinates (stays upright, fixed size)."""
    base = offset_copy(ax.transAxes, fig=ax.figure, x=offset_pt[0],
                       y=offset_pt[1], units="points")
    fig_dpi = ax.figure.dpi
    _, axpix_h = ax.figure.get_size_inches() * fig_dpi
    dy_axes = arrow_len_pt / axpix_h

    x_axes, y_axes = ax_xy
    style = ArrowStyle.Simple(head_length=head_length_pt,
                              head_width=head_width_pt,
                              tail_width=tail_width_pt)
    patch = FancyArrowPatch((x_axes, y_axes), (x_axes, y_axes + dy_axes),
                            arrowstyle=style, mutation_scale=1.0,
                            linewidth=edgewidth, edgecolor=edgecolor,
                            facecolor=color, alpha=alpha, joinstyle="miter",
                            capstyle="butt", zorder=6)
    patch.set_transform(base)
    ax.add_patch(patch)

    if draw_N:
        tx = ax.text(x_axes, y_axes + dy_axes, "N", transform=ax.transAxes,
                     ha="center", va="bottom", fontsize=N_fontsize,
                     color=color, zorder=6)
        tx.set_transform(offset_copy(
            ax.transAxes, fig=ax.figure,
            x=N_offset_pt[0] + offset_pt[0], y=N_offset_pt[1] + offset_pt[1],
            units="points"))


def plot_escape_map(
    img,
    extent,
    pos,
    bearings,
    focal_track,
    *,
    other_color="white", other_alpha=0.7, other_lw=0.8, other_stride=1,
    focal_line_color="black", focal_lw=2.2, focal_stride=1,
    arrow_mode="distance",
    fps=30.0,
    arrow_every_s=20.0,
    arrows_every_m=5.0,
    arrows_start_offset_m=0.0,
    arrow_length=12.0,
    tail_width_pt=1.2, head_width_pt=6.0, head_length_pt=7.5,
    arrow_pre_color=None, arrow_post_color=None, arrow_alpha=0.85,
    arrow_edge_color="black", arrow_edge_width=0.7, arrow_edge_alpha=0.95,
    arrow_zorder=4.0,
    frame_max=None,
    xlim=None, ylim=None, title=None,
    figsize=(10, 10),
    save=True,
    save_path="figures/zebra/escape_map.pdf",
):
    """Draw one escape-response map panel on the aerial mosaic.

    Parameters
    ----------
    img, extent
        RGBA mosaic crop (``H, W, 4`` uint8) and its ``(left, right,
        bottom, top)`` extent in UTM 37N meters, from
        :func:`unwrap.datasets.load_zebra_mosaic`.
    pos
        Posture positions from :func:`unwrap.datasets.load_zebra_postures`
        (columns ``frame``, ``track_number``, ``kp0_x``, ``kp0_y``, ...);
        keypoint 0 (the animal centre) is drawn as the trajectory.
    bearings
        Per-frame bearings from :func:`unwrap.datasets.load_zebra_bearings_raw`
        (columns ``track_number``, ``frame``, ``compass_bearing``,
        ``pre_scare``).
    focal_track
        ``track_number`` to draw black with orientation arrows; all other
        tracks in the mosaic window become the faint background group.
    arrow_mode
        ``"distance"`` places arrows every ``arrows_every_m`` meters of
        straight-line travel; ``"time"`` every ``arrow_every_s`` seconds at
        ``fps`` video frames per second.
    frame_max
        Optional inclusive frame cap (recorded footage extends slightly
        past the usable tracking, e.g. the published panels use 90698).

    Arrow fills are the house pre/post colors unless overridden.
    """
    if arrow_pre_color is None:
        arrow_pre_color = PRE_COLOR
    if arrow_post_color is None:
        arrow_post_color = POST_COLOR

    def _line_segments(x, y):
        x = np.asarray(x)
        y = np.asarray(y)
        keep = np.isfinite(x) & np.isfinite(y)
        x, y = x[keep], y[keep]
        if x.size < 2:
            return np.zeros((0, 2, 2))
        pts = np.column_stack([x, y])
        return np.stack([pts[:-1], pts[1:]], axis=1)

    left, right, bottom, top = extent

    # ---- focal positions, indexed by frame ----
    focal = pos.loc[pos["track_number"] == focal_track].sort_values("frame")
    focal_x = focal["kp0_x"].to_numpy(float)
    focal_y = focal["kp0_y"].to_numpy(float)
    focal_frame = focal["frame"].to_numpy(int)
    if frame_max is not None:
        clip = focal_frame <= int(frame_max)
        focal_x, focal_y, focal_frame = focal_x[clip], focal_y[clip], focal_frame[clip]

    fig, ax = plt.subplots(figsize=figsize)
    ax.imshow(img, extent=(left, right, bottom, top), origin="upper")

    # ---- non-focal trajectories (faint context) ----
    other_stride = max(int(other_stride), 1)
    segs = []
    for tg, g in pos.groupby("track_number", sort=False):
        if tg == focal_track:
            continue
        g = g.sort_values("frame")
        if frame_max is not None:
            g = g[g["frame"] <= int(frame_max)]
        gx = g["kp0_x"].to_numpy(float)[::other_stride]
        gy = g["kp0_y"].to_numpy(float)[::other_stride]
        s = _line_segments(gx, gy)
        if s.size:
            segs.append(s)
    if segs:
        ax.add_collection(LineCollection(
            np.concatenate(segs, axis=0), colors=other_color,
            linewidths=other_lw, alpha=other_alpha, zorder=2))

    # ---- focal trajectory ----
    focal_stride = max(int(focal_stride), 1)
    seg_focal = _line_segments(focal_x[::focal_stride], focal_y[::focal_stride])
    if seg_focal.size:
        ax.add_collection(LineCollection(
            seg_focal, colors=focal_line_color, linewidths=focal_lw, zorder=3))

    # ---- focal bearings ----
    b = bearings.loc[bearings["track_number"] == focal_track]
    b = b[["frame", "compass_bearing", "pre_scare"]].dropna(subset=["compass_bearing"])
    if frame_max is not None:
        b = b[b["frame"] <= int(frame_max)]
    b = b.sort_values("frame").set_index("frame")

    if b.empty:
        _finish_escape_map(ax, extent, xlim, ylim, title, focal_line_color,
                           focal_lw, other_color, arrow_pre_color, arrow_post_color)
        _save_map(fig, save, save_path)
        return fig

    # positions per frame for arrow anchors (frame -> row lookup)
    frame_to_row = {int(fr): r for r, fr in enumerate(focal_frame)}

    if str(arrow_mode).lower() == "time":
        step = max(int(round(arrow_every_s * fps)), 1)
        arrow_frames = b.index.to_numpy()[::step]
    elif str(arrow_mode).lower() == "distance":
        cand = frames_every_distance_greedy(
            focal_x, focal_y, spacing_m=float(arrows_every_m),
            start_offset_m=float(arrows_start_offset_m))
        arrow_frames = b.index.intersection(cand).to_numpy()
    else:
        raise ValueError("arrow_mode must be 'time' or 'distance'")

    if arrow_frames.size:
        bearings_deg = b.loc[arrow_frames, "compass_bearing"].to_numpy(float)
        pre_at = b.loc[arrow_frames, "pre_scare"].to_numpy().astype(bool)
        theta = np.deg2rad((90.0 - bearings_deg) % 360.0)
        ux, vy = np.cos(theta), np.sin(theta)
        n = np.hypot(ux, vy)
        n[n == 0] = 1.0
        ux, vy = ux / n, vy / n

        arrow_style = ArrowStyle.Simple(head_length=head_length_pt,
                                        head_width=head_width_pt,
                                        tail_width=tail_width_pt)
        from matplotlib.colors import to_rgba

        for fr, u, v, is_pre in zip(arrow_frames, ux, vy, pre_at):
            r = frame_to_row.get(int(fr))
            if r is None:
                continue
            x0, y0 = focal_x[r], focal_y[r]
            patch = FancyArrowPatch(
                (x0, y0), (x0 + u * arrow_length, y0 + v * arrow_length),
                arrowstyle=arrow_style, mutation_scale=1.0,
                linewidth=arrow_edge_width, edgecolor=arrow_edge_color,
                facecolor=arrow_pre_color if is_pre else arrow_post_color,
                alpha=arrow_alpha, zorder=arrow_zorder, joinstyle="miter",
                capstyle="butt")
            if arrow_edge_alpha is not None and arrow_edge_alpha != arrow_alpha:
                patch.set_edgecolor(to_rgba(arrow_edge_color, arrow_edge_alpha))
            ax.add_patch(patch)

    _finish_escape_map(ax, extent, xlim, ylim, title, focal_line_color,
                       focal_lw, other_color, arrow_pre_color, arrow_post_color)
    _save_map(fig, save, save_path)
    return fig


def _finish_escape_map(ax, extent, xlim, ylim, title, focal_line_color,
                       focal_lw, other_color, arrow_pre_color, arrow_post_color):
    """Equal aspect, panel limits, legend, scale bar and compass."""
    ax.set_aspect("equal", adjustable="box")
    if title:
        ax.set_title(title)
    ax.set_xlim(sorted(xlim) if xlim is not None else (extent[0], extent[1]))
    ax.set_ylim(sorted(ylim) if ylim is not None else (extent[2], extent[3]))

    handles = [
        Line2D([], [], color=focal_line_color, lw=focal_lw, label="Focal"),
        Line2D([], [], color=str(other_color), lw=2.4, alpha=1.0,
               label="Non-focal"),
        Line2D([], [], color=arrow_pre_color, marker=">", lw=0, ms=8,
               label="Pre-stimulus"),
        Line2D([], [], color=arrow_post_color, marker=">", lw=0, ms=8,
               label="Post-stimulus"),
    ]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(0.02, 0.5),
              frameon=True, framealpha=0.9, facecolor="lightgrey")

    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    ax.set_yticks([])

    add_scale_bar_dual(ax, total_length_m=20, n_segments=2,
                       segment_colors=("white", "black"), segment_height_pt=8,
                       edgewidth=0.8, ax_xy=(0.97, 0.10), align="right",
                       direction="left", label=True, label_above=False)
    add_compass(ax, ax_xy=(0.92, 0.18), arrow_len_pt=42, head_width_pt=12,
                head_length_pt=14, tail_width_pt=2.5, color="black",
                edgecolor="black", edgewidth=1.0, draw_N=True, N_fontsize=12,
                N_offset_pt=(0, 4))
    plt.tight_layout()
