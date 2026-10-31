"""Figure helpers for the zebra_escape example notebook.

The notebook keeps the statistical work (loading the block-averaged data,
fitting the von Mises GLM, and computing the per-draw per-track directions,
concentrations, KDEs, HDIs, and the bootstrap comparisons). This module holds
only the matplotlib "how": the palette, marker scheme, gridspecs, panels,
legends, and the export of the two figures. Each ``plot_*`` function takes
values the notebook has already computed.

Saved figures go to ``figures/zebra/`` next to this file, resolved from the
file location so the output does not depend on the working directory (the
notebook may run from ``notebooks/`` or from the repository root).

Palette: pre-stimulus is orange, post-stimulus is purple (colorblind-safe).
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D

FIG_DIR = Path(__file__).resolve().parent / "figures" / "zebra"

DPI = 300
TICK_LEN = 3
TICK_W = 0.6
AX_LW = 0.6
FONT_SZ = 7

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


def set_style():
    """Consistent figure style across the notebook (colorblind palette,
    Arial, 300 dpi, type-42 PDF outlines, hairline axes)."""
    mpl.rcParams.update({
        "figure.dpi": DPI,
        "savefig.dpi": DPI,
        "font.size": FONT_SZ,
        "axes.labelsize": FONT_SZ,
        "axes.titlesize": FONT_SZ,
        "legend.fontsize": FONT_SZ,
        "xtick.labelsize": FONT_SZ,
        "ytick.labelsize": FONT_SZ,
        "font.family": "Arial",
        "axes.linewidth": AX_LW,
        "xtick.major.width": TICK_W,
        "ytick.major.width": TICK_W,
        "xtick.major.size": TICK_LEN,
        "ytick.major.size": TICK_LEN,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    print(f"Figures will be saved to: {FIG_DIR}")


def _ensure_fig_dir():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    return FIG_DIR


def _restore_ticks(ax):
    ax.tick_params(
        axis="both", which="both",
        direction="out",
        bottom=True, top=False, left=True, right=False,
        length=TICK_LEN, width=TICK_W,
    )
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")


def _despine(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(AX_LW)
    _restore_ticks(ax)


# ---------------------------------------------------------------------------
# Figure 1: 2x2 posterior summary (group/individual, direction/concentration)
# ---------------------------------------------------------------------------
def plot_summary(theta,
                 r_group_pre, r_group_post,
                 r_ind_pre, r_ind_post,
                 mu_group_pre_mean, mu_group_post_mean,
                 R_group_pre_draw, R_group_post_draw,
                 emp_mu_pre, emp_mu_post,
                 R_ind_pre_mean, R_ind_post_mean,
                 R_ind_pre_hdi, R_ind_post_hdi):
    """Two-by-two posterior summary: group-level mean direction (polar KDE),
    group-level alignment (mean R with 95% HDI), individual-level mean
    directions (polar KDE per track plus empirical dots), and individual-level
    concentration (paired means with 95% HDIs). Pre is orange, post is purple.

    All geometry arrays are precomputed by the notebook; ``theta`` is the
    angular grid, the ``r_*`` are the KDE radii (1 = baseline, 2 = max), and
    the draws/HDI arrays index tracks by column and draws by row.
    """
    _ensure_fig_dir()
    nT = R_ind_pre_mean.shape[0]
    r_dot = 2.25
    rmax = 2.5
    FIG_W, FIG_H = 3.2, 4.1
    fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=DPI)
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
    ax_mu_group.plot(theta, r_group_pre, color=PRE_COLOR, lw=1.6)
    ax_mu_group.plot(theta, r_group_post, color=POST_COLOR, lw=1.6)
    ax_mu_group.scatter(mu_group_pre_mean, r_dot, s=S_GROUP_POLAR,
                        color=PRE_COLOR, marker=PRE_MARK, zorder=5)
    ax_mu_group.scatter(mu_group_post_mean, r_dot, s=S_GROUP_POLAR,
                        color=POST_COLOR, marker=POST_MARK, zorder=5)
    ax_mu_group.set_title("Group-level\nmean direction", pad=10)

    # (2) Group alignment (vertical)
    x_pre, x_post = 0, 1
    m_pre = float(np.mean(R_group_pre_draw))
    lo_pre, hi_pre = _linear_hdi(R_group_pre_draw, prob=0.95)
    m_post = float(np.mean(R_group_post_draw))
    lo_post, hi_post = _linear_hdi(R_group_post_draw, prob=0.95)

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
        ax_mu_ind.plot(theta, r_ind_pre[j], color=PRE_COLOR, lw=0.8, alpha=0.5)
        ax_mu_ind.plot(theta, r_ind_post[j], color=POST_COLOR, lw=0.8, alpha=0.5)
    for m_pre, m_post in zip(emp_mu_pre, emp_mu_post):
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
        ax_R_ind.plot([x_pre[j], x_post[j]], [R_ind_pre_mean[j], R_ind_post_mean[j]],
                      color="k", lw=0.8, alpha=1.0, zorder=0)
    for j in range(nT):
        lo, hi = R_ind_pre_hdi[j]
        ax_R_ind.vlines(x_pre[j], lo, hi, color=PRE_COLOR, lw=0.8, alpha=1.0, zorder=1)
        ax_R_ind.scatter(x_pre[j], R_ind_pre_mean[j], s=S_INDIV_POINTS,
                         color=PRE_COLOR, marker=PRE_MARK, zorder=2)
    for j in range(nT):
        lo, hi = R_ind_post_hdi[j]
        ax_R_ind.vlines(x_post[j], lo, hi, color=POST_COLOR, lw=0.8, alpha=1.0, zorder=1)
        ax_R_ind.scatter(x_post[j], R_ind_post_mean[j], s=S_INDIV_POINTS,
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

    plt.tight_layout()

    out_path = _ensure_fig_dir() / "zebra_main_summary.pdf"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", pad_inches=0.01)
    print(f"Saved figure to: {out_path}")
    plt.show()


def _linear_hdi(x, prob=0.95):
    x = np.sort(np.asarray(x, float))
    n = x.size
    if n == 0:
        return np.nan, np.nan
    k = int(np.ceil(prob * n))
    spans = x[k - 1:] - x[:n - k + 1]
    j = int(np.argmin(spans))
    return float(x[j]), float(x[j + k - 1])


# ---------------------------------------------------------------------------
# Figure 2: supplemental PPC — bootstrap (data) vs posterior (model).
# ---------------------------------------------------------------------------
def _circ_mod(x):
    return (x + np.pi) % (2 * np.pi) - np.pi


def _circ_xerr(m, lo, hi):
    w_lo = abs(float(_circ_mod(m - lo)))
    w_hi = abs(float(_circ_mod(hi - m)))
    return np.degrees(w_lo), np.degrees(w_hi)


def plot_ppc(tracks, nT,
             muG_pre_m_b, muG_pre_lo_b, muG_pre_hi_b,
             muG_post_m_b, muG_post_lo_b, muG_post_hi_b,
             muG_pre_m_p, muG_pre_lo_p, muG_pre_hi_p,
             muG_post_m_p, muG_post_lo_p, muG_post_hi_p,
             R_pre_m_b, R_pre_lo_b, R_pre_hi_b,
             R_post_m_b, R_post_lo_b, R_post_hi_b,
             R_pre_m_p, R_pre_lo_p, R_pre_hi_p,
             R_post_m_p, R_post_lo_p, R_post_hi_p,
             mu_tbl, R_tbl):
    """Supplemental PPC figure comparing a nonparametric bootstrap of the data
    against the model posterior, for group mean direction, group alignment,
    and per-track concentration, pre vs post. The notebook precomputes every
    summary (means, 95% HDIs) and the per-track tables; this function only
    lays out the panels."""
    _ensure_fig_dir()
    fig_w = 8.8
    fig_h = max(3.4, 0.28 * nT + 2.2)
    bottom_ratio = max(1.2, 0.55 + 0.05 * nT)

    fig = plt.figure(figsize=(fig_w, fig_h), dpi=DPI)
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

    _plot_group_mu_point(ax_g_mu, muG_pre_m_b, muG_pre_lo_b, muG_pre_hi_b, y_boot[0], PRE_COLOR, PRE_MARK, filled=False)
    _plot_group_mu_point(ax_g_mu, muG_pre_m_p, muG_pre_lo_p, muG_pre_hi_p, y_post[0], PRE_COLOR, PRE_MARK, filled=True)
    _plot_group_mu_point(ax_g_mu, muG_post_m_b, muG_post_lo_b, muG_post_hi_b, y_boot[1], POST_COLOR, POST_MARK, filled=False)
    _plot_group_mu_point(ax_g_mu, muG_post_m_p, muG_post_lo_p, muG_post_hi_p, y_post[1], POST_COLOR, POST_MARK, filled=True)

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
        ("Pre", PRE_COLOR, PRE_MARK, R_pre_m_b, R_pre_lo_b, R_pre_hi_b, R_pre_m_p, R_pre_lo_p, R_pre_hi_p),
        ("Post", POST_COLOR, POST_MARK, R_post_m_b, R_post_lo_b, R_post_hi_b, R_post_m_p, R_post_lo_p, R_post_hi_p),
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

    out_path = _ensure_fig_dir() / "zebra_ppc_bootstrap_vs_posterior.pdf"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", pad_inches=0.01)
    plt.show()
    print(f"Saved figure to: {out_path}")
