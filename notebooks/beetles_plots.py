"""Figure helpers for the beetles example notebook.

The notebook keeps the statistical work (loading the data, fitting the von Mises
GLM with ``bambi``/``unwrap``, and computing every directional quantity). This
module holds only the matplotlib "how": the palette, polar grids, KDE rings,
error bars, layout, and the export of the figures. Each ``plot_*`` function
takes values the notebook has already computed.

Saved figures go to ``figures/beetles/`` next to this file, resolved from the
file location so the output does not depend on the working directory (the
notebook may run from ``notebooks/`` or from the repository root).

Palette: light OFF is grey, light ON is gold (matching the paper figures).
"""
from pathlib import Path
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import arviz as az
from scipy.special import i0, i1, ive
from scipy.stats import circmean, circstd, vonmises
from matplotlib.lines import Line2D

FIG_DIR = Path(__file__).resolve().parent / "figures" / "beetles"

# Light OFF is grey, light ON is gold (matches the paper figures).
OFF = "#9e9e9e"
ON = "#f4a300"
DPI = 300
TICK_LEN = 3
TICK_W = 0.6
AX_LW = 0.6


def set_style():
    """Consistent figure style across the notebook (Arial, 300 dpi, type-42
    PDF outlines, hairline axes)."""
    mpl.rcParams.update({
        "font.family": "Arial",
        "figure.dpi": DPI,
        "savefig.dpi": DPI,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "axes.linewidth": AX_LW,
        "xtick.major.width": TICK_W,
        "ytick.major.width": TICK_W,
        "xtick.major.size": TICK_LEN,
        "ytick.major.size": TICK_LEN,
    })
    print(f"Figures will be saved to: {FIG_DIR}")


def _ensure_fig_dir():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    return FIG_DIR


def _save_svg(fig, name):
    _ensure_fig_dir()
    path = FIG_DIR / name
    fig.savefig(path, format="svg", bbox_inches="tight")
    plt.close(fig)
    return str(path)


# ---------------------------------------------------------------------------
# Small circular-statistics helpers (plotting-adjacent, kept here).
# ---------------------------------------------------------------------------
def circ_median(angles):
    """Circular median by the minimum cumulative angular difference."""
    angles = np.sort(angles)
    diffs = np.abs(np.angle(np.exp(1j * (angles[:, None] - angles[None, :]))))
    sum_diffs = np.sum(diffs, axis=1)
    return angles[np.argmin(sum_diffs)]


def circular_kde(theta_grid, samples, kappa_bw=80.0):
    """Average von Mises KDE of the samples over ``theta_grid`` (bandwidth in
    concentration units)."""
    samples = np.asarray(samples).ravel()
    return np.mean(vonmises.pdf(theta_grid[None, :], kappa=kappa_bw,
                                loc=samples[:, None]), axis=0)


# ---------------------------------------------------------------------------
# Figure 1: polar grid of raw headings per individual.
# ---------------------------------------------------------------------------
def plot_data_polar(df, radial_max=10):
    """Polar grid of raw headings per individual, light OFF (grey) vs light ON
    (gold), each condition summarised by its circular mean (solid), circular
    median (dashed) and +/- one circular standard deviation (shaded band)."""
    unique_individuals = df["Individual"].unique()
    n_rows, n_cols = 5, 2
    fig, axes = plt.subplots(n_rows, n_cols, subplot_kw={"projection": "polar"},
                             figsize=(10, 18))
    axes = axes.flatten()

    for idx, individual in enumerate(unique_individuals):
        ax = axes[idx]
        individual_data = df[df["Individual"] == individual]

        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)

        on = individual_data[individual_data["Light"] == 1].dropna(subset=["HeadingRadians"])
        off = individual_data[individual_data["Light"] == 0].dropna(subset=["HeadingRadians"])

        _scatter_polar(ax, on, ON, "Light On")
        _scatter_polar(ax, off, OFF, "Light Off")

        if not on.empty:
            _summary_polar(ax, on["HeadingRadians"].values, radial_max, ON)
        if not off.empty:
            _summary_polar(ax, off["HeadingRadians"].values, radial_max, OFF)

        ax.set_title(f"Individual {individual}", pad=30)
        ax.set_ylim(0, radial_max + 1)
        ax.legend().set_visible(False)

    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=ON, markersize=8, label="Light On Data"),
        Line2D([0], [0], color=ON, lw=2, label="Light On Mean"),
        Line2D([0], [0], color=ON, lw=2, linestyle="--", label="Light On Median"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=OFF, markersize=8, label="Light Off Data"),
        Line2D([0], [0], color=OFF, lw=2, label="Light Off Mean"),
        Line2D([0], [0], color=OFF, lw=2, linestyle="--", label="Light Off Median"),
    ]
    fig.legend(handles=handles, loc="center right", bbox_to_anchor=(1.10, 0.5),
               fontsize="medium")

    plt.tight_layout()
    plt.show()
    _save_svg(fig, "data_individual_polar_grid.svg")
    plt.show()


def _scatter_polar(ax, data, color, label):
    if not data.empty:
        ax.scatter(data["HeadingRadians"], data["Trial"], color=color, alpha=0.6, label=label)


def _summary_polar(ax, angles, radial_max, color):
    mean_val = circmean(angles, high=np.pi, low=-np.pi)
    median_val = circ_median(angles)
    std_val = circstd(angles, high=np.pi, low=-np.pi, normalize=True)

    ax.plot([mean_val, mean_val], [0, radial_max], color=color, linewidth=2)
    ax.plot([median_val, median_val], [0, radial_max], color=color,
            linewidth=2, linestyle="--")

    theta_vals = np.linspace(mean_val - std_val, mean_val + std_val, 100)
    r_vals = np.linspace(0, radial_max, 100)
    Theta, R = np.meshgrid(theta_vals, r_vals)
    ax.contourf(Theta, R, np.ones_like(Theta), alpha=0.2, colors=[color])


# ---------------------------------------------------------------------------
# Figure 2 (panel b): polar KDE of the population-level mu posterior.
# ---------------------------------------------------------------------------
def plot_polar_mu_kde(idata_treat, df, shifted_modulo, flood_dir=None,
                      flood_halfwidth=np.deg2rad(45), theta_grid_n=361,
                      kappa_bw=80.0):
    """Polar kernel density of the population-level posterior of mu under each
    lighting condition (grey = off, gold = on), with the floodlight direction
    marked by an arc and each individual's condition mean as a dot on the outer
    ring."""
    _ensure_style()

    post = idata_treat.posterior
    mu_off = shifted_modulo(post["Intercept"].values.reshape(-1))
    mu_on = shifted_modulo((post["Intercept"] + post["Light"]).values.reshape(-1))

    theta = np.linspace(-np.pi, np.pi, theta_grid_n)
    kde_off = circular_kde(theta, mu_off, kappa_bw=kappa_bw)
    kde_on = circular_kde(theta, mu_on, kappa_bw=kappa_bw)

    max_kde = max(kde_off.max(), kde_on.max())
    ring_base, ring_amp = 0.5, 0.5
    r_off = ring_base + ring_amp * (kde_off / max_kde)
    r_on = ring_base + ring_amp * (kde_on / max_kde)

    fig = plt.figure(figsize=(3.0, 3.0))
    ax = fig.add_subplot(111, projection="polar")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.grid(False)
    ax.spines["polar"].set_visible(False)

    ax.fill_between(theta, ring_base, r_off, color=OFF, alpha=0.22, linewidth=0)
    ax.fill_between(theta, ring_base, r_on, color=ON, alpha=0.22, linewidth=0)
    ax.plot(theta, r_off, color=OFF, lw=2)
    ax.plot(theta, r_on, color=ON, lw=2)

    r_pts = ring_base + ring_amp * 1.22
    ax.plot(theta, np.full_like(theta, r_pts), color="0.75", lw=0.9, zorder=1)

    if flood_dir is not None:
        th = np.linspace(flood_dir - flood_halfwidth, flood_dir + flood_halfwidth, 200)
        r_on_max = float(np.max(r_on))
        r_mark = r_on_max + 0.4 * (r_pts - r_on_max)

        ax.plot(th, np.full_like(th, r_mark), color="black", lw=1.1, ls="-",
                alpha=0.75, zorder=6)

        cap = 0.03 * r_pts
        ax.plot([flood_dir, flood_dir], [r_mark - cap, r_mark + cap],
                color="black", lw=1.1, alpha=0.75, zorder=7)
        for t_end in (flood_dir - flood_halfwidth, flood_dir + flood_halfwidth):
            ax.plot([t_end, t_end], [r_mark - cap, r_mark + cap],
                    color="black", lw=1.1, alpha=0.75, zorder=7)

    for ind in df["Individual"].unique():
        sub = df[df["Individual"] == ind]
        off_heads = sub.loc[sub["Light"] == 0, "HeadingRadians"].dropna().to_numpy()
        on_heads = sub.loc[sub["Light"] == 1, "HeadingRadians"].dropna().to_numpy()

        if off_heads.size:
            ang = shifted_modulo(circmean(off_heads, high=np.pi, low=-np.pi))
            ax.scatter([ang], [r_pts], s=100, color=OFF, edgecolor="none", alpha=0.95, zorder=3)
        if on_heads.size:
            ang = shifted_modulo(circmean(on_heads, high=np.pi, low=-np.pi))
            ax.scatter([ang], [r_pts], s=100, color=ON, edgecolor="none", alpha=0.95, zorder=3)

    ax.set_rlim(0, r_pts * 1.1)

    fig.tight_layout(pad=0.2)
    _ensure_fig_dir()
    outpath = FIG_DIR / "treat_polar_mu_kde_clean.svg"
    fig.savefig(outpath, format="svg")
    plt.show()
    plt.close(fig)
    print(f"Saved: {outpath}")


def _ensure_style():
    plt.rcParams.update({
        "font.size": 11,
        "font.family": "Arial",
        "pdf.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "svg.fonttype": "none",
    })


# ---------------------------------------------------------------------------
# Figure 3 (panels c-f): within/between concentration + floodlight attraction.
# ---------------------------------------------------------------------------
def _bessel_ratio_mrl(kappa):
    kappa = np.asarray(kappa)
    return ive(1, kappa) / ive(0, kappa)


def _rho_bar(angles, axis=-1):
    angles = np.asarray(angles)
    C = np.mean(np.cos(angles), axis=axis)
    S = np.mean(np.sin(angles), axis=axis)
    return np.sqrt(C**2 + S**2)


def _mean_hdi(samples, hdi_prob=0.89):
    s = np.asarray(samples).ravel()
    lo, hi = az.hdi(s, hdi_prob=hdi_prob)
    return float(s.mean()), float(lo), float(hi)


def _errorbar_mean_hdi(ax, x, samples, color, hdi_prob=0.89, capsize=3, lw=1, ms=6,
                       label=None):
    m, lo, hi = _mean_hdi(samples, hdi_prob=hdi_prob)
    low = m - lo
    high = hi - m
    ax.errorbar([x], [m], yerr=[[low], [high]], fmt="o", color=color, ecolor=color,
                elinewidth=lw, capsize=capsize, markersize=ms, label=label, zorder=5)


def _vonmises_interval_mass(mu, kappa, a, b):
    twopi = 2 * np.pi
    a = (a + np.pi) % twopi - np.pi
    b = (b + np.pi) % twopi - np.pi
    if a <= b:
        return vonmises.cdf(b, kappa=kappa, loc=mu) - vonmises.cdf(a, kappa=kappa, loc=mu)
    return (vonmises.cdf(np.pi, kappa=kappa, loc=mu) - vonmises.cdf(a, kappa=kappa, loc=mu)
            + vonmises.cdf(b, kappa=kappa, loc=mu) - vonmises.cdf(-np.pi, kappa=kappa, loc=mu))


def plot_concentration(idata_treat, df, shifted_modulo, softplus,
                       flood_dir=np.deg2rad(-135), hdi_prob=0.89):
    """Four panels: within-individual concentration (MRL), between-individual
    concentration (MRL), posterior predictive mass within +/-45 of the
    floodlight, and the on-minus-off change in that mass. Each with an
    off/on error bar (mean, 89% HDI)."""
    _ensure_style()
    _ensure_fig_dir()

    post = idata_treat.posterior

    # c: within-individual concentration (MRL)
    kappa_off = softplus(post["kappa_Intercept"].values.reshape(-1))
    kappa_on = softplus((post["kappa_Intercept"] + post["kappa_Light"]).values.reshape(-1))
    mrl_within_off = _bessel_ratio_mrl(kappa_off)
    mrl_within_on = _bessel_ratio_mrl(kappa_on)

    # d: between-individual concentration (MRL)
    mu_I_fixed = post["Intercept"].values
    mu_L_fixed = post["Light"].values
    mu_I_rand = post["1|Individual"].values
    mu_L_rand = post["Light|Individual"].values

    mu_off = shifted_modulo(mu_I_fixed[..., None] + mu_I_rand)
    mu_on = shifted_modulo(mu_I_fixed[..., None] + mu_L_fixed[..., None] +
                           mu_I_rand + mu_L_rand)
    n_chain, n_draw, n_ind = mu_off.shape
    mrl_between_off = _rho_bar(mu_off.reshape(n_chain * n_draw, n_ind), axis=1)
    mrl_between_on = _rho_bar(mu_on.reshape(n_chain * n_draw, n_ind), axis=1)

    # e/f: predictive mass within +/-45 of the floodlight
    mu_pop_off = shifted_modulo(post["Intercept"].values.reshape(-1))
    mu_pop_on = shifted_modulo((post["Intercept"] + post["Light"]).values.reshape(-1))

    halfwidth = np.deg2rad(45)
    a = flood_dir - halfwidth
    b = flood_dir + halfwidth
    mass_off = _vonmises_interval_mass(mu_pop_off, kappa_off, a, b)
    mass_on = _vonmises_interval_mass(mu_pop_on, kappa_on, a, b)
    mass_diff = mass_on - mass_off

    fig, axes = plt.subplots(1, 4, figsize=(8.0, 3.0), sharey=False)

    ax = axes[0]
    ax.set_title("Within-individual\nconcentration")
    _errorbar_mean_hdi(ax, 0, mrl_within_off, OFF, hdi_prob=hdi_prob, capsize=3)
    _errorbar_mean_hdi(ax, 1, mrl_within_on, ON, hdi_prob=hdi_prob, capsize=3)
    ax.set_xticks([0, 1], ["Off", "On"])
    ax.set_ylim(0, 1.0)

    ax = axes[1]
    ax.set_title("Between-individual\nconcentration")
    _errorbar_mean_hdi(ax, 0, mrl_between_off, OFF, hdi_prob=hdi_prob, capsize=3)
    _errorbar_mean_hdi(ax, 1, mrl_between_on, ON, hdi_prob=hdi_prob, capsize=3)
    ax.set_xticks([0, 1], ["Off", "On"])
    ax.set_ylim(0, 1.0)

    ax = axes[2]
    ax.set_title("Floodlight\nattraction")
    _errorbar_mean_hdi(ax, 0, mass_off, OFF, hdi_prob=hdi_prob, capsize=3)
    _errorbar_mean_hdi(ax, 1, mass_on, ON, hdi_prob=hdi_prob, capsize=3)
    ax.set_xticks([0, 1], ["Off", "On"])
    ax.set_ylabel("Predictive mass (+-45 deg)")
    uniform_mass = (2 * halfwidth) / (2 * np.pi)
    ax.axhline(uniform_mass, color="k", lw=1, ls=(0, (5, 5)), alpha=1.0, zorder=1)
    ax.set_ylim(0.0, 1.0)
    ax.text(0.98, uniform_mass - 0.03, "Uniform",
            transform=ax.get_yaxis_transform(), ha="right", va="top",
            fontsize=10, color="k", alpha=1.0)

    ax = axes[3]
    ax.set_title("Floodlight\ntreatment effect")
    _errorbar_mean_hdi(ax, 0, mass_diff, ON, hdi_prob=hdi_prob, capsize=3)
    ax.set_xticks([0], ["On - Off"])
    ax.set_ylabel("Difference in\npredictive mass (+-45 deg)")
    ax.set_ylim(0.0, 1.0)

    yt_mrl = np.linspace(0, 1.0, 6)
    for i in [0, 1]:
        axes[i].set_yticks(yt_mrl)
        axes[i].tick_params(axis="y", which="both", left=True, labelleft=True)
        axes[i].set_ylabel("Mean resultant length")

    for ax in axes:
        ax.margins(x=0.25)
        if ax is not axes[0]:
            ax.tick_params(axis="y", which="both", left=True, labelleft=True)

    fig.tight_layout()
    outpath = FIG_DIR / "treat_concentration_plus_floodlight.svg"
    fig.savefig(outpath, format="svg")
    plt.show()
    plt.close(fig)
    print(f"Saved: {outpath}")


# ---------------------------------------------------------------------------
# Figure 4: supplemental PPC — bootstrap (data) vs posterior (model).
# ---------------------------------------------------------------------------
def _circ_mod(x):
    return (x + np.pi) % (2 * np.pi) - np.pi


def _natural_sort_key(s):
    return [int(t) if t.isdigit() else t.lower() for t in re.split("([0-9]+)", str(s))]


def _to_360(rad):
    return (np.degrees(rad) + 360) % 360


def _circ_mean(theta):
    theta = np.asarray(theta, float)
    z = np.exp(1j * theta).mean()
    return float(_circ_mod(np.angle(z)))


def _linear_hdi(x, prob=0.95):
    x = np.sort(np.asarray(x, float))
    n = x.size
    if n == 0:
        return np.nan, np.nan
    k = int(np.ceil(prob * n))
    spans = x[k - 1:] - x[:n - k + 1]
    j = int(np.argmin(spans))
    return float(x[j]), float(x[j + k - 1])


def _circ_hdi_centered(samples, prob=0.95):
    a = np.asarray(samples, float)
    if a.size == 0:
        return np.nan, np.nan, np.nan
    m = _circ_mean(a)
    d = _circ_mod(a - m)
    lo_d, hi_d = _linear_hdi(d, prob)
    lo = float(_circ_mod(m + lo_d))
    hi = float(_circ_mod(m + hi_d))
    return m, lo, hi


def _circ_xerr(m, lo, hi):
    w_lo = abs(float(_circ_mod(m - lo)))
    w_hi = abs(float(_circ_mod(hi - m)))
    return np.degrees(w_lo), np.degrees(w_hi)


def _A_of_kappa(kappa):
    kappa = np.asarray(kappa, float)
    return i1(kappa) / i0(kappa)


def _resultant_length_angles(a, axis=-1):
    a = np.asarray(a, float)
    return np.sqrt(np.mean(np.cos(a), axis=axis)**2 + np.mean(np.sin(a), axis=axis)**2)


def _stack_1d(post, name):
    if name not in post.data_vars:
        raise KeyError(f"Missing posterior variable '{name}'.")
    return post[name].stack(sample=("chain", "draw")).values.astype(float)


def _restore_ticks(ax):
    ax.tick_params(axis="both", which="both", direction="out", bottom=True, top=False,
                   left=True, right=False, length=TICK_LEN, width=TICK_W)
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")


def plot_ppc(bootG, boot_mu_off, boot_mu_on, boot_R_off, boot_R_on, indivs,
             mu_pop_off_draw, mu_pop_on_draw,
             mu_off_draw, mu_on_draw, R_off_draw, R_on_draw,
             R_group_off_draw, R_group_on_draw):
    """Supplemental PPC figure comparing a nonparametric bootstrap of the data
    against the model posterior, for group mean direction, group alignment, and
    per-individual concentration, off (top-left/right) vs on (bottom)."""
    _ensure_style()
    _ensure_fig_dir()
    nT = len(indivs)

    # ---- group-level summaries (95% HDIs) ----
    muG_off_m_b, muG_off_lo_b, muG_off_hi_b = _circ_hdi_centered(bootG["mu_off"][~np.isnan(bootG["mu_off"])])
    muG_on_m_b, muG_on_lo_b, muG_on_hi_b = _circ_hdi_centered(bootG["mu_on"][~np.isnan(bootG["mu_on"])])
    muG_off_m_p, muG_off_lo_p, muG_off_hi_p = _circ_hdi_centered(mu_pop_off_draw)
    muG_on_m_p, muG_on_lo_p, muG_on_hi_p = _circ_hdi_centered(mu_pop_on_draw)

    R_off_b = bootG["R_off"][~np.isnan(bootG["R_off"])]
    R_on_b = bootG["R_on"][~np.isnan(bootG["R_on"])]
    R_off_m_b, (R_off_lo_b, R_off_hi_b) = np.nanmean(R_off_b), _linear_hdi(R_off_b)
    R_on_m_b, (R_on_lo_b, R_on_hi_b) = np.nanmean(R_on_b), _linear_hdi(R_on_b)
    R_off_m_p, (R_off_lo_p, R_off_hi_p) = np.mean(R_group_off_draw), _linear_hdi(R_group_off_draw)
    R_on_m_p, (R_on_lo_p, R_on_hi_p) = np.mean(R_group_on_draw), _linear_hdi(R_group_on_draw)

    # ---- per-individual rows ----
    indiv_mu_rows, indiv_R_rows = [], []
    for j, t in enumerate(indivs):
        for cond, b_mu, b_r, p_mu, p_r in [("OFF", boot_mu_off, boot_R_off,
                                            mu_off_draw, R_off_draw),
                                           ("ON", boot_mu_on, boot_R_on,
                                            mu_on_draw, R_on_draw)]:
            if b_mu[t].size:
                m, lo, hi = _circ_hdi_centered(b_mu[t])
                indiv_mu_rows.append((cond, t, m, lo, hi, "bootstrap"))
                indiv_R_rows.append((cond, t, float(np.mean(b_r[t])),
                                     *_linear_hdi(b_r[t]), "bootstrap"))
            m, lo, hi = _circ_hdi_centered(p_mu[:, j])
            indiv_mu_rows.append((cond, t, m, lo, hi, "posterior"))
            indiv_R_rows.append((cond, t, float(np.mean(p_r[:, j])),
                                 *_linear_hdi(p_r[:, j]), "posterior"))

    mu_tbl = pd.DataFrame(indiv_mu_rows, columns=["Condition", "Individual", "mean", "low", "high", "Source"])
    R_tbl = pd.DataFrame(indiv_R_rows, columns=["Condition", "Individual", "mean", "low", "high", "Source"])

    # ---- figure layout ----
    OFF_COLOR, ON_COLOR = OFF, ON
    OFF_MARK, ON_MARK = "s", "o"
    fig = plt.figure(figsize=(8.8, max(3.4, 0.28 * nT + 2.2)), dpi=DPI)
    gs = fig.add_gridspec(2, 4, height_ratios=[0.55, max(1.2, 0.55 + 0.05 * nT)],
                          wspace=0.45, hspace=0.38)

    ax_g_mu, ax_g_R = fig.add_subplot(gs[0, 0:2]), fig.add_subplot(gs[0, 2:4])
    ax_mu_off, ax_mu_on = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
    ax_R_off, ax_R_on = fig.add_subplot(gs[1, 2]), fig.add_subplot(gs[1, 3])

    PANEL_FS = 11

    def _panel_label(ax, lab):
        ax.text(0.05, 1.0, lab, transform=ax.transAxes, ha="left", va="top",
                fontsize=PANEL_FS, fontweight="bold", fontfamily="Arial",
                clip_on=False, zorder=10)

    _panel_label(ax_g_mu, "a.")
    _panel_label(ax_g_R, "b.")
    _panel_label(ax_mu_off, "c.")
    _panel_label(ax_mu_on, "d.")
    _panel_label(ax_R_off, "e.")
    _panel_label(ax_R_on, "f.")

    # (A) group mean direction
    y_base = np.array([0.0, 1.0])
    y_boot, y_post = y_base - 0.12, y_base + 0.12

    def _plot_dir_point(ax, m, lo, hi, y, color, marker, filled):
        m_deg = _to_360(m)
        wlo, whi = _circ_xerr(m, lo, hi)
        ax.errorbar(m_deg, y, xerr=np.array([[wlo], [whi]]), fmt=marker, color=color,
                    ecolor=color, capsize=2, ms=3.2, mfc=(color if filled else "white"),
                    mec=color, mew=0.9, lw=(1.2 if filled else 1.0), zorder=3)

    _plot_dir_point(ax_g_mu, muG_off_m_b, muG_off_lo_b, muG_off_hi_b, y_boot[0], OFF_COLOR, OFF_MARK, filled=False)
    _plot_dir_point(ax_g_mu, muG_off_m_p, muG_off_lo_p, muG_off_hi_p, y_post[0], OFF_COLOR, OFF_MARK, filled=True)
    _plot_dir_point(ax_g_mu, muG_on_m_b, muG_on_lo_b, muG_on_hi_b, y_boot[1], ON_COLOR, ON_MARK, filled=False)
    _plot_dir_point(ax_g_mu, muG_on_m_p, muG_on_lo_p, muG_on_hi_p, y_post[1], ON_COLOR, ON_MARK, filled=True)

    ax_g_mu.set_yticks(y_base); ax_g_mu.set_yticklabels(["OFF", "ON"])
    ax_g_mu.set_xlim(0, 360); ax_g_mu.set_xticks([0, 90, 180, 270, 360])
    ax_g_mu.set_xlabel("Group mean direction (deg)"); ax_g_mu.set_title("Group mean direction")
    _restore_ticks(ax_g_mu)

    # (B) group alignment
    rows = [("OFF", OFF_COLOR, OFF_MARK, R_off_m_b, R_off_lo_b, R_off_hi_b, R_off_m_p, R_off_lo_p, R_off_hi_p),
            ("ON", ON_COLOR, ON_MARK, R_on_m_b, R_on_lo_b, R_on_hi_b, R_on_m_p, R_on_lo_p, R_on_hi_p)]
    for i, (_, col, mark, mb, lob, hib, mp, lop, hip) in enumerate(rows):
        if np.isfinite(mb):
            ax_g_R.errorbar(mb, y_boot[i], xerr=np.array([[mb - lob], [hib - mb]]), fmt=mark,
                            color=col, ecolor=col, capsize=2, ms=3.2, mfc="white", mec=col,
                            mew=0.9, lw=1.0, zorder=3)
        if np.isfinite(mp):
            ax_g_R.errorbar(mp, y_post[i], xerr=np.array([[mp - lop], [hip - mp]]), fmt=mark,
                            color=col, ecolor=col, capsize=2, ms=3.2, mfc=col, mec=col,
                            mew=0.9, lw=1.2, zorder=3)
    ax_g_R.set_yticks(y_base); ax_g_R.set_yticklabels(["OFF", "ON"])
    ax_g_R.set_xlim(0, 1.0); ax_g_R.set_xlabel("Mean resultant length"); ax_g_R.set_title("Group alignment")
    _restore_ticks(ax_g_R)

    # individuals (bottom)
    def _plot_ind_mu(ax, cond, col, mark):
        sub = mu_tbl[mu_tbl["Condition"] == cond]
        for j, t in enumerate(indivs):
            for src, dy, fill in [("bootstrap", -0.14, "white"), ("posterior", 0.14, col)]:
                r = sub[(sub["Individual"] == t) & (sub["Source"] == src)]
                if not r.empty:
                    m = r["mean"].iloc[0]
                    wlo, whi = _circ_xerr(m, r["low"].iloc[0], r["high"].iloc[0])
                    ax.errorbar(_to_360(m), j + dy, xerr=np.array([[wlo], [whi]]), fmt=mark,
                                color=col, ecolor=col, capsize=2, ms=2.7, mfc=fill, mec=col,
                                mew=0.8, lw=0.9, zorder=3)

    def _plot_ind_R(ax, cond, col, mark):
        sub = R_tbl[R_tbl["Condition"] == cond]
        for j, t in enumerate(indivs):
            for src, dy, fill in [("bootstrap", -0.14, "white"), ("posterior", 0.14, col)]:
                r = sub[(sub["Individual"] == t) & (sub["Source"] == src)]
                if not r.empty:
                    m = r["mean"].iloc[0]
                    ax.errorbar(m, j + dy, xerr=np.array([[m - r["low"].iloc[0]],
                                                         [r["high"].iloc[0] - m]]),
                                fmt=mark, color=col, ecolor=col, capsize=2, ms=2.7,
                                mfc=fill, mec=col, mew=0.8, lw=0.9, zorder=3)

    _plot_ind_mu(ax_mu_off, "OFF", OFF_COLOR, OFF_MARK)
    _plot_ind_mu(ax_mu_on, "ON", ON_COLOR, ON_MARK)
    _plot_ind_R(ax_R_off, "OFF", OFF_COLOR, OFF_MARK)
    _plot_ind_R(ax_R_on, "ON", ON_COLOR, ON_MARK)

    for ax in (ax_mu_off, ax_mu_on, ax_R_off, ax_R_on):
        ax.set_yticks(np.arange(nT))
        ax.invert_yaxis()
        _restore_ticks(ax)
        for y in np.arange(-0.5, nT, 1.0):
            ax.axhline(y, color="0.90", lw=0.5, zorder=0)

    ax_mu_off.set_yticklabels(indivs)
    for ax in (ax_mu_on, ax_R_off, ax_R_on):
        ax.set_yticklabels([])
    ax_mu_off.set_xlim(0, 360); ax_mu_on.set_xlim(0, 360)
    ax_mu_off.set_xticks([0, 180, 360]); ax_mu_on.set_xticks([0, 180, 360])

    pL, pR = ax_mu_off.get_position(), ax_mu_on.get_position()
    fig.text(0.5 * (pL.x0 + pR.x1), pL.y0 - 0.045, "Individual mean direction (deg)",
             ha="center", va="top")
    pL, pR = ax_R_off.get_position(), ax_R_on.get_position()
    fig.text(0.5 * (pL.x0 + pR.x1), pL.y0 - 0.045, "Individual concentration (R)",
             ha="center", va="top")

    handles = [Line2D([0], [0], marker="o", linestyle="none", mfc="white", mec="k", label="Bootstrap (data)"),
               Line2D([0], [0], marker="o", linestyle="none", mfc="k", mec="k", label="Posterior (model)"),
               Line2D([0], [0], color=OFF_COLOR, marker=OFF_MARK, label="Light OFF"),
               Line2D([0], [0], color=ON_COLOR, marker=ON_MARK, label="Light ON")]
    fig.legend(handles=handles, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.01))

    out_path = FIG_DIR / "beetle_posterior_predictive_check.pdf"
    fig.savefig(out_path, bbox_inches="tight", pad_inches=0.01)
    plt.show()
    print(f"Saved: {out_path}")
