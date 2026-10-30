"""Figure helpers for the dung-beetle security-light notebook.

The notebook keeps the statistical work (the model-ready table from
``beetles_preprocess``, the three bambi von Mises fits, the LOO model
comparison, and the post-hoc concentration analysis). This module holds
the plotting "how": the per-individual data plots, the trace plots, the
population posterior PDFs, and the five-panel summary figure, plus the
post-hoc kappa-mu extractions those figures need.

Each function renders its figure at the display tier (``figure.dpi=100``,
so the inline notebook embed stays small) and writes a publication vector
PDF (type-42 outlines, 300 dpi raster fallback) under ``figures/beetles/``
at the repository root, resolved from this file's location so the output
does not depend on the working directory.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import arviz as az
from matplotlib.lines import Line2D
from scipy.stats import circmean, circstd, vonmises

from unwrap import plots as _plots
from unwrap.link.numpy import shifted_modulo, softplus

#: Publication figures for this study live in <repo>/figures/beetles/.
FIG_DIR = Path(__file__).resolve().parent.parent / "figures" / "beetles"


def set_style():
    """Apply the shared two-resolution house style (see unwrap.plots)."""
    _plots.set_style(FIG_DIR)


def _circ_median(angles):
    """Circular median: the input angle minimising the summed circular distances."""
    angles = np.sort(angles)
    diffs = np.abs(np.angle(np.exp(1j * (angles[:, None] - angles[None, :]))))
    sum_diffs = np.sum(diffs, axis=1)
    return angles[np.argmin(sum_diffs)]


def plot_data_scatter(df):
    """Polar scatter of each individual's headings, light on (blue) vs
    off (red), with per-condition circular mean (solid), median (dashed)
    and +/- circular standard deviation filled. One shared legend for all
    panels (saved as ``data_scatter.pdf``)."""
    set_style()

    unique_individuals = df["Individual"].unique()
    n_rows, n_cols = 5, 2
    fig, axes = plt.subplots(n_rows, n_cols,
                             subplot_kw={"projection": "polar"}, figsize=(10, 18))
    axes = axes.flatten()

    radial_max = 10  # maximum radial value (trials per individual)

    for idx, individual in enumerate(unique_individuals):
        ax = axes[idx]
        individual_data = df[df["Individual"] == individual]

        ax.set_theta_zero_location("N")  # place 0 at the top
        ax.set_theta_direction(-1)  # clockwise

        data_on = individual_data[individual_data["Light"] == 1].dropna(subset=["HeadingRadians"])
        data_off = individual_data[individual_data["Light"] == 0].dropna(subset=["HeadingRadians"])

        for data, color, label in ((data_on, "blue", "Light On"),
                                   (data_off, "red", "Light Off")):
            if not data.empty:
                ax.scatter(data["HeadingRadians"], data["Trial"],
                           color=color, alpha=0.6, label=label)
                angles = data["HeadingRadians"].values
                mean_val = circmean(angles, high=np.pi, low=-np.pi)
                median_val = _circ_median(angles)
                std_val = circstd(angles, high=np.pi, low=-np.pi, normalize=True)
                ax.plot([mean_val, mean_val], [0, radial_max], color=color, linewidth=2)
                ax.plot([median_val, median_val], [0, radial_max],
                        color=color, linewidth=2, linestyle="--")
                theta_vals = np.linspace(mean_val - std_val, mean_val + std_val, 100)
                r_vals = np.linspace(0, radial_max, 100)
                Theta, R = np.meshgrid(theta_vals, r_vals)
                ax.contourf(Theta, R, np.ones_like(Theta), alpha=0.2, colors=[color])

        ax.set_title(f"Individual {individual}", pad=30)
        ax.set_ylim(0, radial_max + 1)
        # Remove per-panel legends in favour of one shared legend.
        ax.legend().set_visible(False)

    legend_handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="blue",
               markersize=8, label="Light On Data"),
        Line2D([0], [0], color="blue", lw=2, label="Light On Mean"),
        Line2D([0], [0], color="blue", lw=2, linestyle="--", label="Light On Median"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="red",
               markersize=8, label="Light Off Data"),
        Line2D([0], [0], color="red", lw=2, label="Light Off Mean"),
        Line2D([0], [0], color="red", lw=2, linestyle="--", label="Light Off Median"),
    ]
    fig.legend(handles=legend_handles, loc="center right",
               bbox_to_anchor=(1.10, 0.5), fontsize="medium")

    fig.tight_layout()
    _plots.save_figure(fig, FIG_DIR, "data_scatter")
    return fig


def plot_trace(idata, name, var_names=None, circ_var_names=None):
    """ArviZ traces for one of the fits, under the house style
    (saved as ``<name>.png``). Pass either the subset of ``var_names``
    or the circular ``circ_var_names`` to plot, as needed.

    ``compact=True`` keeps each (possibly multi-dimensional) parameter in
    one row so these 20+-panel figures stay a reasonable size inline; the
    saved file is a 300 dpi raster, since the density panels contain
    thousands of vector paths whose PDF would be several MB.
    """
    set_style()
    if var_names is not None:
        az.plot_trace(idata, var_names=var_names, compact=True)
    else:
        az.plot_trace(idata, circ_var_names=circ_var_names, compact=True)
    fig = plt.gcf()
    # Defensive: keep the inline embed at the display tier no matter how
    # the plotting library configured the figure it opened.
    fig.set_dpi(_plots.DEFAULT_DISPLAY_DPI)
    fig.tight_layout()
    _plots.save_figure(fig, FIG_DIR, name, pdf=False, png=True)
    return fig


def plot_posterior_pdfs(idata_treat, df):
    """Median population-level posterior von Mises PDFs (with 80% credible
    bands) for light off/on from the treatment model, overlaid with each
    individual's heading circular means (saved as ``posterior_pdfs.pdf``)."""
    set_style()

    # Posterior samples for the population-level parameters (bambi naming).
    mu_intercept = idata_treat.posterior["Intercept"].values  # fixed intercept for mu
    mu_slope = idata_treat.posterior["Light"].values  # fixed slope for mu
    kappa_intercept = idata_treat.posterior["kappa_Intercept"].values  # kappa linear predictor
    kappa_slope = idata_treat.posterior["kappa_Light"].values

    theta = np.linspace(-np.pi, np.pi, 200)
    fig, ax = plt.subplots(figsize=(14, 7))

    # Population-level parameters per condition.
    mu_pop_off_samples_all = shifted_modulo(mu_intercept.flatten())
    mu_pop_on_samples_all = shifted_modulo((mu_intercept + mu_slope).flatten())
    kappa_pop_off_samples_all = softplus(kappa_intercept.flatten())
    kappa_pop_on_samples_all = softplus((kappa_intercept + kappa_slope).flatten())

    # Von Mises PDF of each posterior sample at each theta.
    pdf_off = vonmises.pdf(theta[None, :], kappa=kappa_pop_off_samples_all[:, None],
                           loc=mu_pop_off_samples_all[:, None])
    pdf_on = vonmises.pdf(theta[None, :], kappa=kappa_pop_on_samples_all[:, None],
                          loc=mu_pop_on_samples_all[:, None])

    median_pdf_off = np.median(pdf_off, axis=0)
    lower_pdf_off = np.percentile(pdf_off, 10, axis=0)
    upper_pdf_off = np.percentile(pdf_off, 90, axis=0)
    median_pdf_on = np.median(pdf_on, axis=0)
    lower_pdf_on = np.percentile(pdf_on, 10, axis=0)
    upper_pdf_on = np.percentile(pdf_on, 90, axis=0)

    ax.plot(theta, median_pdf_off, color="red", lw=2, label="Median Posterior PDF Off")
    ax.fill_between(theta, lower_pdf_off, upper_pdf_off, color="red", alpha=0.3,
                    label="80% Credible Interval Off")
    ax.plot(theta, median_pdf_on, color="blue", lw=2, label="Median Posterior PDF On")
    ax.fill_between(theta, lower_pdf_on, upper_pdf_on, color="blue", alpha=0.3,
                    label="80% Credible Interval On")

    # Individual-level circular means.
    for i, individual in enumerate(df["Individual"].unique()):
        individual_data = df[df["Individual"] == individual]
        headings_off = individual_data.loc[individual_data["Light"] == 0, "HeadingRadians"].dropna()
        headings_on = individual_data.loc[individual_data["Light"] == 1, "HeadingRadians"].dropna()
        individual_circmean_off = shifted_modulo(circmean(headings_off))
        individual_circmean_on = shifted_modulo(circmean(headings_on))
        ax.scatter(individual_circmean_off, 0, color="red", alpha=0.6,
                   label="Individual Means Off" if i == 0 else "")
        ax.scatter(individual_circmean_on, 0, color="blue", alpha=0.6,
                   label="Individual Means On" if i == 0 else "")

    ax.set_title("Von Mises PDF Credible Intervals with Individual Circular Means")
    ax.set_xlabel("Heading (radians)")
    ax.set_ylabel("Density")
    ax.legend()
    _plots.save_figure(fig, FIG_DIR, "posterior_pdfs")
    return fig


def kappamu_individual(idata_treat, df):
    """Post-hoc kappa-mu, per individual (the "mean over draws" recipe).

    For each individual the mu posterior draws are wrapped fixed + random
    effects; the circular mean over draws gives each individual's mu
    estimate, and a von Mises fit to the set of individual estimates
    (separately per condition) gives the concentration of individual
    heading biases, kappa-mu. Prints the two fitted concentrations.

    Returns a dict of the effect arrays and individual estimates;
    :func:`kappamu_per_draw` refines this per draw.
    """
    set_style()

    mu_intercept_fixed = idata_treat.posterior["Intercept"].values
    mu_light_fixed = idata_treat.posterior["Light"].values
    mu_intercept_random = idata_treat.posterior["1|Individual"].values
    mu_light_random = idata_treat.posterior["Light|Individual"].values

    # Sum the fixed and random effects for each draw and individual; take
    # the circular mean over draws for the point estimate.
    mu_intercept_draws = shifted_modulo(mu_intercept_fixed[..., np.newaxis] + mu_intercept_random)
    mu_intercept_individual_estimates = shifted_modulo(circmean(mu_intercept_draws, axis=(0, 1)))

    # For the light-on condition.
    mu_light_draws = shifted_modulo(
        mu_intercept_fixed[..., np.newaxis] + mu_light_fixed[..., np.newaxis]
        + mu_intercept_random + mu_light_random
    )
    mu_light_individual_estimates = shifted_modulo(circmean(mu_light_draws, axis=(0, 1)))

    # Fit a von Mises to the set of individual mu estimates (fscale=1);
    # the alternative is inverse_bessel_ratio on the mean resultant length.
    kappa_mu_intercept_fit, _, _ = vonmises.fit(mu_intercept_individual_estimates, fscale=1)
    print("\n--- Kappa_mu from the Distribution of Individual Mu (Intercept) ---")
    print(f"The Kappa_mu (concentration) of the individual mu estimates: {kappa_mu_intercept_fit:.4f}")

    kappa_mu_light_fit, _, _ = vonmises.fit(mu_light_individual_estimates, fscale=1)
    print("\n--- Kappa_mu from the Distribution of Individual Mu (Light) ---")
    print(f"The Kappa_mu (concentration) of the individual mu estimates: {kappa_mu_light_fit:.4f}")

    return {
        "mu_intercept_fixed": mu_intercept_fixed,
        "mu_light_fixed": mu_light_fixed,
        "mu_intercept_random": mu_intercept_random,
        "mu_light_random": mu_light_random,
        "individual_ids": df["Individual"].unique(),
        "mu_intercept_individual_estimates": mu_intercept_individual_estimates,
        "mu_light_individual_estimates": mu_light_individual_estimates,
        "kappa_mu_intercept_fit": kappa_mu_intercept_fit,
        "kappa_mu_light_fit": kappa_mu_light_fit,
    }


def kappamu_per_draw(post_hoc):
    """Post-hoc kappa-mu, per draw (refines :func:`kappamu_individual`).

    For every (chain, draw) sample, the set of that draw's wrapped
    individual mus is fit with a von Mises (per condition); the resulting
    arrays are the posterior distributions of kappa-mu and its condition
    difference. A failed MLE is treated as very high concentration (100).
    Prints the ArviZ summary of the three posteriors.

    ``post_hoc`` is the dict from :func:`kappamu_individual`; returns it
    extended with the per-draw posterior arrays
    (kappamu_intercept/light/delta_posterior).
    """
    set_style()

    mu_intercept_fixed = post_hoc["mu_intercept_fixed"]
    mu_light_fixed = post_hoc["mu_light_fixed"]
    mu_intercept_random = post_hoc["mu_intercept_random"]
    mu_light_random = post_hoc["mu_light_random"]

    n_chain, n_draws, n_individuals = mu_intercept_random.shape
    total_draws = n_chain * n_draws

    kappamu_intercept_posterior = np.empty(total_draws)
    kappamu_light_posterior = np.empty(total_draws)
    kappamu_delta_posterior = np.empty(total_draws)

    # (total_draws, n_individuals) and (total_draws, 1) views.
    mu_i_random_flat = mu_intercept_random.reshape(total_draws, n_individuals)
    mu_L_random_flat = mu_light_random.reshape(total_draws, n_individuals)
    mu_I_fixed_flat = mu_intercept_fixed.reshape(total_draws, 1)
    mu_L_fixed_flat = mu_light_fixed.reshape(total_draws, 1)

    print(f"Calculating second-level kappamu for {total_draws} draws...")

    for i in range(total_draws):
        # Individual-level mu for this draw, per condition.
        mu_intercept_individual_draw = (
            mu_I_fixed_flat[i] + mu_i_random_flat[i]
        )
        mu_light_individual_draw = (
            mu_I_fixed_flat[i] + mu_L_fixed_flat[i]
            + mu_i_random_flat[i] + mu_L_random_flat[i]
        )
        mu_intercept_wrapped = shifted_modulo(mu_intercept_individual_draw)
        mu_light_wrapped = shifted_modulo(mu_light_individual_draw)

        try:
            kappamu_intercept_fit, _, _ = vonmises.fit(mu_intercept_wrapped, fscale=1)
            kappamu_intercept_posterior[i] = kappamu_intercept_fit
        except RuntimeError:
            # The MLE failed (e.g. all mus identical): use a high value.
            kappamu_light_posterior[i] = 100  # high enough
        try:
            kappamu_light_fit, _, _ = vonmises.fit(mu_light_wrapped, fscale=1)
            kappamu_light_posterior[i] = kappamu_light_fit
        except RuntimeError:
            # The MLE failed (e.g. all mus identical): use a high value.
            kappamu_light_posterior[i] = 100  # high enough

        kappamu_delta_posterior[i] = (
            kappamu_light_posterior[i] - kappamu_intercept_posterior[i]
        )  # change in kappa_mu

    # Filter out any NaNs if the optimizer failed on certain draws.
    kappamu_intercept_posterior = kappamu_intercept_posterior[~np.isnan(kappamu_intercept_posterior)]
    kappamu_light_posterior = kappamu_light_posterior[~np.isnan(kappamu_light_posterior)]

    kappamu_posterior_idata = {
        "kappa_mu_intercept": kappamu_intercept_posterior,
        "kappa_mu_light": kappamu_light_posterior,
        "kappa_mu_delta": kappamu_delta_posterior,
    }

    print("\n--- Summary of Kappa_mu (Light Condition) ---")
    print("This represents the posterior distribution of the concentration of individual mu's.")
    print(az.summary(kappamu_posterior_idata, hdi_prob=0.95))

    post_hoc.update({
        "kappamu_intercept_posterior": kappamu_intercept_posterior,
        "kappamu_light_posterior": kappamu_light_posterior,
        "kappamu_delta_posterior": kappamu_delta_posterior,
    })
    return post_hoc


def _individual_heading_means(df, ax, y):
    """Scatter each individual's conditional heading circular means at height ``y``."""
    for i, individual in enumerate(df["Individual"].unique()):
        individual_data = df[df["Individual"] == individual]
        headings_off = individual_data.loc[individual_data["Light"] == 0, "HeadingRadians"].dropna()
        headings_on = individual_data.loc[individual_data["Light"] == 1, "HeadingRadians"].dropna()
        individual_circmean_off = shifted_modulo(circmean(headings_off))
        individual_circmean_on = shifted_modulo(circmean(headings_on))
        ax.scatter(individual_circmean_off, y, color="red", alpha=0.6,
                   label="Individual Means Off" if i == 0 else "")
        ax.scatter(individual_circmean_on, y, color="blue", alpha=0.6,
                   label="Individual Means On" if i == 0 else "")


def plot_summary(df, idata_treat, post_hoc):
    """The five-panel summary of the treatment model: polar KDEs of the
    population mu (plain and driven by kappa-mu), split violins of the
    posterior kappa and kappa-mu, and of the kappa-mu difference
    (post-hoc vs model) (saved as ``post_summary.pdf``).

    ``post_hoc`` is the dict returned by :func:`kappamu_per_draw`.
    """
    set_style()

    # Posterior samples for the population-level parameters (bambi naming).
    mu_intercept = idata_treat.posterior["Intercept"].values  # fixed intercept for mu
    mu_slope = idata_treat.posterior["Light"].values  # fixed slope for mu
    kappa_intercept = idata_treat.posterior["kappa_Intercept"].values  # kappa linear predictor
    kappa_slope = idata_treat.posterior["kappa_Light"].values
    kappamu_intercept = idata_treat.posterior["1|Individual_kappa"].values  # kappa-mu linear predictor
    kappamu_slope = idata_treat.posterior["Light|Individual_kappa"].values

    theta = np.linspace(-np.pi, np.pi, 200)

    mu_pop_off_samples_all = shifted_modulo(mu_intercept.flatten())
    mu_pop_on_samples_all = shifted_modulo((mu_intercept + mu_slope).flatten())
    kappa_pop_off_samples_all = softplus(kappa_intercept.flatten())
    kappa_pop_on_samples_all = softplus((kappa_intercept + kappa_slope).flatten())
    kappamu_pop_off_samples_all = softplus(kappamu_intercept.flatten())
    kappamu_pop_on_samples_all = softplus(kappamu_slope.flatten())

    kappamu_intercept_posterior = post_hoc["kappamu_intercept_posterior"]
    kappamu_light_posterior = post_hoc["kappamu_light_posterior"]
    kappamu_delta_posterior = post_hoc["kappamu_delta_posterior"]

    # ---------------------------
    # Process mu (circular) samples
    # ---------------------------
    # Circular KDE for each condition: evaluate a narrow von Mises kernel
    # at each theta for every sample, then average over samples.
    mu_off_centered = shifted_modulo(mu_pop_off_samples_all)
    mu_on_centered = shifted_modulo(mu_pop_on_samples_all)
    kappa_bandwidth = 500  # kernel concentration

    kde_mu_off = np.mean(
        vonmises.pdf(theta[None, :], kappa=kappa_bandwidth,
                     loc=mu_off_centered[:, None]),
        axis=0
    )
    kde_mu_on = np.mean(
        vonmises.pdf(theta[None, :], kappa=kappa_bandwidth,
                     loc=mu_on_centered[:, None]),
        axis=0
    )
    kde_kappamu_off = np.mean(
        vonmises.pdf(theta[None, :], kappa=np.median(kappamu_pop_off_samples_all),
                     loc=mu_off_centered[:, None]),
        axis=0
    )
    kde_kappamu_on = np.mean(
        vonmises.pdf(theta[None, :], kappa=np.median(kappamu_pop_on_samples_all),
                     loc=mu_on_centered[:, None]),
        axis=0
    )

    # ---------------------------
    # Process kappa / kappa-mu (linear) samples
    # ---------------------------
    df_kappa = pd.DataFrame({
        "Parameter": ["Kappa"] * (len(kappa_pop_off_samples_all) + len(kappa_pop_on_samples_all)),
        "Light": ["Light Off"] * len(kappa_pop_off_samples_all) + ["Light On"] * len(kappa_pop_on_samples_all),
        "Value": np.concatenate([kappa_pop_off_samples_all, kappa_pop_on_samples_all]),
    })
    df_kappamu = pd.DataFrame({
        "Parameter": ["Kappa_mu"] * (len(kappamu_intercept_posterior) + len(kappamu_light_posterior)),
        "Light": ["Light Off"] * len(kappamu_intercept_posterior) + ["Light On"] * len(kappamu_light_posterior),
        "Value": np.concatenate([kappamu_intercept_posterior, kappamu_light_posterior]),
    })

    # ---------------------------
    # The five panels
    # ---------------------------
    fig = plt.figure(figsize=(14, 6))

    # (A) Polar KDE of the population mu.
    ax1 = fig.add_subplot(1, 5, 1, projection="polar")
    ax1.set_theta_zero_location("N")
    ax1.set_theta_direction(-1)
    # Offset the density by 1 so there is white space in the center.
    ax1.plot(theta, kde_mu_off + 1, color="red", lw=2, label="Light Off")
    ax1.plot(theta, kde_mu_on + 1, color="blue", lw=2, label="Light On")
    _individual_heading_means(df, ax1, 2.5)
    ax1.set_title("Population-Level μ (Centered)\nvia Circular KDE", va="bottom", pad=30)
    max_density = max(np.max(kde_mu_off), np.max(kde_mu_on))
    ax1.set_rlim(0, 1 + max_density * 1.5)
    ax1.set_rticks([1, 2, 3], [0, 1, 2])

    # (B) Posterior distribution of kappa.
    ax2 = fig.add_subplot(1, 5, 2)
    sns.violinplot(x="Parameter", y="Value", hue="Light",
                   data=df_kappa, split=True, inner="quartile",
                   palette={"Light Off": "red", "Light On": "blue"},
                   ax=ax2)
    ax2.set_title("Posterior Distribution of κ (Kappa)", pad=30)
    ax2.set_xlabel("")
    ax2.set_ylabel("Kappa")
    ax2.set_xticks([])

    # (C) Polar KDE of the population mu driven by kappa-mu.
    ax3 = fig.add_subplot(1, 5, 3, projection="polar")
    ax3.set_theta_zero_location("N")
    ax3.set_theta_direction(-1)
    ax3.plot(theta, kde_kappamu_off + 1, color="red", lw=2, label="Light Off")
    ax3.plot(theta, kde_kappamu_on + 1, color="blue", lw=2, label="Light On")
    _individual_heading_means(df, ax3, 2.5)
    ax3.set_title("Population-Level μ (Centered)\nvia Circular KDE using $κ_μ$", va="bottom", pad=30)
    max_density = max(np.max(kde_kappamu_off), np.max(kde_kappamu_on))
    ax3.set_rlim(0, 1 + max_density * 1.5)
    ax3.set_rticks([1, 2, 3], [0, 1, 2])

    # (D) Posterior distribution of kappa-mu.
    ax4 = fig.add_subplot(1, 5, 4)
    sns.violinplot(x="Parameter", y="Value", hue="Light",
                   data=df_kappamu, split=True, inner="quartile",
                   palette={"Light Off": "red", "Light On": "blue"},
                   ax=ax4)
    ax4.set_title("Posterior Distribution\n $κ_μ$ (Kappa of mu)", pad=30)
    ax4.set_xlabel("")
    ax4.set_ylabel("Kappa of mu")
    ax4.set_xticks([])

    # (E) Posterior distribution of the kappa-mu difference.
    ax5 = fig.add_subplot(1, 5, 5)
    sns.violinplot(data=pd.DataFrame(kappamu_delta_posterior),
                   split=False, inner="quartile",
                   color="darkgreen",
                   ax=ax5)
    ax5.set_title("Posterior Distribution\n $Diff. κ_μ$ (Kappa of mu)", pad=30)
    ax5.set_xlabel("")
    ax5.set_ylabel("Diff. Kappa of mu")
    ax5.set_xticks([])
    ax5.set_ylim(-5, 15)

    # One common legend for the whole figure.
    for ax in (ax1, ax2, ax3, ax4):
        if ax.get_legend() is not None:
            ax.get_legend().remove()
    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels, title="Light Condition", loc="lower right")

    fig.tight_layout(rect=[0, 0, 0.85, 1])  # leave space for the legend
    _plots.save_figure(fig, FIG_DIR, "post_summary")
    return fig
