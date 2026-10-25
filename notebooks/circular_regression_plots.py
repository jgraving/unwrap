"""Figure helpers for the circular_regression example notebook.

The notebook focuses on the statistical ideas (what we fit and why); the
matplotlib "how" (axes, lines, legends, layout) lives here. The notebook
imports this sibling module and calls the functions below, passing the values
it has already computed.
"""
import numpy as np
import matplotlib.pyplot as plt
import arviz as az


def plot_trace(idata, var_names):
    """Compact ArviZ trace plot of the posterior for the named variables."""
    az.plot_trace(idata, var_names=var_names, compact=True)
    plt.tight_layout()
    plt.show()


def plot_predictions(data, grid, mean_direction, concentration_draws):
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
    plt.show()


def plot_ppc(replicated_R, observed_R):
    """Histogram of replicated-dataset mean resultant lengths against the
    observed value (a circular-concentration-style PPC)."""
    fig, ax = plt.subplots(figsize=(6, 3))
    ax.hist(replicated_R.ravel(), bins=30, alpha=0.7, label="Replicated datasets")
    ax.axvline(observed_R, color="black", label="Observed dataset")
    ax.set(xlabel="Mean resultant length", ylabel="Replicated dataset count")
    ax.legend()
    fig.tight_layout()
    plt.show()
