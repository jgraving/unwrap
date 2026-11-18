"""Reproducible circular-model example for the graphical abstract."""
from functools import lru_cache

import numpy as np
from matplotlib.lines import Line2D
from scipy.special import i0e

from unwrap.circular import circ_mod
from unwrap.link.numpy import softplus


SEED = 159
# Dot centres measured on the original panel's −180° to 180° axis.
OBSERVED_DEG = (-76, -55, -55, -33, -33, -33, -12, -12,
                10, 10, 10, 31, 31, 53, 53, 75)
PRIOR_SD = (3.5, 2.5)
BURN_IN = 600
N_DRAWS = 1000

BLUE = "#0072b2"
GRAY = "#6b7280"
INK = "#231f20"
GRID = "#d4d7d9"
LINK_COLORS = ("#d55e00", "#bd69ac")


def _parameters(eta):
    eta_mu, eta_kappa = np.moveaxis(eta, -1, 0)
    return circ_mod(eta_mu), softplus(eta_kappa)


def _density(theta, mu, kappa):
    kappa = kappa[..., None]
    return np.exp(kappa * (np.cos(theta - mu[..., None]) - 1.0)) / (
        2 * np.pi * i0e(kappa))


@lru_cache(maxsize=1)
def example():
    """Sample an intercept-only von Mises model on the real predictor scale."""
    rng = np.random.default_rng(SEED)
    observed = np.deg2rad(OBSERVED_DEG)
    prior = rng.normal(scale=PRIOR_SD, size=(N_DRAWS, 2))
    resultant = np.exp(1j * observed).sum()

    def log_posterior(eta):
        mu, kappa = _parameters(eta)
        # Subtract κ before evaluating the log normalizer for stability.
        log_likelihood = kappa * (
            (resultant * np.exp(-1j * mu)).real - observed.size
        ) - observed.size * np.log(i0e(kappa))
        return log_likelihood - 0.5 * np.square(eta / PRIOR_SD).sum()

    draws = np.empty((BURN_IN + N_DRAWS, 2))
    eta = np.zeros(2)
    logp = log_posterior(eta)
    for i in range(len(draws)):
        proposal = eta + rng.normal(scale=(0.42, 0.78))
        # A symmetric jump between equivalent angles helps traverse wraps.
        if rng.random() < 0.12:
            proposal[0] += rng.choice((-2 * np.pi, 2 * np.pi))
        candidate = log_posterior(proposal)
        if np.log(rng.random()) < candidate - logp:
            eta, logp = proposal, candidate
        draws[i] = eta
    posterior = draws[BURN_IN:]
    theta = np.linspace(-np.pi, np.pi, 361)
    prior_density = _density(theta, *_parameters(prior))
    posterior_density = _density(theta, *_parameters(posterior))
    return theta, observed, posterior, prior_density, posterior_density


def _page_axes(fig, page, x, y, w, h, **kwargs):
    width, height = page
    return fig.add_axes([x / width, (height - y - h) / height,
                         w / width, h / height], **kwargs)


def draw_model_example(fig, page, family):
    """Draw predictive polar curves and a two-lane sampling diagram."""
    theta, observed, draws, prior, posterior = example()
    label_kw = dict(fontfamily=family, color=INK, fontsize=7.0)
    polar = _page_axes(fig, page, 270, 86, 68, 68, projection="polar")
    polar.set_theta_zero_location("N")
    polar.set_theta_direction(-1)
    polar.set_xticks(np.arange(4) * np.pi / 2)
    polar.set_xticklabels([])
    polar.set_yticks([])
    polar.spines["polar"].set_color(GRID)
    polar.spines["polar"].set_linewidth(0.65)
    polar.xaxis.grid(True, color=GRID, linewidth=0.45)
    polar.set_axisbelow(True)
    upper = 0.0
    for sample, color in ((prior, GRAY), (posterior, BLUE)):
        low, high = np.quantile(sample, (0.1, 0.9), axis=0)
        upper = max(upper, high.max())
        polar.fill_between(theta, low, high, color=color, alpha=0.22,
                           linewidth=0)
        polar.plot(theta, sample.mean(axis=0), color=color, lw=1.05)
    radii = np.empty(observed.size)
    for angle in np.unique(observed):
        indices = np.flatnonzero(observed == angle)
        radii[indices] = upper * (1.24 - 0.095 * np.arange(indices.size))
    polar.set_ylim(0, 1.34 * upper)
    polar.scatter(observed, radii, s=5.0, color=BLUE, edgecolors="none", zorder=5)
    handles = [
        Line2D([], [], color=GRAY, lw=1.2, label="Prior"),
        Line2D([], [], color=BLUE, lw=1.2, label="Posterior"),
        Line2D([], [], color=BLUE, marker="o", ms=2.7, lw=0, label="Observed angles"),
    ]
    key = _page_axes(fig, page, 342, 94, 45, 50)
    key.set_axis_off()
    key.legend(handles=handles, loc="center left", frameon=False,
               title="Predictive density", alignment="left", title_fontsize=6.2,
               borderaxespad=0, handlelength=1.05, handletextpad=0.46,
               labelspacing=0.82, prop={"family": family, "size": 6.05})

    flow = _page_axes(fig, page, 393, 28, 200, 123)
    flow.set(xlim=(0, 200), ylim=(123, 0))
    flow.set_axis_off()
    for x, title in ((30, "MCMC draws"), (94, "Link"), (161, "Parameter support")):
        flow.text(x, 7, title, ha="center", weight="bold", **label_kw)
    flow.text(30, 15, "Unconstrained", ha="center", color=GRAY,
              fontfamily=family, fontsize=6.05)
    flow.plot([0, 198], [60, 60], color=GRID, lw=0.55)
    flow.text(100, 122, "Illustrative fit", ha="center", color=GRAY,
              fontfamily=family, fontsize=6.2)

    indices = np.unique(np.linspace(25, N_DRAWS - 10, 20, dtype=int))
    mapped = _parameters(draws[indices])
    for col, (y, color, eta_label, link, parameter, support) in enumerate((
        (31, LINK_COLORS[0], r"$\eta_\mu$", "Circular modulo",
         r"Direction $\mu$", r"$[-\pi,\,\pi)$"),
        (83, LINK_COLORS[1], r"$\eta_\kappa$", "Softplus",
         r"Concentration $\kappa$", r"$\kappa > 0$"),
    )):
        flow.text(30, y - 8, eta_label, ha="center", color=color, fontsize=8)
        trace = draws[:, col]
        trace_y = y + 17 - 21 * (trace - trace.min()) / np.ptp(trace)
        x = np.linspace(2, 60, len(trace))
        flow.plot(x, trace_y, color=color, alpha=0.52, lw=0.55)
        flow.scatter(x[indices], trace_y[indices], color=color, s=3,
                     edgecolors="none", zorder=4)
        flow.text(94, y, link, ha="center", **label_kw)
        flow.annotate("", (120, y + 8), (69, y + 8),
                      arrowprops={"arrowstyle": "->", "color": color,
                                  "lw": 1.0, "mutation_scale": 7})
        flow.text(160, y - 8, parameter, ha="center", **label_kw)
        flow.text(160, y + 24, support, ha="center", **label_kw)
    # The same highlighted samples appear after each transformation.
    angle = np.linspace(-np.pi, np.pi, 181)
    flow.plot(160 + 10 * np.sin(angle), 35 - 10 * np.cos(angle),
              color=GRID, lw=0.65)
    mu, kappa = mapped
    flow.scatter(160 + 10 * np.sin(mu), 35 - 10 * np.cos(mu),
                 color=LINK_COLORS[0], s=5, alpha=0.64, edgecolors="none")
    flow.annotate("", (186, 94), (134, 94),
                  arrowprops={"arrowstyle": "->", "color": GRID, "lw": 0.75})
    flow.scatter(134 + 45 * kappa / kappa.max(),
                 92 + (np.arange(len(kappa)) % 3 - 1) * 1.5,
                 color=LINK_COLORS[1], s=5, alpha=0.64, edgecolors="none")
    flow.plot([134, 134], [91, 97], color=GRAY, lw=0.6)
