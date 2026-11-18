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
# Two concentration-only moves improve tail mixing at a fixed draw count.
KAPPA_UPDATES = 2
KAPPA_STEP = 1.6

BLUE = "#0072b2"
GRAY = "#6b7280"
INK = "#231f20"
GRID = "#d4d7d9"
CHAIN_COLORS = (("#973e08", "#d55e00", "#e98b43"),
                ("#804181", "#b65fa0", "#d486bc"))
WRAP_BRANCHES = (-1, 0, 1)


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
    prior_seed, *chain_seeds = np.random.SeedSequence(SEED).spawn(
        1 + len(WRAP_BRANCHES))
    rng = np.random.default_rng(prior_seed)
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

    posterior = np.empty((len(WRAP_BRANCHES), N_DRAWS, 2))
    for chain, (branch, seed) in enumerate(zip(WRAP_BRANCHES, chain_seeds)):
        rng = np.random.default_rng(seed)
        eta = np.array([2 * np.pi * branch, 0.0])
        logp = log_posterior(eta)
        for i in range(BURN_IN + N_DRAWS):
            proposal = eta + rng.normal(scale=(0.42, 0.78))
            # Symmetric jumps let the chains visit equivalent wraps.
            if rng.random() < 0.12:
                proposal[0] += rng.choice((-2 * np.pi, 2 * np.pi))
            candidate = log_posterior(proposal)
            if np.log(rng.random()) < candidate - logp:
                eta, logp = proposal, candidate
            for _ in range(KAPPA_UPDATES):
                proposal = eta.copy()
                proposal[1] += rng.normal(scale=KAPPA_STEP)
                candidate = log_posterior(proposal)
                if np.log(rng.random()) < candidate - logp:
                    eta, logp = proposal, candidate
            if i >= BURN_IN:
                posterior[chain, i - BURN_IN] = eta
    theta = np.linspace(-np.pi, np.pi, 361)
    prior_density = _density(theta, *_parameters(prior))
    posterior_density = _density(theta, *_parameters(posterior.reshape(-1, 2)))
    return theta, observed, posterior, prior_density, posterior_density


def _page_axes(fig, page, x, y, w, h, **kwargs):
    width, height = page
    return fig.add_axes([x / width, (height - y - h) / height,
                         w / width, h / height], **kwargs)


def _trace_windows(draws):
    """Choose contiguous trace excerpts in three different periodic modes."""
    windows = []
    for chain, branch in zip(draws, WRAP_BRANCHES):
        wraps = np.floor((chain[:, 0] + np.pi) / (2 * np.pi))
        runs = np.split(chain, np.flatnonzero(np.diff(wraps)) + 1)
        candidates = [run for run in runs
                      if np.floor((run[0, 0] + np.pi) / (2 * np.pi)) == branch]
        windows.append(max(candidates, key=len))
    count = min(60, *(len(window) for window in windows))
    return np.stack([window[:count] for window in windows])


def _predictive_plot(fig, page, family, theta, observed, prior, posterior):
    polar = _page_axes(fig, page, 270, 94, 63, 63, projection="polar")
    polar.set_theta_zero_location("N")
    polar.set_theta_direction(-1)
    polar.set_xticks([])
    polar.set_yticks([])
    polar.spines["polar"].set_color(GRID)
    polar.spines["polar"].set_linewidth(0.6)
    polar.set_title("Predictive density", fontfamily=family, fontsize=7,
                    color=INK, pad=2.5)

    # An additive radius keeps zero density on an inner circle rather than
    # at the origin. Both sets of draws use one radial scale.
    prior_lines = prior[np.linspace(0, len(prior) - 1, 20, dtype=int)]
    posterior_lines = posterior[np.linspace(0, len(posterior) - 1, 24, dtype=int)]
    scale = max(prior_lines.max(), posterior_lines.max())
    polar.set_ylim(0, 2.17)
    polar.plot(theta, np.ones_like(theta), color=GRID, lw=0.55)
    for angle in np.arange(4) * np.pi / 2:
        polar.plot([angle, angle], [1, 2.17], color=GRID, lw=0.45)
    for sample, color, opacity in ((prior_lines, GRAY, 0.32),
                                    (posterior_lines, BLUE, 0.43)):
        polar.plot(theta, 1 + sample.T / scale, color=color, alpha=opacity,
                   lw=0.48)
    middle = np.median(posterior_lines, axis=0)
    radii = 1 + np.interp(observed, theta, middle) / scale
    for angle in np.unique(observed):
        indices = np.flatnonzero(observed == angle)
        radii[indices] += 0.09 * (np.arange(len(indices)) - (len(indices) - 1) / 2)
    polar.scatter(observed, radii, s=5.8, color=BLUE, edgecolors="white",
                  linewidths=0.33, zorder=5)

    handles = [
        Line2D([], [], color=GRAY, lw=0.8, alpha=0.65, label="Prior"),
        Line2D([], [], color=BLUE, lw=0.8, label="Posterior"),
        Line2D([], [], color=BLUE, marker="o", ms=2.4, lw=0, label="Observed"),
    ]
    key = _page_axes(fig, page, 338, 100, 43, 51)
    key.set_axis_off()
    legend = key.legend(handles=handles, loc="center left", frameon=True,
                        edgecolor="#c7c7c7", facecolor="white", framealpha=1,
                        borderaxespad=0, borderpad=0.46,
                        handlelength=1.0, handletextpad=0.48,
                        labelspacing=0.48, prop={"family": family, "size": 6.2})
    legend.get_frame().set_linewidth(0.65)


def _sampling_diagram(fig, page, family, draws):
    flow = _page_axes(fig, page, 389, 26, 207, 126)
    flow.set(xlim=(0, 207), ylim=(126, 0))
    flow.set_axis_off()
    text_style = dict(fontfamily=family, color=INK, fontsize=6.7, ha="center")
    for x, title in ((31, "MCMC draws"), (100, "Link"), (169, "Likelihood parameters")):
        flow.text(x, 8, title, weight="bold", **text_style)
    flow.plot([2, 204], [63, 63], color=GRID, lw=0.55)

    traces = _trace_windows(draws)
    mu_draws, _ = _parameters(traces)
    mu_center = np.angle(np.exp(1j * circ_mod(draws[:, :, 0])).mean())
    mean_points = np.argmin(np.abs(circ_mod(mu_draws - mu_center)), axis=1)
    kappa_traces = traces[:, :, 1]
    targets = [np.quantile(trace, q)
               for trace, q in zip(kappa_traces, (0.20, 0.50, 0.80))]
    kappa_points = np.argmin(np.abs(kappa_traces - np.array(targets)[:, None]), axis=1)
    x_trace = np.linspace(5, 57, traces.shape[1])

    def mark(x, y, color, *, size=10):
        flow.scatter([x], [y], s=size, color=color, edgecolors="white",
                     linewidths=0.35, zorder=5)

    def arrow(x0, x1, y, color="#848484"):
        flow.annotate("", (x1, y), (x0, y), arrowprops={
            "arrowstyle": "->", "color": color, "lw": 0.7, "mutation_scale": 6})

    for i, (top, trace_label, link, support_label, points) in enumerate((
        (0, "μ", "Circular modulo", "Mean direction μ", mean_points),
        (63, "κ", "Softplus", "Concentration κ", kappa_points),
    )):
        flow.text(29, top + 23, "η", **text_style)
        flow.text(33, top + 24.6, trace_label, fontfamily=family, color=INK,
                  fontsize=4.8, ha="center")
        flow.text(100, top + 23, link, **text_style)
        flow.text(169, top + 23, support_label, **text_style)
        arrow(60, 72, top + 42)
        arrow(128, 140, top + 42)

        trace = traces[:, :, i]
        low, high = trace.min(), trace.max()
        y_trace = top + 55 - 24 * (trace - low) / (high - low)
        values = trace[np.arange(len(points)), points]
        for color, y, point in zip(CHAIN_COLORS[i], y_trace, points):
            flow.plot(x_trace, y, color=color, alpha=0.72, lw=0.6)
            mark(x_trace[point], y[point], color)

        left, width, baseline, height = 75, 50, top + 54, 24
        if i == 0:
            for branch in WRAP_BRANCHES:
                start = (2 * branch - 1) * np.pi
                eta = np.linspace(start, start + 2 * np.pi, 100, endpoint=False)
                flow.plot(left + width * (eta + 3 * np.pi) / (6 * np.pi),
                          baseline - height * (circ_mod(eta) + np.pi) / (2 * np.pi),
                          color="#bb6a3a", lw=0.8)
            flow.plot([left, left + width], [baseline - height / 2] * 2,
                      color=GRID, lw=0.5, zorder=0)
            mapped = circ_mod(values)
            for color, eta, mu in zip(CHAIN_COLORS[i], values, mapped):
                mark(left + width * (eta + 3 * np.pi) / (6 * np.pi),
                     baseline - height * (mu + np.pi) / (2 * np.pi), color)
            circle = np.linspace(-np.pi, np.pi, 181)
            flow.plot(169 + 10 * np.sin(circle), 41 - 10 * np.cos(circle),
                      color=GRID, lw=0.7)
            for radius, mu, color in zip((7, 10, 13), mapped, CHAIN_COLORS[i]):
                mark(169 + radius * np.sin(mu), 41 - radius * np.cos(mu), color)
            flow.text(169, 59, "[−π, π)", **text_style)
        else:
            eta = np.linspace(-4, max(7, high), 160)
            support_limit = softplus(eta[-1])
            x_range = (left, left + width)
            y_range = (baseline, baseline - height)
            input_limits = (eta[0], eta[-1])
            support_limits = (0, support_limit)
            flow.plot([left, left + width], [baseline] * 2, color=GRID, lw=0.5)
            flow.plot(np.interp(eta, input_limits, x_range),
                      np.interp(softplus(eta), support_limits, y_range),
                      color="#ad5b9a", lw=0.8)
            mapped = softplus(values)
            for color, raw, kappa in zip(CHAIN_COLORS[i], values, mapped):
                mark(np.interp(raw, input_limits, x_range),
                     np.interp(kappa, support_limits, y_range), color)
            arrow(147, 192, 107, GRID)
            flow.plot([148, 148], [104, 110], color="#848484", lw=0.6)
            for color, kappa in zip(CHAIN_COLORS[i], mapped):
                mark(148 + 41 * kappa / support_limit, 107, color)
            flow.text(169, 121, "κ > 0", **text_style)


def draw_model_example(fig, page, family):
    """Draw individual predictive samples and the linked chain excerpts."""
    theta, observed, draws, prior, posterior = example()
    _predictive_plot(fig, page, family, theta, observed, prior, posterior)
    _sampling_diagram(fig, page, family, draws)
