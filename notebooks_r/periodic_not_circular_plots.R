# Periodic (not circular) beta: figure construction.
#
# Sourced by notebooks_r/periodic_not_circular.ipynb (no standalone
# output). R twin of notebooks/periodic_not_circular_plots.py.

trace_plot <- function(fit, parameters, angular = character()) {
  d <- draw_matrix(fit, parameters)
  long <- do.call(rbind, lapply(parameters, function(p) {
    x <- d[[p]]
    if (p %in% angular) x <- wrap(x - circ_mean(x))
    data.frame(parameter = p, value = x, iteration = d[[".iteration"]],
               chain = factor(d[[".chain"]]))
  }))
  ggplot(long, aes(iteration, value, colour = chain, group = chain)) +
    geom_line(linewidth = 0.25, alpha = 0.6) +
    facet_wrap(~ parameter, scales = "free_y", ncol = 2) +
    labs(x = "retained iteration", y = "draw (angular traces centered)") +
    theme(legend.position = "bottom")
}

curve_table <- function(draws, truth, panel) {
  bounds <- t(apply(draws, 1, hdi))
  data.frame(time, mean = rowMeans(draws), low = bounds[, 1], high = bounds[, 2], truth, panel)
}

# Posterior latent curves: the fitted mean activity and precision phi
# against time with the 95% HDI ribbon and the simulated truth dashed.
latent_curves_plot <- function(time, mu_draws, phi_draws, mu_true, phi_true) {
  latent <- rbind(curve_table(mu_draws, mu_true, "mean activity"),
    curve_table(phi_draws, phi_true, "precision phi"))
p_latent <- ggplot(latent, aes(time, mean)) +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "#E69F00", alpha = 0.2) +
  geom_line() + geom_line(aes(y = truth), linetype = "dashed", colour = "#0072B2") +
  facet_wrap(~ panel, scales = "free_y") +
  labs(x = "Time (hours)", y = "fitted mean and 95% HDI; dashed = simulated truth")
  p_latent
}

# Posterior predictive fit: the observed 24-hour activity with the
# posterior predictive median and 50%/90% HDI bands (the Python twin
# draws the replicates inside its plot function the same way; quick-run
# budgets cap the draw count as before).
ppc_fit <- function(time, activity, mu_draws, phi_draws, quick_run) {
use <- seq_len(ncol(mu_draws))
if (quick_run) use <- use[seq_len(min(length(use), 1000L))]
  G <- length(time)
yrep <- matrix(0, nrow = G, ncol = length(use))
for (j in seq_along(use)) {
  k <- use[j]
  yrep[, j] <- rbeta(G, mu_draws[, k] * phi_draws[, k],
                        (1 - mu_draws[, k]) * phi_draws[, k])
}
b50 <- t(apply(yrep, 1, hdi, prob = 0.5))
b90 <- t(apply(yrep, 1, hdi, prob = 0.9))
ppc <- data.frame(time, median = apply(yrep, 1, median),
  low50 = b50[, 1], high50 = b50[, 2], low90 = b90[, 1], high90 = b90[, 2])
p_ppc <- ggplot(ppc, aes(time, median)) +
  geom_ribbon(aes(ymin = low90, ymax = high90), fill = "#E69F00", alpha = 0.18) +
  geom_ribbon(aes(ymin = low50, ymax = high50), fill = "#E69F00", alpha = 0.35) +
  geom_point(data = data.frame(time, activity), aes(time, activity),
    inherit.aes = FALSE, size = 0.6, alpha = 0.35) +
  geom_line(linewidth = 0.6) + scale_x_continuous(breaks = c(0, 6, 12, 18, 24)) +
  labs(x = "Time (hours)", y = "Activity")
  p_ppc
}