# Circular regression: figure construction.
#
# Sourced by notebooks_r/circular_regression.ipynb (no standalone
# output). R twin of notebooks/circular_regression_plots.py: the
# notebook keeps the fit, the diagnostic values and the posterior
# computation; only the plot construction lives here.

# Posterior traces for the parameters: one line per chain against the
# post-warmup iteration, one panel per parameter (the Python twin's
# arviz trace plot).
trace_plot <- function(fit, parameters) {
  draws_df <- posterior::as_draws_df(fit, variables = parameters)
  trace_df <- do.call(rbind, lapply(seq_along(parameters), function(k)
    data.frame(parameter = parameters[k],
      iteration = as.integer(draws_df[[".iteration"]]),
      value = as.numeric(draws_df[[parameters[k]]]),
      chain = factor(as.integer(draws_df[[".chain"]])))))
  ggplot(trace_df, aes(x = iteration, y = value, colour = chain)) +
    geom_line(alpha = 0.4, linewidth = 0.3) +
    facet_wrap(~ parameter, scales = "free_y") +
    scale_color_discrete(name = "chain") +
    labs(x = "iteration", y = "") +
    theme(legend.position = "bottom")
}

# Post-warmup predictions: the posterior mean direction (averaged as
# unit vectors, so angles near +/-pi are neighbours) and the mean
# concentration across draws, with the observations and the known
# generating curves dashed underneath.
predictions_plot <- function(data, grid_x, mean_direction, concentration) {
  obs_df <- data.frame(grid = data$x, value = data$y,
    panel = "mean direction (rad)")
  pred_df <- data.frame(
    grid = rep(grid_x, 2),
    value = c(mean_direction, concentration),
    panel = c(rep("mean direction (rad)", length(grid_x)),
      rep("concentration (kappa)", length(grid_x))))
  ggplot(pred_df, aes(grid, value)) +
    geom_line(alpha = 0.8) +
    geom_point(data = obs_df, aes(grid, value),
      alpha = 0.35, position = position_jitter(width = 0.02, height = 0.15)) +
    geom_line(data = data.frame(grid = grid_x,
      value = mod_circular(grid_x),
      panel = "mean direction (rad)"),
      aes(grid, value), inherit.aes = FALSE,
      linetype = "dashed") +
    geom_line(data = data.frame(grid = grid_x,
      value = softplus(2 + 0.5 * grid_x),
      panel = "concentration (kappa)"),
      aes(grid, value), inherit.aes = FALSE,
      linetype = "dashed") +
    facet_wrap(~ panel, scales = "free_y") +
    labs(x = expression(x), y = "")
}

# Posterior predictive check: the density of the posterior replicates'
# mean resultant length against the observed one.
ppc_plot <- function(replicated_R, observed_R) {
  ggplot(data.frame(replicated_R), aes(replicated_R)) +
    geom_density(alpha = 0.7, fill = "darkseagreen", colour = "grey40") +
    geom_vline(xintercept = observed_R, colour = "firebrick", linewidth = 1) +
    labs(x = "mean resultant length R", y = "density")
}
