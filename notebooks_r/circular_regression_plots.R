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

bearing_degrees <- function(theta) (pi/2 + mod_circular(theta - pi/2)) * 180/pi

# Draws occupy rows. Center angles before taking quantiles so a band can
# cross the signed-angle boundary at 180 degrees.
direction_band <- function(x, draws) {
  center <- Arg(colMeans(exp(1i * draws)))
  deviations <- mod_circular(sweep(draws, 2, center, "-"))
  offsets <- apply(deviations, 2, quantile, probs = c(0.05, 0.95)) * 180/pi
  value <- bearing_degrees(center)
  data.frame(x, value, low = value + offsets[1, ], high = value + offsets[2, ])
}

linear_band <- function(x, draws) {
  q <- apply(draws, 2, quantile, probs = c(0.05, 0.5, 0.95))
  data.frame(x, value = q[2, ], low = q[1, ], high = q[3, ])
}

cue_facets <- function(p, labels) {
  p + facet_wrap(~ panel, scales = "free_y", labeller = as_labeller(labels)) +
    scale_x_continuous(breaks = seq(-1, 1, by = 0.5)) +
    labs(x = "Cue strength (x)", y = NULL) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(), strip.text = element_text(face = "bold"))
}

predictions_plot <- function(data, grid_x, pred, true_mu, true_kappa) {
  mu <- transform(direction_band(grid_x, pred$mu), panel = "A", truth = bearing_degrees(true_mu))
  kappa <- transform(linear_band(grid_x, pred$kappa), panel = "B", truth = true_kappa)
  curves <- rbind(mu, kappa)
  observed <- data.frame(x = data$x, value = bearing_degrees(data$y), panel = "A")
  limits <- data.frame(x = 0, value = c(-90, 270, 0), panel = c("A", "A", "B"))
  p <- ggplot(curves, aes(x, value)) +
    geom_blank(data = limits) +
    geom_point(data = observed, size = 1, alpha = 0.24, colour = "grey35") +
    geom_ribbon(aes(ymin = low, ymax = high), fill = "#086b72", alpha = 0.24) +
    geom_line(aes(colour = "Fitted"), linewidth = 0.9) +
    geom_line(aes(y = truth, colour = "Generating curve"), linetype = 2, linewidth = 0.7) +
    scale_colour_manual(values = c("Fitted" = "#086b72", "Generating curve" = "#a33b12"),
                        name = NULL) +
    labs(caption = "Grey points: individual headings. Shading: pointwise 90% credible band.")
  cue_facets(p, c(A = "Preferred bearing (degrees)", B = "Concentration (kappa)"))
}

ppc_plot <- function(levels, observed, predicted) {
  mu <- transform(direction_band(levels, Arg(predicted)), panel = "A",
                  observed = bearing_degrees(Arg(observed)))
  alignment <- transform(linear_band(levels, Mod(predicted)), panel = "B", observed = Mod(observed))
  intervals <- rbind(mu, alignment)
  limits <- data.frame(x = 0, value = c(-90, 270, 0, 1), panel = c("A", "A", "B", "B"))
  p <- ggplot(intervals, aes(x, value)) +
    geom_blank(data = limits) +
    geom_linerange(aes(ymin = low, ymax = high), colour = "#086b72", linewidth = 1.2) +
    geom_point(aes(shape = "Predictive center"), size = 2, colour = "#086b72") +
    geom_point(aes(y = observed, shape = "Observed group"), size = 2.6,
               fill = "#a33b12", colour = "white", stroke = 0.3) +
    scale_shape_manual(values = c("Predictive center" = 16, "Observed group" = 23),
                       name = NULL) +
    labs(caption = "Bars enclose 90% of replicated groups at the same cue strength.")
  cue_facets(p, c(A = "Group mean bearing (degrees)", B = "Group mean resultant length"))
}
