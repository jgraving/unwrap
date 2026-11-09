# Zebra escape: figure construction.
#
# Sourced by notebooks_r/zebra_escape.ipynb (no standalone output). R twin
# of notebooks/zebra_plots.py: the notebook keeps the fit, the diagnostic
# values and the KDE grids; only the plot construction lives here.

trace_plot <- function(fit, trace_parameters, circular_parameters) {
  # as_draws_df is wide: explicitly assemble a long trace table.
  draws_df <- posterior::as_draws_df(posterior::subset_draws(
      posterior::as_draws_array(fit), variable = trace_parameters))
  trace_df <- do.call(rbind, lapply(trace_parameters, function(p) {
    value <- as.numeric(draws_df[[p]])
    if (p %in% circular_parameters) value <- mod_circular(value - circ_mean(value))
    data.frame(parameter = p, value = value,
               iteration = draws_df[[".iteration"]],
               chain = factor(draws_df[[".chain"]]))
  }))
  p_trace <- ggplot(trace_df, aes(iteration, value, colour = chain)) +
    geom_line(alpha = 0.6, linewidth = 0.25) +
    facet_wrap(~ parameter, scales = "free_y", ncol = 2) +
    labs(x = "retained iteration", y = "draw (angular traces centered)") +
    theme(legend.position = "bottom")
  p_trace
}

# Polar direction posteriors (twin of zebra_plots.py's polar panels): the
# plotted radius is 1 + kde/max_kde with 0 degrees at the top and the angle
# axis clockwise; empirical circular means sit at fixed radius 2.25.
polar_direction_plot <- function(kde_df, empirical, x_lab, condition_colors,
  track = FALSE) {
  polar_x <- scale_x_continuous(
    limits = c(0, 2 * pi), expand = c(0, 0),
    breaks = (0:7) * pi / 4,
    labels = c("0°", "45°", "90°", "135°", "180°", "225°", "270°", "315°"))
  polar_y <- scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0),
    breaks = c(1.0, 2.0), labels = c("0", "1.0"))
  polar_theme <- theme(legend.position = "bottom",
    axis.line = element_blank(), panel.border = element_blank())
  polar_guides <- guides(shape = "none")
  p <- ggplot(kde_df, aes(theta, r, colour = condition)) +
    geom_line(linewidth = if (track) 0.3 else 0.5,
      alpha = if (track) 0.5 else 1) +
    geom_point(data = empirical, aes(theta, r, colour = condition,
      shape = shape), inherit.aes = FALSE, size = if (track) 0.9 else 1.2) +
    coord_polar(start = pi / 2, direction = -1) + polar_x + polar_y +
    scale_colour_manual(values = condition_colors) +
    labs(x = x_lab, y = "") + polar_theme + polar_guides
  if (track) p <- p + facet_wrap(~ Track)
  p
}

# Between-individual alignment: mean R with 95% HDI per condition.
alignment_plot <- function(alignment_df, condition_colors) {
  p_alignment <- ggplot(alignment_df, aes(condition, mean, colour = condition)) +
    geom_pointrange(aes(ymin = low, ymax = high)) +
    scale_colour_manual(values = condition_colors) +
    labs(x = "", y = "between-individual alignment (95% HDI)")
  p_alignment
}

# Individual posterior summaries: paired mean direction and within-individual
# R(kappa) by track and condition (twin of zebra_plots.py's individual
# panels).
individual_plot <- function(individual_summary, condition_colors) {
  p_individual <- ggplot(individual_summary,
    aes(id, mean, colour = condition, group = condition)) +
    geom_pointrange(aes(ymin = low, ymax = high),
      position = position_dodge(width = 0.5), linewidth = 0.25) +
    facet_wrap(~ quantity, scales = "free_y", ncol = 1) +
    scale_colour_manual(values = condition_colors) +
    labs(x = "individual", y = "posterior mean and 95% HDI") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "bottom")
  p_individual
}

# Supplemental PPC: data bootstrap vs posterior, group level (open = data
# bootstrap, filled = posterior).
group_ppc_plot <- function(ppc_group, condition_colors) {
  p_group_ppc <- ggplot(ppc_group,
    aes(condition, mean, colour = condition, shape = source)) +
    geom_pointrange(aes(ymin = low, ymax = high),
      position = position_dodge(width = 0.35)) +
    facet_wrap(~ quantity, scales = "free_y") +
    scale_colour_manual(values = condition_colors) +
    scale_shape_manual(values = c("data bootstrap" = 1, "posterior" = 16)) +
    labs(x = "", y = "mean and 95% HDI") + theme(legend.position = "bottom")
  p_group_ppc
}

# Supplemental PPC: individual level, data bootstrap rows and the posterior
# rows from individual_plot's summary table.
individual_ppc_plot <- function(ppc_individual) {
  p_ind_ppc <- ggplot(ppc_individual,
    aes(id, mean, colour = source, shape = source)) +
    geom_pointrange(aes(ymin = low, ymax = high),
      position = position_dodge(width = 0.4), linewidth = 0.2) +
    facet_grid(quantity ~ condition, scales = "free_y") +
    scale_shape_manual(values = c("data bootstrap" = 1, "posterior" = 16)) +
    labs(x = "individual", y = "mean and 95% HDI") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
      legend.position = "bottom")
  p_ind_ppc
}

# Conventional replicated-data check: replicated mean resultant length per
# track/condition with the observed R as an open marker.
replicated_r_plot <- function(rep_stats, condition_colors) {
  p_replicates <- ggplot(rep_stats, aes(id, mean, colour = condition)) +
    geom_pointrange(aes(ymin = low, ymax = high)) +
    geom_point(aes(y = observed), shape = 1) +
    facet_wrap(~ condition) + scale_colour_manual(values = condition_colors) +
    labs(x = "individual", y = "replicated R, 95% HDI; open = observed") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
      legend.position = "bottom")
  p_replicates
}
