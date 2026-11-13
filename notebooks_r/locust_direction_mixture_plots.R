# Locust direction/mixture: figure construction.
#
# Sourced by notebooks_r/locust_direction_mixture_model.ipynb (no
# standalone output). R twin of notebooks/locust_direction_plots.py:
# the notebook keeps the data, the preprocessing, the fits (and their
# fit-cache names), the seeds and the printed post-hoc tables; only the
# figure construction lives here (plus the colour and axis-label
# constants and the validation-figure helpers the builders share).

BLUE <- "#1F77B4"
ORANGE <- "#FF7F0E"
RED <- "#D62728"
HYP_COLORS <- c("Optic flow" = BLUE, "Vec. encoding" = ORANGE)
angle_breaks <- c(-pi, -pi / 2, 0, pi / 2, pi)
angle_labels <- c("-pi", "-pi/2", "0", "pi/2", "pi")

model_labels <- c(Bimodal = "Bimodal\n(Vec. encoding)",
                  Unimodal = "Unimodal\n(Optic flow)")

tab20 <- c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C", "#98DF8A",
           "#D62728", "#FF9896", "#9467BD", "#C5B0D5", "#8C564B", "#C49C94",
           "#E377C2", "#F7B6D2", "#7F7F7F", "#C7C7C7", "#BCBD22", "#DBDB8D",
           "#17BECF", "#9EDAE5")

trajectory_panel <- function(limit, points = FALSE) {
  p <- ggplot(segments, aes(x, y, colour = factor(animal_id %% 20L))) +
    geom_segment(aes(xend = xend, yend = yend), linewidth = 0.25) +
    scale_colour_manual(values = tab20) +
    coord_equal(xlim = c(-limit, limit), ylim = c(-limit, limit)) +
    labs(x = "X Coordinate (meters)", y = "Y Coordinate (meters)") +
    theme(legend.position = "none")
  if (points) p <- p + geom_point(data = df,
    aes(X, Y, colour = factor(animal_id %% 20L)), inherit.aes = FALSE, size = 0.6)
  p
}

observed_polar_plot <- function(raw_densities) {
  # Same radius convention as the other polar figures: density zero sits on
  # the ring at radius 1 (the bars rise from it) and the shared peak of all
  # 52 histograms sits on the ring at radius 2; the 0-degree spoke marks the
  # 3 o'clock zero of the mathematical orientation.
  peak_max_hist <- max(raw_densities$density)
  raw_densities$radius <- 1 + raw_densities$density / peak_max_hist
  ggplot(raw_densities) +
    geom_rect(aes(xmin = left, xmax = right, ymin = 1, ymax = radius), fill = "grey35") +
    geom_hline(yintercept = c(1, 2), colour = "grey40", linewidth = 0.25) +
    geom_segment(data = data.frame(x = 0, xend = 0, y = 1, yend = 2.01),
      inherit.aes = FALSE,
      aes(x = x, xend = xend, y = y, yend = yend),
      colour = "grey20", linewidth = 0.25) +
    scale_x_continuous(limits = c(-pi, pi), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 2.05), expand = c(0, 0)) +
    coord_polar(start = pi / 2, direction = -1) +
    facet_wrap(~ individual, ncol = 7) + theme_void(base_size = 8)
}

trace_plot <- function(state, parameters) {
  long <- do.call(rbind, lapply(parameters, function(p) {
    x <- state$quantities[[p]]
    stopifnot(!is.null(x))
    if (grepl("^mu(\\[|_pop($|\\[))", p)) x <- wrap(x - circ_mean(x))
    data.frame(parameter = p, iteration = state$draws[[".iteration"]],
               chain = factor(state$draws[[".chain"]]), value = x)
  }))
  ggplot(long, aes(iteration, value, colour = chain, group = chain)) +
    geom_line(linewidth = 0.25, alpha = 0.65) +
    facet_wrap(~ parameter, scales = "free_y", ncol = 2) +
    labs(x = "retained iteration", y = "draw (angular traces centered)") +
    theme(legend.position = "bottom")
}

elpd_plot <- function(elpd_table) {
  ggplot(elpd_table, aes(elpd, factor(model, levels = c("Bimodal", "Unimodal")),
                                   colour = model, shape = model)) +
    scale_y_discrete(labels = model_labels) +
    geom_errorbar(aes(xmin = elpd - se, xmax = elpd + se), orientation = "y", width = 0.15) +
    geom_point(size = 2.5) + scale_colour_manual(values = c(Bimodal = ORANGE, Unimodal = BLUE)) +
    scale_shape_manual(values = c(Bimodal = 15, Unimodal = 16)) +
    labs(x = "PSIS-LOO ELPD +/-1 SE", y = "") + theme(legend.position = "none")
}

elpd_individual_plot <- function(elpd_long) {
  ggplot(elpd_long,
    aes(value, factor(individual, levels = individual_levels), colour = model, shape = model)) +
    geom_point(position = position_dodge(width = 0.5), size = 1.8) +
    scale_colour_manual(values = c(Bimodal = ORANGE, Unimodal = BLUE)) +
    scale_shape_manual(values = c(Bimodal = 15, Unimodal = 16)) +
    labs(x = "PSIS-LOO ELPD contribution", y = "Individual") + theme(legend.position = "bottom")
}

elpd_delta_plot <- function(elpd_individual) {
  ggplot(elpd_individual,
    aes(delta, factor(individual, levels = individual_levels))) +
    geom_point(size = 1.8) + geom_vline(xintercept = 0, linetype = "dashed") +
    labs(x = "Delta ELPD (Bimodal - Unimodal)", y = "")
}

model_layer <- function() {
  ggplot(model_density, aes(theta, center)) +
    geom_ribbon(aes(ymin = low, ymax = high), fill = RED, alpha = 0.25) +
    geom_line(colour = RED, linewidth = 0.6) +
    scale_x_continuous(breaks = angle_breaks, labels = angle_labels,
                       limits = c(-pi, pi), expand = c(0, 0)) +
    labs(x = "Direction (radians)", y = "Density")
}

# Hand-drawn Model/Histogram/Kernel legend, matching the Python patch+line
# legend handlers. Built with annotate() (no ggplot2 legend-key geoms) so it
# emits no warning streams; placed as a right-hand patchwork column.
val_legend <- function(framed = FALSE) {
  lp <- ggplot(data.frame(y = c(8, 6.5, 5)), aes(x = 1, y = 0)) +
    annotate("rect", xmin = 0.8, xmax = 1.4, ymin = 7.7, ymax = 8.3, fill = RED, alpha = 0.25) +
    annotate("segment", x = 0.8, y = 8, xend = 1.4, yend = 8, colour = RED, linewidth = 0.6) +
    annotate("point", x = 1.1, y = 6.5, size = 2.5, colour = "black") +
    annotate("rect", xmin = 0.8, xmax = 1.4, ymin = 4.7, ymax = 5.3, fill = "black", alpha = 0.15) +
    annotate("segment", x = 0.8, y = 5, xend = 1.4, yend = 5, colour = "black", linetype = "dashed", linewidth = 0.5) +
    scale_y_continuous(breaks = c(8, 6.5, 5), labels = c("Model", "Histogram", "Kernel"),
                       expand = c(0, 0)) +
    scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
    coord_cartesian(xlim = c(0.7, 1.5), ylim = c(0, 10), clip = "off") +
    theme_void() +
    theme(axis.text.y = element_text(size = 9, colour = "black", hjust = 0,
                                     margin = margin(0, 0, 0, 3, "pt")))
  if (framed) {
    lp <- lp + annotate("rect", xmin = 0.62, xmax = 1.58, ymin = 4.45, ymax = 8.55,
                        colour = "grey50", fill = NA, linewidth = 0.3)
  }
  lp
}

validation_cartesian_plot <- function() {
  p_hist <- model_layer() +
    geom_errorbar(data = hist_density, aes(ymin = low, ymax = high), width = 0.06) +
    geom_point(data = hist_density, size = 1.5)
  p_kde <- model_layer() +
    geom_ribbon(data = kde_density, aes(ymin = low, ymax = high), fill = "black", alpha = 0.15) +
    geom_line(data = kde_density, linetype = "dashed", linewidth = 0.5)
  p_validation_cartesian <- (p_hist / p_kde) | val_legend(framed = FALSE)
  p_validation_cartesian + plot_layout(widths = c(4, 1))
}

population_validation_polar <- function() {
  edge_radius <- 0.5
  p_validation_polar <- ggplot(model_density, aes(theta, edge_radius + center)) +
    geom_ribbon(aes(ymin = edge_radius + low, ymax = edge_radius + high), fill = RED, alpha = 0.25) +
    geom_line(colour = RED, linewidth = 0.6) +
    geom_ribbon(data = kde_density, aes(ymin = edge_radius + low, ymax = edge_radius + high),
                fill = "black", alpha = 0.15) +
    geom_line(data = kde_density, linetype = "dashed", linewidth = 0.5) +
    geom_errorbar(data = hist_density, aes(ymin = edge_radius + low, ymax = edge_radius + high), width = 0.02) +
    geom_point(data = hist_density, size = 1.3) +
    geom_hline(yintercept = edge_radius, colour = "black", linewidth = 0.5) +
    scale_x_continuous(limits = c(-pi, pi),
                       breaks = c(-pi / 2, 0, pi / 2, pi),
                       labels = c("-90°", "0°", " 90°", "±180°     "),
                       expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, edge_radius + max(model_density$high, kde_density$high) * 1.2)) +
    coord_polar(start = pi / 2, direction = -1) +
    # The twin shows no axis labels; keep blank radius tick labels like it.
    labs(x = "", y = "") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p_validation_polar <- p_validation_polar | val_legend(framed = TRUE)
  p_validation_polar + plot_layout(widths = c(4, 1))
}

individual_ppc_plot <- function(ppc_curves, ppc_ticks, ppc_arrows, rmax) {
  p_individual_ppc <- ggplot(ppc_curves, aes(theta, 1 + model, group = individual)) +
    geom_hline(yintercept = c(1, 2), colour = "grey40", linewidth = 0.25) +
    geom_vline(xintercept = c(-pi, -pi / 2, 0, pi / 2), colour = "grey70", linetype = "dashed", linewidth = 0.2) +
    geom_ribbon(aes(ymin = 1 + kde_low, ymax = 1 + kde_high), fill = "grey50", alpha = 0.25) +
    geom_line(aes(y = 1 + kde), colour = "grey20", linetype = "dashed", linewidth = 0.4) +
    geom_ribbon(aes(ymin = 1 + model_low, ymax = 1 + model_high), fill = RED, alpha = 0.35) +
    geom_line(colour = RED, linewidth = 0.5) +
    geom_segment(data = ppc_ticks,
      aes(x = theta, xend = theta, y = 1, yend = 1.1), inherit.aes = FALSE, linewidth = 0.2) +
    geom_segment(data = ppc_ticks,
      aes(x = theta, xend = theta, y = 1, yend = pmax(0, 1 + residual)),
      inherit.aes = FALSE, colour = "#17BECF", linewidth = 0.2) +
    geom_segment(data = ppc_arrows,
      aes(x = theta, xend = theta, y = 0, yend = length), inherit.aes = FALSE,
      colour = RED, linewidth = 0.4, arrow = grid::arrow(length = grid::unit(0.05, "inches"))) +
    scale_x_continuous(limits = c(-pi, pi), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, rmax), expand = c(0, 0)) +
    coord_polar(start = pi / 2, direction = -1) +
    facet_wrap(~ individual, ncol = 7) + theme_void(base_size = 8) +
    theme(strip.text = element_text(size = 8), panel.spacing = grid::unit(0.5, "lines"))

  # Figure legend matching the Python twin: one-vertex line per key, so the
  # legend entries come from the manual scales and nothing is drawn on the
  # panels (manual guide overrides would warn on dev ggplot2; this does not).
  ppc_key <- data.frame(theta = rep(0, 4), r = rep(0, 4),
    key = factor(rep(c("Kernel (bootstrap 95%)", "Model (posterior 95%)"), each = 2L),
                 levels = c("Kernel (bootstrap 95%)", "Model (posterior 95%)")))
  p_individual_ppc +
    geom_line(data = ppc_key,
      aes(x = theta, y = r, colour = key, linetype = key, group = key),
      inherit.aes = FALSE, linewidth = 0.4) +
    scale_colour_manual(name = NULL,
      values = c("Kernel (bootstrap 95%)" = "grey20", "Model (posterior 95%)" = RED)) +
    scale_linetype_manual(name = NULL,
      values = c("Kernel (bootstrap 95%)" = "dashed", "Model (posterior 95%)" = "solid")) +
    guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
    theme(legend.position = "top", legend.justification = "right",
          legend.key.width = grid::unit(1.1, "cm"),
          legend.key.height = grid::unit(0.25, "cm"),
          legend.text = element_text(size = 8),
          legend.spacing = grid::unit(0.4, "lines"))
}

population_predictive_plot <- function(population_predictive, observed_hist, observed_kde) {
  ggplot(population_predictive, aes(theta, center)) +
    geom_ribbon(aes(ymin = low, ymax = high, fill = model), alpha = 0.2) +
    geom_line(aes(colour = model), linewidth = 0.55) +
    geom_errorbar(data = observed_hist, aes(ymin = low, ymax = high), width = 0.05) +
    geom_point(data = observed_hist, size = 1) +
    geom_ribbon(data = observed_kde, aes(ymin = low, ymax = high), alpha = 0.12) +
    geom_line(data = observed_kde, linetype = "dashed", linewidth = 0.5) +
    scale_colour_manual(values = c(Bimodal = ORANGE, Unimodal = BLUE)) +
    scale_fill_manual(values = c(Bimodal = ORANGE, Unimodal = BLUE)) +
    scale_x_continuous(breaks = angle_breaks, labels = angle_labels) +
    facet_grid(estimate ~ model) +
    labs(x = "Heading (radians)", y = "Density",
      caption = "Black: observed bootstrap. Colour: posterior predictive replicates.") +
    theme(legend.position = "none")
}

individual_predictive_plot <- function(individual_predictive, ppc_curves) {
  peak_max_ppc <- max(individual_predictive$center)
  individual_predictive$r_low <- 1 + individual_predictive$low / peak_max_ppc
  individual_predictive$r_center <- 1 + individual_predictive$center / peak_max_ppc
  individual_predictive$r_high <- 1 + individual_predictive$high / peak_max_ppc
  ppc_curves$kde_r <- 1 + ppc_curves$kde / peak_max_ppc
  rmax_ppc <- max(2, 1 + 1.15 * max(individual_predictive$high,
    ppc_curves$kde) / peak_max_ppc)
  ggplot(individual_predictive,
    aes(theta, r_center)) +
    geom_hline(yintercept = c(1, 2), colour = "grey40", linewidth = 0.25) +
    geom_segment(data = data.frame(x = 0, xend = 0, y = 1, yend = rmax_ppc),
      inherit.aes = FALSE,
      aes(x = x, xend = xend, y = y, yend = yend),
      colour = "grey20", linewidth = 0.25) +
    geom_ribbon(aes(ymin = r_low, ymax = r_high, fill = model), alpha = 0.14) +
    geom_line(aes(colour = model), linewidth = 0.4) +
    geom_line(data = ppc_curves, aes(y = kde_r),
      linetype = "dashed", linewidth = 0.4) +
    scale_colour_manual(values = c(Bimodal = ORANGE, Unimodal = BLUE)) +
    scale_fill_manual(values = c(Bimodal = ORANGE, Unimodal = BLUE)) +
    scale_x_continuous(limits = c(-pi, pi), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, rmax_ppc), expand = c(0, 0)) +
    coord_polar(start = pi / 2, direction = -1) +
    facet_wrap(~ individual, ncol = 7) + theme_void(base_size = 8) +
    labs(colour = NULL, fill = NULL, caption = "Dashed: observed bootstrap KDE.") +
    theme(legend.position = "bottom")
}

population_dots <- function(table, label, baseline, invert = FALSE,
                          baseline_label = NULL) {
  levels <- if (invert) c("Optic flow", "Vec. encoding")
          else c("Vec. encoding", "Optic flow")
  p <- ggplot(table, aes(center, factor(hypothesis, levels = levels),
          colour = hypothesis, shape = hypothesis)) +
    geom_errorbar(aes(xmin = low, xmax = high), orientation = "y", width = 0.15) +
    geom_point(size = 2.5) + geom_vline(xintercept = baseline, linetype = "dashed") +
    scale_colour_manual(values = HYP_COLORS) +
    scale_shape_manual(values = c("Optic flow" = 16, "Vec. encoding" = 15)) +
    labs(x = label, y = "") + theme(legend.position = "none")
  if (!is.null(baseline_label)) {
    # Python rotates a label next to the dashed baseline; vjust = 1 puts
    # the text body to the left of the anchor once angle = 90.
    p <- p +
      annotate("text",
               x = baseline - diff(range(table$center)) * 0.015, y = 1.5,
               label = baseline_label, angle = 90, hjust = 0.5, vjust = 1,
               size = 2.8)
  }
  p
}

individual_dots <- function(table, label, baseline, limits = NULL) {
  labels <- unique(table$individual)
  table$individual <- factor(table$individual, levels = rev(labels))
  dodge <- position_dodge(width = 0.5)
  main <- ggplot(table, aes(center, individual, colour = hypothesis, shape = hypothesis)) +
    geom_errorbar(aes(xmin = low, xmax = high), orientation = "y", width = 0.2, position = dodge) +
    geom_point(size = 1.6, position = dodge) +
    geom_vline(xintercept = baseline, linetype = "dashed") +
    scale_colour_manual(values = HYP_COLORS) +
    scale_shape_manual(values = c("Optic flow" = 16, "Vec. encoding" = 15)) +
    labs(x = label, y = "Individual", tag = "a.") +
    theme(legend.position = "top",
          legend.background = element_rect(fill = "white", colour = "grey50"))
  if (!is.null(limits)) main <- main + coord_cartesian(xlim = limits)
  probs <- table[!duplicated(table$individual), ]
  probability <- ggplot(probs, aes(probability, individual)) +
    geom_point(size = 1.6) + geom_vline(xintercept = 0.5, linetype = "dashed") +
    geom_vline(xintercept = 1) + scale_x_continuous(limits = c(0, 1.02)) +
    labs(x = "P(Vec. encoding > Optic flow)", y = "", tag = "b.") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  main + probability + plot_layout(widths = c(1.35, 1))
}
