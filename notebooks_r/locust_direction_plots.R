# Locust direction mixture: figure construction.
#
# Sourced by notebooks_r/locust_direction_mixture_model.ipynb (no
# standalone output). R twin of notebooks/locust_direction_plots.py.

BLUE <- "#1F77B4"
ORANGE <- "#FF7F0E"
RED <- "#D62728"
HYP_COLORS <- c("Optic flow" = BLUE, "Vec. encoding" = ORANGE)
angle_breaks <- c(-pi, -pi / 2, 0, pi / 2, pi)
angle_labels <- c("-pi", "-pi/2", "0", "pi/2", "pi")

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

observed_polar_grid <- function(raw_densities) {
  p_observed <- ggplot(raw_densities) +
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

elpd_comparison_plot <- function(elpd_table) {
  p_elpd <- ggplot(elpd_table, aes(elpd, factor(model, levels = c("Bimodal", "Unimodal")),
                                   colour = model, shape = model)) +
    scale_y_discrete(labels = model_labels) +
    geom_errorbar(aes(xmin = elpd - se, xmax = elpd + se), orientation = "y", width = 0.15) +
    geom_point(size = 2.5) + scale_colour_manual(values = c(Bimodal = ORANGE, Unimodal = BLUE)) +
    scale_shape_manual(values = c(Bimodal = 15, Unimodal = 16)) +
    labs(x = "PSIS-LOO ELPD +/-1 SE", y = "") + theme(legend.position = "none")
}

elpd_individual_plot <- function(elpd_individual, n_ind) {
  individual_levels <- rev(as.character(seq_len(num_individuals)))
  elpd_long <- rbind(data.frame(individual = elpd_individual$individual,
                               value = elpd_individual$bimodal, model = "Bimodal"),
                     data.frame(individual = elpd_individual$individual,
                               value = elpd_individual$unimodal, model = "Unimodal"))
  p_elpd_individual <- ggplot(elpd_long,
    aes(value, factor(individual, levels = individual_levels), colour = model, shape = model)) +
    geom_point(position = position_dodge(width = 0.5), size = 1.8) +
    scale_colour_manual(values = c(Bimodal = ORANGE, Unimodal = BLUE)) +
    scale_shape_manual(values = c(Bimodal = 15, Unimodal = 16)) +
    labs(x = "PSIS-LOO ELPD contribution", y = "Individual") + theme(legend.position = "bottom")
  p_delta <- ggplot(elpd_individual,
    aes(delta, factor(individual, levels = individual_levels))) +
    geom_point(size = 1.8) + geom_vline(xintercept = 0, linetype = "dashed") +
    labs(x = "Delta ELPD (Bimodal - Unimodal)", y = "")
  p_elpd_individual + p_delta
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

# Angle ticks use the twin's degree labels (locust_direction_plots.py):
# the mathematical orientation (0 at 3 o'clock, counter-clockwise) is the
# twin default for this x range, and the density-zero ring at radius 0.5
# is drawn as a black circle, as in the twin.
population_polar_plot <- function(model_density, kde_density, hist_density, edge_radius = 0.5) {
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
  p_validation_polar <- p_validation_polar + plot_layout(widths = c(4, 1))
  p_validation_polar
}

individual_ppc_plot <- function(ppc_curves, ppc_ticks, ppc_arrows) {
  rmax <- max(2, 1 + max(ppc_curves$model_high, ppc_curves$kde_high) * 1.15)
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
  
  p_individual_ppc <- p_individual_ppc +
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

