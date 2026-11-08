# Beetles (security light): figure construction.
#
# Sourced by notebooks_r/beetles.ipynb (no standalone output). R twin of
# notebooks/beetles_plots.py: the notebook keeps the fit, the diagnostic
# values, the KDE grids and the bootstrap tables; only the plot
# construction lives here.

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

# Per-individual raw-heading polar grid: each condition's points with the
# circular mean (solid spoke), circular median (dashed spoke) and a band of
# mean +/- one normalized circular standard deviation per beetle panel
# (twin of beetles_plots.py's plot_data_polar).
data_polar_plot <- function(raw_plot, data, radial_max, condition_colors) {
  circ_sd_norm <- function(a) {
    as.numeric(circular::sd.circular(vmu(a), high = pi, low = -pi,
      normalize = TRUE))[1]
  }
  blocks <- do.call(rbind, lapply(split(data, list(data$Individual,
    data$condition), drop = TRUE), function(g) {
    a <- c(circ_mean(g$HeadingRadians),
      as.numeric(circular::median.circular(vmu(g$HeadingRadians)))[1])
    data.frame(Individual = g$Individual[1], condition = g$condition[1],
      theta = a %% (2 * pi), radius = 0, radius_end = radial_max,
      wedge_center = a[1], wedge_sd = circ_sd_norm(g$HeadingRadians),
      statistic = c("circular mean", "circular median"))
  }))
  # Wedge sector mean +/- sd, split into at most two polygons across the
  # (0, 2 pi) seam. Each polygon carries its Individual so every panel only
  # draws the two wedges of its own beetle; without that column ggplot paints
  # all twenty wedges into every panel.
  wedges <- do.call(rbind, lapply(which(blocks$statistic == "circular mean"), function(k) {
    lo <- (blocks$wedge_center[k] - blocks$wedge_sd[k]) %% (2 * pi)
    hi <- (blocks$wedge_center[k] + blocks$wedge_sd[k]) %% (2 * pi)
    chunks <- if (lo < hi) {
      list(seq(lo, hi, length.out = 100L))
    } else {
      list(seq(lo, 2 * pi, length.out = 50L), c(0, seq(0, hi, length.out = 49L)))
    }
    do.call(rbind, lapply(seq_along(chunks), function(j) {
      th <- chunks[[j]]
      data.frame(Individual = blocks$Individual[k],
        theta = c(th[1L], th, th[length(th)]),
        radius = c(0, rep(radial_max, length(th)), 0),
        condition = blocks$condition[k], group = 10L * k + j)
    }))
  }))
  ring_labels <- data.frame(theta = pi / 8,
    radius = seq(2, radial_max, by = 2L),
    lab = as.character(seq(2, radial_max, by = 2L)))
  p_data <- ggplot(raw_plot, aes(theta, radius, colour = condition)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_polygon(data = wedges, aes(theta, radius, group = group,
      fill = condition), inherit.aes = FALSE, alpha = 0.2, colour = NA) +
    geom_segment(data = blocks,
      aes(x = theta, xend = theta, y = radius, yend = radius_end,
          colour = condition, linetype = statistic),
      inherit.aes = FALSE, linewidth = 0.8) +
    scale_x_continuous(limits = c(0, 2 * pi),
      breaks = (0:7) * pi / 4,
      labels = c("0°", "45°", "90°", "135°",
        "180°", "225°", "270°", "315°"),
      expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, radial_max + 1),
      breaks = seq(2, radial_max, by = 2), expand = c(0, 0)) +
    geom_text(data = ring_labels, aes(theta, radius, label = lab),
      inherit.aes = FALSE, size = 3.2) +
    coord_polar() +
    facet_wrap(~ Individual, ncol = 2,
      labeller = as_labeller(setNames(
        paste0("Individual ", seq_len(nlevels(data$Individual))),
        levels(data$Individual)))) +
    scale_colour_manual(values = condition_colors, name = NULL) +
    scale_fill_manual(values = condition_colors, name = NULL) +
    labs(x = "", y = "") +
    theme(legend.position = "right",
      axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p_data
}

# Polar population-mu KDE: radius 0.5 + 0.5*kde/max_kde with the shared
# peak at 1.0, the outer reference ring at 1.11, and the marked
# +/-45-degree floodlight sector arc with end-and-centre ticks.
polar_mu_kde_plot <- function(mu_kde, theta, r_pts, flood_arc, flood_ticks,
  ind_means, condition_colors) {
  p_mu <- ggplot(mu_kde, aes(theta, r, colour = condition)) +
    geom_ribbon(aes(ymin = base, ymax = r, fill = condition),
      alpha = 0.22, colour = NA) +
    geom_line(linewidth = 0.7) +
    geom_path(data = data.frame(theta = theta %% (2 * pi), r = r_pts),
      inherit.aes = FALSE, aes(theta, r), colour = "#BFBFBF", linewidth = 0.32) +
    geom_line(data = flood_arc, inherit.aes = FALSE, aes(theta, r),
      colour = "black", linewidth = 0.4, alpha = 0.75) +
    geom_segment(data = flood_ticks, inherit.aes = FALSE,
      aes(x = theta, xend = theta, y = lo, yend = hi),
      colour = "black", linewidth = 0.4, alpha = 0.75) +
    geom_point(data = ind_means, inherit.aes = FALSE,
      aes(theta, r, colour = condition), size = 3.5, alpha = 0.95) +
    scale_x_continuous(limits = c(0, 2 * pi), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, r_pts * 1.1), expand = c(0, 0)) +
    coord_polar() +
    scale_colour_manual(values = condition_colors) +
    scale_fill_manual(values = condition_colors) +
    labs(x = "", y = "") +
    theme(legend.position = "none", axis.text = element_blank(),
      axis.ticks = element_blank(), panel.grid = element_blank(),
      panel.border = element_blank())
  p_mu
}

# Python-twin treatment figure: one row of four point-range panels -
# within-individual concentration, between-individual concentration,
# floodlight attraction (dashed Uniform reference at the 25% sector mass)
# and the ON-minus-OFF change in that mass. Vertical error bars per
# condition (Off grey / On gold), mean marker + 89% HDI with caps, value on
# a 0-1 y axis.
treatment_hdi_row <- function(population_kappa, group_alignment, flood_mass,
  condition_colors) {
  off_idx <- 1L
  on_idx <- 2L
  cond2_colors <- c(Off = condition_colors[["light OFF"]],
    On = condition_colors[["light ON"]])
  add_hdi <- function(samples, cond, pos) {
    iv <- interval_row(samples, prob = 0.89)
    data.frame(pos = pos, cond = cond, value = iv[["mean"]],
      low = iv[["low"]], high = iv[["high"]], kind = "post")
  }
  # Per-row fills keyed on the actual row keys, so single-condition
  # panels (e.g. the "On - Off" treatment effect) work with a 1-entry cols.
  hdi_fill_values <- function(rows, cols) {
    vals <- vapply(rows$key, function(k) {
      if (startsWith(k, "boot_")) "white" else cols[[sub("^post_", "", k)]]
    }, character(1))
    names(vals) <- rows$key
    vals
  }
  hdi_bar <- function(rows, title, ylab, cols, xlims) {
    rows$key <- paste(rows$kind, rows$cond, sep = "_")
    p <- ggplot(rows, aes(x = pos, y = value, colour = cond, fill = key)) +
      geom_errorbar(aes(ymin = low, ymax = high), width = 0.16,
        linewidth = 0.3) +
      geom_point(shape = 16, size = 1.5) +
      scale_x_continuous(limits = xlims, breaks = unique(rows$pos),
        labels = rows$cond[order(rows$pos)], expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0),
        breaks = seq(0, 1, by = 0.2),
        labels = function(v) format(v, digits = 1)) +
      scale_colour_manual(values = cols, guide = "none") +
      scale_fill_manual(values = hdi_fill_values(rows, cols), guide = "none") +
      labs(x = "", y = ylab, title = title) +
      theme_classic(base_size = 10) +
      theme(plot.title = element_text(hjust = 0.5, size = 11))
    p
  }
  within_df <- rbind(add_hdi(A1_safe(population_kappa[, off_idx]), "Off", 1L),
                     add_hdi(A1_safe(population_kappa[, on_idx]), "On", 2L))
  between_df <- rbind(add_hdi(group_alignment[, off_idx], "Off", 1L),
                      add_hdi(group_alignment[, on_idx], "On", 2L))
  mass_df <- rbind(add_hdi(flood_mass[, off_idx], "Off", 1L),
                   add_hdi(flood_mass[, on_idx], "On", 2L))
  delta_df <- add_hdi(flood_mass[, on_idx] - flood_mass[, off_idx],
    "On - Off", 1L)
  xlim2 <- c(0.75, 2.25)
  xlim1 <- c(0.25, 1.75)
  p_conc_within <- hdi_bar(within_df, "Within-individual\nconcentration",
    "Mean resultant length", cond2_colors, xlim2)
  p_conc_between <- hdi_bar(between_df, "Between-individual\nconcentration",
    "Mean resultant length", cond2_colors, xlim2)
  p_conc_attract <- hdi_bar(mass_df, "Floodlight\nattraction",
    "Predictive mass (+-45 deg)", cond2_colors, xlim2) +
    geom_hline(yintercept = 0.25, linetype = "dashed", linewidth = 0.35) +
    annotate("text", x = 2.22, y = 0.22, label = "Uniform",
      hjust = 1, vjust = 1, size = 3.3, colour = "black")
  p_conc_delta <- hdi_bar(delta_df, "Floodlight\ntreatment effect",
    "Difference in\npredictive mass (+-45 deg)",
    c("On - Off" = cond2_colors[["On"]]), xlim1)
  p_concentration <- patchwork::wrap_plots(p_conc_within, p_conc_between,
  p_conc_attract, p_conc_delta, ncol = 4)
  p_concentration
}

# Density of the posterior ON-minus-OFF change in floodlight-sector mass.
mass_delta_plot <- function(mass_delta) {
  p_mass_delta <- ggplot(mass_delta, aes(delta)) + geom_density(fill = "#D9A521", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(x = "change in floodlight-sector probability (ON - OFF)", y = "density")
  p_mass_delta
}

# Python-twin supplemental PPC figure: group mean direction (a) and group
# alignment (b) by condition on top, per-individual mean direction
# (c = OFF, d = ON) and per-individual concentration (e = OFF, f = ON)
# below. Bootstrap = open marker, posterior = filled; light OFF grey,
# light ON gold; circular intervals wrapped to 0-360 degrees with
# minimal-arc error widths.
ppc_twin_plot <- function(J, ids, individual_summary, boot_individual,
  boot_group_mu, boot_group_R, population_mu, group_alignment,
  condition_labels, condition_colors) {
    off_idx <- 1L
  on_idx <- 2L
  nT <- J
to360 <- function(rad) (rad * 180 / pi + 360) %% 360
  dir_iv <- function(samples) interval_row(samples, angular = TRUE)
  lin_iv <- function(samples) interval_row(samples)
  ppc_cond_colors <- c(Off = condition_colors[["light OFF"]],
    On = condition_colors[["light ON"]])
  ppc_fill_vals <- c(boot_Off = "white", post_Off = ppc_cond_colors[["Off"]],
    boot_On = "white", post_On = ppc_cond_colors[["On"]])
  mu_row <- function(iv, cond, y, kind) {
    m <- iv["mean"]; lo <- iv["low"]; hi <- iv["high"]
    w <- c(abs(mod_circular(m - lo)), abs(mod_circular(hi - m))) * 180 / pi
    x <- to360(m)
    data.frame(cond = cond, y = y, x = x, xlo = x - w[1], xhi = x + w[2],
      kind = kind)
  }
  lin_row <- function(iv, cond, y, kind) {
    data.frame(cond = cond, y = y, x = iv["mean"], xlo = iv["low"],
      xhi = iv["high"], kind = kind)
  }
  add_key <- function(d) {
    d$key <- paste(d$kind, d$cond, sep = "_")
    d
  }
  df_gmu <- rbind(
    mu_row(dir_iv(boot_group_mu[, off_idx]), "Off", -0.12, "boot"),
    mu_row(dir_iv(population_mu[, off_idx]), "Off", 0.12, "post"),
    mu_row(dir_iv(boot_group_mu[, on_idx]), "On", 1 - 0.12, "boot"),
    mu_row(dir_iv(population_mu[, on_idx]), "On", 1 + 0.12, "post"))
  df_gr <- rbind(
    lin_row(lin_iv(boot_group_R[, off_idx]), "Off", -0.12, "boot"),
    lin_row(lin_iv(group_alignment[, off_idx]), "Off", 0.12, "post"),
    lin_row(lin_iv(boot_group_R[, on_idx]), "On", 1 - 0.12, "boot"),
    lin_row(lin_iv(group_alignment[, on_idx]), "On", 1 + 0.12, "post"))
  top_y <- scale_y_continuous(limits = c(-0.35, 1.35), expand = c(0, 0),
    breaks = c(0, 1), labels = c("OFF", "ON"))
  panel_label <- function(txt, x, y) annotate("text", x = x, y = y,
    label = txt, size = 3.6, fontface = "bold", hjust = 0, vjust = 1,
    colour = "black")
  p_g_mu <- ggplot(add_key(df_gmu), aes(x = x, y = y)) +
    geom_segment(aes(xend = xlo, yend = y, colour = cond), linewidth = 0.3) +
    geom_segment(aes(xend = xhi, yend = y, colour = cond), linewidth = 0.3) +
    geom_point(aes(shape = cond, colour = cond, fill = key), size = 1.5) +
    scale_colour_manual(values = ppc_cond_colors, guide = "none") +
    scale_fill_manual(values = ppc_fill_vals, guide = "none") +
    scale_shape_manual(values = c(Off = 22, On = 21), guide = "none") +
    scale_x_continuous(limits = c(0, 360), breaks = c(0, 90, 180, 270, 360)) +
    top_y +
    labs(x = "Group mean direction (deg)", y = "",
      title = "Group mean direction") +
    panel_label("a.", 4, 1.32) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 11.5),
      legend.position = "none", clip = "on")
  p_gr <- ggplot(add_key(df_gr), aes(x = x, y = y)) +
    geom_segment(aes(xend = xlo, yend = y, colour = cond), linewidth = 0.3) +
    geom_segment(aes(xend = xhi, yend = y, colour = cond), linewidth = 0.3) +
    geom_point(aes(shape = cond, colour = cond, fill = key), size = 1.5) +
    scale_colour_manual(values = ppc_cond_colors, guide = "none") +
    scale_fill_manual(values = ppc_fill_vals, guide = "none") +
    scale_shape_manual(values = c(Off = 22, On = 21), guide = "none") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    top_y +
    labs(x = "Mean resultant length", y = "", title = "Group alignment") +
    panel_label("b.", 0.02, 1.32) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 11.5),
      legend.position = "none", clip = "on")
  # Per-individual rows: open bootstrap (y - 0.14) under filled posterior
  # (y + 0.14) for every beetle in every condition.
  ind_mu_off <- ind_mu_on <- NULL
  ind_r_off <- ind_r_on <- NULL
  po_pick <- function(id, cond, qty) {
    x <- individual_summary[individual_summary$id == id &
      individual_summary$condition == cond & individual_summary$quantity == qty, ]
    stopifnot(nrow(x) == 1L)
    c(mean = x$mean[1], low = x$low[1], high = x$high[1])
  }
  for (t in seq_len(J)) {
    id <- as.character(ids[t])
    ind_mu_off <- rbind(ind_mu_off,
      mu_row(dir_iv(boot_individual[[t]]$mu), "Off", t - 0.14, "boot"),
      mu_row(po_pick(id, condition_labels[1], "mean direction (rad)"),
        "Off", t + 0.14, "post"))
    ind_r_off <- rbind(ind_r_off,
      lin_row(lin_iv(boot_individual[[t]]$R), "Off", t - 0.14, "boot"),
      lin_row(po_pick(id, condition_labels[1], "within-individual R(kappa)"),
        "Off", t + 0.14, "post"))
    ind_mu_on <- rbind(ind_mu_on,
      mu_row(dir_iv(boot_individual[[J + t]]$mu), "On", t - 0.14, "boot"),
      mu_row(po_pick(id, condition_labels[2], "mean direction (rad)"),
        "On", t + 0.14, "post"))
    ind_r_on <- rbind(ind_r_on,
      lin_row(lin_iv(boot_individual[[J + t]]$R), "On", t - 0.14, "boot"),
      lin_row(po_pick(id, condition_labels[2], "within-individual R(kappa)"),
        "On", t + 0.14, "post"))
  }
  ind_y <- function(show_labels) scale_y_reverse(
    limits = c(0.25, nT + 0.75), expand = c(0, 0),
    breaks = seq_len(nT),
    labels = if (show_labels) as.character(ids) else NULL)
  ind_grid <- geom_hline(yintercept = seq(0.5, nT - 0.5, by = 1),
    colour = "#e6e6e6", linewidth = 0.15)
  ind_mu_panel <- function(rows, lab, show_labels, xlab) {
    ggplot(add_key(rows), aes(x = x, y = y)) +
      ind_grid +
      geom_segment(aes(xend = xlo, yend = y, colour = cond), linewidth = 0.25) +
      geom_segment(aes(xend = xhi, yend = y, colour = cond), linewidth = 0.25) +
      geom_point(aes(shape = cond, colour = cond, fill = key), size = 1.9) +
      scale_colour_manual(values = ppc_cond_colors, guide = "none") +
      scale_fill_manual(values = ppc_fill_vals, guide = "none") +
      scale_shape_manual(values = c(Off = 22, On = 21), guide = "none") +
      scale_x_continuous(limits = c(0, 360), breaks = c(0, 180, 360)) +
      ind_y(show_labels) +
      labs(x = xlab, y = "") +
      panel_label(lab, 3, 0.3) +
      theme_classic(base_size = 9) +
      theme(legend.position = "none", clip = "on")
  }
  ind_r_panel <- function(rows, lab, show_labels, xlab) {
    ggplot(add_key(rows), aes(x = x, y = y)) +
      ind_grid +
      geom_segment(aes(xend = xlo, yend = y, colour = cond), linewidth = 0.25) +
      geom_segment(aes(xend = xhi, yend = y, colour = cond), linewidth = 0.25) +
      geom_point(aes(shape = cond, colour = cond, fill = key), size = 1.9) +
      scale_colour_manual(values = ppc_cond_colors, guide = "none") +
      scale_fill_manual(values = ppc_fill_vals, guide = "none") +
      scale_shape_manual(values = c(Off = 22, On = 21), guide = "none") +
      scale_x_continuous(expand = c(0.05, 0.05)) +
      ind_y(show_labels) +
      labs(x = xlab, y = "") +
      panel_label(lab, 0.55, 0.3) +
      theme_classic(base_size = 9) +
      theme(legend.position = "none", clip = "on")
  }
  p_c <- ind_mu_panel(ind_mu_off, "c.", TRUE, "")
  p_d <- ind_mu_panel(ind_mu_on, "d.", FALSE, "Individual mean direction (deg)")
  p_e <- ind_r_panel(ind_r_off, "e.", FALSE, "")
  p_f <- ind_r_panel(ind_r_on, "f.", FALSE, "Individual concentration (R)")
  # Manual figure legend row (the Python fig.legend items).
  leg_items <- data.frame(
    x = c(0.7, 2.9, 5.0, 7.0),
    label = c("Bootstrap (data)", "Posterior (model)", "Light OFF", "Light ON"))
  p_leg <- ggplot(leg_items, aes(x = x, y = 1)) +
    geom_point(shape = c(21, 21, 15, 16),
      colour = c("black", "black", ppc_cond_colors[["Off"]],
        ppc_cond_colors[["On"]]),
      fill = c("white", "black", ppc_cond_colors[["Off"]],
        ppc_cond_colors[["On"]]),
      size = 2.2) +
    geom_text(aes(x = x + 0.34, y = 1, label = label),
      hjust = 0, size = 3.6, colour = "black") +
    scale_x_continuous(limits = c(0, 8.4), expand = c(0, 0)) +
    theme_void(base_size = 10)
  ppc_twin <- ((p_g_mu | p_gr) /
      (p_c | p_d | p_e | p_f) / p_leg) +
    patchwork::plot_layout(heights = c(0.55, 1.2, 0.16))
  ppc_twin
}

# Individual posterior summaries: paired mean direction and
# within-individual R(kappa) by individual and condition.
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

# Conventional replicated-data check: replicated mean resultant length
# per individual/condition with the observed R as an open marker.
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
