# Locust density/order: figure construction.
#
# Sourced by notebooks_r/locust_density_order.ipynb (no standalone
# output). R twin of notebooks/locust_density_plots.py: the notebook
# keeps the data, the preprocessing, the fits and their fit-cache
# names, the diagnostics and the post-hoc tables; only the figure
# construction lives here (plus the colour and theme helpers the
# grids need).


cell_color2 <- function(r, c2, n_rows, n_cols, tmin = 0.10, tmax = 0.95,
  blue = c(0.10, 0.25, 1.00), red = c(1.00, 0.15, 0.15), gamma = 1.15) {
  to <- tmin + (tmax - tmin) * (r / max(n_rows - 1, 1)) ^ gamma
  td <- tmin + (tmax - tmin) * (c2 / max(n_cols - 1, 1)) ^ gamma
  # Channelwise max, per-channel (no matrix ops: keeps scalar and
  # vectorised calls warning-free and dimension-stable).
  rr <- pmax(0, pmin(1, pmax(to * red[1], td * blue[1])))
  gg <- pmax(0, pmin(1, pmax(to * red[2], td * blue[2])))
  bb <- pmax(0, pmin(1, pmax(to * red[3], td * blue[3])))
  sprintf("#%02X%02X%02X", round(rr * 255), round(gg * 255), round(bb * 255))
}

cell_title <- function(o, d)
  paste0("Density = ", format(d, nsmall = 1, trim = TRUE), "\nOrder = ", format(o, nsmall = 1))

grid_thm <- function() theme_classic(base_size = 8) +
  theme(strip.text = element_text(size = 8, hjust = 0.5, lineheight = 0.95),
    strip.background = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.28),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = margin(1, 1, 1, 1, "pt"))

nice_ceil <- function(v) {
  if (v <= 0) return(1.0)
  v <- v * 1.10
  e <- floor(log10(v))
  f <- v / 10 ^ e
  for (n in c(1.0, 2.0, 2.5, 5.0, 10.0)) if (f <= n) return(n * 10 ^ e)
  10.0 * 10 ^ e
}

grid_color_key <- function(n_rows, n_cols) {
  key_grid <- expand.grid(rr = 0:(n_rows - 1), cc = 0:(n_cols - 1))
  key_grid$ypos <- n_rows - 1L - key_grid$rr
  key_grid$col <- cell_color2(key_grid$rr, key_grid$cc, n_rows, n_cols)
  ggplot(key_grid, aes(x = cc + 0.5, y = ypos + 0.5, fill = col)) +
    geom_tile() +
    annotate("rect", xmin = 0, xmax = n_cols, ymin = 0, ymax = n_rows,
      fill = NA, colour = "black", linewidth = 0.28) +
    scale_fill_identity(guide = "none") +
    scale_x_continuous(limits = c(-2, n_cols), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-6.5, n_rows + 6.5), expand = c(0, 0)) +
    annotate("text", x = 0.5, y = -0.85, label = "low", size = 2.3) +
    annotate("text", x = n_cols - 0.5, y = -0.85, label = "high", size = 2.3) +
    annotate("text", x = n_cols / 2, y = -1.5, label = "Density", size = 2.3) +
    annotate("text", x = -0.60, y = n_rows - 0.5, label = "low", hjust = 1,
      size = 2.3) +
    annotate("text", x = -0.60, y = 0.5, label = "high", hjust = 1, size = 2.3) +
    annotate("text", x = -1.60, y = n_rows / 2, label = "Order", angle = 90,
      size = 2.3) +
    coord_equal() +
    theme_void(base_size = 7) +
    theme(plot.margin = margin(1.5, 1.5, 1.5, 1.5, "pt"))
}

trajectory_grid_plot <- function(df, n_cols) {
  ggplot(df, aes(X, Y, group = fname, colour = gcol)) +
    geom_path(linewidth = 0.28, alpha = 0.75) +
    annotate("segment", x = 0, xend = 8, y = 0, yend = 0,
      colour = "black", linewidth = 0.35) +
    facet_wrap(~ gtitle, ncol = n_cols, axes = "all") +
    scale_colour_identity(guide = "none") +
    scale_x_continuous(breaks = c(-8, 0, 8), expand = c(0, 0)) +
    scale_y_continuous(breaks = c(-8, 0, 8), expand = c(0, 0)) +
    coord_equal(xlim = c(-8, 8), ylim = c(-8, 8)) +
    labs(x = "", y = "") +
    grid_thm()
}

mean_vector_grid_plot <- function(vectors, circle, zero_seg, zero_tick, n_cols, Lvec) {
  ggplot(vectors, aes(group = fname)) +
    geom_path(data = circle, aes(x, y), inherit.aes = FALSE, colour = "black",
      linewidth = 0.28) +
    geom_segment(aes(x = 0, y = 0, xend = vx, yend = vy, colour = gcol),
      linewidth = 0.39, alpha = 0.8) +
    geom_point(data = data.frame(x = 0, y = 0), aes(x, y), inherit.aes = FALSE,
      size = 0.12, colour = "black") +
    geom_segment(data = zero_seg, aes(x = x, y = y, xend = xend, yend = yend), inherit.aes = FALSE,
      colour = "black", linewidth = 0.32) +
    geom_segment(data = zero_tick, aes(x = x, y = y, xend = xend, yend = yend), inherit.aes = FALSE,
      colour = "black", linewidth = 0.35, lineend = "round") +
    facet_wrap(~ gtitle, ncol = n_cols) +
    scale_colour_identity(guide = "none") +
    scale_x_continuous(limits = c(-Lvec, Lvec), breaks = NULL) +
    scale_y_continuous(limits = c(-Lvec, Lvec), breaks = NULL) +
    coord_equal() +
    labs(x = "", y = "") +
    grid_thm() +
    theme(axis.line = element_blank())
}

stimulus_grid_plot <- function() {
  orders_s <- order_levels
  densities_s <- density_levels
  n_rows_s <- length(orders_s)
  n_cols_s <- length(densities_s)
  n_steps_s <- 60L
  step_L_s <- 0.04
  cell_color_s <- function(r, c2) {
    to <- 0.10 + 0.85 * (r / (n_rows_s - 1)) ^ 1.15
    td <- 0.10 + 0.85 * (c2 / (n_cols_s - 1)) ^ 1.15
    rgb <- pmax(0, pmin(1, pmax(to * c(1.00, 0.15, 0.15),
      td * c(0.10, 0.25, 1.00))))
    sprintf("#%02X%02X%02X", round(rgb[1] * 255), round(rgb[2] * 255),
      round(rgb[3] * 255))
  }
  walks <- vector("list", n_rows_s * n_cols_s)
  wi <- 0L
  for (r in 0:(n_rows_s - 1)) {
    for (c2 in 0:(n_cols_s - 1)) {
      o <- orders_s[r + 1L]
      n_traj <- densities_s[c2 + 1L]
      wi <- wi + 1L
      set.seed(123 + 1000L * r + 10L * c2)
      if (o == 0) {
        x0 <- runif(n_traj, 0.10, 0.90)
        y0 <- runif(n_traj, 0.10, 0.90)
      } else {
        x0 <- rep(0, n_traj)
        y0 <- if (n_traj == 1L) 0.5 else seq(0.05, 0.95, length.out = n_traj)
      }
      k <- inverse_A1(o)
      pts <- matrix(0, nrow = 0, ncol = 2)
      for (i in seq_len(n_traj)) {
        if (o == 1) {
          angles <- rep(0, n_steps_s)
        } else {
          angles <- mod_circular(as.numeric(circular::rvonmises(n_steps_s,
            mu = vmu(0.0), kappa = k, control.circular = vm_ctrl)))
        }
        pts <- rbind(pts, cbind(
          c(x0[i], x0[i] + cumsum(step_L_s * cos(angles))),
          c(y0[i], y0[i] + cumsum(step_L_s * sin(angles)))))
      }
      walks[[wi]] <- data.frame(o = r + 1L, d = c2 + 1L,
        i = rep(seq_len(n_traj), each = n_steps_s + 1L),
        x = pts[, 1], y = pts[, 2])
    }
  }
  schematic <- do.call(rbind, walks)
  schematic$o_lab <- factor(orders_s[schematic$o],
    levels = orders_s, labels = format(orders_s, nsmall = 1))
  schematic$d_lab <- factor(densities_s[schematic$d],
    levels = densities_s, labels = format(densities_s))
  schematic$cell <- mapply(cell_color_s, schematic$o - 1L, schematic$d - 1L)
  schematic$lw <- ifelse(densities_s[schematic$d] >= 16, 0.11, 0.19)

  border_thm <- element_rect(colour = "black", linewidth = 0.28, fill = NA)
  p_panels <- list()
  for (r in 0:(n_rows_s - 1)) {
    for (c2 in 0:(n_cols_s - 1)) {
      sub <- schematic[schematic$o == r + 1L & schematic$d == c2 + 1L, ,
        drop = FALSE]
      p <- ggplot(sub[densities_s[sub$d] < 16, , drop = FALSE],
        aes(x, y, group = i, colour = cell)) +
        geom_path(data = sub[densities_s[sub$d] >= 16, , drop = FALSE],
          aes(x, y, group = i, colour = cell),
          linewidth = 0.11, alpha = 0.85) +
        geom_path(linewidth = 0.19, alpha = 0.85) +
        scale_colour_identity(guide = "none") +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
        theme_void() +
        theme(panel.border = border_thm,
          plot.margin = margin(1.5, 1.5, 1.5, 1.5, "pt"),
          plot.title = element_text(size = 8, hjust = 0.5))
      if (r == 0L) {
        p <- p + labs(title = format(densities_s[c2 + 1L]))
      }
      p_panels[[length(p_panels) + 1L]] <- p
    }
  }
  grid_block <- patchwork::wrap_plots(p_panels, ncol = n_cols_s)
  gutter_df <- data.frame(lbl = format(orders_s, nsmall = 1),
    y = n_rows_s:1L)
  p_gutter <- ggplot(gutter_df) +
    geom_text(aes(x = 0.98, y = y, label = lbl), hjust = 1, size = 2.6) +
    annotate("text", x = 0.30, y = (n_rows_s + 1) / 2,
      label = "Order", angle = 90, size = 2.6) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, n_rows_s + 0.5), expand = c(0, 0)) +
    theme_void() +
    labs(title = "") +
    theme(plot.margin = margin(1.5, 2, 1.5, 2, "pt"))
  key_grid <- expand.grid(rr = 0:(n_rows_s - 1), cc = 0:(n_cols_s - 1))
  key_grid$ypos <- n_rows_s - 1L - key_grid$rr
  key_grid$col <- mapply(cell_color_s, key_grid$rr, key_grid$cc)
  p_key <- ggplot(key_grid, aes(x = cc + 0.5, y = ypos + 0.5, fill = col)) +
    geom_tile() +
    annotate("rect", xmin = 0, xmax = n_cols_s, ymin = 0, ymax = n_rows_s,
      fill = NA, colour = "black", linewidth = 0.28) +
    scale_fill_identity(guide = "none") +
    scale_x_continuous(limits = c(-2, n_cols_s), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-6.5, n_rows_s + 6.5),
      expand = c(0, 0)) +
    annotate("text", x = 0.5, y = -0.75, label = "low", size = 2.3) +
    annotate("text", x = n_cols_s - 0.5, y = -0.75, label = "high",
      size = 2.3) +
    annotate("text", x = n_cols_s / 2, y = -1.35, label = "Density",
      size = 2.3) +
    annotate("text", x = -0.60, y = n_rows_s - 0.5, label = "low",
      hjust = 1, size = 2.3) +
    annotate("text", x = -0.60, y = 0.5, label = "high", hjust = 1,
      size = 2.3) +
    annotate("text", x = -1.75, y = n_rows_s / 2, label = "Order",
      angle = 90, size = 2.3) +
    coord_equal() +
    labs(title = "") +
    theme_void(base_size = 7) +
    theme(plot.margin = margin(1.5, 1.5, 1.5, 1.5, "pt"))
  p_title <- ggplot(data.frame(), aes(x = 0, y = 0)) +
    annotate("text", x = 0.5, y = 0.5, label = "Density (m\u207b\u00b2)",
      size = 3.2) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void()
  patchwork::wrap_plots(
    p_title,
    patchwork::wrap_plots(p_gutter, grid_block, p_key,
      widths = c(0.4, 6, 1)),
    ncol = 1, heights = c(0.07, 1))
}

schematic_plot <- function() {
  sc_orders <- c(1.0, 0.8, 0.6, 0.4, 0.2, 0.0)
  sc_denss <- c(64, 32, 16, 8, 4, 1)
  sc_reds <- c("64" = "#7E0610", "32" = "#AC1117", "16" = "#CA181D",
    "8" = "#E83429", "4" = "#F6583E", "1" = "#FB7A5A")
  sc_purps <- c("1.0" = "#470F84", "0.8" = "#582F93", "0.6" = "#6950A3",
    "0.4" = "#7B74B5", "0.2" = "#928FC3", "0.0" = "#AAA8D0")
  sc_A1 <- function(k)
    besselI(k, 1, expon.scaled = TRUE) / besselI(k, 0, expon.scaled = TRUE)
  sc_kappa <- function(o) {
    if (o <= 0) return(0)
    if (o >= 1) return(Inf)
    uniroot(function(k) sc_A1(k) - o, c(1e-8, 500))$root
  }
  sc_n_traj <- 12L
  sc_n_steps <- 100L
  sc_step <- 0.04
  sc_box <- theme_void() +
    theme(panel.border = element_rect(colour = "black", fill = NA,
      linewidth = 0.28),
      plot.margin = margin(1, 1, 1, 1, "pt"))
  p_dots <- vector("list", 6)
  p_walks <- vector("list", 6)
  for (r in seq_len(6)) {
    d <- sc_denss[r]
    if (d == 1L) {
      xs <- 0.5; ys <- 0.5
    } else if (d == 4L) {
      g <- expand.grid(c(0.25, 0.75), c(0.25, 0.75))
      xs <- g$Var1; ys <- g$Var2
    } else {
      ns <- ceiling(sqrt(d))
      axis_seq <- seq(0.1, 0.9, length.out = ns)
      g <- expand.grid(axis_seq, axis_seq)
      xs <- g$Var1[seq_len(d)]
      ys <- g$Var2[seq_len(d)]
    }
    xs <- pmin(0.925, pmax(0.075, xs + rnorm(length(xs), 0, 0.025)))
    ys <- pmin(0.925, pmax(0.075, ys + rnorm(length(ys), 0, 0.025)))
    p_dots[[r]] <- ggplot(data.frame(x = xs, y = ys), aes(x, y)) +
      geom_point(size = 0.06, colour = sc_reds[as.character(d)]) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1)) + sc_box
    o <- sc_orders[r]
    k <- sc_kappa(o)
    if (o == 0) {
      x0 <- runif(sc_n_traj, 0.1, 0.9)
      y0 <- runif(sc_n_traj, 0.1, 0.9)
    } else {
      x0 <- rep(0.05, sc_n_traj)
      y0 <- seq(0.15, 0.85, length.out = sc_n_traj)
    }
    pts <- data.frame()
    for (i in seq_len(sc_n_traj)) {
      angles <- if (o == 1) rep(0, sc_n_steps)
        else as.numeric(circular::rvonmises(sc_n_steps, mu = vmu(0.0),
          kappa = k, control.circular = vm_ctrl))
      pts <- rbind(pts, data.frame(
        i = i,
        x = c(x0[i], x0[i] + cumsum(sc_step * cos(angles))),
        y = c(y0[i], y0[i] + cumsum(sc_step * sin(angles)))))
    }
    p_walks[[r]] <- ggplot(pts, aes(x, y, group = i)) +
      geom_path(linewidth = 0.25, colour = sc_purps[sprintf("%.1f", o)]) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1)) + sc_box
  }
  gutters_l <- data.frame(y = rev(seq_len(6)), lbl = as.character(sc_denss))
  gutters_r <- data.frame(y = rev(seq_len(6)),
    lbl = format(sc_orders, nsmall = 1))
  p_gut_l <- ggplot(gutters_l) +
    geom_text(aes(x = 0.9, y = y, label = lbl), hjust = 1, size = 2.4) +
    scale_y_continuous(limits = c(0.5, 6.5), expand = c(0, 0)) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void()
  p_gut_r <- ggplot(gutters_r) +
    geom_text(aes(x = 0.1, y = y, label = lbl), hjust = 0, size = 2.4) +
    scale_y_continuous(limits = c(0.5, 6.5), expand = c(0, 0)) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void()
  patchwork::wrap_plots(
    patchwork::wrap_plots(p_gut_l,
      patchwork::wrap_plots(p_dots, ncol = 1),
      patchwork::wrap_plots(p_walks, ncol = 1), p_gut_r,
      widths = c(0.30, 1, 1, 0.38)),
    ncol = 1)
}

wrap_label <- function(s, width = 40)
  vapply(s, function(x) paste(strwrap(x, width = width), collapse = "\n"),
    character(1))

selection_forest_plot <- function(sweep_df) {
  pretty_names <- sweep_df$model
  pretty_names[pretty_names == "Density x Order"] <- "Density \u00d7 Order"
  group_name <- ifelse(grepl("Order", sweep_df$model),
    "Order included", "Order excluded")
  ypos <- nrow(sweep_df):1L
  # Explicit limits so the annotated legend lands in stable data coordinates on both axes.
  data_lo <- min(sweep_df$delta_elpd - sweep_df$se_delta_elpd)
  data_hi <- max(sweep_df$delta_elpd + sweep_df$se_delta_elpd)
  x_lo_sel <- data_lo - 0.08 * (data_hi - data_lo)
  x_hi_sel <- data_hi + 0.08 * (data_hi - data_lo)
  span_sel <- x_hi_sel - x_lo_sel
  ggplot(sweep_df,
    aes(x = delta_elpd, y = ypos, colour = group_name, fill = group_name)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35) +
    geom_errorbar(aes(xmin = delta_elpd - se_delta_elpd,
      xmax = delta_elpd + se_delta_elpd), orientation = "y",
      width = 0.24, linewidth = 0.35) +
    geom_point(aes(shape = group_name), size = 1.7) +
    scale_colour_manual(values = c("Order included" = "#1F77B4",
      "Order excluded" = "#FF7F0E"), guide = "none") +
    scale_fill_manual(values = c("Order included" = "#1F77B4",
      "Order excluded" = "#FF7F0E"), guide = "none") +
    scale_shape_manual(values = c("Order included" = 16,
      "Order excluded" = 15), guide = "none") +
    scale_x_continuous(limits = c(x_lo_sel, x_hi_sel)) +
    scale_y_continuous(limits = c(0.45, nrow(sweep_df) + 0.55),
      breaks = nrow(sweep_df):1L, labels = wrap_label(pretty_names),
      expand = c(0, 0)) +
    labs(x = "ELPD difference", y = "") +
    # Framed upper-left legend (drawn by hand: this ggplot2 build
    # drops legend.position = "upperleft" without an error).
    annotate("rect", xmin = x_lo_sel + 0.02 * span_sel,
      xmax = x_lo_sel + 0.88 * span_sel, ymin = 4.55, ymax = 5.52,
      fill = "white", colour = "black", linewidth = 0.3) +
    annotate("point", x = x_lo_sel + 0.055 * span_sel, y = 5.24,
      shape = 16, colour = "#1F77B4", size = 2.1) +
    annotate("text", x = x_lo_sel + 0.135 * span_sel, y = 5.24,
      label = "Order included", hjust = 0, vjust = 0.5, size = 2.2) +
    annotate("point", x = x_lo_sel + 0.055 * span_sel, y = 4.78,
      shape = 15, colour = "#FF7F0E", size = 2.1) +
    annotate("text", x = x_lo_sel + 0.135 * span_sel, y = 4.78,
      label = "Order excluded", hjust = 0, vjust = 0.5, size = 2.2) +
    theme_classic(base_size = 8)
}

trace_plot <- function(fit, parameters, angular = character()) {
  # as_draws_df is wide: explicitly assemble a long trace table.
  draws_df <- posterior::as_draws_df(posterior::subset_draws(
      posterior::as_draws_array(fit), variable = trace_parameters))
  trace_df <- do.call(rbind, lapply(trace_parameters, function(p) {
    value <- as.numeric(draws_df[[p]])
    if (p %in% angular) value <- mod_circular(value - circ_mean(value))
    data.frame(parameter = p, value = value,
               iteration = draws_df[[".iteration"]],
               chain = factor(draws_df[[".chain"]]))
  }))
  p_trace <- ggplot(trace_df, aes(iteration, value, colour = chain)) +
    geom_line(alpha = 0.6, linewidth = 0.25) +
    facet_wrap(~ parameter, scales = "free_y", ncol = 2) +
    labs(x = "retained iteration", y = "draw (angular traces centered)") +
    theme(legend.position = "bottom")
}

screen_alignment_plot <- function(screen_alignment_df) {
  ggplot(screen_alignment_df, aes(order, mean)) +
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2, fill = "#4477AA") +
    geom_line() + facet_wrap(~ density) +
    labs(x = "virtual order", y = "screening-fit alignment (68% HDI)")
}

alignment_order_plot <- function(alignment_summary) {
  ggplot(alignment_summary, aes(order, mean)) +
    geom_ribbon(aes(ymin = low, ymax = high), fill = "#4477AA", alpha = 0.2) +
    geom_line(colour = "#4477AA") + facet_wrap(~ density) +
    labs(x = "virtual order", y = "alignment A(kappa) cos(mu), 90% HDI")
}

alignment_density_plot <- function(density_summary) {
  ggplot(density_summary, aes(density, mean)) +
    geom_ribbon(aes(ymin = low, ymax = high), fill = "#CC6677", alpha = 0.2) +
    geom_line(colour = "#CC6677") + scale_x_continuous(trans = "log2", breaks = density_levels) +
    facet_wrap(~ order) + labs(x = "virtual density (per square m)", y = "alignment, 90% HDI")
}

elpd0_order_plot <- function(pdf0) {
  ggplot(pdf0, aes(order, mean)) +
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2, fill = "#4477AA") +
    geom_line() + geom_hline(yintercept = 1 / (2 * pi), linetype = "dashed") +
    facet_wrap(~ density) + labs(x = "virtual order", y = "predictive density at 0 (per radian)")
}

elpd0_density_plot <- function(pdf0_density) {
  ggplot(pdf0_density, aes(density, mean)) +
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2, fill = "#CC6677") +
    geom_line() + geom_hline(yintercept = 1 / (2 * pi), linetype = "dashed") +
    scale_x_continuous(trans = "log2", breaks = density_levels) +
    facet_wrap(~ order) + labs(x = "virtual density (per square m)", y = "density at 0 (per radian)")
}

polar_pdfs_plot <- function(pdf_summary, theta) {
  peak_max <- max(pdf_summary$high)
  pdf_summary$r_med <- 1 + pdf_summary$median / peak_max
  pdf_summary$r_low <- 1 + pdf_summary$low / peak_max
  pdf_summary$r_high <- 1 + pdf_summary$high / peak_max
  RMAX_polar <- 1 + 1.03
  order_levels <- sort(unique(pdf_summary$order))
  density_levels <- sort(unique(pdf_summary$density))
  order_colors <- setNames(
    scales::colour_ramp(c("white", "purple"))(
      seq(0.45, 0.95, length.out = length(order_levels))),
    as.character(order_levels))
  density_colors <- setNames(
    scales::colour_ramp(c("white", "red"))(
      seq(0.45, 0.95, length.out = length(density_levels))),
    as.character(density_levels))
  polar_rings <- rbind(
    data.frame(theta = theta %% (2 * pi), y = 1, ring = 1),
    data.frame(theta = theta %% (2 * pi), y = 2, ring = 2))
  polar_deco <- function(p) {
    p + geom_path(data = polar_rings, inherit.aes = FALSE,
      aes(x = theta, y = y, group = ring),
      colour = "grey25", linewidth = 0.2) +
      geom_segment(data = data.frame(x = 0, xend = 0, y = 1,
        yend = RMAX_polar), inherit.aes = FALSE,
      aes(x = x, xend = xend, y = y, yend = yend),
      colour = "grey20", linewidth = 0.3) +
      scale_x_continuous(limits = c(0, 2 * pi), expand = c(0, 0),
        breaks = NULL) +
      scale_y_continuous(limits = c(0, RMAX_polar), expand = c(0, 0),
        breaks = NULL) +
      coord_polar(start = -pi / 2, direction = -1)
  }
  p_row_orders <- polar_deco(ggplot(pdf_summary,
    aes(theta, r_med,
      colour = as.factor(density), group = as.factor(density))) +
    geom_ribbon(aes(ymin = r_low, ymax = r_high, fill = as.factor(density)),
      alpha = 0.5, colour = NA) +
    geom_line(linewidth = 0.2) +
    facet_wrap(~ as.factor(order), ncol = length(order_levels)) +
    scale_colour_manual(values = density_colors) +
    scale_fill_manual(values = density_colors, name = "Density") +
    theme_void(base_size = 8) + theme(legend.position = "right") +
    guides(colour = "none"))
  p_row_density <- polar_deco(ggplot(pdf_summary,
    aes(theta, r_med,
      colour = as.factor(order), group = as.factor(order))) +
    geom_ribbon(aes(ymin = r_low, ymax = r_high, fill = as.factor(order)),
      alpha = 0.5, colour = NA) +
    geom_line(linewidth = 0.2) +
    facet_wrap(~ as.factor(density), ncol = length(order_levels)) +
    scale_colour_manual(values = order_colors) +
    scale_fill_manual(values = order_colors, name = "Order") +
    theme_void(base_size = 8) + theme(legend.position = "right") +
    guides(colour = "none"))
  patchwork::wrap_plots(p_row_orders, p_row_density, ncol = 1)
}

coef_label <- c("Intercept", "Density", "Order", "Density \u00d7 Order")

coefficients_plot <- function(cs) {
  rng_x <- range(c(cs$low, cs$high))
  pad_x <- 0.12 * diff(rng_x)
  x_lo_c <- rng_x[1] - pad_x
  x_hi_c <- rng_x[2] + pad_x
  cspan_c <- x_hi_c - x_lo_c
  ggplot(cs, aes(mean, ypos, colour = comp)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35) +
    geom_errorbar(aes(xmin = low, xmax = high), width = 0.22,
      linewidth = 0.35) +
    geom_point(aes(shape = comp), size = 1.7) +
    scale_colour_manual(values = c(mu = "#8C564B", kappa = "#2CA02C"),
      guide = "none") +
    scale_shape_manual(values = c(mu = 18, kappa = 17), guide = "none") +
    scale_x_continuous(limits = c(x_lo_c, x_hi_c)) +
    scale_y_reverse(limits = c(4.7, -0.6), breaks = 1:4,
      labels = coef_label, expand = c(0, 0)) +
    # Framed upper-right legend, drawn by hand (this ggplot2 build drops
    # "upperleft"/"bottomright" positions without an error).
    annotate("rect", xmin = x_lo_c + 0.62 * cspan_c,
      xmax = x_lo_c + 0.97 * cspan_c, ymin = -0.5, ymax = 0.1,
      fill = "white", colour = "black", linewidth = 0.3) +
    annotate("point", x = x_lo_c + 0.665 * cspan_c, y = -0.35, shape = 18,
      colour = "#8C564B", size = 2.1) +
    annotate("text", x = x_lo_c + 0.73 * cspan_c, y = -0.35, label = "\u03bc",
      hjust = 0, vjust = 0.5, size = 2.6) +
    annotate("point", x = x_lo_c + 0.665 * cspan_c, y = -0.05, shape = 17,
      colour = "#2CA02C", size = 2.1) +
    annotate("text", x = x_lo_c + 0.73 * cspan_c, y = -0.05, label = "\u03ba",
      hjust = 0, vjust = 0.5, size = 2.6) +
    labs(x = "Coefficient value", y = "") +
    theme_classic(base_size = 8)
}

coef_comp_plot <- function(cs) {
  rng_cc <- range(c(cs$low, cs$high))
  span_cc <- diff(rng_cc)
  x_lo_cc <- rng_cc[1] - 0.10 * span_cc
  x_hi_cc <- rng_cc[2] + 0.10 * span_cc
  cc_span <- x_hi_cc - x_lo_cc
  ggplot(cs, aes(mean, ypos, colour = comp)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35) +
    geom_errorbar(aes(xmin = low, xmax = high), width = 0.22,
      linewidth = 0.35) +
    geom_point(aes(shape = comp), size = 1.5) +
    scale_colour_manual(values = c(mu = "#8C564B", kappa = "#2CA02C"),
      guide = "none") +
    scale_shape_manual(values = c(mu = 18, kappa = 17), guide = "none") +
    scale_x_continuous(limits = c(x_lo_cc, x_hi_cc)) +
    scale_y_reverse(limits = c(4.6, -0.6), breaks = 1:4,
      labels = coef_label, expand = c(0, 0)) +
    annotate("rect", xmin = x_lo_cc + 0.60 * cc_span,
      xmax = x_lo_cc + 0.96 * cc_span, ymin = -0.5, ymax = 0.1,
      fill = "white", colour = "black", linewidth = 0.3) +
    annotate("point", x = x_lo_cc + 0.65 * cc_span, y = -0.35, shape = 18,
      colour = "#8C564B", size = 1.9) +
    annotate("text", x = x_lo_cc + 0.73 * cc_span, y = -0.35, label = "\u03bc",
      hjust = 0, vjust = 0.5, size = 2.0) +
    annotate("point", x = x_lo_cc + 0.65 * cc_span, y = -0.05, shape = 17,
      colour = "#2CA02C", size = 1.9) +
    annotate("text", x = x_lo_cc + 0.73 * cc_span, y = -0.05, label = "\u03ba",
      hjust = 0, vjust = 0.5, size = 2.0) +
    labs(x = "Coefficient value", y = "") +
    theme_classic(base_size = 7)
}

model_red <- "#D62728"

fmt_num <- function(v) format(v, digits = 15, trim = TRUE)

individual_ppc_plot <- function(individual_ppc) {
  dens_levels_ppc <- sort(unique(individual_ppc$density))
  order_levels_ppc <- sort(unique(individual_ppc$order))
  nrows_ppc <- length(dens_levels_ppc)
  ncols_ppc <- length(order_levels_ppc)
  model_red <- "#D62728"
  fmt_num <- function(v) format(v, digits = 15, trim = TRUE)
  ppc_panel <- function(sub, show_y, top_row, left_col, last_panel, o, d) {
    sub <- sub[order(as.character(sub$animal_id)), , drop = FALSE]
    aids <- unique(sub$animal_id)
    n_an <- length(aids)
    xpos <- match(sub$animal_id, aids) +
      ifelse(sub$source == "data bootstrap", -0.14, 0.14)
    pdata <- data.frame(x = xpos[sub$source == "data bootstrap"],
      y = sub$mean[sub$source == "data bootstrap"],
      low = sub$low[sub$source == "data bootstrap"],
      high = sub$high[sub$source == "data bootstrap"])
    pmodel <- data.frame(x = xpos[sub$source == "posterior predictive"],
      y = sub$mean[sub$source == "posterior predictive"],
      low = sub$low[sub$source == "posterior predictive"],
      high = sub$high[sub$source == "posterior predictive"])
    p <- ggplot() +
      geom_segment(data = pdata, inherit.aes = FALSE,
        aes(x = x, y = y, xend = x, yend = low),
        colour = "black", linewidth = 0.32) +
      geom_segment(data = pdata, inherit.aes = FALSE,
        aes(x = x, y = y, xend = x, yend = high),
        colour = "black", linewidth = 0.32) +
      geom_segment(data = pmodel, inherit.aes = FALSE,
        aes(x = x, y = y, xend = x, yend = low),
        colour = model_red, linewidth = 0.35) +
      geom_segment(data = pmodel, inherit.aes = FALSE,
        aes(x = x, y = y, xend = x, yend = high),
        colour = model_red, linewidth = 0.35) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
      geom_point(data = pdata, inherit.aes = FALSE, aes(x = x, y = y),
        shape = 21, colour = "black", fill = "white", size = 1.4) +
      geom_point(data = pmodel, inherit.aes = FALSE, aes(x = x, y = y),
        shape = 16, colour = model_red, size = 1.3) +
      scale_x_continuous(
        limits = c(0.25, n_an + if (last_panel) 2.4 else 0.75),
        breaks = seq_len(n_an), labels = aids, expand = c(0, 0)) +
      scale_y_continuous(limits = c(-1.05, 1.05), expand = c(0, 0),
        breaks = seq(-1, 1, by = 0.5)) +
      labs(x = NULL, y = if (left_col) paste0("Density=", fmt_num(d),
        "\nAlignment") else "") +
      theme_classic(base_size = 7) +
      theme(axis.text.x = element_text(size = 5.5, angle = 90, hjust = 1,
        vjust = 0.5),
        axis.title.y = element_text(size = 7),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
    if (!show_y) {
      p <- p + theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
    }
    if (top_row) {
      p <- p + labs(title = paste0("Order=", fmt_num(o))) +
        theme(plot.title = element_text(size = 8, hjust = 0.5))
    }
    if (last_panel) {
      lg_x <- n_an + 0.30
      p <- p +
        annotate("rect", xmin = lg_x, xmax = lg_x + 1.85,
          ymin = -1.00, ymax = -0.56,
          fill = "white", colour = "black", linewidth = 0.3) +
        annotate("point", x = lg_x + 0.30, y = -0.66, shape = 16,
          colour = model_red, size = 1.5) +
        annotate("text", x = lg_x + 0.50, y = -0.66, hjust = 0,
          vjust = 0, size = 3, colour = "black", label = "Model") +
        annotate("point", x = lg_x + 0.30, y = -0.90, shape = 21,
          colour = "black", fill = "white", size = 1.5) +
        annotate("text", x = lg_x + 0.50, y = -0.90, hjust = 0,
          vjust = 0, size = 3, colour = "black", label = "Data")
    }
    p
  }
  ppc_panels <- list()
  for (ri in seq_len(nrows_ppc)) {
    for (ci in seq_len(ncols_ppc)) {
      d <- dens_levels_ppc[ri]
      o <- order_levels_ppc[ci]
      sub <- individual_ppc[individual_ppc$density == d &
        individual_ppc$order == o, , drop = FALSE]
      ppc_panels[[length(ppc_panels) + 1L]] <- ppc_panel(
        sub, show_y = (ci == 1L), top_row = (ri == 1L),
        left_col = (ci == 1L),
        last_panel = (ri == nrows_ppc && ci == ncols_ppc), o, d)
    }
  }
  p_individual_ppc <- do.call(patchwork::wrap_plots,
    c(ppc_panels, list(ncol = ncols_ppc)))
  list(plt = p_individual_ppc, width = max(6.8, 2.2 * ncols_ppc),
    height = max(3.2, 1.8 * nrows_ppc))
}

cell_ppc_plot <- function(cell_ppc) {
  dens_reds <- c("64" = "#980C13", "32" = "#B81419", "16" = "#D32020",
    "8" = "#ED392B", "4" = "#F75B40", "1" = "#FB7A5A")
  order_purples <- c("1" = "#4F1F8B", "0.8" = "#5F3C99", "0.6" = "#6F5BA8",
    "0.4" = "#7F7BB9", "0.2" = "#9490C3", "0" = "#AAA8D0")
  msize <- function(n) (2.6 + 0.12 * sqrt(n)) * 25.4 / 72.27
  pcc_thm <- function(bottom) {
    theme_classic(base_size = 8) +
      theme(axis.line = element_line(colour = "black", linewidth = 0.28),
        panel.grid = element_blank(),
        axis.text.x = if (bottom) element_text(size = 7) else element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"))
  }
  p_cell_left <- lapply(seq_along(dens_row), function(r) {
    d <- as.character(dens_row[r])
    sub <- cell_ppc[cell_ppc$density == dens_row[r], , drop = FALSE]
    sub$xd <- sub$order - 0.02
    sub$xm <- sub$order + 0.02
    sub$col <- dens_reds[d]
    p <- ggplot(sub, aes(mean, mean, colour = col)) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35) +
      geom_segment(data = sub[sub$source == "data bootstrap", , drop = FALSE],
        aes(x = xd, y = mean, xend = xd, yend = low),
        linewidth = 0.32) +
      geom_segment(data = sub[sub$source == "data bootstrap", , drop = FALSE],
        aes(x = xd, y = mean, xend = xd, yend = high),
        linewidth = 0.32) +
      geom_segment(data = sub[sub$source == "data bootstrap", , drop = FALSE],
        aes(x = xd - 0.02, y = low, xend = xd + 0.02, yend = low),
        linewidth = 0.32) +
      geom_segment(data = sub[sub$source == "data bootstrap", , drop = FALSE],
        aes(x = xd - 0.02, y = high, xend = xd + 0.02, yend = high),
        linewidth = 0.32) +
      geom_point(data = sub[sub$source == "data bootstrap", , drop = FALSE],
        aes(x = xd, size = msize(n_an)), shape = 1, fill = "white") +
      geom_segment(data = sub[sub$source == "posterior predictive", , drop = FALSE],
        aes(x = xm, y = mean, xend = xm, yend = low),
        linewidth = 0.35) +
      geom_segment(data = sub[sub$source == "posterior predictive", , drop = FALSE],
        aes(x = xm, y = mean, xend = xm, yend = high),
        linewidth = 0.35) +
      geom_segment(data = sub[sub$source == "posterior predictive", , drop = FALSE],
        aes(x = xm - 0.02, y = low, xend = xm + 0.02, yend = low),
        linewidth = 0.35) +
      geom_segment(data = sub[sub$source == "posterior predictive", , drop = FALSE],
        aes(x = xm - 0.02, y = high, xend = xm + 0.02, yend = high),
        linewidth = 0.35) +
      geom_point(data = sub[sub$source == "posterior predictive", , drop = FALSE],
        aes(x = xm, size = msize(n_an) * sqrt(0.9)), shape = 16) +
      scale_colour_identity(guide = "none") +
      scale_x_continuous(limits = c(-0.05, 1.05),
        breaks = seq(0, 1, by = 0.2),
        labels = format(seq(0, 1, by = 0.2), nsmall = 1)) +
      scale_y_continuous(limits = c(-0.55, 1.05), breaks = c(-0.5, 0, 0.5, 1)) +
      pcc_thm(r == length(dens_row)) +
      labs(title = paste0("Density = ", d),
        x = if (r == length(dens_row)) "Order" else "",
        y = if (r == 1L) "Alignment" else "") +
      scale_size_identity(guide = "none")
    p
  })
  p_cell_right <- lapply(seq_along(order_row), function(r) {
    o <- order_row[r]
    sub <- cell_ppc[cell_ppc$order == o, , drop = FALSE]
    sub$xd <- sub$density / 1.08
    sub$xm <- sub$density * 1.08
    sub$col <- order_purples[format(o, nsmall = 0)]
    p <- ggplot(sub, aes(mean, mean, colour = col)) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35) +
      geom_segment(data = sub[sub$source == "data bootstrap", , drop = FALSE],
        aes(x = xd, y = mean, xend = xd, yend = low),
        linewidth = 0.32) +
      geom_segment(data = sub[sub$source == "data bootstrap", , drop = FALSE],
        aes(x = xd, y = mean, xend = xd, yend = high),
        linewidth = 0.32) +
      geom_segment(data = sub[sub$source == "data bootstrap", , drop = FALSE],
        aes(x = xd / 1.02, y = low, xend = xd * 1.02, yend = low),
        linewidth = 0.32) +
      geom_segment(data = sub[sub$source == "data bootstrap", , drop = FALSE],
        aes(x = xd / 1.02, y = high, xend = xd * 1.02, yend = high),
        linewidth = 0.32) +
      geom_point(data = sub[sub$source == "data bootstrap", , drop = FALSE],
        aes(x = xd, size = msize(n_an)), shape = 1, fill = "white") +
      geom_segment(data = sub[sub$source == "posterior predictive", , drop = FALSE],
        aes(x = xm, y = mean, xend = xm, yend = high),
        linewidth = 0.35) +
      geom_segment(data = sub[sub$source == "posterior predictive", , drop = FALSE],
        aes(x = xm, y = mean, xend = xm, yend = low),
        linewidth = 0.35) +
      geom_segment(data = sub[sub$source == "posterior predictive", , drop = FALSE],
        aes(x = xm / 1.02, y = high, xend = xm * 1.02, yend = high),
        linewidth = 0.35) +
      geom_segment(data = sub[sub$source == "posterior predictive", , drop = FALSE],
        aes(x = xm / 1.02, y = low, xend = xm * 1.02, yend = low),
        linewidth = 0.35) +
      geom_point(data = sub[sub$source == "posterior predictive", , drop = FALSE],
        aes(x = xm, size = msize(n_an) * sqrt(0.9)), shape = 16) +
      scale_colour_identity(guide = "none") +
      scale_x_continuous(
        transform = scales::trans_new("log2", log2, function(x) 2 ^ x),
        limits = c(0.8, 75), breaks = c(1, 2, 4, 8, 16, 32, 64)) +
      scale_y_continuous(limits = c(-0.55, 1.05), breaks = c(-0.5, 0, 0.5, 1)) +
      pcc_thm(r == length(order_row)) +
      labs(title = paste0("Order = ", format(o, nsmall = 0)),
        x = if (r == length(order_row)) "Density" else "",
        y = "") +
      scale_size_identity(guide = "none")
    if (r == length(order_row)) {
      # framed Model/Data legend in the bottom-right panel, drawn by hand in
      # data coordinates (legend positions in corners are dropped silently
      # by this ggplot2 build).
      p <- p +
        annotate("rect", xmin = 32, xmax = 74, ymin = 0.48, ymax = 0.99,
          fill = "white", colour = "black", linewidth = 0.3) +
        annotate("point", x = 40, y = 0.85, shape = 16, colour = "black",
          size = 2.2) +
        annotate("text", x = 48, y = 0.85, label = "Model", hjust = 0,
          vjust = 0.5, size = 2.2) +
        annotate("point", x = 40, y = 0.60, shape = 1, colour = "black",
          fill = "white", size = 2.2) +
        annotate("text", x = 48, y = 0.60, label = "Data", hjust = 0,
          vjust = 0.5, size = 2.2)
    }
    p
  })
  patchwork::wrap_plots(
    patchwork::wrap_plots(p_cell_left, ncol = 1),
    patchwork::wrap_plots(p_cell_right, ncol = 1),
    ncol = 2)
}
