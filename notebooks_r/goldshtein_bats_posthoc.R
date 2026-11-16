# Bats (Goldshtein): post-fit computation.
#
# Sourced by notebooks_r/goldshtein_bats.ipynb (no standalone output). R
# twin of the Python notebook's local posthoc code: the KDE and kappa
# helpers, the response/likelihood matrices, the summary tables, the PPC
# simulation helpers and the curve-smoothing machinery. The notebook
# keeps the data, both fits, the seeds and printed diagnostics, and
# calls the figure builders from goldshtein_bats_plots.R.

theta_grid <- seq(0, 2 * pi, length.out = 360)
KAPPA_UNIFORM <- 1e-6
# Preserve the source post-hoc helpers' near-uniform proxy in display/PPC
# summaries. The fitted mixture likelihood and responsibilities use the
# exact uniform-circle density, as in Python c53 and the responsibility helper.

kde_vm <- function(y, theta = theta_grid, kappa = 30) {
  # Matrix form of the grid-point averages: exp(vm_logpdf(t, y, kappa))
  # = exp(kappa*cos(t - y)) / (2*pi*I0(kappa)); rowMeans takes the mean
  # over the observations at each grid point, as the vapply did.
  rowMeans(exp(kappa * cos(outer(theta, y, `-`))) /
    (2 * pi * besselI(kappa, 0, expon.scaled = FALSE)))
}

estimate_kappa <- function(rbar, n) {
  rbar <- min(max(rbar, 0), 0.999999999)
  if (rbar < 0.53) {
    k <- 2 * rbar + rbar^3 + (5 / 6) * rbar^5
  } else if (rbar < 0.85) {
    k <- -0.4 + 1.39 * rbar + 0.43 / (1 - rbar)
  } else {
    k <- 1 / (rbar^3 - 4 * rbar^2 + 3 * rbar)
  }
  if (n > 1L) {
    if (k < 2) {
      k <- if (k > 0) max(k - 2 / (n * k), 0) else 0
    } else {
      k <- (n - 1)^3 * k / (n^3 + n)
    }
  }
  k
}

response_parameters <- function(fit, newdata, mixture = FALSE,
                                population = FALSE, draw_ids = NULL) {
  # brms returns [draws x grid]; transpose explicitly to [grid x draws].
  get_lp <- function(dpar = NULL, nlpar = NULL, transform = FALSE) {
    t(posterior_linpred(fit, newdata = newdata, dpar = dpar, nlpar = nlpar,
      transform = transform, re_formula = if (population) NA else NULL,
      draw_ids = draw_ids))
  }
  if (population) {
    # Circular tag vectors are nonlinear population coefficients, so a
    # re_formula alone cannot exclude them. Evaluate the relevant nlpars.
    mu <- if (mixture) get_lp(nlpar = "mu0") + get_lp(nlpar = "mutreat") else get_lp(nlpar = "beta0")
  } else {
    mu <- get_lp(dpar = "mu")
  }
  kappa <- get_lp(dpar = "kappa", transform = TRUE)
  p <- if (mixture) get_lp(dpar = "pgoal", transform = TRUE) else matrix(1, nrow(mu), ncol(mu))
  stopifnot(nrow(mu) == nrow(newdata), identical(dim(mu), dim(kappa)),
            identical(dim(mu), dim(p)), all(kappa > 0), all(p >= 0 & p <= 1))
  list(mu = wrap(mu), kappa = kappa, p = p)
}

explicit_loo <- function(fit, mixture = FALSE) {
  d <- posterior::as_draws_df(fit)
  ll <- matrix(0, nrow = nrow(d), ncol = nrow(df_track))
  for (start in seq.int(1L, nrow(df_track), by = 128L)) {
    ii <- seq.int(start, min(start + 127L, nrow(df_track)))
    block <- df_track[ii, , drop = FALSE]
    par <- response_parameters(fit, block, mixture = mixture)
    y <- matrix(rep(block$HeadingRadians, ncol(par$mu)), nrow = length(ii))
    log_goal <- vm_logpdf(y, par$mu, par$kappa)
    if (mixture) {
      a <- log(par$p) + log_goal
      b <- log1p(-par$p) - log(2 * pi)
      maximum <- pmax(a, b)
      log_density <- maximum + log(exp(a - maximum) + exp(b - maximum))
    } else {
      log_density <- log_goal
    }
    ll[, ii] <- t(log_density)
  }
  # Scale each likelihood column before exponentiating; relative_eff is
  # invariant to this scaling. Preserve chain IDs rather than flattening.
  likelihood <- exp(sweep(ll, 2, apply(ll, 2, max), "-"))
  r_eff <- loo::relative_eff(likelihood, chain_id = d[[".chain"]])
  loo::loo(ll, r_eff = r_eff)
}

delta_row <- function(ii, label) {
  delta <- pointwise_delta[ii]
  data.frame(label = label, n = length(ii), delta_elpd = sum(delta),
    se_delta_elpd = if (length(ii) > 1L) sqrt(length(ii) * var(delta)) else NA_real_)
}

mass_vm <- function(mu, kappa, half_width = pi / 4, grid_n = 721L) {
  theta <- seq(-half_width, half_width, length.out = grid_n)
  weights <- c(0.5, rep(1, grid_n - 2L), 0.5) * diff(theta)[1]
  value <- numeric(length(mu))
  for (j in seq_along(theta)) value <- value + weights[j] * exp(vm_logpdf(theta[j], mu, kappa))
  value
}

parameter_summary <- function(par, grid) {
  do.call(rbind, lapply(seq_len(nrow(grid)), function(j) {
    mu <- par$mu[j, ]
    k <- par$kappa[j, ]
    p <- par$p[j, ]
    f_dir <- exp(vm_logpdf(0, mu, k))
    f_unif <- exp(vm_logpdf(0, 0, KAPPA_UNIFORM))
    values <- list(mu = mu, p = p, kappa = k, mrl = A1_safe(k),
      dens0_dir = f_dir, dens0_mix = p * f_dir + (1 - p) * f_unif,
      mass45 = p * mass_vm(mu, k) + (1 - p) * 0.25,
      align0 = p * A1_safe(k) * cos(mu) + (1 - p) * A1_safe(KAPPA_UNIFORM))
    do.call(rbind, lapply(names(values), function(q) {
      data.frame(grid[j, , drop = FALSE], quantity = q,
        as.list(interval(values[[q]], angular = q == "mu")))
    }))
  }))
}

# Posterior directedness with the EXACT uniform-circle background: its
# -log(2*pi) cancels the von Mises density constant, so
# logB = kappa*(cos(y - mu) - 1) - log(besselI(kappa, 0,
# expon.scaled = TRUE)), the same stable log-I0 form as vm_logpdf
# (finite at any kappa).
log_evid_ratio <- function(y, mu, kappa) {
  kappa * (cos(y - mu) - 1) - log(besselI(kappa, 0, expon.scaled = TRUE))
}

hist_pdf <- function(y) hist(wrap(y), breaks = bins, plot = FALSE)$density

sim_one_tag <- function(j, d) {
  n <- length(tag_rows[[j]])
  goal <- runif(n) < ppc_par$p[j, d]
  # One n-draw call per component: the draw stream of each call matches
  # the n-scalar-call stream exactly (probed), at a large constant-factor
  # speedup. The two pools are drawn separately rather than interleaved
  # row by row, so replicated steps here are stream-different (though
  # distribution-identical) to the earlier interleaved version.
  gd <- as.numeric(circular::rvonmises(n, vmu(ppc_par$mu[j, d]),
    ppc_par$kappa[j, d], control.circular = vm_ctrl))
  fd <- as.numeric(circular::rvonmises(n, vmu(0), KAPPA_UNIFORM,
    control.circular = vm_ctrl))
  wrap(ifelse(goal, gd, fd))
}

hier_y <- function(ys) {
  selected <- sample.int(length(ys), length(ys), replace = TRUE)
  unlist(lapply(selected, function(j) {
    y <- ys[[j]]
    y[sample.int(length(y), length(y), replace = TRUE)]
  }), use.names = FALSE)
}

band_table <- function(draws, x, treatment, source, kind, equal_tail = FALSE) {
  # draws are [grid x draws]; summarize each grid row independently.
  bounds <- if (equal_tail) t(apply(draws, 1, quantile, probs = c(0.025, 0.975))) else t(apply(draws, 1, hdi))
  data.frame(x = x, median = apply(draws, 1, median), low = bounds[, 1],
    high = bounds[, 2], treatment = treatment, source = source, kind = kind)
}
# Recency-weighted trailing directedness smoother (issue #142).
# R twin of smooth_directedness() in notebooks/goldshtein_bats_posthoc.py.
#
# scores: [steps x draws] matrix for ONE animal in recorded-step order.
# valid: per-step logical, FALSE = missing/unusable heading (the row keeps
#   its values and resets the trailing history).
# gap_before: per-step logical, TRUE = heading crosses a tracking gap (the
#   row keeps its raw score but never enters the history, and history
#   restarts at the following step).
#
# alpha = 2 / (1 + sqrt(W)); decay = 1 - alpha; a finite trailing
# exponential window with normalized partial weights at segment starts --
# NOT an infinite recursive EMA. window = 1L returns the input exactly.
# Never apply this in order of a plotted coordinate (distance to the
# destination is not monotone in time) and never backward.
smooth_directedness <- function(scores, valid, gap_before, window = 10L) {
  if (is.logical(window)) stop("window must be an integer in 1..20")
  window <- as.integer(window)
  if (length(window) != 1L || is.na(window) || window < 1L || window > 20L) {
    stop("window must be an integer in 1..20")
  }
  if (!is.matrix(scores)) stop("scores must be a [steps x draws] matrix")
  if (length(valid) != nrow(scores) || length(gap_before) != nrow(scores)) {
    stop("expected scores[steps, draws] and aligned step masks")
  }
  valid <- as.logical(valid)
  gap_before <- as.logical(gap_before)
  vs <- scores[valid, , drop = FALSE]
  if (nrow(vs) > 0L && (any(!is.finite(vs)) || any(vs < 0) || any(vs > 1))) {
    stop("valid steps must have finite scores in [0, 1]")
  }

  out <- scores
  if (window == 1L) return(out)

  alpha <- 2 / (1 + sqrt(window))
  decay <- 1 - alpha
  history <- integer(0)
  for (i in seq_len(nrow(scores))) {
    if (!valid[i] || gap_before[i]) {
      # Missing rows retain their missing values; valid gap-crossing rows
      # retain their raw scores, but never enter the history.
      history <- integer(0)
      next
    }
    history <- c(history, i)
    if (length(history) > window) {
      history <- history[seq.int(length(history) - window + 1L, length(history))]
    }
    w <- decay ^ rev(seq_along(history) - 1L)
    out[i, ] <- colSums(scores[history, , drop = FALSE] * (w / sum(w)))
  }
  out
}

summarize_curve <- function(curve, grid, treatment, tag = NA_character_, hdi_prob = 0.95) {
  valid <- rowSums(is.finite(curve)) == ncol(curve)
  if (!any(valid)) return(NULL)
  z <- curve[valid, , drop = FALSE]
  bounds <- t(apply(z, 1, hdi, hdi_prob))
  data.frame(x = grid[valid], mean = rowMeans(z), low = bounds[, 1], high = bounds[, 2],
    treatment = treatment, tag = tag)
}

# Curves for one view (issue #142). draws is [nrow(df_track) x draws]: the
# per-step posterior directedness, one row per recorded step.
#
# The trailing exponential filter smooth_directedness() is the ONLY smoother:
# it runs per animal along RECORDED steps over the animal's full sequence
# (missing rows and tracking gaps reset the history; gap-crossing headings
# keep their raw score but never enter a window). Only then are the
# already-smoothed step scores mapped onto the plotting grid by equal-width
# binning: each joint draw's value in a bin is the plain mean of that
# animal's smoothed steps whose x falls in the bin. Binning is aggregation
# over repeated / nonmonotone x values at the SAME x -- not a decay
# smoother: it fills no bin an animal has no steps in, so unsupported
# ranges and tracking gaps are never bridged and nothing is extrapolated.
# Equal-animal averaging happens per joint draw; grid cells with fewer than
# two contributing animals are dropped; means and intervals are computed
# ACROSS the smoothed joint draws (never of interval endpoints). W defaults
# to 10 (the issue's default); window = 1L plots the raw step scores.
make_curves <- function(kind, draws, hdi_prob = 0.90, window = 10L, bin_width = NULL,
                        min_tags_at_x = 2L) {
  x <- df_track$dist_to_dest_grad_m
  in_window <- rep(TRUE, nrow(df_track))
  if (kind == "steps_from_end") {
    x <- ave(seq_len(nrow(df_track)), df_track$tag,
      FUN = function(ii) (seq_along(ii) - length(ii)) * 60)
    in_window <- x >= -110 * 60 & x <= 0
  } else if (kind != "distance") {
    stop("kind must be 'distance' or 'steps_from_end'")
  }
  if (is.null(bin_width)) bin_width <- if (kind == "steps_from_end") 60 else 250
  if (!(bin_width > 0)) stop("bin_width must be positive")
  gap <- !is.na(df_track$gap_before) & as.logical(df_track$gap_before)
  eligible <- in_window & is.finite(x)
  population <- individual <- list()

  # Per-treatment bin grids (bin centres are the plotting grid), built from
  # eligible rows; smoothing itself ignores the plotting range entirely.
  grids <- list()
  for (t in levels(df_track$treatment)) {
    xs <- x[eligible & df_track$treatment == t]
    if (!length(xs)) next
    lo_edge <- floor(min(xs) / bin_width) * bin_width
    n_bins <- ceiling((max(xs) - lo_edge) / bin_width) + 1L
    grids[[t]] <- list(lo = lo_edge, n = n_bins,
      x_mid = lo_edge + bin_width * (seq_len(n_bins) - 1L),  # bin-start labels
      tot = matrix(0, n_bins, ncol(draws)), cnt = numeric(n_bins))
  }

  # One pass per animal: smooth in recorded-step order, then bin.
  tag_rows <- split(seq_len(nrow(df_track)), droplevels(df_track$tag))
  for (tag in names(tag_rows)) {
    ii <- tag_rows[[tag]]
    t <- as.character(df_track$treatment[ii[1]])
    g <- grids[[t]]
    if (is.null(g)) next
    xs <- x[ii]
    valid <- is.finite(xs) & rowSums(is.finite(draws[ii, , drop = FALSE])) == ncol(draws)
    sm <- smooth_directedness(draws[ii, , drop = FALSE], valid, gap[ii], window)
    b <- floor((xs - g$lo) / bin_width) + 1L
    keep <- which(eligible[ii] & valid & b >= 1L & b <= g$n)
    if (!length(keep)) next
    bg <- factor(b[keep], levels = seq_len(g$n))
    cnt <- as.numeric(table(bg))
    rs <- rowsum(sm[keep, , drop = FALSE], group = bg)
    sums <- matrix(0, g$n, ncol(draws))
    sums[match(rownames(rs), as.character(seq_len(g$n))), ] <- as.matrix(rs)
    have <- cnt > 0
    tag_curve <- matrix(NA_real_, g$n, ncol(draws))
    tag_curve[have, ] <- sums[have, , drop = FALSE] / cnt[have]

    individual[[paste(t, tag)]] <- summarize_curve(tag_curve, g$x_mid, t, tag, hdi_prob)
    g$tot[have, ] <- g$tot[have, , drop = FALSE] + tag_curve[have, , drop = FALSE]
    g$cnt[have] <- g$cnt[have] + 1
    grids[[t]] <- g
  }

  # Equal-animal average per joint draw; minimum-two-animal population rule;
  # intervals across draws.
  for (t in names(grids)) {
    g <- grids[[t]]
    if (!any(g$cnt > 0)) next
    curve <- matrix(NA_real_, g$n, ncol(draws))
    have <- g$cnt > 0
    curve[have, ] <- g$tot[have, , drop = FALSE] / g$cnt[have]
    curve[g$cnt < min_tags_at_x, ] <- NA_real_
    population[[t]] <- summarize_curve(curve, g$x_mid, t, hdi_prob = hdi_prob)
  }
  list(population = do.call(rbind, population), individual = do.call(rbind, individual))
}
