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

# Recency-weighted trailing directedness smoother (issue #142). Replaces the
# Gaussian x-domain smoothing: each animal's per-draw scores are smoothed in
# original chronological/recorded-step order with a finite normalized trailing
# exponential window, then the already-smoothed observations are mapped onto a
# common per-treatment x grid by nearest-cell binning (pure alignment; there is
# no second Gaussian smoother). A tracking gap or a missing/unusable heading
# resets an animal's history, so segments never contaminate one another and the
# filter never runs backward from the end.
#
#   alpha = 2 / (1 + sqrt(W));  decay = 1 - alpha
#   q_smooth[i] = sum(decay^k * q[i-k]) / sum(decay^k)   (k = 0..m-1)
#
# where m is the number of consecutive valid steps available through i, capped
# at the integer window W (1..20, default 10). W=1 returns the input exactly.
# This is a finite normalized trailing window, not an infinite recursive EMA.
#
# draws is [steps x draws] for ONE animal's full recorded sequence (including
# missing and gap rows). valid and gap_before have one entry per recorded step.
# Missing rows keep their NA; a valid gap-crossing row keeps its raw score but
# never enters the history. Call once per animal in chronological order.
smooth_directedness <- function(draws, valid, gap_before, window = 10L) {
  if (!is.numeric(window) || window != floor(window) || window < 1L || window > 20L)
    stop("window must be an integer in 1..20")
  draws <- as.matrix(draws)
  valid <- as.logical(valid)
  gap_before <- as.logical(gap_before)
  n <- nrow(draws)
  if (length(valid) != n || length(gap_before) != n)
    stop("draws rows must match valid and gap_before")
  if (any(!is.finite(draws[valid]))) stop("valid steps must have finite scores")
  out <- draws
  if (window == 1L) return(out)
  alpha <- 2 / (1 + sqrt(window))
  decay <- 1 - alpha
  history <- integer(0)
  for (i in seq_len(n)) {
    if (!valid[i] || gap_before[i]) {
      history <- integer(0)
      next
    }
    history <- c(history, i)
    if (length(history) > window) history <- history[(length(history) - window + 1L):length(history)]
    age <- (length(history) - 1L):0L
    w <- decay^age
    w <- w / sum(w)
    out[i, ] <- as.numeric(t(draws[history, , drop = FALSE]) %*% w)
  }
  out
}

# Nearest-cell index of x on a sorted grid (1-based; x within the grid range).
nearest_cell <- function(x, grid) {
  g <- length(grid)
  gi <- findInterval(x, grid)            # 0..g: grid[gi] <= x < grid[gi+1]
  left <- pmin(pmax(gi, 1L), g)          # nearest cell at/below x
  right <- pmin(pmax(gi + 1L, 1L), g)    # nearest cell at/above x
  idx <- ifelse(abs(x - grid[left]) <= abs(x - grid[right]), left, right)
  pmin(pmax(idx, 1L), g)
}

summarize_curve <- function(curve, grid, treatment, tag = NA_character_, hdi_prob = 0.95) {
  valid <- rowSums(is.finite(curve)) == ncol(curve)
  if (!any(valid)) return(NULL)
  z <- curve[valid, , drop = FALSE]
  bounds <- t(apply(z, 1, hdi, hdi_prob))
  data.frame(x = grid[valid], mean = rowMeans(z), low = bounds[, 1], high = bounds[, 2],
    treatment = treatment, tag = tag)
}

# Curves for one view. draws is [nrow(df_track) x draws] (steps x draws), the
# per-step posterior directedness (default) or directed-component probability.
# Stage 1: smooth each animal's draws in recorded-step order, aligning the
# fitted scores (rows of draws) back onto the full track via track_positions
# (1-based positions of the usable rows in track_w_treat). Stage 2: nearest-cell
# binning onto a common per-treatment grid (pure alignment; equal-animal
# per-draw average; population shown only where >=2 animals contribute).
make_curves <- function(kind, draws, hdi_prob = 0.95, window = 10L) {
  n_full <- nrow(track_w_treat)
  pos <- as.integer(track_positions)
  full_tag <- track_w_treat$tag
  full_gap <- as.logical(track_w_treat$gap_before)
  full_to_use <- integer(n_full)
  full_to_use[pos] <- seq_along(pos)
  smoothed <- matrix(NA_real_, nrow = nrow(df_track), ncol = ncol(draws))
  # New tag block starts where the tag differs from the previous row.
  bounds <- c(which(full_tag[-1] != full_tag[-length(full_tag)]) + 1L, n_full + 1L)
  for (k in seq_along(bounds[-length(bounds)])) {
    t0 <- bounds[k]
    t1 <- bounds[k + 1L] - 1L
    seg <- t0:t1
    use_idx <- full_to_use[seg]
    if (all(use_idx == 0L)) next
    valid <- use_idx > 0L
    scores_full <- matrix(NA_real_, nrow = length(seg), ncol = ncol(draws))
    scores_full[valid, ] <- draws[use_idx[valid], , drop = FALSE]
    sm <- smooth_directedness(scores_full, valid, full_gap[seg], window = window)
    smoothed[use_idx[valid], ] <- sm[valid, , drop = FALSE]
  }
  # x per usable row (unchanged meaning from the Gaussian path).
  x <- df_track$dist_to_dest_grad_m
  if (kind == "steps_from_end") {
    x <- ave(seq_len(nrow(df_track)), df_track$tag,
      FUN = function(ii) (seq_along(ii) - length(ii)) * 60)
  }
  G <- if (quick_run) 400L else 2000L
  population <- individual <- list()
  for (t in levels(df_track$treatment)) {
    eligible <- which(df_track$treatment == t & is.finite(x))
    if (kind == "steps_from_end") eligible <- eligible[x[eligible] >= -110 * 60]
    if (!length(eligible)) next
    grid <- seq(min(x[eligible]), max(x[eligible]), length.out = G)
    total <- matrix(0, G, ncol(draws))
    counts <- numeric(G)
    tags <- split(eligible, droplevels(df_track$tag[eligible]))
    for (tag in names(tags)) {
      ii <- tags[[tag]]
      if (length(ii) < 2L) next
      # Nearest-cell binning (pure alignment; no cross-cell weighting).
      cell <- nearest_cell(x[ii], grid)
      curve <- matrix(NA_real_, G, ncol(draws))
      for (cc in unique(cell)) {
        sel <- cell == cc
        curve[cc, ] <- colMeans(smoothed[ii[sel], , drop = FALSE])
      }
      valid <- grid >= min(x[ii]) & grid <= max(x[ii]) & !is.na(curve[, 1])
      total[valid, ] <- total[valid, , drop = FALSE] + curve[valid, , drop = FALSE]
      counts[valid] <- counts[valid] + 1L
      curve[!valid, ] <- NA_real_
      individual[[paste(t, tag)]] <- summarize_curve(curve, grid, t, tag, hdi_prob)
    }
    curve <- total / pmax(counts, 1)
    curve[counts < 2L, ] <- NA_real_
    population[[t]] <- summarize_curve(curve, grid, t, hdi_prob)
  }
  list(population = do.call(rbind, population), individual = do.call(rbind, individual))
}
