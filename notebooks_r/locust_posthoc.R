# Locust direction mixture: post-fit computation.
#
# Sourced by notebooks_r/locust_direction_mixture_model.ipynb (no
# standalone output). R twin of notebooks/locust_posthoc.py:
# wrapped von Mises / mixture densities, posterior extraction, the
# diagnostic gate, bootstrap KDEs, explicit PSIS-LOO, per-individual
# PPC layers and the hypothesis / alignment tables.

softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

vm_logpdf <- function(y, mu, kappa) {
  kappa * (cos(y - mu) - 1) - log(2 * pi) -
    log(besselI(kappa, 0, expon.scaled = TRUE))
}

log_add <- function(a, b) {
  m <- pmax(a, b)
  m + log(exp(a - m) + exp(b - m))
}

A1_safe <- function(kappa) {
  besselI(kappa, 1, expon.scaled = TRUE) / besselI(kappa, 0, expon.scaled = TRUE)
}

hdi <- function(x, prob = 0.95) {
  x <- sort(as.numeric(x))
  stopifnot(length(x) > 1L, all(is.finite(x)))
  gap <- min(length(x) - 1L, floor(prob * length(x)))
  starts <- seq_len(length(x) - gap)
  j <- starts[which.min(x[starts + gap] - x[starts])]
  c(low = x[j], high = x[j + gap])
}

# Every prediction/density matrix in this notebook is [grid x draws].
# Every prediction/density matrix in this notebook is [grid x draws].
vm_grid <- function(theta, mu, kappa) {
  stopifnot(length(mu) == length(kappa))
  G <- length(theta)
  K <- matrix(rep(kappa, each = G), nrow = G)
  delta <- outer(theta, mu, "-")
  exp(K * (cos(delta) - 1) - log(2 * pi) -
    matrix(rep(log(besselI(kappa, 0, expon.scaled = TRUE)), each = G), nrow = G))
}

mixture_grid <- function(theta, mu, kappa, weight) {
  stopifnot(nrow(mu) == 2L, identical(dim(mu), dim(kappa)),
            ncol(mu) == length(weight))
  W <- matrix(rep(weight, each = length(theta)), nrow = length(theta))
  W * vm_grid(theta, mu[1, ], kappa[1, ]) +
    (1 - W) * vm_grid(theta, mu[2, ], kappa[2, ])
}

clip_mixture_weight <- function(w) {
  w1 <- pmax(w, 1e-12)
  w2 <- pmax(1 - w, 1e-12)
  w1 / (w1 + w2)
}

mixture_at <- function(theta, state, clip = FALSE) {
  W <- state$weight
  stopifnot(is.matrix(W), identical(dim(W), dim(state$muone)))
  if (clip) {
    W <- clip_mixture_weight(W)
    log_w1 <- log(W)
    log_w2 <- log1p(-W)
  } else {
    log_w1 <- -softplus(-state$eta)
    log_w2 <- -softplus(state$eta)
  }
  log_add(log_w1 + vm_logpdf(theta, state$muone, state$kone),
          log_w2 + vm_logpdf(theta, state$mutwo, state$ktwo))
}

circular_kde <- function(theta, y, bw) {
  # The grid is passed as an explicit circular object: density.circular
  # would otherwise coerce the raw numeric grid internally, which emits
  # a default-value coercion warning on every call (thousands of them in
  # the bootstrap cells).
  y <- circular::circular(wrap(y), units = "radians", template = "none",
                          modulo = "asis")
  z <- circular::circular(theta, units = "radians", template = "none",
                          modulo = "asis")
  circular::density.circular(y, z = z, bw = bw, kernel = "vonmises")$y
}

kde_basis <- function(theta, y, bw) {
  # The KDE is linear in observation weights. Repeating a single angle
  # twice keeps density.circular's kernel matrix two dimensional and
  # gives the same kernel as that angle alone.
  vapply(y, function(yi) circular_kde(theta, rep(yi, 2L), bw),
         numeric(length(theta)))
}

density_frame <- function(theta, draws) {
  qs <- t(apply(draws, 1, quantile, probs = c(0.025, 0.5, 0.975)))
  data.frame(theta = theta, low = qs[, 1], center = qs[, 2], high = qs[, 3])
}

extract_state <- function(fit, bimodal = TRUE) {
  draws <- posterior::as_draws_array(fit)
  draw_index <- expand.grid(.iteration = seq_len(dim(draws)[1]),
    .chain = seq_len(dim(draws)[2]), KEEP.OUT.ATTRS = FALSE)
  get <- function(name) {
    stopifnot(name %in% posterior::variables(draws))
    as.numeric(draws[, , name])
  }
  predict_matrix <- function(dpar, transform = FALSE) {
    x <- t(posterior_linpred(fit, newdata = animal_grid,
                            dpar = dpar, transform = transform))
    stopifnot(identical(dim(x), c(num_individuals, nrow(draw_index))))
    x
  }
  m1 <- wrap(predict_matrix("mu"))
  k1 <- predict_matrix("kappa", transform = TRUE)
  p1 <- get("b_popone_Intercept")
  kp1 <- softplus(get("b_kappa_Intercept"))
  if (!bimodal) {
    return(list(draws = draw_index, muone = m1, kone = k1,
      pop_muone = matrix(wrap(p1), nrow = 1L),
      pop_kone = matrix(kp1, nrow = 1L),
      quantities = c(list(mu_pop = wrap(p1), mu_kappa = get("kmu_one"),
        kappa_log_pop = get("b_kappa_Intercept"), kappa_pop = kp1,
        kappa_sigma = get("sd_animal__kappa_Intercept")),
        setNames(lapply(seq_len(num_individuals), function(i) m1[i, ]),
                 paste0("mu[", seq_len(num_individuals), "]")),
        setNames(lapply(seq_len(num_individuals), function(i) k1[i, ]),
                 paste0("kappa[", seq_len(num_individuals), "]")))))
  }
  m2 <- wrap(predict_matrix("mutwo"))
  k2 <- predict_matrix("kappatwo", transform = TRUE)
  eta <- predict_matrix("eta")
  p2 <- get("b_poptwo_Intercept")
  kp2 <- softplus(get("b_kappatwo_Intercept"))
  eta_pop <- get("b_eta_Intercept")
  first_component <- get("first_component")
  stopifnot(all(first_component %in% 1:2),
            all((first_component == 2) == (p1 > p2)))
  swap <- first_component == 2
  # Read the saved Stan index for all summaries and predictive draws.
  reorder_pair <- function(a, b) {
    stopifnot(identical(dim(a), dim(b)), ncol(a) == length(swap))
    first <- a
    second <- b
    first[, swap] <- b[, swap, drop = FALSE]
    second[, swap] <- a[, swap, drop = FALSE]
    list(first, second)
  }
  means <- reorder_pair(m1, m2)
  kappas <- reorder_pair(k1, k2)
  pop_means <- reorder_pair(matrix(p1, nrow = 1L), matrix(p2, nrow = 1L))
  pop_kappas <- reorder_pair(matrix(kp1, nrow = 1L), matrix(kp2, nrow = 1L))
  pool_k <- reorder_pair(matrix(get("kmu_one"), nrow = 1L),
                         matrix(get("kmu_two"), nrow = 1L))
  sigma_k <- reorder_pair(matrix(get("sd_animal__kappa_Intercept"), nrow = 1L),
                          matrix(get("sd_animal__kappatwo_Intercept"), nrow = 1L))
  eta[, swap] <- -eta[, swap, drop = FALSE]
  eta_pop[swap] <- -eta_pop[swap]
  w <- plogis(eta)
  wp <- plogis(eta_pop)
  quantities <- list(logit_sigma = get("sd_animal__eta_Intercept"),
                     weights_pop = wp)
  for (k in 1:2) {
    quantities[[paste0("mu_pop[", k, "]")]] <- as.numeric(pop_means[[k]])
    quantities[[paste0("mu_kappa[", k, "]")]] <- as.numeric(pool_k[[k]])
    quantities[[paste0("kappa_pop[", k, "]")]] <- as.numeric(pop_kappas[[k]])
    quantities[[paste0("kappa_sigma[", k, "]")]] <- as.numeric(sigma_k[[k]])
    for (i in seq_len(num_individuals)) {
      quantities[[paste0("mu[", i, ",", k, "]")]] <- means[[k]][i, ]
      quantities[[paste0("kappa[", i, ",", k, "]")]] <- kappas[[k]][i, ]
    }
  }
  for (i in seq_len(num_individuals)) quantities[[paste0("weights[", i, "]")]] <- w[i, ]
  list(draws = draw_index, muone = means[[1]], mutwo = means[[2]],
       kone = kappas[[1]], ktwo = kappas[[2]], weight = w, eta = eta,
       pop_muone = pop_means[[1]], pop_mutwo = pop_means[[2]],
       pop_kone = pop_kappas[[1]], pop_ktwo = pop_kappas[[2]],
       pop_weight = matrix(wp, nrow = 1L), pop_eta = matrix(eta_pop, nrow = 1L),
       quantities = quantities)
}

diagnose_fit <- function(fit, state) {
  # Same retained-flag probe as the P1 helpers; absence is not zero.
  divergences <- tryCatch({
    samples <- fit$fit@sim$samples
    if (is.null(samples) || !length(samples)) {
      NA_real_
    } else {
      flags <- lapply(samples, function(ch) ch[["divergent"]])
      if (any(vapply(flags, is.null, logical(1)))) NA_real_ else sum(unlist(flags))
    }
  }, error = function(e) NA_real_)
  if (is.na(divergences)) {
    cat("Divergences: not retained in this fit object.\n")
  } else {
    cat("Divergences:", divergences, "\n")
  }
  table <- do.call(rbind, lapply(names(state$quantities), function(p) {
    x <- chain_matrix(state$quantities[[p]], state$draws)
    stopifnot(all(is.finite(x)))
    angular <- grepl("^mu(\\[|_pop($|\\[))", p)
    z <- if (angular) wrap(x - circ_mean(x)) else x
    data.frame(parameter = p,
      mean = if (angular) circ_mean(x) else mean(x),
      sd = if (angular) sqrt(-2 * log(min(1, resultant(x)))) else sd(as.numeric(x)),
      rhat = if (angular) Rhat_unwrap(x, method = "sine") else posterior::rhat(x),
      ess_bulk = posterior::ess_bulk(z), row.names = NULL)
  }))
  # Angular R-hat uses the Rhat_unwrap sine/cosine convention. Every ESS
  # and R-hat calculation receives an [iteration x chain] matrix.
  failed <- table[!is.finite(table$rhat) | table$rhat >= 1.01 |
    !is.finite(table$ess_bulk) | table$ess_bulk <= 400, , drop = FALSE]
  print(failed, row.names = FALSE)
  stopifnot(nrow(failed) == 0L, is.na(divergences) || divergences == 0)
  table
}

explicit_loo <- function(state, bimodal = TRUE, chunk_size = 128L) {
  S <- nrow(state$draws)
  ll <- matrix(NA_real_, nrow = S, ncol = nrow(df))
  for (start in seq.int(1L, nrow(df), by = chunk_size)) {
    ii <- seq.int(start, min(start + chunk_size - 1L, nrow(df)))
    a <- df$animal_id[ii] + 1L
    y <- matrix(rep(df$direction[ii], S), nrow = length(ii))
    l1 <- vm_logpdf(y, state$muone[a, , drop = FALSE], state$kone[a, , drop = FALSE])
    if (bimodal) {
      eta <- state$eta[a, , drop = FALSE]
      l2 <- vm_logpdf(y, state$mutwo[a, , drop = FALSE], state$ktwo[a, , drop = FALSE])
      ll[, ii] <- t(log_add(-softplus(-eta) + l1, -softplus(eta) + l2))
    } else {
      ll[, ii] <- t(l1)
    }
  }
  stopifnot(all(is.finite(ll)))
  # Scale each likelihood column before exponentiating; relative_eff is
  # invariant to this scaling. Preserve chain IDs rather than flattening.
  likelihood <- exp(sweep(ll, 2, apply(ll, 2, max), "-"))
  r_eff <- loo::relative_eff(likelihood, chain_id = state$draws[[".chain"]])
  # Summary columns flow into the ELPD table the notebook prints; the
  # full per-draw pareto-k table is deliberately not printed (long
  # output). Suppress loo()'s redundant high-k nag.
  result <- suppressWarnings(loo::loo(ll, r_eff = r_eff))
  result
}

# Per-individual PPC layers: posterior predictive band, bootstrap KDE
# band, plug-in model curve, inverse-density residual ticks and the
# component-mean arrows. Seeds are the per-animal values used by the
# Python twin's individual checks.
individual_ppc_layers <- function(state, df, rows_by_animal,
                                num_individuals, n_boot_kde) {
individual_theta <- -pi + (0:359) * (2 * pi / 360)
individual_curves <- vector("list", num_individuals)
individual_ticks <- vector("list", num_individuals)
individual_arrows <- vector("list", num_individuals)
for (i in seq_len(num_individuals)) {
  y <- wrap(df$direction[rows_by_animal[[i]]])
  mu <- rbind(state$muone[i, ], state$mutwo[i, ])
  kappa <- rbind(state$kone[i, ], state$ktwo[i, ])
  weight <- state$weight[i, ]
  pdf <- mixture_grid(individual_theta, mu, kappa, clip_mixture_weight(weight))
  bounds <- t(apply(pdf, 1, quantile, probs = c(0.025, 0.975)))
  # Source semantics: PDF of posterior mean parameters, not mean PDF.
  mu_mean <- wrap(rowMeans(mu))
  kappa_mean <- rowMeans(kappa)
  w_mean <- clip_mixture_weight(mean(weight))
  plugin <- as.numeric(mixture_grid(individual_theta,
    matrix(mu_mean, ncol = 1L), matrix(kappa_mean, ncol = 1L), w_mean))
  kernel <- kde_basis(individual_theta, y, bw = 15)
  boot <- matrix(0, nrow = length(individual_theta), ncol = n_boot_kde)
  set.seed(123L + 10000L + i - 1L)
  for (b in seq_len(n_boot_kde)) {
    boot[, b] <- rowMeans(kernel[, sample.int(length(y), length(y), replace = TRUE), drop = FALSE])
  }
  kde <- t(apply(boot, 1, quantile, probs = c(0.025, 0.5, 0.975)))
  curve <- data.frame(individual = i, theta = individual_theta,
    model = plugin, model_low = bounds[, 1], model_high = bounds[, 2],
    kde_low = kde[, 1], kde = kde[, 2], kde_high = kde[, 3])
  # Close the periodic curve explicitly with the first grid value.
  end <- curve[1, ]
  end$theta <- pi
  individual_curves[[i]] <- rbind(curve, end)
  inv_grid <- 1 / plugin
  integral <- sum(diff(individual_theta) *
                  (head(inv_grid, -1L) + tail(inv_grid, -1L)) / 2)
  at_data <- as.numeric(mixture_grid(y, matrix(mu_mean, ncol = 1L),
                                     matrix(kappa_mean, ncol = 1L), w_mean))
  individual_ticks[[i]] <- data.frame(individual = i, theta = y,
    residual = -1 / (at_data * (integral + 1e-12)))
  individual_arrows[[i]] <- data.frame(individual = i, theta = mu_mean,
    length = c(w_mean, 1 - w_mean) * A1_safe(kappa_mean))
}
  levels <- seq_len(num_individuals)
  label <- function(x) {
    x$individual <- factor(x$individual, levels = levels,
                           labels = paste("Individual", levels))
    x
  }
  list(curves = label(do.call(rbind, individual_curves)),
       ticks = label(do.call(rbind, individual_ticks)),
       arrows = label(do.call(rbind, individual_arrows)))
}

# The 52 animal sections of the predictive check are independent, so
# they run on a PSOCK cluster (fresh child R processes) instead of
# in the kernel. Workers are seeded per (model, animal), so results do
# not depend on the core count. (The earlier parallelize attempt used
# parallel::mclapply, whose fork is not safe under IRkernel: workers
# inherit the kernel's live ZMQ socket and the fork corrupts it. This
# also re-seeds the PPC draw streams relative to the original
# shared-stream implementation; the values move by Monte Carlo noise
# only, and the in-notebook gates are on the fit, not on these draws.
# The workers are plain top-level functions, which serialize as
# ~50-byte closures over the globalenv; all data flows through the
# explicit parLapply arguments.)
# The 52 animal sections of the predictive check are independent, so
# they run on a PSOCK cluster (fresh child R processes) instead of
# in the kernel. Workers are seeded per (model, animal), so results do
# not depend on the core count. (The earlier parallelize attempt used
# parallel::mclapply, whose fork is not safe under IRkernel: workers
# inherit the kernel's live ZMQ socket and the fork corrupts it. This
# also re-seeds the PPC draw streams relative to the original
# shared-stream implementation; the values move by Monte Carlo noise
# only, and the in-notebook gates are on the fit, not on these draws.
# The workers are plain top-level functions, which serialize as
# ~50-byte closures over the globalenv; all data flows through the
# explicit parLapply arguments.)
ppc_worker <- function(i, state, draws, n, theta_grid, bins, bin_widths,
                       bin_centers, rng_seed) {
  wrap <- function(x) (x + pi) %% (2 * pi) - pi
  circular_kde <- function(theta, y, bw) {
    yc <- circular::circular(wrap(y), units = "radians", template = "none",
                             modulo = "asis")
    z <- circular::circular(theta, units = "radians", template = "none",
                            modulo = "asis")
    circular::density.circular(yc, z = z, bw = bw,
                               kernel = "vonmises")$y
  }
  bin_heading_density <- function(y) {
    b <- findInterval(wrap(y), bins, rightmost.closed = TRUE)
    tabulate(b, nbins = length(bin_centers)) /
      (length(y) * bin_widths)
  }
  draw_vm <- function(n, mu, kappa) {
    if (!n) return(numeric())
    vmu <- circular::circular(mu, units = "radians", template = "none",
                              modulo = "2pi")
    as.numeric(circular::rvonmises(n, vmu, kappa,
      control.circular = list(units = "radians", template = "none",
                              modulo = "2pi")))
  }
  if (is.null(state$weight)) {
    rep_head <- function(draw) {
      wrap(draw_vm(n, state$muone[i, draw], state$kone[i, draw]))
    }
  } else {
    rep_head <- function(draw) {
      first <- rbinom(1, n, state$weight[i, draw])
      wrap(c(draw_vm(first, state$muone[i, draw], state$kone[i, draw]),
             draw_vm(n - first, state$mutwo[i, draw], state$ktwo[i,
               draw])))
    }
  }
  set.seed(rng_seed)
  y_rep <- lapply(draws, rep_head)
  list(hist = vapply(y_rep, bin_heading_density,
                     numeric(length(bin_centers))),
    population = vapply(y_rep, function(y) {
      circular_kde(theta_grid, y, bw = 30)
    }, numeric(length(theta_grid))),
    individual = vapply(y_rep, function(y) {
      circular_kde(theta_grid, y, bw = 15)
    }, numeric(length(theta_grid))))
}

ppc_task <- function(i, state, draws, n_all, theta_grid, bins,
                     bin_widths, bin_centers, rng_tag) {
  ppc_worker(i, state, draws, n_all[i], theta_grid, bins, bin_widths,
             bin_centers, rng_tag + i)
}

predictive_checks <- function(state, model, model_tag) {
  draws <- sample.int(nrow(state$draws), N_PREDICTIVE)
  n_all <- vapply(rows_by_animal, length, integer(1L))
  ppc_cores <- min(8L, max(2L, parallel::detectCores()))
  cl <- parallel::makePSOCKcluster(ppc_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(cl, quote({
    .libPaths(c(Filter(nzchar, c(Sys.getenv("R_LIBS_USER"),
      Sys.getenv("R_LIBS_SITE"))), .libPaths()))
  }))
  parallel::clusterExport(cl, c("ppc_worker", "ppc_task"),
                          envir = globalenv())
  individuals <- parallel::parLapply(
    cl, seq_len(num_individuals), ppc_task,
    state = state, draws = draws, n_all = n_all, theta_grid = theta_grid,
    bins = bins, bin_widths = bin_width, bin_centers = bin_centers,
    rng_tag = SEED + 90000L + 1000L * model_tag)
  pop_hist <- matrix(0, nrow = length(bin_centers), ncol = length(draws))
  pop_kde <- matrix(0, nrow = length(theta_grid), ncol = length(draws))
  for (ind in individuals) {
    pop_hist <- pop_hist + ind$hist / num_individuals
    pop_kde <- pop_kde + ind$population / num_individuals
  }
  stopifnot(all(abs(colSums(pop_hist) * bin_width - 1) < 1e-10))
  population <- rbind(
    transform(density_frame(bin_centers, pop_hist), estimate = "Histogram"),
    transform(density_frame(theta_grid, pop_kde), estimate = "Von Mises KDE"))
  individual <- do.call(rbind, lapply(seq_along(individuals), function(i) {
    transform(density_frame(theta_grid, individuals[[i]]$individual),
              individual = paste("Individual", i))
  }))
  list(population = transform(population, model = model),
       individual = transform(individual, model = model))
}

hypothesis_arrays <- function(state, clip = FALSE) {
  forward <- mixture_at(0, state, clip)
  left <- mixture_at(-pi / 2, state, clip)
  right <- mixture_at(pi / 2, state, clip)
  list(forward = forward, perpendicular = log_add(left, right) - log(2))
}

interval_row <- function(x, mode = c("density", "elpd", "alignment")) {
  mode <- match.arg(mode)
  if (mode == "elpd") {
    center <- mean(x)
    error <- 1.96 * sd(x) / sqrt(length(x))
    return(c(center = center, low = center - error, high = center + error))
  }
  bounds <- hdi(x)
  if (mode == "density") {
    return(c(center = exp(median(x)), low = exp(unname(bounds[1])), high = exp(unname(bounds[2]))))
  }
  c(center = mean(x), low = unname(bounds[1]), high = unname(bounds[2]))
}

make_hyp_table <- function(arrays, mode, pooled = FALSE) {
  f <- arrays$forward
  p <- arrays$perpendicular
  stopifnot(is.matrix(f), identical(dim(f), dim(p)))
  labels <- as.character(seq_len(nrow(f)))
  if (pooled) {
    f <- rbind(f, colMeans(f))
    p <- rbind(p, colMeans(p))
    labels <- c(labels, "Pooled")
  }
  do.call(rbind, lapply(seq_len(nrow(f)), function(i) {
    z <- rbind(interval_row(f[i, ], mode), interval_row(p[i, ], mode))
    data.frame(individual = labels[i], hypothesis = c("Optic flow", "Vec. encoding"),
      center = z[, 1], low = z[, 2], high = z[, 3],
      probability = mean(p[i, ] > f[i, ]), row.names = NULL)
  }))
}

alignment_arrays <- function(state) {
  w <- state$weight
  stopifnot(is.matrix(w), identical(dim(w), dim(state$muone)))
  list(forward = w * cos(state$muone) + (1 - w) * cos(state$mutwo),
       perpendicular = w * cos(2 * (state$muone - pi / 2)) +
         (1 - w) * cos(2 * (state$mutwo - pi / 2)))
}

