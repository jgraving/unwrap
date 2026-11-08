# Beetles (security light): post-fit computation.
#
# Sourced by notebooks_r/beetles.ipynb (no standalone output). R twin of
# the Python notebook's local posthoc code: the circular and interval math,
# the convergence diagnostic gate, the per-draw direction posteriors and
# the data-bootstrap PPC tables. The notebook keeps the fit, prints the
# returned summary table and the seeds, and calls the figure builders from
# beetles_plots.R.

circ_mean <- function(a) atan2(mean(sin(a)), mean(cos(a)))
resultant <- function(a) Mod(mean(exp(1i * a)))
circ_sd <- function(a, m = circ_mean(a)) {
  sqrt(-2 * log(min(1, Mod(mean(exp(1i * (a - m)))))))
}
# Scaled Bessel functions avoid overflow when concentrations are large.
A1_safe <- function(k) {
  out <- besselI(k, 1, expon.scaled = TRUE) /
    besselI(k, 0, expon.scaled = TRUE)
  if (is.matrix(k)) dim(out) <- dim(k)
  out
}
vm_logpdf <- function(y, mu, kappa) {
  kappa * (cos(y - mu) - 1) - log(2 * pi) -
    log(besselI(kappa, 0, expon.scaled = TRUE))
}
hdi <- function(x, prob = 0.95) {
  x <- sort(as.numeric(x))
  stopifnot(length(x) > 1L, all(is.finite(x)))
  gap <- min(length(x) - 1L, floor(prob * length(x)))
  starts <- seq_len(length(x) - gap)
  j <- starts[which.min(x[starts + gap] - x[starts])]
  c(low = x[j], high = x[j + gap])
}
interval_row <- function(x, angular = FALSE, prob = 0.95) {
  m <- if (angular) circ_mean(x) else mean(x)
  # Unwrap angular intervals around their circular mean before drawing;
  # display endpoints can lie just beyond +/-pi to avoid a false long arc.
  z <- if (angular) m + mod_circular(x - m) else x
  c(mean = m, hdi(z, prob))
}
inverse_A1 <- function(r) {
  if (r <= 0) {
    return(0)
  }
  if (r >= 1) {
    return(Inf)
  }
  upper <- max(2, 1 / (1 - r))
  uniroot(function(k) A1_safe(k) - r, c(0, upper), tol = 1e-9)$root
}
# Exact tutorial sampler convention: scalar mu/kappa, then wrap the draw.
vmu <- function(m) circular::circular(
    m, units = "radians", template = "none", modulo = "2pi")
vm_ctrl <- list(units = "radians", template = "none", modulo = "2pi")
simulate_angles <- function(mu, kappa) {
  stopifnot(is.matrix(mu), identical(dim(mu), dim(kappa)))
  replicated <- matrix(0, nrow = nrow(mu), ncol = ncol(mu))
  for (d in seq_len(nrow(mu))) {
    for (i in seq_len(ncol(mu))) {
      replicated[d, i] <- as.numeric(circular::rvonmises(1,
        mu = vmu(mu[d, i]), kappa = kappa[d, i],
        control.circular = vm_ctrl))
    }
  }
  mod_circular(replicated)
}

diagnose_fit <- function(fit, circular_parameters = character()) {
  # NUTS divergences. Not every brms/cmdstanr combination keeps the
  # per-iteration divergence flag in the fit object, so probe for it and
  # report what the fit object actually carries (the R-hat and ESS values
  # below are the primary convergence checks either way).
  divergences <- tryCatch({
    div_flags <- lapply(fit$fit@sim$samples, function(ch) ch[["divergent"]])
    if (any(vapply(div_flags, is.null, logical(1)))) NA_real_
    else sum(unlist(div_flags))
  }, error = function(e) NA_real_)
  if (is.na(divergences)) {
      cat("Divergences: not retained in this fit object (see R-hat/ESS below).\n")
  } else {
      cat("Divergences:", divergences, "\n")
  }
  available <- posterior::variables(posterior::as_draws_array(fit))
  parameters <- available[grepl("^(b_|sd_|r_|kmu_)", available)]
  draws_array <- posterior::subset_draws(posterior::as_draws_array(fit),
                                         variable = parameters)
  out <- do.call(rbind, lapply(parameters, function(p) {
    # Preserve iteration x chain shape: flattening loses between-chain
    # information and would make the R-hat/ESS checks misleading.
    x <- as.matrix(draws_array[, , p])
    angular <- p %in% circular_parameters
    z <- if (angular) mod_circular(x - circ_mean(x)) else x
    stopifnot(all(is.finite(x)))
    data.frame(parameter = p,
      mean = if (angular) circ_mean(x) else mean(x),
      sd = if (angular) circ_sd(x) else sd(as.numeric(x)),
      rhat = if (angular) Rhat_unwrap(x, method = "sine") else posterior::rhat(x),
      ess_bulk = posterior::ess_bulk(z),
      row.names = NULL, check.names = FALSE)
  }))
  out$rhat <- round(out$rhat, 4)
  out$ess_bulk <- round(out$ess_bulk, 1)
  out
}

# Per-draw individual-condition direction posteriors: the mu/kappa draw
# matrices over the per-individual x 2 conditions grid, the
# between-individual alignment and the individual summary table.
posterior_directions <- function(fit, data, condition_labels) {
  ids <- levels(data$Individual)
  prediction_data <- expand.grid(Individual = ids, Light = c(0, 1),
                                 KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  prediction_data$Individual <- factor(prediction_data$Individual, levels = ids)
  prediction_data$condition <- factor(prediction_data$Light,
                                       levels = c(0, 1), labels = condition_labels)
  # Draw x individual-condition matrices. Evaluate each fitted distributional
  # predictor explicitly; no arithmetic averaging across the angular seam.
  mu <- mod_circular(posterior_linpred(fit, newdata = prediction_data,
                                       transform = FALSE, dpar = "mu"))
  kappa <- posterior_linpred(fit, newdata = prediction_data,
                             transform = TRUE, dpar = "kappa")
  stopifnot(identical(dim(mu), dim(kappa)), all(kappa > 0))
  J <- length(ids)
  pre <- seq_len(J)
  post <- J + pre
  group_alignment <- cbind(
    Mod(rowMeans(exp(1i * mu[, pre, drop = FALSE]))),
    Mod(rowMeans(exp(1i * mu[, post, drop = FALSE]))))
  individual_summary <- do.call(rbind, lapply(seq_len(nrow(prediction_data)), function(j) {
    rbind(data.frame(id = as.character(prediction_data$Individual[j]),
      condition = prediction_data$condition[j], quantity = "mean direction (rad)",
      as.list(interval_row(mu[, j], angular = TRUE))),
      data.frame(id = as.character(prediction_data$Individual[j]),
      condition = prediction_data$condition[j], quantity = "within-individual R(kappa)",
      as.list(interval_row(A1_safe(kappa[, j])))))
  }))
  list(ids = ids, prediction_data = prediction_data, mu = mu,
    kappa = kappa, J = J, pre = pre, post = post,
    group_alignment = group_alignment,
    individual_summary = individual_summary)
}

# Supplemental PPC: data bootstraps (group and per-individual) and the
# fitted population posteriors that the twin figure compares against.
ppc_bootstrap_stats <- function(fit, prediction_data, ys, J,
  B_group, B_indiv, condition_labels, group_alignment) {
  bootstrap_pair <- function(y, B) {
    z <- replicate(B, mean(exp(1i * y[sample.int(length(y), length(y), replace = TRUE)])))
    list(mu = Arg(z), R = Mod(z))
  }
  boot_individual <- lapply(ys, bootstrap_pair, B = B_indiv)
  boot_group_mu <- boot_group_R <- matrix(0, nrow = B_group, ncol = 2)
  for (b in seq_len(B_group)) {
    selected <- sample.int(J, J, replace = TRUE)
    for (c in 1:2) {
      jj <- selected + (c - 1L) * J
      directions <- vapply(jj, function(j) {
        y <- ys[[j]]
        circ_mean(y[sample.int(length(y), length(y), replace = TRUE)])
      }, numeric(1))
      boot_group_mu[b, c] <- circ_mean(directions)
      boot_group_R[b, c] <- resultant(directions)
    }
  }
  fixed <- as.data.frame(posterior::as_draws_df(fit))
  population_mu <- cbind(mod_circular(fixed$b_beta0_Intercept),
    mod_circular(fixed$b_beta0_Intercept + fixed$b_beta0_Light))
  list(boot_individual = boot_individual, boot_group_mu = boot_group_mu,
    boot_group_R = boot_group_R, population_mu = population_mu)
}

# Conventional replicated-data check: each individual's own conditional
# likelihood resampled on the observed sample size, summarised as the
# replicated mean resultant length next to the observed one.
replicated_stats <- function(mu, kappa, ys, draw_ids, prediction_data) {
  rep_stats <- do.call(rbind, lapply(seq_along(ys), function(j) {
    rep_mu <- matrix(mu[draw_ids, j], nrow = length(draw_ids), ncol = length(ys[[j]]))
    rep_kappa <- matrix(kappa[draw_ids, j], nrow = length(draw_ids), ncol = length(ys[[j]]))
    yrep <- simulate_angles(rep_mu, rep_kappa)
    data.frame(id = as.character(prediction_data$Individual[j]),
      condition = prediction_data$condition[j],
      as.list(interval_row(apply(yrep, 1, resultant))),
      observed = resultant(ys[[j]]))
  }))
}
