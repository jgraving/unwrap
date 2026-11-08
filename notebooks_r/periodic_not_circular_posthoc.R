# Periodic (not circular) beta: post-fit computation.
#
# Sourced by notebooks_r/periodic_not_circular.ipynb (no standalone
# output). R twin of the Python notebook's local helpers: the circular
# and interval math, the von Mises draw wrapper and the convergence
# diagnostic gate (R-hat/ESS with stopifnot thresholds; the notebook
# prints the returned summary table itself).

wrap <- function(x) atan2(sin(x), cos(x))

circ_mean <- function(x) atan2(mean(sin(x)), mean(cos(x)))

resultant <- function(x) Mod(mean(exp(1i * x)))

softplus_stable <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

log1p_exp <- softplus_stable

A1_safe <- function(k) {
  z <- besselI(k, 1, expon.scaled = TRUE) / besselI(k, 0, expon.scaled = TRUE)
  if (!is.null(dim(k))) dim(z) <- dim(k)
  z
}

vm_logpdf <- function(y, mu, kappa) {
  kappa * (cos(y - mu) - 1) - log(2 * pi) - log(besselI(kappa, 0, expon.scaled = TRUE))
}

hdi <- function(x, prob = 0.95) {
  x <- sort(as.numeric(x))
  stopifnot(length(x) > 1L, all(is.finite(x)))
  gap <- min(length(x) - 1L, floor(prob * length(x)))
  ii <- seq_len(length(x) - gap)
  j <- ii[which.min(x[ii + gap] - x[ii])]
  c(low = x[j], high = x[j + gap])
}

interval <- function(x, angular = FALSE, prob = 0.95) {
  m <- if (angular) circ_mean(x) else mean(x)
  z <- if (angular) m + wrap(x - m) else x
  c(mean = m, hdi(z, prob))
}

vmu <- function(m) circular::circular(
    m, units = "radians", template = "none", modulo = "2pi")

vm_ctrl <- list(units = "radians", template = "none", modulo = "2pi")

draw_vm <- function(mu, kappa) {
  wrap(as.numeric(circular::rvonmises(1, mu = vmu(mu), kappa = kappa,
                                     control.circular = vm_ctrl)))
}

draw_matrix <- function(fit, parameters) {
  as.data.frame(posterior::as_draws_df(posterior::subset_draws(
    posterior::as_draws_array(fit), variable = parameters)))
}

diagnose_fit <- function(fit, angular = character(), parameters = NULL) {
  # Probe exactly the flags retained by this brms/backend combination.
  divergences <- tryCatch({
    samples <- fit$fit@sim$samples
    if (is.null(samples) || !length(samples)) {
      NA_real_
    } else {
      div_flags <- lapply(samples, function(ch) ch[["divergent"]])
      if (any(vapply(div_flags, is.null, logical(1)))) NA_real_
      else sum(unlist(div_flags))
    }
  }, error = function(e) NA_real_)
  if (is.na(divergences)) {
    cat("Divergences: not retained in this fit object (see R-hat/ESS below).\n")
  } else {
    cat("Divergences:", divergences, "\n")
    stopifnot(divergences == 0)
  }
  a <- posterior::as_draws_array(fit)
  if (is.null(parameters)) {
    parameters <- posterior::variables(a)
    parameters <- parameters[grepl("^(b_|sd_|r_|kmu_)", parameters)]
  }
  stopifnot(length(parameters) > 0L)
  summary_frame <- do.call(rbind, lapply(parameters, function(p) {
    x <- as.matrix(a[, , p])
    # as_draws_array keeps the draws_array class, whose [ method does not
    # drop the parameter slice, so x is retained as one long column
    # (draws x chains flattened); the rhat/ESS calls below are written
    # for that layout (as in the P1 notebooks).
    stopifnot(all(is.finite(x)))
    is_angle <- p %in% angular
    z <- if (is_angle) wrap(x - circ_mean(x)) else x
    data.frame(parameter = p,
      mean = if (is_angle) circ_mean(x) else mean(x),
      sd = if (is_angle) sqrt(-2 * log(min(1, resultant(x)))) else sd(as.numeric(x)),
      rhat = if (is_angle) Rhat_unwrap(x, method = "sine") else posterior::rhat(x),
      ess_bulk = posterior::ess_bulk(z), row.names = NULL)
  }))
  stopifnot(all(is.finite(summary_frame$rhat)), all(summary_frame$rhat < 1.01),
    all(is.finite(summary_frame$ess_bulk)), all(summary_frame$ess_bulk > 400))
  summary_frame
}
