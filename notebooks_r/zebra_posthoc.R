# Zebra escape: post-fit computation.
#
# Sourced by notebooks_r/zebra_escape.ipynb (no standalone output).
# R twin of notebooks/zebra_posthoc.py: the per-draw direction posteriors
# and the data-bootstrap / replicated-data PPC statistics. The circular
# and interval math, the convergence diagnostic gate and the shared
# sampling helpers come from R/unwrap_helpers.R (the notebook sources it
# before this script); only the zebra-specific post-fit analysis lives
# here. The notebook keeps the fit, prints the returned summary table and
# the seeds, and calls the figure builders from zebra_plots.R.

# Per-draw track-condition direction posteriors: the mu/kappa draw
# matrices over the 18 tracks x 2 conditions grid, the between-individual
# alignment and the individual summary table (twin of zebra_posthoc.py's
# posterior_directions/summary_stats input computation).
posterior_directions <- function(fit, data, condition_labels) {
  ids <- levels(data$Track)
  prediction_data <- expand.grid(Track = ids, PostStimulus = c(0, 1),
                                 KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  prediction_data$Track <- factor(prediction_data$Track, levels = ids)
  prediction_data$condition <- factor(prediction_data$PostStimulus,
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
    rbind(data.frame(id = as.character(prediction_data$Track[j]),
      condition = prediction_data$condition[j], quantity = "mean direction (rad)",
      as.list(interval_row(mu[, j], angular = TRUE))),
      data.frame(id = as.character(prediction_data$Track[j]),
      condition = prediction_data$condition[j], quantity = "within-individual R(kappa)",
      as.list(interval_row(A1_safe(kappa[, j])))))
  }))
  list(prediction_data = prediction_data, mu = mu, kappa = kappa,
    J = J, pre = pre, post = post,
    group_alignment = group_alignment,
    individual_summary = individual_summary)
}

# Supplemental PPC: data bootstraps (group and individual) against the
# fitted population and per-individual posteriors, assembled as the two
# point-range frame layouts the plot functions draw.
ppc_bootstrap_stats <- function(fit, prediction_data, ys, J,
  B_group, B_indiv, condition_labels, group_alignment, individual_summary) {
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
    mod_circular(fixed$b_beta0_Intercept + fixed$b_beta0_PostStimulus))
  ppc_group <- do.call(rbind, lapply(1:2, function(c) {
    rbind(
      data.frame(condition = condition_labels[c], source = "data bootstrap",
        quantity = "group mean direction (rad)",
        as.list(interval_row(boot_group_mu[, c], TRUE))),
      data.frame(condition = condition_labels[c], source = "posterior",
        quantity = "group mean direction (rad)",
        as.list(interval_row(population_mu[, c], TRUE))),
      data.frame(condition = condition_labels[c], source = "data bootstrap",
        quantity = "between-individual alignment",
        as.list(interval_row(boot_group_R[, c]))),
      data.frame(condition = condition_labels[c], source = "posterior",
        quantity = "between-individual alignment",
        as.list(interval_row(group_alignment[, c]))))
  }))
  ppc_individual <- do.call(rbind, lapply(seq_along(ys), function(j) {
    rbind(data.frame(id = as.character(prediction_data$Track[j]),
      condition = prediction_data$condition[j], source = "data bootstrap",
      quantity = "mean direction (rad)",
      as.list(interval_row(boot_individual[[j]]$mu, TRUE))),
    data.frame(id = as.character(prediction_data$Track[j]),
      condition = prediction_data$condition[j], source = "data bootstrap",
      quantity = "within-individual R(kappa)",
      as.list(interval_row(boot_individual[[j]]$R))),
    transform(individual_summary[individual_summary$id ==
      as.character(prediction_data$Track[j]) &
      individual_summary$condition == prediction_data$condition[j], ],
      source = "posterior"))
  }))
  list(ppc_group = ppc_group, ppc_individual = ppc_individual)
}

# Conventional replicated-data check: each track's own conditional
replicated_stats <- function(mu, kappa, ys, draw_ids, prediction_data) {
  rep_stats <- do.call(rbind, lapply(seq_along(ys), function(j) {
    rep_mu <- matrix(mu[draw_ids, j], nrow = length(draw_ids), ncol = length(ys[[j]]))
    rep_kappa <- matrix(kappa[draw_ids, j], nrow = length(draw_ids), ncol = length(ys[[j]]))
    yrep <- simulate_angles(rep_mu, rep_kappa)
    data.frame(id = as.character(prediction_data$Track[j]),
      condition = prediction_data$condition[j],
      as.list(interval_row(apply(yrep, 1, resultant))),
      observed = resultant(ys[[j]]))
  }))
}
