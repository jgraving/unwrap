# Locust density/order: post-fit computation.
#
# Sourced by notebooks_r/locust_density_order.ipynb (no standalone
# output). R twin of notebooks/locust_density_posthoc.py: the PSIS-LOO
# helper, the screening alignment stats, the five-model sweep
# delta-ELPD table, the final-fit grid/polar/elpd0/coefficient
# summaries and the bootstrap + posterior predictive simulation
# helpers. The notebook keeps the data, the preprocessing, the fits
# (and their fit-cache names), the seeds and the printed diagnostics,
# and calls the figure builders from
# locust_density_order_plots.R.


explicit_loo <- function(fit, newdata, response = "y") {
  draw_index <- posterior::as_draws_df(fit)
  log_lik <- matrix(0, nrow = nrow(draw_index), ncol = nrow(newdata))
  # Keep only the full log-likelihood matrix, evaluating response parameters
  # in observation blocks to limit memory in the 18,641-step final model.
  for (start in seq.int(1L, nrow(newdata), by = 256L)) {
    ii <- seq.int(start, min(start + 255L, nrow(newdata)))
    block <- newdata[ii, , drop = FALSE]
    mu_ll <- mod_circular(posterior_linpred(fit, newdata = block,
                                            dpar = "mu", transform = FALSE))
    k_ll <- posterior_linpred(fit, newdata = block,
                              dpar = "kappa", transform = TRUE)
    log_lik[, ii] <- vm_logpdf(matrix(block[[response]], nrow = nrow(mu_ll),
      ncol = nrow(block), byrow = TRUE), mu_ll, k_ll)
  }
  # Column scaling leaves relative efficiencies unchanged and avoids
  # underflow when converting log likelihoods to likelihoods.
  likelihood <- exp(sweep(log_lik, 2, apply(log_lik, 2, max), "-"))
  r_eff <- loo::relative_eff(likelihood, chain_id = draw_index[[".chain"]])
  loo::loo(log_lik, r_eff = r_eff)
}

candidate_rhs <- c("Intercept-only" = "1", "Density" = "t_density",
  "Order" = "t_order", "Density + Order" = "t_density + t_order",
  "Density x Order" = "t_density * t_order")

sweep_delta_table <- function(sweep_loo) {
  elpd <- vapply(sweep_loo, function(l) l$estimates["elpd_loo", "Estimate"], numeric(1))
  best_name <- names(which.max(elpd))
  best_pointwise <- sweep_loo[[best_name]]$pointwise[, "elpd_loo"]
  sweep_df <- do.call(rbind, lapply(names(sweep_loo), function(name) {
    l <- sweep_loo[[name]]
    delta <- l$pointwise[, "elpd_loo"] - best_pointwise
    data.frame(model = name, elpd = elpd[[name]],
      delta_elpd = sum(delta), se_delta_elpd = sqrt(length(delta) * var(delta)),
      pareto_k_max = max(l$diagnostics$pareto_k), best_model = best_name)
  }))
  sweep_df[order(sweep_df$delta_elpd, decreasing = TRUE), ]
}

screen_alignment_stats <- function(screen_fit) {
  screen_grid <- expand.grid(order = seq(0, 1, length.out = 50), density = density_levels)
  screen_grid$t_order <- 2 * screen_grid$order - 1
  screen_grid$t_density <- 2 * log2(screen_grid$density) / log2(64) - 1
  screen_mu <- mod_circular(posterior_linpred(screen_fit, newdata = screen_grid,
                                              dpar = "mu", transform = FALSE))
  screen_k <- posterior_linpred(screen_fit, newdata = screen_grid,
                                dpar = "kappa", transform = TRUE)
  screen_alignment <- A1_safe(screen_k) * cos(screen_mu)
  do.call(rbind, lapply(seq_len(nrow(screen_grid)), function(j) {
    data.frame(screen_grid[j, c("order", "density")],
      as.list(interval_row(screen_alignment[, j], prob = 0.68)))
  }))
}

summarize_log_density <- function(m, k, newdata) {
  do.call(rbind, lapply(seq_len(nrow(newdata)), function(j) {
    lp <- vm_logpdf(0, m[, j], k[, j])
    # Python summarizes in log space, then exponentiates the interval.
    s <- interval_row(lp, prob = 0.9)
    data.frame(newdata[j, c("order", "density")], as.list(exp(s)))
  }))
}

final_grid_stats <- function(fit, data) {
  grid <- expand.grid(order = seq(0, 1, length.out = 50), density = density_levels)
  grid$t_order <- 2 * grid$order - 1
  grid$t_density <- 2 * log2(grid$density) / log2(64) - 1
  grid$animal_id <- data$animal_id[1]
  mu_grid <- mod_circular(posterior_linpred(fit, newdata = grid,
    dpar = "mu", transform = FALSE, re_formula = NA))
  k_grid <- posterior_linpred(fit, newdata = grid,
    dpar = "kappa", transform = TRUE, re_formula = NA)
  alignment <- A1_safe(k_grid) * cos(mu_grid)
  alignment_summary <- do.call(rbind, lapply(seq_len(nrow(grid)), function(j) {
    data.frame(grid[j, c("order", "density")], as.list(interval_row(alignment[, j], prob = 0.9)))
  }))
  density_grid <- expand.grid(density = 2^seq(0, 6, length.out = 50), order = order_levels)
  density_grid$t_order <- 2 * density_grid$order - 1
  density_grid$t_density <- 2 * log2(density_grid$density) / log2(64) - 1
  density_grid$animal_id <- data$animal_id[1]
  mu_density <- mod_circular(posterior_linpred(fit, newdata = density_grid,
    dpar = "mu", transform = FALSE, re_formula = NA))
  k_density <- posterior_linpred(fit, newdata = density_grid,
    dpar = "kappa", transform = TRUE, re_formula = NA)
  density_summary <- do.call(rbind, lapply(seq_len(nrow(density_grid)), function(j) {
    data.frame(density_grid[j, c("order", "density")],
      as.list(interval_row(A1_safe(k_density[, j]) * cos(mu_density[, j]), prob = 0.9)))
  }))
  list(grid = grid, mu_grid = mu_grid, k_grid = k_grid, density_grid = density_grid,
  mu_density = mu_density, k_density = k_density, alignment_summary = alignment_summary,
  density_summary = density_summary)
}

polar_pdf_stats <- function(fit, data) {
  cell_grid <- expand.grid(order = order_levels, density = density_levels)
  cell_grid$t_order <- 2 * cell_grid$order - 1
  cell_grid$t_density <- 2 * log2(cell_grid$density) / log2(64) - 1
  cell_grid$animal_id <- data$animal_id[1]
  m_cell <- mod_circular(posterior_linpred(fit, newdata = cell_grid,
    dpar = "mu", transform = FALSE, re_formula = NA))
  k_cell <- posterior_linpred(fit, newdata = cell_grid,
    dpar = "kappa", transform = TRUE, re_formula = NA)
  theta <- seq(-pi, pi, length.out = 900)
  pdf_summary <- do.call(rbind, lapply(seq_len(nrow(cell_grid)), function(j) {
    do.call(rbind, lapply(theta, function(t) {
      p <- exp(vm_logpdf(t, m_cell[, j], k_cell[, j]))
      data.frame(order = cell_grid$order[j], density = cell_grid$density[j],
        theta = t %% (2 * pi), median = median(p), as.list(hdi(p, 0.9)))
    }))
  }))
  list(theta = theta, pdf_summary = pdf_summary)
}

coefficient_summary <- function(fit, trace_parameters) {
  coefficient_draws <- as.data.frame(posterior::as_draws_df(fit))
  mu_pars <- trace_parameters[!grepl("b_kappa_", trace_parameters)]
  k_pars <- paste0("b_kappa_", sub("^b_", "", mu_pars))
  stopifnot(length(mu_pars) == 4L, all(grepl("b_kappa_", k_pars)))
  cs <- do.call(rbind, lapply(seq_along(mu_pars), function(j) {
    rbind(data.frame(term = j, comp = "mu",
      as.list(interval_row(coefficient_draws[[mu_pars[j]]],
        angular = mu_pars[j] == "b_Intercept"))),
      data.frame(term = j, comp = "kappa",
        as.list(interval_row(coefficient_draws[[k_pars[j]]]))))
  }))
  cs$ypos <- cs$term + ifelse(cs$comp == "kappa", 0.16, -0.16)
  cs
}

individual_ppc_simulation <- function(fit, data) {
  draw_ids <- sample.int(posterior::ndraws(posterior::as_draws_array(fit)), N_REPLICATES)
  # Subset before allocating draw x observation matrices (~18,641 columns).
  mu_obs <- mod_circular(posterior_linpred(fit, dpar = "mu", transform = FALSE, draw_ids = draw_ids))
  k_obs <- posterior_linpred(fit, dpar = "kappa", transform = TRUE, draw_ids = draw_ids)
  yrep <- simulate_angles(mu_obs, k_obs)
  animal_rows <- split(seq_len(nrow(data)), data$animal_id)
  animal_metadata <- data[vapply(animal_rows, function(ii) ii[1], integer(1)),
                          c("animal_id", "order", "density")]
  individual_ppc <- do.call(rbind, lapply(seq_along(animal_rows), function(j) {
    ii <- animal_rows[[j]]
    y <- data$y[ii]
    empirical <- replicate(B_INDIV, mean(cos(y[sample.int(length(y), length(y), replace = TRUE)])))
    predictive <- rowMeans(cos(yrep[, ii, drop = FALSE]))
    rbind(data.frame(animal_metadata[j, ], source = "data bootstrap",
      as.list(interval_row(empirical))),
      data.frame(animal_metadata[j, ], source = "posterior predictive",
      as.list(interval_row(predictive))))
  }))
  list(yrep = yrep, animal_rows = animal_rows, animal_metadata = animal_metadata,
  individual_ppc = individual_ppc)
}

cell_ppc_stats <- function(animal_metadata, data, animal_rows, yrep) {
  cell_rows <- split(seq_len(nrow(animal_metadata)),
    interaction(animal_metadata$order, animal_metadata$density, drop = TRUE))
  do.call(rbind, lapply(cell_rows, function(animals) {
    empirical <- replicate(B_DATA, {
      selected <- animals[sample.int(length(animals), length(animals), replace = TRUE)]
      mean(vapply(selected, function(j) {
        y <- data$y[animal_rows[[j]]]
        mean(cos(y[sample.int(length(y), length(y), replace = TRUE)]))
      }, numeric(1)))
    })
    predictive <- replicate(B_MODEL, {
      d <- sample.int(nrow(yrep), 1)
      selected <- animals[sample.int(length(animals), length(animals), replace = TRUE)]
      mean(vapply(selected, function(j) {
        y <- yrep[d, animal_rows[[j]]]
        mean(cos(y[sample.int(length(y), length(y), replace = TRUE)]))
      }, numeric(1)))
    })
    g <- animal_metadata[animals[1], c("order", "density")]
    rbind(data.frame(g, source = "data bootstrap", n_an = length(animals),
      as.list(interval_row(empirical))),
      data.frame(g, source = "posterior predictive", n_an = length(animals),
      as.list(interval_row(predictive))))
  }))
}
