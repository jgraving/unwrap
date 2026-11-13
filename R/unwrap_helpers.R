# Shared R helpers and model components for the unwrap R examples.
#
# Every R notebook sources this file in its Setup cell, right after
# R/unwrap_family.R and R/unwrap_functions.R and before its own
# notebooks_r/<stem>_posthoc.R / <stem>_plots.R scripts, so a stem script
# can still specialise any definition here for a single example.
#
# Layers:
#   1. Data and figure plumbing: ensure_data, show_figure, save_figure.
#   2. Circular statistics: wrap, mod_circular, the unwrap helpers,
#      circ_mean, circ_sd, resultant, A1_safe, vm_logpdf, vmu,
#      simulate_angles, hdi, interval_row, inverse_A1, softplus,
#      Rhat_unwrap.
#   3. Reusable model components: the vector von Mises Stan accumulator
#      (stan_unwrap_vector_lpdf), the von Mises / lognormal
#      posterior-target prior builders (vm_target_prior,
#      vm_vector_target_prior, lognormal_target_prior), the cached brms
#      fit runner (unwrap_brms_fit over the unwrap_*_fit_cache
#      primitives), and the shared convergence diagnostic gate
#      (diagnose_fit).
#
# Provenance: each definition below is the copy the examples carried in
# duplicate. show_figure sat in all eight notebooks, save_figure in seven,
# ensure_data in four, and the circular statistics in the notebook setup
# cells and again in the per-stem *_posthoc.R scripts. Where those copies
# had drifted in cosmetics only (argument names, line layout, a missing
# dimension-preserving guard), the copy here is the most general body and
# the stems keep their existing call sites. The fit-cache primitives now
# live only here: the byte-identical duplicate in R/unwrap_functions.R has
# been pruned, so this file is their single source.

# --- Data and figure plumbing -------------------------------------------------

# The notebook defines data_root in its Setup before this is first called
# (each example keeps its own UNWRAP_DATA_ROOT override there).
ensure_data <- function(subdir) {
  target <- file.path(data_root, subdir)
  if (!dir.exists(target)) {
    if (!dir.exists(data_root)) {
      status <- system2("git", c("clone", "--depth", "1",
        "--filter=blob:none", "--sparse",
        "https://github.com/jgraving/unwrap-data.git", shQuote(data_root)))
      if (status != 0L) stop("Could not fetch the published study data.")
    }
    # Sparse checkout downloads only this study; later notebooks may add
    # their own subdirectory to the same checkout without replacing it.
    status <- system2("git", c("-C", shQuote(data_root),
                              "sparse-checkout", "add", subdir))
    if (status != 0L) stop("Could not materialize the study subdirectory.")
  }
  stopifnot(dir.exists(target))
  target
}

# The notebook defines FIG_DIR in its Setup before this is first called.
# Without X11/cairo (most server R installs) the default device cannot
# rasterise, so each figure renders to a PNG through ragg; a vector PDF is
# saved alongside as the portable copy. In a Jupyter (ir) kernel the PNG
# is displayed inline through IRdisplay, exactly as the notebook copies
# do. Outside Jupyter there is no display channel, so the plot object is
# returned instead and a native-R renderer (RStudio, Quarto) draws it -
# the twin of the fallback tools/ipynb_to_qmd.py has always spliced into
# the generated .qmd files.
show_figure <- function(plot, name, width, height, res = 110) {
  file <- file.path(FIG_DIR, paste0(name, ".png"))
  ragg::agg_png(file, width = width, height = height, units = "in",
                res = res)
  print(plot)
  invisible(dev.off())
  tryCatch({
    IRdisplay::display_png(file = file)
    invisible(NULL)
  }, error = function(e) plot)
}

# All figure call sites in the current examples pass width and height
# explicitly; the defaults only serve a standalone script calling bare.
save_figure <- function(plot, name, width = 8, height = 4) {
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), plot,
         width = width, height = height)
  show_figure(plot, name, width, height)
}

# --- Circular statistics -------------------------------------------------------

# Canonical modulo wrap (-pi, pi]; byte-identical to the copies in the
# bats, locust-mixture and policy examples, and value-identical to the
# older atan2(sin(x), cos(x)) form on all reals (wrap PR #126).
wrap <- function(x) (x + pi) %% (2 * pi) - pi

# the radian equivalent
mod_circular = function(x)
{
  atan2(y = sin(x),
        x = cos(x))
}

# the unwrapped (no discontinuities)
unwrap_circular = function(x)
{
  mux = mean.circular(x = circular(x = x, template = 'none'))
  centx = atan2(y = sin(x - mux),
                x = cos(x  - mux))
  unwrx = centx + mux
}

# the unwrapped (no discontinuities)
unwrap_circular_median = function(x)
{
  mux = median.circular(x = circular(x = x, template = 'none'))[[1]]
  centx = atan2(y = sin(x - mux),
                x = cos(x  - mux))
  unwrx = centx + mux
}

#calculate rhat for circular variables with extreme ranges
Rhat_unwrap = function(x,
                       method = 'mean')
{
  switch(EXPR = method,
         mean = rhat(unwrap_circular(x)),
         median = rhat(unwrap_circular_median(x)),
         sine = max(rhat(sin(x)), rhat(cos(x))),
         rhat(unwrap_circular(x)) )

}

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

# Stable inverse softplus link. Mathematically identical to the legacy
# log(exp(x) + 1) in R/unwrap_functions.R but finite for large |x|; the
# stems that source this file pick up this copy.
softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

# --- Reusable model components -------------------------------------------------

# Vector-valued von Mises log pdf (Stan): applies a scalar von Mises
# density to every element of a direction-coefficient vector and sums
# the log densities. This is the accumulator the `target += ...` prior
# terms below apply to a whole coefficient vector. It composes with the
# family's unwrap_von_mises_lpdf from R/unwrap_family.R, so supply it to
# brm() together with stan_unwrap_fun. The Stan text is the exact block
# the examples carried inline.
stan_unwrap_vector_lpdf <- stanvar(scode = "
  real unwrap_von_mises_vect_lpdf(vector y, real mu, real kappa) {
    real lp = 0;
    for (i in 1:num_elements(y)) {
      lp += unwrap_von_mises_lpdf(y[i] | mu, kappa);
    }
    return lp;
  }
", block = "functions")

# Von Mises posterior-target prior for one direction coefficient. The
# coefficient is an unconstrained real with a periodic von Mises prior:
# its log density is added straight to the target instead of using a
# bounded uniform (uniform boxes pin NUTS chains against the +-pi walls).
vm_target_prior <- function(coef, center, kappa) {
  set_prior(paste0("target += unwrap_von_mises_lpdf(", coef,
                   " | ", center, ", ", kappa, ")"), check = FALSE)
}

# Same, applied to a whole direction-coefficient vector (it needs the
# stan_unwrap_vector_lpdf Stan function supplied via stanvars).
vm_vector_target_prior <- function(coef, center, kappa) {
  set_prior(paste0("target += unwrap_von_mises_vect_lpdf(", coef,
                   " | ", center, ", ", kappa, ")"), check = FALSE)
}

# LogNormal posterior-target prior for a declared positive parameter.
lognormal_target_prior <- function(param, mu, sigma) {
  set_prior(paste0("target += lognormal_lpdf(", param,
                   " | ", mu, ", ", sigma, ")"), check = FALSE)
}

# Development fit cache.
# When UNWRAP_R_FIT_CACHE=1, a notebook's brms fits are reused across
# local runs instead of resampled: before fitting, a cell asks
# unwrap_load_fit_cache() for a saved snapshot; when none is valid it
# fits as usual and stores the result with unwrap_save_fit_cache().
# A snapshot is only reused when name, seed, iter, and warmup all match
# (and the fixed-effects size still matches), so a stale cache can
# never sneak into a run.
# name is caller-provided, so it must identify the MODEL as well as its
# budget: a model change (priors/formula/family) under an unchanged name
# would silently reuse the old snapshot. Bump the name (e.g. a version
# suffix like "fit_v2") whenever the model changes.
# With UNWRAP_R_FIT_CACHE unset (the default,
# and how committed notebooks are validated) these helpers are silent
# no-ops and nothing is read from or written to disk.
unwrap_fit_cache_enabled <- function() {
  identical(Sys.getenv("UNWRAP_R_FIT_CACHE"), "1")
}

unwrap_fit_cache_path <- function(name) {
  if (!unwrap_fit_cache_enabled()) {
    return(NULL)
  }
  dir <- Sys.getenv("UNWRAP_R_FIT_CACHE_DIR")
  if (!nzchar(dir)) {
    # Stable across R sessions (tempdir() is per-session and would be
    # wiped on kernel exit, defeating the cache). Override with
    # UNWRAP_R_FIT_CACHE_DIR for a shared or scratch location.
    home <- Sys.getenv("HOME")
    dir <- if (nzchar(home)) file.path(home, ".unwrap_fit_cache")
    else file.path(tempdir(), "unwrap_fit_cache")
  }
  file.path(dir, paste0(name, ".rds"))
}

unwrap_load_fit_cache <- function(name, seed, iter, warmup) {
  path <- unwrap_fit_cache_path(name)
  if (is.null(path) || !file.exists(path)) {
    return(NULL)
  }
  entry <- tryCatch(readRDS(path), error = function(e) NULL)
  meta <- if (is.list(entry)) attr(entry, "unwrap_fit_cache_meta") else NULL
  ok <- is.list(meta) &&
    identical(as.character(meta$name), name) &&
    identical(as.numeric(meta$seed), as.numeric(seed)) &&
    identical(as.numeric(meta$iter), as.numeric(iter)) &&
    identical(as.numeric(meta$warmup), as.numeric(warmup))
  if (ok && !is.na(meta$n_fixed)) {
    # `==` (not identical): length() is integer, the stored count numeric.
    ok <- tryCatch(
      length(fixef(entry)) == as.numeric(meta$n_fixed),
      error = function(e) FALSE)
  }
  if (!ok) {
    message("Unwrap fit cache: '", name, "' is stale; refitting.")
    return(NULL)
  }
  message("Unwrap fit cache: reusing '", name, "' (seed ", meta$seed,
          ", iter ", meta$iter, ", warmup ", meta$warmup, ").")
  entry
}

unwrap_save_fit_cache <- function(name, fit, seed, iter, warmup) {
  path <- unwrap_fit_cache_path(name)
  if (is.null(path)) {
    return(invisible(NULL))
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  n_fixed <- tryCatch(length(fixef(fit)), error = function(e) NA_real_)
  meta <- list(name = name, seed = seed, iter = iter, warmup = warmup,
               n_fixed = n_fixed,
               saved_utc = format(Sys.time(), usetz = TRUE))
  if (is.list(fit)) {
    attr(fit, "unwrap_fit_cache_meta") <- meta
  }
  tmp <- paste0(path, ".tmp")
  saveRDS(fit, tmp)
  file.rename(tmp, path)
  invisible(NULL)
}

# Standard load-or-fit-for this example: reuse a valid cached snapshot
# when the fit cache is enabled; otherwise run the (quoted) brm() call
# and store its result. fit_call is only evaluated on a cache miss.
unwrap_brms_fit <- function(name, seed, iter, warmup, fit_call) {
  fit <- unwrap_load_fit_cache(name, seed, iter, warmup)
  if (is.null(fit)) {
    fit <- fit_call()
    unwrap_save_fit_cache(name, fit, seed, iter, warmup)
  }
  fit
}

# Convergence diagnostic gate shared by the examples. NUTS divergences.
# Not every brms/cmdstanr combination keeps the per-iteration divergence
# flag in the fit object, so probe for it and report what the fit object
# actually carries (the R-hat and ESS values below are the primary
# convergence checks either way).
diagnose_fit <- function(fit, circular_parameters = character()) {
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
