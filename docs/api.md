# Python API guide

Import helpers from their defining modules. Unwrap provides components for Bambi
and PyMC; it does not define a separate top-level model-fitting class.

## Links

See [Circular links](links.md) for the numerical functions, backend registries,
inputs, aliases, and boundary conventions. Source docstrings provide the same
conventions next to each function.

## Prior helpers: `unwrap.prior`

- `CircularUniformPrior()` creates a Bambi prior backed by a von Mises
  distribution with near-zero concentration. It approximates a circular uniform.
- `SoftplusNormalPrior(mu, sigma, *, mu_is_on_positive_scale=True)` creates a
  Normal prior whose location is converted to the pre-softplus scale when the
  flag is true, and attaches PyMC's `LogExpM1` transform. `sigma` is on the Normal
  scale. This transform controls PyMC's internal sampling coordinates; it does
  not itself turn Normal draws into a positive-valued prior. This helper is not
  used in the introductory tutorial; review its intended support before use.

The tutorial uses ordinary Normal priors for the concentration's linear
predictor and the softplus inverse link for the concentration itself.

## Distributions: `unwrap.distributions`

`CircularUniform(name, ...)` constructs a PyMC von Mises variable with a mean of
zero and concentration at machine epsilon.

`VonMisesMixture(name, w, mu, kappa, **kwargs)` wraps a PyMC mixture of von Mises
components. Weights must sum to one, means are radians, and concentrations are
positive. `.dist(w, mu, kappa, **kwargs)` constructs an unnamed distribution.

`BimodalVonMises(name, mu, kappa, p, **kwargs)` is a two-component convenience
wrapper with weights `p` and `1-p`; the final axis of `mu` and `kappa` must have
length two. `.dist(...)` is also provided. Lists, NumPy arrays, and symbolic tensors are accepted. Known invalid
component dimensions are rejected at construction; unknown symbolic dimensions
are checked at evaluation. Leading dimensions follow PyMC's broadcasting rules.

## Synthetic data: `unwrap.data`

`generate_circular_regression_data(num_samples, x_range, mean_intercept,
mean_slope, conc_intercept, conc_slope, mean_link, conc_link, seed=42)` generates
one predictor and a von Mises angular response. Supply link names from the
registry. The returned pandas DataFrame has columns `x`, `y`, `mean`, and
`concentration`; the last two are generating parameters for comparison with a fit.

## Concentration conversion: `unwrap.bessel`

`bessel_ratio(kappa)` maps a nonnegative scalar concentration to the von Mises
mean resultant length, using exponentially scaled Bessel functions. Positive
infinity maps to one; negative inputs and NaN are rejected.

`inverse_bessel_ratio(R, tolerance=1e-10, max_iterations=1000)` inverts this
mapping for `0 <= R < 1`. Zero maps to zero. Tolerance is relative to the
concentration bracket; failure to converge raises `RuntimeError`. The mapping
becomes poorly conditioned near one, where floating-point precision limits
recoverable concentration accuracy.
