# Unwrap

Unwrap supplies circular links, priors, and distribution helpers for Bayesian
models of individual angular observations. Build formula-based models with
Bambi or use the components directly with PyMC.

Start with [Python setup](getting-started.md), then run the
[circular regression tutorial](https://github.com/jgraving/unwrap/blob/main/notebooks/circular_regression.ipynb).
It generates its own data, fits mean direction and concentration, and checks
posterior predictive behaviour.

## Choose your entry point

| Task | Guide |
| --- | --- |
| Install from GitHub with uv or pip | [Python setup](getting-started.md) |
| Understand wrapping, concentration links, and aliases | [Circular links](links.md) |
| Find package components and their limits | [Python API guide](api.md) |
| Explore the other language implementations | [R and MATLAB](languages.md) |
| Run tests and build these docs | [Development](development.md) |

The package is a work in progress. The Python test suite and tutorial provide a
starting point for reproducible use; the R and MATLAB scripts have separate
setup and validation requirements. Cross-language equivalence is not established.
