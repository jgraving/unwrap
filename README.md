<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/assets/unwrap-logo-dark.svg">
    <img src="docs/assets/unwrap-logo.svg" alt="unwrap — Bayesian circular regression" width="560" height="193">
  </picture>
</p>

`unwrap` helps you explain where animals tend to point, how consistent their
directions are, and how these patterns change with experimental conditions.
It brings circular data into **Bayesian generalized linear models (GLMs)**:
regression models that relate an outcome, such as heading, to predictors such
as treatment, neighbour density, or a directional cue.

For repeated observations of animals or groups, **generalized linear mixed
models (GLMMs)** add individual or group differences to those relationships.
These multilevel models can separate variation within animals from variation
between animals while estimating population-level responses.

The Python package targets **PyMC** for flexible model construction and Bayesian
inference, and **Bambi**, a high-level, formula-based interface built on PyMC
for GLMs and GLMMs. `unwrap` supplies circular links and prior helpers for these existing
APIs. The repository also includes examples in **R** and **MATLAB**.

**Start here:** [Read the guide](https://jgraving.github.io/unwrap/) ·
[Try a worked example](notebooks/circular_regression.ipynb) ·
[Explore the biological analyses](notebooks/README.md)

## What can you learn from a circular model?

A change in average direction and a change in directional consistency are
different biological responses. The models here estimate both:

- **Mean direction:** the direction around which observations tend to cluster.
- **Directional concentration:** how tightly observations cluster around that
  direction. Higher concentration means more consistent directions.
- **Individual differences:** repeated observations can be grouped by animal,
  separating differences between animals from variation within an animal.

The framework uses the **von Mises distribution**, a circular counterpart of
the normal distribution that treats 359° and 1° as neighbouring directions.
Predictors can influence mean direction, concentration, or both. For example,
animals might turn toward a cue without becoming more consistent, or become
more consistent without changing their mean direction. Modelling both
parameters is called **distributional regression**; it extends the usual GLM
approach beyond modelling only the mean.

The [model guide](docs/models.md) connects these biological questions to GLM
and GLMM formulas, links, and interpretation. The
[mathematical guide](docs/math/von-mises.md)
builds from the von Mises likelihood through GLMs, GLMMs, and mixture
extensions, with equations and biological interpretations at each level.
For the full treatment, see the
[accompanying manuscript](https://github.com/jgraving/circstats/blob/main/elife-template.tex).

## Python installation

Install from GitHub in a virtual environment:

```bash
python -m venv .venv
# macOS / Linux
source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1

python -m pip install --upgrade pip
python -m pip install git+https://github.com/jgraving/unwrap.git
```

Python 3.10 or newer is required. Install using the GitHub URL shown above.
See the [setup guide](docs/getting-started.md) to run a notebook in this environment.

### Install with uv

With [uv installed](https://docs.astral.sh/uv/getting-started/installation/):

```bash
uv venv --python 3.12
uv pip install git+https://github.com/jgraving/unwrap.git
uv run --no-project python -c "import unwrap"
```

To add `unwrap` to an existing uv project instead:

```bash
uv add git+https://github.com/jgraving/unwrap.git
```

For local development, clone the repository and install in editable mode:

```bash
git clone https://github.com/jgraving/unwrap.git
cd unwrap
python -m pip install -e .
```

## Quick start: circular regression in Python

This example generates synthetic angles and fits a circular GLM in which both
mean direction (`mu`) and concentration (`kappa`) depend on a predictor `x`.
In a biological study, `x` might be a standardised measure of neighbour density
and `y` the animal's heading. Bambi's formulas `y ~ x` and `kappa ~ x` specify
the two relationships; `unwrap` supplies their links and a circular prior.

```python
import numpy as np
import bambi as bmb
import arviz as az

from unwrap.data import generate_circular_regression_data
from unwrap.link.bambi import LINKS
from unwrap.prior import CircularUniformPrior

data = generate_circular_regression_data(
    num_samples=200,
    x_range=(-1, 1),
    mean_intercept=0,
    mean_slope=1,
    conc_intercept=2,
    conc_slope=0.5,
    mean_link="circular_modulo",
    conc_link="softplus",
    seed=42,
)

model = bmb.Model(
    bmb.Formula("y ~ x", "kappa ~ x"),
    data,
    family="vonmises",
    link={"mu": LINKS["circular_modulo"], "kappa": LINKS["softplus"]},
    priors={
        "Intercept": CircularUniformPrior(),
        "x": bmb.Prior("Normal", mu=0, sigma=0.5 * np.pi),
        "kappa": {
            "Intercept": bmb.Prior("Normal", mu=2, sigma=0.5),
            "x": bmb.Prior("Normal", mu=0, sigma=1),
        },
    },
    center_predictors=False,
)

idata = model.fit(draws=1000, tune=1000, chains=4, random_seed=42)

print(
    az.summary(
        idata,
        var_names=["Intercept", "x", "kappa_Intercept", "kappa_x"],
        circ_var_names=["Intercept"],
    )
)
```

Angles are in **radians**. The `circular_modulo` link wraps the mean predictor
to `[-pi, pi)`; the `softplus` link maps the concentration predictor to positive
values. Larger `kappa` means a more concentrated distribution of angles.
These **inverse links** turn the intercept-plus-slope predictions into valid
directions and concentrations. The mean slope describes angular change in
radians per unit of `x` before wrapping; the concentration slope is on the
pre-softplus scale, so interpret its size through predicted `kappa` values.
`CircularUniformPrior()` approximates a uniform distribution on the circle
using a von Mises distribution with near-zero concentration.

Priors express plausible parameter values before fitting; those above are
illustrative and should be adapted to the scale of your predictors. Inspect
sampling diagnostics and use posterior predictive checks, comparing observations
simulated from the fitted model with your data, before interpreting results.
Use circular summaries for angular parameters, as shown for the intercept.

See [the circular regression notebook](notebooks/circular_regression.ipynb)
for a longer example with plots.

## Examples and documentation

Start with [circular regression](notebooks/circular_regression.ipynb), which
generates a small dataset and walks through fitting, interpreting, and checking
a model. The [notebook guide](notebooks/README.md) describes the biological
analyses and their data requirements.

The [documentation](https://jgraving.github.io/unwrap/) includes installation,
[circular GLMs and GLMMs](docs/models.md), biological interpretation, circular
links, and a Python API reference. Readers familiar with GLMs can go directly
to the [link conventions](docs/links.md) and
[API guide](docs/api.md).

## R and MATLAB

The language implementations currently have separate interfaces and setup
requirements. Installing the Python package does not install the R or MATLAB
dependencies, and equivalent behaviour across implementations has not yet been
established by a shared automated test suite.

### R

Start with [`R/brms_example.R`](R/brms_example.R), which uses
[`R/unwrap_family.R`](R/unwrap_family.R) to define a custom von Mises family
for brms. The example loads `brms` and `circular` and uses
[`R/vonmises_data.csv`](R/vonmises_data.csv). Its relative paths assume that
the working directory is `R/`.

Additional functions and analysis are in
[`R/unwrap_functions.R`](R/unwrap_functions.R) and
[the R Markdown analysis](R/Circular%20Dataset%20Structures.Rmd).

### MATLAB

Start with [the MATLAB setup notes](Matlab/ReadmeMATLAB.txt) and
[`Matlab/UnwrapStan.m`](Matlab/UnwrapStan.m). The Stan example requires
MatlabStan, MatlabProcessManager, and CmdStan, with paths and version settings
adjusted for your installation.

[`Matlab/MatlabOnly/`](Matlab/MatlabOnly/) contains an alternative example
using MATLAB's `hmcSampler`.

## Contributing

See the [development guide](docs/development.md) for tests, documentation builds,
and the scope of automated checks.

## License

[Apache License 2.0](LICENSE.md).
