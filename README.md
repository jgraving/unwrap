# Unwrap

`unwrap` provides tools for Bayesian circular regression: modelling angular
responses with von Mises distributions, circular links for mean direction, and
positive links for concentration. The Python package works with **Bambi** and
**PyMC**. This repository also contains **R** and **MATLAB** implementations and
examples.

<p align="center">
  <img src="assets/graphic-abstract.png" alt="Unwrap visual abstract" width="900">
</p>

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

The declared dependencies are Bambi, PyMC, and ArviZ; see
[requirements.txt](requirements.txt) for their version constraints. Use a Python
version supported by the versions of those dependencies you install. A tested
Python/dependency compatibility matrix is planned.

For local development, clone the repository and install in editable mode:

```bash
git clone https://github.com/jgraving/unwrap.git
cd unwrap
python -m pip install -e .
```

## Quick start: circular regression in Python

This example generates synthetic angles and fits a model in which both mean
direction (`mu`) and concentration (`kappa`) depend on a predictor `x`.
It uses Bambi's model interface with links and a prior supplied by `unwrap`.

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
`CircularUniformPrior()` approximates a uniform distribution on the circle
using a von Mises distribution with near-zero concentration.

The priors above are illustrative and should be adapted to the scale of your
predictors. Inspect sampling diagnostics and posterior predictive behaviour
before interpreting a fitted model. Use circular summaries for angular
parameters, as shown for the intercept.

See [the circular regression notebook](notebooks/circular_regression.ipynb)
for a longer example with plots.

## Code and examples

| Location | Contents |
| --- | --- |
| [`unwrap/link/`](unwrap/link/) | NumPy and PyTensor/PyMC link functions, plus Bambi link objects |
| [`unwrap/prior.py`](unwrap/prior.py) | Bambi prior helpers |
| [`unwrap/distributions.py`](unwrap/distributions.py) | Circular uniform helper, von Mises mixtures, and a sufficient-statistics likelihood |
| [`unwrap/data.py`](unwrap/data.py) | Synthetic circular regression data |
| [`unwrap/bessel.py`](unwrap/bessel.py) | Helpers relating von Mises concentration to mean resultant length |
| [`notebooks/`](notebooks/) | Python examples and analyses |
| [`notebooks/datasets/`](notebooks/datasets/) | Data files included for the notebooks |
| [`R/`](R/) | R functions, a brms custom family, and examples |
| [`Matlab/`](Matlab/) | MATLAB examples using Stan and MATLAB sampling tools |

Python notebooks include:

- [Circular regression](notebooks/circular_regression.ipynb)
- [Beetles](notebooks/beetles.ipynb)
- [Locust density and order](notebooks/locust_density_order.ipynb)
- [Locust direction mixture model](notebooks/locust_direction_mixture_model.ipynb)
- [Zebra escape](notebooks/zebra_escape.ipynb)

These are working analysis notebooks; their setup and dependencies are not yet
standardized. Check the installation cells and data paths before running them.
Notebook organization and reproducible setup instructions are planned.

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

## Development

Planned work includes:

- Organizing and documenting the code and its public interfaces.
- Adding automated tests and continuous integration.
- Building user and API documentation.
- Organizing example notebooks with reproducible setup instructions.
- Comparing the language implementations on shared reference examples.

Changes should be proposed as small, focused pull requests with a clear purpose
and a description of how they were checked. Keep code moves, behaviour changes,
and notebook rewrites in separate PRs where practical so each change is easy to
review and revert.

## License

[Apache License 2.0](LICENSE.md).
