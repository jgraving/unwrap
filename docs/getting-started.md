# Python setup

Use Python 3.10 or newer and install from GitHub in a virtual environment.
This keeps the notebook's dependencies together.

## Install and open the tutorial

From a terminal:

```bash
git clone https://github.com/jgraving/unwrap.git
cd unwrap
python -m venv .venv
```

Activate the environment on macOS or Linux:

```bash
source .venv/bin/activate
```

On Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
```

Then install the package and open the tutorial:

```bash
python -m pip install -e ".[notebook]"
python -m jupyterlab notebooks/circular_regression.ipynb
```

Select the kernel for this environment and run the cells in order. The tutorial
creates its own data, so no download is needed. It fits a circular generalized
linear model (GLM) using Bambi's formula interface, with `unwrap` links and
priors and PyMC performing Bayesian inference. Both mean direction and
directional concentration depend on a predictor. The tutorial shows how to
display uncertainty and compare model predictions with observations.

The [model guide](models.md) explains GLM formulas and their extension to
generalized linear mixed models (GLMMs) for repeated observations of animals
or groups.

For an existing analysis environment, install just the package with:

```bash
python -m pip install git+https://github.com/jgraving/unwrap.git
```

## Optional: use uv

If you already use [uv](https://docs.astral.sh/uv/), you can install and run the
same tutorial from the repository root with:

```bash
uv venv --python 3.12
uv pip install -e ".[notebook]"
uv run --no-project jupyter lab notebooks/circular_regression.ipynb
```

## Prepare your angles

Python model inputs use **radians**. Convert degrees with `np.deg2rad` and keep
the reference direction and clockwise or anticlockwise convention consistent.
For the von Mises examples here, `circular_modulo` wraps angles to `[-pi, pi)`.

Summarise angles on the circle: the average of 179° and −179° is near 180°, not
0°. The tutorial demonstrates circular summaries and plots. The
[link guide](links.md) explains the mathematical conventions.

## Run a biological example

The [example guide](examples.md) lists the data and extra dependencies for each
research notebook. Zebra and locust mixture observations are bundled with the
repository. The locust density analysis also needs original trajectory files
that are not included.

The introductory notebook is executed automatically in CI. Research notebooks
are not yet checked end to end. To make your own analysis reproducible, record
the repository commit and installed package versions with your results.
