# Python setup

Python 3.10 or newer is required. CI checks Python 3.10–3.14 on Linux and Python
3.12 on macOS and Windows; dependency versions vary with the Python interpreter.
The package currently uses PyMC 5 and ArviZ 0.x APIs.

## Install from GitHub

With [uv](https://docs.astral.sh/uv/getting-started/installation/):

```bash
uv venv --python 3.12
uv pip install git+https://github.com/jgraving/unwrap.git
uv run --no-project python -c "import unwrap"
```

In an existing uv project, use `uv add git+https://github.com/jgraving/unwrap.git`.
For pip, activate a virtual environment and run:

```bash
python -m pip install git+https://github.com/jgraving/unwrap.git
```

Use the GitHub URL: the `unwrap` name on PyPI currently refers to a different
project. For a reproducible analysis, record the Git commit and resolved package
versions, or pin the dependency to a commit in your project's lockfile.

## Run the introductory notebook

```bash
git clone https://github.com/jgraving/unwrap.git
cd unwrap
uv venv --python 3.12
uv pip install -e ".[notebook]"
uv run --no-project jupyter lab notebooks/circular_regression.ipynb
```

Select the kernel belonging to this environment and run all cells in order.
The tutorial uses seeded synthetic data and prints dependency versions. Its
saved outputs are cleared; CI executes it and retains an executed copy as an
artifact. Other analysis notebooks are not yet part of that execution check.

## Research notebooks

See the [notebook inventory](https://github.com/jgraving/unwrap/blob/main/notebooks/README.md)
for each analysis's purpose, bundled or missing data, additional dependencies,
and current setup issues. Only the introductory notebook is executed in CI.

## Angular conventions

Use radians for Python model inputs. Convert degrees with `np.deg2rad` and choose
a consistent reference direction and rotation convention. For the von Mises
workflow here, wrap observations with `circular_modulo` to `[-pi, pi)`.
An angular parameter needs circular summaries; averaging -179 and 179 degrees
as ordinary numbers would give a misleading direction.
