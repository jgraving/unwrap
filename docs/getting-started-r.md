# R setup

Use a recent R release (developed on R 4.5; anything from 4.2 up should
work) and follow the [brms installation
instructions](https://paulbuerkner.com/brms/) for the one-time setup of
the Stan/C++ toolchain. After that the notebooks are self-contained:
the first code cell of each notebook checks for the packages it needs
and installs any that are missing from CRAN.

## One-time setup

In an R session:

```r
install.packages(
  c("brms", "posterior", "loo", "circular", "ggplot2", "ragg", "arrow"),
  repos = "https://cloud.r-project.org"
)
# As of 2026, cmdstanr is distributed by the Stan team on r-universe
# (it is no longer on CRAN); the second repo lets its CRAN dependencies
# resolve.
install.packages(
  "cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", "https://cloud.r-project.org")
)
cmdstanr::install_cmdstan()  # downloads CmdStan once (~50 MB)
```

`brms` pulls in the rest of the Bayesian toolchain (the Stan interface,
`posterior`, ...). `cmdstanr` is the Stan backend the models use, and
`circular` provides the von Mises sampling and angle utilities.
`ragg` renders inline PNG figures on headless machines; `arrow` reads
the published Parquet data used by the locust-density and zebra examples.
The data-loading cells require Git and fetch each study's data on first
use into a repository-relative `unwrap-data/` checkout. Set
`UNWRAP_DATA_ROOT` to an existing data checkout to reuse local files.

To run the notebooks in a Jupyter frontend (JupyterLab, VS Code, or
Positron), register the R kernel once:

```r
install.packages("IRkernel", repos = "https://cloud.r-project.org")
IRkernel::installspec()
```

## Install and open the tutorial

The repository hosts the R examples under `notebooks_r/` (the Python
versions live under `notebooks/`).

```bash
git clone https://github.com/jgraving/unwrap.git
cd unwrap
jupyter lab notebooks_r/circular_regression.ipynb
```

Select the **R** kernel and run the cells in order. The tutorial creates
its own data, so no download is needed. It fits a circular generalized
linear model (GLM) using `brms`'s formula interface and the repository's
custom von Mises family (`R/unwrap_family.R`), with the same circular
mean link (the modulo link) and softplus concentration link as the
Python tutorial uses with Bambi. Both mean direction and directional
concentration depend on a predictor. The tutorial shows how to display
uncertainty and compare model predictions with observations.

The [model guide](models.md) explains GLM formulas and their extension
to generalized linear mixed models (GLMMs) for repeated observations of
animals or groups. The [R and MATLAB](languages.md) page covers the
standalone R scripts (the tutorial notebook runs the same machinery as
`R/brms_example.R`, in notebook form) and the MATLAB examples.

## Research examples

After the tutorial, open `notebooks_r/beetles.ipynb`,
`notebooks_r/locust_density_order.ipynb`, or
`notebooks_r/zebra_escape.ipynb` with the R kernel. Their study titles
and biological questions match the [biological examples](examples.md).
These use the same observations, formulas, and prior scales as their
Python counterparts, with circular individual effects represented using
the repository's nonlinear brms machinery where needed.

The beetle notebook loads per-trial exit bearings. The zebra notebook
loads the published 10-second block averages. The locust-density
notebook downloads approximately 870 MB of tracking parts, reproduces
the trajectory rediscretization, fits five individual-level screening
models, and fits the final step-level hierarchy; allow substantial time
and memory for that complete workflow. No research notebook is executed
while building the documentation site.

Figures render inline as PNGs and are saved as vector PDFs under
`figures/<notebook_stem>/` (gitignored). No GitHub workflow executes
these notebooks; they are verified end to end by the local notebook
test (see the [development](development.md) page):

```bash
UNWRAP_TEST_NOTEBOOKS=1 uv run --no-project python -m pytest tests/test_notebooks.py -v
```

Set `UNWRAP_R_QUICK=1` for the short two-chain smoke budget; unset it
for full sampling. The grouped models use longer adaptation and
retained-draw budgets than the introductory GLM in both modes, and
always keep the full datasets and model structures. Review R-hat, bulk
ESS, the divergence report, and predictive checks before interpreting
or publishing results.

