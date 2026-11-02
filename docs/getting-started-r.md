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
  c("brms", "cmdstanr", "posterior", "loo", "circular", "ggplot2"),
  repos = "https://cloud.r-project.org",
)
cmdstanr::install_cmdstan()  # downloads CmdStan once (~50 MB)
```

`brms` pulls in the rest of the Bayesian toolchain (the Stan interface,
`posterior`, ...). `cmdstanr` is the Stan backend the models use, and
`circular` provides the von Mises sampling and angle utilities.

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

Figures render inline and are saved as vector files under
`figures/circular_regression/` (gitignored). The tutorial is also the
only R notebook executed by the CI pipeline
(`.github/workflows/notebooks.yml`), with shortened sampling.
