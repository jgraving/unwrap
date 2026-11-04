# R setup

The R notebooks are verified against a pinned package environment that
ships with the repository in the form of an
[renv](https://rstudio.github.io/renv/) lockfile (`renv.lock`). The
one-time setup below restores exactly those versions — checksum verified,
from CRAN and the Stan team's r-universe — so you do not have to install
and hope. Developed on R 4.5 (anything from 4.2 up should work); the
lockfile records R 4.5.3.

## One-time setup

The compiled R packages (`brms`, `sf`, `ragg`, `stringi`, ...) install
from CRAN source, so the machine needs a working C/C++/Fortran toolchain
plus the development libraries those packages link against (GDAL/PROJ/GEOS
for `sf`, freetype/fontconfig for `ragg`, ICU for `stringi`). On
Ubuntu/Debian:

```bash
sudo apt-get install -y build-essential gfortran \
  libgdal-dev libproj-dev libgeos-dev libudunits2-dev \
  libfreetype-dev libfontconfig1-dev libpng-dev libicu-dev
```

On macOS, the Homebrew equivalents (`gdal proj geos udunits freetype
icu4c`) fill the same roles; on Windows, [Rtools 40](https://cloud.r-project.org/utils/Rtools-40/)
plus the same native libraries (the [brms installation
instructions](https://paulbuerkner.com/brms/) walk through the Windows
toolchain in detail).

Then, from the repository root, in a terminal:

```bash
Rscript -e 'install.packages("renv", repos = "https://cloud.r-project.org")'
Rscript -e 'renv::restore()'
Rscript -e 'renv::activate()'
Rscript -e 'cmdstanr::install_cmdstan()'
Rscript -e 'IRkernel::installspec()'
```

- `renv::restore()` installs the 118 pinned packages — the rest of the
  stack from CRAN, and `cmdstanr`, `loo`, `bayesplot`, and `posterior`
  from the Stan team's `stan-dev.r-universe.dev` — into a project-local
  library under `renv/`; every tarball is SHA-256 checked against the
  lockfile. The first run downloads about 110 MB of source, then
  compiles (allow roughly 20–60 minutes on a typical laptop, less on a
  fast machine); later re-runs are fast.
- `renv::activate()` writes a machine-local `.Rprofile` (gitignored) that
  puts any R session started in the repository — including the Jupyter
  R kernel — on that project library. RStudio also picks the project up
  automatically when you open `unwrap.Rproj`.
- `cmdstanr::install_cmdstan()` downloads the CmdStan compiler once
  (~50 MB) into your home directory.
- `IRkernel::installspec()` registers the **R** kernel so
  `jupyter lab` / VS Code / Positron can run the notebooks.

`brms` is the modeling front end (its formula interface and the
repository's custom von Mises family, `R/unwrap_family.R`); `circular`
provides the von Mises sampling and angle utilities; `ragg` renders
inline PNG figures on headless machines; `arrow` reads the published
Parquet data used by the locust-density and zebra examples. The
data-loading cells require Git and fetch each study's data on first use
into a repository-relative `unwrap-data/` checkout. Set
`UNWRAP_DATA_ROOT` to an existing data checkout to reuse local files.

Start the Jupyter frontend **from the repository root** so the R kernel
inherits the activated project:

```bash
cd unwrap
jupyter lab
```

The notebook setup cells remain self-contained as a fallback: run outside
the activated project (no `.Rprofile`), they install any missing packages
from CRAN. Inside an activated project the pinned set is already there,
so those cells are no-ops.

### Keeping the lockfile current

Maintainers: after deliberately upgrading a package, re-verify the
notebooks and regenerate the lockfile from the repository root with
`Rscript -e 'renv::snapshot()'`, committing the updated `renv.lock`. Do
not commit the machine-local state under `renv/` (see `renv/.gitignore`)
or the generated `.Rprofile`.

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

Set `UNWRAP_R_QUICK=1` for the short smoke budget (reduced draws); unset
it for full sampling. The grouped models use longer adaptation and
retained-draw budgets than the introductory GLM in both modes, and
always keep the full datasets and model structures. Review R-hat, bulk
ESS, the divergence report, and predictive checks before interpreting
or publishing results.
