# R setup

There are two ways to run the R examples, and both share the same
one-time environment:

- **RStudio + Quarto (recommended)**: open the `.qmd` file and press
  *Run All*. It runs with plain R — no Jupyter at all.
- **Jupyter**: the same example as a `.ipynb` with the R (`ir`) kernel
  selected. The documentation site renders from these.

The environment itself is pinned in `renv.lock`, which plays the role
`requirements.txt` plays for Python: package versions and sources are
recorded for the environment used to verify the examples. And
`renv::restore()` plays the role `uv sync` plays: one command installs
exactly those versions, into a project-local library, so you do not
have to install anything by hand and hope.

Use R 4.5.3, the version recorded in `renv.lock` and tested by CI.

## One-time environment

### 1. System packages (the only OS-specific step)

The pinned R packages link against a handful of native libraries
(GDAL/PROJ/GEOS for `sf`, freetype/fontconfig for `ragg`, ICU for
`stringi`, libuv for `fs`, zmq for the Jupyter R kernel). The tested
Ubuntu and Windows restores use
prebuilt binaries for the pinned CRAN packages; Ubuntu also builds the
four r-universe packages from source. On macOS the restore uses a mix of
binaries and source builds, so the listed headers and libraries support
packages that compile locally.

**Linux (Ubuntu 24.04, the tested runner):**

```bash
sudo apt-get update
sudo apt-get install -y build-essential gfortran \
  libgdal-dev libproj-dev libgeos-dev libudunits2-dev \
  libfreetype-dev libfontconfig1-dev libpng-dev libicu-dev \
  libuv1-dev libzmq3-dev
```

**macOS:** install R 4.5.3 from CRAN ([Apple silicon](https://cran.r-project.org/bin/macosx/big-sur-arm64/base/R-4.5.3-arm64.pkg)
or [Intel](https://cran.r-project.org/bin/macosx/big-sur-x86_64/base/R-4.5.3-x86_64.pkg)),
then install the system libraries with Homebrew:

```bash
xcode-select --install
brew install gdal proj geos udunits freetype icu4c fontconfig libpng libuv
```

**Windows:** install [R 4.5.3](https://cran.r-project.org/bin/windows/base/old/4.5.3/)
and [Rtools 4.5](https://cran.r-project.org/bin/windows/Rtools/). The tested
package restore used prebuilt Windows binaries, but the CmdStan step below
builds C++ code and needs the toolchain.

### 2. The pinned R packages (identical on every OS)

From the **repository root**, in a terminal:

```bash
Rscript -e 'ru <- "https://stan-dev.r-universe.dev"; repos <- c(CRAN = "https://packagemanager.posit.co/cran/2026-11-03", setNames(ru, ru)); install.packages("renv", repos = repos); renv::activate(); renv::load(); renv::restore(repos = repos)'
```

Both URLs need to be in the same command (not just in `renv.lock`):
passing `repos` to `renv::restore()` replaces the repositories recorded
in the lockfile, and four pinned entries (`cmdstanr`, `loo`,
`posterior`, `bayesplot`) record the Stan r-universe repository URL as their repository ID. The
`CRAN` entry routes the other 129 entries to the RSPM snapshot above.

- `renv::restore()` restores the 133 package records — 129 from the
  `CRAN` entry and four from `stan-dev.r-universe.dev` — into a
  project-local library under `renv/`. The `CRAN` entry points to Posit's
  dated **RSPM** snapshot. CI installed its pinned CRAN packages as
  binaries on Ubuntu 24.04, Windows, and macOS; Ubuntu built the four
  r-universe packages from source, and macOS built `pbdZMQ` and `IRkernel`
  from source. The `Hash` fields in `renv.lock` identify packages for
  renv's cache; they are not archive checksums. Allow more time when
  source builds are needed; later re-runs can reuse cached packages.
- `renv::activate()` writes a machine-local `.Rprofile` (gitignored)
  that puts any R session started in the repository — including the
  Jupyter R kernel — on that project library. RStudio also activates
  the project automatically when you open `unwrap.Rproj`.

### 3. CmdStan (one-time, ~50 MB)

```bash
Rscript -e 'cmdstanr::install_cmdstan(version = "2.40.0")'
```

Downloads the Stan compiler into your home directory. (The examples'
setup cells check for it and will download it automatically on first
use if it is missing.)

### 4. Jupyter mode only

```bash
Rscript -e 'IRkernel::installspec()'
```

Registers the **R** kernel so `jupyter lab` / VS Code / Positron can
run the `.ipynb` examples. Start the Jupyter frontend **from the
repository root** so the R kernel inherits the activated project:

```bash
cd unwrap
jupyter lab
```

## Running the examples

The R examples live under `notebooks_r/`; each one exists both as a
Quarto twin (`.qmd`, opened in RStudio) and as a Jupyter notebook
(`.ipynb`, opened in Jupyter). The `.ipynb` is the source of truth; the
`.qmd` twins are generated from it (`tools/ipynb_to_qmd.py`), and each
example's documentation page carries a "Download Quarto" button.

```bash
git clone https://github.com/jgraving/unwrap.git
cd unwrap
jupyter lab   # or: open unwrap.Rproj in RStudio
```

- **RStudio/Quarto**: open `notebooks_r/<example>.qmd` and run all
  cells.
- **Jupyter**: open `notebooks_r/<example>.ipynb`, select the **R**
  kernel, and run the cells in order.

The `brms` package is the modeling front end (its formula interface and
the repository's custom von Mises family, `R/unwrap_family.R`);
`circular` provides the von Mises sampling and angle utilities; `ragg`
renders figures on headless machines; `arrow` reads the published
Parquet data used by the locust-density and zebra examples. The
data-loading cells require Git and fetch each study's data on first use
into a repository-relative `unwrap-data/` checkout. Set
`UNWRAP_DATA_ROOT` to an existing data checkout to reuse local files.

The first code cell of each example re-checks the environment and
restores missing pieces from CRAN, so on a machine set up as above it
is a no-op — and running an example outside the activated project
(without the `.Rprofile`) still works, just not at the pinned versions.

## The tutorial

Start with `circular_regression`: it creates its own data, so no
download is needed. It fits a circular generalized linear model (GLM)
using `brms`' formula interface and the repository's custom von Mises
family (`R/unwrap_family.R`), with the same circular mean link (the
modulo link) and softplus concentration link as the Python tutorial
uses with Bambi. Both mean direction and directional concentration
depend on a predictor. The tutorial shows how to display uncertainty
and compare model predictions with observations. It is the cheapest R
example to run end to end.

The [model guide](models.md) explains GLM formulas and their extension
to generalized linear mixed models (GLMMs) for repeated observations of
animals or groups. The [R and MATLAB](languages.md) page covers the
standalone R scripts (the tutorial example runs the same machinery as
`R/brms_example.R`, in notebook form) and the MATLAB examples.

## Research examples

After the tutorial, try `beetles`, `locust_density_order`, or
`zebra_escape` (RStudio/Quarto or Jupyter). Their study titles and
biological questions match the [biological examples](examples.md).
These use the same observations, formulas, and prior scales as their
Python counterparts, with circular individual effects represented using
the repository's nonlinear brms machinery where needed.

The beetle example loads per-trial exit bearings. The zebra example
loads the published 10-second block averages. The locust-density
example downloads approximately 870 MB of tracking parts, reproduces
the trajectory rediscretization, fits five individual-level screening
models, and fits the final step-level hierarchy; allow substantial time
and memory for that complete workflow. No research example is executed
while building the documentation site.

Figures render inline and are saved as vector PDFs under
`figures/<example_stem>/` (gitignored). No GitHub workflow executes
these examples; they are verified end to end by the local notebook
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
