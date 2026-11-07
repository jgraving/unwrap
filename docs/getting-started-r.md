# R setup

There are two ways to run the R examples, and both share the same
one-time environment:

- **RStudio + Quarto (recommended)**: open the `.qmd` file and press
  *Run All*. It runs with plain R — no Jupyter at all.
- **Jupyter**: the same example as an `.ipynb` with the **R** kernel.

The environment is pinned in `renv.lock`: it records the exact version
of every package this project uses, and the short R script in step 3
installs precisely those versions into a project-local library,
separate from your system packages. Re-running it is a fast no-op once
everything is in place. Use R **4.5.3**, the version recorded in
`renv.lock`.

## Install and open the tutorial

### 1. Get the repository

```bash
git clone https://github.com/jgraving/unwrap.git
cd unwrap
```

Every step below runs from this directory.

### 2. System packages (the only OS-specific step)

A few pinned packages link against native libraries, and on some
operating systems a few are built from source; install those
prerequisites first. The C++ toolchain is also used to build CmdStan
in step 4.

**Linux (tested on Ubuntu 24.04):**

```bash
sudo apt-get update
sudo apt-get install -y build-essential gfortran \
  libgdal-dev libproj-dev libgeos-dev libudunits2-dev \
  libfreetype-dev libfontconfig1-dev libpng-dev libicu-dev \
  libuv1-dev libzmq3-dev
```

**macOS:** install [R 4.5.3](https://cran.r-project.org/bin/macosx/big-sur-arm64/base/R-4.5.3-arm64.pkg)
(Apple silicon; the [Intel build](https://cran.r-project.org/bin/macosx/big-sur-x86_64/base/R-4.5.3-x86_64.pkg)
follows the same steps), then:

```bash
xcode-select --install
brew install gdal proj geos udunits freetype icu4c fontconfig \
  libpng libuv
```

**Windows:** install [R 4.5.3](https://cran.r-project.org/bin/windows/base/old/4.5.3/)
and [Rtools 4.5](https://cran.r-project.org/bin/windows/Rtools/). On
Windows the pinned CRAN packages install as prebuilt binaries; Rtools
serves the CmdStan build in step 4 and any source-build fallback.

### 3. The pinned R packages

In any R session (RStudio's console or `R` in a terminal), from the
repository root:

```r
# Four Stan development packages (cmdstanr, loo, posterior, bayesplot),
# from the Stan team's r-universe:
ru <- "https://stan-dev.r-universe.dev"
# The rest, from a dated CRAN snapshot that serves prebuilt binaries
# on Linux and Windows:
cran <- "https://packagemanager.posit.co/cran/2026-11-03"
repos <- c(CRAN = cran, setNames(ru, ru))

install.packages("renv", repos = repos)
renv::activate()      # record the project library in a local .Rprofile
renv::load()
renv::restore(repos = repos)
```

The final line installs the 136 pinned packages (132 from CRAN, four
from r-universe) into a project-local `renv/` library. The CRAN
packages arrive as prebuilt binaries on Linux and Windows; the
r-universe packages — and a few others on macOS — are built from
source, so a first install can take 20–60 minutes. Re-running the
script is fast.

`renv::activate()` writes a machine-local `.Rprofile` (gitignored), so
any R session started in the repository — RStudio, scripts, or the
Jupyter R kernel — uses that project library.

### 4. CmdStan

```r
cmdstanr::install_cmdstan(version = "2.40.0")
```

Fetches the Stan compiler into your home folder (~50 MB) and builds it
with the C++ toolchain from step 2.

### 5. The Jupyter R kernel (only if you use Jupyter)

```r
IRkernel::installspec()
```

Skip this if you only use RStudio/Quarto. Registers the **R** kernel
so `jupyter lab`, VS Code, or Positron can run the `.ipynb` examples.

### Open the tutorial

- **RStudio/Quarto:** open `notebooks_r/circular_regression.qmd` and
  run all cells.
- **Jupyter:** from the repository root run `jupyter lab`, open
  `notebooks_r/circular_regression.ipynb`, select the **R** kernel,
  and run the cells in order.

The first code cell checks the environment and stops, with a pointer
back to this page, if the one-time setup has not been run — a failed
setup fails fast instead of slowly at the first fit.

The tutorial creates its own data, so no download is needed. It fits
a circular generalized linear model (GLM) with `brms`' formula
interface and the repository's custom von Mises family
(`R/unwrap_family.R`): a modulo mean link for the direction and a
softplus link for the concentration, with both depending on a
predictor.
The [model guide](models.md) explains GLM formulas and their
extension to generalized linear mixed models (GLMMs) for repeated
observations of animals or groups.

On a machine that is already set up, a fresh clone plus step 3
re-checks the pinned versions; steps 4 and 5 run once per machine.

## Prepare your angles

The models take angles in **radians** (the example data's angle
columns already are). If you work in degrees, convert with
`degrees * pi / 180` and keep one reference direction and one
clockwise or counterclockwise convention throughout the analysis.

Remember that the angles wrap: the average of 179° and −179° is near
180°, not 0°. The tutorial shows the circular summaries and plots
that handle this. The [link guide](links.md) explains the
mathematical conventions.

## Run a study example

The [example guide](examples.md) organizes the notebooks into the
circular regression tutorial (simulated toy data, not a biological
study), five biological studies, and two method examples (periodic
but not circular, and the von Mises navigation policy). After the
tutorial, try `beetles`, `locust_density_order`, or `zebra_escape`,
in RStudio or Jupyter.

The animal-behavior notebooks load each study through a data cell that
fetches the raw files from the `jgraving/unwrap-data` repository on
first use (needs Git) into a `unwrap-data/` checkout next to the
repository, then resolves from that checkout. Set `UNWRAP_DATA_ROOT`
to an existing checkout to reuse local files.

`locust_density_order` downloads roughly 870 MB of tracking parts and
fits five individual-level screening models before the final
hierarchy; allow substantial time and memory. The notebooks follow a
shorter draw budget when `UNWRAP_R_QUICK=1` is set; unset, they run
the full budgets. Review R-hat, bulk ESS, the divergence report, and
the predictive checks before interpreting results.

Every notebook's stored outputs come from an end-to-end run of that
notebook (see the [development](development.md) page).

## Skip the refit with the saved fits

The fit cells are seeded, and the fits they produce are published in
the `jgraving/unwrap-data` repository under `fit_cache/` (one RDS
file per model and draw budget). The data checkout the notebooks
create is sparse -- one study at a time -- so add that directory
first (a full clone already has it):

    git -C <data-root> sparse-checkout add fit_cache

Then start RStudio or Jupyter with

    UNWRAP_R_FIT_CACHE=1
    UNWRAP_R_FIT_CACHE_DIR=<data-root>/fit_cache

where `<data-root>` is the data checkout (by default the
`unwrap-data/` directory inside this repository, or the path you
hand to `UNWRAP_DATA_ROOT`). A snapshot is reused only when its
model name, seed, draw budget, and warmup all match the fit cell,
and the fixed-effects size is unchanged; otherwise the cell says the
cache entry is stale and refits as usual. Leave the cache unset for
`UNWRAP_R_QUICK` runs: a quick-budget fit saves under the same name
and would replace the full-budget snapshot. The stored notebook
outputs always come from runs with the cache off -- the helpers are
no-ops while it is -- so the committed runs exercise the full
sampling path.
