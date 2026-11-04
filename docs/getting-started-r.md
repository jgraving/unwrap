# R setup

There are two ways to run the R examples, and both share the same
one-time environment:

- **RStudio + Quarto (recommended)**: open the `.qmd` file and press
  *Run All*. It runs with plain R — no Jupyter at all.
- **Jupyter**: the same example as a `.ipynb` with the R (`ir`) kernel
  selected. The documentation site renders from these.

The environment itself is pinned in `renv.lock`, which plays the role
`requirements.txt` plays for Python: every package version is fixed
(checksummed) to exactly what the examples were verified against. And
`renv::restore()` plays the role `uv sync` plays: one command installs
exactly those versions, into a project-local library, so you do not
have to install anything by hand and hope.

Developed on R 4.5 (anything from 4.2 up should work); the lockfile
records R 4.5.3.

## One-time environment

### 1. System packages (the only OS-specific step)

The pinned R packages link against a handful of native libraries
(GDAL/PROJ/GEOS for `sf`, freetype/fontconfig for `ragg`, ICU for
`stringi`, libuv for `fs`, zmq for the Jupyter R kernel). On Linux and
Windows the restore in step 2 uses **prebuilt binaries**, so these are
runtime libraries (the `-dev` headers are only needed if a package ever
builds from source). On macOS the compiled packages build from source,
so they are build dependencies.

**Linux (Ubuntu/Debian):**

```bash
sudo apt-get install -y build-essential gfortran \
  libgdal-dev libproj-dev libgeos-dev libudunits2-dev \
  libfreetype-dev libfontconfig1-dev libpng-dev libicu-dev \
  libuv1-dev libzmq3-dev
```

**macOS** (Homebrew; `r` installs R itself):

```bash
xcode-select --install
brew install r gdal proj geos udunits freetype icu4c fontconfig libpng libuv
```

**Windows:** install a current R release from
<https://cloud.r-project.org>. That is all you need up front: every
pinned package ships a prebuilt Windows binary at the repositories the
step-2 command uses, so the restore is plain downloads and no compiler
is required. If a package ever falls back to a source build, install
Rtools ([Rtools 40](https://cloud.r-project.org/utils/Rtools-40/),
or in R: `install.packages("rtools40"); rtools40::install_rtools()`).

### 2. The pinned R packages (identical on every OS)

From the **repository root**, in a terminal:

```bash
Rscript -e 'ru <- "https://stan-dev.r-universe.dev"; repos <- c(CRAN = "https://packagemanager.posit.co/cran/2026-11-03", setNames(ru, ru)); install.packages("renv", repos = repos); renv::restore(repos = repos); renv::activate()'
```

Both URLs need to be in the same command (not just in `renv.lock`):
passing `repos` to `renv::restore()` replaces the repositories recorded
in the lockfile, and four of the pinned packages (`cmdstanr`, `loo`,
`posterior`, `bayesplot`) are dev versions that only exist on the Stan
teams r-universe, while the `CRAN` entry is what routes the other 129
packages to the RSPM snapshot above.

- `renv::restore()` installs the 133 pinned packages — the rest of the
  stack from the `CRAN` repository in the command, and `cmdstanr`,
  `loo`, `bayesplot`, and `posterior` from the Stan team's
  `stan-dev.r-universe.dev` — into a project-local library under `renv/`.
  The `CRAN` entry is Posit's free **RSPM** mirror of CRAN, pinned to a
  snapshot verified to carry prebuilt binaries for exactly the locked
  versions: on Linux and Windows that makes the restore a download-only,
  compile-free install (a few minutes). Every package the lockfile
  resolves is integrity-checked — source tarballs against the SHA-256
  hashes recorded in `renv.lock`, prebuilt binaries against the hashes
  the repositories publish for them. On macOS there is no free binary
  source, so the compiled packages build from source — allow roughly
  20–60 minutes on a typical laptop for the first run; later re-runs are
  fast everywhere.
- `renv::activate()` writes a machine-local `.Rprofile` (gitignored)
  that puts any R session started in the repository — including the
  Jupyter R kernel — on that project library. RStudio also activates
  the project automatically when you open `unwrap.Rproj`.
- If the snapshot URL above ever 404s (it is a permanent snapshot, so
  it normally will not), use any RSPM snapshot dated on or after the
  latest package `Date` field in `renv.lock` (for example, the date the
  lockfile was last changed in this repository's history).

### 3. CmdStan (one-time, ~50 MB)

```bash
Rscript -e 'cmdstanr::install_cmdstan()'
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
