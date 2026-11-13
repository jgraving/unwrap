# R and MATLAB

Every study in the
[biological examples](examples.md) exists as a Python/R pair - the R notebook
lives next to it under `notebooks_r/` - and the repository ships MATLAB
examples of the same models. Python installation does not install the other
language stacks, and each setup below is self-contained.

## R: brms and CmdStan

The R examples fit circular observations with
[brms](https://paulbuerkner.com/brms/) using the repository's custom von
Mises family
([`R/unwrap_family.R`](https://github.com/jgraving/unwrap/blob/main/R/unwrap_family.R)),
with **cmdstanr** and CmdStan as the Stan backend. The circular structure is
carried by the same links as in the Python models: the modulo link keeps
the mean direction on the circle, and softplus keeps the concentration
positive. Circular coefficients take von Mises
priors, which can be nearly uniform (small concentration) where no prior
direction is intended; individual angular offsets use the same construction.
"Uninformed" is therefore never a flat uniform over an interval - it is a
von Mises prior with a near-zero concentration.

[`notebooks_r/circular_regression.ipynb`](notebooks_r/circular_regression.ipynb)
is the R version of the
[Python tutorial](notebooks/circular_regression.ipynb), and the remaining
notebooks mirror the
[biological examples](examples.md) and the method examples one for one.

To use one, install the one-time environment from the
[R setup page](getting-started-r.md) (the pinned package environment,
CmdStan, and, for the Jupyter route, the R kernel), then open the
example's `.qmd` in RStudio or its `.ipynb` in Jupyter (select the
**R (`ir`) kernel**) and run the cells in order.
The first code cell checks the environment and stops, with a pointer to
the setup page, if anything is missing. The study fits are long, so the
committed outputs come from a complete end-to-end run of each notebook
(see the
[notebook test](development.md#notebook-execution-local-only) on the
development page).

### Standalone scripts

- [`R/brms_example.R`](https://github.com/jgraving/unwrap/blob/main/R/brms_example.R)
  is a small self-contained brms fit: it loads `circular` and `brms`,
  sources `unwrap_family.R`, reads the included `vonmises_data.csv`, and fits
  `y ~ x` and `kappa ~ x` with the custom von Mises family. Its relative
  paths require `R/` as the working directory:

  ```r
  install.packages(c("brms", "circular"))  # one-time package setup
  setwd("R")
  source("brms_example.R")
  ```

- [`R/unwrap_functions.R`](https://github.com/jgraving/unwrap/blob/main/R/unwrap_functions.R)
  adds circular summaries and analysis helpers on top of the family (it also
  loads `cmdstanr`); it is not required by the small example script.
- [`docs/Circular Dataset Structures.Rmd`](https://github.com/jgraving/unwrap/blob/main/docs/Circular%20Dataset%20Structures.Rmd)
  is the longer early analysis of circular dataset structures. Its
  dependencies and paths still need a reproducibility pass, so treat it as a
  reference for the intended approach rather than a ready-to-run script.

## MATLAB

MATLAB support covers the same models by two routes; neither is integrated
with the Python or R environments.

### Route 1: MATLAB + Stan (MatlabStan)

Read [`Matlab/ReadmeMATLAB.txt`](https://github.com/jgraving/unwrap/blob/main/Matlab/ReadmeMATLAB.txt)
and [`Matlab/UnwrapStan.m`](https://github.com/jgraving/unwrap/blob/main/Matlab/UnwrapStan.m).
This route needs MatlabStan, MatlabProcessManager, and CmdStan, with the
configured CmdStan path and version matched in the integration files.
`UnwrapStan.m` also carries a hard-coded `stan_home` argument; point it at
your local CmdStan installation before running. The script includes a small
example dataset and plots posterior samples. Because MatlabStan provides no
tools to unwrap circular estimates, inspect the periodic-chain behavior when
interpreting its output.

### Route 2: native MATLAB sampling

[`Matlab/MatlabOnly.m`](https://github.com/jgraving/unwrap/blob/main/Matlab/MatlabOnly.m)
fits a small circular regression with `hmcSampler` and finite-difference
gradients; observations and helpers live in the one script.

[`Matlab/CircularDatasetStructures.m`](https://github.com/jgraving/unwrap/blob/main/Matlab/CircularDatasetStructures.m)
is the full dataset analysis: effects on direction, concentration, and
individual differences, with analytic gradients and helpers in
`Matlab/unwrap_functions/`. It reads the bundled CSV files in
`Matlab/circular_datasets/` (exported from the R simulations by
`R/Export Circular Datasets.R`), and resolves input and helper paths from the
script's location, so run it from MATLAB as-is. Set `quickTest = true`
before running for short chains that check execution, or `saveFigures = true`
to save figures.

Both native examples require the Statistics and Machine Learning Toolbox and
do not require Stan. The dataset analysis reports classic R-hat before and
after unwrapping circular parameters and uses PSIS-LOO, refitting for high
Pareto-k observations up to its configured limit. Short test chains are not
sufficient to assess convergence or support scientific conclusions.

## Checking that two implementations agree

If you extend these models in another language - or port an existing
analysis - pin down what "the same model" means before comparing outputs:
the same individual observations, angular units, reference direction, design
matrices, centering, likelihood, links, and priors. Compare log densities and
predictions at fixed parameter values first, and only then compare posterior
summaries, within Monte Carlo uncertainty. Matching names or similar plots do
not establish that two implementations fit the same model.
