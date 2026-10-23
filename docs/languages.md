# R and MATLAB implementations

Keep the language implementations together while their interfaces and scope are
being clarified. Python installation does not install the other language stacks.
The setup notes here are based on the checked-in scripts; R and MATLAB execution
is not yet covered by this repository's CI.

## R: brms and Stan

Start with [`R/brms_example.R`](https://github.com/jgraving/unwrap/blob/main/R/brms_example.R).
It loads `circular` and `brms`, sources `unwrap_family.R`, and reads the included
`vonmises_data.csv`. Follow the [brms installation instructions](https://paulbuerkner.com/brms/)
for the Stan/C++ toolchain, then from an R session in the repository root:

```r
install.packages(c("brms", "circular"))  # one-time package setup
setwd("R")
source("brms_example.R")
```

The example fits `y ~ x` and `kappa ~ x` using a custom von Mises family. Its
relative paths require `R/` as the working directory. `unwrap_functions.R`
contains additional circular summaries and analysis helpers and also loads
`cmdstanr`; it is not required by the small `brms_example.R` script.

The longer analysis is in
[`Circular Dataset Structures.Rmd`](https://github.com/jgraving/unwrap/blob/main/R/Circular%20Dataset%20Structures.Rmd).
Those analysis dependencies and paths still need a dedicated reproducibility pass.

## MATLAB: Stan example

Read [`Matlab/ReadmeMATLAB.txt`](https://github.com/jgraving/unwrap/blob/main/Matlab/ReadmeMATLAB.txt)
and [`Matlab/UnwrapStan.m`](https://github.com/jgraving/unwrap/blob/main/Matlab/UnwrapStan.m).
The current route needs MatlabStan, MatlabProcessManager, and CmdStan.
The setup notes describe matching the configured CmdStan path and version in
the integration files. `UnwrapStan.m` itself also has a hard-coded `stan_home`
argument that must be changed to your local installation.

The script includes a small example dataset and plots posterior samples.
It notes that MatlabStan does not provide tools to unwrap estimates, so inspect
periodic chain behaviour when interpreting its output.

## MATLAB: native sampling

[`Matlab/MatlabOnly.m`](https://github.com/jgraving/unwrap/blob/main/Matlab/MatlabOnly.m)
fits a small circular regression using `hmcSampler` and finite-difference
gradients. It includes the observations and helper functions in one script.

[`Matlab/CircularDatasetStructures.m`](https://github.com/jgraving/unwrap/blob/main/Matlab/CircularDatasetStructures.m)
provides examples of effects on direction, concentration, and individual
differences, using analytic gradients and helpers in `Matlab/unwrap_functions/`.
It reads the bundled CSV files in `Matlab/circular_datasets/`, exported from the
R simulations by `R/Export Circular Datasets.R`. To run the analysis, open the
script in MATLAB and run it; input and helper paths are resolved from the
script's location. Set `quickTest = true` before running for short chains that
check execution, or `saveFigures = true` to save figures.

Both native examples require the Statistics and Machine Learning Toolbox and
do not require Stan. The dataset analysis reports classic R-hat before and
after unwrapping circular parameters and uses PSIS-LOO, with refits for high
Pareto-k observations up to its configured limit. Short test chains are not
sufficient to assess convergence or support scientific conclusions.

## Establishing equivalence later

A shared example should specify the same individual observations, angular units,
reference direction, design matrices, centering, likelihood, links, and priors.
Compare log densities and predictions at fixed parameters before comparing
posterior summaries within Monte Carlo uncertainty. Matching names or similar
plots do not establish that the implementations fit the same model.
