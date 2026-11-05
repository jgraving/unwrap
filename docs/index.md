---
title: unwrap
---

<div class="unwrap-brand">
  <img class="unwrap-logo-light" src="assets/unwrap-logo.svg" alt="unwrap — Bayesian circular regression" width="560" height="193">
  <img class="unwrap-logo-dark" src="assets/unwrap-logo-dark.svg" alt="unwrap — Bayesian circular regression" width="560" height="193">
</div>

`unwrap` helps you ask how angles and directions change with experimental
conditions. Do animals turn toward a cue? Do their headings become more
consistent? Do individuals respond differently?

These questions can be addressed with **Bayesian generalized linear models
(GLMs)**, which relate an outcome such as heading to predictors such as
treatment or neighbor density. The circular GLMs here use a von Mises
distribution, a circular counterpart of the normal distribution, so that 359° and 1° are
treated as nearby directions. Predictors can explain **mean direction** and
**directional concentration** together: where animals point and how tightly
their headings cluster.

For repeated observations, **generalized linear mixed models (GLMMs)** add
animal- or group-specific differences to the regression. These hierarchical,
or multilevel, models estimate population-level relationships alongside
variation among individuals.

The Python package targets **PyMC** for model construction and Bayesian
inference, and **Bambi** for a high-level formula interface to GLMs and GLMMs.
`unwrap` supplies circular links and prior helpers that work with their APIs.

## Start with a biological question

| Your question | What the model describes |
| --- | --- |
| Does a treatment change where animals point? | A change in mean direction |
| Are headings more consistent under one condition? | A change in directional concentration |
| Do individuals differ in their responses? | Animal-specific intercepts or treatment effects in a GLMM |
| Are there two preferred directions? | A mixture of directional distributions |
| Does behavior follow a daily rhythm rather than a compass heading? | A periodic (not circular) regression on time of day |
| What should optimal goal-directed steering do? | The von Mises navigation policy |

Install the Python side (the [Python setup page](getting-started.md)) or
take the [R setup](getting-started-r.md), then do the circular
regression tutorial in [Python](notebooks/circular_regression.ipynb) or
[R](notebooks_r/circular_regression.ipynb). It generates its own
observations and shows how to fit a model, interpret its predictions, and
check how well it describes the data.
The [model guide](models.md) explains how these questions become GLM and GLMM
formulas, with biological interpretations alongside the statistical terms.

Every analysis on this site - the tutorial and all the
[biological examples](examples.md) - runs in both **Python** (PyMC and
Bambi) and **R** (brms on CmdStan); the notebooks are Python/R pairs that
use the same observations, formulas, and prior scales. The
[R setup page](getting-started-r.md) covers the one-time R environment.

## Read a fitted model

**Mean direction** tells you where observations tend to cluster. **Concentration**
tells you how tightly they cluster: near zero, directions are nearly uniform;
larger values mean greater directional consistency. Concentration is not the
same as certainty about an estimated mean. Both estimates have uncertainty.

Bayesian models combine **priors**, which describe plausible parameter values
before fitting, with the **likelihood**, which describes how those parameters
could generate the observed angles. The result is a **posterior distribution**
describing parameter uncertainty given the data and model. The examples show
predicted directions and differences between conditions, together with uncertainty
intervals. Before interpreting a result, check both the sampling diagnostics
and whether observations simulated from the fitted model resemble the measured
behavior. This comparison is a **posterior predictive check**.

A GLMM can also account for repeated measurements within animals. This helps
separate consistency within an animal from alignment between animals: animals
can each hold a steady heading while pointing in different directions.

## Go further

- [Examples](examples.md): the tutorial, the biological studies, and the method examples.
- [Circular GLMs and GLMMs](models.md): connect biological questions to model formulas.
- [Mathematical guide](math/von-mises.md): build from the von Mises likelihood to GLMs, GLMMs, and mixtures.
- [Circular links](links.md): how predictors map to directions and concentration.
- [Python API guide](api.md): components for Bambi and PyMC models.
- [R and MATLAB](languages.md): the R versions of the tutorial and every study, the standalone R scripts, and the MATLAB examples.

The introductory tutorial is executed automatically in CI, and every
notebook is rerun end to end before a change that touches it is merged.
The example guide describes each analysis and its inputs.

## Citation

> Graving, J. M., Sayin, S., Costelloe, B. R., Goldshtein, A., and Foster, J. J.
> *Unwrapping Circular Statistics with Bayesian Multilevel Models.*
