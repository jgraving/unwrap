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

A circular model respects the fact that 359° and 1° are close together.
It can estimate **mean direction** and **directional concentration** together,
while accounting for repeated observations from the same animal.

## Start with a biological question

| Your question | What the model describes |
| --- | --- |
| Does a treatment change where animals point? | A change in mean direction |
| Are headings more consistent under one condition? | A change in directional concentration |
| Do individuals differ in their responses? | Individual differences in a hierarchical model |
| Are there two preferred directions? | A mixture of directional distributions |

[Install `unwrap`](getting-started.md), then follow the
[circular regression tutorial](https://github.com/jgraving/unwrap/blob/main/notebooks/circular_regression.ipynb).
It generates its own observations and shows how to fit a model, interpret its
predictions, and check how well it describes the data.

## Read a fitted model

**Mean direction** tells you where observations tend to cluster. **Concentration**
tells you how tightly they cluster: near zero, directions are nearly uniform;
larger values mean greater directional consistency. Concentration is not the
same as certainty about an estimated mean. Both estimates have uncertainty.

Bayesian models describe that uncertainty with distributions of plausible
parameter values. The examples use these distributions to show predicted
directions and differences between conditions, together with uncertainty
intervals. Before interpreting a result, check both the sampling diagnostics
and whether simulated observations resemble the measured behaviour.

A model can also account for repeated measurements within animals. This helps
separate consistency within an animal from alignment between animals.

## Go further

- [Biological examples](examples.md): choose an analysis and check its inputs.
- [Circular links](links.md): how predictors map to directions and concentration.
- [Python API guide](api.md): components for Bambi and PyMC models.
- [R and MATLAB](languages.md): setup and examples for the other implementations.
- [Accompanying manuscript](https://github.com/jgraving/circstats/blob/main/elife-template.tex):
  the statistical framework and biological applications.

The introductory tutorial is checked automatically. Research notebooks have
additional data and validation requirements, described in the example guide.
