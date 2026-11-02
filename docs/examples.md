# Biological examples

Start with the
[circular regression tutorial](notebooks/circular_regression.ipynb)
for a complete introduction using simulated observations. It fits a circular
generalized linear model (GLM) with predictors for both mean direction and
concentration, using Bambi's formula interface and PyMC for inference.

The research examples extend the same framework to biological experiments.
Repeated observations can be represented with generalized linear mixed models
(GLMMs), which estimate animal-specific variation alongside population-level
effects. The [model guide](models.md) connects these terms to biological
questions and model formulas.

The [notebook guide](https://github.com/jgraving/unwrap/blob/main/notebooks/README.md)
lists all examples, their inputs, and what each validated run covers.

## Zebra: reorientation and directional consistency

The [zebra analysis](notebooks/zebra_escape.ipynb)
compares body orientation before and after a human approach. It distinguishes
changes in the group's mean direction, alignment between animals, and
consistency within animals.

## Homing bats: alignment with the destination

The [homing bat analysis](notebooks/goldshtein_bats.ipynb)
examines whether kilometer-scale bat flights align with the roost or the lake
destination, and how strongly. It fits a von Mises–uniform mixture model to
goal-relative bearings and compares it with a unimodal baseline.

## Locusts: preferred directions and responses to virtual neighbours

The [mixture analysis](notebooks/locust_direction_mixture_model.ipynb)
compares a single preferred direction with two preferred directions in a
virtual reality experiment.

The [density and order analysis](notebooks/locust_density_order.ipynb)
examines responses to the number and directional alignment of virtual
neighbours, with a model-selection sweep over the density and order
specifications.

## Dung beetles: orientation cues

The [beetle analysis](notebooks/beetles.ipynb)
examines heading responses under manipulated orientation cues.

## Before running one of these notebooks

Each notebook on this page was rerun end to end from this checkout on the
current stack, and the outputs shown come from that validated run. The
[Python setup guide](getting-started.md) covers the environment. The
animal-behaviour notebooks load their observations through the
`unwrap.datasets` accessors, which fetch the raw data from the
[`jgraving/unwrap-data`](https://github.com/jgraving/unwrap-data) repository
on first use (cloned to `~/.cache/unwrap-data`) and then resolve to the
cache without network access.
