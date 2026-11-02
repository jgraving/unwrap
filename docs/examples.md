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

## Grevy's zebra: escape response to an approaching human observer

A focal group of Grevy's zebra (*Equus grevyi*) was recorded with an overhead
drone and automated tracking at the Mpala Conservancy in Laikipia, Kenya: the
first flights capture baseline (grazing) body orientation, and the final
flight contains a disturbance event in which a human observer approaches the
group. The [zebra analysis](notebooks/zebra_escape.ipynb) compares
body orientation in the **pre-stimulus baseline period** (pre-approach) with
the **post-stimulus period** (post-approach), using circularly averaged
orientation over 10 s windows. The hierarchical model separates three
complementary aspects of the escape response: a group-level shift in the
mean direction, alignment between individuals (how the individual mean
directions cluster), and within-individual consistency (how tightly each
zebra's observations cluster around its own mean direction).

## Homing bats: goal-directedness and directional precision under sensory deprivation

The [homing bat analysis](notebooks/goldshtein_bats.ipynb) revisits
kilometer-scale homing flights of echolocating bats from the experiments of
Goldshtein et al.: free-flying bats translocated within their home range are
tracked on their return to a known goal, the roost or the Agamon–Hula lake
region, under a control condition or one of several sensory-deprivation
treatments (temporary occlusion of vision, experimental disruption of
magnetoreception and olfaction). The original study found that bats typically
do home successfully even under sensory deprivation; this reanalysis asks
whether that deprivation nonetheless produces detectable changes in
goal-relative orientation. Because the treatments could change how straight a
bat flies, the reanalysis avoids feature-based selection on post-treatment
quantities (for example a "homing" threshold) and instead fits a von
Mises–uniform mixture to the full set of goal-relative bearings, in which
each step arises either from a directed von Mises component (goal-directed
"exploitation") or an undirected uniform component (search "exploration"),
compared with a unimodal von Mises baseline.

## Marching locusts: order and density, optic flow, and vectorial encoding

The [order and density analysis](notebooks/locust_density_order.ipynb)
revisits a virtual-reality experiment in which marching locusts walk on a
motion-compensating ball surrounded by a panoramic projection of moving
virtual conspecifics, parameterized by *order* (how aligned the virtual
locusts' motion is, 0 to 1) and *density* (virtual locusts per m$^2$,
1–64). The question is whether the focal locust's directional alignment is
driven primarily by the directional structure of the neighbor motion (order)
or by the number of neighbors (density); a PSIS-LOO model-selection sweep
over the order and density specifications supports an order-dominant
interpretation - order captures the bulk of the improvement in out-of-sample
performance, and adding density provides at most a secondary gain.

The [mixture analysis](notebooks/locust_direction_mixture_model.ipynb) tests
a second virtual-reality paradigm, in which moving virtual conspecifics are
placed symmetrically at a lateral distance of 8 or 50 cm. A hierarchical
bimodal von Mises mixture fitted to the full distribution of headings
discriminates the two mechanistic accounts: optic flow (headings
concentrated near the forward march direction, consistent with optomotor
alignment to wide-field motion) versus vectorial encoding (a choice-like
commitment to one of the two lateral swarms, with leftward and rightward
turns).

## Dung beetles: artificial light reorients orientation

The [beetle analysis](notebooks/beetles.ipynb) revisits a two-condition
arena experiment in the nocturnal ball-rolling dung beetle *Escarabaeus
satyrus* from the experiments of Foster et al.: each beetle is placed at the
center of a circular arena and allowed to roll to the boundary, with its exit
bearing recorded repeatedly under natural nighttime illumination (control)
and with an artificial floodlight switched on (treatment). The model fits the
treatment effect on both the mean direction and the concentration, and the
analysis quantifies, in a single coherent fit, the group-level reorientation
toward the floodlight, the change in within-individual concentration, the
post-hoc change in the between-individual clustering of the individuals' mean
directions ($\kappa_\mu$), and "floodlight attraction" as the posterior
predictive mass within $\pm45^{\circ}$ of the floodlight direction.

## Before running one of these notebooks

Each notebook on this page was rerun end to end from this checkout on the
current stack, and the outputs shown come from that validated run. The
[Python setup guide](getting-started.md) covers the environment. The
animal-behavior notebooks load their observations through the
`unwrap.datasets` accessors, which fetch the raw data from the
[`jgraving/unwrap-data`](https://github.com/jgraving/unwrap-data) repository
on first use (cloned to `~/.cache/unwrap-data`) and then resolve to the
cache without network access.
