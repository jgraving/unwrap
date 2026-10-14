# Biological examples

Start with the
[circular regression tutorial](https://github.com/jgraving/unwrap/blob/main/notebooks/circular_regression.ipynb)
for a complete introduction using simulated observations.

The [notebook guide](https://github.com/jgraving/unwrap/blob/main/notebooks/README.md)
lists all examples, their inputs, and their current execution status.

## Zebra: reorientation and directional consistency

The [zebra analysis](https://github.com/jgraving/unwrap/blob/main/notebooks/zebra_escape.ipynb)
compares body orientation before and after a human approach. It distinguishes
changes in the group's mean direction, alignment between animals, and
consistency within animals. The preprocessed observations are included.

## Locusts: preferred directions and responses to virtual neighbours

The [mixture analysis](https://github.com/jgraving/unwrap/blob/main/notebooks/locust_direction_mixture_model.ipynb)
compares a single preferred direction with two preferred directions in a
virtual reality experiment. The included data are preprocessed heading steps.

The [density and order analysis](https://github.com/jgraving/unwrap/blob/main/notebooks/locust_density_order.ipynb)
examines responses to the number and directional alignment of virtual
neighbours. Its raw trajectory preprocessing requires files not included in
this repository.

## Dung beetles: orientation cues

The [beetle analysis](https://github.com/jgraving/unwrap/blob/main/notebooks/beetles.ipynb)
examines heading responses under manipulated orientation cues. Its data file
is included in the repository.

## Before using a research notebook

These notebooks retain research outputs and are being prepared for easier
reuse. They have not yet been validated end to end with the current package
environment. Check the notebook guide for setup requirements and limitations;
use the introductory tutorial as the checked starting point.
