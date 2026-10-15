# 1. The von Mises distribution

An animal's headings may cluster around one direction, with some variation
from observation to observation. The von Mises distribution describes this
pattern on the circle. It is the observation model used throughout the
[GLM](glms.md), [GLMM](glmms.md), and [mixture](mixtures.md) extensions below.

The [conceptual guide](../models.md) introduces the same ideas through
biological questions and Bambi formulas.

| Level | Biological question | Model |
| --- | --- | --- |
| 1. [Von Mises](von-mises.md) | Where do headings cluster, and how tightly? | One circular distribution |
| 2. [GLMs](glms.md) | How do treatments or continuous predictors change heading? | Regression for mean direction and concentration |
| 3. [GLMMs](glmms.md) | How do animals differ around the population response? | Partially pooled animal or group effects |
| 4. [Mixtures](mixtures.md) | Are headings consistent with more than one process or preferred direction? | Weighted combinations of circular distributions |

## A likelihood for angular observations

For a heading $y$, mean direction $\mu$, and concentration $\kappa$, the density is

$$
\mathcal{M}(y\mid\mu,\kappa)
=\frac{\exp\!\left[\kappa\cos(y-\mu)\right]}{2\pi I_0(\kappa)},
\qquad \kappa\geq 0.
$$

| Symbol | Meaning |
| --- | --- |
| $y$ | An observed angle in radians |
| $\mu$ | Mean direction, around which observations cluster |
| $\kappa$ | Concentration: larger values give tighter clustering |
| $I_0(\kappa)$ | The modified Bessel function of the first kind, order zero |
| $2\pi I_0(\kappa)$ | The normalising factor, also written $Z(\kappa)$; it makes the density integrate to one |

The cosine measures agreement between an observation and the mean direction.
It treats directions a full turn apart as identical. Concentration controls
how strongly observations cluster around that mean. At $\kappa=0$, the density
is uniform, $1/(2\pi)$, and no mean direction is preferred. Positive values
describe progressively stronger directional consistency.

This is a **probability density** over angles. A probability for an angular
range is found by integrating the density over that range; the height at an
exact direction is not a probability.

## Angular conventions

All model equations use radians. The Python circular-modulo function represents
angles in $[-\pi,\pi)$. The alternative $(-\pi,\pi]$ convention is equivalent:
the endpoints describe the same direction. Keep the angular origin
and clockwise or anticlockwise convention consistent across observations and
predictions.

## Concentration and circular variation

For a von Mises distribution, the mean resultant length is

$$
R(\kappa)=\frac{I_1(\kappa)}{I_0(\kappa)},
\qquad V=1-R(\kappa).
$$

$I_1$ is the corresponding Bessel function of order one. $R$ describes
directional agreement on a scale from zero toward one; $V$ is circular
variance. These describe variation in the distribution of observations,
rather than uncertainty about an estimated mean. The
[Bessel helpers](../api.md#concentration-conversion-unwrapbessel) convert
between von Mises concentration and mean resultant length.

**Next:** [Relate mean direction and concentration to predictors with a GLM](glms.md).
