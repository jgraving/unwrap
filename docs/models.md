# Circular GLMs and GLMMs

Suppose you want to know whether a directional cue changes how an animal moves.
It could change the animal's preferred heading, make its headings more
consistent, or do both. A circular regression models the observed headings
directly and estimates these responses separately, with uncertainty.

For the equations, follow the mathematical guide in increasing complexity:
[von Mises distribution](math/von-mises.md), [GLMs](math/glms.md),
[GLMMs](math/glmms.md), and [mixture extensions](math/mixtures.md).
Each section defines the notation alongside its biological interpretation.

## From a biological question to a GLM

A **generalized linear model (GLM)** relates a measured outcome to explanatory
variables, called **predictors**. Here the outcome is an angle, and predictors
could include cue presence, neighbour density, or their interaction: whether
the response to the cue depends on density.

The **von Mises likelihood** describes the distribution of headings around a
mean direction, `mu`, with concentration `kappa`. Larger `kappa` means tighter
clustering; near zero, headings are nearly uniform around the circle. This
respects the closeness of 359° and 1° and keeps changes in preferred direction
distinct from changes in directional consistency.

Each parameter can have its own **linear predictor**, a sum of an intercept
and predictor effects. The **intercept** sets the baseline when predictors are
zero or at their reference category; a **slope** describes change per unit of
a continuous predictor. For a single predictor `x`, the model used in the
[README quick start](https://github.com/jgraving/unwrap#quick-start-circular-regression-in-python)
has this structure:

```text
mean predictor          = mean intercept + mean slope * x
concentration predictor = concentration intercept + concentration slope * x

mu    = circular_modulo(mean predictor)
kappa = softplus(concentration predictor)
y     follows a von Mises distribution with parameters mu and kappa
```

The **inverse links** convert the linear predictors into valid distribution
parameters. `circular_modulo` wraps a predicted direction onto the circle;
`softplus` maps a predictor to positive concentration. Angles use radians.
The [link guide](links.md) gives the functions and numerical conventions.

A basic GLM puts predictors on the mean while keeping other distribution
parameters constant. Allowing predictors to influence concentration as well
is **distributional regression**. We use *circular GLM* for this broader
framework, in which mean direction and concentration can each have predictors.

## Express the model with Bambi or PyMC

`unwrap` works within the Python Bayesian modelling ecosystem. **PyMC** lets
you construct models explicitly from distributions, predictors, and priors,
then estimate their posterior distributions. **Bambi** builds PyMC models
from regression formulas, providing a high-level API for GLMs and GLMMs.
`unwrap` adds circular link objects and prior helpers to these interfaces.

For example, the following Bambi formula describes two regressions fitted
together to the same angular observations:

```python
import bambi as bmb

formula = bmb.Formula("heading ~ cue", "kappa ~ cue")
```

`heading` and `cue` are columns in your data: headings in radians and a
treatment indicator or category. `heading ~ cue` asks whether the cue changes
mean direction; `kappa ~ cue` asks whether it changes directional consistency.
Pass the formula to `bmb.Model` with `family="vonmises"`, `unwrap`'s
`circular_modulo` and `softplus` links, and suitable priors, as in the README.
The formula specifies the relationships; it does not fit a model on its own.

The [API guide](api.md) maps these model components to imports. See Bambi's
[model building guide](https://bambinos.github.io/bambi/notebooks/how_bambi_works.html)
for formula syntax and distributional models, and the
[PyMC documentation](https://www.pymc.io/projects/docs/en/stable/)
for explicit model construction.

## Repeated observations: GLMMs

If you measure each animal several times, the model should represent which
observations belong to the same animal. A **generalized linear mixed model
(GLMM)** extends a GLM with animal- or group-specific effects. These are also
called **hierarchical** or **multilevel** models.

Population-level effects describe the average relationship with a predictor.
Animal-specific effects describe departures from that relationship, such as
a different baseline heading or response to the cue. These are often called
**fixed effects** and **random effects**, respectively.

For data with an `animal` identifier, a Bambi formula could be:

```python
formula = bmb.Formula(
    "heading ~ cue + (1 + cue | animal)",
    "kappa ~ cue + (1 | animal)",
)
```

`(1 + cue | animal)` allows baseline direction and the cue response to vary
among animals. `(1 | animal)` allows baseline concentration to vary among
animals. Together, they distinguish how consistently each animal points from
how closely different animals agree. This design assumes the cue varies
within animals so their individual responses can be estimated; choose the
group structure to match the experiment.

Animal effects share population distributions, so estimates for an animal
with few observations are informed by the others. This is **partial pooling**:
it retains individual differences while estimating the population response.
Grouping by animal represents shared individual differences; closely spaced
trajectory observations may also need a model of dependence through time.

## Interpret effects on the biological scale

With circular modulo, a mean-direction slope describes an angular change in
radians per unit of a continuous predictor before wrapping. For a categorical
predictor, coefficients describe angular shifts relative to its reference
category. With softplus, concentration coefficients act before the link, so
their size is not a constant change in `kappa`. Compare predicted directions
and concentrations under the conditions of interest to interpret both.

Concentration describes variation among observations, not certainty about an
estimated mean. A population can also contain individually consistent animals
with different mean headings. Use predictions at both animal and population
levels to distinguish these patterns.

Bayesian fitting combines priors with the observed angles to estimate a
posterior distribution. Check sampling diagnostics, summarise angles on the
circle, and compare simulated headings from the fitted model with the data
using **posterior predictive checks**. The
[circular regression tutorial](https://github.com/jgraving/unwrap/blob/main/notebooks/circular_regression.ipynb)
demonstrates fitting and checking a GLM; the [biological examples](examples.md)
explore repeated observations and richer directional patterns.
