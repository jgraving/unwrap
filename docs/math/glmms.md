# 3. Circular GLMMs

Repeated observations can reveal whether each animal keeps a consistent
heading and whether different animals agree. A generalized linear mixed
model (GLMM), also called a hierarchical or multilevel model, extends the
[GLM](glms.md) with animal- or group-specific departures from population effects.

## Add individual effects to the regression coefficients

Let $i$ identify an animal and $n$ an observation from that animal. The
coefficient vectors are

$$
\boldsymbol{\beta}_i=\boldsymbol{\beta}+\mathbf{b}_i,
\qquad
\boldsymbol{\alpha}_i=\boldsymbol{\alpha}+\mathbf{a}_i.
$$

$\boldsymbol{\beta}$ and $\boldsymbol{\alpha}$ are population-level
coefficients. $\mathbf{b}_i$ and $\mathbf{a}_i$ are animal-specific offsets for
mean direction and concentration. An offset can change a baseline
intercept, a treatment response, or a continuous slope. The observation model is

$$
\begin{aligned}
\eta_{\mu,in}&=\mathbf{x}_{in}^\top(\boldsymbol{\beta}+\mathbf{b}_i),\\
\eta_{\kappa,in}&=\mathbf{x}_{in}^\top(\boldsymbol{\alpha}+\mathbf{a}_i),\\
\mu_{in}&=g_\mu^{-1}(\eta_{\mu,in}),\\
\kappa_{in}&=g_\kappa^{-1}(\eta_{\kappa,in}),\\
y_{in}\mid\mathbf{x}_{in},\boldsymbol{\beta}_i,\boldsymbol{\alpha}_i
&\sim\mathcal{M}(\mu_{in},\kappa_{in}).
\end{aligned}
$$

The inverse links are the same circular-modulo and softplus functions used in
the GLM. The sampling unit can also be a trial, site, or social group; choose
the grouping and varying effects to match the experiment.

## Partial pooling on the predictor scale

For offsets on an unconstrained scale before the links, a common construction is

$$
\mathbf{b}_i\sim\mathcal{N}(\mathbf{0},\Sigma_\beta),
\qquad
\mathbf{a}_i\sim\mathcal{N}(\mathbf{0},\Sigma_\alpha).
$$

The covariance matrices $\Sigma_\beta$ and $\Sigma_\alpha$ describe the spread
of animal effects and any modelled association among their intercepts and
slopes. Their unknown scales receive **hyperpriors**, priors on the parameters
of these shared distributions. Learning the spread from all animals induces
**partial pooling**: individual estimates share information while retaining
differences supported by the observations.

These multivariate Normal priors are one predictor-scale construction.
An angular offset can instead use the circular construction below.

## Circular offsets: between-animal alignment

For a genuinely angular effect, such as an animal's baseline directional
offset, a circular alternative is

$$
b_{i,0}\sim\mathcal{M}(0,\kappa_\mu),
\qquad
\mu_i(\mathbf{x}=\mathbf{x}_0)
=g_\mu^{-1}(\beta_0+b_{i,0}),
$$

where $\mathbf{x}_0=(1,0,\ldots,0)^\top$ is the baseline design vector.
This replaces the Normal prior for that angular offset. The remaining
offsets can still use appropriate predictor-scale distributions.

$\kappa_\mu$ measures how tightly animal directional offsets cluster around
the population direction. Give it a hyperprior when estimating this alignment
from the data. It is different from the observation concentration
$\kappa_{in}$:

| Quantity | Biological interpretation |
| --- | --- |
| $\kappa_{in}$ | Consistency of headings within an animal at the specified predictors |
| $\kappa_\mu$ | Alignment of animal mean-direction offsets around the population value |

Thus, animals can each hold a steady heading (high observation concentration)
while pointing in different directions (low between-animal alignment).
Angular treatment offsets can be modelled similarly when supported by the
experimental design.

Bambi's standard group-effect formulas, such as `(1 + cue | animal)`, specify
group effects on the predictor scale. They do not by themselves replace an
offset prior with a von Mises distribution. Circular offset constructions can
be specified explicitly in PyMC to match the intended model.

## Choose the level of prediction

The same fitted model supports different biological summaries:

| Target | Calculation |
| --- | --- |
| An observed animal | Use its population coefficients plus estimated offsets |
| The population relationship at zero offsets | Set $\mathbf{a}_i=\mathbf{b}_i=\mathbf{0}$ before applying the links |
| A new animal from the population | Draw new offsets from the fitted population distributions and propagate their uncertainty |

Setting offsets to zero describes a reference animal. It generally differs
from the prediction averaged over animal variation, because the inverse
links are nonlinear. In particular, observation concentration and between-animal
dispersion both shape the distribution of headings pooled across animals.
Posterior predictive checks at both levels help distinguish these patterns.

**Next:** [Extend the observation model with mixtures](mixtures.md).
