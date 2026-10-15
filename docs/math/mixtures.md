# 4. Mixture extensions

A single von Mises distribution has one preferred direction. Some data
instead combine directed and undirected headings, or contain competing
preferred directions. Mixtures extend the [GLM](glms.md) and [GLMM](glmms.md)
framework by combining observation densities with estimated weights.

These equations follow *Methods: Mixture likelihoods for circular outcomes*
and the *Mixture von Mises Likelihoods* and *Mixture weights* appendices in the
[paper](https://github.com/jgraving/circstats/blob/438d2f2e2880defb792bc8d4e2ecefbc7337470d/elife-template.tex).
Here $\mathcal{M}$ denotes the [von Mises density](von-mises.md), $n$ indexes
observations, and $i$ identifies an animal.

## Directed and undirected headings: von Mises–uniform

The von Mises–uniform mixture (vMU) combines a directional component and a
uniform component:

$$
p(y_n\mid\Theta)
=\omega_n\,\mathcal{M}(y_n\mid\mu_n,\kappa_n)
+(1-\omega_n)\frac{1}{2\pi}.
$$

$\omega_n$ is the probability that an observation comes from the directed
component; $1-\omega_n$ is the probability of the uniform component.
$\kappa_n$ describes consistency **within the directed component**. The model
can therefore separate a change in the prevalence of directed headings from
a change in their precision, rather than making one concentration absorb both.

Use the same circular-modulo and softplus links for $\mu_n$ and $\kappa_n$.
For the weight, define a linear predictor with coefficients $\boldsymbol{\gamma}$:

$$
\begin{aligned}
\eta_{\omega,n}&=\mathbf{x}_n^\top\boldsymbol{\gamma},\\
\omega_n&=\operatorname{logit}^{-1}(\eta_{\omega,n})
=\frac{1}{1+\exp(-\eta_{\omega,n})}.
\end{aligned}
$$

The inverse-logit, also called the sigmoid, maps an unconstrained prediction
to a probability between zero and one. Its forward link is
$\operatorname{logit}(\omega)=\log[\omega/(1-\omega)]$.
Predictors can change mean direction, directed-component concentration, and
the directed fraction separately.

### Which component could have generated an observation?

Let $w_n=1$ denote the directed component and $w_n=0$ the uniform component.
For a particular parameter draw $\Theta$,

$$
P(w_n=1\mid y_n,\Theta)
=\frac{\omega_n\mathcal{M}(y_n\mid\mu_n,\kappa_n)}
{\omega_n\mathcal{M}(y_n\mid\mu_n,\kappa_n)+(1-\omega_n)/(2\pi)}.
$$

This **responsibility** is a probability of component membership under the
model. Evaluate it for each posterior draw and summarise across draws to
include parameter uncertainty. A component's biological interpretation
depends on the experiment; membership alone does not establish a behavioural
mechanism. If $\kappa_n$ is near zero, the directed component also becomes
nearly uniform, making its weight difficult to distinguish from concentration.

## Two competing directions: a two-component von Mises mixture

For two directional components,

$$
\begin{aligned}
p(y_n\mid\Theta)
={}&\omega_n\,\mathcal{M}(y_n\mid\mu_{n1},\kappa_{n1})\\
&+(1-\omega_n)\,\mathcal{M}(y_n\mid\mu_{n2},\kappa_{n2}).
\end{aligned}
$$

Each component has a mean direction and concentration, with its own predictor
effects if needed. The weight uses the same inverse-logit mapping as above.
Two components need not produce two visibly separate peaks if they overlap;
they describe candidate directional processes, and the fitted density should
be checked against the observations. The mean of the whole mixture can lie
between its preferred directions or be undefined when their effects cancel.
Inspect the predictive density and component-specific predictions as well.

## Add animal differences to a mixture

Mixture weights can vary among animals alongside component means and
concentrations. For example, a hierarchical weight model is

$$
\begin{aligned}
\eta_{\omega,in}&=\mathbf{x}_{in}^\top\boldsymbol{\gamma}+c_i,\\
c_i&\sim\mathcal{N}(0,\tau_\omega^2),\\
\omega_{in}&=\operatorname{logit}^{-1}(\eta_{\omega,in}).
\end{aligned}
$$

$c_i$ is an animal-specific offset and $\tau_\omega$ its population standard
deviation, with a suitable hyperprior. Component means and concentrations
can have the same partially pooled coefficients described for [GLMMs](glmms.md).

Repeated observations help distinguish animals that each mix between two
directions from animals that consistently choose different directions.
Priors on the weight offsets matter: narrow distributions favour similar
weights among animals, while wider distributions allow weights close to
zero or one. Check what those priors imply on the probability scale.
This mixture specifies component probabilities for observations; temporal
persistence or transitions between states require additional structure.

## More than two components

For $K$ directional components, the general model is

$$
p(y_n\mid\Theta)
=\sum_{k=1}^{K}\pi_{nk}\,\mathcal{M}(y_n\mid\mu_{nk},\kappa_{nk}),
\qquad \sum_{k=1}^{K}\pi_{nk}=1.
$$

The weights $\pi_{nk}$ are positive probabilities. Each component can use
regressions for its mean and concentration, with the same links as before:

$$
\begin{aligned}
\mu_{nk}&=g_\mu^{-1}(\mathbf{x}_n^\top\boldsymbol{\beta}_k),\\
\kappa_{nk}&=g_\kappa^{-1}(\mathbf{x}_n^\top\boldsymbol{\alpha}_k).
\end{aligned}
$$

The **softmax** inverse link maps weight predictors to probabilities that sum
to one:

$$
\eta_{\pi,nk}=\mathbf{x}_n^\top\boldsymbol{\gamma}_k,
\qquad
\pi_{nk}=\frac{\exp(\eta_{\pi,nk})}
{\sum_{j=1}^{K}\exp(\eta_{\pi,nj})}.
$$

Adding the same constant to every weight predictor leaves all weights
unchanged. Fix one reference predictor, such as $\eta_{\pi,nK}=0$, or impose
a sum-to-zero constraint to remove this redundancy. This is separate from
the component-label issue below. For numerical stability, subtract the
largest weight predictor before evaluating the exponentials; the weights
are unchanged. The paper's worked applications focus on
two components; additional components require more data and careful evaluation.

## Component labels and interpretation

With exchangeable von Mises components, swapping component labels leaves the
mixture density unchanged. This **label switching** can make component-wise
posterior averages misleading. Use scientifically justified identifying
structure or appropriate posterior relabelling, and check the predictive
density. Ordering angles requires an explicit circular reference convention.

Fixing one component to uniform removes the vMU model's label symmetry, but
does not prevent weak identification when its components look similar.
For either mixture, compare predictions with observed angular distributions
at both animal and population levels before assigning biological meaning.

## Python implementation

Specify mixture structure explicitly in PyMC. The
[`VonMisesMixture` helper](../api.md#distributions-unwrapdistributions) combines
von Mises components; component regressions, weight links, priors, and
identifying structure belong in the model specification. The vMU equation
also requires a uniform component. These mixture equations are not selected
by Bambi's single-component `family="vonmises"` setting.

See the [biological examples](../examples.md) for the locust mixture analysis
and its data and execution requirements.
