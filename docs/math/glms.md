# 2. Circular GLMs

A treatment might turn animals toward a cue, make their headings more
consistent, or do both. A circular generalized linear model (GLM) describes
these as separate relationships with predictors, using the
[von Mises distribution](von-mises.md) for the observed angles.

The equations follow *Models: A Generalized Linear Model for Circular Data*
in the [paper](https://github.com/jgraving/circstats/blob/438d2f2e2880defb792bc8d4e2ecefbc7337470d/elife-template.tex).

## Two regression equations, one angular outcome

Let $n=1,\ldots,N$ index observations. The predictor vector
$\mathbf{x}_n=(1,x_{n1},\ldots,x_{nM})^\top$ includes an intercept and $M$
predictors, such as treatment, density, or an interaction. Define

$$
\begin{aligned}
\eta_{\mu,n}
&=\mathbf{x}_n^\top\boldsymbol{\beta}
=\beta_0+\sum_{m=1}^{M}\beta_m x_{nm},\\
\eta_{\kappa,n}
&=\mathbf{x}_n^\top\boldsymbol{\alpha}
=\alpha_0+\sum_{m=1}^{M}\alpha_m x_{nm}.
\end{aligned}
$$

$\boldsymbol{\beta}$ contains effects on the mean-direction predictor;
$\boldsymbol{\alpha}$ contains effects on the concentration predictor.
The $\eta$ terms are **linear predictors**: intercepts plus weighted sums of
the explanatory variables, before conversion to distribution parameters.
Each equation can contain a different set of predictors; concentration can
also be held constant.

## Map predictions onto the circle and positive values

The paper's inverse links give

$$
\begin{aligned}
\mu_n
&=g_\mu^{-1}(\eta_{\mu,n})
=\bigl((\eta_{\mu,n}+\pi)\bmod 2\pi\bigr)-\pi,\\
\kappa_n
&=g_\kappa^{-1}(\eta_{\kappa,n})
=\log\!\left(1+\exp(\eta_{\kappa,n})\right),\\
y_n\mid\mathbf{x}_n,\boldsymbol{\beta},\boldsymbol{\alpha}
&\sim\mathcal{M}(\mu_n,\kappa_n).
\end{aligned}
$$

The circular-modulo operation wraps the predicted angle to $[-\pi,\pi)$.
Softplus maps the concentration predictor to positive values. These are
conventionally called inverse links; circular wrapping is many-to-one and
has no global mathematical inverse. The [link reference](../links.md)
describes the numerical functions and boundary behaviour.

Including predictors for both parameters is **distributional regression**:
the experiment can explain directional consistency as well as mean direction.

## Example: a cue can affect direction and consistency separately

For cue absence $x=0$ and presence $x=1$,

$$
\begin{aligned}
\mu(x)&=g_\mu^{-1}(\beta_0+\beta_1x),\\
\kappa(x)&=g_\kappa^{-1}(\alpha_0+\alpha_1x).
\end{aligned}
$$

$\beta_1$ is the treatment's angular shift in radians before wrapping. For a
continuous predictor, it is the change per predictor unit. $\alpha_1$ is a
change on the pre-softplus scale, so its biological size is best expressed as
the resulting difference in predicted concentrations. For either effect,
derive predictions and contrasts within each posterior draw to preserve
uncertainty.

With the specified links, the corresponding Bambi formula is
`bmb.Formula("heading ~ cue", "kappa ~ cue")`. The
[conceptual guide](../models.md#express-the-model-with-bambi-or-pymc) connects
this formula to Bambi and PyMC, and the
[README quick start](https://github.com/jgraving/unwrap#quick-start-circular-regression-in-python)
shows a complete model specification with priors.

## Bayesian fitting

Let $\Theta$ collect the unknown parameters and
$\mathcal{D}=\{(\mathbf{x}_n,y_n)\}_{n=1}^N$ be the data. Assuming observations
are conditionally independent given the model parameters,

$$
\begin{aligned}
p(\mathcal{D}\mid\Theta)
&=\prod_{n=1}^{N}\mathcal{M}(y_n\mid\mu_n,\kappa_n),\\
p(\Theta\mid\mathcal{D})
&\propto p(\mathcal{D}\mid\Theta)\,p(\Theta).
\end{aligned}
$$

Priors $p(\Theta)$ describe plausible parameters before fitting; the posterior
describes uncertainty after combining those priors with the observations.
Check priors by simulating headings before fitting, then check sampling
diagnostics and posterior predictions afterward. For repeated observations,
the [GLMM extension](glmms.md) adds shared animal or group effects; additional
temporal dependence may require further structure.

**Next:** [Represent individual differences with a GLMM](glmms.md).
