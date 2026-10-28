# Circular convergence diagnostics

`unwrap.diagnostics.rhat` and `summary` can check a whole posterior while
treating named angular variables as circular. Pass their names in
`circ_var_names`; other parameters retain ordinary ArviZ diagnostics. Without
circular names, the existing ArviZ behavior is preserved.

For circular variables, the default `circ_method="circular"` selects a joint
within/between-chain covariance diagnostic. This is a different statistic from
modern rank-normalized R-hat. The summary also reports separate standard ArviZ
diagnostics for `sin(theta)` and `cos(theta)`, so both projections remain visible.
`rhat_components` provides those two rank R-hats directly.

The joint circular method compares within-chain and between-chain covariance of
the joint unit-vector representation

$$
u(\theta)=\begin{pmatrix}\cos\theta\\\sin\theta\end{pmatrix}.
$$

Directions separated by a full turn have the same representation. A change
of angular origin rotates every vector together, so there is no arbitrary
cut at which the representation jumps. Both coordinates enter one covariance
calculation, rather than two separate scalar R-hats.

## Usage and relationship to ArviZ

```python
import arviz as az
import numpy as np
from unwrap.diagnostics import rhat, rhat_components, summary

rng = np.random.default_rng(42)
theta = rng.vonmises(mu=0.5, kappa=8, size=(4, 1000))

# Standard ArviZ rank-normalized split R-hat; behavior is unchanged.
standard = rhat(theta)
np.testing.assert_equal(standard, az.rhat(theta))

# Report each angular projection using ArviZ's standard rank method.
components = rhat_components(theta)
sin_rhat = components["sin"]
cos_rhat = components["cos"]

# Explicitly select the joint circular covariance diagnostic.
circular = rhat(theta, method="circular")
```

`rhat` extends the ArviZ 0.x signature with circular-variable selection:

```python
rhat(
    data, *, var_names=None, method="rank", dask_kwargs=None,
    circ_var_names=None, circ_method="circular",
)
rhat_components(data, *, var_names=None, dask_kwargs=None)
```

| Input | Circular-method interpretation and output |
| --- | --- |
| Two-dimensional NumPy array | Axes are `(chain, draw)`; returns one scalar |
| xarray Dataset | Uses named `chain` and `draw` dimensions; returns a Dataset preserving other variable dimensions and coordinates |
| xarray DataArray | Converted to a Dataset for the same named-dimension calculation |
| ArviZ InferenceData or other convertible input | Uses the posterior group through ArviZ's conversion and selection machinery |

For a von Mises posterior, name angular parameters and include all other
parameters in the same call. Use the actual names in your posterior:

```python
circular_names = ["Intercept", "animal_offset", "mu"]
checks = rhat(idata, circ_var_names=circular_names)
table = summary(idata, circ_var_names=circular_names)

# Use a particular rank-normalized projection for the selected circular values.
sin_checks = rhat(idata, circ_var_names=circular_names, circ_method="sin")
cos_checks = rhat(idata, circ_var_names=circular_names, circ_method="cos")

# Or request both projection R-hats directly for only the angular variables.
direction_components = rhat_components(idata, var_names=circular_names)
```

Here, a continuous-predictor slope or concentration `kappa` remains in the
results with ordinary diagnostics. The call does not apply circular methods to
every parameter merely because the likelihood is von Mises.

`circ_var_names` contains exact variable names. Names are validated against the
full relevant dataset; a name excluded by `var_names` is then ignored. Unknown
names raise an error. `var_names` retains ArviZ selection conventions, and
`dask_kwargs` is forwarded for the R-hat calculations. For an unnamed NumPy
angular array, `rhat(theta, method="circular")` remains the direct route to the
joint method. Selecting `method="circular"` applies it to every selected
variable, so select only angles in that mode. Do not combine that mode with
`circ_var_names`; use `circ_method="circular"` for a mixed dataset.

For mixed data, `method` controls ordinary variables and `circ_method` controls
the named circular variables. The latter accepts `"circular"`, `"sin"`, or
`"cos"`; the projection choices use ArviZ's rank-normalized split R-hat.
The wrapper does not infer angular variables or convert degrees. All angular
inputs must be radians.

## A summary for mixed circular and ordinary parameters

`summary` accepts the same arguments as `arviz.summary`, including
`var_names`, `filter_vars`, `group`, `coords`, `kind`, `fmt`, `round_to`,
`hdi_prob`, `skipna`, and custom `stat_funcs`. It adds
`circ_method="circular"` and `include_circular_diagnostics=True`:

```python
table = summary(
    idata,
    circ_var_names=["Intercept", "animal_offset", "mu"],
    circ_method="circular",
    include_circular_diagnostics=True,
    round_to=3,
)

# The same choices also work for a labeled result or diagnostics-only table.
labeled = summary(idata, circ_var_names=["mu"], fmt="xarray")
diagnostics = summary(idata, circ_var_names=["mu"], kind="diagnostics")
```

Without circular names, the call delegates to ArviZ unchanged. With circular
names, ordinary rows retain ArviZ's statistics and diagnostics. Circular rows
use circular mean, circular standard deviation, and circular HDI endpoints.
A circular interval that crosses the angular cut can have a numerically larger
first endpoint than second endpoint. Circular summaries require
`stat_focus="mean"`; median-focused circular summaries are not supported.

The diagnostic columns have explicit meanings:

| Column on a circular row | Meaning |
| --- | --- |
| `r_hat` | The method selected by `circ_method`: joint circular covariance by default, or rank R-hat of sine or cosine |
| `mcse_mean` | ArviZ's circular-mean Monte Carlo error estimate, retained for both `kind="all"` and `kind="diagnostics"` |
| `mcse_sd`, `ess_bulk`, `ess_tail` | NaN: ordinary raw-angle estimates are not presented as circular diagnostics |
| `r_hat_sin`, `r_hat_cos` | Separate ArviZ rank-normalized split R-hats of the projections |
| `ess_bulk_sin`, `ess_bulk_cos`, `ess_tail_sin`, `ess_tail_cos` | ArviZ effective sample sizes for the corresponding projections |
| `mcse_mean_sin`, `mcse_mean_cos`, `mcse_sd_sin`, `mcse_sd_cos` | ArviZ Monte Carlo errors for the mean or standard deviation of the corresponding projection |

Projection columns are included by default for `kind="all"` and
`kind="diagnostics"`; their ordinary-parameter rows contain NaN. Set
`include_circular_diagnostics=False` to omit those extra columns. `kind="stats"`
contains only descriptive statistics. Custom statistics retain ArviZ's behavior:
their functions receive the original draws, and `extend=False` requests only
those custom statistics.

These projection errors describe estimates of quantities such as
$E[\sin\theta]$, not errors in the circular mean direction. Likewise, tail ESS
describes the tails of a sine or cosine projection; the circle has no intrinsic
lower and upper tails. The circular `mcse_mean` is ArviZ's batch-based estimate,
not a new joint covariance or rank-based estimator. A mean direction can itself
be unstable or uninformative when the resultant length is near zero.

ArviZ 0.20 and 0.23 already accept `circ_var_names` for descriptive summaries,
but leave `r_hat`, `ess_bulk`, `ess_tail`, and `mcse_sd` based on the raw angles.
They also use circular `mcse_mean` for `kind="all"` but raw-angle `mcse_mean`
for `kind="diagnostics"`. This wrapper makes the circular mean error consistent
between the two modes and replaces or leaves undefined the other raw-angle
diagnostics as described above. Thus adding a full turn to a chain need not
produce an apparent diagnostic failure solely because its stored angles use a
different numerical interval.

The default circular `r_hat` is **not** the same statistic as the rank R-hat
on ordinary rows or in the projection columns. Its derivation and limitations
appear below; a universal threshold of 1.01 has not been calibrated for it.
The [validation experiments](validation/circular_rhat.md) compare the methods
on agreement, disagreement, and known blind spots.

## Reporting sine and cosine R-hat

`rhat_components` applies the two trigonometric transformations to the selected
posterior variables, then calls `arviz.rhat(..., method="rank")` for each.
The result is a dictionary with exactly two keys, `"sin"` and `"cos"`.
Each value follows ArviZ's scalar or Dataset output behavior. Variable selection
happens before the transformations, and the same `dask_kwargs` is forwarded
to both calculations. There is no averaging, maximum, or automatic pass/fail
decision in the public result.

The rank method includes split-chain rank normalization and folding. This
retains ArviZ's handling of invalid input and degenerate projections. For
example, a coordinate that is constant across all draws can produce NaN;
`rhat_components` does not replace that result with one, silently discard the
coordinate, or substitute the other result.

Report both values. An arithmetic average can dilute a large R-hat in one
coordinate with a value near one in the other. A maximum preserves the larger
value, but is still an aggregation with its own interpretation. The package
does not recommend a separately calibrated aggregate cutoff for the two
projections.

The pair $(\cos\theta,\sin\theta)$ is continuous around the circle and unchanged
by full turns. However, separate rank R-hats on its coordinates are **not
rotation invariant**: changing the angular origin changes the projections,
and nonlinear rank transformations do not commute with that rotation. This
is a difference from the joint covariance method derived below.

The two scalar diagnostics also do not test equality of the full joint
angular distributions. For example, equal-weight distributions on
$\{\pi/4,-3\pi/4\}$ and $\{3\pi/4,-\pi/4\}$ have the same sine and cosine
marginals but different pairings of those coordinates. Both component
diagnostics can therefore look satisfactory while the angular distributions
differ. Inspect per-chain angular distributions alongside their diagnostic
values, especially for multimodal posteriors.

## Which model quantities are circular?

Choose diagnostics according to the quantity whose sampling matters, rather
than the variable name alone.

| Quantity | Appropriate use |
| --- | --- |
| Population intercept defining a mean direction at the reference predictor values | Apply angular component or joint circular diagnostics to its directional interpretation |
| Angular random intercept or angular offset | Apply angular diagnostics when it represents a rotation modulo $2\pi$ |
| Subject-specific direction, such as $\operatorname{wrap}(\alpha+u_j)$ | Compose each posterior draw first, then diagnose that subject's direction |
| Predicted mean direction $\mu_j(x)$ | Apply angular diagnostics at the predictor values of scientific interest |
| Slope for a continuous predictor | Use ordinary coefficient diagnostics; adding $2\pi$ generally changes predictions at intermediate predictor values |
| Concentration $\kappa$, random-effect scale, or other non-angular parameter | Use ordinary ArviZ diagnostics |

It is useful to check both angular offsets and the resulting subject directions:
the latter incorporates their posterior dependence with the population
direction. A real-valued intercept or offset with a nonperiodic prior can also
be assessed on its original coefficient scale when that is a quantity of
interest. Angular diagnostics assess its direction and do not change its
prior or sampling parameterization.

For example, if an xarray posterior contains scalar `Intercept` draws and
angular `animal_offset` draws indexed by animal, one can form and assess the
resulting directions without wrapping them first:

```python
subject_direction = (
    idata.posterior["Intercept"] + idata.posterior["animal_offset"]
).rename("subject_direction")
subject_checks = rhat_components(subject_direction)

# Keep the original parameters and add the composed direction to one summary.
posterior = idata.posterior.assign(subject_direction=subject_direction)
table = summary(
    posterior,
    circ_var_names=["Intercept", "animal_offset", "subject_direction"],
)
```

Sine and cosine perform the periodic identification. Neither function aligns
different chains separately, changes component labels, or infers which
random effects are angular.

## Joint method: circular within-chain and between-chain dispersion

For $m$ equally long chains with $n$ draws, write their complex mean vectors as

$$
z_j=\frac1n\sum_{t=1}^n e^{i\theta_{jt}},\qquad
r_j=|z_j|,\qquad
\bar z=\frac1m\sum_{j=1}^m z_j.
$$

The usual circular variance is one minus the mean resultant length. It has
the following exact within-group and between-group decomposition:

$$
T=1-|\bar z|,\qquad
W_{\rm circ}=\frac1m\sum_j(1-r_j),\qquad
B_{\rm circ}=\frac1m\sum_j r_j-|\bar z|,
$$

$$
T=W_{\rm circ}+B_{\rm circ}.
$$

The triangle inequality makes $B_{\rm circ}$ nonnegative. Thus it is possible
to compare within-chain and between-chain circular variation directly, as
in circular analysis of variance. See the directional decomposition of
[Anderson and Wu](https://uwaterloo.ca/business-and-industrial-statistics-research-group/sites/default/files/uploads/documents/rr-92-10.pdf).

This particular between-group term assesses alignment of mean directions.
If all chain mean vectors point in the same direction, then
$|\bar z|=m^{-1}\sum_j r_j$ and $B_{\rm circ}=0$, even if their lengths, and
therefore concentrations, differ. A ratio based on this decomposition alone
does not detect that difference.

The unit-vector embedding also gives a squared-chordal decomposition:

$$
T_2=1-|\bar z|^2,\qquad
W_2=\frac1m\sum_j(1-|z_j|^2),\qquad
B_2=\frac1m\sum_j|z_j-\bar z|^2,
$$

$$
T_2=W_2+B_2.
$$

Here the chain mean vectors lie inside the unit disk; they are not normalized
back onto the unit circle. Consequently $B_2$ detects differences in both
their lengths and directions. This identity follows by expanding the squared
Euclidean distances and is the ordinary variance decomposition of the embedded
vectors. These dispersion identities explain the construction, but the
implemented method uses full covariance matrices with sample-size corrections,
not the scalar ratio $T_2/W_2$.

## Joint method: the implemented statistic

The circular method first splits each original chain in half. If there are
$M$ original chains and $N$ draws, the calculation uses
$m=2M$ chains of $n=\lfloor N/2\rfloor$ draws. The first and last $n$ draws
are retained; an odd middle draw is omitted. Splitting makes changes between
the early and late parts of a chain contribute to between-chain disagreement.
It does not remove warmup; supply post-warmup draws.

For the split chains define

$$
\bar u_j=\frac1n\sum_t u_{jt},\qquad
\bar u=\frac1m\sum_j\bar u_j,
$$

$$
W=\frac1m\sum_{j=1}^m\left[
  \frac1{n-1}\sum_{t=1}^n
  (u_{jt}-\bar u_j)(u_{jt}-\bar u_j)^\mathsf{T}
\right],
$$

$$
S=\frac1{m-1}\sum_{j=1}^m
(\bar u_j-\bar u)(\bar u_j-\bar u)^\mathsf{T}.
$$

Both sample covariances use `ddof=1`. $W$ is the average within-chain
covariance; $S$ is the covariance of the chain means. In the conventional
notation $B=nS$. The corrected covariance estimate used here is

$$
\widehat V=\frac{n-1}{n}W+\frac{m+1}{m}S.
$$

For positive-definite $W$, let $\lambda_{\max}(S,W)$ denote the largest
generalized eigenvalue, satisfying $Sa=\lambda Wa$. The returned value is

$$
\widehat R_{\rm circular}
=\sqrt{\frac{n-1}{n}+\frac{m+1}{m}\lambda_{\max}(S,W)}
=\sqrt{\max_{a\ne0}\frac{a^\mathsf{T}\widehat V a}{a^\mathsf{T}Wa}}.
$$

It checks the most discrepant linear projection of the two-dimensional
embedding. Every such projection has the form
$a_1\cos\theta+a_2\sin\theta=A\cos(\theta-\alpha)$.

This follows the multivariate covariance construction in
[Brooks and Gelman (1998), Section 4.1](https://sites.stat.columbia.edu/gelman/research/published/brooksgelman2.pdf).
Their published $\widehat R^p$ is the **variance ratio without the square root**.
This implementation returns $\sqrt{\widehat R^p}$ to put the result on a
standard-deviation scale. The squared return value recovers that published
ratio. The finite-chain factor is $(m+1)/m$, where $m$ counts the split chains,
and is not a factor based on the two embedding coordinates. No scalar
degrees-of-freedom correction is applied.

This is a split, classical multivariate covariance diagnostic. It is not
rank-normalized or folded, and does not claim the robustness properties of
ArviZ's modern `method="rank"`. Ranking the two coordinates separately would
change the joint geometry and lose the rotation invariance described below.

## Joint method: invariance and the concentrated-angle limit

Adding any integer multiple of $2\pi$ to any draw leaves $u(\theta)$
unchanged. Under a common rotation or a reversal of angular orientation,
$u$ changes to $Qu$ for an orthogonal matrix $Q$. The covariance matrices
then change to $QWQ^\mathsf{T}$ and $QSQ^\mathsf{T}$. Their generalized
eigenvalues are unchanged, proving that the diagnostic is independent of
the angular origin, the chosen interval, and clockwise versus anticlockwise
orientation in exact arithmetic.

The same covariance argument permits a common numerical rotation, translation,
or nonsingular rescaling of the embedding. It does not permit independently
aligning each chain: that could remove actual differences between chains.
Period and rotation invariance are approximate in floating-point arithmetic,
especially for enormous angles whose directional information is already
limited by numerical precision.

For concentrated draws, write $\theta=\phi+\delta$ with small $\delta$ and
rotate the embedding by $-\phi$:

$$
u(\delta)=
\begin{pmatrix}1-\delta^2/2+O(\delta^4)\\
\delta+O(\delta^3)\end{pmatrix}.
$$

The tangent projection is approximately $\delta$. Its within/between ratio
therefore approaches the ordinary scalar covariance PSRF for the local angular
deviations, with the same finite-chain correction. The joint statistic maximizes
over every projection, so it is at least as large as that tangent-projection
value. It need not equal it: the radial coordinate also contains information
about $\delta^2$, and can respond to different concentrations. This relationship
is to a classical scalar PSRF, not rank-normalized R-hat.

## Joint method: invalid and degenerate inputs

The circular method requires at least two original chains and four original
draws per chain. Insufficient samples or nonfinite values give NaN with a
diagnostic warning. There is no automatic removal of missing draws because
that would alter the chains' time ordering and effective lengths.

If all draws represent one constant direction, including representations
separated by whole turns, the result is NaN: no within-chain variation is
available to assess mixing. If every chain is constant but the chains represent
different directions, the result is infinity. Other singular within-chain
covariance matrices give NaN with a warning. Singular cases are not silently
regularized, converted to a pseudoinverse result, or reported as convergence.

Very concentrated draws can make the radial variance much smaller than the
tangent variance and leave $W$ numerically singular. This is a limitation of
the covariance calculation, not proof that the model has failed to converge.
Conversely, a mean resultant near zero does not itself make the calculation
undefined: a broad circular distribution can have a well-conditioned $W$
without having a meaningful single mean direction.

Finite-sample values slightly below one are possible and are not clamped.

## Joint method: interpretation and limits

An elevated value indicates that at least one angular projection has greater
between-chain disagreement than the within-chain covariance accounts for.
Agreement about a direction is assessed modulo a full turn. The usual
rank-R-hat threshold of 1.01 has not been calibrated as a universal cutoff for
this circular method, and a value close to one does not establish convergence.

The method uses differences between **mean vectors**. Covariance determines
how those differences are scaled; it is not a separate test that all chains'
covariance matrices or full angular distributions agree. For example, chains
exploring equal-weight antipodal modes at $0$ and $\pi$ can have the same zero
mean vector as chains exploring modes at $\pi/2$ and $3\pi/2$. The pooled
within-chain covariance can be nonsingular while the diagnostic is close to
one, despite different angular distributions. More generally, distributions
can share their first trigonometric moment and still differ substantially.

Use the circular result alongside per-chain angular distributions, traces,
ordinary sampler diagnostics, and diagnostics of the scientific quantities
of interest. For mixtures, consistent component identities require separate
attention; a circular diagnostic does not resolve label switching. Existing
ArviZ methods remain available unchanged through this wrapper.

## References

- Brooks, S. P. and Gelman, A. (1998). *General Methods for Monitoring
  Convergence of Iterative Simulations*. Journal of Computational and Graphical
  Statistics, 7(4), 434–455. [Author-hosted paper](https://sites.stat.columbia.edu/gelman/research/published/brooksgelman2.pdf).
- Anderson, C. M. and Wu, C. F. J. (1992). *Decomposition of Variation in
  Directional Data*. IIQP Research Report RR-92-10.
  [University of Waterloo report](https://uwaterloo.ca/business-and-industrial-statistics-research-group/sites/default/files/uploads/documents/rr-92-10.pdf).
- Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., and Bürkner, P.-C.
  (2021). *Rank-normalization, folding, and localization: An improved R-hat
  for assessing convergence of MCMC*. Bayesian Analysis, 16(2), 667–718.
  [Paper](https://doi.org/10.1214/20-BA1221).
- [ArviZ 0.x `rhat` API](https://python.arviz.org/en/v0.20.0/api/generated/arviz.rhat.html).
- [ArviZ 0.x `summary` API](https://python.arviz.org/en/v0.20.0/api/generated/arviz.summary.html).
