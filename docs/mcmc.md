# From angles to inference

If a beetle heads 20° further clockwise than usual, did its preferred direction
change, or do its headings vary more from trial to trial? Mean direction
$\mu$ and concentration $\kappa$ describe different features of the same
observations. The explorer fits both and shows which combinations the data
support.

Select **Single distribution** and press **Play sampling**. The dots are
simulated headings. The dashed curve is the density that generated them.
Four Hamiltonian Monte Carlo chains explore the posterior given just the
observations and the model's priors. Each chain uses 120 warmup iterations
to tune its step size, then retains 1,000 draws.

[Open the explorer at full width](assets/mcmc/explorer.html)

<iframe src="../assets/mcmc/explorer.html" title="Interactive circular MCMC explorer with six models" loading="lazy" style="width:100%;height:760px;border:0"></iframe>

## Read the displays

- **Posterior on the circle.** Each sampled pair $(\mu,\kappa)$ defines a
  density for a future heading. The solid curve averages those densities.
  At each angle, the shaded band spans the middle 90% of the sampled
  densities. The band describes uncertainty in the height of the curve.
- **Posterior geometry.** A point in the left panel gives the mean and
  concentration predictors for one draw. The right panel shows the same
  draw after wrapping the angle and mapping concentration to a positive
  value. Chains one complete turn apart on the left can describe the
  same heading on the right.
- **Chains.** A trace follows one parameter through successive iterations.
  A rejected proposal repeats the preceding state. Shading marks warmup.
  Compare the chains in both coordinate systems. A jump from +179° to
  −179° in a constrained trace is a 2° change on the circle. Open **Model
  details** for acceptance, divergences and retained draws.

## Try three comparisons

1. In **Single distribution**, reduce **Observations**, then play again.
   Compare the uncertainty about the mean direction with the spread of
   the observed headings. Increase **Concentration** and look for
   a tighter cluster in the simulated observations.
2. In **Continuous predictor**, enable **Concentration varies with x** and
   set a positive concentration slope. Follow the two fitted curves. Move
   **Predictor x** from −1 to +1 to see where the peak moves and whether
   it sharpens. Changing this viewing position reuses the same posterior
   draws.
3. In **Treatments and individuals**, choose an individual and compare
   its treatment estimates with the bold group baseline. Both mean
   direction and concentration have individual offsets whose spreads
   are learned across animals. The **Bimodal mixture** and **von
   Mises–uniform** examples let you investigate a different source of
   spread, in which observations can come from more than one component.

**Restart chains** replays the same simulated experiment. **New data**
generates another sample under the selected settings. For the derivation
of the likelihood and guidance on fitting these models, see the
[Bayesian primer](bayes-basics.md), [circular GLMs](math/glms.md) and
[circular GLMMs](math/glmms.md).
