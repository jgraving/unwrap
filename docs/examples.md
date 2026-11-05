# Examples

**Every example is available in Python and in R.** The Python notebooks
live under `notebooks/`, built on **PyMC** (with **Bambi** for the formula
interface), and the site's "Open in Colab" button launches them in the
cloud; the R notebooks live under `notebooks_r/` next to their Python
counterparts, fitted with [brms](https://paulbuerkner.com/brms/) on
CmdStan. The two sides of each pair use the same observations, formulas,
and prior scales. The [Python setup page](getting-started.md) covers the
Python environment and the [R setup page](getting-started-r.md) the R one;
[R and MATLAB](languages.md) covers the standalone R scripts and the MATLAB
examples.

The [notebook guide](https://github.com/jgraving/unwrap/blob/main/notebooks/README.md)
lists all examples, their inputs, and what each validated run covers.

## The tutorial

<p class="unwrap-lang-row">
  <a class="unwrap-lang unwrap-lang-py" href="../notebooks/circular_regression/">Python</a>
  <a class="unwrap-lang unwrap-lang-r" href="../notebooks_r/circular_regression/">R</a>
</p>

The circular regression tutorial
using simulated observations: a math example on toy data, not a
biological study, and the cheapest way to see the whole workflow. It
fits a circular generalized linear model (GLM) with a formula for both
the mean direction and the concentration. The tutorial creates its own
data, so no download is needed.

## Biological examples

The biological examples extend the same framework to real experiments.
Repeated observations can be represented with generalized linear mixed
models (GLMMs), which estimate animal-specific variation alongside
population-level effects. The [model guide](models.md) connects these
terms to biological questions and model formulas.

### Grevy's zebra body orientation during an escape response

<p class="unwrap-lang-row">
  <a class="unwrap-lang unwrap-lang-py" href="../notebooks/zebra_escape/">Python</a>
  <a class="unwrap-lang unwrap-lang-r" href="../notebooks_r/zebra_escape/">R</a>
</p>

A focal group of Grevy's zebra (*Equus grevyi*) was recorded with an overhead
drone and automated tracking at the Mpala Conservancy in Laikipia, Kenya: the
body orientations come from the drone recordings of B. Koger, A. Deshpande,
J. T. Kerby, J. M. Graving, B. R. Costelloe, and I. D. Couzin, "Quantifying
the movement, behaviour and environmental context of group-living animals
using drones and computer vision," *Journal of Animal Ecology* **92**(6),
1415-1433 (2023); [doi:10.1111/1365-2656.13904](https://doi.org/10.1111/1365-2656.13904).
The first flights capture baseline (grazing) body orientation, and the final
flight contains a disturbance event in which a human observer approaches the
group. The zebra analysis
body orientation in the **pre-stimulus baseline period** (pre-approach) with
the **post-stimulus period** (post-approach), using circularly averaged
orientation over 10 s windows. The hierarchical model separates three
complementary aspects of the escape response: a group-level shift in the
mean direction, alignment between individuals (how the individual mean
directions cluster), and within-individual consistency (how tightly each
zebra's observations cluster around its own mean direction).

### Decomposing goal-directedness and directional precision in sensory-deprived echolocating bats during homing

<p class="unwrap-lang-row">
  <a class="unwrap-lang unwrap-lang-py" href="../notebooks/goldshtein_bats/">Python</a>
  <a class="unwrap-lang unwrap-lang-r" href="../notebooks_r/goldshtein_bats/">R</a>
</p>

The homing bat analysis
kilometer-scale homing flights of echolocating bats from the experiments of
A. Goldshtein et al., "Acoustic cognitive map-based navigation in
echolocating bats," *Science* **386**(6721), 561-567 (2024);
[doi:10.1126/science.adn6269](https://doi.org/10.1126/science.adn6269):
free-flying bats translocated within their home range are tracked on their
return to a known goal, the roost or the Agamon–Hula lake region, under a
control condition or one of several sensory-deprivation treatments (temporary
occlusion of vision, experimental disruption of magnetoreception and
olfaction). The original study found that bats typically
do home successfully even under sensory deprivation; this reanalysis asks
whether that deprivation nonetheless produces detectable changes in
goal-relative orientation. Because the treatments could change how straight a
bat flies, the reanalysis avoids feature-based selection on post-treatment
quantities (for example a "homing" threshold) and instead fits a von
Mises–uniform mixture to the full set of goal-relative bearings, in which
each step arises either from a directed von Mises component (goal-directed
"exploitation") or an undirected uniform component (search "exploration"),
compared with a unimodal von Mises baseline.

### Order and density shape locust directional alignment in virtual reality

<p class="unwrap-lang-row">
  <a class="unwrap-lang unwrap-lang-py" href="../notebooks/locust_density_order/">Python</a>
  <a class="unwrap-lang unwrap-lang-r" href="../notebooks_r/locust_density_order/">R</a>
</p>

The order and density analysis
a virtual-reality experiment from S. Sayin et al., "The behavioral
mechanisms governing collective motion in swarming locusts," *Science*
**387**(6737), 995-1000 (2025);
[doi:10.1126/science.adq7832](https://doi.org/10.1126/science.adq7832),
in which marching locusts walk on a motion-compensating ball surrounded by a
panoramic projection of moving virtual conspecifics, parameterized by *order*
(how aligned the virtual locusts' motion is, 0 to 1) and *density* (virtual
locusts per m$^2$, 1–64). The question is whether the focal locust's directional alignment is
driven primarily by the directional structure of the neighbor motion (order)
or by the number of neighbors (density); a PSIS-LOO model-selection sweep
over the order and density specifications supports an order-dominant
interpretation - order captures the bulk of the improvement in out-of-sample
performance, and adding density provides at most a secondary gain.

### Testing optic-flow vs. vectorial encoding in marching locusts

<p class="unwrap-lang-row">
  <a class="unwrap-lang unwrap-lang-py" href="../notebooks/locust_direction_mixture_model/">Python</a>
  <a class="unwrap-lang unwrap-lang-r" href="../notebooks_r/locust_direction_mixture_model/">R</a>
</p>

The mixture analysis
a second paradigm from the same experiments, in which moving virtual
conspecifics are placed symmetrically at a lateral distance of 8 or 50 cm. A hierarchical
bimodal von Mises mixture fitted to the full distribution of headings
discriminates the two mechanistic accounts: optic flow (headings
concentrated near the forward march direction, consistent with optomotor
alignment to wide-field motion) versus vectorial encoding (a choice-like
commitment to one of the two lateral swarms, with leftward and rightward
turns).

### Artificial light reorients nocturnal dung beetles

<p class="unwrap-lang-row">
  <a class="unwrap-lang unwrap-lang-py" href="../notebooks/beetles/">Python</a>
  <a class="unwrap-lang unwrap-lang-r" href="../notebooks_r/beetles/">R</a>
</p>

The beetle analysis a two-condition
arena experiment in the nocturnal ball-rolling dung beetle *Escarabaeus
satyrus* from the experiments of J. J. Foster et al., "Light pollution
forces a change in dung beetle orientation behavior," *Current Biology*
**31**(17), 3935-3942 (2021);
[doi:10.1016/j.cub.2021.06.038](https://doi.org/10.1016/j.cub.2021.06.038):
each beetle is placed at the center of a circular arena and allowed to roll
to the boundary, with its exit bearing recorded repeatedly under natural
nighttime illumination (control) and with an artificial floodlight switched
on (treatment). The model fits the
treatment effect on both the mean direction and the concentration, and the
analysis quantifies, in a single coherent fit, the group-level reorientation
toward the floodlight, the change in within-individual concentration, the
post-hoc change in the between-individual clustering of the individuals' mean
directions ($\kappa_\mu$), and "floodlight attraction" as the posterior
predictive mass within $\pm45^{\circ}$ of the floodlight direction.

## Method examples

Two further notebooks round out the framework: one shows when a circular
model is *not* the right tool, and one works the von Mises distribution
backwards from likelihood to behavior.

### Periodic but not circular: a 24-hour activity rhythm

<p class="unwrap-lang-row">
  <a class="unwrap-lang unwrap-lang-py" href="../notebooks/periodic_not_circular/">Python</a>
  <a class="unwrap-lang unwrap-lang-r" href="../notebooks_r/periodic_not_circular/">R</a>
</p>

The periodicity example
*periodic* and *circular*: time of day repeats every 24 hours, but the
quantity of interest - activity - is a level bounded in (0, 1), not an angle.
The periodicity belongs to the *predictor*, so the model is a Beta GLM whose
mean (logit-linked) and precision (softplus-linked) each vary sinusoidally
across the day, rather than a von Mises likelihood on a circle: the wrap at
midnight is a periodic continuity, not the angular equivalence that makes
bearings circular. The notebook re-simulates a seeded 24-hour rhythm, fits
the model, and recovers the injected daily mean and precision curves with
posterior-predictive uncertainty bands.

### Von Mises navigation policy

<p class="unwrap-lang-row">
  <a class="unwrap-lang unwrap-lang-py" href="../notebooks/von_mises_policy/">Python</a>
  <a class="unwrap-lang unwrap-lang-r" href="../notebooks_r/von_mises_policy/">R</a>
</p>

The policy example
the opposite direction - as a *navigation policy* rather than a likelihood.
With utility defined as the negated distance to a goal, the heading
distribution that maximizes utility in constant-length steps is exactly a von
Mises distribution centered on the goal direction, with the concentration
setting how strongly each step concentrates toward it. The notebook evaluates
the policy density at each concentration and samples trajectories under it:
at zero concentration the policy is uniform and the walk is a
constant-step random walk with no net progress, while rising concentration
progressively aligns the trajectories with the goal. It is the visual
companion to the likelihood use in the [beetle](notebooks/beetles.ipynb) and
[zebra](notebooks/zebra_escape.ipynb) analyses, where fitted mean direction
and concentration describe observed directional control.

## Before running one of these notebooks

Each notebook on this page was rerun end to end from this checkout on the
current stack, and the outputs shown come from that validated run. The
[Python setup guide](getting-started.md) covers the environment; the R
notebooks follow the [R setup guide](getting-started-r.md) and load the same
published data through a Git-managed local checkout. The
animal-behavior notebooks load their observations through the
`unwrap.datasets` accessors, which fetch the raw data from the
[`jgraving/unwrap-data`](https://github.com/jgraving/unwrap-data) repository
on first use (cloned to `~/.cache/unwrap-data`) and then resolve to the
cache without network access.
