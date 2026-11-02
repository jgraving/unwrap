# Circular models: examples and biological analyses

Start with **[circular regression](circular_regression.ipynb)**. It generates
its own data and walks through estimating mean direction and directional
concentration, displaying uncertainty, and checking model predictions.

The research notebooks apply these ideas to animal behaviour. Each has
been rerun end to end from this checkout on the current stack, and
their saved outputs come from that validated run; each row above says
what that run validated.

| Notebook | Biological question or purpose | Inputs and readiness |
| --- | --- | --- |
| [Circular regression](circular_regression.ipynb) | How does a predictor affect mean direction and directional consistency? | Simulated data; executed automatically in CI |
| [Grevy's zebra body orientation during an escape response](zebra_escape.ipynb) | How do group orientation, alignment between animals, and consistency within animals change after an approach? | Block bearings pulled from `unwrap-data` (`load_zebra_bearings`); fully rerun and validated from this checkout (MCMC diagnostics pass; kappa random effects fitted non-centered; bootstrap-vs-posterior PPC, with the per-observation parameters from `compute_deterministics` cross-checked against the stored log-likelihood, and the `compute_deterministics` workaround for the bambi `predict(kind="response_params")` kappa bug, bambinos/bambi#1084) |
| [Decomposing goal-directedness and directional precision in sensory-deprived echolocating bats during homing](goldshtein_bats.ipynb) | Do kilometer-scale bat flights align with the destination (goal-homing), and how strongly? | Tracks, metadata and the Agamon–Hula lake outline pulled from `unwrap-data` (`load_bats_tracks`, `load_bats_info`, `load_bats_lake_outline`); fully rerun and validated from this checkout (fixed 60 m re-discretization with gap segmentation; roost-vs-lake destination assignment from the harmonic lake field and the radial-to-roost bearing; vMU von Mises–uniform mixture fit vs the vM baseline on the goal-relative bearings; PSIS-LOO total / by treatment / per animal, best model the mixture; population and tag-level mixture summaries: mixture weight, goal-directed concentration, directed-component density at the goal and the ±45° mass around it; hierarchical-bootstrap PPC by treatment plus the per-animal bootstrap-vs-posterior grid; per-step goal-component responsibility with kernel-smoothed approach curves; the vM per-observation log-likelihood taken from the PyMC graph, working around the bambi vonmises kappa link bug, bambinos/bambi#1084) |
| [Order and density shape locust directional alignment in virtual reality](locust_density_order.ipynb) | How do the density and alignment of virtual neighbours relate to movement direction? | Raw tracking parts pulled from `unwrap-data` (`load_density_order`); fully rerun and validated from this checkout (MCMC diagnostics pass; individual-level PSIS-LOO sweep selecting the Density × Order specification; final random-intercept fit validated with alignment and predictive-density statistics plus cell- and individual-level predictive checks; per-observation log-likelihood attached from the PyMC graph via `compute_deterministics`, working around the bambi `predict(kind="response_params")` kappa bug, bambinos/bambi#1084) |
| [Testing optic-flow vs. vectorial encoding in marching locusts](locust_direction_mixture_model.ipynb) | Do locust headings support forward alignment or two preferred directions in virtual reality? | Raw trajectories pulled from `unwrap-data` (`load_directional` + `rediscretize_directional`); fully rerun and validated from this checkout (MCMC diagnostics pass; animal ids renumbered densely for the sparse file-derived ids; PSIS-LOO comparison of the bimodal mixture against the unimodal baseline, population and individual predictive checks, and the optic-flow vs vectorial-encoding hypothesis tests) |
| [Artificial light reorients nocturnal dung beetles](beetles.ipynb) | How do manipulated orientation cues affect heading? | Security-light trials pulled from `unwrap-data` (`load_beetle_security_light`); fully rerun and validated from this checkout (seeded 2000+2000 von Mises Bambi fit of light ON/OFF heading with non-centered per-individual concentration random effects, zero divergences; circular summaries of the mean and concentration parameters; per-individual mean directions with the between-individual concentration kappa_mu of the individual means (intercept / light / difference) and the per-draw summary; the PPC figure with the per-observation log-likelihood cross-checked against the PyMC graph, working around the bambi vonmises kappa link bug, bambinos/bambi#1084) |
| [Periodic but not circular](periodic_not_circular.ipynb) | Illustrate the distinction between periodic functions and circular responses | Simulated 24-hour activity data (seeded); fully rerun and validated from this checkout (seeded 1000+1000 NUTS fit of the beta-process activity rhythm with harmonic time-of-day predictors on the mean (logit link) and the concentration (softplus link), zero divergences; the figure overlays the posterior-predictive median and 50%/90% HDI on the simulated rhythm) |
| [Von Mises navigation policy](von_mises_policy.ipynb) | Explore a von Mises directional policy | Simulated trajectories (seeded); fully rerun and validated from this checkout (exact von Mises densities as radial modulations and seeded Monte-Carlo steering-policy trajectories, one panel per concentration 0 / 0.5 / 2 / 5; a single two-resolution figure) |
## Setup

Follow the [Python setup guide](../docs/getting-started.md). The animal-behaviour
notebooks also need the optional notebook stack in the same environment:

```bash
python -m pip install seaborn pyarrow
```

The animal-behaviour notebooks load their observations through the
`unwrap.datasets` accessors, which fetch raw data from the
[`jgraving/unwrap-data`](https://github.com/jgraving/unwrap-data) repository on
first use (cloned to `~/.cache/unwrap-data`) or from a local checkout passed as
`data_root=`. After the initial clone, the loaders resolve to cached files with
no network access.

The research notebooks still contain setup and path assumptions from the
original analyses (including package version checks). Read those cells before
running a notebook in an existing analysis environment.

The locust direction notebook runs its own trajectory rediscretization
(`load_directional` then `rediscretize_directional`); the density notebook has
that preprocessing plus a model sweep, and its comparison table is generated
during that workflow. The zebra notebook reads a preprocessed table that is
published in `unwrap-data` alongside its raw source.

## Interpretation

A change in **mean direction** describes reorientation. A change in
**concentration** describes greater or lower directional consistency. In a
hierarchical model, consistency within animals and alignment between animals
are distinct quantities. Read the notebook's figure descriptions to see which
level is being summarised.

For the statistical framework and biological context, see the
[model guide](../docs/models.md) and
[mathematical guide](../docs/math/von-mises.md).

