# Circular models: examples and biological analyses

Start with **[circular regression](circular_regression.ipynb)**. It generates
its own data and walks through estimating mean direction and directional
concentration, displaying uncertainty, and checking model predictions.

The research notebooks apply these ideas to animal behaviour. They retain
saved results from the research workflow; those outputs are not evidence that
a fresh run succeeds with today's dependencies.

| Notebook | Biological question or purpose | Inputs and readiness |
| --- | --- | --- |
| [Circular regression](circular_regression.ipynb) | How does a predictor affect mean direction and directional consistency? | Simulated data; executed automatically in CI |
| [Zebra body orientation](zebra_escape.ipynb) | How do group orientation, alignment between animals, and consistency within animals change after an approach? | Block bearings pulled from `unwrap-data` (`load_zebra_bearings`); fully rerun and validated (MCMC diagnostics pass; bootstrap-vs-posterior PPC, with the `predict(kind="response_params")` kappa fix in place) |
| [Homing bats](goldshtein_bats.ipynb) | Do kilometer-scale bat flights align with the destination (goal-homing), and how strongly? | Tracks, metadata and the Agamon–Hula lake outline pulled from `unwrap-data` (`load_bats_tracks`, `load_bats_info`, `load_bats_lake_outline`); full workflow runs from this checkout |
| [Locust direction mixtures](locust_direction_mixture_model.ipynb) | Do locust headings support forward alignment or two preferred directions in virtual reality? | Raw trajectories from `unwrap-data` (`load_directional`) + in-notebook rediscretization; full rerun pending |
| [Locust density and order](locust_density_order.ipynb) | How do the density and alignment of virtual neighbours relate to movement direction? | Tracking parts from `unwrap-data` (`load_density_order`); full workflow rerun pending |
| [Dung beetle orientation](beetles.ipynb) | How do manipulated orientation cues affect heading? | Security-light trials pulled from `unwrap-data` (`load_beetle_security_light`); full rerun pending |
| [Periodic and circular functions](periodic_not_circular.ipynb) | Illustrate the distinction between periodic functions and circular responses | Simulated data; full workflow runs from this checkout |
| [Von Mises policy](von_mises_policy.ipynb) | Explore a von Mises directional policy | Simulated trajectories; full workflow runs from this checkout |

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
during that workflow. The zebra and beetle notebooks read preprocessed tables
that are published in `unwrap-data` alongside their raw sources.

## Interpretation

A change in **mean direction** describes reorientation. A change in
**concentration** describes greater or lower directional consistency. In a
hierarchical model, consistency within animals and alignment between animals
are distinct quantities. Read the notebook's figure descriptions to see which
level is being summarised.

For the statistical framework and biological context, see the
[model guide](../docs/models.md) and
[mathematical guide](../docs/math/von-mises.md).

