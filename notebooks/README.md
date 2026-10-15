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
| [Zebra body orientation](zebra_escape.ipynb) | How do group orientation, alignment between animals, and consistency within animals change after an approach? | Bundled observations; setup cleanup under review; full rerun and prediction validation pending |
| [Locust direction mixtures](locust_direction_mixture_model.ipynb) | Do locust headings support forward alignment or two preferred directions in virtual reality? | Bundled heading steps; data path and setup still need cleanup |
| [Locust density and order](locust_density_order.ipynb) | How do the density and alignment of virtual neighbours relate to movement direction? | Original `XY_parquet_parts` trajectories are not bundled; the full workflow cannot run from this checkout alone |
| [Dung beetle orientation](beetles.ipynb) | How do manipulated orientation cues affect heading? | Bundled `datasets/SecurityLightBergsig.txt`; full rerun pending |
| [Periodic and circular functions](periodic_not_circular.ipynb) | Illustrate the distinction between periodic functions and circular responses | Exploratory illustration; full rerun pending |
| [Von Mises policy](von_mises_policy.ipynb) | Explore a von Mises directional policy | Exploratory illustration; full rerun pending |

## Setup

Follow the [Python setup guide](../docs/getting-started.md). For zebra, also
install `seaborn` and `pyarrow` in the same environment:

```bash
python -m pip install seaborn pyarrow
```

A separate zebra cleanup PR documents its setup and interpretation. Until
that is merged, its input path needs the `datasets/` prefix when running from
`notebooks/`. Research notebooks still contain setup and path assumptions from
the original analyses, including package installation cells. Read those cells
before running a notebook in an existing analysis environment.

The locust density notebook includes trajectory preprocessing and a model
sweep; its comparison table is generated during that workflow. The bundled
locust mixture data contain preprocessed heading steps, not the full raw
trajectory dataset.

## Interpretation

A change in **mean direction** describes reorientation. A change in
**concentration** describes greater or lower directional consistency. In a
hierarchical model, consistency within animals and alignment between animals
are distinct quantities. Read the notebook's figure descriptions to see which
level is being summarised.

For the statistical framework and biological context, see the
[model guide](../docs/models.md) and
[mathematical guide](../docs/math/von-mises.md).

