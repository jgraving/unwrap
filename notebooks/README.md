# Notebook guide

Start with [Circular regression](circular_regression.ipynb), a worked example
using simulated observations. It introduces model fitting, summaries of mean
direction and concentration, and comparisons between predictions and data.

## Animal navigation and collective movement

These examples use circular regression to connect directional behaviour to
biological questions. The descriptions follow the case studies in the
[accompanying manuscript](https://github.com/jgraving/circstats/blob/main/elife-template.tex).

| Notebook | Question |
| --- | --- |
| [Artificial light and dung beetle orientation](beetles.ipynb) | Examine how a floodlight changes mean direction, within-individual directional consistency, and alignment between individuals. |
| [Grevy’s zebra escape response](zebra_escape.ipynb) | Compare body orientation before and after a human approach, separating group reorientation, between-individual alignment, and within-individual directional consistency. |
| [Order, density, and locust alignment](locust_density_order.ipynb) | Ask how the directional coherence and density of a virtual swarm influence a focal locust’s alignment with its neighbours. |
| [Directional decisions in marching locusts](locust_direction_mixture_model.ipynb) | Compare optic-flow and vectorial-encoding hypotheses through their predictions of forward marching or lateral turns toward virtual neighbour groups. |

## Further illustrations

- [Periodic but not circular](periodic_not_circular.ipynb): Model a response that varies over a daily cycle but is not itself an angle.
- [Goal-directed movement](von_mises_policy.ipynb): Explore simulated trajectories in which headings are drawn from a von Mises distribution.

## Run a notebook

From the repository root:

```bash
python -m venv .venv
# macOS / Linux
source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -e ".[notebook]"
python -m jupyterlab notebooks/circular_regression.ipynb
```

Alternatively, with `uv` (optional):

```bash
uv venv --python 3.12
uv pip install -e ".[notebook]"
uv run --no-project jupyter lab notebooks/circular_regression.ipynb
```

Select the kernel for your environment and run the cells in order. The research
examples fit models to experimental observations and can take longer to run
than the introductory tutorial. Later figures depend on the model fits from
earlier cells; opening a notebook with saved figures does not restore a fitted
model.

## Data and additional packages

The beetle, zebra, and locust direction examples include their analysis data.
The introductory tutorial and further illustrations generate simulated data.

| Example | Included data | Additional packages |
| --- | --- | --- |
| Dung beetle orientation | `datasets/SecurityLightBergsig.txt` | `seaborn` |
| Grevy’s zebra escape response | `datasets/zebra_bearings_blocks.parquet` | `seaborn`, `pyarrow` |
| Locust directional decisions | `datasets/bimodal_locust.parquet` | `seaborn`, `joblib`, `numba`, `tqdm`, `pyarrow` |
| Locust order and density | The preprocessing section requires external trajectories. | `pyarrow`, `joblib`, `tqdm` |

Install the additional packages needed for your chosen example in the same
environment, for example `python -m pip install seaborn pyarrow` for the zebra
analysis.

Start the research notebooks in `notebooks/` and follow any setup instructions
at the top of your chosen notebook. If a data-loading cell uses only a filename,
add the `datasets/` prefix to read the included file. Skip package-installation
cells if you have already installed Unwrap in your environment, and create any
figure output directories requested by the notebook.

The full locust trajectories are tens of GB and are intentionally not distributed
here. Shared locust data consist of rediscretized movement vectors produced
after preprocessing. The order-and-density notebook includes a preprocessing
section that requires the external `XY_parquet_parts/` trajectory directory;
that section cannot be run from the included data alone. Its later comparisons
also use a table produced by an earlier model-fitting step.
