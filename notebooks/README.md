# Notebook guide

Start with [circular_regression.ipynb](circular_regression.ipynb). It uses seeded
synthetic data, explains the model and diagnostics, and runs end to end in CI.
The other notebooks are research analyses or exploratory illustrations. Their
saved outputs do not establish that they run from a fresh kernel today.

## Inventory

| Notebook | Purpose | Data and current setup |
| --- | --- | --- |
| [Circular regression](circular_regression.ipynb) | Introductory Bambi regression for mean direction and concentration, with posterior predictive checks. | Generates its own data. Executed in CI. |
| [Beetles](beetles.ipynb) | Directional responses under lighting treatments; alternative mean/concentration models. | Uses bundled `datasets/SecurityLightBergsig.txt` (tab-separated). Its relative path assumes a kernel working directory of `notebooks/`. Contains package uninstall/install cells. |
| [Zebra escape](zebra_escape.ipynb) | Circular regression of zebra escape bearings using block-level data. | Bundled file: `datasets/zebra_bearings_blocks.parquet`. The notebook currently reads `zebra_bearings_blocks.parquet` without the `datasets/` prefix. Contains package uninstall/install cells. |
| [Locust direction mixture](locust_direction_mixture_model.ipynb) | Direct PyMC mixture and single-component models of locust directions, with model comparisons and figures. | Bundled file: `datasets/bimodal_locust.parquet`. The notebook currently reads `bimodal_locust.parquet` without the `datasets/` prefix. |
| [Locust density and order](locust_density_order.ipynb) | Trajectory processing and circular regression against density and collective order, including model sweeps. | Its preprocessing section expects an `XY_parquet_parts/` directory of full trajectory files. These are intentionally not bundled: the full locust datasets are tens of GB, and only rediscretized vectors are shared after preprocessing. Later cells read `mu_kappa_matched_model_sweep_pointwise_delta.parquet`, produced by an earlier model-sweep cell. |
| [Periodic, not circular](periodic_not_circular.ipynb) | Simulated Beta regression with periodic time predictors and a non-angular response. | Generates its own data; not executed in CI. |
| [Von Mises policy](von_mises_policy.ipynb) | Illustrates goal-directed trajectories sampled from von Mises heading distributions. | Generates its own data and writes figures; not executed in CI. |

## Environment and working directory

From the repository root, install the notebook environment and launch the tested
introduction:

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

Select the kernel belonging to that environment. The introductory notebook's
setup is also described in the [installation guide](../docs/getting-started.md).

The research notebooks import additional packages. The following is an inventory
of direct imports beyond the core modelling stack, not a validated environment
specification for reproducing every analysis:

| Notebook | Additional imports or file-reader requirements |
| --- | --- |
| Beetles | `seaborn` |
| Zebra escape | `seaborn`, a pandas Parquet engine such as `pyarrow` |
| Locust direction mixture | `seaborn`, `joblib`, `numba`, `tqdm`, a Parquet engine such as `pyarrow` |
| Locust density and order | `pyarrow`, `joblib`, `tqdm` |

Relative paths are resolved from the kernel's working directory. Check it with
`%pwd`; for the bundled research data, `notebooks/` is the intended starting
point. The two missing `datasets/` prefixes listed above still need correction.
Some figure cells write to `../figures/`; check that their output directories
exist. Several illustrations request Arial, so typography can vary by system.

The beetle and zebra notebooks still contain cells that uninstall or reinstall
Unwrap. Review those cells before running them in an existing environment.

## Locust data scope

The full locust trajectory datasets are tens of GB and are intentionally kept
outside this example repository. The shared data are rediscretized vectors
produced by preprocessing, rather than the full trajectories. The bundled
`datasets/bimodal_locust.parquet` is the input used by the direction-mixture
notebook.

The density-and-order notebook still includes preprocessing from the external
`XY_parquet_parts/` files. A fresh checkout therefore cannot run that section
without those inputs. Its cleanup should clearly separate trajectory
preprocessing from analysis of shared, preprocessed inputs and identify the
exact input required by each analysis; it should not require adding the full
trajectory archive to GitHub.

## Scope of the next cleanup

Work through one research notebook per PR: clarify setup, correct data paths,
and document inputs and outputs before changing execution order or removing
cells. Preserve model formulas, priors, mixture constraints, and scientific
results during that organizational work. Changes to the analyses themselves
need separate review and validation.
