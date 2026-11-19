"""Load verified bat tracks and posterior summaries for the revised layout."""
from functools import partial
from unittest.mock import patch

import geopandas as gpd
import numpy as np
from pyproj import Transformer
from scipy.special import expit
from scipy.stats import vonmises

import goldshtein_bats_publication_figure as publication
from unwrap import datasets
from unwrap.link.numpy import circular_modulo, softplus


def _read_context(data_root):
    # The notebook helpers resolve a lake path at import time. Supply the
    # pinned path to those imports so a figure build does not need a checkout
    # in the default dataset cache.
    configured_lake_path = partial(datasets.bats_lake_outline_path, data_root=data_root)
    with patch.object(datasets, "bats_lake_outline_path", configured_lake_path):
        import goldshtein_bats_posthoc as posthoc
        import goldshtein_bats_preprocess as preprocess
        import goldshtein_bats_plots as plots

    lake_path = configured_lake_path()
    with (patch.object(publication, "_import_notebook_modules",
                       return_value=(preprocess, posthoc, plots)),
          patch.object(plots, "LAKE_GEOJSON", str(lake_path))):
        ctx = publication.build_context(data_root)
    return ctx, posthoc, plots, lake_path


def _population_parameters(posterior, levels, posthoc):
    """Return the pinned vMU draws with individual offsets set to zero."""
    if levels != [1, 2, 3, 4]:
        raise ValueError(f"Unexpected treatment order for the bat fit: {levels}")

    def population_draws(intercept, effect):
        base = posthoc._stack_samples_first(posterior[intercept])
        offsets = posthoc._stack_samples_first(posterior[effect])
        if base.shape != (8000,) or offsets.shape != (8000, len(levels) - 1):
            raise ValueError(f"Unexpected shape for {intercept} or {effect} in the pinned fit")
        return base[:, None] + np.column_stack([np.zeros(base.size), offsets])

    mu = circular_modulo(population_draws("mu0", "mu_treat"))
    kappa = softplus(population_draws("eta_kappa0", "kappa_treat"))
    weight = expit(population_draws("logit_p0", "beta_treat"))
    return mu, kappa, weight


def _goal_window_mass(tolerance, mu, kappa, weight, kappa_uniform):
    """Integrate the predictive mixture from -tolerance to +tolerance."""
    directed_mass = (vonmises.cdf(tolerance, kappa, loc=mu)
                     - vonmises.cdf(-tolerance, kappa, loc=mu))
    uniform_mass = vonmises.cdf(tolerance, kappa_uniform) - vonmises.cdf(-tolerance, kappa_uniform)
    return np.clip(weight * directed_mass + (1 - weight) * uniform_mass, 0, 1)


def _mixture_predictions(posterior, levels, counts, posthoc):
    """Summarize the pinned population-level vMU predictive densities."""
    mu, kappa, weight = _population_parameters(posterior, levels, posthoc)

    def summary(draws):
        return (float(np.mean(draws)), *map(float, posthoc.hdi_1d(draws, prob=0.95)))

    # Use the near-uniform background specified by the notebook's predictive
    # summaries for both the density and the window probabilities.
    kappa_uniform = posthoc.KAPPA_UNIFORM
    bearings = np.deg2rad(np.arange(-180, 181))
    goal_index = bearings.size // 2
    uniform_density = posthoc.vonmises_pdf_scalar(bearings, 0.0, kappa_uniform)
    density_at_goal = np.empty_like(weight)
    directed_at_goal = np.empty_like(weight)
    density_profiles = {}
    for i, level in enumerate(levels):
        directed_density = posthoc.vonmises_pdf_scalar(
            bearings, mu[:, i, None], kappa[:, i, None]
        )
        p = weight[:, i, None]
        density = p * directed_density + (1 - p) * uniform_density
        density_at_goal[:, i] = density[:, goal_index]
        directed_at_goal[:, i] = directed_density[:, goal_index]
        lo, hi = posthoc.hdi_along_samples(density, prob=0.95)
        density_profiles[level] = {"mean": density.mean(axis=0), "lo": lo, "hi": hi}

    half_width = posthoc.MASS_HALF_WIDTH_DEG
    mass = _goal_window_mass(np.deg2rad(half_width), mu, kappa, weight, kappa_uniform)
    by_treatment = mass - mass[:, :1]
    pooled = np.average(by_treatment[:, 1:], axis=1,
                        weights=[counts[level] for level in levels[1:]])
    goal_window_groups = {level: summary(draws) for level, draws in zip(levels, mass.T, strict=True)}
    goal_window_contrasts = {level: summary(by_treatment[:, i]) for i, level in enumerate(levels[1:], 1)}
    return {
        "density_bearings": bearings,
        "density_profiles": density_profiles,
        "density_at_zero_by_treatment": {
            level: summary(draws) for level, draws in zip(levels, density_at_goal.T, strict=True)
        },
        "full_contrasts": {
            level: summary(density_at_goal[:, i] - density_at_goal[:, 0])
            for i, level in enumerate(levels[1:], 1)
        },
        "directed_contrasts": {
            level: summary(directed_at_goal[:, i] - directed_at_goal[:, 0])
            for i, level in enumerate(levels[1:], 1)
        },
        "goal_window_half_width_degrees": half_width,
        "goal_window_groups": goal_window_groups,
        "goal_window_contrasts": goal_window_contrasts,
        "pooled_goal_window_contrast": summary(pooled),
    }


def load_layout_data(data_root):
    """Validate the input pins and assemble the map, model, and heading panels."""
    root = publication.prepare_inputs(data_root)
    return assemble_layout_data(*_read_context(root))


def assemble_layout_data(ctx, posthoc, plots, lake_path):
    """Build the composite's panel data from the restored notebook context."""
    observed = ctx["df_track"]
    tracks = ctx["track_w_treat"]
    if not (observed.index.is_unique and tracks.index.is_unique and observed.index.isin(tracks.index).all()):
        raise ValueError("Modeling rows do not have unique source track identifiers")
    coordinates = ["x_m", "y_m"]
    source_rows = tracks.loc[observed.index]
    if (not np.array_equal(observed["tag"].to_numpy(), source_rows["tag"].to_numpy()) or
            not np.array_equal(observed[coordinates].to_numpy(), source_rows[coordinates].to_numpy())):
        raise ValueError("Modeling rows disagree with source track tags or coordinates")

    levels = [int(code) for code in ctx["treatment_order"]]
    counts = {int(code): int(count) for code, count in
              observed.groupby("treatment", observed=True)["tag"].nunique().items()}
    if counts != {1: 29, 2: 14, 3: 15, 4: 9}:
        raise ValueError(f"Unexpected modeled bat counts: {counts}")

    model_keys = {"von Mises": "vm", "von Mises-Uniform": "vmu"}
    scores = {model_keys[row.model]: {"elpd": float(row.elpd), "se": float(row.se_elpd)}
              for row in ctx["res_total"].itertuples()}
    lake = gpd.read_file(lake_path).to_crs(publication.METRIC_EPSG).union_all()
    transform = Transformer.from_crs(4326, publication.METRIC_EPSG, always_xy=True)
    sites = [transform.transform(lon, lat) for lat, lon in plots.TRANS_POINTS]
    roost = transform.transform(publication.ROOST_LON, publication.ROOST_LAT)
    return {
        "tracks": tracks[["tag", "treatment", "dest", *coordinates, "gap_before",
                          "step_bearing_deg", "bearing_to_dest_deg"]],
        "evidence": observed["evi_post_mean"],
        "levels": levels,
        "curves": ctx["out_steps_k"],
        "lake": lake,
        "sites": sites,
        "roost": roost,
        "counts": counts,
        "scores": scores,
        **_mixture_predictions(ctx["idata_vmu"].posterior, levels, counts, posthoc),
    }
