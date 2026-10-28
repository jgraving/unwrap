"""Reproduce circular diagnostic comparisons and two small PyMC fits.

Run from the repository root::

    python examples/circular_rhat_validation.py --output docs/validation

The 1.01 exceedance rates are descriptive comparisons, not a calibration of the
new circular statistic or of the exploratory maximum/mean of two R-hats.
"""

import argparse
import importlib.metadata
import json
from pathlib import Path

import arviz as az
import numpy as np
import pymc as pm

from unwrap.diagnostics import rhat, rhat_components, summary

SCENARIOS = (
    "iid_same_target",
    "whole_turn_offsets",
    "near_branch_cut",
    "different_mean_directions",
    "different_concentrations",
    "within_chain_shift",
    "stationary_ar1",
    "different_antipodal_pairs",
    "same_projection_marginals",
    "small_direction_difference",
)
METRICS = (
    "joint_circular",
    "rank_sin",
    "rank_cos",
    "rank_component_max",
    "rank_component_mean",
    "rank_wrapped_angle",
)


def wrap(values):
    return (values + np.pi) % (2 * np.pi) - np.pi


def antipodal_pairs(rng, draws, rotations):
    """Each half-chain has exactly balanced antipodal pairs, with nonzero jitter."""
    chains = []
    for rotation in rotations:
        halves = []
        for _ in range(2):
            jitter = rng.normal(0, 0.06, draws // 4)
            half = np.concatenate([jitter, jitter + np.pi]) + rotation
            rng.shuffle(half)
            halves.append(half)
        chains.append(np.concatenate(halves))
    return np.asarray(chains)


def simulate(scenario, draws, seed):
    rng = np.random.default_rng(seed)
    if scenario in {"iid_same_target", "whole_turn_offsets", "near_branch_cut"}:
        values = rng.vonmises(0.4, 4, (4, draws))
        if scenario == "whole_turn_offsets":
            values += 2 * np.pi * np.arange(4)[:, None]
        elif scenario == "near_branch_cut":
            values += np.pi - 0.4
        return values
    if scenario == "different_mean_directions":
        return rng.vonmises(np.array([0, 0, 1.5, 1.5])[:, None], 8, (4, draws))
    if scenario == "different_concentrations":
        return rng.vonmises(0.4, np.array([2, 2, 20, 20])[:, None], (4, draws))
    if scenario == "within_chain_shift":
        values = rng.vonmises(0.4, 8, (4, draws))
        values[:, draws // 2 :] += 1.5
        return values
    if scenario == "stationary_ar1":
        rho, sigma = 0.95, 0.4
        values = np.empty((4, draws))
        values[:, 0] = rng.normal(0, sigma, 4)
        noise = rng.normal(0, sigma * np.sqrt(1 - rho**2), (4, draws - 1))
        for draw in range(1, draws):
            values[:, draw] = rho * values[:, draw - 1] + noise[:, draw - 1]
        return wrap(values + 0.4)
    if scenario == "different_antipodal_pairs":
        return antipodal_pairs(rng, draws, (0, 0, np.pi / 2, np.pi / 2))
    if scenario == "same_projection_marginals":
        return antipodal_pairs(rng, draws, (np.pi / 4, np.pi / 4, 3 * np.pi / 4, 3 * np.pi / 4))
    if scenario == "small_direction_difference":
        return rng.normal(np.array([-0.08, -0.08, 0.08, 0.08])[:, None], 0.5, (4, draws))
    raise ValueError(scenario)


def diagnose(values):
    components = rhat_components(values)
    sine, cosine = float(components["sin"]), float(components["cos"])
    return {
        "joint_circular": float(rhat(values, method="circular")),
        "rank_sin": sine,
        "rank_cos": cosine,
        "rank_component_max": max(sine, cosine),
        "rank_component_mean": (sine + cosine) / 2,
        "rank_wrapped_angle": float(az.rhat(wrap(values))),
    }


def summarize(rows):
    summaries = []
    for scenario in SCENARIOS:
        for draws in sorted({row["draws_per_chain"] for row in rows}):
            selected = [
                row
                for row in rows
                if row["scenario"] == scenario and row["draws_per_chain"] == draws
            ]
            metrics = {}
            for metric in METRICS:
                values = np.array([row[metric] for row in selected])
                metrics[metric] = {
                    "q05": float(np.quantile(values, 0.05)),
                    "median": float(np.median(values)),
                    "q95": float(np.quantile(values, 0.95)),
                    "fraction_above_1_01": float(np.mean(values > 1.01)),
                }
            summaries.append({"scenario": scenario, "draws_per_chain": draws, "metrics": metrics})
    return summaries


def validate_mixed_api(idata):
    """Check a real mixed posterior against independently computed ArviZ outputs."""
    angular = ["Intercept", "mu"]
    mixed_rhat = rhat(idata, circ_var_names=angular)
    mixed_summary = summary(idata, circ_var_names=angular, round_to="none")
    ordinary_rhat = az.rhat(idata)
    ordinary_summary = az.summary(idata, round_to="none")
    arviz_circular = az.summary(idata, circ_var_names=angular, round_to="none")
    joint = rhat(idata, var_names=angular, method="circular")

    np.testing.assert_array_equal(mixed_rhat.kappa, ordinary_rhat.kappa)
    np.testing.assert_array_equal(
        mixed_summary.loc["kappa", ordinary_summary.columns], ordinary_summary.loc["kappa"]
    )
    for name in angular:
        np.testing.assert_array_equal(mixed_rhat[name], joint[name])
        np.testing.assert_array_equal(mixed_summary.loc[name, "r_hat"], joint[name])
        np.testing.assert_array_equal(
            arviz_circular.loc[name, ["r_hat", "ess_bulk", "ess_tail"]],
            ordinary_summary.loc[name, ["r_hat", "ess_bulk", "ess_tail"]],
        )
        assert mixed_summary.loc[name, ["ess_bulk", "ess_tail", "mcse_sd"]].isna().all()
    projection_columns = []
    for suffix, projection in (("sin", np.sin), ("cos", np.cos)):
        projected = az.summary(
            projection(idata.posterior[angular]), kind="diagnostics", round_to="none"
        )
        columns = [f"{metric}_{suffix}" for metric in projected.columns]
        projection_columns.extend(columns)
        np.testing.assert_array_equal(
            mixed_summary.loc[angular, columns].values, projected.loc[angular].values
        )
    assert mixed_summary.loc["kappa", projection_columns].isna().all()

    records = []
    for name in angular + ["kappa"]:
        row = {
            "variable": name,
            "is_circular": name in angular,
            "raw_arviz_rhat": float(ordinary_rhat[name]),
            "arviz_circ_var_names_rhat": float(arviz_circular.loc[name, "r_hat"]),
            "mixed_rhat": float(mixed_rhat[name]),
        }
        # JSON null denotes an intentionally unavailable circular scalar diagnostic.
        row["mixed_summary"] = {
            metric: None if np.isnan(value) else float(value)
            for metric, value in mixed_summary.loc[name].items()
        }
        records.append(row)
    return {
        "circ_var_names": angular,
        "ordinary_variables": ["kappa"],
        "ordinary_rhat_and_summary_match_arviz_exactly": True,
        "arviz_circ_var_names_keeps_raw_rhat_and_ess": True,
        "angular_projection_diagnostics_match_arviz_exactly": True,
        "angular_raw_ess_and_mcse_sd_are_null": True,
        "rows": records,
    }


def fit_intercept(scenario, seed, draws=500, tune=500):
    rng = np.random.default_rng(seed)
    if scenario == "period_separated_intercepts":
        y = rng.vonmises(0.6, 8, 160)
        prior_mu, prior_sigma = 0, 20
        initial = [0.6 + 2 * np.pi * turn for turn in (-1, 0, 1, 2)]
    else:
        errors = rng.vonmises(0, 8, 80)
        y = wrap(np.pi + np.concatenate([errors, -errors]))
        prior_mu, prior_sigma = np.pi, 0.5
        initial = [np.pi + shift for shift in (-0.1, 0.1, -0.2, 0.2)]
    with pm.Model():
        intercept = pm.Normal("Intercept", prior_mu, prior_sigma)
        pm.Deterministic("mu", (intercept + np.pi) % (2 * np.pi) - np.pi)
        kappa = pm.Exponential("kappa", lam=1 / 8)
        pm.VonMises("y", mu=intercept, kappa=kappa, observed=y)
        idata = pm.sample(
            draws=draws,
            tune=tune,
            chains=4,
            cores=1,
            init="adapt_diag",
            initvals=[{"Intercept": value, "kappa": 8} for value in initial],
            random_seed=seed,
            target_accept=0.9,
            progressbar=False,
            compute_convergence_checks=False,
        )
    values = idata.posterior.Intercept.values
    return {
        "scenario": scenario,
        "seed": seed,
        "chains": 4,
        "draws_per_chain": draws,
        "warmup_per_chain": tune,
        "raw_rank_rhat": float(az.rhat(values)),
        "divergences": int(idata.sample_stats.diverging.sum()),
        "nuts_steps": int(idata.sample_stats.n_steps.sum()),
        **diagnose(values),
        "mixed_api": validate_mixed_api(idata),
    }


def write_report(result, output):
    maximum_draws = max(row["draws_per_chain"] for row in result["summaries"])
    lines = [
        "# Circular diagnostic validation",
        "",
        "This reproducible study compares the joint covariance diagnostic with separate",
        "ArviZ rank R-hats on sine and cosine. The maximum and arithmetic mean of those",
        "two R-hats are included only as exploratory comparisons; neither is a public",
        "aggregate API or a separately calibrated convergence criterion.",
        "",
        "See [the mathematical and API guide](../diagnostics.md) and",
        "[the script](https://github.com/jgraving/unwrap/blob/feat/circular-rhat/examples/circular_rhat_validation.py).",
        "",
        "## Design",
        "",
        f"The simulation has {len(result['rows'])} four-chain draw sets: ten scenarios,",
        f"{result['replicates']} seeds, and lengths {result['draw_lengths']} per chain.",
        "Most controls are independent draws. The stationary AR(1) control uses",
        "autocorrelation 0.95 and starts in its stationary distribution. Antipodal",
        "controls balance opposite directions within each half-chain, exposing a",
        "known first-moment limitation without singular within-chain covariance.",
        "All data and sampler seeds are recorded in the JSON results.",
        "",
        "The 1.01 exceedance rates in the JSON are descriptive. In particular, a",
        "large value for finite autocorrelated draws from one stationary target is",
        "evidence of limited finite-run agreement, not proof of a different target.",
        "This study does not calibrate a universal threshold or prove convergence.",
        "",
        f"## Simulation medians at {maximum_draws:,} draws per chain",
        "",
        "| Scenario | Joint circular | Rank sine | Rank cosine | Component max | Component mean | Rank wrapped angle |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in result["summaries"]:
        if row["draws_per_chain"] == maximum_draws:
            metrics = " | ".join(f"{row['metrics'][metric]['median']:.4f}" for metric in METRICS)
            lines.append(f"| {row['scenario']} | {metrics} |")
    lines += [
        "",
        "Whole-turn offsets leave the angular diagnostics unchanged. Joint circular",
        "R-hat detects the tested mean-direction, concentration, and within-chain",
        "shift differences. Its first-moment limitation is visible for balanced",
        "antipodal pairs: different angular distributions can share the same mean",
        "unit vector. The sine/cosine checks detect some such cases but can themselves",
        "miss distributions with matching projection marginals. Neither approach",
        "dominates across these constructed examples.",
        "",
        "The small-direction-difference case illustrates why reporting both",
        "components is useful: averaging a larger sine R-hat with a near-unit",
        "cosine R-hat dilutes the visible disagreement.",
        "",
        "## Rotation comparison",
        "",
        "The same draws from the matching-projection-marginals case are rotated",
        "together. The joint covariance statistic is invariant; separate ranked",
        "projections can change because a rotation changes the scalar projections.",
        "",
        "| Rotation | Joint circular | Rank sine | Rank cosine |",
        "| --- | ---: | ---: | ---: |",
    ]
    for row in result["rotation_check"]:
        lines.append(
            f"| {row['rotation_radians']:.4f} | {row['joint_circular']:.4f} | "
            f"{row['rank_sin']:.4f} | {row['rank_cos']:.4f} |"
        )
    lines += ["", "## Actual PyMC sampling", ""]
    if result["mcmc"]:
        lines += [
            "Two intercept-only von Mises models use four NUTS chains, 500 warmup",
            "iterations and 500 retained draws per chain. Concentration is inferred",
            "with an Exponential prior of mean 8. One initializes the intercept in different whole-turn modes under",
            "a broad Normal(0, 20) prior; the other uses a Normal(π, 0.5) prior and",
            "symmetric data centered on the angular cut. These are smoke experiments",
            "for diagnostic behavior, not a broad sampler benchmark.",
            "",
            "| Model | Raw rank R-hat | Joint circular | Rank sine | Rank cosine | Divergences |",
            "| --- | ---: | ---: | ---: | ---: | ---: |",
        ]
        for row in result["mcmc"]:
            lines.append(
                f"| {row['scenario']} | {row['raw_rank_rhat']:.4f} | "
                f"{row['joint_circular']:.4f} | {row['rank_sin']:.4f} | "
                f"{row['rank_cos']:.4f} | {row['divergences']} |"
            )
        lines += [
            "",
            "The period-separated real intercepts remain in different modes. Their",
            "directions agree after accounting for periodicity. These are different",
            "quantities to diagnose, and the raw coefficient diagnostic is retained.",
            "",
            "## Whole-posterior API checks",
            "",
            "Each fit retains the real `Intercept`, its wrapped fitted direction `mu`,",
            "and the positive concentration `kappa`. Both angular variables can be",
            "declared in one call while `kappa` retains ordinary ArviZ behavior:",
            "",
            "```python",
            "from unwrap.diagnostics import rhat, summary",
            "",
            'angular = ["Intercept", "mu"]',
            "rhats = rhat(idata, circ_var_names=angular)",
            "table = summary(idata, circ_var_names=angular)",
            "```",
            "",
            "Assertions on both actual posteriors verify exact equality with ArviZ",
            "for ordinary `kappa` R-hat and summary cells. Angular `r_hat` cells equal",
            "the joint circular diagnostic; the additional `_sin` and `_cos` R-hat,",
            "bulk/tail ESS, and mean/SD MCSE columns equal ArviZ diagnostics of the",
            "corresponding projections. Projected MCSEs describe the projections,",
            "not an angular mean or SD. The raw-angle ESS and `mcse_sd` cells are",
            "blank for angular rows; projection columns are blank for ordinary rows.",
            "",
            "In the tested ArviZ 0.x version, its own `circ_var_names` option changes",
            "circular descriptive summaries but leaves R-hat and ESS computed on",
            "the raw angles. The wrapper replaces those R-hat cells and reports",
            "explicitly labeled projection diagnostics. The JSON records the",
            "comparison and every mixed summary cell, with blanks represented as null.",
            "",
            "| Model | Variable | ArviZ rank R-hat | Mixed R-hat | Rank sine | Rank cosine |",
            "| --- | --- | ---: | ---: | ---: | ---: |",
        ]
        for fit in result["mcmc"]:
            for row in fit["mixed_api"]["rows"]:
                values = row["mixed_summary"]
                sine = "—" if values["r_hat_sin"] is None else f"{values['r_hat_sin']:.4f}"
                cosine = "—" if values["r_hat_cos"] is None else f"{values['r_hat_cos']:.4f}"
                lines.append(
                    f"| {fit['scenario']} | {row['variable']} | "
                    f"{row['raw_arviz_rhat']:.4f} | {row['mixed_rhat']:.4f} | {sine} | {cosine} |"
                )
    else:
        lines.append("Sampling was skipped for this run with `--skip-mcmc`.")
    lines += [
        "",
        "## Reproduce",
        "",
        "```bash",
        "python examples/circular_rhat_validation.py --output docs/validation",
        "```",
        "",
        "Use `--replicates` to change the simulation repetitions or `--skip-mcmc`",
        "to run only the synthetic controls. The JSON records every run and dependency",
        "version. No simulated data are used to alter the diagnostic implementation.",
    ]
    (output / "circular_rhat.md").write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("docs/validation"))
    parser.add_argument("--replicates", type=int, default=20)
    parser.add_argument("--skip-mcmc", action="store_true")
    args = parser.parse_args()
    if args.replicates < 1:
        parser.error("--replicates must be positive")
    args.output.mkdir(parents=True, exist_ok=True)
    result = {
        "versions": {
            name: importlib.metadata.version(name)
            for name in ("unwrap", "arviz", "numpy", "scipy", "pymc", "pytensor")
        },
        "replicates": args.replicates,
        "draw_lengths": [200, 1000, 4000],
        "rows": [],
    }
    for draws in result["draw_lengths"]:
        for replicate in range(args.replicates):
            seed = 51000 + replicate
            for scenario in SCENARIOS:
                result["rows"].append(
                    {
                        "scenario": scenario,
                        "draws_per_chain": draws,
                        "seed": seed,
                        **diagnose(simulate(scenario, draws, seed)),
                    }
                )
        print(f"Finished simulations at {draws} draws per chain", flush=True)
    result["summaries"] = summarize(result["rows"])
    rotation_draws = simulate("same_projection_marginals", 4000, 51000)
    result["rotation_check"] = [
        {"rotation_radians": rotation, **diagnose(rotation_draws + rotation)}
        for rotation in (0.0, np.pi / 4, np.pi / 2)
    ]
    result["mcmc"] = (
        []
        if args.skip_mcmc
        else [
            fit_intercept("period_separated_intercepts", 62000),
            fit_intercept("at_branch_cut", 62001),
        ]
    )
    (args.output / "circular_rhat.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )
    write_report(result, args.output)


if __name__ == "__main__":
    main()
