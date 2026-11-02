"""Build notebooks/goldshtein_bats.ipynb (study-PR layout).

Rewrites the 96-cell annotated base notebook into thin cells: markdown +
data + model + evaluation stay inline; the heavy preprocessing, post-hoc
statistics and figure code live in the sibling scripts
(goldshtein_bats_preprocess.py, goldshtein_bats_posthoc.py,
goldshtein_bats_plots.py). Every figure cell is a comment + a single
``plots.plot_*(save_path=...)`` call + ``plt.show()``.

Base cell map (c## = base notebook index):
  keep verbatim: all markdown; c03 c04 c09 c12 c14 c16 c17 c45 c46 c47 c48
                 c57
  drop (superseded by the scripts): c19 c25 c27 c62 c68 c78 c88 c92
  replace: c01 setup; c06 c08 c11 preprocess calls; c32 c37 c41 thin
           posthoc compute; c43 vM model (fit call ported to the installed
           bambi: no idata_kwargs, log-likelihood attached afterwards via
           pm.compute_log_likelihood + posthoc.graph_free_rvs, house
           divergence print added); c55 vMU model (base verbatim except the
           house divergence print after the fit call); c54 sigmoid
           percentiles; figure slots (see FIGURE_CELLS below).
"""
import json
import shutil
from pathlib import Path

NB_PATH = Path("/tmp/opencode/wt-unified/notebooks/goldshtein_bats.ipynb")
BAK = Path("/tmp/opencode/goldshtein_bats_base.ipynb.bak")
if not BAK.exists():
    shutil.copy2(NB_PATH, BAK)

base = json.load(open(BAK))
base_cells = {f"c{i:02d}": c for i, c in enumerate(base["cells"])}


def code(text):
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": text.splitlines(keepends=True),
    }


def keep_md(i):
    c = base_cells[f"c{i:02d}"]
    assert c["cell_type"] == "markdown", f"c{i:02d} is not markdown"
    out = dict(c)
    out["metadata"] = {}
    return out


def keep_code(i):
    c = base_cells[f"c{i:02d}"]
    assert c["cell_type"] == "code", f"c{i:02d} is not code"
    out = dict(c)
    out["metadata"] = {}
    out["execution_count"] = None
    out["outputs"] = []
    return out


# --------------------------------------------------------------------------
# c01 -> house setup
SETUP = """\
%matplotlib inline
from importlib.metadata import version
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import arviz as az

import unwrap
from unwrap.datasets import load_bats_info, load_bats_tracks, bats_lake_outline_path

# The preprocessing, post-hoc and figure code lives in sibling scripts next
# to this notebook (goldshtein_bats_preprocess.py,
# goldshtein_bats_posthoc.py, goldshtein_bats_plots.py); make them
# importable regardless of the kernel's working directory (the kernel may
# run from notebooks/ or from the repository root).
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(unwrap.__file__)), "notebooks"))

import goldshtein_bats_preprocess as pre
import goldshtein_bats_posthoc as posthoc
import goldshtein_bats_plots as plots

plots.set_style()

SEED = 123
{package: version(package) for package in (
    "unwrap", "pymc", "arviz", "bambi",
    "numpy", "pandas", "scipy", "matplotlib",
)}
"""

# --------------------------------------------------------------------------
# c06
C06 = """\
# Resample every tag's lat/lon trajectory to a constant 60 m step length;
# gaps (a jump > gap_factor * spacing_m = 600 m, or a time gap > 2 min)
# split the track into contiguous segments flagged by ``gap_before``.
rediscretized_ll_df = pre.rediscretize_latlon_parallel(
    raw_track_df,
    spacing_m=60.0,
    n_workers=16,
    include_last=False,
    gap_factor=10.0,
    gap_time_s=120.0,
    time_unit="days",
)

print(rediscretized_ll_df.head())
print("rows:", len(rediscretized_ll_df), "tags:", rediscretized_ll_df["tag"].nunique())
print("n segments total:", rediscretized_ll_df["segment_id"].nunique() if "segment_id" in rediscretized_ll_df else None)
print("gap markers:", int(rediscretized_ll_df["gap_before"].sum()) if "gap_before" in rediscretized_ll_df else None)

step_stats = pre.check_step_lengths_latlon(rediscretized_ll_df)
print(step_stats.describe(include="all"))
"""

# --------------------------------------------------------------------------
# c08
C08 = """\
# One metadata row per tag (treatment + destination label), merged back
# onto the resampled steps so every step knows its treatment and goal.
rediscretized_ll_w_meta_df = pre.merge_tag_metadata(rediscretized_ll_df, meta_df)

print("resampled rows:", len(rediscretized_ll_w_meta_df))
print("missing treatment rows:",
      rediscretized_ll_w_meta_df["treatment"].isna().sum(),
      " (tags missing:",
      rediscretized_ll_w_meta_df.loc[rediscretized_ll_w_meta_df["treatment"].isna(), "tag"].nunique(),
      ")"
)
"""

# --------------------------------------------------------------------------
# c11
C11 = """\
# Distance + bearing to each animal's destination (lake: negative
# gradient of a harmonic distance field; roost: bearing to the roost
# point), plus the goal-relative angle ``flightAngDiffDest`` used as the
# circular response (0 deg = aimed at the goal).
df_harm = pre.add_dist_and_angle_to_destination(
    rediscretized_ll_w_meta_df,
    LAKE_GEOJSON=bats_lake_outline_path(),
    METRIC_EPSG=32636,
    ROOST_LAT=33.112373,
    ROOST_LON=35.584135,
    lake_grid_res_m=20.0,
    lake_domain_margin_m=None,
    lake_jacobi_tol=1e-5,
    lake_field="harmonic",
)
df_harm[
    ["tag", "dest",
     "dist_to_lake_m", "dist_to_lake_grad_m",
     "dist_to_dest_m", "dist_to_dest_grad_m",
     "step_len_m", "step_bearing_deg",
     "bearing_to_dest_deg", "flightAngDiffDest_deg"]
].head(10)

rediscretized_ll_w_meta_df2 = df_harm
"""

# --------------------------------------------------------------------------
# thin posthoc compute cells
C32 = """\
# Per-tag descriptive concentration: mean resultant length -> kappa_hat
# (Fisher), full trajectory, grouped by treatment.
kappa_df = posthoc.kappa_by_tag_full_trajectory(track_w_treat)
kappa_df.head()
"""

C37 = """\
# Per-tag alignment of the animal's mean bearing to the goal direction
# (0 deg), full trajectory, grouped by treatment.
align_df = posthoc.alignment_by_tag_full_trajectory(track_w_treat)
align_df.head()
"""

C41 = """\
# Modeling table for the vM fit: the usable rows, the wrapped
# goal-relative heading (HeadingRadians in (-pi, pi]), categorical
# tag/treatment and the standardized log distance (prints the design
# summary).
df_track = posthoc.build_modeling_table(track_w_treat)
"""

C54 = """\
# Sigmoid-space percentiles of a Normal on the logit scale:
# P = sigmoid(Normal(logit(0.9), 1)) — how a posterior mixture weight p
# with that prior/uncertainty reads on the [0, 1] scale.
from scipy.special import logit

qs, pqs = posthoc.sigmoid_normal_percentiles(logit(0.9), 1.0)
for q, p in zip(qs, pqs):
    print(f"{q:>5.2%}: {p:.6f}")
"""

# --------------------------------------------------------------------------
# c43 -> vM Bambi model (kept inline; the fit call is ported to the
# installed bambi, which takes no idata_kwargs: the per-observation
# log-likelihood is attached in place afterwards).
_src = "".join(base_cells["c43"]["source"])
_fit_old = '''idata_treatment = model.fit(tune=2000, draws=2000, target_accept=0.9, idata_kwargs={"log_likelihood": True}
)'''
_fit_new = '''idata_treatment = model.fit(tune=2000, draws=2000, target_accept=0.9, random_seed=SEED, progressbar=False)

# Store the per-observation log-likelihood needed for the PSIS-LOO
# comparison below (the pymc 6 replacement for
# idata_kwargs={"log_likelihood": True}): bambi's own
# Model.compute_log_likelihood applies the vonmises *mean* link to every
# submodel, so it would wrap the kappa linear predictor into an angle
# (bambinos/bambi#1084); the PyMC graph evaluation uses each submodel's
# actual link (posthoc.graph_free_rvs supplies the non-centered *_offset
# free RVs that bambi does not store).
idata_treatment = pm.compute_log_likelihood(
    posthoc.graph_free_rvs(idata_treatment, model),
    model=model.backend.model,
    progressbar=False,
)

# Report the divergence count first: R-hat is only meaningful when it is 0.
print("Divergences:", int(idata_treatment.sample_stats.diverging.sum()))'''
assert _src.count(_fit_old) == 1, "c43 fit call not found verbatim"
C43 = _src.replace(_fit_old, _fit_new)
# the graph evaluation needs pymc in the cell's namespace
C43 = C43.replace("import arviz as az\nimport matplotlib.pyplot as plt",
                  "import arviz as az\nimport pymc as pm\nimport matplotlib.pyplot as plt")

# the vMU cell keeps the base verbatim EXCEPT the house divergence print
# (zebra/density pattern) added right after the fit call
_MODEL55_SRC = "".join(keep_code(55)["source"])
_VMU_FIT_OLD = '''idata_mix_full, model_full = fit_vmum_model_full(
    include_treatment_in_p=True,
    include_treatment_in_mu=True,
    include_treatment_in_kappa=True,
    seed=123,
)'''
_VMU_FIT_NEW = '''idata_mix_full, model_full = fit_vmum_model_full(
    include_treatment_in_p=True,
    include_treatment_in_mu=True,
    include_treatment_in_kappa=True,
    seed=123,
)

# Report the divergence count first: R-hat is only meaningful when it is 0.
print("Divergences:", int(idata_mix_full.sample_stats.diverging.sum()))'''
assert _MODEL55_SRC.count(_VMU_FIT_OLD) == 1, "c55 fit call not found verbatim"
C55 = code(_MODEL55_SRC.replace(_VMU_FIT_OLD, _VMU_FIT_NEW))

# --------------------------------------------------------------------------
# figure cells. Each: comment + one plots.* call + plt.show().
F21 = """\
# Full-track destination panels by treatment x destination (basemap at
# zoom 15 with retina tiles, no repeated treatment legend).
plots.plot_destination_panels(
    track_w_treat,
    meta_df,
    basemap_zoom=15,
    basemap_retina=True,
    show_treatment_legend=False,
    save_path="figures/bats/rediscretized_destination_panels.pdf",
)
plt.show()
"""

F23 = """\
# Compact version: one panel per treatment / destination combination,
# small figure with larger relative text (trajectory layers rasterized
# inside the saved PDF).
plots.plot_destination_panels_compact(
    track_w_treat,
    meta_df,
    fig_w=2.5,
    fig_h=3.8,
    title_fs=9,
    legend_fs=6,
    save_path="figures/bats/rediscretized_destination_panels_compact.pdf",
)
plt.show()
"""

F29 = """\
# Goal-approach direction fields feeding bearing_to_dest: the
# smoothed-distance field over the lake panel and the radial-to-roost
# field over the roost panel (scale/compass UI off for the export).
plots.plot_destination_flow_fields(
    track_w_treat,
    meta_df,
    lake_field="smoothed_distance",
    lake_grid_res_m=200.0,
    lake_smooth_sigma_pix=1.0,
    roost_grid_res_m=200.0,
    quiver_scale_lake=30,
    quiver_scale_roost=30,
    show_ui_scale_compass=False,
    save_path="figures/bats/destination_flow_fields_smoothed_distance.pdf",
)
plt.show()
"""

F31 = """\
# Polar KDE of the goal-relative bearings by treatment x destination
# (von Mises kernel, kappa = 15).
plots.plot_per_tag_polar_kde(
    track_w_treat,
    kappa=15,
    r0=1.0,
    r_scale=1.0,
    save_path="figures/bats/polar_vonmises_kde_all.pdf",
)
plt.show()
"""

F34 = """\
# Estimated population concentration kappa_hat by treatment (median with
# percentile bands over the per-tag estimates).
plots.plot_kappa_hat(kappa_df,
                     save_path="figures/bats/kappa_hat_by_treatment_median_percentiles.pdf")
plt.show()
"""

F36 = """\
# Mean resultant length R by treatment (median with percentile bands
# over the per-tag estimates).
plots.plot_rbar(kappa_df,
                save_path="figures/bats/rbar_by_treatment_median_percentiles.pdf")
plt.show()
"""

F39 = """\
# Alignment of each animal's mean bearing to the goal direction by
# treatment (median with percentile bands).
plots.plot_alignment(align_df,
                     save_path="figures/bats/alignment0_by_treatment_median_percentiles.pdf")
plt.show()
"""

F50 = r"""
# Treatment coefficients of the vM fit: the mu (mean direction) and
# kappa (concentration) submodels with 95% HDI, control intercept on
# top.
plots.plot_treatment_coefficients_mu_kappa(idata_treatment, df_track,
                                           save_path="figures/bats/treatment_coefficients_mu_kappa.pdf")
plt.show()
"""

F52 = r"""
# Posterior predictive density at the goal direction (0 deg) by
# treatment for the vM baseline (mean + 95% HDI over posterior draws).
plots.plot_treatment_pred_density_at0_mu_kappa(idata_treatment, df_track,
                                               save_path="figures/bats/treatment_pred_density_at0_mu_kappa.pdf")
plt.show()
"""

C59A = """\
# PSIS-LOO ELPD of both fitted models, best first (per-observation
# leave-one-out log-likelihoods from each model's idata).
models = {
    "von Mises-Uniform": idata_mix_full,
    "von Mises": idata_treatment,
}

res_total, best_model, loo_cache = posthoc.compute_psisloo_elpd_summary(
    models=models,
    df_fit=df_track,
    var_name=None,
    pareto_thresh=0.7,
    order_best_first=True,
)

print("Best model:", best_model)
print(res_total)
"""

F59 = """\
# Total PSIS-LOO ELPD by model (vMU vs vM), best first.
plots.plot_total_elpd(
    res_total,
    title=None,
    xlabel="ELPD",
    vline_at_best=False,
    save_path="figures/bats/psisloo_total_elpd_all_models.pdf",
)
plt.show()
"""

C61A = """\
# Delta ELPD (vMU - vM), total and stratified by treatment, with the SE
# of the total.
delta_df = posthoc.elpd_diff_total_and_by_treatment(
    idata_vmu=idata_mix_full,      # von Mises-Uniform
    idata_vm=idata_treatment,      # von Mises
    df_fit=df_track,
    treat_col="treatment",
    treat_name_col="treatmentName",
    var_name=None,
)

treatment_order = (
    df_track["treatment"].cat.categories.tolist()
    if hasattr(df_track["treatment"], "cat")
    else sorted(delta_df.loc[delta_df["row_type"] == "treatment", "treatment"].dropna().unique())
)
delta_df
"""

F61 = """\
# Delta ELPD (vMU - vM) total row on top of the within-treatment rows.
plots.plot_delta_total_and_by_treatment(
    delta_df,
    treatment_order=treatment_order,
    title=None,
    xlabel="ELPD difference",
    vline_at_zero=True,
    save_path="figures/bats/psisloo_delta_total_and_by_treatment_vMU_minus_vM.pdf",
)
plt.show()
"""

F64 = """\
# Combined two-panel figure: total ELPD by model (top) next to the
# within-treatment delta ELPD bracket (bottom).
plots.plot_psisloo_combined(res_total, delta_df, treatment_order,
                            save_path="figures/bats/psisloo_combined.pdf")
plt.show()
"""

C66A = """\
# Pointwise delta ELPD (vMU - vM) per animal, grouped by treatment —
# the vMU advantage can be concentrated in a subset of animals.
res_ind = posthoc.delta_by_individual_treatment(
    idata_vmu=idata_mix_full,
    idata_vm=idata_treatment,
    df_fit=df_track,
    individual_col="tag",
    treat_col="treatment",
    treat_name_col="treatmentName",
    var_name=None,
)

treatment_order = (
    df_track["treatment"].cat.categories.tolist()
    if hasattr(df_track["treatment"], "cat")
    else sorted(res_ind["treatment"].unique())
)
"""

F66 = """\
# Pointwise delta ELPD (vMU - vM) by animal and by treatment (animals
# sorted by their own delta ELPD within each panel).
plots.plot_delta_subplots_by_treatment(
    res_ind,
    treatment_order=treatment_order,
    save_path="figures/bats/psisloo_delta_by_individual_by_treatment_vMU_minus_vM.pdf",
)
plt.show()
"""

F69 = """\
# Per-treatment population polar densities of the goal-relative
# bearings for the vMU fit: the directed + undirected mixture.
plots.plot_mixture_polar_pdfs(idata_mix_full, df_track,
                              save_path="figures/bats/vmum_mixture_full_polar_pdfs_population_only.pdf")
plt.show()
"""

C71A = """\
# Population-level vMU summaries by treatment: mixture weight p,
# goal-directed mean direction and concentration, and the predictive
# density / mass in the goal window (posterior mean + 95% HDI).
res_full = posthoc.summarize_pop_p_kappa_mu_density0_mass45_alignment0_by_treatment(
    idata_mix_full=idata_mix_full,
    df_fit=df_track,
    treat_col="treatment",
    treat_name_col="treatmentName",
    hdi_prob=0.95,
    max_draws=2000,
    kappa_uniform=posthoc.KAPPA_UNIFORM,
    mass_half_width_deg=45.0,
)

# Shared x-limit so the mixture and directed-density panels are
# directly comparable.
dens_xmax = float(max(res_full["dens0_mix_hi"].max(), res_full["dens0_dir_hi"].max()) * 1.05)
dens_xlim = (0.0, dens_xmax)
"""

F71_1 = """\
# Mixture weight p (probability a step is goal-directed) by treatment;
# dashed line at the uniform 0.5.
plots.plot_mean_hdi_horizontal(
    res_full, "p_mean", "p_lo", "p_hi", "Mixture weight",
    xlim=(0, 1), vline_x=0.5,
    show_treatment_labels=False, show_treatment_ticks=False,
    save_path="figures/bats/fullmodel_mixture_weight_by_treatment_mean_hdi.pdf",
)
plt.show()
"""

F71_2 = r"""
# Goal-directed concentration kappa by treatment.
plots.plot_mean_hdi_horizontal(
    res_full, "kappa_mean", "kappa_lo", "kappa_hi", "Concentration " + rf"($\kappa$)",
    show_treatment_labels=False, show_treatment_ticks=False,
    save_path="figures/bats/fullmodel_kappa_by_treatment_mean_hdi.pdf",
)
plt.show()
"""

F71_3 = """\
# Predictive density at the goal direction (full mixture) by treatment;
# dashed line at the exact uniform density 1/(2 pi).
plots.plot_mean_hdi_horizontal(
    res_full, "dens0_mix_mean", "dens0_mix_lo", "dens0_mix_hi",
    "Predictive density (mixture)",
    xlim=dens_xlim, vline_x=1.0 / (2.0 * np.pi),
    show_treatment_labels=False, show_treatment_ticks=False,
    save_path="figures/bats/fullmodel_predictive_density_at_0_mixture_by_treatment_mean_hdi.pdf",
)
plt.show()
"""

F71_4 = """\
# Predictive density at the goal direction (directed component only) by
# treatment; the uniform line stays as the no-preference reference.
plots.plot_mean_hdi_horizontal(
    res_full, "dens0_dir_mean", "dens0_dir_lo", "dens0_dir_hi",
    "Predictive density (directed component)",
    xlim=dens_xlim, vline_x=1.0 / (2.0 * np.pi),
    show_treatment_labels=False, show_treatment_ticks=False,
    save_path="figures/bats/fullmodel_predictive_density_at_0_directed_component_by_treatment_mean_hdi.pdf",
)
plt.show()
"""

F71_5 = """\
# Predictive mass within +/-45 deg of the goal direction by treatment;
# dashed line at the uniform 0.25.
plots.plot_mean_hdi_horizontal(
    res_full, "mass45_mean", "mass45_lo", "mass45_hi", "Predictive mass (±45°)",
    xlim=(0, 1), vline_x=0.25,
    show_treatment_labels=False, show_treatment_ticks=False,
    save_path="figures/bats/fullmodel_predictive_mass_pm45deg_by_treatment_mean_hdi.pdf",
)
plt.show()
"""

F71_6 = """\
# Alignment to the goal direction E[cos(theta)] under the mixture by
# treatment.
plots.plot_mean_hdi_horizontal(
    res_full, "align0_mean", "align0_lo", "align0_hi", "Alignment",
    xlim=(0, 1),
    show_treatment_labels=False, show_treatment_ticks=False,
    save_path="figures/bats/fullmodel_alignment_to_zero_by_treatment_mean_hdi.pdf",
)
plt.show()
"""

F73 = """\
# Combined summary figure: delta ELPD (top) next to the population vMU
# parameters in a 2x2 block (mean angle, mixture weight, MRL,
# predictive density at the goal).
plots.plot_psisloo_plus_params_2x2(
    delta_df, treatment_order, idata_mix_full, df_track,
    save_path="figures/bats/psisloo_plus_params_2x2.pdf",
)
plt.show()
"""

F75 = """\
# Matched two-panel comparison of the predictive density at the goal
# direction: vMU directed component (top) vs vM (bottom); the uniform
# dashed line is the no-preference reference.
plots.plot_predictive_density_vmu_directed_vs_vm(
    res_full, df_track, idata_treatment,
    save_path="figures/bats/predictive_density_two_model_comparison_matched_vmu_directed_vs_vm.pdf",
)
plt.show()
"""

F76 = """\
# Three-panel summary: total ELPD by model, delta ELPD by treatment,
# and the vMU mixture predictive density at the goal direction.
plots.plot_psisloo_plus_predictive_density_3panel(
    res_total, delta_df, treatment_order, res_full,
    save_path="figures/bats/psisloo_plus_predictive_density_3panel.pdf",
)
plt.show()
"""

C79A = """\
# Per-animal (tag-level) vMU summaries from the model terms: mixture
# weight p, goal-directed kappa and density at the goal, alignment to
# 0 deg — plus the population-only table that draws the reference band.
res_ind, treatment_order = posthoc.summarize_taglevel_from_model_terms(
    idata_mix_full=idata_mix_full,
    df_fit=df_track,
    treat_col="treatment",
    treat_name_col="treatmentName",
    individual_col="tag",
    hdi_prob=0.95,
    max_draws=2000,
    kappa_uniform=posthoc.KAPPA_UNIFORM,
)

pop_df = posthoc.pop_summary_by_treatment_from_model_terms(
    idata_mix_full=idata_mix_full,
    df_fit=df_track,
    treat_col="treatment",
    treat_name_col="treatmentName",
    hdi_prob=0.95,
    max_draws=2000,
    kappa_uniform=posthoc.KAPPA_UNIFORM,
)

# reference values for the dashed lines
uniform_density_at_0 = 1.0 / (2.0 * np.pi)
mixture_weight_ref = 0.5
"""

F79_1 = """\
# Per-animal mixture weight by treatment (gray population band + solid
# population mean, dashed line at 0.5; animals sorted within each
# panel).
plots.plot_individual_subplots_by_treatment(
    res_ind, treatment_order,
    value_col="p_mean", lo_col="p_lo", hi_col="p_hi",
    xlabel="Mixture weight",
    pop_df=pop_df, pop_mean_col="p_mean", pop_lo_col="p_lo", pop_hi_col="p_hi",
    xlim=(0, 1),
    ref_vline_x=mixture_weight_ref,
    save_path="figures/bats/fullmodel_taglevel_mixture_weight_by_treatment_mean_hdi.pdf",
)
plt.show()
"""

F79_2 = r"""
# Per-animal goal-directed concentration kappa by treatment.
plots.plot_individual_subplots_by_treatment(
    res_ind, treatment_order,
    value_col="kappa_mean", lo_col="kappa_lo", hi_col="kappa_hi",
    xlabel="Concentration " + rf"($\kappa$)",
    pop_df=pop_df, pop_mean_col="kappa_mean", pop_lo_col="kappa_lo", pop_hi_col="kappa_hi",
    save_path="figures/bats/fullmodel_taglevel_kappa_by_treatment_mean_hdi.pdf",
)
plt.show()
"""

F79_3 = """\
# Per-animal predictive density at the goal direction by treatment
# (dashed line at the uniform density).
plots.plot_individual_subplots_by_treatment(
    res_ind, treatment_order,
    value_col="dens0_mean", lo_col="dens0_lo", hi_col="dens0_hi",
    xlabel="Predictive density",
    pop_df=pop_df, pop_mean_col="dens0_mean", pop_lo_col="dens0_lo", pop_hi_col="dens0_hi",
    ref_vline_x=uniform_density_at_0,
    save_path="figures/bats/fullmodel_taglevel_predictive_density_at_0_by_treatment_mean_hdi.pdf",
)
plt.show()
"""

F79_4 = """\
# Per-animal alignment to the goal direction by treatment.
plots.plot_individual_subplots_by_treatment(
    res_ind, treatment_order,
    value_col="align0_mean", lo_col="align0_lo", hi_col="align0_hi",
    xlabel="Alignment",
    pop_df=pop_df, pop_mean_col="align0_mean", pop_lo_col="align0_lo", pop_hi_col="align0_hi",
    xlim=(0, 1),
    save_path="figures/bats/fullmodel_taglevel_alignment_to_zero_by_treatment_mean_hdi.pdf",
)
plt.show()
"""

C80A = """\
# Population-level PPC by treatment: hierarchical bootstrap of the
# observed goal-relative headings vs the same summaries simulated from
# the vMU posterior (tag random effects included); the heavy machinery
# runs in goldshtein_bats_posthoc.
N_KDE_POINTS = 360       # evaluation grid for the KDE
N_HIST_BINS = 16
N_BOOT_OBS = 1000        # hierarchical bootstrap replicates (observed)
N_PPC_DRAWS = 4000       # posterior draws for the PPC
PPC_BOOT_PER_DRAW = 1
KDE_KAPPA = 30.0
N_JOBS = max(1, (os.cpu_count() or 2) - 1)
RNG_SEED = 123

df_ppc = df_track.reset_index(drop=True).copy()  # index must be 0..N-1 for PPC slicing
treat_levels, groups = posthoc.build_groups(
    df_ppc, y_col="HeadingRadians", tag_col="tag", treat_col="treatment"
)
treatment_label_map = (
    df_ppc.groupby("treatment", observed=False)["treatmentName"].first().to_dict()
)

angles = np.linspace(-np.pi, np.pi, N_KDE_POINTS, endpoint=False)
bin_edges = np.linspace(-np.pi, np.pi, N_HIST_BINS + 1)
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

obs_boot = posthoc.bootstrap_observed(
    groups=groups, angles=angles, bin_edges=bin_edges,
    n_boot=N_BOOT_OBS, kde_kappa=KDE_KAPPA, n_jobs=N_JOBS, seed=RNG_SEED,
)

ppc = posthoc.posterior_predictive_summaries(
    idata=idata_mix_full, df=df_ppc, treat_levels=treat_levels,
    angles=angles, bin_edges=bin_edges, kde_kappa=KDE_KAPPA,
    n_draws=N_PPC_DRAWS, n_boot_within_draw=PPC_BOOT_PER_DRAW,
    n_jobs=N_JOBS, seed=RNG_SEED,
)
"""

F80 = """\
# Population PPC by treatment: observed hierarchical-bootstrap KDE vs
# posterior predictive KDE of the goal-relative bearings.
plots.plot_ppc_by_treatment(
    treat_levels, angles, bin_centers, obs_boot, ppc,
    treatment_label_map=treatment_label_map,
    save_path="figures/bats/bat_population_ppc_by_treatment.pdf",
)
plt.show()
"""

C82 = """\
# Individual PPC: per-animal bootstrap (observed) vs posterior KDE of
# the goal-directed component; the grid is built in posthoc and drawn
# in the plotting script.
ppc_grid = posthoc.individual_ppc(idata_mix_full, df_track, n_cols=8)
"""

F84 = """\
# Individual PPC grid: per-animal data bootstrap vs posterior
# predictive KDE of the goal-directed component.
plots.plot_individual_ppc_grid(ppc_grid,
                               save_path="figures/bats/bat_individual_ppc_bootstrap_vs_ppc_kde.pdf")
plt.show()
"""

C85 = """\
# Posterior probability that each step is goal-directed (p_goal):
# responsibilities of the two mixture components per step, summarized
# over 2000 draws with the exact uniform-circle likelihood; the p
# columns are appended to df_track and joined onto the tracks in the
# figures below.
df_track["_step"] = df_track.groupby("tag", observed=False).cumcount()
df_track, _ = posthoc.append_pgoal_post_to_df(
    df_track, idata_mix_full, y_col="HeadingRadians", max_draws=2000, hdi_prob=0.95
)

print(df_track[["tag", "treatment", "_step", "p_goal_post_mean", "p_goal_post_lo", "p_goal_post_hi"]].head())
"""

F87 = """\
# Tracks colored by the per-step goal-directed probability (viridis),
# one panel per treatment.
plots.plot_tracks_pgoal_viridis(
    track_w_treat,
    p_col="p_goal_post_mean",
    df_prob=df_track,
    prob_tag_col="tag",
    prob_step_col="_step",
    prob_p_col="p_goal_post_mean",
    gap_col="gap_before",
    gap_linestyle="--",
    gap_alpha=0.80,
    seg_lw=1.0,
    wspace=0.10,
    hspace=0.10,
    fig_w=8.0,
    fig_h=6,
    save_path="figures/bats/treatment_panels_tracks_pgoal.pdf",
)
plt.show()
"""

F89 = """\
# Compact version of the p_goal-colored tracks (smaller figure, bigger
# colorbar text, track layers rasterized inside the saved PDF).
plots.plot_tracks_pgoal_small_bigfont(
    track_w_treat,
    p_col="p_goal_post_mean",
    df_prob=df_track,
    prob_tag_col="tag",
    prob_step_col="_step",
    prob_p_col="p_goal_post_mean",
    gap_col="gap_before",
    gap_linestyle="--",
    gap_alpha=0.75,
    seg_lw=0.8,
    wspace=0.03,
    hspace=0.02,
    save_path="figures/bats/treatment_panels_tracks_pgoal_small_bigfont.pdf",
)
plt.show()
"""

C91A = """\
# Population-level p_goal along the trajectory, kernel-smoothed per tag
# then averaged per draw (per-tag truncation to the observed range,
# post-smoothing across x, masking where few tags contribute).
# df_track_with_p, resp and the curves below are also reused by the
# compact and per-animal cells further down.
STEP_LENGTH_M = 60.0
DIST_COL = "dist_to_dest_grad_m"
LAST_N_STEPS = 110

df_track_with_p, resp = posthoc.append_pgoal_post_to_df(
    df_track, idata_mix_full, y_col="HeadingRadians", max_draws=2000, hdi_prob=0.95
)
print(df_track_with_p[["p_goal_post_mean", "p_goal_post_lo", "p_goal_post_hi"]].head())

treat_levels = (
    list(df_track_with_p["treatment"].cat.categories)
    if isinstance(df_track_with_p["treatment"].dtype, pd.CategoricalDtype)
    else sorted(df_track_with_p["treatment"].dropna().unique().tolist())
)

# A) Distance-to-destination
_, out_dist_k = posthoc.population_curve_kernel_per_tag(
    df_track_with_p, resp,
    x_kind="distance",
    dist_col=DIST_COL,
    n_grid=2000,
    bandwidth=120.0,
    chunk_size=4000,
    hdi_prob=0.95,
    truncate_to_data_range=True,
    truncate_by_weight=True,
    den_frac=0.01,
    pop_post_smooth=True,
    pop_smooth_bandwidth=240.0,
    min_tags_at_x=2,
)

# B) Steps-from-end (x in meters)
_, out_steps_k = posthoc.population_curve_kernel_per_tag(
    df_track_with_p, resp,
    x_kind="steps_from_end",
    last_n_steps=LAST_N_STEPS,
    steps_x_units="m",
    step_length_m=STEP_LENGTH_M,
    n_grid=2000,
    bandwidth=120.0,
    chunk_size=4000,
    hdi_prob=0.95,
    truncate_to_data_range=True,
    truncate_by_weight=True,
    den_frac=0.01,
    pop_post_smooth=True,
    pop_smooth_bandwidth=240.0,
    min_tags_at_x=2,
)
"""

F91_1 = """\
# Population mean p_goal vs distance to the destination (x-axis
# inverted: approaching the goal moves right).
plots.plot_population_curve_smooth(
    out_dist_k, treat_levels,
    xlabel="Distance to destination (m)",
    ylabel="Mean directed-component probability",
    invert_x=True,
    save_path="figures/bats/p_goal_pop_kernel_distance_TRUNC_SMOOTH.pdf",
)
plt.show()
"""

F91_2 = """\
# Population mean p_goal vs cumulative path distance from the end of
# the flight (negative = steps before the end).
plots.plot_population_curve_smooth(
    out_steps_k, treat_levels,
    xlabel="Cumulative path distance from end (m)",
    ylabel="Mean directed-component probability",
    invert_x=False,
    xlim=(-LAST_N_STEPS * STEP_LENGTH_M, 0),
    x_tick_abs=True,
    save_path="figures/bats/p_goal_pop_kernel_distFromEndM_TRUNC_SMOOTH.pdf",
)
plt.show()
"""

F93_1 = r"""
# Compact (paper-width) version: p_goal vs distance to the destination
# (same curves as above).
plots.plot_population_curve_smooth_compact(
    out_dist_k, treat_levels,
    xlabel="Distance to destination (m)",
    ylabel="Directed-component\nprobability",
    invert_x=True,
    save_path="figures/bats/p_goal_pop_kernel_distance_TRUNC_SMOOTH_compact.pdf",
)
plt.show()
"""

F93_2 = r"""
# Compact version: p_goal vs cumulative path distance from the end of
# the flight (same curves as above).
plots.plot_population_curve_smooth_compact(
    out_steps_k, treat_levels,
    xlabel="Cumulative path distance from end (m)",
    ylabel="Directed-component\nprobability",
    invert_x=False,
    xlim=(-LAST_N_STEPS * STEP_LENGTH_M, 0),
    x_tick_abs=True,
    save_path="figures/bats/p_goal_pop_kernel_distFromEndM_TRUNC_SMOOTH_compact.pdf",
)
plt.show()
"""

C95A = """\
# Per-animal p_goal curves (kernel-smoothed per tag, truncated to each
# tag's data range): the faceted panels show when each bat commits to
# goal-directed travel.
_, out_tags_dist = posthoc.tag_curves_kernel_by_treatment(
    df_track_with_p, resp,
    x_kind="distance",
    dist_col=DIST_COL,
    n_grid=2000,
    bandwidth=120.0,
    chunk_size=4000,
    hdi_prob=0.95,
    truncate_to_data_range=True,
    truncate_by_weight=True,
    den_frac=0.01,
)

_, out_tags_steps = posthoc.tag_curves_kernel_by_treatment(
    df_track_with_p, resp,
    x_kind="steps_from_end",
    last_n_steps=LAST_N_STEPS,
    steps_x_units="m",
    step_length_m=STEP_LENGTH_M,
    n_grid=2000,
    bandwidth=120.0,
    chunk_size=4000,
    hdi_prob=0.95,
    truncate_to_data_range=True,
    truncate_by_weight=True,
    den_frac=0.01,
)
"""

F95_1 = """\
# Per-animal p_goal vs distance to the destination (four treatments
# per page; one saved PDF per page).
plots.plot_tag_curves_faceted_4panel(
    out_tags_dist, treat_levels,
    xlabel="Distance to destination (m)",
    ylabel="Directed-component probability",
    invert_x=True,
    fig_h=6,
    fig_w=8,
    save_path="figures/bats/p_goal_tag_kernel_distance_faceted_TRUNC.pdf",
)
plt.show()
"""

F95_2 = """\
# Per-animal p_goal vs cumulative path distance from the end of the
# flight (one saved PDF per page).
plots.plot_tag_curves_faceted_4panel(
    out_tags_steps, treat_levels,
    xlabel="Cumulative path distance from end (m)",
    ylabel="Directed-component probability",
    invert_x=False,
    xlim=(-LAST_N_STEPS * STEP_LENGTH_M, 0),
    x_tick_abs=True,
    fig_h=6,
    fig_w=8,
    save_path="figures/bats/p_goal_tag_kernel_distFromEndM_faceted_TRUNC.pdf",
)
plt.show()
"""

# --------------------------------------------------------------------------
# assemble (base order; dropped cells vanish, slots may expand)
CELLS = [
    keep_md(0),
    code(SETUP),
    keep_md(2),
    keep_code(3),
    keep_code(4),
    keep_md(5),
    code(C06),
    keep_md(7),
    code(C08),
    keep_code(9),
    keep_md(10),
    code(C11),
    keep_code(12),
    keep_md(13),
    keep_code(14),
    keep_md(15),
    keep_code(16),
    keep_code(17),
    keep_md(18),
    # c19 dropped: trajectory-panel helpers live in goldshtein_bats_plots
    keep_md(20),
    code(F21),
    keep_md(22),
    code(F23),
    keep_md(24),
    # c25 dropped: superseded by plots.plot_destination_panels_compact
    keep_md(26),
    # c27 dropped: superseded by plots.plot_destination_panels
    keep_md(28),
    code(F29),
    keep_md(30),
    code(F31),
    code(C32),
    keep_md(33),
    code(F34),
    keep_md(35),
    code(F36),
    code(C37),
    keep_md(38),
    code(F39),
    keep_md(40),
    code(C41),
    keep_md(42),
    code(C43),
    keep_md(44),
    keep_code(45),
    keep_code(46),
    keep_code(47),
    keep_code(48),
    keep_md(49),
    code(F50),
    keep_md(51),
    code(F52),
    keep_md(53),
    code(C54),
    C55,
    keep_md(56),
    keep_code(57),
    keep_md(58),
    code(C59A),
    code(F59),
    keep_md(60),
    code(C61A),
    code(F61),
    # c62 dropped: superseded by the posthoc PSIS-LOO machinery above
    keep_md(63),
    code(F64),
    keep_md(65),
    code(C66A),
    code(F66),
    keep_md(67),
    # c68 dropped: superseded by plots.plot_mixture_polar_pdfs
    code(F69),
    keep_md(70),
    code(C71A),
    code(F71_1),
    code(F71_2),
    code(F71_3),
    code(F71_4),
    code(F71_5),
    code(F71_6),
    keep_md(72),
    code(F73),
    keep_md(74),
    code(F75),
    code(F76),
    keep_md(77),
    # c78 dropped: tag-level summaries live in goldshtein_bats_posthoc
    code(C79A),
    code(F79_1),
    code(F79_2),
    code(F79_3),
    code(F79_4),
    code(C80A),
    code(F80),
    keep_md(81),
    code(C82),
    keep_md(83),
    code(F84),
    code(C85),
    keep_md(86),
    code(F87),
    # c88 dropped: superseded by plots.plot_tracks_pgoal_small_bigfont
    code(F89),
    keep_md(90),
    code(C91A),
    code(F91_1),
    code(F91_2),
    # c92 dropped: superseded by plots.plot_population_curve_smooth_compact
    code(F93_1),
    code(F93_2),
    keep_md(94),
    code(C95A),
    code(F95_1),
    code(F95_2),
]

# --------------------------------------------------------------------------
# model cells keep top-level defs by design (the vM Bambi and vMU PyMC
# specifications are model code, not plot helpers)
MODEL_SOURCES = {C43, "".join(C55["source"])}

# --------------------------------------------------------------------------
# audit
FIG_STEMS = sorted({
    "rediscretized_destination_panels",
    "rediscretized_destination_panels_compact",
    "destination_flow_fields_smoothed_distance",
    "polar_vonmises_kde_all",
    "kappa_hat_by_treatment_median_percentiles",
    "rbar_by_treatment_median_percentiles",
    "alignment0_by_treatment_median_percentiles",
    "treatment_coefficients_mu_kappa",
    "treatment_pred_density_at0_mu_kappa",
    "psisloo_total_elpd_all_models",
    "psisloo_delta_total_and_by_treatment_vMU_minus_vM",
    "psisloo_combined",
    "psisloo_delta_by_individual_by_treatment_vMU_minus_vM",
    "vmum_mixture_full_polar_pdfs_population_only",
    "fullmodel_mixture_weight_by_treatment_mean_hdi",
    "fullmodel_kappa_by_treatment_mean_hdi",
    "fullmodel_predictive_density_at_0_mixture_by_treatment_mean_hdi",
    "fullmodel_predictive_density_at_0_directed_component_by_treatment_mean_hdi",
    "fullmodel_predictive_mass_pm45deg_by_treatment_mean_hdi",
    "fullmodel_alignment_to_zero_by_treatment_mean_hdi",
    "psisloo_plus_params_2x2",
    "predictive_density_two_model_comparison_matched_vmu_directed_vs_vm",
    "psisloo_plus_predictive_density_3panel",
    "fullmodel_taglevel_mixture_weight_by_treatment_mean_hdi",
    "fullmodel_taglevel_kappa_by_treatment_mean_hdi",
    "fullmodel_taglevel_predictive_density_at_0_by_treatment_mean_hdi",
    "fullmodel_taglevel_alignment_to_zero_by_treatment_mean_hdi",
    "bat_population_ppc_by_treatment",
    "bat_individual_ppc_bootstrap_vs_ppc_kde",
    "treatment_panels_tracks_pgoal",
    "treatment_panels_tracks_pgoal_small_bigfont",
    "p_goal_pop_kernel_distance_TRUNC_SMOOTH",
    "p_goal_pop_kernel_distFromEndM_TRUNC_SMOOTH",
    "p_goal_pop_kernel_distance_TRUNC_SMOOTH_compact",
    "p_goal_pop_kernel_distFromEndM_TRUNC_SMOOTH_compact",
    "p_goal_tag_kernel_distance_faceted_TRUNC",
    "p_goal_tag_kernel_distFromEndM_faceted_TRUNC",
})

n_code = n_md = 0
for i, c in enumerate(CELLS):
    if c["cell_type"] != "code":
        n_md += 1
        continue
    n_code += 1
    text = "".join(c["source"])
    if text in MODEL_SOURCES:
        continue
    for bad in (
        "def ",
        "sns.",
        "contourf",
        "plt.savefig",
        "get_cmap",
        "out_pdf",
        "out_svg",
        "savepath",
        "save_svg_also",
        "az.hdi",
        "pm.sample",
    ):
        assert bad not in text, f"cell {i}: forbidden {bad!r}"
    plt_refs = [l for l in text.splitlines() if "plt." in l]
    if "save_path=" in text:
        # a figure cell: comment + one plots.* call + plt.show()
        assert plt_refs and plt_refs[-1].strip() == "plt.show()", f"cell {i}: figure cell must end with plt.show()"
        non_comment = [l for l in text.splitlines() if l.strip() and not l.lstrip().startswith("#")]
        assert any(l.strip() == "plt.show()" for l in non_comment), f"cell {i}: needs plt.show()"
        assert text.count("plots.") >= 1, f"cell {i}: figure cell must call plots.*"
        assert not any("plt." in l and l.strip() != "plt.show()" for l in non_comment), f"cell {i}: stray plt use"
    else:
        assert not plt_refs, f"cell {i}: non-figure code cell must not touch plt"
    # syntax gate (IPython magics like %matplotlib inline are not python)
    if not any(l.lstrip().startswith("%") for l in text.splitlines()):
        compile(text, f"<cell {i}>", "exec")

n_saves = sum(1 for c in CELLS if c["cell_type"] == "code" and "save_path=" in "".join(c["source"]))
assert n_saves == len(FIG_STEMS), f"figure cells {n_saves} != expected {len(FIG_STEMS)}"
for stem in FIG_STEMS:
    hits = [i for i, c in enumerate(CELLS)
            if c["cell_type"] == "code" and f'figures/bats/{stem}.pdf' in "".join(c["source"])]
    assert len(hits) == 1, f"stem {stem}: {len(hits)} save cells (expected 1)"

# markdown must be byte-identical to the base annotation (same order,
# nothing reworded)
n_md_base = sum(1 for c in base["cells"] if c["cell_type"] == "markdown")
assert n_md == n_md_base == 38, f"markdown count {n_md} != base {n_md_base}"
assert all(
    "".join(a["source"]) == "".join(b["source"])
    for a, b in zip((c for c in CELLS if c["cell_type"] == "markdown"),
                    (c for c in base["cells"] if c["cell_type"] == "markdown"))
), "markdown drifted from the base annotation"

# --------------------------------------------------------------------------
out = {**base}  # keeps nbformat, nbformat_minor, metadata
out["cells"] = CELLS
json.dump(out, open(NB_PATH, "w"), indent=1, ensure_ascii=False)
print(f"OK: wrote {NB_PATH} — {len(CELLS)} cells ({n_code} code, {n_md} markdown, {n_saves} figures)")
