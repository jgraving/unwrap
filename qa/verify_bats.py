"""Gate for the re-baked goldshtein_bats.ipynb (study PR layout: 107 cells).

Cell map: c01 setup (version dict incl. pymc/bambi + set_style echo to
figures/bats), c06 rediscretize (tags: 67), c38 modeling table, c40 vM
Bambi model (top-level model code, whitelisted), c42 az.summary, c43/c44
diagnostics filters (expected Empty), c45 full summary, c51 sigmoid-space
percentiles, c52 vMU PyMC model (whitelisted), c54 convergence flags
(RHAT/ESS none), c56 PSIS-LOO (Best model: von Mises-Uniform), figure
cells: 37 (100 dpi inline + 300 dpi PDF under figures/bats/; the two
faceted p_goal cells may save per-page _p{i} stems, unreachable with the
4-treatment data), both model fits print Divergences: 0.

Source audit (mirrors the builder): no ``def ``/``sns.``/``contourf``/
``savefig``/residue tokens outside the two model cells; figure cells are
comment + one plots.plot_*(save_path=...) + plt.show() and no other plt
use; every other code cell never touches plt.
"""
import base64
import io
import json
import re
from pathlib import Path

from PIL import Image

ROOT = Path("/tmp/opencode/wt-unified")
nb = json.load(open(ROOT / "notebooks/goldshtein_bats.ipynb"))

MODEL_CELLS = {40, 52}
FIG_CELLS = {20, 22, 26, 28, 31, 33, 36, 47, 49, 57, 60, 62, 65, 67,
             70, 71, 72, 73, 74, 75, 77, 79, 80, 83, 84, 85, 86, 88, 92,
             95, 96, 99, 100, 101, 102, 105, 106}
# expected saved stem (no .pdf); per-page variants <stem>_p{i}.pdf are
# accepted for every cell
EXPECT_SAVED = {
    20: "rediscretized_destination_panels",
    22: "rediscretized_destination_panels_compact",
    26: "destination_flow_fields_smoothed_distance",
    28: "polar_vonmises_kde_all",
    31: "kappa_hat_by_treatment_median_percentiles",
    33: "rbar_by_treatment_median_percentiles",
    36: "alignment0_by_treatment_median_percentiles",
    47: "treatment_coefficients_mu_kappa",
    49: "treatment_pred_density_at0_mu_kappa",
    57: "psisloo_total_elpd_all_models",
    60: "psisloo_delta_total_and_by_treatment_vMU_minus_vM",
    62: "psisloo_combined",
    65: "psisloo_delta_by_individual_by_treatment_vMU_minus_vM",
    67: "vmum_mixture_full_polar_pdfs_population_only",
    70: "fullmodel_mixture_weight_by_treatment_mean_hdi",
    71: "fullmodel_kappa_by_treatment_mean_hdi",
    72: "fullmodel_predictive_density_at_0_mixture_by_treatment_mean_hdi",
    73: "fullmodel_predictive_density_at_0_directed_component_by_treatment_mean_hdi",
    74: "fullmodel_predictive_mass_pm45deg_by_treatment_mean_hdi",
    75: "fullmodel_alignment_to_zero_by_treatment_mean_hdi",
    77: "psisloo_plus_params_2x2",
    79: "predictive_density_two_model_comparison_matched_vmu_directed_vs_vm",
    80: "psisloo_plus_predictive_density_3panel",
    83: "fullmodel_taglevel_mixture_weight_by_treatment_mean_hdi",
    84: "fullmodel_taglevel_kappa_by_treatment_mean_hdi",
    85: "fullmodel_taglevel_predictive_density_at_0_by_treatment_mean_hdi",
    86: "fullmodel_taglevel_alignment_to_zero_by_treatment_mean_hdi",
    88: "bat_population_ppc_by_treatment",
    92: "bat_individual_ppc_bootstrap_vs_ppc_kde",
    95: "treatment_panels_tracks_pgoal",
    96: "treatment_panels_tracks_pgoal_small_bigfont",
    99: "p_goal_pop_kernel_distance_TRUNC_SMOOTH",
    100: "p_goal_pop_kernel_distFromEndM_TRUNC_SMOOTH",
    101: "p_goal_pop_kernel_distance_TRUNC_SMOOTH_compact",
    102: "p_goal_pop_kernel_distFromEndM_TRUNC_SMOOTH_compact",
    105: "p_goal_tag_kernel_distance_faceted_TRUNC",
    106: "p_goal_tag_kernel_distFromEndM_faceted_TRUNC",
}
PER_PAGE_OK = {105, 106}

BANNED = ("def ", "sns.", "contourf", "plt.savefig", "savefig", "get_cmap",
          "out_pdf", "out_svg", "savepath", "save_svg_also", "az.hdi", "pm.sample")

errs = []
version_line = None
style_line = None
tags67 = None
best_model = None
rhok = None
essk = None
empty43 = None
empty44 = None
summary45 = None
div_lines = []
pgoal_cols = {}

fig_cells = {}
saved_files = {}   # cell -> [repo-relative saved paths]
source_bad = []

for i, c in enumerate(nb["cells"]):
    if c["cell_type"] != "code":
        continue
    text = "".join(c.get("source", []))

    # ---- source audit -------------------------------------------------
    if i not in MODEL_CELLS:
        for bad in BANNED:
            if bad in text:
                source_bad.append((i, bad))
        plt_refs = [l for l in text.splitlines() if "plt." in l]
        if "save_path=" in text:
            non_comment = [l for l in text.splitlines()
                           if l.strip() and not l.lstrip().startswith("#")]
            if not (plt_refs and plt_refs[-1].strip() == "plt.show()"):
                source_bad.append((i, "figure cell must end with plt.show()"))
            if not any(l.strip() == "plt.show()" for l in non_comment):
                source_bad.append((i, "missing bare plt.show() line"))
            if text.count("plots.") < 1:
                source_bad.append((i, "figure cell must call plots.*"))
            if any("plt." in l and l.strip() != "plt.show()" for l in non_comment):
                source_bad.append((i, "stray plt use in figure cell"))
        else:
            if i != 1 and plt_refs:
                source_bad.append((i, "non-figure code cell touches plt"))

    # ---- outputs ------------------------------------------------------
    cell_imgs = []
    for o in c.get("outputs", []):
        if o.get("output_type") == "error":
            errs.append((i, o.get("ename", ""), o.get("evalue", "")))
        if o.get("output_type") == "stream":
            txt = "".join(o.get("text", []))
            for l in txt.splitlines():
                if "Figures will be saved to:" in l:
                    style_line = (i, l)
                if l.startswith("Saved figure to:"):
                    saved_files.setdefault(i, []).append(l.split("Saved figure to:")[1].strip())
                if l.startswith("Best model:"):
                    best_model = (i, l)
                if "RHAT > 1.01:" in l:
                    rhok = (i, l)
                if "ESS bulk or tail < 400:" in l:
                    essk = (i, l)
                if i == 6 and "tags: 67" in l and tags67 is None:
                    tags67 = i
                if i in (93, 98) and "p_goal_post_mean" in l and i not in pgoal_cols:
                    pgoal_cols[i] = True
                if i in (40, 52) and l.startswith("Divergences:"):
                    div_lines.append((i, l))
        if o.get("output_type") in ("execute_result", "display_data"):
            data = o.get("data", {})
            txt = ""
            for key in ("text/plain", "text/html"):
                v = data.get(key, "")
                txt += (" " if txt else "") + (v if isinstance(v, str) else "".join(v))
            if "pymc" in txt and "bambi" in txt and version_line is None:
                version_line = (i, txt.strip()[:110])
            if i == 43 and empty43 is None:
                empty43 = "Empty" if txt.strip().startswith("Empty DataFrame") else txt.strip()[:40]
            if i == 44 and empty44 is None:
                empty44 = "Empty" if txt.strip().startswith("Empty DataFrame") else txt.strip()[:40]
            if i == 45 and summary45 is None:
                summary45 = ("has r_hat" if "r_hat" in txt else "NO r_hat",
                             len(txt.strip().splitlines()))
            img = data.get("image/png")
            if img:
                im = Image.open(io.BytesIO(base64.b64decode(img)))
                cell_imgs.append((im.size, im.info.get("dpi")))
    if cell_imgs:
        fig_cells[i] = cell_imgs
        print(f"c{i:02d}: inline={len(cell_imgs)} saved={saved_files.get(i)}")

print("errors:", errs or "none")
print("source audit:", source_bad or "clean")
print("version dict c01:", version_line)
print("set_style echo:", style_line)
print("Best model:", best_model)
print("RHAT line (c54):", rhok)
print("ESS line (c54):", essk)
print("c43 empty:", empty43, "| c44 empty:", empty44)
print("c45 summary:", summary45)
print("p_goal columns head (c93/c98):", sorted(pgoal_cols))
print("divergences (c40/c52):", div_lines)

ok = True
if errs:
    ok = False
if source_bad:
    ok = False
if not version_line or version_line[0] != 1:
    ok = False
if not style_line or style_line[0] != 1 or "figures/bats" not in style_line[1]:
    ok = False
if tags67 != 6:
    ok = False
if not best_model or best_model[0] != 56 or "von Mises-Uniform" not in best_model[1]:
    ok = False
if not rhok or rhok[0] != 54 or "none" not in rhok[1]:
    ok = False
if not essk or essk[0] != 54 or "none" not in essk[1]:
    ok = False
if empty43 != "Empty" or empty44 != "Empty":
    ok = False
if not summary45 or summary45[0] != "has r_hat" or summary45[1] < 10:
    ok = False
if sorted(pgoal_cols) != [93, 98]:
    ok = False
# both model fits must report zero divergences (house diagnostic)
div_map = {i: (1 if "Divergences: 0" in l else 0) for i, l in div_lines}
if div_map.get(40, 0) != 1 or div_map.get(52, 0) != 1:
    print("BAD divergence prints (expect 'Divergences: 0' in c40 and c52):", div_lines)
    ok = False

if set(fig_cells) != FIG_CELLS:
    print("BAD figure-cell set:", sorted(fig_cells))
    ok = False
for i, imgs in sorted(fig_cells.items()):
    if len(imgs) != 1:
        print(f"BAD: c{i} has {len(imgs)} inline images")
        ok = False
    for size, dpi in imgs:
        if dpi and abs(dpi[0] - 100) > 0.5:
            print(f"BAD dpi c{i}: {dpi} size={size}")
            ok = False

fig_dir = ROOT / "figures" / "bats"
on_disk = sorted(p.name for p in fig_dir.glob("*.pdf")) if fig_dir.exists() else []
for i, stem in EXPECT_SAVED.items():
    lines = saved_files.get(i, [])
    if not lines:
        print(f"MISSING 'Saved figure to' line in c{i}")
        ok = False
    matched_files = []
    for rel in lines:
        if not rel.startswith("figures/bats/"):
            print(f"BAD save path {rel} in c{i} (must be repo-relative under figures/bats/)")
            ok = False
            continue
        name = Path(rel).name           # <stem>.pdf or <stem>_p{i}.pdf
        m = re.match(r"^(.*)\.pdf$", name)
        stem_got = m.group(1)
        per_page = bool(re.match(rf"^{re.escape(stem)}_p\d+$", stem_got))
        if stem_got != stem and not (i in PER_PAGE_OK and per_page):
            print(f"BAD saved stem in c{i}: {name} (expected {stem}[, _p{{i}}])")
            ok = False
        p = ROOT / rel
        if not p.exists():
            print("MISSING file:", p)
            ok = False
            continue
        valid = p.read_bytes()[:5] == b"%PDF-" and p.stat().st_size >= 1000
        print(f"  c{i:02d} saved {name}: {p.stat().st_size} B valid={valid}")
        if not valid:
            ok = False
        matched_files.append(name)
    if not matched_files:
        print(f"NO valid saved PDF for c{i} ({stem})")
        ok = False

# every PDF on disk must be accounted for by a figure cell
expected_on_disk = set()
for i, lines in saved_files.items():
    for rel in lines:
        if rel.startswith("figures/bats/"):
            expected_on_disk.add(Path(rel).name)
stray = sorted(set(on_disk) - expected_on_disk)
if stray:
    print("STAY PDFs on disk not claimed by any cell:", stray)
    ok = False

print("BATS GATE:", "PASS" if ok else "FAIL")
