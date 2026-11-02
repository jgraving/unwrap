"""Gate for the re-baked beetles.ipynb (study PR layout: 25 cells).

Cell map: c02 setup (version dict incl. bambi/pymc, SEED = 42), c04 set
style (echo to figures/beetles), c06 E2E data load (load_beetle_model_table
from unwrap-data; the modelling table head), c08/c15/c16/c24 figure cells
(100 dpi inline + 300 dpi PDF under figures/beetles/), c10 the von Mises
Bambi fit (whitelisted model cell; 2000+2000, seeded, prints Divergences
first), c12/c13 circular_summary diagnostics, c19/c21 kappa_mu summaries,
c24 the PPC figure with the PyMC-graph log-likelihood cross-check.

Source audit: no ``def ``/matplotlib-"how" tokens outside the model cell;
figure cells are comment(s) + one plots.plot_*(save_path=...) +
plt.show() and nothing else; every other code cell never touches plt
(the setup cell's import is exempt).
"""
import base64
import io
import json
from pathlib import Path

from PIL import Image

ROOT = Path("/tmp/opencode/wt-beetles-study")
nb = json.load(open(ROOT / "notebooks/beetles.ipynb"))

MODEL_CELLS = {10}
SETUP_CELLS = {2}
FIG_CELLS = {8, 15, 16, 24}
EXPECT_SAVED = {
    8: "data_individual_polar_grid",
    15: "treat_polar_mu_kde_clean",
    16: "treat_concentration_plus_floodlight",
    24: "beetle_posterior_predictive_check",
}

BANNED = ("def ", "sns.", "contourf", "plt.savefig", "savefig", "get_cmap",
          "out_pdf", "out_svg", "savepath", "save_svg_also", "az.hdi", "pm.sample")

errs = []
version_line = None
style_line = None
head_cols = None
fit_lines = []
div_lines = []
summary12 = None
summary13 = None
kappamu19 = None
kappamu21_stream = None
kappamu21_display = None
ppc_check_line = None

fig_cells = {}
saved_files = {}
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
            if i not in SETUP_CELLS and plt_refs:
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
                    style_line = (i, l.strip())
                if l.startswith("Saved figure to:"):
                    saved_files.setdefault(i, []).append(l.split("Saved figure to:")[1].strip())
                if i == 10 and "Sampling 4 chains for 2_000 tune and 2_000 draw iterations" in l:
                    fit_lines.append((i, l.strip()))
                if i == 10 and l.startswith("Divergences:"):
                    div_lines.append((i, l.strip()))
                if i == 19 and "Kappa_mu" in l and kappamu19 is None:
                    kappamu19 = (i, l.strip())
                if i == 21 and "Calculating second-level kappamu" in l:
                    kappamu21_stream = (i, l.strip())
                if i == 24 and "log-likelihood cross-check passed" in l:
                    ppc_check_line = (i, l.strip())
        if o.get("output_type") in ("execute_result", "display_data"):
            data = o.get("data", {})
            txt = ""
            for key in ("text/plain", "text/html"):
                v = data.get(key, "")
                txt += (" " if txt else "") + (v if isinstance(v, str) else "".join(v))
            if "pymc" in txt and "bambi" in txt and version_line is None:
                version_line = (i, txt.strip()[:110])
            if i == 6 and head_cols is None:
                head_cols = ("HeadingRadians" in txt and "Individual" in txt,
                             txt.strip().splitlines()[0][:60] if txt.strip() else "")
            if i == 12 and summary12 is None:
                summary12 = "has r_hat" if "r_hat" in txt else "NO r_hat"
            if i == 13 and summary13 is None:
                summary13 = "has r_hat" if "r_hat" in txt else "NO r_hat"
            if i == 21 and kappamu21_display is None and "kappa_mu" in txt:
                kappamu21_display = ("kappa_mu_intercept" in txt and "kappa_mu_delta" in txt,
                                     txt.strip().splitlines()[0][:60])
            img = data.get("image/png")
            if img:
                im = Image.open(io.BytesIO(base64.b64decode(img)))
                cell_imgs.append((im.size, im.info.get("dpi")))
    if cell_imgs:
        fig_cells[i] = cell_imgs
        print(f"c{i:02d}: inline={len(cell_imgs)} saved={saved_files.get(i)}")

print("errors:", errs or "none")
print("source audit:", source_bad or "clean")
print("version dict c02:", version_line)
print("set_style echo c04:", style_line)
print("c06 table head ok:", head_cols)
print("c10 fit line:", fit_lines)
print("c10 divergences:", div_lines)
print("c12 circular_summary:", summary12, "| c13:", summary13)
print("c19 kappa_mu:", kappamu19)
print("c21 stream:", kappamu21_stream, "| display:", kappamu21_display)
print("c24 cross-check:", ppc_check_line)

ok = True
if errs:
    ok = False
if source_bad:
    ok = False
if not version_line or version_line[0] != 2:
    ok = False
if not style_line or style_line[0] != 4 or "figures/beetles" not in style_line[1]:
    ok = False
if not head_cols or not head_cols[0]:
    print("BAD c06 modelling-table head (expect Light/Individual/Trial/HeadingRadians columns)")
    ok = False
if not fit_lines:
    print("BAD: c10 fit line (expect 'Sampling 4 chains for 2_000 tune and 2_000 draw iterations')")
    ok = False
div_ok = any(i == 10 and l.startswith("Divergences: 0") for i, l in div_lines)
if not div_ok:
    print("BAD divergence print (expect 'Divergences: 0' in c10):", div_lines)
    ok = False
if summary12 != "has r_hat" or summary13 != "has r_hat":
    ok = False
if kappamu19 is None:
    print("BAD: no Kappa_mu line in c19")
    ok = False
if kappamu21_stream is None or not kappamu21_display or not kappamu21_display[0]:
    print("BAD: c21 second-level kappa_mu summary missing")
    ok = False
if ppc_check_line is None:
    print("BAD: c24 PyMC-graph log-likelihood cross-check line missing")
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

fig_dir = ROOT / "figures" / "beetles"
on_disk_pdf = sorted(p.name for p in fig_dir.glob("*.pdf")) if fig_dir.exists() else []
on_disk_svg = sorted(p.name for p in fig_dir.glob("*.svg")) if fig_dir.exists() else []
for i, stem in EXPECT_SAVED.items():
    lines = saved_files.get(i, [])
    if not lines:
        print(f"MISSING 'Saved figure to' line in c{i}")
        ok = False
    matched_files = []
    for rel in lines:
        if not rel.startswith("figures/beetles/"):
            print(f"BAD save path {rel} in c{i} (must be repo-relative under figures/beetles/)")
            ok = False
            continue
        name = Path(rel).name
        stem_got = name.rsplit(".", 1)[0]
        if stem_got != stem:
            print(f"BAD saved name in c{i}: {name} (expected {stem}.pdf)")
            ok = False
        p = ROOT / rel
        if not p.exists():
            print("MISSING file:", p)
            ok = False
            continue
        valid = p.suffix == ".pdf" and p.read_bytes()[:5] == b"%PDF-" and p.stat().st_size >= 1000
        print(f"  c{i:02d} saved {name}: {p.stat().st_size} B valid={valid}")
        if not valid:
            ok = False
        matched_files.append(name)
    if not matched_files:
        print(f"NO valid saved PDF for c{i} ({stem})")
        ok = False

# every figure file on disk must be claimed by exactly one figure cell, and
# no .svg may remain (the house save convention is 300 dpi PDF)
expected_on_disk = set()
for i, lines in saved_files.items():
    for rel in lines:
        if rel.startswith("figures/beetles/"):
            expected_on_disk.add(Path(rel).name)
stray = sorted(set(on_disk_pdf) - expected_on_disk)
if stray:
    print("STRAY PDFs on disk not claimed by any cell:", stray)
    ok = False
if on_disk_svg:
    print("STRAY SVGs on disk (house convention is PDF):", on_disk_svg)
    ok = False

print("BEETLES GATE:", "PASS" if ok else "FAIL")
