"""Gate for the re-baked locust_density_order.ipynb (study PR layout:
43 cells).

Cell map: c01 setup (version dict + style echo), c09 nunique (440),
c11 trajectory grid (Sorted orders print), c17 screening fit
(Divergences: 0), c21 five-model sweep (5 Divergences prints + Best model:
Density × Order), c22 sweep table print, c25 fit table (num_samples: 18641),
c26 final fit (response params from PyMC graph OK), c27 LOO
(p_loo / N / p_loo), c28 fraction k>0.7, c31/c32 diagnostics filters
(expected empty, as in the original bake), figure cells: c11-c14, c19, c23,
c35-c42 (14 figures, each 100 dpi inline + 300 dpi PDF under figures/locust/).
"""
import base64
import io
import json
from pathlib import Path

from PIL import Image

ROOT = Path("/tmp/opencode/wt-unified")
nb = json.load(open(ROOT / "notebooks/locust_density_order.ipynb"))

FIG_CELLS = {11, 12, 13, 14, 19, 23, 35, 36, 37, 38, 39, 40, 41, 42}
EXPECT_SAVED = {
    11: "trajectories_order_by_density_grid.pdf",
    12: "mean_vectors_order_by_density_grid.pdf",
    13: "order_density_trajectories.pdf",
    14: "order_density_pairwise_trajectories_only.pdf",
    19: "order_density_alignment_screening.pdf",
    23: "order_density_model_selection.pdf",
    35: "order_density_alignment.pdf",
    36: "order_density_elpd0.pdf",
    37: "order_density_polar_pdfs_grid.pdf",
    38: "order_density_coefficients.pdf",
    39: "order_density_model_selection_and_coefficients.pdf",
    40: "order_density_alignment_ppc_bootstrap.pdf",
    41: "order_density_cell_ppc_dots.pdf",
    42: "order_density_individual_cell_ppc.pdf",
}

errs = []
version_line = None
style_line = None
sorted_orders = None
nunique = None
div_lines = []
fitting_lines = []
best_model = None
sweep_print = None
num_samples = None
graph_ok = None
loo_lines = []
frac_k = None
summary_txt = None
diag_ess = None
diag_rhat = None
fig_cells = {}
saved_lines = {}

for i, c in enumerate(nb["cells"]):
    if c["cell_type"] != "code":
        continue
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
                    saved_lines.setdefault(i, []).append(l)
                if l.startswith("Divergences:"):
                    div_lines.append((i, l))
                if l.startswith("=== Fitting "):
                    fitting_lines.append((i, l))
                if l.startswith("Sorted orders:"):
                    sorted_orders = (i, l)
                if l.startswith("Best model:"):
                    best_model = (i, l)
                if l.startswith("num_samples:"):
                    num_samples = (i, l)
                if "response params from PyMC graph OK" in l:
                    graph_ok = (i, l)
                if l.startswith("p_loo (effective params):"):
                    loo_lines.append((i, l))
                if l.startswith("N / p_loo:"):
                    loo_lines.append((i, l))
                if l.startswith("fraction k>0.7:"):
                    frac_k = (i, l)
                if "pretty_name" in l and i == 22:
                    sweep_print = (i, l.strip()[:60]) if sweep_print is None else sweep_print
        if o.get("output_type") in ("execute_result", "display_data"):
            data = o.get("data", {})
            txt = ""
            for key in ("text/plain", "text/html"):
                v = data.get(key, "")
                txt += (" " if txt else "") + (v if isinstance(v, str) else "".join(v))
            if "pymc" in txt and version_line is None:
                version_line = (i, txt.strip().replace("\n", " | "))
            if i == 9 and txt.strip() == "440" and nunique is None:
                nunique = i
            if i == 22 and sweep_print is None and "pretty_name" in txt:
                sweep_print = (i, txt.strip().splitlines()[0][:80])
            if i == 30 and summary_txt is None:
                summary_txt = (i, "has r_hat" if "r_hat" in txt else "NO r_hat",
                               len(txt.strip().splitlines()))
            if i == 31 and diag_ess is None:
                diag_ess = (i, "Empty" if txt.strip().startswith("Empty DataFrame") else "non-empty")
            if i == 32 and diag_rhat is None:
                diag_rhat = (i, "Empty" if txt.strip().startswith("Empty DataFrame") else "non-empty")
            img = data.get("image/png")
            if img:
                im = Image.open(io.BytesIO(base64.b64decode(img)))
                cell_imgs.append((im.size, im.info.get("dpi")))
    if cell_imgs:
        fig_cells[i] = cell_imgs
        print(f"c{i:02d}: inline_images={len(cell_imgs)} {cell_imgs} saved={saved_lines.get(i)}")

print("errors:", errs or "none")
print("version dict cell:", version_line[0] if version_line else None, "|",
      (version_line[1][:110] if version_line else ""))
print("set_style echo:", style_line)
print("Sorted orders:", sorted_orders)
print("nunique cell:", nunique, "(expect 9 -> 440)")
print("divergences prints:", div_lines)
print("=== Fitting lines:", len(fitting_lines), fitting_lines[:2], "...", fitting_lines[-1:] if len(fitting_lines) > 2 else "")
print("best model:", best_model)
print("sweep table print:", sweep_print)
print("num_samples:", num_samples)
print("graph OK:", graph_ok)
print("loo lines:", loo_lines)
print("fraction k>0.7:", frac_k)
print("summary (c30):", summary_txt)
print("diagnostics ess (c31):", diag_ess, "(expect Empty)")
print("diagnostics rhat (c32):", diag_rhat, "(expect Empty)")

ok = True
if errs:
    ok = False
if not version_line or version_line[0] != 1:
    ok = False
if not style_line or style_line[0] != 1:
    ok = False
if not style_line or "figures/locust" not in style_line[1]:
    ok = False
if nunique != 9:
    ok = False
if sorted_orders is None or sorted_orders[0] != 11:
    ok = False
# divergences: 1 screening (c17) + 5 sweep (c21), all zero
screenings = [(i, l) for i, l in div_lines if i == 17]
sweeps = [(i, l) for i, l in div_lines if i == 21]
if len(screenings) != 1 or "Divergences: 0" not in screenings[0][1]:
    ok = False
if len(sweeps) != 5 or any("Divergences: 0" not in l for _, l in sweeps):
    ok = False
if len(fitting_lines) != 5:
    ok = False
if not best_model or best_model[0] != 21 or "Density × Order" not in best_model[1]:
    ok = False
if not sweep_print or sweep_print[0] != 22:
    ok = False
if not num_samples or num_samples[0] != 25 or "18641" not in num_samples[1]:
    ok = False
if not graph_ok or graph_ok[0] != 26:
    ok = False
loo_cells = sorted({i for i, _ in loo_lines})
if loo_cells != [27] or len(loo_lines) != 2:
    ok = False
if not frac_k or frac_k[0] != 28:
    ok = False
if not summary_txt or summary_txt[1] != "has r_hat" or summary_txt[2] < 10:
    ok = False
if diag_ess is None or diag_ess[1] != "Empty":
    ok = False
if diag_rhat is None or diag_rhat[1] != "Empty":
    ok = False

if set(fig_cells) != FIG_CELLS:
    print("BAD figure-cell set:", sorted(fig_cells))
    ok = False
for i, imgs in sorted(fig_cells.items()):
    if len(imgs) != 1:
        print(f"BAD: c{i} has {len(imgs)} images")
        ok = False
    for size, dpi in imgs:
        if dpi and abs(dpi[0] - 100) > 0.5:
            print(f"BAD dpi c{i}: {dpi} size={size}")
            ok = False

for i, name in EXPECT_SAVED.items():
    if not any(name in s for s in saved_lines.get(i, [])):
        print(f"MISSING 'Saved figure to' for {name} in c{i}")
        ok = False
    p = ROOT / "figures" / "locust" / name
    if not p.exists():
        print("MISSING file:", p)
        ok = False
        continue
    head = p.read_bytes()[:5]
    valid = head == b"%PDF-"
    print(f"saved {name}: {p.stat().st_size} B valid={valid}")
    if not valid or p.stat().st_size < 1000:
        ok = False

print("LOCUST DENSITY GATE:", "PASS" if ok else "FAIL")
