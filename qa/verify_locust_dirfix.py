"""Gate for the re-baked locust_direction_mixture_model.ipynb (study PR
layout: 33 cells).

Cell map: c01 setup (version dict + style echo), c02 data, c05 nunique
(52), c10 bimodal fit / c14 unimodal fit (Divergences prints),
c11 bimodal diagnostics filter (expected empty), c15 unimodal diagnostics
filter (expected NON-empty - low-ESS kappa_log rows, as in the original
bake), figure cells: c07, c16, c17, c20, c21, c23, c26-c31 (12 figures,
each 100 dpi inline + 300 dpi PDF under figures/locust/).
"""
import base64
import sys
import io
import json
from pathlib import Path

from PIL import Image

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else next(
    p for p in [Path.cwd(), *Path.cwd().parents] if (p / ".git").exists())
nb = json.load(open(ROOT / "notebooks/locust_direction_mixture_model.ipynb"))

FIG_CELLS = {7, 16, 17, 20, 21, 23, 26, 27, 28, 29, 30, 31}
EXPECT_SAVED = {
    7: "locust_directional_trajectories.pdf",
    16: "locust_bimodal_unimodal_elpd.pdf",
    17: "locust_elpd_by_individual_ab.pdf",
    20: "population_directional_validation_cartesian.pdf",
    21: "population_directional_validation.pdf",
    23: "locust_directional_individual_ppc.pdf",
    26: "locust_directional_hypothesis_horizontal.pdf",
    27: "locust_directional_hypothesis_horizontal_elpd.pdf",
    28: "locust_directional_hypothesis_individual_ab.pdf",
    29: "locust_directional_hypothesis_individual_elpd_ab.pdf",
    30: "locust_directional_alignment_horizontal.pdf",
    31: "locust_directional_alignment_individual_ab.pdf",
}

errs = []
version_line = None
style_line = None
data_lines = []
nunique = None
div_lines = []
diag_bi = None
diag_uni = None
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
        if o.get("output_type") in ("execute_result", "display_data"):
            data = o.get("data", {})
            txt = ""
            for key in ("text/plain", "text/html"):
                v = data.get(key, "")
                txt += (" " if txt else "") + (v if isinstance(v, str) else "".join(v))
            if "pymc" in txt and version_line is None:
                version_line = (i, txt.strip().replace("\n", " | "))
            if i == 2 and "animal_id" in txt and data_lines == []:
                data_lines = [(i, txt.strip().splitlines()[0][:60])]
            if i == 5 and txt.strip() == "52" and nunique is None:
                nunique = i
            if i == 11 and diag_bi is None:
                diag_bi = (i, "Empty" if txt.strip().startswith("Empty DataFrame") else "non-empty",
                           len(txt.strip().splitlines()))
            if i == 15 and diag_uni is None:
                diag_uni = (i, "Empty" if txt.strip().startswith("Empty DataFrame") else "non-empty",
                            len(txt.strip().splitlines()))
            img = data.get("image/png")
            if img:
                im = Image.open(io.BytesIO(base64.b64decode(img)))
                cell_imgs.append((im.size, im.info.get("dpi")))
    if cell_imgs:
        fig_cells[i] = cell_imgs
        print(f"c{i:02d}: inline_images={len(cell_imgs)} {cell_imgs} saved={saved_lines.get(i)}")

print("errors:", errs or "none")
print("version dict cell:", version_line[0] if version_line else None, "|",
      (version_line[1][:130] if version_line else ""))
print("set_style echo:", style_line)
print("data head shown:", data_lines[:1])
print("nunique cell:", nunique, "(expect 5 -> 52)")
print("divergences prints:", div_lines)
print("diagnostics bimodal (c11):", diag_bi, "(expect Empty)")
print("diagnostics unimodal (c15):", diag_uni, "(expect non-empty, as in the original bake)")

ok = True
if errs:
    ok = False
if not version_line or version_line[0] != 1:
    ok = False
if not style_line or style_line[0] != 1:
    ok = False
if "figures/locust" not in (style_line and style_line[1] or ""):
    ok = False
if not data_lines or data_lines[0][0] != 2:
    ok = False
if nunique != 5:
    ok = False
if [i for i, _ in div_lines] != [10, 14] or any("Divergences: 0" not in l for _, l in div_lines):
    ok = False
if diag_bi is None or diag_bi[1] != "Empty":
    ok = False
if diag_uni is None or diag_uni[1] != "non-empty":
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

print("LOCUST GATE:", "PASS" if ok else "FAIL")
