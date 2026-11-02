"""Gate for the re-baked periodic_not_circular.ipynb (house layout: 11 cells).

Cell map: c02 setup (version dict incl. pymc/arviz, SEED = 42, sibling
import of the plotting script), c03 set style (echo to
figures/periodic_not_circular), c05 seeded simulation of the 24-hour
beta process, c07 the linear-model + sigma-function spec (top-level model
code, whitelisted), c08 the seeded fit (whitelisted; prints
Divergences: first), c10 the single figure cell (100 dpi inline + 300 dpi
PDF under figures/periodic_not_circular/).

Source audit: no savefig/banned tokens outside the two model cells; the
figure cell is a comment + one plots.plot_*(save_path=...) call +
plt.show() and nothing else; every other code cell never touches plt
(the setup cell's import is exempt).
"""
import base64
import io
import json
from pathlib import Path

from PIL import Image

ROOT = Path("/tmp/opencode/wt-periodic-house")
nb = json.load(open(ROOT / "notebooks/periodic_not_circular.ipynb"))

MODEL_CELLS = {7, 8}
SETUP_CELLS = {2}
FIG_CELLS = {10}
EXPECT_SAVED = {10: "periodic_beta_fit"}

BANNED = ("def ", "sns.", "contourf", "plt.savefig", "savefig", "get_cmap",
          "out_pdf", "out_svg", "savepath", "save_svg_also", "pm.sample")

errs = []
version_line = None
style_line = None
fit_lines = []
div_lines = []

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
                if i == 8 and "NUTS[" in l:
                    fit_lines.append((i, l.strip()))
                if i == 8 and l.startswith("Divergences:"):
                    div_lines.append((i, l.strip()))
        if o.get("output_type") in ("execute_result", "display_data"):
            data = o.get("data", {})
            txt = data.get("text/plain", "")
            if isinstance(txt, list):
                txt = "".join(txt)
            if "pymc" in txt and "arviz" in txt and version_line is None:
                version_line = (i, txt.strip()[:110])
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
print("set_style echo c03:", style_line)
print("c08 fit line:", fit_lines)
print("c08 divergences:", div_lines)

ok = True
if errs:
    ok = False
if source_bad:
    ok = False
if not version_line or version_line[0] != 2:
    ok = False
if not style_line or style_line[0] != 3 or "figures/periodic_not_circular" not in style_line[1]:
    ok = False
if not fit_lines:
    print("BAD: c08 sampler start line (expect 'NUTS[sampler]: [vars]' — progressbar=False suppresses the Sampling line)")
    ok = False
if not any(i == 8 and l.startswith("Divergences: 0") for i, l in div_lines):
    print("BAD divergence print (expect 'Divergences: 0' in c08):", div_lines)
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

fig_dir = ROOT / "figures" / "periodic_not_circular"
on_disk_pdf = sorted(p.name for p in fig_dir.glob("*.pdf")) if fig_dir.exists() else []
on_disk_svg = sorted(p.name for p in fig_dir.glob("*.svg")) if fig_dir.exists() else []
for i, stem in EXPECT_SAVED.items():
    lines = saved_files.get(i, [])
    if not lines:
        print(f"MISSING 'Saved figure to' line in c{i}")
        ok = False
    matched_files = []
    for rel in lines:
        if not rel.startswith("figures/periodic_not_circular/"):
            print(f"BAD save path {rel} in c{i}")
            ok = False
            continue
        name = Path(rel).name
        if name.rsplit(".", 1)[0] != stem or not name.endswith(".pdf"):
            print(f"BAD saved name in c{i}: {name} (expected {stem}.pdf)")
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

expected_on_disk = set()
for lines in saved_files.values():
    for rel in lines:
        if rel.startswith("figures/periodic_not_circular/"):
            expected_on_disk.add(Path(rel).name)
stray = sorted(set(on_disk_pdf) - expected_on_disk)
if stray:
    print("STRAY PDFs on disk not claimed by any cell:", stray)
    ok = False
if on_disk_svg:
    print("STRAY SVGs on disk (house convention is PDF):", on_disk_svg)
    ok = False

print("PERIODIC GATE:", "PASS" if ok else "FAIL")
