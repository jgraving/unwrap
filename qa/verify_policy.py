"""Gate for the re-baked von_mises_policy.ipynb (house layout: 11 cells).

Cell map: c02 setup (version dict, sibling import of the plotting
script), c03 set style (echo to figures/von_mises_policy), c05 policy +
configuration (top-level defs kept inline by design: this is a teaching
notebook, not a fitted study), c07 Monte-Carlo trajectories (pure numpy),
c09 the single figure cell (100 dpi inline + 300 dpi PDF under
figures/von_mises_policy/).

Source audit: the figure cell is a comment + one plots.plot_*(save_path=...)
call + plt.show() and nothing else; no other code cell touches plt.
"""
import base64
import io
import json
from pathlib import Path

from PIL import Image

ROOT = Path("/tmp/opencode/wt-policy")
nb = json.load(open(ROOT / "notebooks/von_mises_policy.ipynb"))

SETUP_CELLS = {2}
FIG_CELLS = {9}
EXPECT_SAVED = {9: "vonmises_navigation"}

errs = []
version_line = None
style_line = None

fig_cells = {}
saved_files = {}
source_bad = []

for i, c in enumerate(nb["cells"]):
    if c["cell_type"] != "code":
        continue
    text = "".join(c.get("source", []))

    # ---- source audit -------------------------------------------------
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
    for bad in ("plt.savefig", "savefig", "out_pdf", "out_svg", "save_svg_also"):
        if bad in text:
            source_bad.append((i, bad))

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
        if o.get("output_type") in ("execute_result", "display_data"):
            data = o.get("data", {})
            txt = data.get("text/plain", "")
            if isinstance(txt, list):
                txt = "".join(txt)
            if "unwrap" in txt and "numpy" in txt and version_line is None:
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

ok = True
if errs:
    ok = False
if source_bad:
    ok = False
if not version_line or version_line[0] != 2:
    ok = False
if not style_line or style_line[0] != 3 or "figures/von_mises_policy" not in style_line[1]:
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

fig_dir = ROOT / "figures" / "von_mises_policy"
on_disk_pdf = sorted(p.name for p in fig_dir.glob("*.pdf")) if fig_dir.exists() else []
on_disk_svg = sorted(p.name for p in fig_dir.glob("*.svg")) if fig_dir.exists() else []
for i, stem in EXPECT_SAVED.items():
    lines = saved_files.get(i, [])
    if not lines:
        print(f"MISSING 'Saved figure to' line in c{i}")
        ok = False
    matched_files = []
    for rel in lines:
        if not rel.startswith("figures/von_mises_policy/"):
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
        if rel.startswith("figures/von_mises_policy/"):
            expected_on_disk.add(Path(rel).name)
stray = sorted(set(on_disk_pdf) - expected_on_disk)
if stray:
    print("STRAY PDFs on disk not claimed by any cell:", stray)
    ok = False
if on_disk_svg:
    print("STRAY SVGs on disk (house convention is PDF):", on_disk_svg)
    ok = False

print("POLICY GATE:", "PASS" if ok else "FAIL")
