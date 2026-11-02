"""Gate for the re-baked zebra_escape.ipynb (study PR layout: 17 cells).

Cell map: c02 setup/version, c03 style, c05 data, c07 nunique,
c09 model+fit (Divergences print), c10 circular_summary (silent),
c12 diagnostics filter, c14 plot_summary, c16 plot_ppc.
"""
import base64
import sys
import io
import json
from pathlib import Path

from PIL import Image

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else next(
    p for p in [Path.cwd(), *Path.cwd().parents] if (p / ".git").exists())
nb = json.load(open(ROOT / "notebooks/zebra_escape.ipynb"))

errs = []
version_line = None
style_line = None
data_lines = []
nunique = None
div_line = None
diag_result = None
fig_cells = {}
saved_lines = {}
xcheck_line = None

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
                    div_line = (i, l)
                if "response params from PyMC graph" in l:
                    xcheck_line = (i, l)
        if o.get("output_type") in ("execute_result", "display_data"):
            data = o.get("data", {})
            txt = data.get("text/plain", "")
            txt = txt if isinstance(txt, str) else "".join(txt)
            if "pymc" in txt and version_line is None:
                version_line = (i, txt.strip().replace("\n", " | "))
            if i == 5 and "HeadingRadians" in txt and data_lines == []:
                data_lines = [(i, txt.strip().splitlines()[0][:60])]
            if txt.strip() == "18" and nunique is None:
                nunique = i
            if i == 12 and diag_result is None:
                diag_result = (i, txt.strip()[:120], len(txt.strip().splitlines()))
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
print("nunique cell:", nunique, "(expect 7)")
print("divergences print:", div_line)
print("diagnostics result (c12):", diag_result)
print("cross-check line:", xcheck_line)

ok = True
if errs:
    ok = False
if not version_line or version_line[0] != 2:
    ok = False
if not style_line or style_line[0] != 3:
    ok = False
if not data_lines or data_lines[0][0] != 5:
    ok = False
if nunique != 7:
    ok = False
if not div_line or div_line[0] != 9 or "Divergences: 0" not in div_line[1]:
    ok = False
if diag_result is None or diag_result[0] != 12:
    ok = False
if not xcheck_line or xcheck_line[0] != 16 or "OK" not in xcheck_line[1]:
    ok = False

if set(fig_cells) != {14, 16}:
    ok = False
for i, imgs in sorted(fig_cells.items()):
    if len(imgs) != 1:
        print(f"BAD: c{i} has {len(imgs)} images")
        ok = False
    for size, dpi in imgs:
        if dpi and abs(dpi[0] - 100) > 0.5:
            print(f"BAD dpi c{i}: {dpi} size={size}")
            ok = False

exp = {14: "zebra_main_summary.pdf", 16: "zebra_ppc_bootstrap_vs_posterior.pdf"}
for i, name in exp.items():
    if not any(name in s for s in saved_lines.get(i, [])):
        print(f"MISSING 'Saved figure to' for {name} in c{i}")
        ok = False
    p = ROOT / "figures" / "zebra" / name
    if not p.exists():
        print("MISSING file:", p)
        ok = False
        continue
    head = p.read_bytes()[:5]
    valid = head == b"%PDF-"
    print(f"saved {name}: {p.stat().st_size} B valid={valid}")
    if not valid:
        ok = False

print("ZEBRA GATE:", "PASS" if ok else "FAIL")
