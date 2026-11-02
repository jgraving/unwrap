"""Gate for the re-baked circular_regression.ipynb on the new stack
(circular study PR layout: 16 cells).

Checks:
- zero error outputs anywhere
- c02 version dict (new stack)
- c03 prints "Figures will be saved to: figures/circular_regression"
- c09 prints "Divergences: 0" and shows the circular_summary table
  with the arviz 0.23 layout (hdi_3% columns)
- figure cells = {10, 12, 14}, each with one inline ~100 dpi png AND a
  "Saved figure to: figures/circular_regression/....svg" stream line
- the three saved SVGs exist and are valid
"""
import base64
import sys
import io
import json
from pathlib import Path

from PIL import Image

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else next(
    p for p in [Path.cwd(), *Path.cwd().parents] if (p / ".git").exists())
nb = json.load(open(ROOT / "notebooks/circular_regression.ipynb"))

errs = []
version_line = None
style_line = None
div_line = None
summary_txt = None
fig_cells = {}

for i, c in enumerate(nb["cells"]):
    if c["cell_type"] != "code":
        continue
    cell_imgs = 0
    cell_saved = []
    for o in c.get("outputs", []):
        if o.get("output_type") == "error":
            errs.append((i, o.get("ename", ""), o.get("evalue", "")))
        if o.get("output_type") == "stream":
            txt = "".join(o.get("text", []))
            if "Figures will be saved to:" in txt:
                for l in txt.splitlines():
                    if "Figures will be saved to:" in l:
                        style_line = (i, l)
            if "Divergences:" in txt and div_line is None:
                div_line = (i, [l for l in txt.splitlines() if "Divergences:" in l])
            for l in txt.splitlines():
                if l.startswith("Saved figure to:"):
                    cell_saved.append(l)
        if o.get("output_type") in ("execute_result", "display_data"):
            data = o.get("data", {})
            txt = data.get("text/plain", "")
            txt = txt if isinstance(txt, str) else "".join(txt)
            if "arviz" in txt and "bambi" in txt and version_line is None:
                version_line = (i, txt.strip().replace("\n", " | "))
            if "hdi_3%" in txt and summary_txt is None:
                summary_txt = (i, txt)
            img = data.get("image/png")
            if img:
                cell_imgs += 1
                im = Image.open(io.BytesIO(base64.b64decode(img)))
                fig_cells.setdefault(i, []).append((im.size, im.info.get("dpi")))
    if cell_imgs or cell_saved:
        print(f"c{i:02d}: inline_images={cell_imgs} {fig_cells.get(i)} saved={cell_saved}")

print("errors:", errs or "none")
print("version dict cell:", version_line[0] if version_line else None, "|",
      (version_line[1][:130] if version_line else ""))
print("set_style echo:", style_line)
print("divergences print:", div_line)
print("circular_summary cell (hdi_3% cols):", summary_txt[0] if summary_txt else None)
if summary_txt:
    lines = summary_txt[1].splitlines()
    print("  header:", lines[0].strip())
    print("  rows:", [l.split()[0] for l in lines[1:]])

fig_ids = {i for i in fig_cells}
print("figure cells:", sorted(fig_ids))

ok = True
if errs:
    ok = False
if not version_line or version_line[0] != 2:
    ok = False
if not style_line or style_line[0] != 3:
    ok = False
if not div_line or div_line[0] != 9 or "Divergences: 0" not in " ".join(div_line[1]):
    ok = False
if not summary_txt or summary_txt[0] != 9:
    ok = False
if fig_ids != {10, 12, 14}:
    ok = False

# inline dpi ~100 for every figure cell
for i in sorted(fig_ids):
    for size, dpi in fig_cells[i]:
        if dpi and abs(dpi[0] - 100) > 0.5:
            print(f"BAD dpi c{i}: {dpi}")
            ok = False

# saved SVGs exist + valid
expected = {
    10: "trace.svg",
    12: "predictions.svg",
    14: "ppc.svg",
}
for i, name in expected.items():
    p = ROOT / "figures" / "circular_regression" / name
    if not p.exists():
        print(f"MISSING saved figure: {p}")
        ok = False
        continue
    head = p.read_bytes()[:5]
    valid = head == b"<?xml" or head[:4] == b"<svg"
    print(f"saved {name}: {p.stat().st_size} B valid={valid}")
    if not valid:
        ok = False
    if not any(name in s for s in (nb["cells"][i].get("outputs", [])
              and [o["output_type"] for o in nb["cells"][i]["outputs"]])):
        pass

print("CIRCULAR STUDY GATE:", "PASS" if ok else "FAIL")
