#!/usr/bin/env python3
"""Generate Quarto (.qmd) twins from the R Jupyter notebooks.

The `.ipynb` files are the source of truth: they are baked, tested, and
rendered by the documentation site. The `.qmd` twins are generated views
for RStudio/Quarto users, who can run them with plain R and no Jupyter
kernel. Run from the repository root after any R notebook change:

    python tools/ipynb_to_qmd.py

Add `--check` to regenerate to memory and compare against the committed
twins instead of writing (exit 1 on drift); CI uses this to keep the
twins in sync with the notebooks.

Per-cell transformation (the R code is otherwise untouched):
- the `IRdisplay::display_png(file = file)` line inside `show_figure` is
  dropped and the trailing `invisible(NULL)` returns the plot object,
  which Quarto renders natively;
- `"IRdisplay"` entries in the setup cells' package vectors are dropped
  (the package is only needed to push figures over the Jupyter display
  channel);
- the setup prose and a repo-path comment are reworded so the twins
  contain no Jupyter-specific residue;
- code cells become `r` engine chunks carrying the notebook cell id
  (`#| cell-id:`), so re-runs keep stable Quarto cell identifiers.
"""
import difflib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
NB_DIR = ROOT / "notebooks_r"

# show_figure has one shared shape across all notebooks_r notebooks.
OLD_SHOW_TAIL = (
    "  invisible(dev.off())\n"
    "  IRdisplay::display_png(file = file)\n"
    "  invisible(NULL)\n"
)
NEW_SHOW_TAIL = (
    "  invisible(dev.off())\n"
    "  # Quarto renders the returned plot object natively; there is no\n"
    "  # display channel to write to outside Jupyter.\n"
    "  plot\n"
)

# Native-R figures must not depend on Cairo SVG: headless R builds (the kind
# people fit models on with a bare R install) often lack it. ragg is already
# in the pinned environment — the Jupyter notebooks use ragg::agg_png for the
# same reason — so Quarto's figure device is pointed at it.
DEVICE_CHUNK = [
    "```{r}",
    "#| label: figure-device",
    "#| dev: ragg_png",
    "#| echo: false",
    "#| warning: false",
    "# The Jupyter notebook shows inline figures through its display",
    "# channel, which native R has no equivalent of, so Quarto captures",
    "# the returned plot objects itself. The device is ragg (headless-safe)",
    "# rather than the Cairo SVG default, which many headless R builds lack.",
    'knitr::opts_chunk$set(dev = "ragg_png")',
    "```",
    "",
]

# "IRdisplay" in the setup cells' package vectors. Three shapes occur:
# mid-vector on a shared line, the element on its own line (trailing
# comma, line wrap), and the final element before the closing paren.
# The lookbehind keeps a qualified form (e.g. `pkg::IRdisplay`) untouched;
# no code uses that form.
IRDISPLAY_VECTOR = [
    (re.compile(r'(?<![\w.])"IRdisplay",\s*\n'), ""),
    (re.compile(r'(?<![\w.])"IRdisplay",\s*'), ""),
    (re.compile(r'(?<![\w.]),\s*"IRdisplay"\s*(?=\))'), ""),
]

# Setup-prose rewordings (the Setup cell is the second markdown cell in
# every notebook, right after the title; a pattern missing in a file is
# fine, the notebooks vary).
SETUP_MD_REPLACEMENTS = [
    (
        "page once, select the\n**R (`ir`)** kernel, and run all cells in order.",
        "page once, and run all cells in order.",
    ),
    (
        "Select the **R (`ir`)** kernel after following the [R setup](../getting-started-r.md)\npage",
        "Follow the [R setup](../getting-started-r.md)\npage",
    ),
    (
        "Select the R (`ir`) kernel and follow the [R setup](../getting-started-r.md)\n"
        "page for the environment and kernel prerequisites",
        "Follow the [R setup](../getting-started-r.md)\n"
        "page for the environment and kernel prerequisites",
    ),
    (
        "CmdStan, and the Jupyter R kernel\nagainst the pinned `renv.lock` versions.",
        'CmdStan, and the Jupyter R kernel (used by the `.ipynb` twin; this\n'
        "file runs without any kernel)\nagainst the pinned `renv.lock` versions.",
    ),
]

# Code comment that only makes sense in the Jupyter side.
WORKDIR_COMMENT = (
    "make them reachable regardless of the\n# kernel's working directory.",
    "make them reachable regardless of the\n# current working directory.",
)

# Phrases that must never survive in a generated twin.
RESIDUE_PHRASES = [
    "IRdisplay",
    "select the **R (`ir`)** kernel",
    "Select the R (`ir`) kernel",
    "Select the **R (`ir`)** kernel",
    "kernel's working directory",
]


def markdown_title(md_text, fallback):
    m = re.search(r"^#\s+(.+)$", md_text, flags=re.MULTILINE)
    return (m.group(1).strip() if m else fallback).replace('"', "\\'")


def check_residue(text: str, name: str) -> None:
    hits = [p for p in RESIDUE_PHRASES if p in text]
    if hits:
        sys.exit(f"{name}: Jupyter residue in generated twin: {hits}")


def convert(nb_path: Path) -> str:
    nb = json.loads(nb_path.read_text())
    stem = nb_path.stem
    lines = []
    lines += [
        "---",
        'engine: r',
        f'title: "{markdown_title("".join(nb["cells"][0]["source"]), stem)}"',
        "format: html",
        "---",
        "",
        f"> Generated by `tools/ipynb_to_qmd.py` from `{stem}.ipynb`, the source",
        "> of truth (its committed outputs are what the documentation site",
        "> renders). Do not edit this file by hand; re-run the generator.",
        "> Runs in RStudio/Quarto with plain R (no Jupyter kernel required).",
        "",
    ]
    lines += DEVICE_CHUNK
    for cell in nb["cells"]:
        src = "".join(cell["source"])
        if cell.get("cell_type") == "markdown":
            # Every pattern below occurs in at most one cell, and only in
            # the one it targets, so applying them to all markdown cells is
            # equivalent to targeting the Setup cell without finding it.
            for old, new in SETUP_MD_REPLACEMENTS:
                src = src.replace(old, new)
            lines += [src, ""]
        else:
            if OLD_SHOW_TAIL in src:
                if src.count(OLD_SHOW_TAIL) != 1:
                    sys.exit(f"{nb_path.name}: expected 1 show_figure tail, "
                             f"found {src.count(OLD_SHOW_TAIL)}")
                src = src.replace(OLD_SHOW_TAIL, NEW_SHOW_TAIL)
            for rx, rep in IRDISPLAY_VECTOR:
                src = rx.sub(rep, src)
            src = src.replace(*WORKDIR_COMMENT)
            if "IRdisplay" in src or WORKDIR_COMMENT[0] in src:
                sys.exit(f"{nb_path.name}: cell {cell.get('id')}: "
                         f"untransformed Jupyter residue")
            lines += ["```{r}", f'#| cell-id: "{cell.get("id", stem)}"', src, "```", ""]
    text = "\n".join(lines).rstrip("\n") + "\n"
    check_residue(text, nb_path.name)
    return text


def main() -> None:
    check_mode = "--check" in sys.argv
    ipynbs = sorted(NB_DIR.glob("*.ipynb"))
    strays = sorted(p.name for p in NB_DIR.glob("*.qmd")
                    if p.stem not in {ip.stem for ip in ipynbs})
    if strays:
        sys.exit(f"refusing to proceed: qmd files without a matching ipynb: "
                 f"{strays}")
    drift = []
    for nb in ipynbs:
        text = convert(nb)
        out = NB_DIR / (nb.stem + ".qmd")
        if check_mode:
            committed = out.read_text() if out.exists() else None
            if committed != text:
                drift.append((nb.stem, committed, text))
            continue
        out.write_text(text)
        print(f"wrote {out.relative_to(ROOT)}")
    left = sorted(p.name for p in NB_DIR.glob("*.qmd"))
    expected = sorted(p.stem + ".qmd" for p in ipynbs)
    if left != expected:
        sys.exit(f"stray qmd files: {set(left) ^ set(expected)}")
    if check_mode:
        if not drift:
            print(f"all {len(ipynbs)} qmd twins in sync")
            return
        for stem, committed, text in drift:
            print(f"DRIFT {stem}.qmd", file=sys.stderr)
            base = " " * len(stem)
            diff = difflib.unified_diff(
                (committed or "").splitlines(keepends=True),
                text.splitlines(keepends=True),
                fromfile=f"a/{stem}.qmd", tofile=f"b/{stem}.qmd")
            sys.stderr.writelines(diff)
        sys.exit(f"{len(drift)} qmd twin(s) out of sync with the notebooks; "
                 f"run 'python tools/ipynb_to_qmd.py' to regenerate")
    print(f"{len(ipynbs)} qmd twins written")


if __name__ == "__main__":
    main()
