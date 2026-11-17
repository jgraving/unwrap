#!/usr/bin/env python3
"""Extract original JPEGs and a versioned notebook's screening result table.

No page screenshots or resampling: image streams retain their native pixels.
PyMuPDF is required only for this optional source-preparation command.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
PAPER_REV = "9bef710eb44fd8d20a61a7119a12e4933d02ce6d"
NOTEBOOK_REV = "075cb7274f8d80f284346dbe30bdc3e0f20dfd89"
PDF_SHA = "48609fdc9185e20227ae0b1c20dd377301d5aa962c7d57078ef1bf5873d0bf21"
NOTEBOOK_SHA = "073d9272f24ce6b648f92550ac25156b269c4a52599b58be00a05bf874e971a7"
IMAGES = {11: "vr_arena", 12: "real_locust", 13: "virtual_locust_lateral",
          14: "virtual_locust_front", 15: "virtual_locust_rear",
          20: "real_swarm", 21: "virtual_swarm"}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def extract_selection(data):
    notebook = json.loads(data)
    matches = []
    for index, cell in enumerate(notebook["cells"]):
        source = "".join(cell["source"])
        if 'print(sweep_df[[' not in source:
            continue
        texts = ["".join(o.get("text", [])) for o in cell.get("outputs", [])]
        rows = []
        for line in "\n".join(texts).splitlines():
            # Names contain spaces; the next four columns are numerical.
            m = re.match(r"^\s*(.+?)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(.+?)\s*$", line)
            if m:
                name, elpd, delta, se, pareto, best = m.groups()
                rows.append(dict(name=name, elpd=float(elpd), delta=float(delta),
                                 se=float(se), pareto_k_max=float(pareto), best=best,
                                 has_order="Order" in name))
        if len(rows) == 5:
            matches.append((index, rows))
    if len(matches) != 1:
        raise ValueError("Expected one five-model screening table; refusing to guess")
    index, rows = matches[0]
    return {"source_repository": "jgraving/unwrap", "source_revision": NOTEBOOK_REV,
            "source_path": "notebooks/locust_density_order.ipynb",
            "source_sha256": digest(data), "cell_index": index,
            "basis": "Saved Python screening output; 440 recording mean directions; printed to six decimals",
            "recompute": "locust_density_posthoc.model_sweep(screen_df, draws=1000, seed=123)",
            "models": rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf", required=True, type=Path)
    parser.add_argument("--notebook", required=True, type=Path,
                        help="Notebook from the pinned source revision")
    parser.add_argument("--output-dir", type=Path,
                        default=ROOT / "notebooks/assets/locust_density_order")
    args = parser.parse_args()
    data = args.pdf.read_bytes()
    if digest(data) != PDF_SHA:
        raise ValueError("PDF differs from the pinned manuscript figure")
    notebook = args.notebook.read_bytes()
    if digest(notebook) != NOTEBOOK_SHA:
        raise ValueError("Notebook differs from the pinned source revision")
    import fitz
    args.output_dir.mkdir(parents=True, exist_ok=True)
    doc = fitz.open(stream=data, filetype="pdf")
    images = {}
    for xref, name in IMAGES.items():
        image = doc.extract_image(xref)
        path = args.output_dir / (name + "." + image["ext"])
        path.write_bytes(image["image"])
        images[path.name] = dict(xref=xref, width=image["width"], height=image["height"],
                                 sha256=digest(image["image"]))
    manifest = {"source_repository": "jgraving/circstats", "source_revision": PAPER_REV,
                "source_path": "figures/locust_density_order_figure.pdf",
                "source_sha256": PDF_SHA, "page": 1, "images": images,
                "credit": "Centre for the Advanced Study of Collective Behaviour; manuscript states used with permission. No additional redistribution licence is inferred."}
    for name, obj in (("images.json", manifest),
                      ("selection.json", extract_selection(notebook))):
        (args.output_dir / name).write_text(json.dumps(obj, indent=2) + "\n")
    print(f"Extracted {len(images)} native JPEGs and the screening table to {args.output_dir}")


if __name__ == "__main__":
    main()
