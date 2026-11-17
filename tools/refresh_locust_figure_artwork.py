"""Refresh artwork exports in the cloud without recomputing model summaries.

Use after typography/artwork edits; statistical paths and numerical summaries
are retained. Full figure regeneration remains locust_density_figure.py.
"""
import argparse
import base64
from copy import deepcopy
import hashlib
import io
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

import fitz
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import cairosvg

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "notebooks"))
sys.path.insert(0, str(ROOT / "tools"))
import locust_density_figure as figure
import prepare_locust_figure_assets as prepare
import upgrade_locust_figure_images as upgrade_tool

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preview", action="store_true", help="Emit review PNGs as base64")
    args = parser.parse_args(argv)
    reference = figure.ASSETS / "reference.pdf"
    assert figure.sha256(reference) == prepare.PDF_SHA
    # Re-extract with the original compact density unit preserved.
    with fitz.open(reference) as original:
        font_bytes = original.extract_font(10)[3]
        source_svg, changed = prepare.extract_original_panels(original)
    assert changed == 0
    source_svg = upgrade_tool.upgrade_svg_payloads(source_svg)
    (figure.ASSETS / "original_panels.svg").write_text(source_svg)
    manifest_path = figure.ASSETS / "images.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["original_panels"]["sha256"] = hashlib.sha256(source_svg.encode()).hexdigest()
    manifest["original_panels"]["extractor"] = "PyMuPDF " + fitz.VersionBind
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    stem = ROOT / "assets/locust-density-publication/locust_density_publication"
    before_svg = stem.with_suffix(".svg").read_text()
    root = ET.fromstring(before_svg)
    parents = {child: parent for parent in root.iter() for child in parent}
    # These named groups are the previous h label and generated e unit only.
    for identifier in ("text_33", "text_65"):
        node = next((item for item in root.iter() if item.get("id") == identifier), None)
        if node is not None and any("Cmex10" in value
                                    for child in node.iter()
                                    for value in child.attrib.values()):
            parents[node].remove(node)
    previous_label = next((item for item in root.iter()
                           if item.get("id") == "predictive-density-label"), None)
    if previous_label is not None:
        parents[previous_label].remove(previous_label)
    original_node = next(item for item in root if item.tag == "{"+figure.SVG+"}svg")
    index = list(root).index(original_node)
    source = ET.fromstring(source_svg)
    figure.polish_original_panels(source)
    source.set("x", "0"); source.set("y", "0")
    root.remove(original_node); root.insert(index, source)
    with tempfile.TemporaryDirectory() as directory:
        font_path = Path(directory) / "original-arial.ttf"
        font_path.write_bytes(font_bytes)
        prop = font_manager.FontProperties(fname=str(font_path))
        fig = plt.figure(figsize=(750/72, 525/72))
        fig.patch.set_visible(False)
        figure.add_compact_unit_label(fig, 380.576, 525-(401.301+89.75/2),
                                      "Predictive density", "rad", -1,
                                      fontsize=9, rotation=90, prop=prop)
        buf = io.StringIO()
        fig.savefig(buf, format="svg", metadata={"Date": None})
        plt.close(fig)
        label_root = ET.fromstring(buf.getvalue())
    label = next(node for node in label_root.iter()
                 if node.get("id") == "predictive-density-label")
    root.append(deepcopy(label))
    ET.register_namespace("", figure.SVG)
    svg = "\n".join(line.rstrip() for line in ET.tostring(root, encoding="unicode").splitlines()) + "\n"
    stem.with_suffix(".svg").write_text(svg)
    cairosvg.svg2pdf(bytestring=svg.encode(), write_to=str(stem.with_suffix(".pdf")))
    figure.restore_pdf_source_images(stem.with_suffix(".pdf"))
    with fitz.open(stem.with_suffix(".pdf")) as output:
        output[0].get_pixmap(matrix=fitz.Matrix(4, 4)).save(stem.with_suffix(".png"))
        if args.preview:
            preview = output[0].get_pixmap(matrix=fitz.Matrix(2, 2))
            print("FIGURE_PREVIEW_START")
            print(base64.b64encode(preview.tobytes("png")).decode())
            print("FIGURE_PREVIEW_END")
            for name, rect in (("units", (225, 199, 267, 232)),
                               ("ylabel", (368, 389, 387, 503))):
                pix = output[0].get_pixmap(matrix=fitz.Matrix(5, 5), clip=fitz.Rect(rect))
                print("PREVIEW_" + name + "_START")
                print(base64.b64encode(pix.tobytes("png")).decode())
                print("PREVIEW_" + name + "_END")
    html = stem.with_suffix(".html").read_text()
    html = re.sub(r"<svg\b[\s\S]*</svg>", lambda m: svg.rstrip(), html, count=1)
    stem.with_suffix(".html").write_text(html)
    metadata = json.loads(stem.with_suffix(".json").read_text())
    metadata["code_revision"] = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    for relative in list(metadata["source_sha256"]):
        metadata["source_sha256"][relative] = figure.sha256(ROOT / relative)
    for path in (reference, Path(__file__)):
        relative = path.relative_to(ROOT).as_posix()
        metadata["source_sha256"][relative] = figure.sha256(path)
    metadata["artwork_refresh"] = {
        "source": "Original circstats PDF",
        "source_font": "Embedded ArialMT (xref 10), used for the h label",
        "font_sha256": hashlib.sha256(font_bytes).hexdigest(),
        "superscript_scale": figure.SUPERSCRIPT_SCALE,
        "superscript_rise_em": figure.SUPERSCRIPT_RISE,
        "density_unit": "Retained original vector artwork",
        "pdf_photos": "Checked-in high-resolution JPEG streams (Keynote export)",
        "png_pixels": [3000, 2100],
        "numerical_summaries": "Unchanged"
    }
    metadata["outputs_sha256"] = {s: figure.sha256(stem.with_suffix(s))
                                 for s in (".svg", ".pdf", ".png", ".html")}
    stem.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n")
    print("ARTWORK_REFRESH", json.dumps(metadata["artwork_refresh"]))

if __name__ == "__main__":
    main()
