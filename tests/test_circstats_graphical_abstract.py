"""Hermetic checks for the circstats graphical abstract generator.

The generator replays the 526 vector drawings and 236 text spans of the
Inkscape-polished reference PDF, so the tests pin the asset hashes,
verify that the pinned/skeleton partition covers every drawing exactly
once, re-check the reference PDF's geometry with PyMuPDF, and assert
that two renders in the same environment are byte-identical.  No
network access or licensed fonts are needed (DejaVu/Liberation resolve
the stacks either way).
"""
import hashlib
import io
import json
import re
import struct
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "notebooks"))
import circstats_graphical_abstract as ga  # noqa: E402

SVG_NS = "http://www.w3.org/2000/svg"

ARTWORK_SHA = ("d2d6991cb522ac180dea76bf130511049a3562cded5b2b501cb8"
               "22685a9fd612")
REFERENCE_SHA = ("5874a6da22930ff3b0e53f37c3ffd0074c2e88ddb25c"
                 "dd4b909f7187a9e4e55c")
FONTS = {"ArialMT", "Arial-BoldMT", "Helvetica",
         "Menlo-Regular", "Menlo-Italic", "STIXTwoMath-Regular"}


def _pdf_spans(page):
    return [s for b in page.get_text("dict")["blocks"]
            if b["type"] == 0
            for l in b["lines"] for s in l["spans"]]


def test_committed_assets_match_manifest():
    assert ga.sha256(ga.ASSET_DIR / "artwork.json") == ARTWORK_SHA
    assert ga.sha256(ga.REFERENCE_PDF) == REFERENCE_SHA


def test_asset_internal_consistency():
    asset = ga.load_asset()
    assert tuple(asset["source"]["page_pt"]) == (599.0, 342.0)
    assert len(asset["drawings"]) == 526
    assert sorted(asset["drawings"], key=int) == [str(i) for i in range(526)]
    assert len(asset["text_spans"]) == 236
    assert len(asset["units"]) == 12
    assert sum(len(u["indices"]) for u in asset["units"]) > 0
    assert set(asset["regions"]) == {
        "trajectory", "code_window", "eta_mu_trace", "eta_k_trace",
        "mu_posterior", "k_posterior", "posterior_cluster",
        "observed_dots", "d_line_panel", "d_ellipse", "c_blobs_fans"}


def test_pinned_split_is_a_partition_of_all_drawings():
    asset = json.loads(ga.ARTWORK_PATH.read_text())
    claimed = [i for u in asset["units"] for i in u["indices"]]
    claimed += [i for idxs in asset["regions"].values() for i in idxs]
    claimed += list(asset["skeleton_indices"])
    assert len(claimed) == len(set(claimed)) == 526
    assert set(claimed) == set(range(526))


def test_reference_pdf_is_pure_vector_with_matching_counts():
    import pymupdf

    doc = pymupdf.open(ga.REFERENCE_PDF)
    page = doc.load_page(0)
    assert (page.rect.width, page.rect.height) == (599.0, 342.0)
    assert page.get_images(full=True) == []
    assert len(page.get_drawings()) == 526
    assert len(_pdf_spans(page)) == 236
    assert {s["font"] for s in _pdf_spans(page)} <= FONTS
    doc.close()


def test_svg_is_deterministic_vector_with_editable_text():
    asset = ga.load_asset()

    def svg_bytes():
        fig, _ax, _stats = ga.draw_figure(asset)
        buf = io.BytesIO()
        fig.savefig(buf, format="svg",
                    metadata={"Date": None})
        plt.close(fig)
        return buf.getvalue()

    a, b = svg_bytes(), svg_bytes()
    assert a == b

    root = ET.fromstring(a)
    assert root.tag == f"{{{SVG_NS}}}svg"
    assert root.get("width") == "599pt"
    assert root.get("height") == "342pt"
    assert root.get("viewBox") == "0 0 599 342"
    assert not root.findall(f".//{{{SVG_NS}}}image")
    assert len(root.findall(f".//{{{SVG_NS}}}path")) >= 526
    # svg.fonttype = "none": the 149 non-placeholder spans stay <text>
    assert len(root.findall(f".//{{{SVG_NS}}}text")) == 149


def test_stabilize_pdf_pins_dates_and_id():
    data = (b"%PDF-1.7\n/CreationDate (D:20261118030400Z)"
            b"/ModDate (D:20261118030400Z) /ID [<AA> <BB>]")
    out = ga._stabilize_pdf(data)
    assert b"20261118030400" not in out
    assert b"/CreationDate (D:20000101000000Z)" in out
    assert b"/ModDate (D:20000101000000Z)" in out
    # the /ID slots are the SHA-256 of the date-pinned document, so a
    # second (or a writer-stamped) run hashes to the same content
    m = re.search(rb"/ID \[<([0-9a-f]{64})> <([0-9a-f]{64})>\]", out)
    assert m is not None
    assert m.group(1) == m.group(2)
    pinned = (b"%PDF-1.7\n/CreationDate (D:20000101000000Z)"
              b"/ModDate (D:20000101000000Z) /ID [<AA> <BB>]")
    assert m.group(1) == hashlib.sha256(pinned).hexdigest().encode()


def test_render_writes_all_exports_with_consistent_provenance(tmp_path):
    prov = ga.render(tmp_path)
    base = tmp_path / ga.FIGURE_NAME
    for ext in (".svg", ".pdf", ".png", ".html", ".json"):
        assert base.with_suffix(ext).exists()
    assert prov["page_pt"] == [599.0, 342.0]
    assert prov["inputs"]["artwork_sha256"] == ARTWORK_SHA
    assert prov["inputs"]["reference_pdf_sha256"] == REFERENCE_SHA
    assert prov["counts"] == {"drawings": 526, "text_drawn": 149,
                              "text_skipped": 87}
    assert len(prov["fonts"]) == len(ga.FONT_STACKS)
    assert all(prov["fonts"].values())
    for name, digest in prov["exports"].items():
        ext = name.split("_sha256")[0]
        assert digest == ga.sha256(base.with_suffix("." + ext))
    assert base.with_suffix(".pdf").read_bytes()[:5] == b"%PDF-"
    png = base.with_suffix(".png").read_bytes()
    assert png[:8] == b"\x89PNG\r\n\x1a\n"
    assert struct.unpack(">II", png[16:24]) == (3 * 599, 3 * 342)
    html = base.with_suffix(".html").read_text()
    assert html.startswith("<!DOCTYPE html>")
    assert ga.CAPTION in html


def test_pdf_round_trip_stays_vector_and_keeps_text(tmp_path):
    prov = ga.render(tmp_path)
    import pymupdf

    pdf = tmp_path / f"{ga.FIGURE_NAME}.pdf"
    doc = pymupdf.open(pdf)
    page = doc.load_page(0)
    assert (page.rect.width, page.rect.height) == (599.0, 342.0)
    assert page.get_images(full=True) == []
    # the 149 emitted text runs may merge or split on re-read, so check
    # the floor and the actual content
    assert len(_pdf_spans(page)) >= 149
    text = page.get_text()
    for needle in ("Experimental design", "Model summaries",
                   "Predictive density", "Softplus"):
        assert needle in text
    doc.close()
    assert prov["exports"]["pdf_sha256"] == ga.sha256(pdf)
