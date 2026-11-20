"""Hermetic checks for the Goldschtein bat figure generator.

The new generator rebuilds every notebook panel from the notebooks' own
plotting code, so nothing that needs the pinned data, the fits or a full
build is exercised here (the CI workflow dispatch job renders the full
pipeline twice and compares). These tests cover what can be checked
from the committed artwork alone: the asset pins, the artwork
structure, the panel-b curation and chrome spans, the measured layout
constants against the asset's image placements, the embed scales
against the reference strokes, ``prepare_inputs`` tamper-proofing, the
placement/affine math, the SVG helpers, the chrome layer rendering and
the PDF trailer stabiliser. No downloads, no sampling, no licensed
fonts (the chrome renders on the default font stack).
"""
import base64
import hashlib
import io
import json
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "notebooks"))
import goldshtein_bats_publication_figure as figure  # noqa: E402

SVG_NS = "http://www.w3.org/2000/svg"
XLINK_NS = "http://www.w3.org/1999/xlink"

ASSET_FILES = {
    "artwork.json",
    "reference.pdf",
    "image_23.jpeg", "image_279.jpeg", "image_45.jpeg",
}


# ---------------------------------------------------------------------------
# Committed asset pins and artwork structure
# ---------------------------------------------------------------------------
def test_asset_pins_cover_all_committed_files():
    assert set(figure.ASSET_SHA256) == ASSET_FILES
    assert all(len(h) == 64 for h in figure.ASSET_SHA256.values())
    manifest = figure.check_assets()
    assert set(manifest) == ASSET_FILES
    assert all(info["sha256"] == figure.ASSET_SHA256[name]
               for name, info in manifest.items())


def test_artwork_structure():
    asset = figure.load_artwork()
    assert asset["schema"] == 1
    assert len(asset["drawings"]) == figure.N_EXPECTED_DRAWINGS
    placements = [pl for im in asset["images"].values()
                  for pl in im["placements"]]
    assert len(asset["images"]) == figure.N_EXPECTED_IMAGES
    assert len(placements) == figure.N_EXPECTED_PLACEMENTS
    for pl in placements:
        for key in ("rect_ydown", "transform", "paths_before",
                    "fill_alpha"):
            assert key in pl, (pl, key)
    # the bat photo (xref 279) is the only placed image with an alpha
    photo = asset["images"]["279"]["placements"][0]
    assert photo["fill_alpha"] == pytest.approx(figure.PANEL_B_PHOTO_ALPHA)
    assert all(pl["fill_alpha"] == 1.0
               for k, im in asset["images"].items()
               for pl in im["placements"] if k != "279")
    spans = asset["text_spans"]
    assert len(spans) == figure.N_EXPECTED_SPANS
    # every chrome span carries a measured fill alpha; the 'Placeholder'
    # watermark is the only one that is not 1.0
    assert sum(1 for s in spans if "fill_alpha" in s) == 60
    assert [s["text"] for s in spans if s.get("fill_alpha", 1.0) != 1.0] == \
        ["Placeholder"]


def test_panel_b_curation():
    asset = figure.load_artwork()
    chosen = figure.panel_b_drawings(asset)
    assert len(chosen) == figure.N_EXPECTED_PANEL_B_DRAWINGS
    for d in chosen.values():
        x0, y0, x1, y1 = figure._drawing_bbox(d)
        cx, cy = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
        assert figure.PANEL_B_X_MIN <= cx <= figure.PANEL_B_X_MAX
        assert figure.PANEL_B_PHOTO[3] - 2.0 <= cy <= figure.PANEL_B_Y_MAX


def test_chrome_spans():
    asset = figure.load_artwork()
    spans = figure.chrome_spans(asset)
    wanted = set(figure.PANEL_LETTERS) | {"Treatment groups",
                                          "Placeholder"}
    assert {s["text"] for s in spans} == wanted
    assert len(spans) == len(wanted)
    letters = [s for s in spans if s["text"] in figure.PANEL_LETTERS]
    assert len(letters) == 8 and all(s["flags"] & 16 for s in letters)
    placeholder = next(s for s in spans if s["text"] == "Placeholder")
    assert placeholder.get("fill_alpha", 1.0) == pytest.approx(
        figure.PANEL_B_PHOTO_ALPHA)


# ---------------------------------------------------------------------------
# Measured layout against the asset placements and reference strokes
# ---------------------------------------------------------------------------
def _ydown_box_of_placement(pl):
    # the json stores the box as (x0, y_yup_bottom, x1, y_yup_top)
    x0, yb, x1, yt = pl["rect_ydown"]
    return (x0, figure.PAGE_H - yt, x1, figure.PAGE_H - yb)


def _sorted(boxes):
    return sorted(tuple(round(v, 2) for v in b) for b in boxes)


def test_layout_constants_match_asset_placements():
    asset = figure.load_artwork()

    def boxes_of(xref):
        return [_ydown_box_of_placement(pl)
                for pl in asset["images"][xref]["placements"]]

    assert _sorted(boxes_of("45")) == _sorted(figure.PANEL_A_AXES)
    assert _sorted(boxes_of("23")) == _sorted(figure.PANEL_H_AXES)
    assert _sorted([figure.PANEL_H_CBAR_BAR]) == _sorted(boxes_of("36"))
    assert _sorted([figure.PANEL_B_PHOTO]) == _sorted(boxes_of("279"))


def test_scales_match_reference_strokes():
    asset = figure.load_artwork()
    dr = asset["drawings"]
    # reference stroke width == code stroke width x embed scale
    assert dr["349"]["width"] == pytest.approx(1.0 * figure.SCALE_A,
                                               rel=1e-9)   # a border
    assert dr["77"]["width"] == pytest.approx(0.8 * figure.SCALE_C,
                                              rel=1e-9)    # c spine
    assert dr["146"]["width"] == pytest.approx(0.8 * figure.SCALE_F,
                                               rel=1e-9)   # f ring
    assert dr["306"]["width"] == pytest.approx(0.8 * figure.SCALE_G,
                                               rel=1e-9)   # g spine
    assert dr["193"]["width"] == pytest.approx(1.5 * figure.SCALE_H,
                                               rel=1e-9)   # h border
    # c/d/e share the Inkscape scale of the c spines
    assert figure.SCALE_C == figure.SCALE_D == figure.SCALE_E


# ---------------------------------------------------------------------------
# prepare_inputs tamper-proofing (no network)
# ---------------------------------------------------------------------------
def test_prepare_inputs_verifies_and_never_overwrites(tmp_path, monkeypatch):
    root = tmp_path / "ud"
    hashes = {}
    for i, rel in enumerate(figure.INPUT_HASHES):
        p = root / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        data = f"pinned-input-{i}".encode()
        p.write_bytes(data)
        hashes[rel] = hashlib.sha256(data).hexdigest()
    monkeypatch.setattr(figure, "INPUT_HASHES", hashes)

    assert figure.prepare_inputs(root) == root

    rel0 = list(hashes)[0]
    (root / rel0).write_bytes(b"tampered")
    with pytest.raises(ValueError, match="Checksum mismatch"):
        figure.prepare_inputs(root)
    (root / rel0).write_bytes(f"pinned-input-{list(hashes).index(rel0)}".encode())

    rel1 = list(hashes)[1]
    (root / rel1).unlink()
    with pytest.raises(FileNotFoundError):
        figure.prepare_inputs(root)


# ---------------------------------------------------------------------------
# Placement and affine math
# ---------------------------------------------------------------------------
def test_translate_maps_natural_box_to_measured():
    fig = plt.figure(figsize=(2.0, 1.0), dpi=100)     # 144 x 72 pt
    ax = fig.add_axes([0.25, 0.20, 0.5, 0.6])         # top at 20% of height
    scale = 1.5
    measured = (100.0, 50.0, 190.0, 250.0)            # page pt, y down
    tx, ty, w, h = figure._translate(fig, ax, measured, scale)
    nx0 = 0.25 * 144.0
    ny0_down = 0.20 * 72.0                            # anchor top, y down
    assert (tx, ty) == pytest.approx((measured[0] - scale * nx0,
                                      measured[1] - scale * ny0_down))
    assert (w, h) == pytest.approx((scale * 144.0, scale * 72.0))
    plt.close(fig)


def test_bar_box_error_zero_at_exact_placement():
    box = (50.0, 10.0, 120.0, 15.0)   # figure pt, y up, in a 100-pt figure
    h_pt = 100.0
    tx, ty_down, scale = 200.0, 300.0, 1.2
    measured_ydown = (tx + scale * 50.0,
                      ty_down + scale * (h_pt - 15.0),
                      tx + scale * 120.0,
                      ty_down + scale * (h_pt - 10.0))
    assert figure._bar_box_error(box, tx, ty_down, scale, h_pt,
                                 measured_ydown) < 1e-9


def test_affine_helpers():
    t1 = (1.0, 0.0, 0.0, 1.0, 10.0, 0.0)              # x += 10
    t2 = (2.0, 0.0, 0.0, 2.0, 0.0, 0.0)               # scale x2
    assert figure._affine_compose(t2, t1) == pytest.approx(
        (2.0, 0.0, 0.0, 2.0, 20.0, 0.0))
    assert figure._affine_box(figure._affine_compose(t2, t1),
                              (0.0, 0.0, 1.0, 1.0)) == \
        pytest.approx((20.0, 0.0, 22.0, 2.0))


# ---------------------------------------------------------------------------
# SVG helpers
# ---------------------------------------------------------------------------
def test_prefix_svg_ids():
    svg = (f'<svg xmlns="{SVG_NS}"><defs>'
           '<linearGradient id="grad"><stop offset="0"/></linearGradient>'
           '</defs><g id="g"><rect fill="url(#grad)"/></g></svg>')
    root = ET.fromstring(svg)
    figure._prefix_svg_ids(root, "p--")
    out = ET.fromstring(ET.tostring(root))
    assert out.find(f".//{{{SVG_NS}}}linearGradient").get("id") == "p--grad"
    assert out.find(f".//{{{SVG_NS}}}g").get("id") == "p--g"
    assert (out.find(f".//{{{SVG_NS}}}rect").get("fill")
            == "url(#p--grad)")


def test_swap_panel_rasters():
    ns = "{%s}" % SVG_NS
    href = "{%s}href" % XLINK_NS
    g = ET.Element(ns + "g")
    imgs = []
    for i, (height, mime) in enumerate(
            ((300.0, "data:image/png;base64,AAA"),
             (300.0, "data:image/png;base64,BBB"),
             (10.0, "data:image/jpeg;base64,CCC"))):
        el = ET.SubElement(g, ns + "image",
                           {"x": "1", "y": str(2 + i), "width": "100",
                            "height": str(height)})
        el.set(href, mime)
        imgs.append(el)
    content = [g]
    figure._swap_panel_rasters(content, 2,
                               figure.ASSET_DIR / "image_279.jpeg")
    want = ("data:image/jpeg;base64," + base64.b64encode(
        (figure.ASSET_DIR / "image_279.jpeg").read_bytes()
    ).decode("ascii"))
    assert imgs[0].get(href) == want
    assert imgs[1].get(href) == want
    assert imgs[2].get(href) == "data:image/jpeg;base64,CCC"   # untouched
    assert imgs[0].get("preserveAspectRatio") == "none"
    with pytest.raises(ValueError, match="expected 3 basemap rasters"):
        figure._swap_panel_rasters(content, 3,
                                   figure.ASSET_DIR / "image_279.jpeg")


# ---------------------------------------------------------------------------
# _stabilize_pdf
# ---------------------------------------------------------------------------
def test_stabilize_pdf_is_a_function_of_content():
    trailer = (b"/ID[<0123456789ABCDEF0123456789ABCDEF>"
               b"<FEDCBA9876543210FEDCBA9876543210>]>>")
    data = (b"%PDF-1.7\n/CreationDate (D:20260101120000Z)\n"
            b"/ModDate (D:20260101120000Z)\n"
            b"trailer <</Size 2/Root 1 0 R" + trailer
            + b">\nstartxref\n0\n%%EOF")
    other = data.replace(b"FEDCBA9876543210FEDCBA9876543210",
                         b"00112233445566778899001122334455")
    other = other.replace(b"(D:20260101120000Z)", b"(D:19990101000000Z)")
    out, out_other = figure._stabilize_pdf(data), figure._stabilize_pdf(other)
    # the random MuPDF ID and the wall-clock dates leave no trace
    assert out == out_other
    assert b"(D:20260101120000Z)" not in out
    m = re.search(rb"/ID \[<([0-9a-f]{64})> <([0-9a-f]{64})>\]", out)
    assert m and m.group(1) == m.group(2)
    # the digest is the SHA-256 of the document with a normalised ID
    base = re.sub(rb"/ID\s*\[[^\]]*\]", b"/ID [<0> <0>]",
                  re.sub(rb"/CreationDate \(.*?\)",
                         b"/CreationDate (D:20000101000000Z)",
                         re.sub(rb"/ModDate \(.*?\)",
                                b"/ModDate (D:20000101000000Z)", other)))
    assert m.group(1) == hashlib.sha256(base).hexdigest().encode("ascii")


# ---------------------------------------------------------------------------
# Chrome layer (pinned artwork, default font stack)
# ---------------------------------------------------------------------------
def test_chrome_layer_renders():
    asset = figure.load_artwork()
    ctx = {"panel_b_drawings": figure.panel_b_drawings(asset),
           "chrome_spans": figure.chrome_spans(asset)}
    fig, stats = figure.draw_chrome(ctx)
    assert stats["drawings"] == figure.N_EXPECTED_PANEL_B_DRAWINGS
    assert stats["text_drawn"] == 10
    assert stats["text_skipped"] == 0
    # the chrome layer must be transparent: the page white backing and
    # the photograph are drawn by the composition itself, and a white
    # chrome patch would hide the photograph in the PDF export
    assert fig.patch.get_facecolor()[3] == 0.0
    buf = io.BytesIO()
    fig.savefig(buf, format="svg", metadata={"Date": None})
    root = ET.fromstring(buf.getvalue())
    figure_g = root.find(f".//{{{SVG_NS}}}g[@id='figure_1']")
    patch = figure_g.find(f"{{{SVG_NS}}}g[@id='patch_1']")
    # the figure background path renders fill: none when transparent
    assert patch is None or \
        "fill: none" in (patch[0].get("style") or "")
    assert len(buf.getvalue()) > 1000
    plt.close(fig)


def test_svg_part_split():
    svg = (f'<svg xmlns="{SVG_NS}"><metadata>x</metadata>'
           f'<defs><linearGradient id="g"/></defs>'
           f'<rect id="r1"/><g id="g1"/></svg>')
    defs, content = figure._split_svg(svg.encode())
    assert defs is not None
    assert [c.tag for c in content] == [f"{{{SVG_NS}}}rect",
                                        f"{{{SVG_NS}}}g"]
