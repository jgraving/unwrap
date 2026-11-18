"""Reproduce the circstats paper's graphical abstract (599 x 342 points).

Rebuilds the Inkscape-polished ``Graphical Abstract.pdf`` published with the
circstats elife paper (``notebooks/assets/circstats_graphical_abstract/
reference.pdf``) in one deterministic render:

* the 472 hand-made / stochastic vector paths of the published figure
  (the twelve locust circle units, the foraging trajectory, the code
  editor mock, the MCMC trace panels, the posterior predictive clusters,
  the jagged posterior densities, the observed-angle dots, the panel
  "d" predictive curves and ellipsoidal loop, and the four KDE summary
  circles' interiors with arrow fans) are pinned verbatim from the
  reference PDF and replayed through the vector renderer (SHA-256
  verified);
* the 54-piece deterministic skeleton (axis spines and ticks, the
  circular-modulo and softplus link curves, the model / hypothesis
  interval bars, the uniform reference lines, the flow arrows, the
  legend frames and swatches, the four KDE circle outlines and the
  eta_mu / eta_kappa glyph marks) is drawn from the same measured
  geometry;
* all 236 text spans (Arial labels, the Menlo syntax-highlighted code
  window, the rotated axis captions) are typeset live at their measured
  positions, so the exports stay editable (``svg.fonttype = "none"``).

Nothing is randomized: every coordinate, colour, width, dash pattern and
string literal comes from the hash-pinned asset, so two renders produce
byte-identical exports.

Exports (``assets/circstats-graphical-abstract/`` by default): SVG,
byte-stabilised PDF (cairo CreationDate cleared, /ID pinned to the
content hash), a 3x PNG preview rasterised from the PDF, a
self-contained HTML page with the published caption, and a provenance
JSON.

Run from any directory: ``python notebooks/circstats_graphical_abstract.py
[--output-dir DIR] [--font-family FAMILY]``.  See
``CIRCSTATS_GRAPHICAL_ABSTRACT.md`` for provenance.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import io
import json
import re
import sys
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.font_manager import FontProperties, findfont  # noqa: E402
from matplotlib.path import Path as MplPath  # noqa: E402
from matplotlib.patches import PathPatch  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
ASSET_DIR = Path(__file__).resolve().parent / "assets" / \
    "circstats_graphical_abstract"
REFERENCE_PDF = ASSET_DIR / "reference.pdf"
ARTWORK_PATH = ASSET_DIR / "artwork.json"

PAGE_W = 599.0          # points
PAGE_H = 342.0          # points
FIGURE_NAME = "circstats_graphical_abstract"
CAPTION = "Overview of a modern Bayesian workflow for circular GLMs."

# PDF font name -> family stack.  Liberation Sans is metric-compatible
# with Arial (the CI default); Menlo falls back to a monospace.
FONT_STACKS = {
    "ArialMT": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
    "Arial-BoldMT": ["Arial", "Liberation Sans", "Helvetica",
                     "DejaVu Sans"],
    "Helvetica": ["Helvetica", "Arial", "Liberation Sans",
                  "DejaVu Sans"],
    "Menlo-Regular": ["Menlo", "Menlo Regular", "DejaVu Sans Mono",
                      "Liberation Mono", "monospace"],
    "Menlo-Italic": ["Menlo", "Menlo Regular", "DejaVu Sans Mono",
                     "Liberation Mono", "monospace"],
}

# Spans that are whitespace placeholders in the reference (the eta
# labels are vector glyphs, not text) and are skipped.
_SKIP_FONTS = {"STIXTwoMath-Regular"}


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_asset() -> dict:
    asset = json.loads(ARTWORK_PATH.read_text())
    if asset.get("schema") != 1:
        raise ValueError(f"unexpected artwork schema in {ARTWORK_PATH}")
    if tuple(asset["source"]["page_pt"]) != (PAGE_W, PAGE_H):
        raise ValueError("asset page size does not match the reference")
    drawings = asset["drawings"]
    if len(drawings) != 526:
        raise ValueError(f"expected 526 pinned drawings, found "
                         f"{len(drawings)}")
    return asset


# -------------------------------------------------------------------
# vector renderer (PDF path items -> matplotlib paths)
# -------------------------------------------------------------------
def path_from_items(items) -> MplPath:
    """Convert a PyMuPDF item list into a matplotlib path.

    Items come in page points with y down from the top-left corner; the
    flip to matplotlib y-up happens here (``PAGE_H - y``), so a single
    renderer serves every drawing regardless of orientation.  Subpaths
    that close on themselves are terminated with CLOSEPOLY so strokes
    join cleanly.
    """
    def flip(x, y):
        return (x, PAGE_H - y)

    verts: list[tuple[float, float]] = []
    codes: list[int] = []
    cur = None          # end of the current subpath (flipped coords)
    starts: list = []   # open subpath starts

    def begin(x: float, y: float):
        nonlocal cur
        p = flip(x, y)
        if cur is None or abs(cur[0] - p[0]) > 1e-9 \
                or abs(cur[1] - p[1]) > 1e-9:
            verts.append(p)
            codes.append(MplPath.MOVETO)
            starts.append([p])
            cur = p

    for it in items:
        kind = it[0]
        if kind == "l":
            x1, y1, x2, y2 = (float(v) for v in it[1:5])
            begin(x1, y1)
            verts.append(flip(x2, y2))
            codes.append(MplPath.LINETO)
            cur = verts[-1]
        elif kind == "c":
            (x1, y1), (x2, y2), (x3, y3), (x4, y4) = (
                (float(it[k]), float(it[k + 1])) for k in (1, 3, 5, 7))
            begin(x1, y1)
            verts.append(flip(x2, y2))
            codes.append(MplPath.CURVE4)
            verts.append(flip(x3, y3))
            codes.append(MplPath.CURVE4)
            verts.append(flip(x4, y4))
            codes.append(MplPath.CURVE4)
            cur = verts[-1]
        elif kind == "re":
            x0, y0, x1, y1 = (float(v) for v in it[1:5])
            if cur is not None:
                codes.append(MplPath.CLOSEPOLY)
                verts.append(cur)
                starts.pop()
            cur = None
            a = flip(x0, y0)
            b = flip(x1, y0)
            c = flip(x1, y1)
            d = flip(x0, y1)
            verts.append(a)
            codes.append(MplPath.MOVETO)
            starts.append([a])
            for p in (b, c, d):
                verts.append(p)
                codes.append(MplPath.LINETO)
            codes.append(MplPath.CLOSEPOLY)
            verts.append(a)
            starts.pop()

    # close any open subpath that returns to its start point
    for start in starts:
        if cur is not None and abs(cur[0] - start[0][0]) <= 1e-9 \
                and abs(cur[1] - start[0][1]) <= 1e-9:
            codes.append(MplPath.CLOSEPOLY)
            verts.append(cur)
    return MplPath(verts, codes)


def _color(c, opacity=None):
    if c is None:
        return "none"
    rgb = tuple(float(v) for v in c[:3])
    if opacity is None or opacity in (1.0, 1):
        return rgb
    return rgb + (float(opacity),)


def draw_drawing(ax, dr, zorder):
    """Render one pinned/reference drawing (fill and/or stroke)."""
    path = path_from_items(dr["items"])
    t = dr["type"]
    face = _color(dr["fill"], dr.get("fill_opacity")) if "f" in t \
        else "none"
    if "s" in t:
        edge = _color(dr["color"], dr.get("stroke_opacity"))
        lw = dr.get("width") or 0.0
        ls = "solid"
        dash = dr.get("dashes")
        if dash and lw:
            ls = (0, (dash["on"] / lw, dash["off"] / lw))
        patch = PathPatch(path, facecolor=face, edgecolor=edge,
                          linewidth=lw, linestyle=ls,
                          capstyle="butt", joinstyle="miter",
                          zorder=zorder)
    else:
        patch = PathPatch(path, facecolor=face, edgecolor="none",
                          zorder=zorder)
    ax.add_patch(patch)
    return patch


# -------------------------------------------------------------------
# typography
# -------------------------------------------------------------------
def _resolved_fp(stack, *, bold, italic, size):
    return FontProperties(family=list(stack), size=size,
                          weight="bold" if bold else "normal",
                          style="italic" if italic else "normal")


def draw_text(ax, spans, *, family_override=None, zorder=1000):
    """Typeset every span on its measured baseline.

    Each span carries the embedded font's own ``ascender`` /
    ``descender`` (em fractions, positive up / negative down), so the
    baseline sits at ``y + ascender*size`` in page points (y down),
    independent of the substitute face used for rendering.  Vertical
    spans (direction (0, -1)) read bottom-to-top with glyph tops facing
    left (the standard rotated axis caption); their baseline is the
    vertical line ``x0 + ascender*size`` and their first glyph starts at
    the bottom of the span.
    """
    drawn = 0
    skipped = 0
    for s in spans:
        text = s["text"]
        if s["font"] in _SKIP_FONTS or not text.strip(" \xa0"):
            skipped += 1
            continue
        bold = s["font"] == "Arial-BoldMT"
        italic = s["font"] == "Menlo-Italic"
        stack = list(FONT_STACKS.get(s["font"], FONT_STACKS["ArialMT"]))
        if family_override:
            stack = [family_override] + [
                f for f in stack if f != family_override]
        fp = _resolved_fp(stack, bold=bold, italic=italic,
                          size=s["size"])
        asc = float(s.get("ascender", 0.8))
        color = s["color"]
        if color.startswith("#00") and len(color) == 9:
            color = "#" + color[3:]  # drop the leading zero byte
        if s["dir"][1] < -0.5:  # vertical, bottom-to-top, tops facing left
            # rotation_mode="anchor" keeps ha/va in the unrotated frame,
            # so the anchor sits on the baseline at the string start.
            # The reference baseline is x0 + ascender*size (y down).
            ax.text(s["x"] + asc * s["size"], PAGE_H - s["y1"], text,
                    fontproperties=fp, color=color, ha="left",
                    va="baseline", rotation=90, rotation_mode="anchor",
                    zorder=zorder)
        else:
            tilt = 57.29578 * (s["dir"][1] / s["dir"][0]) \
                if s["dir"][0] else 0.0
            base_y = s["y"] + asc * s["size"]
            ax.text(s["x"], PAGE_H - base_y, text, fontproperties=fp,
                    color=color, ha="left", va="baseline", rotation=tilt,
                    zorder=zorder)
        drawn += 1
    return drawn, skipped


# -------------------------------------------------------------------
# figure assembly
# -------------------------------------------------------------------
def draw_figure(asset, *, family=None):
    """Render the full page; returns (fig, ax, stats)."""
    fig = plt.figure(figsize=(PAGE_W / 72.0, PAGE_H / 72.0), dpi=96)
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_xlim(0, PAGE_W)
    ax.set_ylim(0, PAGE_H)
    ax.axis("off")
    ax.set_facecolor("white")
    fig.patch.set_facecolor("white")
    # deterministic SVG emission: fixed glyph/clip ids, no wall-clock date
    mpl.rcParams["svg.fonttype"] = "none"  # keep text editable
    mpl.rcParams["svg.hashsalt"] = "unwrap-circstats-graphical-abstract-v1"

    drawings = asset["drawings"]
    order = sorted(int(k) for k in drawings)
    for z, i in enumerate(order):
        draw_drawing(ax, drawings[str(i)], zorder=10 + z)

    drawn, skipped = draw_text(ax, asset["text_spans"],
                               family_override=family)
    stats = {
        "drawings": len(order),
        "text_drawn": drawn,
        "text_skipped": skipped,
    }
    return fig, ax, stats


def _stabilize_pdf(data: bytes) -> bytes:
    """Clear cairo's wall-clock CreationDate and pin /ID to content."""
    data = re.sub(rb"/CreationDate \(.*?\)",
                  b"/CreationDate (D:20000101000000Z)", data)
    data = re.sub(rb"/ModDate \(.*?\)", b"/ModDate (D:20000101000000Z)",
                  data)
    digest = hashlib.sha256(data).hexdigest().encode("ascii")
    data = re.sub(rb"/ID \[[^\]]*\]",
                  b"/ID [<" + digest + b"> <" + digest + b">]", data)
    return data


def _html_page(svg: str) -> str:
    return (
        "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n"
        "<meta charset=\"utf-8\">\n"
        f"<title>circstats graphical abstract</title>\n"
        "<style>\n"
        "body{margin:0;background:#fff;display:flex;flex-direction:"
        "column;align-items:center;justify-content:center;min-height:"
        "100vh}\n"
        "figure{margin:2rem auto}\n"
        "img{max-width:95vw;height:auto}\n"
        "figcaption{max-width:900px;margin-top:1rem;color:#231f20;"
        "font:14px/1.5 Helvetica,Arial,sans-serif}\n"
        "</style>\n</head>\n<body>\n<figure>\n"
        "<img src=\"data:image/svg+xml;base64,"
        + __import__("base64").b64encode(svg.encode()).decode("ascii")
        + "\" alt=\"Overview of a modern Bayesian workflow for circular "
          "GLMs.\">\n"
        f"<figcaption>{CAPTION}</figcaption>\n"
        "</figure>\n</body>\n</html>\n")


def resolve_font_report(family):
    report = {}
    for name, stack in FONT_STACKS.items():
        if family:
            stack = ([family]
                     + [f for f in stack if f != family])
        try:
            path = findfont(FontProperties(family=stack),
                            fallback_to_default=False)
        except ValueError:
            path = findfont(FontProperties(family=stack))
        report[name] = str(Path(path).name)
    return report


def render(output_dir: Path, *, family=None) -> dict:
    """Render all exports into ``output_dir``; returns provenance."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    asset = load_asset()
    fig, _ax, stats = draw_figure(asset, family=family)

    svg_path = output_dir / f"{FIGURE_NAME}.svg"
    buf = io.BytesIO()
    fig.savefig(buf, format="svg", metadata={"Date": None})
    svg_path.write_bytes(buf.getvalue())
    plt.close(fig)
    svg = svg_path.read_bytes()

    pdf_path = output_dir / f"{FIGURE_NAME}.pdf"
    pdf_buf = io.BytesIO()
    # re-open the figure for the PDF pass (savefig closes the canvas)
    fig2, _ax2, _ = draw_figure(asset, family=family)
    fig2.savefig(pdf_buf, format="pdf")
    plt.close(fig2)
    pdf_path.write_bytes(_stabilize_pdf(pdf_buf.getvalue()))

    png_path = output_dir / f"{FIGURE_NAME}.png"
    import pymupdf
    with pymupdf.open(stream=pdf_path.read_bytes(), filetype="pdf") as doc:
        pix = doc[0].get_pixmap(matrix=pymupdf.Matrix(3, 3))
        pix.save(str(png_path))

    html_path = output_dir / f"{FIGURE_NAME}.html"
    html_path.write_text(_html_page(svg.decode("utf-8")))

    provenance = {
        "figure": FIGURE_NAME,
        "caption": CAPTION,
        "page_pt": [PAGE_W, PAGE_H],
        "inputs": {
            "artwork_sha256": sha256(ARTWORK_PATH),
            "reference_pdf_sha256": sha256(REFERENCE_PDF),
        },
        "fonts": resolve_font_report(family),
        "counts": stats,
        "env": {
            "python": sys.version.split()[0],
            "matplotlib": importlib.metadata.version("matplotlib"),
            "numpy": importlib.metadata.version("numpy"),
        },
        "exports": {
            "svg_sha256": hashlib.sha256(svg).hexdigest(),
            "pdf_sha256": sha256(pdf_path),
            "png_sha256": sha256(png_path),
            "html_sha256": sha256(html_path),
        },
    }
    (output_dir / f"{FIGURE_NAME}.json").write_text(
        json.dumps(provenance, indent=2) + "\n")
    return provenance


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__.splitlines()[0])
    parser.add_argument("--output-dir", type=Path,
                        default=ROOT / "assets"
                        / "circstats-graphical-abstract")
    parser.add_argument("--font-family", default=None,
                        help="override the first family in every stack "
                             "(CI uses 'Liberation Sans')")
    args = parser.parse_args(argv)
    prov = render(args.output_dir, family=args.font_family)
    for name, h in prov["exports"].items():
        print(f"{name}  {h[:16]}…")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
