"""Rebuild the circstats graphical abstract with an executable model example.

The illustrations, surrounding panels, and code window use vector artwork
and typography from the paper. Panel B's plots are regenerated from a
seeded circular model and the package's parameter transforms.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from functools import lru_cache
import hashlib
import importlib.metadata
import io
import json
import sys
from pathlib import Path

import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import to_rgb  # noqa: E402
from matplotlib.font_manager import (FontProperties, findSystemFonts, findfont,
                                     fontManager)  # noqa: E402
from matplotlib.path import Path as MplPath  # noqa: E402
from matplotlib.patches import PathPatch  # noqa: E402
from matplotlib.textpath import TextToPath  # noqa: E402
from matplotlib.transforms import Affine2D, Bbox, TransformedBbox  # noqa: E402

import circstats_abstract_model as model  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
ASSET_DIR = Path(__file__).resolve().parent / "assets" / \
    "circstats_graphical_abstract"
REFERENCE_PDF = ASSET_DIR / "reference.pdf"
ARTWORK_PATH = ASSET_DIR / "artwork.json"

PAGE_W = 599.0          # points
PAGE_H = 342.0          # points
FIGURE_NAME = "circstats_graphical_abstract"
CAPTION = (
    "Overview of a modern Bayesian workflow for circular GLMs. In panel B, "
    "the faint polar curves are von Mises densities at individual prior and "
    "posterior draws of μ and κ. The "
    "blue curve is the posterior mean and the dashed curve marks the "
    "uniform reference. Zero density starts on the inner ring. "
    f"The 16 illustrative observations use {model.OBSERVED_SPREAD:.0%} "
    "of the angles’ original "
    "offsets from their circular mean. Translucent dots inside the ring "
    "show jittered copies of these observations. Each lane shows 400 "
    "consecutive draws from three illustrative MCMC chains. Sample dots "
    "follow the chains through the links to the likelihood parameters. "
    "Direction traces use three equivalent turns and softplus maps "
    "concentration draws to positive values."
)
# The old Cartesian plots in B and an overlapping copy of D's inner ring.
_REPLACED_PATHS = set(range(11, 38)) | set(range(119, 240)) | {255}
_WORKFLOW_ARROWHEADS = {1, 3, 6, 8, 10}
_SUMMARY_POLAR_RINGS = {260, 261}
_SUMMARY_POLAR_SPOKES = {256: (0, -1), 257: (1, 0), 258: (0, 1), 259: (-1, 0)}
# The measured center and radii of the two gray rings in panel D, in points.
_SUMMARY_POLAR_CENTER = (167.492, 262.062)
_SUMMARY_POLAR_RADII = (16.998, 33.996)
# Model-comparison and hypothesis axes, and their superseded source spans.
_SUMMARY_AXIS_LABELS = {107: "Predictive density difference",
                        273: "Predictive density"}
_SUMMARY_AXIS_LABEL_SPANS = {53, 54, 230}
_SUMMARY_AXIS_Y = (218.0, 317.0)  # common span enclosing both sets of intervals
_ESTIMATE_DOTS = {115, 116, 117, 118, 281, 282, 283}
_ESTIMATE_DOT_SCALE = 0.68
_PANEL_STROKE_WIDTHS = {
    **dict.fromkeys(set(range(109, 115)) | set(range(275, 281)), 1.05),  # intervals
    **dict.fromkeys(_ESTIMATE_DOTS, 0.40),
    **dict.fromkeys((44, 45, 48, 49), 1.10),  # population arrows
    **dict.fromkeys(set(range(52, 77)) | set(range(79, 107)), 0.64),  # individual arrows
}
_EDITOR_SOURCE = (268.81, 24.07)
_EDITOR_POSITION = (254.0, 20.5)
_EDITOR_ZOOM = 1.19
_EDITOR_WIDTH = 107.37
_EDITOR_BAR_HEIGHT = 8.05

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


def _reference_drawings():
    """Resolve the PDF graphics state in drawing order, including group alpha."""
    import pymupdf

    with pymupdf.open(REFERENCE_PDF) as doc:
        sequence = doc[0].get_drawings(extended=True)

    parents = []
    for item in sequence:
        parents = [parent for parent in parents if parent["level"] < item["level"]]
        if item["type"] in {"group", "clip"}:
            parents.append(item)
            continue
        opacity = 1.0
        clip_box = None
        clip_path = None
        for parent in parents:
            if parent["type"] == "group":
                # The reference's translucent groups contain one path each.
                opacity *= parent["opacity"]
            else:
                clip_box = parent["scissor"]
                if any(segment[0] != "re" for segment in parent["items"]):
                    items = [[segment[0], *(v for point in segment[1:] for v in point)]
                             for segment in parent["items"]]
                    clip_path = path_from_items(items)
        yield item["type"], opacity, clip_box, clip_path


def draw_drawing(ax, dr, zorder, *, opacity, clip_box, clip_path,
                 joinstyle="miter", transform=None, center_scale=1.0):
    """Render one pinned/reference drawing (fill and/or stroke)."""
    data_transform = ax.transData if transform is None else transform + ax.transData
    path = path_from_items(dr["items"])
    if center_scale != 1.0:
        bounds = path.get_extents()
        x, y = bounds.p0 + bounds.size / 2
        path = (Affine2D().translate(-x, -y).scale(center_scale).translate(x, y)
                .transform_path(path))
    t = dr["type"]
    face = _color(dr["fill"], opacity * dr["fill_opacity"]) if "f" in t \
        else "none"
    if "s" in t:
        edge = _color(dr["color"], opacity * dr["stroke_opacity"])
        lw = dr.get("width") or 0.0
        ls = "solid"
        dash = dr.get("dashes")
        if dash and lw:
            ls = (0, (dash["on"] / lw, dash["off"] / lw))
        patch = PathPatch(path, facecolor=face, edgecolor=edge,
                          linewidth=lw, linestyle=ls,
                          capstyle="butt", joinstyle=joinstyle,
                          transform=data_transform, zorder=zorder)
    else:
        patch = PathPatch(path, facecolor=face, edgecolor="none",
                          transform=data_transform, zorder=zorder)
    ax.add_patch(patch)
    if clip_box is not None:
        x0, y0, x1, y1 = clip_box
        box = Bbox.from_extents(x0, PAGE_H - y1, x1, PAGE_H - y0)
        patch.set_clip_box(TransformedBbox(box, data_transform))
    if clip_path is not None:
        patch.set_clip_path(clip_path, data_transform)
    return patch


def _summary_polar_guide(index, drawing):
    """Use the shared polar guide style between the two measured rings."""
    if index in _SUMMARY_POLAR_RINGS:
        return {**drawing, "color": to_rgb(model.POLAR_GUIDE_COLOR),
                "width": model.POLAR_RING_WIDTH}
    dx, dy = _SUMMARY_POLAR_SPOKES[index]
    x, y = _SUMMARY_POLAR_CENTER
    inner, outer = _SUMMARY_POLAR_RADII
    items = [("l", x + dx * inner, y + dy * inner,
              x + dx * outer, y + dy * outer)]
    return {**drawing, "color": to_rgb(model.POLAR_GUIDE_COLOR),
            "width": model.POLAR_SPOKE_WIDTH, "items": items}


# -------------------------------------------------------------------
# typography
# -------------------------------------------------------------------
@lru_cache(maxsize=None)
def _font_name(stack, bold, italic):
    for name in stack:
        fp = FontProperties(family=[name], weight="bold" if bold else "normal",
                            style="italic" if italic else "normal")
        try:
            findfont(fp, fallback_to_default=False)
        except ValueError:
            continue
        return name
    return FontProperties().get_name()


def _resolved_fp(stack, *, bold, italic, size):
    family = _font_name(tuple(stack), bold, italic)
    return FontProperties(family=[family], size=size,
                          weight="bold" if bold else "normal",
                          style="italic" if italic else "normal")


def _label_size(span):
    size = span["size"]
    if size >= 17:
        return 18.0
    if size >= 10:
        return 11.0
    if size < 6.5:
        return model.LEGEND_SIZE
    return model.LABEL_SIZE


def _editor_point(x, y):
    left, top = _EDITOR_POSITION
    return (left + (x - _EDITOR_SOURCE[0]) * _EDITOR_ZOOM,
            top + (y - _EDITOR_SOURCE[1]) * _EDITOR_ZOOM)


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
    fit_column = 0
    for s in spans:
        text = s["text"]
        if s["font"] in _SKIP_FONTS or not text.strip(" \xa0"):
            skipped += 1
            continue
        bold = s["font"] == "Arial-BoldMT"
        italic = s["font"] == "Menlo-Italic"
        stack = list(FONT_STACKS.get(s["font"], FONT_STACKS["ArialMT"]))
        if family_override and not s["font"].startswith("Menlo"):
            stack = [family_override] + [
                f for f in stack if f != family_override]
        editor = 268 <= s["x"] <= 376 and 23 <= s["y"] <= 81
        size = s["size"] * _EDITOR_ZOOM if editor else _label_size(s)
        fp = _resolved_fp(stack, bold=bold, italic=italic, size=size)
        asc = float(s.get("ascender", 0.8))
        color = s["color"]
        if color.startswith("#00") and len(color) == 9:
            color = "#" + color[3:]  # drop the leading zero byte
        if text == "model.py — Python":
            left, top = _EDITOR_POSITION
            fp.set_size(4.2)
            ax.text(left + _EDITOR_WIDTH * _EDITOR_ZOOM / 2,
                    PAGE_H - (top + _EDITOR_BAR_HEIGHT * _EDITOR_ZOOM / 2),
                    "model.py", fontproperties=fp, ha="center", va="center",
                    color=color, zorder=zorder)
        elif s["dir"][1] < -0.5:  # vertical, bottom-to-top, tops facing left
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
            x, y = _editor_point(s["x"], base_y) if editor else (s["x"], base_y)
            if editor and s["x"] >= 280 and s["y"] >= 73:
                # Place the fit tokens on one monospace grid with their
                # original colors and the original italic keyword argument.
                if fit_column == 0:
                    fit_x = x
                    fit_advance = TextToPath().get_text_width_height_descent(
                        "0", fp, ismath=False)[0]
                text = text.replace("ﬁ", "fi")
                columns = len(text)
                if text == "=" and fit_column == len("idata"):
                    fit_column += 1
                    columns += 1
                x = fit_x + fit_column * fit_advance
                fit_column += columns
            ax.text(x, PAGE_H - y, text, fontproperties=fp,
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
    editor_shapes = set(asset["regions"]["code_window"])
    editor_transform = (Affine2D()
                        .translate(-_EDITOR_SOURCE[0], _EDITOR_SOURCE[1] - PAGE_H)
                        .scale(_EDITOR_ZOOM)
                        .translate(_EDITOR_POSITION[0], PAGE_H - _EDITOR_POSITION[1]))
    translucent = 0
    retained = (drawings[str(i)] for i in range(len(drawings)))
    for i, ((draw_type, opacity, box, path), dr) in enumerate(
        zip(_reference_drawings(), retained, strict=True)
    ):
        if draw_type != dr["type"]:
            raise ValueError(f"Artwork and reference drawing {i} differ")
        if i not in _REPLACED_PATHS:
            if i in _SUMMARY_POLAR_RINGS or i in _SUMMARY_POLAR_SPOKES:
                dr = _summary_polar_guide(i, dr)
            if i in _SUMMARY_AXIS_LABELS:
                _, x, _, _, _ = dr["items"][0]
                y0, y1 = _SUMMARY_AXIS_Y
                dr = {**dr, "items": [("l", x, y0, x, y1)]}
            if i in editor_shapes and dr["width"] is not None:
                dr = {**dr, "width": dr["width"] * _EDITOR_ZOOM}
            if i in _PANEL_STROKE_WIDTHS and dr["width"] is not None:
                dr = {**dr, "width": _PANEL_STROKE_WIDTHS[i]}
            draw_drawing(ax, dr, zorder=10 + i, opacity=opacity,
                         clip_box=box, clip_path=path,
                         joinstyle="round" if i in _WORKFLOW_ARROWHEADS else "miter",
                         transform=editor_transform if i in editor_shapes else None,
                         center_scale=_ESTIMATE_DOT_SCALE if i in _ESTIMATE_DOTS else 1.0)
            translucent += opacity < 1.0

    spans = [span for i, span in enumerate(asset["text_spans"])
             if i not in _SUMMARY_AXIS_LABEL_SPANS
             and not (span["x"] >= 250 and 23 <= span["y"] < 158
                      and not (268 <= span["x"] <= 376 and span["y"] <= 81))]
    drawn, skipped = draw_text(ax, spans, family_override=family)
    proportional = (family, *FONT_STACKS["ArialMT"]) if family else tuple(
        FONT_STACKS["ArialMT"])
    resolved_family = _font_name(proportional, False, False)
    y0, y1 = _SUMMARY_AXIS_Y
    for index, label in _SUMMARY_AXIS_LABELS.items():
        _, x, _, _, _ = drawings[str(index)]["items"][0]
        ax.text(x - 10.5, PAGE_H - (y0 + y1) / 2, label, rotation=90,
                ha="center", va="center", fontfamily=resolved_family,
                fontsize=model.LABEL_SIZE, color=model.INK, zorder=1000)
    model.draw_model_example(fig, (PAGE_W, PAGE_H), resolved_family)
    stats = {
        "reference_paths": len(drawings) - len(_REPLACED_PATHS),
        "reference_paths_with_group_alpha": translucent,
        "reference_text_drawn": drawn,
        "reference_text_skipped": skipped,
    }
    return fig, ax, stats


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
        if family and not name.startswith("Menlo"):
            stack = ([family]
                     + [f for f in stack if f != family])
        try:
            path = findfont(_resolved_fp(stack, bold=name == "Arial-BoldMT",
                                         italic=name == "Menlo-Italic", size=10),
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
    svg = svg_path.read_bytes()

    pdf_path = output_dir / f"{FIGURE_NAME}.pdf"
    fig.savefig(pdf_path, metadata={
        "CreationDate": datetime(2000, 1, 1, tzinfo=timezone.utc),
    })
    plt.close(fig)

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
        "illustrative_fit": {
            "reference_observed_degrees": model.REFERENCE_OBSERVED_DEG,
            "observed_offset_factor": model.OBSERVED_SPREAD,
            "observed_degrees": np.rad2deg(model.example()[1]).tolist(),
            "method": "random-walk Metropolis on the unconstrained scale",
            "prior_mean": model.PRIOR_MEAN,
            "prior_sd": model.PRIOR_SD,
            "seed": model.SEED,
            "warmup": model.BURN_IN,
            "posterior_draws_per_chain": model.N_DRAWS,
            "displayed_draws_per_chain": model.TRACE_DRAWS,
            "marked_draws_per_chain": model.MARKED_DRAWS_PER_CHAIN,
            "display_dots_per_observation": model.DISPLAY_DOTS_PER_ANGLE,
            "displayed_prior_predictive_draws": model.PRIOR_LINES,
            "displayed_posterior_predictive_draws": model.POSTERIOR_LINES,
            "chains": len(model.WRAP_BRANCHES),
            "initial_direction_turns": model.WRAP_BRANCHES,
            "direction_proposal_sd": model.MU_STEP,
            "concentration_updates_per_step": model.KAPPA_UPDATES,
            "concentration_proposal_sd": model.KAPPA_STEP,
        },
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
                        help="override the first family in proportional text "
                             "(CI uses 'Liberation Sans')")
    parser.add_argument("--font-dir", type=Path,
                        help="register local proportional fonts for this render")
    args = parser.parse_args(argv)
    if args.font_dir is not None:
        for path in findSystemFonts(fontpaths=[str(args.font_dir)]):
            fontManager.addfont(path)
        _font_name.cache_clear()
    prov = render(args.output_dir, family=args.font_family)
    for name, h in prov["exports"].items():
        print(f"{name}  {h[:16]}…")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
