#!/usr/bin/env python3
"""Replace the low-resolution a-e photo payloads with high-resolution rasters.

Source of the high-resolution pixels: the author's Keynote export of the
manuscript figure, ``locust_density_order_figure.pdf`` (pinned SHA-256
below). For every photo the SVG display box is kept identical: the ancestor
transforms, clip paths and every other attribute are untouched. Only the
base64 payload and the ``width``/``height`` attributes of the seven photo
``<image>`` nodes change, plus a per-image ``transform="scale(...)"`` equal
to the manuscript/stored pixel ratio that re-anchors the larger raster
inside the original box. Panel geometry, the calibrated 10 cm scale bar and
its backing label therefore stay exactly where they were, verified by
``tests/test_locust_density_figure.py``.

The rasters are resized to at most ``--target-dpi`` effective pixels per
display point (never upscaled beyond native), encoded as baseline JPEG at
``--quality`` without ICC profiles (the source profiles are sRGB-family;
measured per-channel bias before dropping them is below 0.5/255), and both
the checked-in ``.jpeg`` files and the in-SVG payloads are replaced with the
same bytes. ``images.json`` records the new dimensions and hashes alongside
the manuscript hashes, so ``upgrade_svg_payloads`` can re-apply the swap to
a freshly extracted SVG — which ``tools/refresh_locust_figure_artwork.py``
does after every artwork re-extraction.
"""
import argparse
import base64
import hashlib
import io
import json
import re
from pathlib import Path
import xml.etree.ElementTree as ET

from PIL import Image

SVG = "http://www.w3.org/2000/svg"
XLINK = "http://www.w3.org/1999/xlink"
ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "notebooks/assets/locust_density_order"
EXPORT_SHA = "4b953238e985098c1d87c82e7ce6c8f25de6374aefc528900a4ee1cd6e471bf2"
IMAGES = {11: "vr_arena", 12: "real_locust", 13: "virtual_locust_lateral",
          14: "virtual_locust_front", 15: "virtual_locust_rear",
          20: "real_swarm", 21: "virtual_swarm"}
# Vertical row ranges verified against the original low-resolution crops by
# brute-force alignment: the exported rasters carry a few boundary rows the
# manuscript crops omitted. Best-match mean-absolute differences were <= 0.7.
ROW_CROPS = {"virtual_locust_lateral": (2, 1079), "virtual_locust_rear": (1, 1079)}
IMAGE_LINE = re.compile(
    r'^(?P<head><ns0:image [^>]*?)ns1:href="data:image/(?P<kind>[a-z]+);base64,'
    r'(?P<b64>[A-Za-z0-9+/= ]+)"(?P<tail>[^>]*/>)$')
TRANSFORM = re.compile(r' transform="scale\(([^)]*)\)"')


def digest(data):
    if isinstance(data, Path):
        data = data.read_bytes()
    return hashlib.sha256(data).hexdigest()


def manifest():
    return json.loads((ASSETS / "images.json").read_text())


def upgrade_svg_payloads(svg_text):
    """Swap manuscript photo payloads for the checked-in high-resolution files.

    Deterministic and idempotent: each replacement photo payload equals the
    checked-in ``.jpeg`` bytes, with ``width``/``height`` from the manifest
    and ``transform="scale(manuscript/stored)"`` preserving the original
    display box. Already-upgraded lines pass through unchanged.
    """
    entries = manifest()["images"]
    by_manuscript = {item["manuscript_sha256"]: (name, item)
                     for name, item in entries.items() if "manuscript_sha256" in item}
    by_stored = {item["sha256"]: name for name, item in entries.items()}
    replaced = kept = 0
    lines = svg_text.splitlines()
    for index, line in enumerate(lines):
        match = IMAGE_LINE.match(line)
        if not match or match.group("kind") != "jpeg":
            continue
        payload = base64.b64decode(match.group("b64"))
        sha = digest(payload)
        if sha in by_stored:
            kept += 1
            continue
        entry = by_manuscript.get(sha)
        if entry is None:
            width = int(re.search(r'width="(\d+)"', match.group("head")).group(1))
            if width > 24:
                raise ValueError(f"unknown photo payload at line {index + 1}")
            continue  # Scale-label backing mask is not a photo.
        name, item = entry
        stored = (ASSETS / name).read_bytes()
        head = re.sub(r'width="\d+"', f'width="{item["width"]}"', match.group("head"), count=1)
        head = re.sub(r'height="\d+"', f'height="{item["height"]}"', head, count=1)
        scale_x = item["manuscript_width"] / item["width"]
        scale_y = item["manuscript_height"] / item["height"]
        head = f'{head.rstrip()} transform="scale({scale_x!r},{scale_y!r})"'
        b64 = base64.b64encode(stored).decode("ascii")
        lines[index] = f'{head} ns1:href="data:image/jpeg;base64,{b64}"{match.group("tail")}'
        replaced += 1
    if replaced + kept != len(by_stored):
        raise ValueError(f"expected {len(by_stored)} photo payloads, saw {replaced + kept}")
    return "\n".join(lines) + "\n"


def load_rasters(path):
    export = path.read_bytes()
    if digest(export) != EXPORT_SHA:
        raise SystemExit(f"export checksum {digest(export)} != pinned {EXPORT_SHA}")
    import fitz  # PyMuPDF; required only for this optional preparation command
    doc = fitz.open(stream=export, filetype="pdf")
    rasters = {}
    for xref, name in IMAGES.items():
        info = doc.extract_image(xref)
        image = Image.open(io.BytesIO(info["image"])).convert("RGB")
        if name in ROW_CROPS:
            top, bottom = ROW_CROPS[name]
            image = image.crop((0, top, image.width, bottom))
        rasters[name] = image
    return rasters


def photo_geometry(svg_path):
    """{name: (manuscript_width, manuscript_height, display_w_pt, display_h_pt)}."""
    root = ET.parse(svg_path).getroot()
    parents = {child: parent for parent in root.iter() for child in parent}
    entries = manifest()["images"]
    by_hash = {item.get("manuscript_sha256", item["sha256"]): name
               for name, item in entries.items()}
    result = {}
    for node in root.iter(f"{{{SVG}}}image"):
        payload = base64.b64decode(
            (node.get(f"{{{XLINK}}}href") or "").split("base64,", 1)[1].replace(" ", ""))
        name = by_hash.get(digest(payload))
        if name is None:
            continue  # Scale-bar backing alpha mask keeps its raster.
        sx = sy = 1.0
        ancestor = node
        while ancestor in parents:
            ancestor = parents[ancestor]
            values = re.fullmatch(
                r"matrix\(([^)]*)\)", ancestor.get("transform") or "")
            if values:
                a, b, c, d = (float(v) for v in values.group(1).split(",")[:4])
                sx *= (a * a + b * b) ** 0.5
                sy *= (c * c + d * d) ** 0.5
        width, height = int(node.get("width")), int(node.get("height"))
        result[name.removesuffix(".jpeg")] = (width, height, width * sx, height * sy)
    return result


def encode(raster, display_w_pt, target_dpi, quality):
    """JPEG pixels at <= target_dpi effective dpi, capped at native, exact aspect."""
    native_w, native_h = raster.size
    target_w = min(native_w, int(round(display_w_pt * target_dpi / 72.0)))
    target_h = min(native_h, int(round(target_w * native_h / native_w)))
    resized = raster.resize((target_w, target_h), Image.LANCZOS)
    buffer = io.BytesIO()
    resized.save(buffer, "JPEG", quality=quality, optimize=True)
    return buffer.getvalue(), (target_w, target_h)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--export-pdf", type=Path,
                        help="High-resolution Keynote export of the manuscript figure")
    parser.add_argument("--target-dpi", type=float, default=400.0,
                        help="Maximum effective pixels per display point (default 400)")
    parser.add_argument("--quality", type=int, default=88, help="JPEG quality (default 88)")
    parser.add_argument("--upgrade-only", action="store_true",
                        help="Re-apply the checked-in payloads to the asset SVG without re-encoding")
    args = parser.parse_args()

    svg_path = ASSETS / "original_panels.svg"
    manifest_path = ASSETS / "images.json"
    if args.upgrade_only:
        svg_path.write_text(upgrade_svg_payloads(svg_path.read_text()))
        entries = manifest()["images"]
        (ASSETS / "images.json").write_text(
            json.dumps(json.loads(manifest_path.read_text()), indent=2) + "\n")
        print(f"payloads re-applied for {len(entries)} images; "
              f"original_panels.svg sha256 {digest(svg_path)[:16]}...")
        return

    rasters = load_rasters(args.export_pdf)
    geometry = photo_geometry(svg_path)
    info = json.loads(manifest_path.read_text())
    for name, (old_w, old_h, disp_w, disp_h) in sorted(geometry.items()):
        data, (new_w, new_h) = encode(rasters[name], disp_w, args.target_dpi, args.quality)
        entry = info["images"][f"{name}.jpeg"]
        if "manuscript_sha256" not in entry:
            entry["manuscript_width"], entry["manuscript_height"] = old_w, old_h
            entry["manuscript_sha256"] = entry["sha256"]
        entry.update(width=new_w, height=new_h, sha256=digest(data))
        (ASSETS / f"{name}.jpeg").write_bytes(data)
        print(f"{name:24s} {old_w}x{old_h} -> {new_w}x{new_h} "
              f"({disp_w:.1f}x{disp_h:.1f}pt, "
              f"{new_w * 72 / disp_w:.0f} dpi, {len(data) / 1e3:.0f} KB)")

    manifest_path.write_text(json.dumps(info, indent=2) + "\n")
    svg_path.write_text(upgrade_svg_payloads(svg_path.read_text()))
    info["high_resolution_export"] = {
        "path": "locust_density_order_figure.pdf",
        "sha256": EXPORT_SHA,
        "creator": "Keynote, macOS 15.8.3 (24H241) Quartz PDFContext",
        "page_points": [750, 525],
        "policy": (f"rasters resized to at most {args.target_dpi:g} effective dpi "
                   "(never above native), "
                   f"JPEG quality {args.quality}, sRGB without ICC profile, rows "
                   "cropped to match the manuscript crops for the two virtual "
                   "views; per-image transform scale preserves each display box"),
    }
    info["original_panels"]["sha256"] = digest(svg_path)
    manifest_path.write_text(json.dumps(info, indent=2) + "\n")
    print(f"original_panels.svg sha256 {info['original_panels']['sha256'][:16]}...")


if __name__ == "__main__":
    main()
