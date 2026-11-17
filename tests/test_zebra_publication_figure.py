"""Hermetic checks for the zebra publication figure generator.

Everything renders on synthetic fits and tiny raster fixtures with
DejaVu Sans, so no ``unwrap-data`` checkout or licensed fonts are needed.
The layout test compares the measured layout constants against the
committed reference PDF, and the restore test verifies that the original
JPEG streams survive the SVG/PDF round trip.
"""
import hashlib
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
import xarray as xr

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "notebooks"))
import zebra_publication_figure as figure  # noqa: E402
import zebra_posthoc as posthog  # noqa: E402

SVG_NS = "http://www.w3.org/2000/svg"
TRACKS = [f"observation088-{i}" for i in range(4)]
FACTORS = {name: ("chain", "draw", "Track__factor_dim")
           for name in ("1|Track", "PostStimulus|Track",
                        "kappa_1|Track", "kappa_PostStimulus|Track")}


def make_idata(seed=2718, n_chains=2, n_draws=24):
    rng = np.random.default_rng(seed)
    post = {name: xr.DataArray(rng.normal(0.0, 0.2, (n_chains, n_draws)),
                               dims=("chain", "draw"))
            for name in ("Intercept", "PostStimulus",
                         "kappa_Intercept", "kappa_PostStimulus")}
    for name, dims in FACTORS.items():
        post[name] = xr.DataArray(
            rng.normal(0.0, 0.2, (n_chains, n_draws, len(TRACKS))),
            dims=dims,
            coords={"Track__factor_dim": TRACKS})
    ds = xr.Dataset(post)
    # the statistics path only reads ``idata.posterior``; a plain
    # namespace keeps the fixture off arviz's constructor surface, which
    # shifted between 0.x and 1.x
    import types
    return types.SimpleNamespace(posterior=ds)


def make_df_blocks(seed=71, heading_offset=None):
    rng = np.random.default_rng(seed)
    rows = []
    for t in TRACKS:
        for cond in (0, 1):
            base = 1.0 + cond * 2.2
            for _ in range(4):
                rows.append({"Track": t, "Observation": "observation088",
                             "Condition": "pre-approach" if cond == 0
                             else "post-approach",
                             "HeadingRadians": (rng.normal(base, 0.25)
                                                + np.pi) % (2 * np.pi) - np.pi,
                             "PostStimulus": cond})
    df = pd.DataFrame(rows)
    if heading_offset is not None:
        df["HeadingRadians"] = df["HeadingRadians"] + heading_offset
    return df


@pytest.fixture
def idata():
    return make_idata()


@pytest.fixture
def df_blocks():
    return make_df_blocks()


def small_jpeg(path, w=60, h=40, seed=3):
    from PIL import Image
    rng = np.random.default_rng(seed)
    a = rng.integers(0, 255, (h, w, 3), dtype="uint8")
    Image.fromarray(a).save(path)
    return Path(path)


def make_context(tmp_path, idata, df_blocks):
    """A synthetic drawing context with the same shape as build_context's.

    Tiny mosaic stand-in, a wandering focal track plus two others, and one
    JPEG per raster slot, so draw_figure exercises every panel without
    unwrap-data.
    """
    rng = np.random.default_rng(11)
    img = rng.integers(0, 255, (400, 160, 4), dtype="uint8")
    extent = (267525.0, 267680.0, 37300.0, 37770.0)

    frames = np.arange(0, 901, 2)
    rows = []
    for k, (name, fx, fy) in enumerate([
            (figure.FOCAL_TRACK,
             lambda f: 267580 + 40 * np.sin(f / 40.0),
             lambda f: 37535 + 120 * np.cos(f / 60.0)),
            (f"{figure.OBS}-1",
             lambda f: 267600 + 25 * np.cos(f / 55.0),
             lambda f: 37510 + 90 * np.sin(f / 70.0)),
            (f"{figure.OBS}-2",
             lambda f: 267620 + 30 * np.sin(f / 51.0 + 1.0),
             lambda f: 37550 - 80 * np.cos(f / 80.0))]):
        for f in frames:
            rows.append({"frame": int(f), "track_number": name,
                         "kp0_x": float(fx(f)), "kp0_y": float(fy(f))})
    df_pos = pd.DataFrame(rows)

    bear_rows = []
    for f in frames:
        bear_rows.append({"observation": figure.OBS,
                          "track_number": figure.FOCAL_TRACK,
                          "frame": int(f),
                          "compass_bearing": float((90.0 + 20 * np.sin(f / 40.0)) % 360.0),
                          "pre_scare": bool(f < 400),
                          "row_index": int(f)})
    df_bearings = pd.DataFrame(bear_rows)

    photo = small_jpeg(tmp_path / "photo.jpg", 120, 80, seed=4)
    insets = {key: small_jpeg(tmp_path / f"inset_{key}.jpg", 60, 60,
                              seed=5 + i)
              for i, key in enumerate(("pre", "post"))}
    anchors = {
        "pre": {"frame": 100, "pos": [267580.0, 37585.0],
                "arrow_bearing_deg": 93.2, "n_window": 30},
        "post": {"frame": 800, "pos": [267610.0, 37480.0],
                 "arrow_bearing_deg": 214.1, "n_window": 30},
    }
    numbers = figure.summary_numbers(idata, df_blocks)
    return {"img": img, "extent": extent, "df_pos": df_pos,
            "df_bearings": df_bearings, "numbers": numbers,
            "anchors": anchors, "photo": photo, "inset_assets": insets,
            "fit_meta": {}, "divergences": 0}


# ---------------------------------------------------------------------------
# Statistics path (summary_numbers on the saved-fit structure, synthetic)
# ---------------------------------------------------------------------------
def test_summary_statistics_shapes(idata, df_blocks):
    stats = figure.summary_numbers(idata, df_blocks)
    prep = stats["stats"]
    assert stats["nT"] == len(TRACKS)
    assert stats["tracks"] == list(TRACKS)
    assert prep["theta"].shape == (200,)
    assert prep["r_group_pre"].shape == (200,)
    assert prep["r_group_post"].shape == (200,)
    assert prep["r_ind_pre"].shape == (len(TRACKS), 200)
    assert prep["r_ind_post"].shape == (len(TRACKS), 200)
    # unit-scaled KDE radii: 1 at zero density, up to 2 at the max
    for key in ("r_group_pre", "r_group_post"):
        assert prep[key].min() >= 1.0 - 1e-12
        assert prep[key].max() <= 2.0 + 1e-12
    # the radial scale is shared across all group and individual KDEs
    maxima = [prep[key].max() for key in ("r_group_pre", "r_group_post",
                                          "r_ind_pre", "r_ind_post")]
    assert max(maxima) == pytest.approx(2.0, abs=1e-6)
    for cond in ("pre", "post"):
        lo = stats["group_alignment"][f"R_group_{cond}_lo"]
        hi = stats["group_alignment"][f"R_group_{cond}_hi"]
        mean = stats["group_alignment"][f"R_group_{cond}_mean"]
        assert 0.0 <= lo <= hi <= 1.0
        assert 0.0 <= mean <= 1.0
    assert len(stats["individual_concentration"]["R_hdi_post"]) == len(TRACKS)
    for hdi in (stats["individual_concentration"]["R_hdi_pre"],
                stats["individual_concentration"]["R_hdi_post"]):
        assert all(lo <= hi for lo, hi in hdi)
    for values in (stats["empirical_mean_direction_deg"]["pre"],
                   stats["empirical_mean_direction_deg"]["post"]):
        assert len(values) == len(TRACKS)


def _assert_direction_lists_equal(a, b):
    assert len(a) == len(b)
    for v_a, v_b in zip(a, b):
        if v_a is None or v_b is None:
            assert v_a is None and v_b is None
        else:
            assert v_a == pytest.approx(v_b, abs=1e-9)


def test_summary_statistics_cross_check(idata, df_blocks):
    """Group alignment and mean direction agree with direct computation."""
    stats = figure.summary_numbers(idata, df_blocks)
    directions = posthog.posterior_directions(idata, df_blocks)
    mu, idx = directions["mu"], directions["idx_pre"]
    n_samples = idata.posterior.sizes["chain"] * idata.posterior.sizes["draw"]
    assert mu.shape == (n_samples, 2 * len(TRACKS))
    assert float(np.nanmax(np.abs(mu))) <= np.pi + 1e-12
    r_pre = np.abs(np.exp(1j * mu[:, idx]).mean(axis=1))
    assert stats["stats"]["R_group_pre_mean"] == pytest.approx(
        float(r_pre.mean()), abs=1e-12)
    from scipy.stats import circmean
    pre_draw = circmean(mu[:, idx], high=np.pi, low=-np.pi, axis=1)
    group_mean = circmean(pre_draw, high=np.pi, low=-np.pi)
    expected_axis = float(np.rad2deg(group_mean) % 360.0)
    assert stats["group_mean_direction_deg"]["pre"] == pytest.approx(
        expected_axis, abs=1e-9)


def test_summary_statistics_wrap_invariance(df_blocks):
    """Adding 2*pi to every heading leaves the empirical summary unchanged."""
    base = make_idata()
    a = figure.summary_numbers(base, df_blocks)
    b = figure.summary_numbers(base, make_df_blocks(heading_offset=2 * np.pi))
    for cond in ("pre", "post"):
        _assert_direction_lists_equal(
            a["empirical_mean_direction_deg"][cond],
            b["empirical_mean_direction_deg"][cond])
    assert a["paired_permanent_turn_deg"] == pytest.approx(
        b["paired_permanent_turn_deg"], abs=1e-9)
    assert a["group_alignment"] == pytest.approx(b["group_alignment"])


def test_paired_turn_is_circular(idata, df_blocks):
    """The paired turn is the circular mean of the per-track turns."""
    stats = figure.summary_numbers(idata, df_blocks)
    turns = []
    for m_pre, m_post in zip(stats["stats"]["emp_mu_pre"],
                             stats["stats"]["emp_mu_post"]):
        if not (np.isnan(m_pre) or np.isnan(m_post)):
            turns.append(np.angle(np.exp(1j * (m_post - m_pre))))
    from unwrap.circular import circ_mean
    expected = float(np.rad2deg(circ_mean(np.asarray(turns))) % 360.0)
    assert stats["paired_permanent_turn_deg"] == pytest.approx(expected,
                                                               abs=1e-9)


# ---------------------------------------------------------------------------
# Vector export (synthetic context, DejaVu Sans)
# ---------------------------------------------------------------------------
def test_vector_export_deterministic_and_text_free(tmp_path, idata,
                                                   df_blocks):
    ctx = make_context(tmp_path, idata, df_blocks)
    svg = figure.draw_figure(ctx, family="DejaVu Sans")
    assert svg == figure.draw_figure(ctx, family="DejaVu Sans")
    root = ET.fromstring(svg)
    # mosaic, photo and the two inset crops: four raster embeds, all text
    # outlined to paths
    assert len(root.findall(f".//{{{SVG_NS}}}image")) == 4
    assert not root.findall(f".//{{{SVG_NS}}}text")
    desc = root.find(f"{{{SVG_NS}}}desc")
    assert desc is not None and desc.text == figure.CAPTION
    title = root.find(f"{{{SVG_NS}}}title")
    assert title is not None and title.text == "Zebra escape response"
    assert root.get("role") == "img"


def test_write_outputs_round_trip(tmp_path, idata, df_blocks):
    from PIL import Image
    ctx = make_context(tmp_path, idata, df_blocks)
    svg = figure.draw_figure(ctx, family="DejaVu Sans")
    stem = figure.write_outputs(tmp_path / "out", ctx, svg,
                                {"figure": "test"}, family="DejaVu Sans")
    assert stem.with_suffix(".pdf").read_bytes().startswith(b"%PDF-")
    assert stem.with_suffix(".png").read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    svg_text = stem.with_suffix(".svg").read_text()
    assert "<svg " in svg_text and svg_text.rstrip().endswith("</svg>")
    assert '<title id="figure-title">Zebra escape response</title>' in svg_text
    assert stem.with_suffix(".html").read_text().startswith("<!doctype html>")
    meta = json.loads(stem.with_suffix(".json").read_text())
    for suffix in (".svg", ".pdf", ".png", ".html"):
        assert meta["outputs_sha256"][suffix] == figure.sha256(
            stem.with_suffix(suffix))
    with Image.open(stem.with_suffix(".png")) as img:
        assert img.size == (4 * 900, 4 * 1080)


def test_restore_pdf_source_images_keeps_jpeg_streams(tmp_path):
    import cairosvg
    import fitz

    jpeg = small_jpeg(tmp_path / "src.jpg", 72, 48, seed=9)
    rect = (10.0, 10.0, 70.0, 50.0)  # page points, y down
    fig = plt.figure(figsize=(1.0, 1.0))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 72)
    ax.set_ylim(72, 0)
    ax.axis("off")
    ax.imshow(plt.imread(str(jpeg)), extent=(rect[0], rect[2],
                                             rect[3], rect[1]),
              aspect="auto")
    svg_path = tmp_path / "fig.svg"
    fig.savefig(svg_path, format="svg",
                metadata={"Date": None})
    plt.close(fig)

    pdf_path = tmp_path / "fig.pdf"
    cairosvg.svg2pdf(bytestring=svg_path.read_bytes(), write_to=str(pdf_path))

    def decode_first_image(path):
        with fitz.open(path) as doc:
            row = doc[0].get_images(full=True)[0]
            filter_ = doc.xref_get_key(row[0], "Filter")
            pix = fitz.Pixmap(doc, row[0])
            data = bytes(pix.samples)
            info = doc[0].get_image_info()[0]
            bbox = [round(v, 2) for v in info["bbox"]]
        return filter_, row[2], row[3], data, bbox

    filter_before, w_before, h_before, _, bbox_before = \
        decode_first_image(pdf_path)
    assert filter_before[1] != "/DCTDecode"  # cairosvg re-encodes as flate
    # cairosvg resampled the embed to display resolution, not 72 x 48
    assert (w_before, h_before) != (72, 48)

    figure.restore_pdf_source_images(pdf_path, [(jpeg, rect)])

    filter_after, w_after, h_after, data_after, bbox_after = \
        decode_first_image(pdf_path)
    assert filter_after == ("name", "/DCTDecode")
    assert (w_after, h_after) == (72, 48)
    # decoding the restored object matches a fresh embed-decode of the
    # same source file (same stream, same decoder -> identical pixels)
    with fitz.open() as src_doc:
        src_page = src_doc.new_page()
        src_page.insert_image(fitz.Rect(0, 0, 72, 48), filename=str(jpeg),
                              keep_proportion=False)
        src_row = src_page.get_images(full=True)[0]
        src_pix = fitz.Pixmap(src_doc, src_row[0])
        src_data = bytes(src_pix.samples)
    assert data_after == src_data
    # placement is untouched (get_image_info's bbox carries a sub-point
    # transform overshoot, so compare against the placed rect with slack)
    placed = fitz.Rect(*rect)
    assert bbox_after == bbox_before
    assert all(abs(a - b) <= 0.5 for a, b in zip(bbox_after,
                                                 (placed.x0, placed.y0,
                                                  placed.x1, placed.y1)))


def test_restore_pdf_source_images_refuses_ambiguous_placement(tmp_path):
    """Two exported rasters near one source rect must raise."""
    import cairosvg

    # The 1.0-inch figure is a 72 x 72 point page and each axes covers
    # the whole of it in page coordinates (no scaling, so the image
    # placement rectangles are exact). The two images must live on two
    # separate axes: this matplotlib's SVG backend emits at most one
    # <image> element per axes in a fresh figure (with or without
    # distinct content), while the real composite runs a canvas draw
    # before adding its overlay images and is exercised by
    # test_write_outputs_round_trip instead. The 72 x 48 JPEGs match
    # the 60 x 40 pt boxes (aspect 1.5), so cairosvg does not
    # letterbox the embeds.
    fig = plt.figure(figsize=(1.0, 1.0))
    ax_a = fig.add_axes([0, 0, 1, 1], zorder=1)
    ax_b = fig.add_axes([0, 0, 1, 1], zorder=2)
    for ax in (ax_a, ax_b):
        ax.set_xlim(0, 72)
        ax.set_ylim(72, 0)
        ax.set_facecolor("none")
        ax.axis("off")
    jpeg = small_jpeg(tmp_path / "src.jpg", 72, 48)
    rect = (10.0, 10.0, 70.0, 50.0)
    # One raster exactly on the source rect; the other's centre is 2.1
    # pt away -- inside the 5 pt destination tolerance, yet far enough
    # to keep the 1 pt scratch embedding unambiguous
    ax_a.imshow(plt.imread(str(small_jpeg(tmp_path / "a.jpg", 72, 48,
                                          seed=20))),
                extent=(rect[0], rect[2], rect[3], rect[1]),
                aspect="auto")
    ax_b.imshow(plt.imread(str(small_jpeg(tmp_path / "b.jpg", 72, 48,
                                          seed=21))),
                extent=(rect[0] + 1.5, rect[2] + 1.5,
                        rect[3] + 1.5, rect[1] + 1.5),
                aspect="auto")
    svg_path = tmp_path / "fig.svg"
    fig.savefig(svg_path, format="svg", metadata={"Date": None})
    plt.close(fig)
    pdf_path = tmp_path / "fig.pdf"
    cairosvg.svg2pdf(bytestring=svg_path.read_bytes(), write_to=str(pdf_path))
    with pytest.raises(ValueError, match="exactly one exported raster"):
        figure.restore_pdf_source_images(pdf_path, [(jpeg, rect)])


def test_stabilize_pdf_id_pins_the_trailer_id():
    trailer = (b"/ID[<0123456789ABCDEF0123456789ABCDEF>"
               b"<FEDCBA9876543210FEDCBA9876543210>]>>")
    data = (b"%PDF-1.7\ntrailer\n<</Size 2/Root 1 0 R" + trailer
            + b"\nstartxref\n0\n%%EOF")
    out = figure._stabilize_pdf_id(data)
    assert len(out) == len(data)
    # the result does not depend on the ID MuPDF happened to stamp
    other = data.replace(b"FEDCBA9876543210FEDCBA9876543210",
                         b"00112233445566778899001122334455")
    assert figure._stabilize_pdf_id(other) == out
    # both slots carry the SHA-256 of the document with the ID zeroed,
    # truncated to the 16-byte slot width
    m = figure._PDF_ID_RE.search(out)
    digest = out[m.start(1):m.end(1)]
    base = (out[:m.start(1)]
            + b"0" * (m.end(2) - m.start(1))
            + out[m.end(2):])
    assert digest == (hashlib.sha256(base).hexdigest()[:32]
                      .upper().encode("ascii"))
    assert m.group(1) == m.group(2)
    for bad in (b"%PDF-1.7\nno trailer here",
                data + b"trailer/ID[<ABCDEF0123456789ABCDEF0123456789>"
                      b"<ABCDEF0123456789ABCDEF0123456789>]>"):
        with pytest.raises(ValueError, match="one trailer /ID"):
            figure._stabilize_pdf_id(bad)


# ---------------------------------------------------------------------------
# Pinned inputs and committed assets
# ---------------------------------------------------------------------------
def test_input_pins_cover_the_documented_inputs():
    assert figure.DATA_REVISION == "86c41b5f43eb61191405b9bc85baf730a492a06f"
    expected = (
        ["fit_cache/zebra_escape_fit.nc",
         "fit_cache/zebra_escape_fit.meta.json",
         "zebra_bearings/zebra_bearings_blocks.parquet",
         "zebra_bearings/mosaic/manifest.json"]
        + [f"zebra_bearings/mosaic/tile_{i:03d}.png" for i in range(8)]
        + [f"zebra_bearings/postures/observation088-{i}.parquet"
           for i in range(18)]
        + [f"zebra_bearings/raw/observation088-{i}.parquet"
           for i in range(18)])
    assert set(figure.INPUT_HASHES) == set(expected)
    assert all(len(h) == 64 for h in figure.INPUT_HASHES.values())


def test_prepare_inputs_verifies_and_tamper_proofs(tmp_path, monkeypatch):
    payload = b"pinned zebra input\n"
    monkeypatch.setattr(
        figure, "INPUT_HASHES",
        {"zebra_bearings/x.dat": hashlib.sha256(payload).hexdigest()})
    root = tmp_path / "data"
    (root / "zebra_bearings").mkdir(parents=True)
    (root / "zebra_bearings" / "x.dat").write_bytes(payload)
    assert figure.prepare_inputs(root) == root

    # existing files are verified, never overwritten
    (root / "zebra_bearings" / "x.dat").write_bytes(b"tampered")
    with pytest.raises(ValueError, match="Checksum mismatch"):
        figure.prepare_inputs(root)

    # missing files require an explicit download
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(FileNotFoundError):
        figure.prepare_inputs(empty, download=False)


def test_committed_assets_match_manifest():
    manifest = figure.check_assets()
    for name in ("drone_inset_pre.jpg", "drone_inset_post.jpg",
                 "zebra_photo.jpg"):
        assert name in manifest["images"]
    assert figure.sha256(figure.ASSETS / "reference.pdf") == \
        manifest["source_sha256"]


# ---------------------------------------------------------------------------
# Layout constants vs the committed reference PDF
# ---------------------------------------------------------------------------
def test_reference_layout_matches_constants():
    import fitz

    with fitz.open(figure.ASSETS / "reference.pdf") as doc:
        page = doc[0]
        assert page.rect == fitz.Rect(0, 0, figure.PAGE_W, figure.PAGE_H)
        rects = [r for row in page.get_images(full=True)
                 for r in page.get_image_rects(row[0])]
        # the mosaic fills panel a
        assert any(abs(r.x0 - figure.MAP_RECT[0]) < 1 and
                   abs(r.y0 - figure.MAP_RECT[1]) < 1 and
                   abs(r.x1 - figure.MAP_RECT[2]) < 1 and
                   abs(r.y1 - figure.MAP_RECT[3]) < 1
                   for r in rects)
        # the photo fills panel b
        assert any(abs(r.x0 - figure.PHOTO_RECT[0]) < 1 and
                   abs(r.y0 - figure.PHOTO_RECT[1]) < 1 and
                   abs(r.x1 - figure.PHOTO_RECT[2]) < 1 and
                   abs(r.y1 - figure.PHOTO_RECT[3]) < 1
                   for r in rects)
        # panel letters: 32 pt, at the measured positions
        spans = [s for b in page.get_text("dict")["blocks"]
                 if b.get("type") == 0
                 for l in b["lines"] for s in l["spans"]]
        letters = {s["text"].strip():
                   (s["bbox"], s["size"]) for s in spans
                   if s["text"].strip() in ("a.", "b.", "c.", "d.",
                                            "e.", "f.")}
        assert set(letters) == {"a.", "b.", "c.", "d.", "e.", "f."}
        for name, (bbox, size) in letters.items():
            letter = name.rstrip(".")
            character, x, y, _ = next((l for l in figure.LETTERS
                                       if l[0] == letter))
            assert character == letter
            assert size == pytest.approx(figure.FS["letter"], abs=0.3)
            assert abs(bbox[0] - x) < 1.5
            assert abs((bbox[1] + bbox[3]) / 2 - y) < 3.0
