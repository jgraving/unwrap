"""Offline checks for the reference artwork and generated model panel."""
import hashlib
import io
import json
import struct
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib
import numpy as np
import pytest
from scipy.integrate import trapezoid

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
    pymupdf = pytest.importorskip("pymupdf")

    doc = pymupdf.open(ga.REFERENCE_PDF)
    page = doc.load_page(0)
    assert (page.rect.width, page.rect.height) == (599.0, 342.0)
    assert page.get_images(full=True) == []
    assert len(page.get_drawings()) == 526
    assert len(_pdf_spans(page)) == 236
    assert {s["font"] for s in _pdf_spans(page)} <= FONTS
    doc.close()


def test_svg_is_deterministic_vector_with_editable_text():
    pytest.importorskip("pymupdf")
    asset = ga.load_asset()

    def svg_bytes():
        fig, ax, _stats = ga.draw_figure(asset)
        original_paths = {patch.get_zorder() - 10: patch for patch in ax.patches}
        for index in ga._WORKFLOW_ARROWHEADS:
            assert original_paths[index].get_joinstyle() == "round"
        title, = [t for t in ax.texts if t.get_text() == "model.py"]
        assert title.get_ha() == "center"
        assert title.get_position()[0] == (
            ga._EDITOR_POSITION[0] + ga._EDITOR_WIDTH * ga._EDITOR_ZOOM / 2)
        assert all(t.get_fontsize() > 3.02 for t in ax.texts if t.get_text() == "import")
        polar = fig.axes[1]
        uniform, = [line for line in polar.lines if line.get_label() == "Uniform"]
        mean, = [line for line in polar.lines if line.get_label() == "Posterior mean"]
        assert uniform.get_linestyle() == "--"
        assert np.ptp(uniform.get_ydata()) == 0
        assert np.ptp(mean.get_ydata()) > 0
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        legend = fig.axes[2].get_legend()
        assert [t.get_text() for t in legend.get_texts()] == [
            "Prior", "Posterior", "Observed"]
        legend_bounds = legend.get_window_extent(renderer)
        predictor, = [t for t in fig.axes[3].texts
                      if t.get_text() == "Concentration predictor"]
        assert not legend_bounds.overlaps(predictor.get_window_extent(renderer))
        predictive_ring = fig.axes[1].spines["polar"]
        for index in ga._SUMMARY_POLAR_RINGS:
            ring = original_paths[index]
            assert ring.get_edgecolor() == predictive_ring.get_edgecolor()
            assert ring.get_linewidth() == predictive_ring.get_linewidth()
        for index, label in ga._SUMMARY_AXIS_LABELS.items():
            _, x, _, _, _ = asset["drawings"][str(index)]["items"][0]
            y0, y1 = ga._SUMMARY_AXIS_Y
            np.testing.assert_allclose(original_paths[index].get_path().vertices,
                                       [(x, ga.PAGE_H - y0), (x, ga.PAGE_H - y1)])
            label_position = (x - 10.5, ga.PAGE_H - (y0 + y1) / 2)
            caption, = [t for t in ax.texts if t.get_position() == label_position]
            assert caption.get_text() == label
            assert caption.get_rotation() == 90
            assert caption.get_fontsize() == ga.model.LABEL_SIZE
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
    assert len(root.findall(f".//{{{SVG_NS}}}path")) >= 377
    text_nodes = root.findall(f".//{{{SVG_NS}}}text")
    labels = {"".join(node.itertext()) for node in text_nodes}
    assert {"MCMC draws", "Direction predictor", "Concentration predictor",
            "Circular modulo", "Softplus", "Predictive density",
            "Observed", "Prior", "Posterior", "Likelihood parameters"} <= labels
    assert labels.isdisjoint({"Posterior draw", "Unconstrained", "Illustrative fit"})
    fit_start, = [node for node in text_nodes if node.text == "idata"]
    start_x = float(fit_start.get("x"))
    fit_tokens = [node for node in text_nodes
                  if node.get("y") == fit_start.get("y")
                  and float(node.get("x")) >= start_x]
    assert [node.text for node in fit_tokens] == [
        "idata", "=", "model", ".", "fit", "(", "draws", "=", "1000", ")"]
    offsets = np.array([float(node.get("x")) - start_x for node in fit_tokens])
    np.testing.assert_allclose(offsets, offsets[-1] / 28 * np.array(
        [0, 6, 8, 13, 14, 17, 18, 23, 24, 28]), atol=1e-6)
    styles = {node.text: node.get("style") for node in fit_tokens}
    assert "fill: #007197" in styles["fit"]
    assert "fill: #0087a8" in styles["draws"]
    assert "font-style: italic" in styles["draws"]
    assert "fill: #b15c00" in styles["1000"]
    assert b"stroke-opacity: 0.28" in a
    assert b"fill-opacity: 0.48" in a


def test_seeded_example_recomputes_normalized_predictive_densities():
    theta, observed, draws, prior, posterior = ga.model.example()
    ga.model.example.cache_clear()
    for previous, current in zip((theta, observed, draws, prior, posterior),
                                 ga.model.example(), strict=True):
        np.testing.assert_array_equal(previous, current)
    assert observed.size == 16
    display_angles, display_radii = ga.model._display_observations(observed)
    assert len(display_angles) == len(display_radii) == 64
    assert np.all(display_radii < 1)
    assert draws.shape == (3, ga.model.N_DRAWS, 2)
    traces = ga.model._display_traces(draws)
    assert traces.shape == (3, ga.model.TRACE_DRAWS, 2)
    original = draws[:, :ga.model.TRACE_DRAWS]
    np.testing.assert_array_equal(traces[:, :, 1], original[:, :, 1])
    np.testing.assert_allclose(ga.model.circ_mod(traces[:, :, 0]),
                               ga.model.circ_mod(original[:, :, 0]), atol=1e-12)
    wraps = np.floor((traces[:, :, 0] + np.pi) / (2 * np.pi))
    np.testing.assert_array_equal(wraps, np.broadcast_to(
        np.array(ga.model.WRAP_BRANCHES)[:, None], wraps.shape))
    mu, kappa = ga.model._parameters(draws)
    assert ((-np.pi <= mu) & (mu < np.pi)).all()
    assert (kappa > 0).all()
    np.testing.assert_allclose(mu, ga.model._parameters(draws + [2 * np.pi, 0.0])[0],
                               atol=1e-12)
    for samples in (prior, posterior):
        np.testing.assert_allclose(trapezoid(samples, x=theta), 1.0, atol=1e-10)
    assert posterior.mean(axis=0).max() > prior.mean(axis=0).max()


def test_render_writes_all_exports_with_consistent_provenance(tmp_path):
    pytest.importorskip("pymupdf")
    prov = ga.render(tmp_path)
    base = tmp_path / ga.FIGURE_NAME
    for ext in (".svg", ".pdf", ".png", ".html", ".json"):
        assert base.with_suffix(ext).exists()
    assert prov["page_pt"] == [599.0, 342.0]
    assert prov["inputs"]["artwork_sha256"] == ARTWORK_SHA
    assert prov["inputs"]["reference_pdf_sha256"] == REFERENCE_SHA
    assert prov["counts"]["reference_paths"] == 377
    assert prov["counts"]["reference_paths_with_group_alpha"] > 0
    assert prov["illustrative_fit"]["seed"] == ga.model.SEED
    assert prov["illustrative_fit"]["chains"] == 3
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
    pymupdf = pytest.importorskip("pymupdf")
    prov = ga.render(tmp_path)

    pdf = tmp_path / f"{ga.FIGURE_NAME}.pdf"
    doc = pymupdf.open(pdf)
    page = doc.load_page(0)
    assert (page.rect.width, page.rect.height) == (599.0, 342.0)
    assert page.get_images(full=True) == []
    assert doc.metadata["creationDate"] == "D:20000101000000Z"
    # The reference bands and individual polar traces retain their alpha.
    drawings = page.get_drawings()
    for key, target in (("fill_opacity", 0.25), ("stroke_opacity", 0.28)):
        assert any(abs(drawing[key] - target) < 1e-4 for drawing in drawings
                   if drawing.get(key) is not None)
    text = page.get_text()
    for needle in ("Experimental design", "Model summaries",
                   "Predictive density", "Softplus"):
        assert needle in text
    doc.close()
    assert prov["exports"]["pdf_sha256"] == ga.sha256(pdf)
