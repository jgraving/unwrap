"""Figure plumbing for the research notebooks."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

import unwrap.plots as plots


def test_set_style_two_resolution_tiers():
    rc = plots.set_style("/tmp/figures_test/x", extra_rc={"font.size": 8})
    # inline embed vs written file: low vs publication
    assert rc["figure.dpi"] == plots.DEFAULT_DISPLAY_DPI == 100
    assert rc["savefig.dpi"] == plots.PUBLICATION_DPI == 300
    # vector PDFs stay crisp at any zoom
    assert rc["pdf.fonttype"] == 42
    plt.rcdefaults()


def test_save_figure_writes_vector_pdf(tmp_path):
    fig, ax = plt.subplots()
    ax.plot([0.0, 1.0], [0.0, 1.0])
    out = plots.save_figure(fig, str(tmp_path / "figures" / "study"), "f", pdf=True)
    assert len(out) == 1
    assert out[0].name == "f.pdf"
    assert out[0].exists() and out[0].stat().st_size > 200


def test_save_figure_pdfs_and_svg_together(tmp_path):
    fig, ax = plt.subplots()
    ax.plot([0.0, 1.0], [0.0, 1.0])
    out = plots.save_figure(fig, str(tmp_path / "figures" / "study"), "f",
                            pdf=True, svg=True)
    assert [p.name for p in out] == ["f.pdf", "f.svg"]
    assert all(p.exists() and p.stat().st_size > 200 for p in out)


def test_repo_rel_for_files_and_dirs():
    import os

    base = os.path.realpath("/tmp")
    d = os.path.join(base, "figures", "bats")
    os.makedirs(d, exist_ok=True)
    assert plots.repo_rel(d) == "figures/bats"
    f = os.path.join(d, "kappa_hat.pdf")
    assert plots.repo_rel(f) == "figures/bats/kappa_hat.pdf"


def test_save_by_path_routes_by_suffix(tmp_path):
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    fig_dir = str(tmp_path / "figures" / "study")
    out_pdf = plots.save_by_path(fig, fig_dir, "figures/study/x.pdf")
    out_svg = plots.save_by_path(fig, fig_dir, "figures/study/y.svg")
    plt.close(fig)
    pdf = tmp_path / "figures" / "study" / "x.pdf"
    svg = tmp_path / "figures" / "study" / "y.svg"
    assert [str(p) for p in out_pdf] == [str(pdf)]
    assert [str(p) for p in out_svg] == [str(svg)]
    assert pdf.exists() and pdf.read_bytes()[:5] == b"%PDF-"
    assert svg.exists() and svg.read_text()[:5] == "<?xml"
    # nothing else in the study dir
    assert sorted(p.name for p in (tmp_path / "figures" / "study").iterdir()) == ["x.pdf", "y.svg"]
