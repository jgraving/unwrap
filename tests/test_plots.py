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


def test_repo_rel_for_files_and_dirs():
    import os

    base = os.path.realpath("/tmp")
    d = os.path.join(base, "figures", "bats")
    os.makedirs(d, exist_ok=True)
    assert plots.repo_rel(d) == "figures/bats"
    f = os.path.join(d, "kappa_hat.pdf")
    assert plots.repo_rel(f) == "figures/bats/kappa_hat.pdf"
