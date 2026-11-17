"""Scientific and file-integrity checks for the cached-fit figure command."""
import hashlib
from pathlib import Path
import sys
from types import SimpleNamespace
from unittest.mock import patch
import xml.etree.ElementTree as ET

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
import xarray as xr

sys.path.insert(0, str(Path(__file__).parents[1] / "notebooks"))
import beetles_figure as figure
import beetles_plots


@pytest.fixture
def posterior():
    shape = (2, 40)
    mu = np.linspace(-np.pi - 0.4, np.pi + 0.4, 80).reshape(shape)
    ds = xr.Dataset({
        "Intercept": (("chain", "draw"), mu),
        "Light": (("chain", "draw"), 0.8 + mu / 10),
        "kappa_Intercept": (("chain", "draw"), 2 + mu / 5),
        "kappa_Light": (("chain", "draw"), 1 - mu / 20),
        "1|Individual": (("chain", "draw", "individual"),
                         np.broadcast_to([0.2, -0.3, 1.4], (*shape, 3))),
        "Light|Individual": (("chain", "draw", "individual"),
                             np.broadcast_to([-0.1, 0.1, -1.2], (*shape, 3))),
    })
    return ds


def test_matches_existing_notebook_estimands(posterior):
    """Compare every draw with what the original plotting function receives."""
    captured = []
    with patch.object(beetles_plots, "_errorbar_mean_hdi", side_effect=lambda ax, x, s, *a, **k: captured.append(s)), \
         patch.object(beetles_plots, "_ensure_style"):
        old = beetles_plots.plot_concentration(
            SimpleNamespace(posterior=posterior), None,
            figure.circular_modulo, figure.softplus, save=False,
        )
    plt.close(old)
    arrays = figure.posterior_arrays(posterior)
    expected = [*arrays["within"], *arrays["between"], *arrays["mass"], arrays["difference"]]
    assert len(captured) == len(expected) == 7
    for before, after in zip(captured, expected):
        np.testing.assert_allclose(before, after, atol=1e-14)


def test_named_dimensions_and_circular_wrap(posterior):
    arrays = figure.posterior_arrays(posterior)
    transposed = figure.posterior_arrays(posterior.transpose("individual", "draw", "chain", missing_dims="ignore"))
    shifted = posterior.copy(deep=True)
    shifted["Intercept"] += 2 * np.pi
    wrapped = figure.posterior_arrays(shifted)
    for key in arrays:
        np.testing.assert_allclose(arrays[key], transposed[key], atol=1e-14)
        np.testing.assert_allclose(arrays[key], wrapped[key], atol=1e-14)
    np.testing.assert_array_equal(arrays["difference"], arrays["mass"][1] - arrays["mass"][0])


def test_sector_center_and_independent_quadrature(posterior):
    """Check the recentered, seam-crossing sector without using the CDF helper.

    The fine grid is independent of the CDF implementation; the coarse grid
    matches the R notebook's existing 181-point trapezoid integration.
    """
    from scipy.integrate import trapezoid
    from scipy.stats import vonmises

    assert figure.STUDY["sector_center_degrees"] == 202.5
    assert figure.STUDY["sector_halfwidth_degrees"] == 45
    arrays = figure.posterior_arrays(posterior)
    kappa = np.stack([
        figure.softplus(posterior["kappa_Intercept"].values.ravel()),
        figure.softplus((posterior["kappa_Intercept"] + posterior["kappa_Light"]).values.ravel()),
    ])
    for points, atol in ((20001, 1e-8), (181, 3e-5)):
        theta = np.linspace(np.deg2rad(157.5), np.deg2rad(247.5), points)
        density = vonmises.pdf(theta[:, None, None], kappa, loc=arrays["mu"])
        mass = trapezoid(density, theta, axis=0)
        np.testing.assert_allclose(arrays["mass"], mass, atol=atol, rtol=0)


def test_negative_difference_interval_is_visible(posterior):
    shifted = posterior.copy(deep=True)
    shifted["Intercept"][:] = figure.FLOOD_DIR
    shifted["Light"][:] = np.pi
    arrays = figure.posterior_arrays(shifted)
    assert arrays["difference"].max() < 0
    fig = plt.figure()
    ax = figure.draw_summary(fig, 0, arrays["difference"], "Change", difference=True)
    assert ax.get_ylim()[0] <= arrays["difference"].min()
    plt.close(fig)
    with patch.object(beetles_plots, "_ensure_style"):
        fig = beetles_plots.plot_concentration(
            SimpleNamespace(posterior=shifted), None,
            figure.circular_modulo, figure.softplus, save=False,
        )
    assert fig.axes[-1].get_ylim()[0] <= min(figure._mean_hdi(arrays["difference"], .89))
    plt.close(fig)


def test_uniform_distribution_and_identical_conditions(posterior):
    uniform = posterior.copy(deep=True)
    uniform["kappa_Intercept"][:] = -1000  # softplus underflows to zero
    uniform["kappa_Light"][:] = 0
    arrays = figure.posterior_arrays(uniform)
    np.testing.assert_allclose(arrays["mass"], 0.25, atol=1e-14)
    np.testing.assert_array_equal(arrays["within"], 0)
    np.testing.assert_allclose(arrays["difference"], 0, atol=1e-14)
    same = posterior.copy(deep=True)
    for key in ("Light", "kappa_Light", "Light|Individual"):
        same[key][:] = 0
    arrays = figure.posterior_arrays(same)
    for key in ("mu", "within", "between", "mass"):
        np.testing.assert_array_equal(arrays[key][0], arrays[key][1])


def test_rejects_incomplete_or_nonfinite_posterior(posterior):
    with pytest.raises(ValueError, match="missing"):
        figure.posterior_arrays(posterior.drop_vars("Light"))
    posterior["Light"][0, 0] = np.nan
    with pytest.raises(ValueError, match="Non-finite"):
        figure.posterior_arrays(posterior)


def test_checksum_failure_never_overwrites_data(tmp_path, monkeypatch):
    path = tmp_path / "data.txt"
    path.write_bytes(b"local modification")
    monkeypatch.setattr(figure, "INPUT_HASHES", {"data.txt": hashlib.sha256(b"expected").hexdigest()})
    with patch.object(figure.urllib.request, "urlopen") as download:
        with pytest.raises(ValueError, match="Checksum mismatch"):
            figure.prepare_inputs(tmp_path, download=True)
        download.assert_not_called()
    assert path.read_bytes() == b"local modification"
    path.unlink()
    with pytest.raises(FileNotFoundError, match="--download"):
        figure.prepare_inputs(tmp_path)


def test_missing_font_is_explicit():
    with pytest.raises(ValueError, match="unavailable"):
        figure.resolve_fonts(family="AnIntentionallyMissingFigureFont")


def test_vector_export_and_reproducible_svg(posterior, tmp_path):
    pytest.importorskip("cairosvg")
    data = pd.DataFrame({"Individual": [1, 2, 3] * 2, "Light": [0] * 3 + [1] * 3,
                         "HeadingRadians": [0.1, 2, -2, -2.1, -2.3, -2.4]})
    arrays = figure.posterior_arrays(posterior)
    svg = figure.build_svg(arrays, data, family="DejaVu Sans")
    assert svg == figure.build_svg(arrays, data, family="DejaVu Sans")
    root = ET.fromstring(svg)
    assert not root.findall(".//{http://www.w3.org/2000/svg}image")
    assert not root.findall(".//{http://www.w3.org/2000/svg}text")
    assert root.find("{http://www.w3.org/2000/svg}desc").text == figure.CAPTION
    stem = figure.write_outputs(tmp_path, arrays, data, {}, family="DejaVu Sans")
    assert stem.with_suffix(".pdf").read_bytes().startswith(b"%PDF-")
    assert stem.with_suffix(".png").read_bytes().startswith(b"\x89PNG")
    assert stem.with_suffix(".json").exists()
    assert "observed circular mean" in stem.with_suffix(".html").read_text()
