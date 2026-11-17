"""Scientific and export checks for the saved-fit publication figure."""
import json
from pathlib import Path
import sys
from unittest.mock import patch
import xml.etree.ElementTree as ET

import numpy as np
import pytest
from scipy.stats import vonmises
import xarray as xr

sys.path.insert(0, str(Path(__file__).parents[1] / "notebooks"))
import locust_density_figure as figure


@pytest.fixture
def posterior():
    rng = np.random.default_rng(13)
    return xr.Dataset({prefix + term: (("chain", "draw"), rng.normal(scale=.3, size=(2, 30)))
                       for prefix in ("", "kappa_") for term in figure.TERMS})


def test_reference_matches_bambi_pymc_graph(posterior):
    """Independent graph check, with unequal means and nonzero interaction mean.

    No fitting: transport the stored-scale draws back to Bambi's centred
    PyMC coordinates, then evaluate each component's actual link.
    """
    import bambi as bmb
    import pandas as pd
    import pymc as pm
    from unwrap.link.bambi import LINKS

    order = np.array([0., .2, .4, .4, .8, 1., 1.])
    density = np.array([1., 4., 4., 8., 32., 64., 64.])
    df = pd.DataFrame({"y": np.arange(7) / 10, "t_order": 2*order-1,
                       "t_density": np.log2(density)/3-1})
    model = bmb.Model(bmb.Formula("y ~ t_density*t_order", "kappa ~ t_density*t_order"),
                      df, family="vonmises", center_predictors=True,
                      link={"mu": LINKS["circular_modulo"], "kappa": LINKS["softplus"]})
    model.build()
    stored = xr.DataTree.from_dict({"posterior": posterior})
    internal = model._re_center_intercept(stored)
    graph = pm.compute_deterministics(internal.posterior.to_dataset(), model=model.backend.model,
                                      var_names=["mu", "kappa"], progressbar=False)
    mu, kappa = figure.reference_parameters(posterior, order, density)
    np.testing.assert_allclose(mu, graph.mu.values.reshape(mu.shape), atol=1e-12)
    np.testing.assert_allclose(kappa, graph.kappa.values.reshape(kappa.shape), atol=1e-12)


def test_stable_density_and_intervals(posterior):
    orders = [0., .3, 1.]; densities = [1., 8., 64.]
    mu, kappa = figure.reference_parameters(posterior, orders, densities)
    expected = vonmises.logpdf(0, kappa, loc=mu)
    result = figure.reference_summary(posterior, orders, densities)
    np.testing.assert_allclose(result["mean"], np.exp(expected.mean(axis=0)), atol=1e-13)
    import arviz as az
    for i in range(3):
        lo, hi = az.hdi(expected[:, i], prob=.90)
        np.testing.assert_allclose([result["low"][i], result["high"][i]], np.exp([lo, hi]), atol=1e-13)
    wrapped = posterior.copy(deep=True); wrapped["Intercept"] += 2*np.pi
    for key in result:
        np.testing.assert_allclose(result[key], figure.reference_summary(wrapped, orders, densities)[key], atol=1e-13)


def test_uniform_and_dimension_order(posterior):
    uniform = posterior.copy(deep=True)
    for term in figure.TERMS:
        uniform["kappa_"+term][:] = -1000 if term == "Intercept" else 0
    result = figure.reference_summary(uniform, [0, 1], [1, 64])
    for values in result.values():
        np.testing.assert_allclose(values, 1/(2*np.pi))
    expected = figure.reference_parameters(posterior, [0, 1], [1, 64])
    actual = figure.reference_parameters(posterior.transpose("draw", "chain"), [0, 1], [1, 64])
    for x, y in zip(expected, actual): np.testing.assert_array_equal(x, y)
    with pytest.raises(ValueError, match="outside"):
        figure.reference_parameters(posterior, 0, 128)


def test_input_integrity_never_replaces_existing_file(tmp_path, monkeypatch):
    path = tmp_path / "example.nc"; path.write_bytes(b"modified")
    monkeypatch.setattr(figure, "INPUT_HASHES", {"example.nc": "expected"})
    with patch.object(figure.urllib.request, "urlopen") as download:
        with pytest.raises(ValueError, match="Checksum mismatch"):
            figure.prepare_inputs(tmp_path, download=True)
        download.assert_not_called()
    assert path.read_bytes() == b"modified"


def test_assets_selection_and_vector_export(posterior, tmp_path):
    pytest.importorskip("cairosvg")
    manifest = json.loads((figure.ASSETS / "images.json").read_text())
    assert len(manifest["images"]) == 7
    for name, item in manifest["images"].items(): assert figure.sha256(figure.ASSETS / name) == item["sha256"]
    selection = json.loads((figure.ASSETS / "selection.json").read_text())
    assert len(selection["models"]) == 5
    assert selection["models"][0]["delta"] == 0
    assert max(row["pareto_k_max"] for row in selection["models"]) < .7
    stats = figure.statistics(posterior)
    svg = figure.draw_figure(stats, selection, family="DejaVu Sans")
    assert svg == figure.draw_figure(stats, selection, family="DejaVu Sans")
    root = ET.fromstring(svg)
    # Seven source photos/renderings plus the original scale-label backing
    # and its alpha mask; statistical panels contain no raster images.
    assert len(root.findall(".//{http://www.w3.org/2000/svg}image")) == 9
    original = manifest["original_panels"]
    assert figure.sha256(figure.ASSETS / original["file"]) == original["sha256"]
    assert original["retained_panel_pixel_changes"] == 0
    assert original["scale_bar_label"] == "10 cm"
    ids = [node.get("id") for node in root.iter() if node.get("id")]
    assert len(ids) == len(set(ids))
    assert not root.findall(".//{http://www.w3.org/2000/svg}text")
    assert root.find("{http://www.w3.org/2000/svg}desc").text == figure.CAPTION
    stem = figure.write_outputs(tmp_path, stats, selection, {}, family="DejaVu Sans")
    assert stem.with_suffix(".pdf").read_bytes().startswith(b"%PDF-")
    assert stem.with_suffix(".png").read_bytes().startswith(b"\x89PNG")
