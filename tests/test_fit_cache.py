"""Round-trip and validation tests for the netCDF load-or-save fit cache
(issue #130, bats part).

``fit_cache`` itself only needs stdlib + xarray; the snapshot format is a
netCDF-4 file, so these tests skip when no netCDF-4 engine (h5netcdf or
netCDF4) is installed in the CI dependency set.
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("xarray")
if not (importlib.util.find_spec("h5netcdf") or importlib.util.find_spec("netCDF4")):
    pytest.skip("no netCDF-4 save/load engine (h5netcdf or netCDF4) installed",
                allow_module_level=True)

_NOTEBOOKS = Path(__file__).resolve().parents[1] / "notebooks"
_spec = importlib.util.spec_from_file_location("fit_cache_under_test",
                                               _NOTEBOOKS / "fit_cache.py")
fc = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = fc
_spec.loader.exec_module(fc)


def _tiny_idata():
    """A minimal xarray InferenceData-shaped DataTree (all groups the
    bats snapshots carry, tiny draws, one string variable on purpose)."""
    import xarray

    rng = np.random.default_rng(0)
    post = xarray.Dataset(
        {
            "mu0": (("chain", "draw"), rng.normal(size=(2, 3))),
            "mu": (("chain", "draw", "mu_dim_0"),
                   rng.normal(size=(2, 3, 4))),
            "mu_treat": (("chain", "draw", "mu_treat_dim_0"),
                         rng.normal(size=(2, 3, 2))),
        },
        coords={"chain": [0, 1], "draw": [0, 1, 2]},
    )
    ss = xarray.Dataset(
        {"diverging": (("chain", "draw"), np.zeros((2, 3), dtype=bool)),
         "step_size": (("chain", "draw"), rng.uniform(size=(2, 3)))},
        coords={"chain": [0, 1], "draw": [0, 1, 2]},
    )
    od = xarray.Dataset(
        {"Heading_obs": (("obs",), rng.normal(size=4))},
        coords={"obs": [0, 1, 2, 3]},
    )
    cd = xarray.Dataset(
        {"tag": (("n",), np.array(["a", "b", "c"], dtype=object))},
        coords={"n": [0, 1, 2]},
    )
    return xarray.DataTree(children={
        "posterior": xarray.DataTree(name="posterior", dataset=post),
        "sample_stats": xarray.DataTree(name="sample_stats", dataset=ss),
        "constant_data": xarray.DataTree(name="constant_data", dataset=cd),
        "observed_data": xarray.DataTree(name="observed_data", dataset=od),
    })


def _cache_env(tmp_path, monkeypatch, *, enabled=True):
    monkeypatch.setenv(fc.ENV_DIR, str(tmp_path))
    if enabled:
        monkeypatch.setenv(fc.ENV_ENABLE, "1")
    else:
        monkeypatch.delenv(fc.ENV_ENABLE, raising=False)


def test_save_load_netcdf_roundtrip(tmp_path, monkeypatch):
    _cache_env(tmp_path, monkeypatch)
    meta = {"model": "tiny", "seed": 123}
    fc.save_snapshot("tiny", _tiny_idata(), meta=meta)
    assert (tmp_path / "tiny.nc").exists()
    out = fc.load_snapshot("tiny", expect=meta)
    src = _tiny_idata()
    for group in ("posterior", "sample_stats", "constant_data", "observed_data"):
        a, b = out[group], src[group]
        assert sorted(map(str, a.data_vars)) == sorted(map(str, b.data_vars))
        for v in a.data_vars:
            np.testing.assert_array_equal(a[v].values, b[v].values)
    # the string variable kept its values
    np.testing.assert_array_equal(
        out.constant_data["tag"].values.astype(str),
        src.constant_data["tag"].values.astype(str))


def test_sha_mismatch_is_loud(tmp_path, monkeypatch):
    import json

    _cache_env(tmp_path, monkeypatch)
    fc.save_snapshot("tiny", _tiny_idata(), meta={"seed": 123})
    meta_path = tmp_path / "tiny.meta.json"
    info = json.loads(meta_path.read_text())
    info["sha256"] = "0" * len(info["sha256"])
    meta_path.write_text(json.dumps(info))
    with pytest.raises(fc.CacheMismatch, match="sha256"):
        fc.load_snapshot("tiny")


def test_layout_mismatch_is_loud(tmp_path, monkeypatch):
    import json

    _cache_env(tmp_path, monkeypatch)
    fc.save_snapshot("tiny", _tiny_idata(), meta={"seed": 123})
    meta_path = tmp_path / "tiny.meta.json"
    info = json.loads(meta_path.read_text())
    info["posterior_vars"] = ["not-a-var"]
    meta_path.write_text(json.dumps(info))
    with pytest.raises(fc.CacheMismatch, match="posterior layout"):
        fc.load_snapshot("tiny")


def test_meta_expect_mismatch_is_loud(tmp_path, monkeypatch):
    _cache_env(tmp_path, monkeypatch)
    fc.save_snapshot("tiny", _tiny_idata(), meta={"seed": 123})
    with pytest.raises(fc.CacheMismatch, match="fit metadata mismatch"):
        fc.load_snapshot("tiny", expect={"seed": 456})


def test_load_or_fit_hit_and_miss(tmp_path, monkeypatch, capsys):
    _cache_env(tmp_path, monkeypatch)
    calls = {"n": 0}

    def fit_fn():
        calls["n"] += 1
        return _tiny_idata()

    out = fc.load_or_fit("tiny", fit_fn, meta={"seed": 123}, label="tiny fit")
    assert calls["n"] == 1
    assert "FIT-CACHE miss" in capsys.readouterr().out
    assert out is not None  # the full fit_fn return value
    out2 = fc.load_or_fit("tiny", fit_fn, meta={"seed": 123}, label="tiny fit")
    assert calls["n"] == 1  # cached, no refit
    logged = capsys.readouterr().out
    assert "FIT-CACHE hit" in logged and "tiny.nc" in logged
    np.testing.assert_array_equal(
        out2.posterior["mu0"].values, out.posterior["mu0"].values)


def test_disabled_is_a_noop(tmp_path, monkeypatch, capsys):
    _cache_env(tmp_path, monkeypatch, enabled=False)
    out = fc.load_or_fit("tiny", _tiny_idata, meta={"seed": 123})
    assert capsys.readouterr().out == ""
    assert not (tmp_path / "tiny.nc").exists()
    assert out is not None


def test_slim_snapshot_roundtrips(tmp_path, monkeypatch):
    _cache_env(tmp_path, monkeypatch)
    slim = fc.slim_snapshot(_tiny_idata(), drop_groups=("sample_stats",),
                            drop_posterior=("mu",))
    assert "sample_stats" not in slim.children
    assert "mu" not in slim.posterior.data_vars
    fc.save_snapshot("slim", slim, meta={"seed": 123, "stripped": ["mu"]})
    out = fc.load_snapshot("slim", expect={"seed": 123, "stripped": ["mu"]})
    assert "sample_stats" not in out.children
    assert "mu" not in out.posterior.data_vars
    np.testing.assert_array_equal(
        out.posterior["mu_treat"].values,
        _tiny_idata().posterior["mu_treat"].values)
