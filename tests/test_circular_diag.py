"""Circular r_hat/ESS in ``circular_summary`` (issue #117).

``circ_diag="circular"`` mirrors the R twins' ``diagnose_fit`` treatment
of angular variables: the ``Rhat_unwrap`` sine convention for ``r_hat``
and ESS on circular-mean-centred wrapped deviations. These tests pin the
Python side against that contract (including agreement with the R
implementation itself when Rscript is available) and keep the default
arviz-0.23 verbatim behaviour.
"""
from __future__ import annotations

import numpy as np
import pytest

from unwrap.circular import circular_summary, wrap_pi


def _idata(arr):
    import arviz as az
    return az.from_dict({"posterior": {"mu": np.asarray(arr, float)}})


def _seam_fixture():
    """Chains ONE circularly, read as two modes by the linear diagnostics.

    Chains 1-2 straddle the +-pi seam around the circle point pi; chains
    3-4 sit just inside the same cluster. Circularly a single cluster;
    linearly the straddling chains store a mean near 0 and read as a
    second mode far from the others.
    """
    rng = np.random.default_rng(11)
    c12 = wrap_pi(rng.vonmises(0.0, 25.0, (2, 1500)) + np.pi)
    c34 = wrap_pi(rng.vonmises(0.0, 200.0, (2, 1500)) + (np.pi - 0.3))
    return np.vstack([c12, c34])


def _ar_fixture():
    """Autocorrelated chains (HMC-like) diffusing through the +-pi seam."""
    rng = np.random.default_rng(4)
    z = np.empty((4, 1500))
    z[:, 0] = 0.05 * rng.normal(size=4)
    for t in range(1, z.shape[1]):
        z[:, t] = 0.9 * z[:, t - 1] + 0.05 * rng.normal(size=4)
    return wrap_pi(z + (np.pi - 0.02))


def test_circ_diag_circular_differs_from_linear_at_the_seam():
    idata = _idata(_seam_fixture())
    lin = circular_summary(idata, circ_var_names=["mu"])
    cir = circular_summary(idata, circ_var_names=["mu"],
                           circ_diag="circular")
    # Same draws: the diagnostics must genuinely differ, and the linear
    # r_hat reads the fake seam split while the circular one is calmer.
    assert lin.loc["mu", "r_hat"] > 1.01
    assert cir.loc["mu", "r_hat"] < lin.loc["mu", "r_hat"]
    assert (cir.loc["mu", "ess_bulk"] != lin.loc["mu", "ess_bulk"])


def test_circ_diag_is_rotation_invariant():
    arr = _seam_fixture()
    a = circular_summary(_idata(arr), circ_var_names=["mu"],
                         circ_diag="circular")
    b = circular_summary(_idata(wrap_pi(arr + 0.9)), circ_var_names=["mu"],
                         circ_diag="circular")
    # The centred deviations are identical under rotation, and the sine
    # max over the two projections is numerically origin-invariant.
    assert a.loc["mu", "ess_bulk"] == pytest.approx(b.loc["mu", "ess_bulk"],
                                                    abs=1e-2)
    assert a.loc["mu", "r_hat"] == pytest.approx(b.loc["mu", "r_hat"],
                                                 abs=2e-3)


def test_circ_diag_default_is_linear_and_validates():
    idata = _idata(_seam_fixture())
    default = circular_summary(idata, circ_var_names=["mu"])
    linear = circular_summary(idata, circ_var_names=["mu"],
                              circ_diag="linear")
    assert default.equals(linear)
    # Default r_hat is the plain arviz value on the raw angles.
    import arviz as az
    assert default.loc["mu", "r_hat"] == pytest.approx(
        round(float(az.rhat(_seam_fixture())), 3), abs=1e-3)
    with pytest.raises(ValueError, match="circ_diag"):
        circular_summary(idata, circ_var_names=["mu"], circ_diag="sine")


def _find_rscript():
    import shutil
    cand = shutil.which("Rscript")
    return cand if cand and __import__("pathlib").Path(cand).exists() else None


def test_circ_diag_matches_r_twin(tmp_path):
    """Python circular r_hat/ESS == the R twins' formulas on one fixture."""
    rscript = _find_rscript()
    if rscript is None:
        pytest.skip("Rscript not available")
    import pathlib
    import subprocess

    arr = _ar_fixture()  # [chain, draw] in Python
    x_iter_chain = np.asarray(arr).T  # [iteration, chain] for R
    csv = tmp_path / "draws.csv"
    np.savetxt(csv, x_iter_chain, delimiter=",")

    # Extract Rhat_unwrap verbatim from the shared R helpers.
    helpers = (pathlib.Path(__file__).resolve().parents[1] / "R" /
               "unwrap_helpers.R").read_text()
    start = helpers.find("Rhat_unwrap = function")
    assert start != -1, "Rhat_unwrap not found in R/unwrap_helpers.R"
    ob = helpers.find("{", start)
    depth, i = 0, ob
    while True:
        if helpers[i] == "{":
            depth += 1
        elif helpers[i] == "}":
            depth -= 1
            if depth == 0:
                break
        i += 1
    rhat_src = helpers[start:i + 1]
    (tmp_path / "rhat.R").write_text(
        "library(posterior)\n"
        "wrap <- function(x) (x + pi) %% (2 * pi) - pi\n"
        "circ_mean <- function(a) atan2(mean(sin(a)), mean(cos(a)))\n"
        + rhat_src + "\n")

    runner = tmp_path / "run.R"
    runner.write_text(
        "args <- commandArgs(trailingOnly = TRUE)\n"
        "source(args[1])\n"
        "x <- as.matrix(read.csv(args[2], header = FALSE))\n"
        "writeLines(sprintf('%.10f', Rhat_unwrap(x, method = 'sine'))[1])\n"
        "writeLines(sprintf('%.10f', "
        "posterior::ess_bulk(wrap(x - circ_mean(x))))[1])\n")
    proc = subprocess.run([rscript, str(runner), str(tmp_path / "rhat.R"),
                           str(csv)], cwd=str(tmp_path), capture_output=True,
                          text=True)
    if proc.returncode != 0:
        pytest.skip(f"Rscript failed: {proc.stderr.strip()[:200]}")
    r_rhat, r_ess = (float(v) for v in proc.stdout.split())

    got = circular_summary(_idata(arr), circ_var_names=["mu"],
                           circ_diag="circular")
    # Same formulas; arviz and posterior implement split-rhat and Geyer
    # ESS independently, so allow implementation noise, not semantics.
    assert got.loc["mu", "r_hat"] == pytest.approx(r_rhat, abs=5e-4)
    assert got.loc["mu", "ess_bulk"] == pytest.approx(r_ess, rel=0.25)
