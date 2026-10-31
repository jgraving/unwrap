"""Small end-to-end checks; short chains test execution, not convergence."""
import arviz as az
import bambi as bmb
import numpy as np
from numpy.testing import assert_allclose
from pandas.testing import assert_frame_equal

from unwrap.circular import circular_summary
from unwrap.data import generate_circular_regression_data
from unwrap.link.bambi import LINKS
from unwrap.prior import CircularUniformPrior


def make_data():
    return generate_circular_regression_data(
        num_samples=30, x_range=(-1, 1), mean_intercept=0, mean_slope=1,
        conc_intercept=2, conc_slope=0.5, mean_link="circular_modulo",
        conc_link="softplus", seed=42,
    )


def test_generated_data_reproducible_and_on_correct_scale():
    data = make_data()
    assert_frame_equal(data, make_data())
    assert list(data.columns) == ["x", "y", "mean", "concentration"]
    assert len(data) == 30
    assert np.isfinite(data.to_numpy()).all()
    assert data.x.between(-1, 1).all()
    assert data.y.between(-np.pi, np.pi).all()
    assert_allclose(data["mean"], data.x)
    assert_allclose(data.concentration, np.logaddexp(0, 2 + 0.5 * data.x))


def test_bambi_fit_and_circular_summary():
    model = bmb.Model(
        bmb.Formula("y ~ x", "kappa ~ x"), make_data(), family="vonmises",
        link={"mu": LINKS["circular_modulo"], "kappa": LINKS["softplus"]},
        priors={
            "Intercept": CircularUniformPrior(),
            "x": bmb.Prior("Normal", mu=0, sigma=0.5 * np.pi),
            "kappa": {
                "Intercept": bmb.Prior("Normal", mu=2, sigma=0.5),
                "x": bmb.Prior("Normal", mu=0, sigma=1),
            },
        },
        center_predictors=False,
    )
    idata = model.fit(draws=20, tune=20, chains=1, cores=1, random_seed=42,
                      progressbar=False, compute_convergence_checks=False)
    names = ["Intercept", "x", "kappa_Intercept", "kappa_x"]
    assert idata.posterior.sizes["draw"] == 20
    for name in names:
        assert np.isfinite(idata.posterior[name]).all()
    summary = circular_summary(idata, var_names=names, circ_var_names=["Intercept"])
    assert np.isfinite(summary["mean"]).all()
