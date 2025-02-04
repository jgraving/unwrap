import numpy as np
import pandas as pd
import pymc as pm
from typing import Tuple, Union
from unwrap.link.pymc import LINKS


def generate_circular_regression_data(
    num_samples: int,
    x_range: Tuple[Union[int, float], Union[int, float]],
    mean_intercept: float,
    mean_slope: float,
    conc_intercept: float,
    conc_slope: float,
    mean_link: str,
    conc_link: str,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Generates synthetic data for circular regression using PyMC with fixed parameter values.

    Parameters:
    - num_samples (int): Number of samples to generate.
    - x_range (tuple): Range (min, max) for the input variable x.
    - mean_intercept (float): Fixed intercept for the mean.
    - mean_slope (float): Fixed slope for the mean.
    - conc_intercept (float): Fixed intercept for the concentration parameter.
    - conc_slope (float): Fixed slope for the concentration parameter.
    - mean_link (str): Name of the link function for the mean (from LINKS).
    - conc_link (str): Name of the link function for the concentration (from LINKS).
    - seed (int): Random seed for reproducibility.

    Returns:
    - pd.DataFrame: DataFrame containing generated x, y, mean, and concentration values.
    """

    with pm.Model() as model:
        x = pm.Uniform("x", lower=x_range[0], upper=x_range[1], shape=num_samples)

        mean = pm.Deterministic(
            "mean", LINKS[mean_link].linkinv(mean_intercept + mean_slope * x)
        )
        concentration = pm.Deterministic(
            "concentration", LINKS[conc_link].linkinv(conc_intercept + conc_slope * x)
        )

        y = pm.VonMises("y", mu=mean, kappa=concentration, shape=num_samples)

        trace = pm.sample_prior_predictive(samples=1, random_seed=seed)

    x_samples = trace.prior["x"].values.flatten()
    y_samples = trace.prior["y"].values.flatten()
    mean_values = trace.prior["mean"].values.flatten()
    concentration_values = trace.prior["concentration"].values.flatten()

    # Combine into a DataFrame
    df = pd.DataFrame(
        {
            "x": x_samples,
            "y": y_samples,
            "mean": mean_values,
            "concentration": concentration_values,
        }
    )

    return df


if __name__ == "__main__":
    data_df = generate_circular_regression_data(
        num_samples=1000,
        x_range=(-1, 1),
        mean_intercept=np.pi,
        mean_slope=np.pi,
        conc_intercept=1,
        conc_slope=0.5,
        mean_link="shifted_modulo",
        conc_link="softplus",
    )
    print(data_df.head())
