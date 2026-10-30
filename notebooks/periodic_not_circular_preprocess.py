"""Data simulation for notebooks/periodic_not_circular.ipynb.

The study fits a Beta model of 24-hour activity whose mean (``mu``) and
scale (``phi``) follow fixed circadian sin/cos shapes. The data are
simulated deterministically (seed 42, 1000 points) — bit-identical to the
old inline simulation cell.
"""

import numpy as np


def generate_periodic_beta_activity(num_points=1000, seed=42):
    """Simulate one day of activity data for the periodic Beta model.

    Parameters
    ----------
    num_points : int
        Grid points on the day (default 1000, as in the notebook).
    seed : int
        Seed of the global NumPy RNG (default 42, as in the notebook).

    Returns
    -------
    (time, activity)
        ``time`` in hours on ``[0, 24]`` (``num_points`` points) and the
        simulated ``activity`` curve, with values in ``(0, 1)``.
    """
    np.random.seed(seed)
    time = np.linspace(0, 24, num_points)

    mu_alpha_true = 0.5
    mu_sin_true = 2.0
    mu_cos_true = 1.5
    phi_alpha_true = 2.0
    phi_sin_true = 1.0
    phi_cos_true = 0.5

    linear_mu_true = mu_alpha_true + mu_sin_true*np.sin(2*np.pi*time/24) + mu_cos_true*np.cos(2*np.pi*time/24)
    mu_true = 1 / (1 + np.exp(-linear_mu_true))

    linear_phi_true = phi_alpha_true + phi_sin_true*np.sin(2*np.pi*time/24) + phi_cos_true*np.cos(2*np.pi*time/24)
    phi_true = np.log1p(np.exp(linear_phi_true))

    alpha_beta_true = mu_true * phi_true
    beta_beta_true = (1 - mu_true) * phi_true

    activity = np.random.beta(alpha_beta_true, beta_beta_true)
    return time, activity
