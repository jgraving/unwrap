import bambi as bmb
from pymc.distributions.transforms import LogExpM1 as Softplus  # softplus transform
from unwrap.link import inverse_softplus
from unwrap.distributions import CircularUniform


def CircularUniformPrior():
    """
    Defines a circular uniform prior using the Von Mises distribution with a near-zero concentration (kappa).

    This effectively approximates a uniform distribution on the circle by setting the concentration
    parameter (kappa) of the Von Mises distribution to a very small value.  The location parameter
    (mu) is set to 0, which can be considered an arbitrary reference point on the circle.

    Returns
    -------
    bambi.Prior
        A Bambi Prior object representing the circular uniform prior.  Specifically, a Von Mises
        distribution with mu=0 and kappa set to the smallest possible positive float value.
        `auto_scale` is set to False to prevent Bambi from automatically scaling the prior.
    """
    return bmb.Prior("CircularUniform", dist=CircularUniform)


def SoftplusNormalPrior(mu, sigma, *, mu_is_on_positive_scale=True):
    """
    Normal prior on the *pre-softplus* scale, with a Softplus transform applied so the
    resulting parameter lives on (0, ∞).

    Parameters
    ----------
    mu : float
        If mu_is_on_positive_scale=True (default), this is the desired location on the
        positive (post-softplus) scale (e.g., target kappa).
        If False, this is interpreted as the location on the pre-softplus scale.
    sigma : float
        Std dev on the pre-softplus scale.
    mu_is_on_positive_scale : bool
        Whether `mu` is specified on the positive scale (default) or pre-softplus scale.

    Returns
    -------
    bambi.Prior
        A Bambi prior: Normal(mu=..., sigma=...) with transform=Softplus().
    """
    mu0 = inverse_softplus(mu) if mu_is_on_positive_scale else mu
    return bmb.Prior("Normal", mu=mu0, sigma=sigma, transform=Softplus())
