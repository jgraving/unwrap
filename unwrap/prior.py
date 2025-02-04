import bambi as bmb
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
