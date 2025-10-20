import pymc as pm
import numpy as np
import pytensor.tensor as pt
from pymc.distributions.dist_math import check_parameters

# Define a very small concentration value (effectively uniform on the circle)
EPS = np.finfo(float).eps


def CircularUniform(name, *args, dims=None, **kwargs):
    """
    Defines a circular uniform prior using the Von Mises distribution with a near-zero concentration (kappa).

    This effectively approximates a uniform distribution on the circle by setting the concentration
    parameter (kappa) of the Von Mises distribution to a very small value. The location parameter
    (mu) is set to 0, which can be considered an arbitrary reference point on the circle.
    """
    return pm.VonMises(name, mu=0, kappa=EPS, *args, dims=dims, **kwargs)


class VonMisesMixture:
    r"""
    Von Mises mixture log-likelihood.

    .. math::

        f(x \mid w, \mu, \kappa) = \sum_{i = 1}^{n} w_i \, \text{VM}(x \mid \mu_i, \kappa_i)

    where the Von Mises density is given by

    .. math::

        \text{VM}(x \mid \mu, \kappa) = \frac{\exp\bigl(\kappa \cos(x - \mu)\bigr)}{2\pi I_0(\kappa)}

    Parameters
    ----------
    w : tensor_like of float
        The mixture weights, with each element in [0, 1] and sum(w) = 1.
    mu : tensor_like of float
        The component means (angles, in radians).
    kappa : tensor_like of float
        The concentration parameters for the components.
    """

    def __new__(cls, name, w, mu, kappa, **kwargs):
        # This creates a new Mixture distribution in a model context.
        return pm.Mixture(name, w, pm.VonMises.dist(mu=mu, kappa=kappa), **kwargs)

    @classmethod
    def dist(cls, w, mu, kappa, **kwargs):
        # This creates a free (standalone) random variable.
        return pm.Mixture.dist(w, pm.VonMises.dist(mu=mu, kappa=kappa), **kwargs)


class BimodalVonMises:
    r"""
    A distribution representing a mixture of two Von Mises distributions, used to model circular data with bimodal characteristics.

    This distribution is parameterized by two means (mu) and two concentration parameters (kappa), allowing for flexible modeling of directional data with two modes.

    The bimodal Von Mises mixture log-likelihood is given by:

    .. math::

        f(x \mid p, \mu, \kappa) = p \text{VM}(x \mid \mu_1, \kappa_1) + (1 - p) \text{VM}(x \mid \mu_2, \kappa_2)

    where the Von Mises density is defined as:

    .. math::

        \text{VM}(x \mid \mu, \kappa) = \frac{\exp(\kappa \cos(x - \mu))}{2\pi I_0(\kappa)}

    Parameters
    ----------
    p : tensor_like of float
        The probability of the first component (0 < p < 1).
    mu : tensor_like of float, 2D tensor
        The component means (angles, in radians). Must have last dimension of size 2.
    kappa : tensor_like of float, 2D tensor
        The concentration parameters for the components. Must have last dimension of size 2.
    kwargs : dict, optional
        Additional keyword arguments passed to the Mixture distribution.
    """

    def __new__(cls, name, mu, kappa, p, **kwargs):
        if mu.shape[-1] != 2 or kappa.shape[-1] != 2:
            raise ValueError("mu and kappa must have last dimension of size 2.")

        w = pt.stack([p, 1 - p], axis=-1)

        return VonMisesMixture(name, w, mu, kappa, **kwargs)

    @classmethod
    def dist(cls, mu, kappa, p, **kwargs):
        if mu.shape[-1] != 2 or kappa.shape[-1] != 2:
            raise ValueError("mu and kappa must have last dimension of size 2.")

        w = pt.stack([p, 1 - p], axis=-1)

        return VonMisesMixture.dist(w, mu, kappa, **kwargs)


def vmss_logp(theta, mu, kappa, R, n):
    r"""
    Elementwise log-likelihood for the **sufficient-statistics von Mises** model
    with **observed** per-unit mean angles ``theta`` and known summaries ``R`` and ``n``.

    Data (observed / known)
    -----------------------
    theta : tensor-like, shape (N,)
        Per-unit circular **mean direction** (radians) in (-π, π] — this is the
        observed value passed to the likelihood (i.e., `observed=theta`).
    R : tensor-like, shape (N,)
        Per-unit **mean resultant length** in [0, 1].
    n : tensor-like, shape (N,)
        Per-unit sample sizes (≥ 1).

    Parameters (to be estimated)
    ----------------------------
    mu : tensor-like, shape (N,)
        Per-unit mean direction (radians), typically a single linear predictor, e.g. `mu = X @ beta`.
        No modulo needed; periodicity is handled by `cos(mu - theta)`.
    kappa : tensor-like, shape (N,)
        Per-unit concentration (> 0), usually via a positive link such as `softplus(X @ eta)`.

    Model & exact joint log-likelihood
    ----------------------------------
    If step angles θ_ij are iid VM(μ_i, κ_i), then the joint log-likelihood condensed by
    sufficient statistics equals:

        log L_i = κ_i * (n_i * R_i) * cos(μ_i - θ_i)  -  n_i * log(2π)  -  n_i * log I0(κ_i),

    where R_i is the **mean** resultant length and θ_i is the circular mean.

    Returns
    -------
    logp : tensor, shape (N,)
        Elementwise log-likelihood; PyMC will sum across the observed axis.
    """
    # match PyMC VonMises style: res → support mask → check_parameters
    res = kappa * n * R * pt.cos(mu - theta) - n * (
        pt.log(2.0 * np.pi) + pt.log(pt.i0(kappa))
    )

    # support: theta ∈ [-π, π]
    res = pt.switch(
        pt.bitwise_and(pt.ge(theta, -np.pi), pt.le(theta, np.pi)),
        res,
        -np.inf,
    )

    # parameter & data checks
    return check_parameters(
        res,
        (kappa > 0) & (R >= 0) & (R <= 1) & (n >= 1),
        msg="kappa > 0, 0 <= R <= 1, n >= 1",
    )


def VonMisesSS(name, mu, kappa, R, n, **kwargs):
    r"""
    Sufficient Statistics von Mises with observed per-unit mean angles, scaled by sample size `n`
    and mean resultant vector length `R`.

    Call with `observed=theta`, and pass `R` (mean resultant length) and `n` as known arrays.

    Parameters
    ----------
    name : str
        Random variable name.
    mu : tensor-like, shape (N,)
        Per-unit mean direction.
    kappa : tensor-like, shape (N,)
        Per-unit concentration (> 0).
    R : tensor-like, shape (N,)
        Per-unit **mean** resultant length in [0, 1].
    n : tensor-like, shape (N,)
        Per-unit step counts (≥ 1).
    **kwargs :
        Forwarded to `pm.CustomDist` (`observed=theta`).

    Example
    -------
    >>> beta  = pm.Normal("beta", 0.0, 1.5, shape=P)
    >>> eta   = pm.Normal("eta",  0.0, 5.0, shape=P)
    >>> mu    = X @ beta
    >>> kappa = pm.math.softplus(X @ eta)
    >>> vm = VonMisesSS(
    ...     "vm_ss",
    ...     mu=mu, kappa=kappa,
    ...     R=R, n=n,
    ...     observed=theta,
    ...     shape=n.shape,
    ... )
    """
    # Note: vmss_logp signature must match (value, *params)
    return pm.CustomDist(
        name,
        mu,
        kappa,
        R,
        n,
        logp=vmss_logp,
        **kwargs,
    )
