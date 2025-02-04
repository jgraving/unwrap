import pymc as pm
import numpy as np

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
