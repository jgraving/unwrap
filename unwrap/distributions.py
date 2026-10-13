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
