import numpy as np
from numbers import Integral
from scipy.special import i0e, i1e

def bessel_ratio(kappa: float) -> float:
    """
    Computes the ratio of modified Bessel functions of the first kind, I1(kappa) / I0(kappa).

    This ratio is equal to the mean resultant vector length (R) for a von Mises distribution
    with concentration parameter kappa.

    Args:
        kappa: Scalar concentration parameter of the von Mises distribution (kappa >= 0).

    Returns:
        The mean resultant vector length R. Returns 1.0 if kappa is infinite.
    """
    if np.isnan(kappa) or kappa < 0:
        raise ValueError("kappa must be nonnegative and not NaN.")
    if np.isposinf(kappa):
        return 1.0
    # Exponential scaling cancels in the ratio and prevents inf / inf.
    return i1e(kappa) / i0e(kappa)


def inverse_bessel_ratio(R: float, tolerance: float = 1e-10, max_iterations: int = 1000) -> float:
    """
    Numerically solves for the concentration parameter kappa of a von Mises distribution,
    given the mean resultant vector length R.  This is the inverse of the bessel_ratio function.

    Uses a bisection method to find the root of the equation bessel_ratio(kappa) - R = 0.

    Args:
        R: The mean resultant vector length. Must be in the range [0, 1).
        tolerance: Positive relative tolerance for the concentration bracket.
            Accuracy near R=1 is ultimately limited by floating-point precision.
        max_iterations: The maximum number of iterations for the bisection method.

    Returns:
        The estimated value of kappa.

    Raises:
        ValueError: If R is outside [0, 1) or solver controls are invalid.
        RuntimeError: If the root cannot be bracketed or the bisection method does not converge
                     within the maximum number of iterations.
    """
    if not (0 <= R < 1):
        raise ValueError("R must be in [0, 1).")

    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("tolerance must be finite and positive.")
    if not isinstance(max_iterations, Integral) or max_iterations < 1:
        raise ValueError("max_iterations must be a positive integer.")
    if R == 0.0:
        return 0.0

    def objective(kappa: float) -> float:
        return bessel_ratio(kappa) - R

    # A scale-aware initial bracket also resolves very small positive R.
    a, b = 0.0, max(4.0 * R, 1.0 / (1.0 - R))
    while objective(b) < 0:
        b *= 2
        if not np.isfinite(b):
            raise RuntimeError("Could not bracket the root at floating-point precision.")

    for _ in range(max_iterations):
        mid = a + 0.5 * (b - a)
        if mid == a or mid == b or b - a <= tolerance * mid:
            return mid
        if objective(mid) > 0:
            b = mid
        else:
            a = mid
    raise RuntimeError("Concentration inversion did not converge within max_iterations.")



def test_kappa_inversion():
    """
    Tests the kappa inversion by comparing true and estimated kappa values.
    Prints a table showing the true kappa, the corresponding R value, the estimated kappa
    (using the inverse function), and the error.
    """
    kappa_values = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 500.0, np.inf]
    print("     kappa_true        R_true        kappa_est       R_est         error         ")
    print("-----------------------------------------------------------------------------")
    for k_true in kappa_values:
        R_true = bessel_ratio(k_true)

        k_num = inverse_bessel_ratio(R_true) if R_true < 1.0 else np.inf

        err_num = k_num - k_true
        R_num = bessel_ratio(k_num)

        print(f"{k_true:12.4f}  {R_true:12.6f}  {k_num:12.6f}  {R_num:12.6f}  {err_num:12.3e} ")


if __name__ == "__main__":
    test_kappa_inversion()