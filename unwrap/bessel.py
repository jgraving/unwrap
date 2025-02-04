import numpy as np
from scipy.special import i0, i1

def bessel_ratio(kappa: float) -> float:
    """
    Computes the ratio of modified Bessel functions of the first kind, I1(kappa) / I0(kappa).

    This ratio is equal to the mean resultant vector length (R) for a von Mises distribution
    with concentration parameter kappa.

    Args:
        kappa: Concentration parameter of the von Mises distribution (kappa >= 0).

    Returns:
        The mean resultant vector length R. Returns 1.0 if kappa is infinite.
    """
    return i1(kappa) / i0(kappa) if kappa < np.inf else 1.0


def inverse_bessel_ratio(R: float, tolerance: float = 1e-10, max_iterations: int = 1000) -> float:
    """
    Numerically solves for the concentration parameter kappa of a von Mises distribution,
    given the mean resultant vector length R.  This is the inverse of the bessel_ratio function.

    Uses a bisection method to find the root of the equation bessel_ratio(kappa) - R = 0.

    Args:
        R: The mean resultant vector length. Must be in the range [0, 1).
        tolerance: The desired accuracy for the solution.
        max_iterations: The maximum number of iterations for the bisection method.

    Returns:
        The estimated value of kappa.

    Raises:
        ValueError: If R is not in the valid range [0, 1).
        RuntimeError: If the root cannot be bracketed or the bisection method does not converge
                     within the maximum number of iterations.
    """
    if not (0 <= R < 1):
        raise ValueError("R must be in [0, 1).")

    if np.isclose(R, 0.0):
        return 0.0

    def objective(kappa: float) -> float:
        return bessel_ratio(kappa) - R

    a, b = 0.0, 10.0  # Initial bracket
    while objective(b) < 0:
        b *= 2
        if b > 1e12:
            raise RuntimeError("Could not bracket the root; R extremely close to 1.")

    for _ in range(max_iterations):
        mid = 0.5 * (a + b)
        fmid = objective(mid)
        if abs(fmid) < tolerance or b - a < tolerance:  # Convergence check
            return mid
        if fmid > 0:
            b = mid
        else:
            a = mid

    return 0.5 * (a + b)  # Return midpoint if max iterations reached (may not have converged)



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