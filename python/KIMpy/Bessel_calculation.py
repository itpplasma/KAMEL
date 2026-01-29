from math import factorial

import numpy as np
from scipy import integrate as spint
from scipy.special import gamma as Gamma
from scipy.special import iv as Bessel

# def calc_needed_bessel_of_mphi(mphi, eval_b):

# BesselProd0 = 0 + 0j
# BesselProd1 = 0 + 0j
# BesselProd0_store = 1.0 + 0j
# BesselProd1_store = 1.0 + 0j
# k = 0
# while (np.abs(BesselProd0 - BesselProd0_store) > 1e-10 or np.abs(BesselProd1 - BesselProd1_store) > 1e-10):
# BesselProd0_store = 1.0 * BesselProd0
# BesselProd1_store = 1.0 * BesselProd1
# BesselProd0 += (0.5 * eval_b)**(2*k) / (factorial(k) * Gamma(k+mphi+1) + 0j) * np.exp(-eval_b)
# BesselProd1 += (0.5 * eval_b)**(2*k) / (factorial(k) * Gamma(k+mphi) + 0j) * np.exp(-eval_b)
# k +=1
# if k > 50:
# break
# BesselProd0 *= (0.5 * eval_b)**mphi
# BesselProd1 *= (0.5 * eval_b)**(mphi-1)

# return BesselProd0, BesselProd1


def calc_needed_bessel_of_mphi(mphi, eval_b, mode="Integral"):
    if mode == "Integral":
        # If integration method quad is used in cxroots eval_b is a scalar
        if np.isscalar(eval_b):
            BesselProd0, _ = spint.quad(
                lambda theta: bessel_integrand(theta, eval_b, mphi), 0, np.pi, complex_func=True
            )
            BesselProd1, _ = spint.quad(
                lambda theta: bessel_integrand(theta, eval_b, mphi - 1), 0, np.pi, complex_func=True
            )
        # If integration method romb is used in cxroots eval_b may be vector -> this needs to be handled
        else:
            intfun = lambda eval_b: spint.quad(
                lambda theta: bessel_integrand(theta, eval_b, mphi), 0, np.pi, complex_func=True
            )[0]
            vec_int = np.vectorize(intfun)
            BesselProd0 = vec_int(eval_b)
            intfun = lambda eval_b: spint.quad(
                lambda theta: bessel_integrand(theta, eval_b, mphi - 1), 0, np.pi, complex_func=True
            )[0]
            vec_int = np.vectorize(intfun)
            BesselProd1 = vec_int(eval_b)
    elif mode == "Asymptotic expansion":
        if np.abs(eval_b) > 30:
            if np.abs(np.angle(eval_b)) < np.pi / 2:
                BesselProd0, BesselProd1 = BesselParabolicExpansion(eval_b)
            else:
                print("Argument b: %f" % (np.abs(np.angle(eval_b))))
        else:
            BesselProd0 = Bessel(mphi, eval_b) * np.exp(-eval_b)
            BesselProd1 = Bessel(mphi - 1, eval_b) * np.exp(-eval_b)
    return BesselProd0, BesselProd1


def bessel_integrand(theta: float, b: complex, m_phi: int) -> complex:
    exp_term = np.exp(b * (np.cos(theta) - 1))  # Complex exponential
    cos_term = np.cos(m_phi * theta)  # Real cosine (can also be complex if needed)

    return (1 / np.pi) * exp_term * cos_term


def BesselParabolicExpansion(z):
    f = (
        lambda z, n: 1
        / np.sqrt(2 * np.pi * z)
        * (
            1
            + (4 * n**2 - 1) / (8 * z)
            + ((4 * n**2 - 1) * (4 * n**2 - 9)) / (factorial(2) * (8 * z) ** 2)
            - ((4 * n**2 - 1) * (4 * n**2 - 9) * (4 * n**2 - 25) / (factorial(3) * (8 * z) ** 3))
        )
    )
    B0 = f(z, 0)
    B1 = f(z, -1)
    return B0, B1
