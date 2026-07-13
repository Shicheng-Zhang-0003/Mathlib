#!/usr/bin/env python3
"""
MathLib V1.0.5 Minimax Polynomial Generator
Uses Chebyshev economization (the industry standard for near-minimax approximations)
to derive the optimal polynomial coefficients for fast math functions.

Usage: pip install numpy && python3 compute_minimax.py
"""
import numpy as np
from numpy.polynomial import Chebyshev

def compute_minimax(f, a, b, degree):
    """Fits a Chebyshev polynomial and converts it to the standard power basis."""
    x = np.linspace(a, b, 10000)
    y = f(x)
    cheb = Chebyshev.fit(x, y, deg=degree, domain=[a, b])
    poly = cheb.convert()
    return poly.coef

if __name__ == "__main__":
    print("Computing true Minimax coefficients for 2^x - 1 on [0, 1]...")
    coeffs = compute_minimax(lambda x: 2**x - 1, 0.0, 1.0, 5)
    print("C Array: {", ", ".join(f"{c:.16e}" for c in coeffs[1:]), "}")

    print("\nComputing true Minimax coefficients for log2(1+x) on [0, 1]...")
    coeffs_log = compute_minimax(lambda x: np.log2(1 + x), 0.0, 1.0, 3)
    print("C Array: {", ", ".join(f"{c:.16e}" for c in coeffs_log[1:]), "}")
