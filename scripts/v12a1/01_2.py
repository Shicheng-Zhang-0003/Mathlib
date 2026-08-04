#!/usr/bin/env python3
"""
01c_minimax_chebyshev.py

Run from the folder that CONTAINS the v12A1 working folder.

Generates near-minimax polynomial coefficients using Chebyshev
interpolation (interpolate at Chebyshev nodes, convert to power
basis via cheb2poly, compose with domain map).

This replaces the broken np.polyfit approach (least-squares, not
minimax) with proper Chebyshev interpolation (near-minimax).

Requires: pip install numpy mpmath

Usage:
    python3 01c_minimax_chebyshev.py
    python3 01c_minimax_chebyshev.py --force
"""
from __future__ import annotations
import struct
import sys
from datetime import date
from pathlib import Path

try:
    import numpy as np
    from numpy.polynomial.chebyshev import chebfit, cheb2poly
    from numpy.polynomial import Polynomial
except ImportError:
    print("ERROR: pip install numpy"); sys.exit(1)

try:
    from mpmath import mp, mpf, sin, cos, exp, log, pi, ln
    mp.dps = 80
except ImportError:
    print("ERROR: pip install mpmath"); sys.exit(1)

MARKER = "MATHLIB_V12A1_MINIMAX_COEFFS"
TODAY = date.today().isoformat()
PI4 = float(mpf(str(float(pi))) / 4)
LN2H = float(mpf(str(float(ln(2)))) / 2)
ZMAX = float((mpf(2).sqrt() - 1) / (mpf(2).sqrt() + 1))

# ---------------------------------------------------------------------------
# Chebyshev interpolation -> power basis
# ---------------------------------------------------------------------------
def cheb_interp_to_power(f, a, b, deg):
    """
    Chebyshev interpolation on [a, b], converted to power basis
    for direct evaluation of x in [a, b].

    Uses n = deg+1 Chebyshev nodes (exact interpolation).
    Returns ascending-order power-basis coefficients.
    """
    n = deg + 1
    k = np.arange(n)
    # Chebyshev nodes of the first kind in [-1, 1]
    u_nodes = np.cos((2*k + 1) * np.pi / (2*n))
    # Map to [a, b]
    t_nodes = a + (b - a) * (u_nodes + 1) / 2
    # Evaluate function at nodes
    y_nodes = np.array([float(f(mpf(str(float(t))))) for t in t_nodes])
    # Fit Chebyshev series on [-1, 1] (exact interpolation at Chebyshev nodes)
    cheb_coeffs = chebfit(u_nodes, y_nodes, deg)
    # Convert to power basis on [-1, 1]
    power_u = cheb2poly(cheb_coeffs)
    # Compose with linear map u = alpha*t + beta to get power basis in t
    alpha = 2.0 / (b - a)
    beta = -(a + b) / (b - a)
    poly_u = Polynomial(power_u)
    linear = Polynomial([beta, alpha])  # u = beta + alpha*t
    poly_t = poly_u(linear)
    return poly_t.coef.tolist()

def cheb_interp_odd(f, a, b, deg):
    """
    Fit odd polynomial: f(x) = x * P(x^2).
    Fits P(t) on t in [a^2, b^2] via Chebyshev interpolation.
    Returns ascending-order power-basis coefficients of P.
    """
    a2 = a * a
    b2 = b * b
    n = deg + 1
    k = np.arange(n)
    u_nodes = np.cos((2*k + 1) * np.pi / (2*n))
    t_nodes = a2 + (b2 - a2) * (u_nodes + 1) / 2
    x_nodes = np.sqrt(t_nodes)
    y_nodes = np.array([float(f(mpf(str(float(x)))) / x) for x in x_nodes])
    cheb_coeffs = chebfit(u_nodes, y_nodes, deg)
    power_u = cheb2poly(cheb_coeffs)
    alpha = 2.0 / (b2 - a2)
    beta = -(a2 + b2) / (b2 - a2)
    poly_u = Polynomial(power_u)
    linear = Polynomial([beta, alpha])
    poly_t = poly_u(linear)
    return poly_t.coef.tolist()

# ---------------------------------------------------------------------------
# ULP validation
# ---------------------------------------------------------------------------
def d2i(d):
    b = struct.unpack('Q', struct.pack('d', d))[0]
    return (0x8000000000000000 - b) if (b >> 63) else b

def ulp(a, b):
    return abs(d2i(a) - d2i(b))

def horner(c, x):
    r = c[-1]
    for i in range(len(c)-2, -1, -1):
        r = r * x + c[i]
    return r

def max_ulp_odd(c, f, pts):
    mx = 0
    for x in pts:
        g = x * horner(c, x*x)
        e = float(f(mpf(str(x))))
        if e == 0.0 and g == 0.0: continue
        u = ulp(g, e)
        if u > mx: mx = u
    return mx

def max_ulp_even(c, f, pts):
    mx = 0
    for x in pts:
        g = horner(c, x*x)
        e = float(f(mpf(str(x))))
        u = ulp(g, e)
        if u > mx: mx = u
    return mx

def max_ulp_plain(c, f, pts):
    mx = 0
    for x in pts:
        g = horner(c, x)
        e = float(f(mpf(str(x))))
        if e == 0.0 and g == 0.0: continue
        u = ulp(g, e)
        if u > mx: mx = u
    return mx

# ---------------------------------------------------------------------------
# Search for degree that achieves <= 5 ULP
# ---------------------------------------------------------------------------
def search_deg(name, fit_fn, f, pts, eval_fn, deg_start, deg_end):
    best_c, best_u, best_d = None, 10**18, 0
    for deg in range(deg_start, deg_end + 1):
        c = fit_fn(f, 0 if name in ("sin","log") else -LN2H if name == "exp" else 0,
                   PI4 if name in ("sin","cos") else LN2H if name == "exp" else ZMAX,
                   deg)
        u = eval_fn(c, f, pts)
        print(f"    {name} deg={deg}: {u} ULP")
        if u < best_u:
            best_c, best_u, best_d = c, u, deg
        if u <= 5:
            return c, u, deg
    print(f"    WARNING: {name} best={best_u} ULP at deg={best_d}")
    return best_c, best_u, best_d

# ---------------------------------------------------------------------------
# C header
# ---------------------------------------------------------------------------
def c_array(name, c):
    lines = [f"static const double {name}[] = {{"]
    for i, v in enumerate(c):
        lines.append(f"    {v:.17e}{',' if i < len(c)-1 else ''}")
    lines.append("};")
    return "\n".join(lines)

def header(sc, cc, ec, lc, sd, cd, ed, ld):
    return "\n".join([
        "#ifndef ML_INTERNAL_MINIMAX_COEFFS_H",
        "#define ML_INTERNAL_MINIMAX_COEFFS_H", "",
        f"/* {MARKER} */",
        f"/* Auto-generated by 01c_minimax_chebyshev.py on {TODAY} */",
        "/* DO NOT EDIT MANUALLY. Regenerate with the pipeline script. */", "",
        f"/* sin(x) = x*P(x^2) on [-pi/4,pi/4], deg P = {sd} */",
        c_array("minimax_sin_coeffs", sc), "",
        f"/* cos(x) = P(x^2) on [-pi/4,pi/4], deg P = {cd} */",
        c_array("minimax_cos_coeffs", cc), "",
        f"/* exp(x) = P(x) on [-ln2/2,ln2/2], deg = {ed} */",
        c_array("minimax_exp_coeffs", ec), "",
        f"/* log((1+z)/(1-z)) = z*P(z^2) on [0,{ZMAX:.17e}], deg P = {ld} */",
        c_array("minimax_log_coeffs", lc), "",
        "#endif /* ML_INTERNAL_MINIMAX_COEFFS_H */", "",
    ])

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    force = "--force" in sys.argv[1:]
    root = Path.cwd()
    v12 = root / "v12A1"
    if not v12.is_dir():
        if (root / "src" / "internal").is_dir():
            v12 = root; root = root.parent
        else:
            fail("Run from folder containing v12A1/")

    hp = v12 / "src" / "internal" / "minimax_coeffs.h"
    rp = v12 / "docs" / "V12A1_MINIMAX_REPORT.md"

    if hp.exists() and not force and MARKER in hp.read_text():
        print(f"  [skip] {hp}: already generated. Use --force.")
        return 0

    print("=========================================================")
    print("  MATHLIB v12A1: MINIMAX ULP FIX (01c Chebyshev)")
    print("=========================================================")

    pts_sin = np.linspace(-PI4, PI4, 10000).tolist()
    pts_cos = np.linspace(-PI4, PI4, 10000).tolist()
    pts_exp = np.linspace(-LN2H, LN2H, 10000).tolist()
    pts_log = np.linspace(1e-15, ZMAX, 10000).tolist()

    print("\n[sin] searching for <= 5 ULP ...")
    sc, su, sd = search_deg("sin", cheb_interp_odd, sin, pts_sin, max_ulp_odd, 9, 24)

    print("\n[cos] searching for <= 5 ULP ...")
    cc, cu, cd = search_deg("cos", cheb_interp_odd, cos, pts_cos, max_ulp_even, 9, 24)

    print("\n[exp] searching for <= 5 ULP ...")
    ec, eu, ed = search_deg("exp", cheb_interp_to_power, exp, pts_exp, max_ulp_plain, 13, 24)

    print("\n[log] searching for <= 5 ULP ...")
    def log_atanh(z): return log((1+z)/(1-z))
    lc, lu, ld = search_deg("log", cheb_interp_odd, log_atanh, pts_log, max_ulp_odd, 10, 24)

    print(f"\n{'='*57}")
    print(f"  sin: deg={sd}, {su} ULP")
    print(f"  cos: deg={cd}, {cu} ULP")
    print(f"  exp: deg={ed}, {eu} ULP")
    print(f"  log: deg={ld}, {lu} ULP")

    all_pass = all(u <= 5 for u in [su, cu, eu, lu])
    if not all_pass:
        print("\n  WARNING: Not all functions reached <= 5 ULP.")
        print("  Coefficients written anyway. Manual review needed.")

    print(f"\n{'='*57}")
    hp.parent.mkdir(parents=True, exist_ok=True)
    hp.write_text(header(sc, cc, ec, lc, sd, cd, ed, ld), encoding="utf-8", newline="\n")
    print(f"  [write] {hp}")

    rp.parent.mkdir(parents=True, exist_ok=True)
    rp.write_text("\n".join([
        "# v12A1 Minimax Coefficient Validation Report", "",
        f"Generated: {TODAY}",
        "Method: Chebyshev interpolation (chebfit + cheb2poly + composition)",
        "Ground truth: mpmath at 80 decimal places", "",
        "## Results", "",
        "| Function | Degree | Max ULP | Status |",
        "|----------|--------|---------|--------|",
        f"| sin | {sd} | {su} | {'PASS' if su<=5 else 'WARN'} |",
        f"| cos | {cd} | {cu} | {'PASS' if cu<=5 else 'WARN'} |",
        f"| exp | {ed} | {eu} | {'PASS' if eu<=5 else 'WARN'} |",
        f"| log | {ld} | {lu} | {'PASS' if lu<=5 else 'WARN'} |",
        "", "## Acceptance Criteria", "",
        "All functions must achieve <= 5 ULP on their reduced domain.", "",
    ]), encoding="utf-8", newline="\n")
    print(f"  [write] {rp}")

    try:
        import shutil
        src = Path(__file__).resolve()
        dst = v12 / "scripts" / "v12a1" / src.name
        if src != dst:
            dst.parent.mkdir(parents=True, exist_ok=True)
            if not dst.exists() or force:
                shutil.copy2(src, dst)
                print(f"  [archive] {dst}")
    except NameError:
        pass

    print(f"\n{'='*57}")
    if all_pass:
        print("  ALL FUNCTIONS PASS (<= 5 ULP)")
    else:
        print("  WARNING: Some functions exceed 5 ULP. Review needed.")
    print(f"\n  Next: 05_trig_minimax.py (swap coefficients into trig.c)")
    print("=========================================================")
    return 0 if all_pass else 1

if __name__ == "__main__":
    raise SystemExit(main())
