#!/usr/bin/env python3
"""
01_minimax_pipeline.py  (v3 — polyfit on Chebyshev nodes)

Run from the folder that CONTAINS the v12A1 working folder.

FIX HISTORY:
  v1: Used Chebyshev.fit().convert().coef — coefficients were for the
      mapped variable u in [-1,1], not for x. Garbage.
  v2: Tried composing with the linear domain map via Polynomial.__call__.
      Numpy's domain/window metadata corrupted the composition. Garbage.
  v3: Uses np.polyfit on Chebyshev-distributed nodes. Returns plain
      power-basis coefficients for direct Horner evaluation. No domain,
      no window, no composition. Just math.

Requires:
    pip install numpy mpmath

Outputs:
    v12A1/src/internal/minimax_coeffs.h
    v12A1/docs/V12A1_MINIMAX_REPORT.md

Usage:
    python3 01_minimax_pipeline.py --force
"""

from __future__ import annotations

import sys
from datetime import date
from pathlib import Path

try:
    import numpy as np
except ImportError:
    print("ERROR: numpy not installed. Run: pip install numpy")
    sys.exit(1)

try:
    from mpmath import mp, mpf, sin, cos, exp, log, pi, ln
    mp.dps = 80
except ImportError:
    print("ERROR: mpmath not installed. Run: pip install mpmath")
    sys.exit(1)

MARKER = "MATHLIB_V12A1_MINIMAX_COEFFS"
TODAY = date.today().isoformat()

PI_OVER_4 = float(mpf(str(float(pi))) / 4)
LN2_OVER_2 = float(mpf(str(float(ln(2)))) / 2)
Z_MAX = float((mpf(2).sqrt() - 1) / (mpf(2).sqrt() + 1))


# ---------------------------------------------------------------------------
# Fitting via np.polyfit on Chebyshev nodes
# ---------------------------------------------------------------------------
def chebyshev_nodes(a, b, n):
    """n Chebyshev nodes on [a, b], mapped from cos formula."""
    k = np.arange(n)
    return 0.5 * (a + b) + 0.5 * (b - a) * np.cos((2*k + 1) * np.pi / (2*n))


def chebyshev_fit(f, a, b, degree, n_nodes=500):
    """
    Near-minimax polynomial fit using Chebyshev nodes + least squares.
    Returns ascending-order coefficients for direct Horner evaluation:
        P(x) = c[0] + c[1]*x + c[2]*x^2 + ...
    """
    nodes = chebyshev_nodes(a, b, n_nodes)
    y = np.array([float(f(mpf(str(float(xi))))) for xi in nodes])
    # np.polyfit returns DESCENDING order: [c_n, ..., c_1, c_0]
    coeffs_desc = np.polyfit(nodes, y, degree)
    # Reverse to ascending: [c_0, c_1, ..., c_n]
    return coeffs_desc[::-1].tolist()


def fit_odd(f, a, b, degree, n_nodes=500):
    """
    Fit f(x) = x * P(x^2) on [0, b].
    Returns ascending-order coefficients of P for direct evaluation.
    """
    t_max = b * b
    nodes_t = chebyshev_nodes(1e-30, t_max, n_nodes)
    nodes_x = np.sqrt(nodes_t)
    y_over_x = np.array([
        float(f(mpf(str(float(xi)))) / xi) for xi in nodes_x
    ])
    coeffs_desc = np.polyfit(nodes_t, y_over_x, degree)
    return coeffs_desc[::-1].tolist()


def fit_even(f, a, b, degree, n_nodes=500):
    """
    Fit f(x) = P(x^2) on [0, b] (f is even).
    Returns ascending-order coefficients of P for direct evaluation.
    """
    t_max = b * b
    nodes_t = chebyshev_nodes(1e-30, t_max, n_nodes)
    nodes_x = np.sqrt(nodes_t)
    y = np.array([float(f(mpf(str(float(xi))))) for xi in nodes_x])
    coeffs_desc = np.polyfit(nodes_t, y, degree)
    return coeffs_desc[::-1].tolist()


# ---------------------------------------------------------------------------
# ULP validation
# ---------------------------------------------------------------------------
def double_to_int(d):
    import struct
    bits = struct.unpack('Q', struct.pack('d', d))[0]
    if bits >> 63:
        bits = 0x8000000000000000 - bits
    return bits


def ulp_distance(a, b):
    return abs(double_to_int(a) - double_to_int(b))


def eval_poly(coeffs, x):
    """Horner: c[0] + c[1]*x + c[2]*x^2 + ..."""
    result = coeffs[-1]
    for i in range(len(coeffs) - 2, -1, -1):
        result = result * x + coeffs[i]
    return result


def eval_odd_poly(coeffs, x):
    """x * P(x^2)"""
    return x * eval_poly(coeffs, x * x)


def eval_even_poly(coeffs, x):
    """P(x^2)"""
    return eval_poly(coeffs, x * x)


def validate(name, coeffs, f_mpmath, test_points, eval_fn):
    max_ulp = 0
    worst_x = 0.0
    for x in test_points:
        got = eval_fn(coeffs, x)
        expected = float(f_mpmath(mpf(str(x))))
        if expected == 0.0 and got == 0.0:
            continue
        ulps = ulp_distance(got, expected)
        if ulps > max_ulp:
            max_ulp = ulps
            worst_x = x
    status = "PASS" if max_ulp <= 5 else "WARN"
    line = f"| {name} | {max_ulp} | {worst_x:.17e} | {status} |"
    return max_ulp, line


# ---------------------------------------------------------------------------
# C header + report generation
# ---------------------------------------------------------------------------
def format_c_array(name, coeffs):
    lines = [f"static const double {name}[] = {{"]
    for i, c in enumerate(coeffs):
        comma = "," if i < len(coeffs) - 1 else ""
        lines.append(f"    {c:.17e}{comma}")
    lines.append("};")
    return "\n".join(lines)


def generate_header(sin_c, cos_c, exp_c, log_c):
    parts = [
        "#ifndef ML_INTERNAL_MINIMAX_COEFFS_H",
        "#define ML_INTERNAL_MINIMAX_COEFFS_H",
        "",
        f"/* {MARKER} */",
        f"/* Auto-generated by 01_minimax_pipeline.py on {TODAY} */",
        "/* DO NOT EDIT MANUALLY. Regenerate with the pipeline script. */",
        "",
        "/*",
        " * sin(x) on [-pi/4, pi/4]",
        " * Odd polynomial: sin(x) = x * P(x^2)",
        f" * Degree of P: {len(sin_c) - 1}",
        " */",
        format_c_array("minimax_sin_coeffs", sin_c),
        "",
        "/*",
        " * cos(x) on [-pi/4, pi/4]",
        " * Even polynomial: cos(x) = P(x^2)",
        f" * Degree of P: {len(cos_c) - 1}",
        " */",
        format_c_array("minimax_cos_coeffs", cos_c),
        "",
        "/*",
        " * exp(x) on [-ln2/2, ln2/2]",
        " * Full polynomial: exp(x) = P(x)",
        f" * Degree: {len(exp_c) - 1}",
        " */",
        format_c_array("minimax_exp_coeffs", exp_c),
        "",
        "/*",
        f" * log((1+z)/(1-z)) on z in [0, {Z_MAX:.17e}]",
        " * Odd polynomial: f(z) = z * P(z^2)",
        f" * Degree of P: {len(log_c) - 1}",
        " */",
        format_c_array("minimax_log_coeffs", log_c),
        "",
        "#endif /* ML_INTERNAL_MINIMAX_COEFFS_H */",
        "",
    ]
    return "\n".join(parts)


def generate_report(results):
    lines = [
        "# v12A1 Minimax Coefficient Validation Report",
        "",
        f"Generated: {TODAY}",
        "Method: np.polyfit on Chebyshev-distributed nodes",
        "Ground truth: mpmath at 80 decimal places",
        "",
        "## Results",
        "",
        "| Function | Max ULP | Worst Input | Status |",
        "|----------|---------|-------------|--------|",
    ]
    lines.extend(results)
    lines += ["", "## Acceptance Criteria", "",
              "All functions must achieve <= 5 ULP on their reduced domain.", ""]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def fail(msg):
    print("ERROR: " + msg)
    sys.exit(1)


def locate_v12a1():
    root = Path.cwd()
    candidate = root / "v12A1"
    if candidate.is_dir():
        return root, candidate
    if (root / "src" / "internal").is_dir() and (root / "include" / "mathlib").is_dir():
        return root.parent, root
    fail("Run from the folder that CONTAINS v12A1/, or from inside v12A1/ itself.")


def write_text(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(content)
    print(f"  [write] {path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    force = "--force" in sys.argv[1:]
    root, v12 = locate_v12a1()

    print("=========================================================")
    print("  MATHLIB v12A1: MINIMAX PIPELINE (v3 — polyfit)")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Date:  {TODAY}")
    print("---------------------------------------------------------")

    header_path = v12 / "src" / "internal" / "minimax_coeffs.h"
    report_path = v12 / "docs" / "V12A1_MINIMAX_REPORT.md"

    if header_path.exists() and not force:
        if MARKER in header_path.read_text(encoding="utf-8"):
            print(f"  [skip] {header_path}: already generated. Use --force.")
            return 0

    # --- sin: odd poly, sin(x) = x * P(x^2) on [0, pi/4] ---
    print("\n[1/4] sin(x) on [-pi/4, pi/4] ...")
    sin_c = fit_odd(sin, 0, PI_OVER_4, degree=9)
    print(f"       Degree of P(x^2): {len(sin_c) - 1}")
    print(f"       c[0] = {sin_c[0]:.17e}  (expect ~1.0)")
    print(f"       c[1] = {sin_c[1]:.17e}  (expect ~-0.1667)")

    # --- cos: even poly, cos(x) = P(x^2) on [0, pi/4] ---
    print("\n[2/4] cos(x) on [-pi/4, pi/4] ...")
    cos_c = fit_even(cos, 0, PI_OVER_4, degree=9)
    print(f"       Degree of P(x^2): {len(cos_c) - 1}")
    print(f"       c[0] = {cos_c[0]:.17e}  (expect ~1.0)")
    print(f"       c[1] = {cos_c[1]:.17e}  (expect ~-0.5)")

    # --- exp: full poly on [-ln2/2, ln2/2] ---
    print("\n[3/4] exp(x) on [-ln2/2, ln2/2] ...")
    exp_c = chebyshev_fit(exp, -LN2_OVER_2, LN2_OVER_2, degree=13)
    print(f"       Degree: {len(exp_c) - 1}")
    print(f"       c[0] = {exp_c[0]:.17e}  (expect ~1.0)")
    print(f"       c[1] = {exp_c[1]:.17e}  (expect ~1.0)")

    # --- log: odd poly, log((1+z)/(1-z)) = z * P(z^2) ---
    print("\n[4/4] log((1+z)/(1-z)) on [0, z_max] ...")
    def log_atanh(z):
        return log((1 + z) / (1 - z))
    log_c = fit_odd(log_atanh, 0, Z_MAX, degree=10)
    print(f"       Degree of P(z^2): {len(log_c) - 1}")
    print(f"       c[0] = {log_c[0]:.17e}  (expect ~2.0)")
    print(f"       c[1] = {log_c[1]:.17e}  (expect ~0.6667)")

    # --- Sanity gate: abort if c[0] is wildly wrong ---
    print("\n---------------------------------------------------------")
    print("  SANITY CHECK")
    print("---------------------------------------------------------")
    ok = True
    for name, c0, expect in [("sin", sin_c[0], 1.0), ("cos", cos_c[0], 1.0),
                              ("exp", exp_c[0], 1.0), ("log", log_c[0], 2.0)]:
        if abs(c0 - expect) > 0.01:
            print(f"  FAIL: {name} c[0] = {c0:.6f}, expected ~{expect}")
            ok = False
        else:
            print(f"  OK:   {name} c[0] = {c0:.6f}")
    if not ok:
        print("\n  Coefficients are garbage. Aborting before writing files.")
        return 1

    # --- Validate against mpmath ---
    print("\n---------------------------------------------------------")
    print("  VALIDATION (mpmath ground truth, 80 dps)")
    print("---------------------------------------------------------")

    test_sin = np.linspace(-PI_OVER_4, PI_OVER_4, 10000).tolist()
    test_cos = np.linspace(-PI_OVER_4, PI_OVER_4, 10000).tolist()
    test_exp = np.linspace(-LN2_OVER_2, LN2_OVER_2, 10000).tolist()
    test_log = np.linspace(1e-15, Z_MAX, 10000).tolist()

    results = []
    _, line = validate("sin", sin_c, sin, test_sin, eval_odd_poly)
    results.append(line)
    _, line = validate("cos", cos_c, cos, test_cos, eval_even_poly)
    results.append(line)
    _, line = validate("exp", exp_c, exp, test_exp, eval_poly)
    results.append(line)
    _, line = validate("log", log_c, log_atanh, test_log, eval_odd_poly)
    results.append(line)

    for r in results:
        print(f"  {r}")

    # --- Write ---
    print("\n---------------------------------------------------------")
    print("  WRITING OUTPUTS")
    print("---------------------------------------------------------")
    write_text(header_path, generate_header(sin_c, cos_c, exp_c, log_c))
    write_text(report_path, generate_report(results))

    # --- Archive ---
    try:
        import shutil
        source = Path(__file__).resolve()
        dest = v12 / "scripts" / "v12a1" / source.name
        if source != dest:
            dest.parent.mkdir(parents=True, exist_ok=True)
            if not dest.exists() or force:
                shutil.copy2(source, dest)
                print(f"  [archive] {dest}")
    except NameError:
        pass

    print("\n---------------------------------------------------------")
    all_pass = all("PASS" in r for r in results)
    if all_pass:
        print("  ALL FUNCTIONS PASS (<= 5 ULP)")
    else:
        print("  WARNING: Some functions exceed 5 ULP.")
    print("")
    print("  Next: 02_error_free_cleanup.py")
    print("=========================================================")
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
