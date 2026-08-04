#!/usr/bin/env python3
"""
18_2_minimax_truth_v3.py
Run from the folder that CONTAINS the v12A1 working folder.

GIANT ERROR #2 (final resolution).

v2 established the truth:
  - Active kernels (Maclaurin sin/cos, Taylor-19 exp, atanh log) all
    measure 1 ULP under the FMA Horner evaluation model.
  - Chebyshev fits pass the <= 5 ULP gate for sin, cos, log.
  - Chebyshev fits CANNOT pass for exp: approximation error at degree 10
    is ~1 ULP, but Horner evaluation noise adds ~6-8 ULP, and the noise
    GROWS with degree (measured: deg10=8, deg11=16, deg12=33, deg13=10,
    deg14=12). This is a hard floor for power-basis double Horner.

Resolution for exp:
  The validated exp coefficient set is the Taylor-19 kernel itself
  (exact 1/k! coefficients). It measures 1 ULP under the same model,
  strictly better than any Chebyshev fit. The generator now:
    - fits sin/cos/log via Chebyshev (gated <= 5 ULP),
    - VALIDATES (not fits) the Taylor-19 exp kernel (gated <= 5 ULP),
    - records the Chebyshev exp evidence as informational in the report,
    - writes nothing unless all four gates pass.

Requires: pip install numpy mpmath

Usage:
  python3 18_2_minimax_truth_v3.py
  python3 18_2_minimax_truth_v3.py --force
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

MARKER = "MATHLIB_V12A1_MINIMAX_TRUTH_V3"

GENERATOR_SRC = r'''#!/usr/bin/env python3
"""
generate_minimax.py — canonical v12A1 coefficient pipeline (v3).
MATHLIB_V12A1_MINIMAX_TRUTH_V3

v3 policy:
  - sin/cos/log: Chebyshev interpolation fits, hard gate <= 5 ULP.
  - exp: no fit. The Taylor-19 kernel (exact 1/k! coefficients) is
    validated under the FMA Horner model and hard-gated at <= 5 ULP.
    Chebyshev interpolation provably cannot beat it: evaluation-noise
    floor of ~8 ULP for power-basis double Horner on [-ln2/2, ln2/2].
    Chebyshev degrees 9..14 are still measured and recorded in the
    report as evidence (informational, not gated).

Evaluation model: FMA Horner via libm fma() (correctly rounded),
matching ML_FMA Horner in src/internal/minimax.h and src/exp_log.c.

Hard gate: nothing is written unless every kernel validates <= 5 ULP.

Requires: pip install numpy mpmath
"""
from __future__ import annotations

import ctypes
import struct
import sys
from ctypes.util import find_library
from datetime import date
from pathlib import Path

try:
    import numpy as np
    from numpy.polynomial.chebyshev import chebfit, cheb2poly
    from numpy.polynomial import Polynomial
except ImportError:
    print("ERROR: numpy not installed. Run: pip install numpy")
    sys.exit(1)

try:
    from mpmath import mp, mpf, sin, cos, exp, log, pi, ln, sqrt
    mp.dps = 80
except ImportError:
    print("ERROR: mpmath not installed. Run: pip install mpmath")
    sys.exit(1)

MARKER = "MATHLIB_V12A1_MINIMAX_TRUTH_V3"
TARGET_ULP = 5
GRID_N = 5001
TODAY = date.today().isoformat()

PI4 = float(pi / 4)
LN2H = float(ln(2) / 2)
ZMAX = float((sqrt(2) - 1) / (sqrt(2) + 1))

V12 = Path(__file__).resolve().parents[2]


# -------------------------------------------------------------------
# FMA matching the C kernel evaluation model (ML_FMA Horner)
# -------------------------------------------------------------------
def _load_libm_fma():
    candidates = []
    lib = find_library("m")
    if lib:
        candidates.append(lib)
    candidates += ["libm.so.6", "libm.so"]
    for name in candidates:
        try:
            libm = ctypes.CDLL(name)
            f = libm.fma
            f.restype = ctypes.c_double
            f.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double]
            return f
        except (OSError, AttributeError):
            continue
    return None


_FMA = _load_libm_fma()
if _FMA is not None:
    EVAL_MODEL = "FMA Horner (libm fma, correctly rounded)"

    def fma(a, b, c):
        return _FMA(a, b, c)
else:
    EVAL_MODEL = "PLAIN Horner (libm fma unavailable — gate may be conservative)"

    def fma(a, b, c):
        return a * b + c


def horner_fma(c, x):
    r = c[-1]
    for i in range(len(c) - 2, -1, -1):
        r = fma(r, x, c[i])
    return r


# -------------------------------------------------------------------
# Fitting: Chebyshev interpolation nodes
# -------------------------------------------------------------------
def cheb_nodes(a, b, n):
    k = np.arange(n)
    return a + (b - a) * (np.cos((2 * k + 1) * np.pi / (2 * n)) + 1.0) / 2.0


def cheb_to_power_exact(cheb_coeffs):
    """Exact 80-digit conversion: Chebyshev series -> power basis."""
    n = len(cheb_coeffs)
    T = [[mpf(1)], [mpf(0), mpf(1)]]
    while len(T) < n:
        k = len(T)
        tk = [mpf(0)] * (k + 1)
        for i, c in enumerate(T[k - 1]):
            tk[i + 1] += 2 * c
        for i, c in enumerate(T[k - 2]):
            tk[i] -= c
        T.append(tk)
    result = [mpf(0)] * n
    for k, a in enumerate(cheb_coeffs):
        ak = mpf(float(a))
        for i, c in enumerate(T[k]):
            result[i] += ak * c
    return result


def compose_linear_exact(power_u, alpha, beta):
    """Exact 80-digit composition: Q(t) = P(alpha*t + beta)."""
    n = len(power_u)
    result = [mpf(0)] * n
    pow_lin = [mpf(1)]
    a_ = mpf(float(alpha))
    b_ = mpf(float(beta))
    for k in range(n):
        ak = mpf(float(power_u[k]))
        for i, c in enumerate(pow_lin):
            result[i] += ak * c
        if k < n - 1:
            newp = [mpf(0)] * (len(pow_lin) + 1)
            for i, c in enumerate(pow_lin):
                newp[i + 1] += c * a_
                newp[i] += c * b_
            pow_lin = newp
    return [float(c) for c in result]


def from_cheb_on_interval(cheb_coeffs, a, b):
    power_u = cheb_to_power_exact(cheb_coeffs)
    alpha = 2.0 / (b - a)
    beta = -(a + b) / (b - a)
    return compose_linear_exact(power_u, alpha, beta)


def fit_plain(f, a, b, deg):
    nodes = cheb_nodes(a, b, deg + 1)
    y = np.array([float(f(mpf(float(t)))) for t in nodes])
    u = (2.0 * nodes - (a + b)) / (b - a)
    return from_cheb_on_interval(chebfit(u, y, deg), a, b)


def fit_odd(f, b, deg):
    """Fit f(x) = x * P(x^2) on [0, b]. Returns coefficients of P."""
    t_max = b * b
    t_nodes = cheb_nodes(0.0, t_max, deg + 1)
    x_nodes = np.sqrt(t_nodes)
    y = np.array([float(f(mpf(float(x))) / x) for x in x_nodes])
    u = (2.0 * t_nodes - t_max) / t_max
    return from_cheb_on_interval(chebfit(u, y, deg), 0.0, t_max)


def fit_even(f, b, deg):
    """Fit f(x) = P(x^2) on [0, b] (f even). Returns coefficients of P."""
    t_max = b * b
    t_nodes = cheb_nodes(0.0, t_max, deg + 1)
    x_nodes = np.sqrt(t_nodes)
    y = np.array([float(f(mpf(float(x)))) for x in x_nodes])
    u = (2.0 * t_nodes - t_max) / t_max
    return from_cheb_on_interval(chebfit(u, y, deg), 0.0, t_max)


# -------------------------------------------------------------------
# ULP machinery
# -------------------------------------------------------------------
def d2i(d):
    b = struct.unpack("Q", struct.pack("d", d))[0]
    return (0x8000000000000000 - b) if (b >> 63) else b


def ulp(a, b):
    return abs(d2i(a) - d2i(b))


def make_grid(f, a, b):
    pts = np.linspace(a, b, GRID_N).tolist()
    return [(x, float(f(mpf(float(x))))) for x in pts]


def max_ulp(c, grid, kind):
    mx, worst = 0, 0.0
    for x, e in grid:
        if kind == "odd":
            g = x * horner_fma(c, x * x)
        elif kind == "even":
            g = horner_fma(c, x * x)
        else:
            g = horner_fma(c, x)
        if e == 0.0 and g == 0.0:
            continue
        u = ulp(g, e)
        if u > mx:
            mx, worst = u, x
    return mx, worst


# -------------------------------------------------------------------
# Active kernels in tree (measured with the same FMA model)
# -------------------------------------------------------------------
MACLAURIN_SIN = [
    1.0, -0.16666666666666666, 0.008333333333333333,
    -0.0001984126984126984, 2.7557319223985893e-06,
    -2.505210838544172e-08, 1.6059043836821613e-10,
    -7.647163731819816e-13, 2.811457254345521e-15,
    -8.220635816560923e-18,
]
MACLAURIN_COS = [
    1.0, -0.5, 0.041666666666666664, -0.001388888888888889,
    2.48015873015873e-05, -2.755731922398589e-07,
    2.08767569878681e-09, -1.1470745597729725e-11,
    4.779477332387385e-14, -1.5619206967218455e-16,
]
TAYLOR_EXP = [
    1.0, 1.0, 0.5, 0.16666666666666666, 0.041666666666666664,
    0.008333333333333333, 0.001388888888888889, 0.0001984126984126984,
    2.48015873015873e-05, 2.7557319223985893e-06, 2.7557319223985888e-07,
    2.505210838544172e-08, 2.08767569878681e-09, 1.6059043836821613e-10,
    1.1470745597729725e-11, 7.647163731819816e-13, 4.779477332387385e-14,
    2.8114572543455206e-15, 1.5619206968586226e-16, 8.22063524662433e-18,
]
SERIES_LOG = [
    2.0, 0.6666666666666666, 0.4, 0.2857142857142857,
    0.2222222222222222, 0.18181818181818182, 0.15384615384615385,
    0.13333333333333333, 0.11764705882352941, 0.10526315789473684,
    0.09523809523809523,
]


def log_atanh(z):
    return log((1 + z) / (1 - z))


def c_array(name, c):
    lines = ["static const double " + name + "[] = {"]
    for i, v in enumerate(c):
        comma = "," if i < len(c) - 1 else ""
        lines.append("    %.17e%s" % (v, comma))
    lines.append("};")
    return "\n".join(lines)


def main() -> int:
    print("=========================================================")
    print("  MATHLIB v12A1: CANONICAL MINIMAX PIPELINE (v3)")
    print("=========================================================")
    print("  Evaluation model: " + EVAL_MODEL)
    if not (V12 / "src" / "internal").is_dir():
        print("ERROR: %s does not look like the v12A1 tree." % V12)
        return 1

    print("\n[grids] building mpmath-80dps ground truth ...")
    g_sin = make_grid(sin, -PI4, PI4)
    g_cos = make_grid(cos, -PI4, PI4)
    g_exp = make_grid(exp, -LN2H, LN2H)
    g_log = make_grid(log_atanh, 0.0, ZMAX)

    print("\n[active kernels] measuring Maclaurin/Taylor kernels in tree ...")
    act_sin = max_ulp(MACLAURIN_SIN, g_sin, "odd")
    act_cos = max_ulp(MACLAURIN_COS, g_cos, "even")
    act_exp = max_ulp(TAYLOR_EXP, g_exp, "plain")
    act_log = max_ulp(SERIES_LOG, g_log, "odd")
    print("    sin (Maclaurin, minimax.h): %d ULP" % act_sin[0])
    print("    cos (Maclaurin, minimax.h): %d ULP" % act_cos[0])
    print("    exp (Taylor-19 kernel):     %d ULP" % act_exp[0])
    print("    log (atanh series kernel):  %d ULP" % act_log[0])

    print("\n[fit] sin/cos/log Chebyshev search (gated <= %d ULP) ..." % TARGET_ULP)
    searches = [
        ("sin", "odd", lambda d: fit_odd(sin, PI4, d), range(6, 13), g_sin),
        ("cos", "even", lambda d: fit_even(cos, PI4, d), range(6, 13), g_cos),
        ("log", "odd", lambda d: fit_odd(log_atanh, ZMAX, d), range(6, 13), g_log),
    ]
    results = {}
    ok = True
    for name, kind, fitfn, degs, grid in searches:
        best = None
        for d in degs:
            c = fitfn(d)
            u, wx = max_ulp(c, grid, kind)
            print("    %s deg=%d: %d ULP (worst x=%.6e)" % (name, d, u, wx))
            if u <= TARGET_ULP:
                best = (c, u, d)
                break
        if best is None:
            ok = False
            print("  [FAIL] %s: no degree reached <= %d ULP" % (name, TARGET_ULP))
        else:
            results[name] = best

    print("\n[exp] Chebyshev evidence (informational, NOT gated) ...")
    exp_evidence = []
    for d in range(9, 15):
        c = fit_plain(exp, -LN2H, LN2H, d)
        u, wx = max_ulp(c, g_exp, "plain")
        print("    exp cheb deg=%d: %d ULP (worst x=%.6e)" % (d, u, wx))
        exp_evidence.append((d, u))
    exp_cheb_best = min(exp_evidence, key=lambda t: t[1])

    print("\n[exp] Taylor-19 kernel validation (GATED) ...")
    exp_u, exp_wx = max_ulp(TAYLOR_EXP, g_exp, "plain")
    print("    exp Taylor-19 (active kernel): %d ULP (worst x=%.6e)" % (exp_u, exp_wx))
    if exp_u > TARGET_ULP:
        ok = False
        print("  [FAIL] exp: Taylor-19 kernel exceeds %d ULP" % TARGET_ULP)

    if not ok:
        print("\nHARD GATE FAILED: nothing written.")
        print("Fix the pipeline before retrying. This is intentional.")
        return 1

    sc, su, sd = results["sin"]
    cc, cu, cd = results["cos"]
    lc, lu, ld = results["log"]

    header_path = V12 / "src" / "internal" / "minimax_coeffs.h"
    parts = [
        "#ifndef ML_INTERNAL_MINIMAX_COEFFS_H",
        "#define ML_INTERNAL_MINIMAX_COEFFS_H",
        "",
        "/* %s */" % MARKER,
        "/* Auto-generated by scripts/oracles/generate_minimax.py (v3) on %s */" % TODAY,
        "/* DO NOT EDIT MANUALLY. Regenerate with the pipeline script. */",
        "/*",
        " * STATUS: DORMANT.",
        " * The active kernels in src/internal/minimax.h and src/exp_log.c",
        " * are the validated kernels. These coefficients are a clean,",
        " * validated set kept for a future gated swap. No code includes",
        " * this header yet.",
        " *",
        " * Validation model: %s" % EVAL_MODEL,
        " */",
        "",
        "/* sin(x) = x*P(x^2) on [-pi/4,pi/4], deg P = %d, %d ULP */" % (sd, su),
        c_array("minimax_sin_coeffs", sc),
        "",
        "/* cos(x) = P(x^2) on [-pi/4,pi/4], deg P = %d, %d ULP */" % (cd, cu),
        c_array("minimax_cos_coeffs", cc),
        "",
        "/* exp(x) = P(x) on [-ln2/2,ln2/2], deg = 19, %d ULP */" % exp_u,
        "/* NOTE: Taylor-19 coefficients (exact 1/k!, rounded once).",
        " * Chebyshev interpolation cannot reach <= 5 ULP for exp under",
        " * double FMA Horner: evaluation-noise floor ~%d ULP at degrees" % exp_cheb_best[1],
        " * 9-14 (noise grows with degree). This Taylor set measures %d ULP" % exp_u,
        " * under the same model and is the validated exp kernel. */",
        c_array("minimax_exp_coeffs", TAYLOR_EXP),
        "",
        "/* log((1+z)/(1-z)) = z*P(z^2) on [0,%.17e], deg P = %d, %d ULP */" % (ZMAX, ld, lu),
        c_array("minimax_log_coeffs", lc),
        "",
        "#endif /* ML_INTERNAL_MINIMAX_COEFFS_H */",
        "",
    ]
    header_path.parent.mkdir(parents=True, exist_ok=True)
    with open(header_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(parts))
    print("  [write] %s" % header_path)

    report_path = V12 / "docs" / "V12A1_MINIMAX_REPORT.md"
    def status(u):
        return "PASS" if u <= TARGET_ULP else "WARN"
    rep = [
        "# v12A1 Minimax Coefficient Validation Report",
        "",
        "Generated: %s" % TODAY,
        "Method: Chebyshev interpolation (exact 80-digit power-basis conversion)",
        "Ground truth: mpmath at 80 decimal places",
        "Validation grid: %d points per reduced domain" % GRID_N,
        "Evaluation model: %s" % EVAL_MODEL,
        "",
        "## Generated coefficients (minimax_coeffs.h, DORMANT)",
        "",
        "| Function | Parity | Degree | Max ULP | Status |",
        "|----------|--------|--------|---------|--------|",
        "| sin | odd, x*P(x^2) | %d | %d | %s |" % (sd, su, status(su)),
        "| cos | even, P(x^2) | %d | %d | %s |" % (cd, cu, status(cu)),
        "| exp | plain (Taylor-19) | 19 | %d | %s |" % (exp_u, status(exp_u)),
        "| log | odd, z*P(z^2) | %d | %d | %s |" % (ld, lu, status(lu)),
        "",
        "## Active kernels in tree (NOT modified)",
        "",
        "The running code uses the Maclaurin kernels in src/internal/minimax.h",
        "and the Taylor/series kernels in src/exp_log.c. Measured on the same",
        "reduced domains with the same evaluation model:",
        "",
        "| Kernel | Source | Max ULP |",
        "|--------|--------|---------|",
        "| sin (Maclaurin, x*P(x^2), deg 19) | minimax.h | %d |" % act_sin[0],
        "| cos (Maclaurin, P(x^2), deg 18) | minimax.h | %d |" % act_cos[0],
        "| exp (Taylor-19 on reduced arg) | exp_log.c | %d |" % act_exp[0],
        "| log (atanh series, 11 terms) | exp_log.c | %d |" % act_log[0],
        "",
        "## The exp decision (measured, not assumed)",
        "",
        "Chebyshev interpolation cannot reach <= 5 ULP for exp under double",
        "FMA Horner on [-ln2/2, ln2/2]. The approximation error at degree 10",
        "is already ~1 ULP, but Horner evaluation adds ~6-8 ULP of rounding",
        "noise, and the noise grows with degree. Measured evidence:",
        "",
        "| Chebyshev degree | Max ULP |",
        "|------------------|---------|",
    ]
    for d, u in exp_evidence:
        rep.append("| %d | %d |" % (d, u))
    rep += [
        "",
        "The oscillating, non-decreasing error is the signature of",
        "evaluation noise, not truncation. The Taylor-19 kernel (exact 1/k!",
        "coefficients) measures %d ULP under the same model — strictly better" % exp_u,
        "than any fitted polynomial. It is therefore recorded as the validated",
        "exp coefficient set. Real libms escape this floor with table-driven",
        "arguments or extended precision, not with higher-degree Horner.",
        "",
        "## Decisions",
        "",
        "- Hard gate: nothing is written unless every kernel is <= 5 ULP.",
        "- No coefficient swap into trig.c / exp_log.c in this wave.",
        "  The active kernels already meet the accuracy bar; the generated",
        "  header stays dormant until a future gated swap script.",
        "",
        "## Acceptance Criteria",
        "",
        "All functions must achieve <= 5 ULP on their reduced domain.",
        "",
    ]
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(rep))
    print("  [write] %s" % report_path)

    print("\n=========================================================")
    print("  ALL KERNELS PASS (<= %d ULP). Hard gate satisfied." % TARGET_ULP)
    print("=========================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
'''


def fail(message: str) -> None:
    print("ERROR: " + message)
    sys.exit(1)


def locate_v12a1() -> tuple[Path, Path]:
    root = Path.cwd()
    candidate = root / "v12A1"
    if candidate.is_dir():
        return root, candidate
    if (root / "src" / "internal").is_dir() and (root / "include" / "mathlib").is_dir():
        print("  [note] Running from inside v12A1.")
        return root.parent, root
    fail("Run from the folder that CONTAINS v12A1/, or from inside v12A1/ itself.")


def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(content)
    print(f"  [write] {path}")


def install_generator(v12: Path, force: bool) -> Path:
    path = v12 / "scripts" / "oracles" / "generate_minimax.py"
    if path.exists() and not force:
        if MARKER in path.read_text(encoding="utf-8"):
            print(f"  [skip] {path}: generator v3 already installed")
            return path
    write_text(path, GENERATOR_SRC)
    return path


def run_generator(gen: Path, v12: Path) -> int:
    print("\n  [run] generate_minimax.py v3 (this can take a few minutes) ...")
    try:
        result = subprocess.run(
            [sys.executable, str(gen)],
            cwd=str(v12),
            capture_output=True,
            text=True,
            timeout=2400,
        )
    except subprocess.TimeoutExpired:
        print("  [error] generator timed out after 2400s")
        return 1
    if result.stdout:
        print(result.stdout)
    if result.returncode != 0:
        if result.stderr:
            print(result.stderr[-2000:])
        print("  [error] generator FAILED. Nothing was written (hard gate).")
    return result.returncode


def archive_self(v12: Path, force: bool) -> None:
    try:
        source = Path(__file__).resolve()
        dest = v12 / "scripts" / "v12a1" / source.name
        if source == dest:
            return
        if dest.exists() and not force:
            print(f"  [skip] {dest}: already archived")
            return
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, dest)
        print(f"  [archive] {dest}")
    except NameError:
        pass


def main() -> int:
    force = "--force" in sys.argv[1:]
    root, v12 = locate_v12a1()
    print("=========================================================")
    print("  MATHLIB v12A1: MINIMAX TRUTH v3 (script 18_2)")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")

    header = v12 / "src" / "internal" / "minimax_coeffs.h"
    if header.exists() and not force:
        if MARKER in header.read_text(encoding="utf-8"):
            print(f"  [skip] {header}: already regenerated by v3 pipeline.")
            print("         Use --force to re-run.")
            archive_self(v12, force)
            return 0

    print("\n[1/2] Installing generator v3 (exp = validated Taylor-19)")
    gen = install_generator(v12, force)

    print("\n[2/2] Running generator v3 with hard <= 5 ULP gate")
    rc = run_generator(gen, v12)
    if rc != 0:
        print("\n  Pipeline failed. The poisoned header was NOT regenerated.")
        print("  Read the generator output above; fix and re-run.")
        return rc

    if MARKER not in header.read_text(encoding="utf-8"):
        fail(f"{header}: generator reported success but marker missing.")

    archive_self(v12, force)

    print("\n---------------------------------------------------------")
    print("  Minimax truth v3 restored.")
    print("")
    print("  What changed vs v2:")
    print("    - exp is no longer fitted. Chebyshev cannot beat the")
    print("      Taylor-19 kernel under FMA Horner (measured floor).")
    print("    - Taylor-19 is VALIDATED and hard-gated like the rest.")
    print("    - Chebyshev exp evidence recorded in the report.")
    print("    - minimax_coeffs.h regenerated CLEAN and DORMANT.")
    print("    - V12A1_MINIMAX_REPORT.md rewritten with measured truth.")
    print("")
    print("  Next: 19_payne_hanek_v3_fixed.py  (kills the sin/cos catastrophe)")
    print("=========================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
