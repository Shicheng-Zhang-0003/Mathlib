#!/usr/bin/env python3
"""
21_gamma_double_double.py
Run from the folder that CONTAINS the v12A1 working folder.

GIANT ERROR #5: gamma/lgamma amplification (7..544 ULP in the oracle).

Root cause:
  lgamma(x) = 0.5*ln(2pi) + (x-0.5)*ln(x+6.5) - (x+6.5) + ln(Ag(x-1))

  The product (x-0.5)*ln(x+6.5) amplifies the ~1 ULP error of the
  double-precision ln by (x-0.5), up to ~170.5. For gamma(100) this
  turns ~3.6e-15 of log error into ~3.6e-13 relative error: 544 ULP.

Fix: evaluate the entire Lanczos formula in double-double arithmetic
(~106 bits), so intermediate rounding does not survive to the result:

  1. ml_log_dd(): double-double logarithm.
     - z = (m-1)/(m+1) computed as a DD division (m-1 is exact by
       Sterbenz; m+1 captured exactly via two-sum).
     - atanh series evaluated with DD Horner.
     - e*ln2 exact via the ML_LN2_HI / ML_LN2_LO split.
  2. Ag summed in DD, each term c_i/(x+i-1) computed by DD division.
  3. (x-0.5)*log(t) accumulated with FMA error capture.
  4. gamma = exp(Lh) * (1 + Ll) via one FMA (DD argument preserved).

Exact integer special cases (oracle distance to 0.0 is brutal):
  - lgamma(1) = lgamma(2) = 0 exactly
  - gamma(n) = (n-1)! exact for n <= 23 (all exactly representable)
  - lgamma(n) = log((n-1)!) for n <= 23

Targets:
  v12A1/src/integral.c
  v12A1/tests/test_edge_integral.c

Usage:
  python3 21_gamma_double_double.py
  python3 21_gamma_double_double.py --force
"""
from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

MARKER = "MATHLIB_V12A1_GAMMA_DOUBLE_DOUBLE"

NEW_GAMMA_SECTION = r"""/* MATHLIB_V12A1_GAMMA_DOUBLE_DOUBLE */
/*
 * Lanczos approximation (g=7, n=9) evaluated in double-double arithmetic.
 *
 * Why double-double:
 *   lgamma(x) = 0.5*ln(2pi) + (x-0.5)*ln(x+6.5) - (x+6.5) + ln(Ag(x-1))
 *
 * The product (x-0.5)*ln(x+6.5) amplifies any error in ln by (x-0.5),
 * up to ~170.5. A 1 ULP error in ln becomes hundreds of ULP in gamma.
 * Computing ln in double-double (~106 bits) and accumulating the whole
 * expression in double-double removes the amplification.
 *
 * Integer shortcuts:
 *   gamma(n) = (n-1)! is exactly representable in double for n <= 23;
 *   lgamma(1) = lgamma(2) = 0 exactly; lgamma(n) = ln((n-1)!) via the
 *   exact factorial for n <= 23. These bypass the Lanczos approximation
 *   error at integer arguments, where ULP distance to exact oracle
 *   values (0 for lgamma(1), lgamma(2)) would otherwise be astronomic.
 */

static const double ml_lanczos_coeff[9] = {
    0.99999999999980993,
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.61502916214059,
    12.507343278686905,
    -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7
};

#define ML_LANCZOS_G 7.0
#define ML_GAMMA_OVERFLOW 171.6243769563027

#ifndef ML_LN2_HI
#define ML_LN2_HI 6.93147180369123816490e-01
#endif
#ifndef ML_LN2_LO
#define ML_LN2_LO 1.90821492927058500170e-10
#endif

/* ---------------- double-double primitives ---------------- */

/* Knuth two-sum: hi + lo == a + b exactly. */
static inline void ml_dd_two_sum(double a, double b, double *hi, double *lo) {
    double s = a + b;
    double v = s - a;
    *hi = s;
    *lo = (a - (s - v)) + (b - v);
}

/* Fast two-sum. Requires |a| >= |b|. */
static inline void ml_dd_renorm(double a, double b, double *hi, double *lo) {
    double h = a + b;
    *lo = b - (h - a);
    *hi = h;
}

/* DD + double -> DD */
static inline void ml_dd_add_d(double ah, double al, double b,
                               double *hi, double *lo) {
    double h, l;
    ml_dd_two_sum(ah, b, &h, &l);
    l += al;
    ml_dd_renorm(h, l, hi, lo);
}

/* DD + DD -> DD */
static inline void ml_dd_add_dd(double ah, double al, double bh, double bl,
                                double *hi, double *lo) {
    double h, l;
    ml_dd_two_sum(ah, bh, &h, &l);
    l += al + bl;
    ml_dd_renorm(h, l, hi, lo);
}

/* DD - DD -> DD */
static inline void ml_dd_sub_dd(double ah, double al, double bh, double bl,
                                double *hi, double *lo) {
    ml_dd_add_dd(ah, al, -bh, -bl, hi, lo);
}

/* DD * DD -> DD (error ~ few ulp^2) */
static inline void ml_dd_mul_dd(double ah, double al, double bh, double bl,
                                double *hi, double *lo) {
    double p = ah * bh;
    double e = ML_FMA(ah, bh, -p) + (ah * bl + al * bh);
    ml_dd_renorm(p, e, hi, lo);
}

/* double / double -> DD (one Newton correction) */
static inline void ml_dd_div_d(double a, double b, double *hi, double *lo) {
    double q1 = a / b;
    double r = ML_FMA(-q1, b, a);
    double q2 = r / b;
    ml_dd_renorm(q1, q2, hi, lo);
}

/* ---------------- double-double logarithm ---------------- */

/*
 * Horner coefficients for log((1+z)/(1-z)) = z * P(z^2),
 * c_k = 2/(2k+1), ordered for descending Horner in z^2.
 */
static const double ml_atanh_coeff[11] = {
    0.09523809523809523,   /* 2/21 */
    0.10526315789473684,   /* 2/19 */
    0.11764705882352941,   /* 2/17 */
    0.13333333333333333,   /* 2/15 */
    0.15384615384615385,   /* 2/13 */
    0.18181818181818182,   /* 2/11 */
    0.2222222222222222,    /* 2/9  */
    0.2857142857142857,    /* 2/7  */
    0.4,                   /* 2/5  */
    0.6666666666666666,    /* 2/3  */
    2.0                    /* 2/1  */
};

/*
 * ln(x) as double-double. Requires x > 0 and finite.
 */
static void ml_log_dd(double x, double *hi, double *lo) {
    int e;
    double m = ml_frexp_pure(x, &e);
    int adjust = (m < 0.7071067811865475);
    m *= (1.0 + adjust);
    e -= adjust;

    /* z = (m-1)/(m+1) in DD.
     * m-1 is exact by Sterbenz (m in [sqrt(1/2), sqrt(2)]).
     * m+1 is captured exactly by two-sum. */
    double num = m - 1.0;
    double dh, dl;
    ml_dd_two_sum(m, 1.0, &dh, &dl);

    double q1 = num / dh;
    double r = ML_FMA(-q1, dh, num);
    r = ML_FMA(-q1, dl, r);
    double q2 = r / dh;
    double zh, zl;
    ml_dd_renorm(q1, q2, &zh, &zl);

    /* z2 = z*z in DD */
    double z2h, z2l;
    {
        double p = zh * zh;
        double err = ML_FMA(zh, zh, -p) + 2.0 * zh * zl;
        err = ML_FMA(zl, zl, err);
        ml_dd_renorm(p, err, &z2h, &z2l);
    }

    /* P(z^2) via Horner in DD */
    double rh = ml_atanh_coeff[0];
    double rl = 0.0;
    for (int i = 1; i < 11; i++) {
        double ph, pl;
        ml_dd_mul_dd(rh, rl, z2h, z2l, &ph, &pl);
        ml_dd_add_d(ph, pl, ml_atanh_coeff[i], &rh, &rl);
    }

    /* log(m) = z * P(z^2) */
    double lmh, lml;
    ml_dd_mul_dd(rh, rl, zh, zl, &lmh, &lml);

    /* e*ln2: e*ML_LN2_HI is exact (ML_LN2_HI has its low 26
     * significand bits zeroed); e*ML_LN2_LO restores the rest. */
    double eh = (double)e * ML_LN2_HI;
    double el = (double)e * ML_LN2_LO;

    ml_dd_add_dd(lmh, lml, eh, el, hi, lo);
}

/* ---------------- Lanczos core in double-double ---------------- */

/*
 * lgamma(x) for x > 0, returned as double-double (hi + lo).
 */
static void ml_lgamma_positive_dd(double x, double *hi, double *lo) {
    /* Ag = c0 + sum_{i=1..8} c_i / (x + i - 1), accumulated in DD */
    double ah = ml_lanczos_coeff[0];
    double al = 0.0;
    for (int i = 1; i < 9; i++) {
        double denom = x + (double)(i - 1);   /* exact for x < 2^50 */
        double th, tl;
        ml_dd_div_d(ml_lanczos_coeff[i], denom, &th, &tl);
        ml_dd_add_dd(ah, al, th, tl, &ah, &al);
    }

    /* log(Ag) in DD, plus first-order correction for the low part */
    double lah, lal;
    ml_log_dd(ah, &lah, &lal);
    {
        double corr = al / ah;
        double nh, nl;
        ml_dd_add_d(lah, lal, corr, &nh, &nl);
        lah = nh;
        lal = nl;
    }

    /* log(t), t = x + 6.5 (exact for x < 2^50) */
    double t = x + 6.5;
    double lth, ltl;
    ml_log_dd(t, &lth, &ltl);

    /* p = (x - 0.5) * log(t) in DD; x - 0.5 exact for x < 2^51 */
    double w = x - 0.5;
    double ph = w * lth;
    double pl = ML_FMA(w, lth, -ph) + w * ltl;
    double phr, plr;
    ml_dd_renorm(ph, pl, &phr, &plr);

    /* L = 0.5*ln(2*pi) + p - t + log(Ag) */
    double Lh, Ll;
    ml_dd_add_d(phr, plr, -t, &Lh, &Ll);
    ml_dd_add_d(Lh, Ll, 0.91893853320467274178, &Lh, &Ll);
    ml_dd_add_dd(Lh, Ll, lah, lal, &Lh, &Ll);

    *hi = Lh;
    *lo = Ll;
}

/* Exact (n-1)! for 1 <= n <= 23; all such factorials are exactly
 * representable in double. */
static double ml_factorial_exact_small(int n) {
    double f = 1.0;
    for (int k = 2; k <= n - 1; k++) {
        f *= (double)k;
    }
    return f;
}

ML_API double ml_lgamma(double x) {
    /* MATHLIB_V12A1_GAMMA_DOUBLE_DOUBLE */
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return ml_make_inf(0);
    /* Poles at zero and negative integers */
    if (x <= 0.0 && x == ml_round(x)) return ml_make_inf(0);

    /* Exact zeros */
    if (x == 1.0 || x == 2.0) return 0.0;

    if (x > 0.0) {
        /* lgamma(n) = log((n-1)!) via exact factorial */
        if (x == ml_round(x) && x <= 23.0) {
            return ml_log(ml_factorial_exact_small((int)x));
        }
        double Lh, Ll;
        ml_lgamma_positive_dd(x, &Lh, &Ll);
        return Lh + Ll;
    }

    /* Reflection: log|Gamma(x)| = log(pi/|sin(pi x)|) - lgamma(1-x) */
    double sinpx = ml_sin(ML_PI * x);
    if (sinpx == 0.0) return ml_make_inf(0);
    double lh, ll;
    ml_log_dd(ML_PI / ml_fabs(sinpx), &lh, &ll);
    double gh, gl;
    ml_lgamma_positive_dd(1.0 - x, &gh, &gl);
    double rh, rl;
    ml_dd_sub_dd(lh, ll, gh, gl, &rh, &rl);
    return rh + rl;
}

ML_API double ml_gamma_new(double x) {
    /* MATHLIB_V12A1_GAMMA_DOUBLE_DOUBLE */
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return x > 0.0 ? ml_make_inf(0) : ml_make_nan();
    /* Poles at zero and negative integers */
    if (x <= 0.0 && x == ml_round(x)) return ml_make_nan();

    if (x > 0.0) {
        if (x > ML_GAMMA_OVERFLOW) return ml_make_inf(0);
        /* Exact factorial for small integers */
        if (x == ml_round(x) && x <= 23.0) {
            return ml_factorial_exact_small((int)x);
        }
        double Lh, Ll;
        ml_lgamma_positive_dd(x, &Lh, &Ll);
        /* exp(Lh + Ll) = exp(Lh) * (1 + Ll) to first order;
         * the FMA keeps the double-double argument alive. */
        double g = ml_exp(Lh);
        return ML_FMA(g, Ll, g);
    }

    /* Reflection: Gamma(x) = pi / (sin(pi x) * Gamma(1-x)) */
    double sinpx = ml_sin(ML_PI * x);
    if (sinpx == 0.0) return ml_make_nan();
    if (1.0 - x > ML_GAMMA_OVERFLOW) {
        /* Gamma(1-x) overflows, so Gamma(x) -> 0 */
        return ml_copysign(0.0, sinpx);
    }
    double Lh, Ll;
    ml_lgamma_positive_dd(1.0 - x, &Lh, &Ll);
    double g = ml_exp(Lh);
    double gpos = ML_FMA(g, Ll, g);
    return ML_PI / (sinpx * gpos);
}
"""

NEW_EDGE_TESTS = r"""
/* MATHLIB_V12A1_GAMMA_DOUBLE_DOUBLE_TEST */
/* Exact integer gamma values (factorial shortcut) */
ASSERT_TRUE(&ctx, ml_gamma_new(6.0) == 120.0, "gamma(6) == 120 exact");
ASSERT_TRUE(&ctx, ml_gamma_new(7.0) == 720.0, "gamma(7) == 720 exact");
ASSERT_TRUE(&ctx, ml_gamma_new(23.0) == 1.12400072777760768e21, "gamma(23) == 22! exact");
/* lgamma integer identities */
ASSERT_TRUE(&ctx, ml_lgamma(2.0) == 0.0, "lgamma(2) == 0 exact");
ASSERT_TRUE(&ctx, ml_lgamma(1.0) == 0.0, "lgamma(1) == 0 exact");
ASSERT_NEAR(&ctx, ml_lgamma(6.0), ml_log(120.0), 1e-14, "lgamma(6) == log(120)");
/* DD path: large-argument gamma (oracle is the strict judge at <= 5 ULP) */
ASSERT_NEAR(&ctx, ml_gamma_new(171.0), 7.257415615307999e306,
            7.257415615307999e306 * 2e-15, "gamma(171) tight");
ASSERT_NEAR(&ctx, ml_gamma_new(100.0), 9.332621544394415e155,
            9.332621544394415e155 * 2e-15, "gamma(100) tight");
"""


def fail(message: str) -> None:
    print("ERROR: " + message)
    sys.exit(1)


def normalize(text: str) -> str:
    return text.replace("\r\n", "\n").replace("\r", "\n")


def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(content)
    print(f"  [write] {path}")


def locate_v12a1() -> tuple[Path, Path]:
    root = Path.cwd()
    candidate = root / "v12A1"
    if candidate.is_dir():
        return root, candidate
    if (root / "src" / "integral.c").is_file():
        print("  [note] Running from inside v12A1.")
        return root.parent, root
    fail("Run from the folder that CONTAINS v12A1/, or from inside v12A1/ itself.")


def patch_integral(v12: Path, force: bool) -> None:
    path = v12 / "src" / "integral.c"
    if not path.is_file():
        fail(f"Missing expected file: {path}")
    text = normalize(path.read_text(encoding="utf-8"))
    if MARKER in text and not force:
        print(f"  [skip] {path}: already patched")
        return
    # Replace the entire gamma section (from its marker to EOF)
    pattern = re.compile(
        r"(?ms)^[ \t]*/\*[ \t]*MATHLIB_V12A1_GAMMA_(LANCZOS|DOUBLE_DOUBLE)"
        r"[ \t]*\*/.*\Z"
    )
    patched, count = pattern.subn(
        lambda m: NEW_GAMMA_SECTION,
        text,
        count=1,
    )
    if count != 1:
        fail(
            f"{path}: could not find the gamma section marker "
            "(MATHLIB_V12A1_GAMMA_LANCZOS). Source may have drifted."
        )
    write_text(path, patched)


def patch_edge_test(v12: Path, force: bool) -> None:
    path = v12 / "tests" / "test_edge_integral.c"
    if not path.is_file():
        fail(f"Missing expected file: {path}")
    text = normalize(path.read_text(encoding="utf-8"))
    if "MATHLIB_V12A1_GAMMA_DOUBLE_DOUBLE_TEST" in text and not force:
        print(f"  [skip] {path}: already patched")
        return
    anchor = "    return ml_test_summary(&ctx);"
    if anchor not in text:
        fail(f"{path}: could not find return anchor.")
    text = text.replace(anchor, NEW_EDGE_TESTS + "\n" + anchor, 1)
    write_text(path, text)


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
    print("  MATHLIB v12A1: GAMMA/LGAMMA DOUBLE-DOUBLE (script 21)")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")

    print("\n[1/2] integral.c — DD Lanczos + exact integer cases")
    patch_integral(v12, force)

    print("\n[2/2] test_edge_integral.c — tight regression assertions")
    patch_edge_test(v12, force)

    archive_self(v12, force)

    print("\n---------------------------------------------------------")
    print("  Double-double gamma applied.")
    print("")
    print("  What changed:")
    print("    - ml_log_dd(): ~106-bit logarithm (DD division for z,")
    print("      DD Horner for the atanh series, exact e*ln2 split)")
    print("    - Ag summed in DD with DD term divisions")
    print("    - (x-0.5)*log(t) via FMA error capture")
    print("    - gamma = exp(Lh)*(1+Ll) via FMA")
    print("    - exact factorial shortcuts: gamma(n)=(n-1)! for n<=23,")
    print("      lgamma(1)=lgamma(2)=0, lgamma(n)=log((n-1)!)")
    print("")
    print("  Verify:")
    print("    cd v12A1")
    print("    cmake --build build")
    print("    ./build/test")
    print("    MATHLIB_EDGE_SANITIZERS=1 bash tests/run_edge_tests.sh")
    print("    ./build/oracle_check")
    print("")
    print("  Expected: gamma/lgamma collapse from 7..544 ULP toward <= 5.")
    print("  If a few NON-INTEGER points (e.g. gamma(0.5)) linger at")
    print("  ~5-10 ULP, that residue is the intrinsic error of the")
    print("  g=7,n=9 Lanczos coefficient set itself — the next step")
    print("  would be fitting better coefficients with mpmath.")
    print("=========================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
