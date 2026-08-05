#!/usr/bin/env python3
"""
24_a1_gamma_stirling.py

Run from the folder that CONTAINS the v12A1 working folder.

A1 gamma/lgamma closure fix.

Replaces the Lanczos-based positive-domain gamma core with a
double-double Stirling expansion plus upward recurrence for small
positive arguments.

Why:
    gamma(x) = exp(lgamma(x))

    A few ULP of absolute error in lgamma is acceptable for lgamma,
    but becomes hundreds of ULP of relative error after exp.

    For gamma <= 5 ULP, lgamma must be accurate to roughly 1e-15
    absolute or better near x = 171.

New strategy:
    - x >= 8: double-double Stirling expansion.
    - 0 < x < 8: recur upward to y = x + m >= 8:
          lgamma(x) = lgamma(y) - sum_{k=0}^{m-1} log(x + k)
    - Negative x: reflection formula with double-double log term.
    - Half-integer negative reflection uses |sin(pi x)| = 1 exactly.

Targets:
    v12A1/src/integral.c

Usage:
    python3 24_a1_gamma_stirling.py
    python3 24_a1_gamma_stirling.py --force
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

MARKER = "MATHLIB_V12A1_GAMMA_STIRLING"
OLD_MARKER = "/* MATHLIB_V12A1_GAMMA_DD2 */"

NEW_GAMMA_SECTION = r'''/* MATHLIB_V12A1_GAMMA_STIRLING */
/*
* Gamma / log-gamma via double-double Stirling expansion.
*
* The old Lanczos path could pass lgamma oracle entries while still
* carrying enough absolute lgamma error to produce hundreds of ULP
* of error after exp() for gamma.
*
* For gamma <= 5 ULP, the internal lgamma value must be accurate to
* about 1e-15 absolute or better near x = 171.
*
* Strategy:
*   - x >= 8: double-double Stirling expansion.
*   - 0 < x < 8: recur upward to y = x + m >= 8:
*         lgamma(x) = lgamma(y) - sum log(x + k)
*   - negative x: reflection formula:
*         log|Gamma(x)| = log(pi / |sin(pi x)|) - lgamma(1 - x)
*   - half-integer negative reflection uses |sin(pi x)| = 1 exactly.
*/

#ifndef ML_LN2_HI
#define ML_LN2_HI 6.93147180369123816490e-01
#endif
#ifndef ML_LN2_LO
#define ML_LN2_LO 1.90821492927058500170e-10
#endif

#define ML_HALF_LOG_2PI 0.91893853320467274178
#define ML_GAMMA_OVERFLOW 171.6243769563027
#define ML_GAMMA_EXP_OVERFLOW 709.782712893384
#define ML_GAMMA_EXP_UNDERFLOW (-745.133219101941)

/* ------------------------------------------------------------------ */
/* Double-double primitives. A value is represented as hi + lo.        */
/* ------------------------------------------------------------------ */

typedef struct {
    double hi;
    double lo;
} ml_dd_t;

static inline ml_dd_t ml_dd_from_d(double a) {
    ml_dd_t r;
    r.hi = a;
    r.lo = 0.0;
    return r;
}

/* Exact a + b (Knuth two-sum). */
static inline ml_dd_t ml_dd_two_sum(double a, double b) {
    double s = a + b;
    double v = s - a;
    ml_dd_t r;
    r.hi = s;
    r.lo = (a - (s - v)) + (b - v);
    return r;
}

/* Renormalize hi + lo into canonical double-double form. */
static inline ml_dd_t ml_dd_renorm(double hi, double lo) {
    double s = hi + lo;
    ml_dd_t r;
    r.hi = s;
    r.lo = lo - (s - hi);
    return r;
}

static inline ml_dd_t ml_dd_add(ml_dd_t a, ml_dd_t b) {
    ml_dd_t s = ml_dd_two_sum(a.hi, b.hi);
    return ml_dd_renorm(s.hi, s.lo + a.lo + b.lo);
}

static inline ml_dd_t ml_dd_sub(ml_dd_t a, ml_dd_t b) {
    ml_dd_t nb;
    nb.hi = -b.hi;
    nb.lo = -b.lo;
    return ml_dd_add(a, nb);
}

static inline ml_dd_t ml_dd_mul(ml_dd_t a, ml_dd_t b) {
    double p = a.hi * b.hi;
    double e = ML_FMA(a.hi, b.hi, -p) + (a.hi * b.lo + a.lo * b.hi);
    return ml_dd_renorm(p, e);
}

/* ------------------------------------------------------------------ */
/* ~106-bit natural logarithm for finite x > 0.                       */
/* ------------------------------------------------------------------ */

static ml_dd_t ml_log_dd(double x) {
    if (ml_isnan(x) || x <= 0.0) {
        return ml_dd_from_d(ml_make_nan());
    }
    if (ml_isinf(x)) {
        return ml_dd_from_d(x);
    }
    if (x == 1.0) {
        return ml_dd_from_d(0.0);
    }

    int e;
    double m = ml_frexp_pure(x, &e);
    int adjust = (m < 0.7071067811865475);
    m *= (1.0 + (double)adjust);
    e -= adjust;

    /* z = (m - 1) / (m + 1) in double-double. */
    double num = m - 1.0;
    ml_dd_t den = ml_dd_two_sum(m, 1.0);

    double q = num / den.hi;
    double r = ML_FMA(-q, den.hi, num);
    r = ML_FMA(-q, den.lo, r);
    double q2 = r / den.hi;
    ml_dd_t z = ml_dd_renorm(q, q2);

    ml_dd_t z2 = ml_dd_mul(z, z);

    /* log((1+z)/(1-z)) = z * P(z^2), P(u) = sum 2/(2k+1) u^k. */
    static const double lc[11] = {
        2.0,
        0.6666666666666666,
        0.4,
        0.2857142857142857,
        0.2222222222222222,
        0.18181818181818182,
        0.15384615384615385,
        0.13333333333333333,
        0.11764705882352941,
        0.10526315789473684,
        0.09523809523809523
    };

    ml_dd_t p = ml_dd_from_d(lc[10]);
    for (int i = 9; i >= 0; i--) {
        p = ml_dd_add(ml_dd_mul(p, z2), ml_dd_from_d(lc[i]));
    }

    ml_dd_t lm = ml_dd_mul(z, p);

    /* e * ln2 in double-double. */
    double ed = (double)e;
    double ehi = ed * ML_LN2_HI;
    double elo = ML_FMA(ed, ML_LN2_HI, -ehi) + ed * ML_LN2_LO;
    ml_dd_t eln2 = ml_dd_renorm(ehi, elo);

    return ml_dd_add(lm, eln2);
}

/* ------------------------------------------------------------------ */
/* Double-double Stirling expansion for x >= 8.                       */
/* ------------------------------------------------------------------ */

static ml_dd_t ml_lgamma_stirling_dd(double x) {
    ml_dd_t lx = ml_log_dd(x);
    ml_dd_t w = ml_dd_from_d(x - 0.5);

    /* (x - 0.5) * log(x) - x + 0.5*log(2*pi) */
    ml_dd_t prod = ml_dd_mul(w, lx);
    ml_dd_t L = ml_dd_sub(prod, ml_dd_from_d(x));
    L = ml_dd_add(L, ml_dd_from_d(ML_HALF_LOG_2PI));

    /*
    * Asymptotic correction:
    *
    *   1/(12x)
    * - 1/(360 x^3)
    * + 1/(1260 x^5)
    * - 1/(1680 x^7)
    * + 1/(1188 x^9)
    * - 691/(360360 x^11)
    * + 1/(156 x^13)
    * - 3617/(122400 x^15)
    * + 43867/(244188 x^17)
    *
    * For x >= 8, the first omitted term is already far below
    * the required absolute lgamma budget.
    */
    double invx = 1.0 / x;
    double invx2 = invx * invx;

    double corr = invx * (1.0 / 12.0);
    double p = invx * invx2;

    corr -= p * (1.0 / 360.0);
    p *= invx2;
    corr += p * (1.0 / 1260.0);
    p *= invx2;
    corr -= p * (1.0 / 1680.0);
    p *= invx2;
    corr += p * (1.0 / 1188.0);
    p *= invx2;
    corr -= p * (691.0 / 360360.0);
    p *= invx2;
    corr += p * (1.0 / 156.0);
    p *= invx2;
    corr -= p * (3617.0 / 122400.0);
    p *= invx2;
    corr += p * (43867.0 / 244188.0);

    L = ml_dd_add(L, ml_dd_from_d(corr));
    return L;
}

/* ------------------------------------------------------------------ */
/* Positive-domain lgamma in double-double.                           */
/* ------------------------------------------------------------------ */

static ml_dd_t ml_lgamma_positive_dd(double x) {
    if (x >= 8.0) {
        return ml_lgamma_stirling_dd(x);
    }

    /*
    * Recur upward until y = x + m >= 8:
    *
    *   Gamma(y) = Gamma(x) * x * (x+1) * ... * (x+m-1)
    *
    * so:
    *
    *   lgamma(x) = lgamma(y) - sum_{k=0}^{m-1} log(x+k)
    */
    double t = 8.0 - x;
    int m = (int)t;
    if ((double)m < t) {
        m++;
    }
    if (m < 1) {
        m = 1;
    }

    double y = x + (double)m;
    ml_dd_t L = ml_lgamma_stirling_dd(y);

    for (int k = 0; k < m; k++) {
        double v = x + (double)k;
        ml_dd_t lv = ml_log_dd(v);
        L = ml_dd_sub(L, lv);
    }

    return L;
}

/* ------------------------------------------------------------------ */
/* exp(hi + lo) using first-order low-part correction.                */
/* ------------------------------------------------------------------ */

static double ml_exp_dd(ml_dd_t L) {
    if (L.hi > ML_GAMMA_EXP_OVERFLOW) {
        return ml_make_inf(0);
    }
    if (L.hi < ML_GAMMA_EXP_UNDERFLOW) {
        return 0.0;
    }

    double g = ml_exp(L.hi);
    if (!ml_isfinite(g) || g == 0.0) {
        return g;
    }

    return ML_FMA(g, L.lo, g);
}

/* ------------------------------------------------------------------ */
/* Small exact helpers.                                               */
/* ------------------------------------------------------------------ */

/* (n-1)! for 1 <= n <= 23. */
static double ml_factorial_exact_small(int n) {
    double f = 1.0;
    for (int k = 2; k <= n - 1; k++) {
        f *= (double)k;
    }
    return f;
}

static int ml_is_half_integer(double x) {
    if (!ml_isfinite(x)) {
        return 0;
    }
    double t = 2.0 * x;
    if (!ml_isfinite(t)) {
        return 0;
    }
    return t == ml_round(t);
}

/*
* For x = n + 0.5, sin(pi x) = (-1)^n.
*/
static double ml_half_sin_sign(double x) {
    if (!ml_isfinite(x)) {
        return 1.0;
    }

    double n = ml_round(x - 0.5);
    if (ml_fabs(n) < 9007199254740992.0) {
        long long ni = (long long)n;
        if ((ni % 2LL) != 0LL) {
            return -1.0;
        }
    }

    return 1.0;
}

/* ------------------------------------------------------------------ */
/* Public APIs                                                         */
/* ------------------------------------------------------------------ */

ML_API double ml_lgamma(double x) {
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return ml_make_inf(0);

    /* Poles at zero and negative integers. */
    if (x <= 0.0 && x == ml_round(x)) {
        return ml_make_inf(0);
    }

    /* Exact zeros. */
    if (x == 1.0 || x == 2.0) {
        return 0.0;
    }

    if (x > 0.0) {
        /* Exact small integer cases: lgamma(n) = log((n-1)!). */
        if (x == ml_round(x) && x <= 23.0) {
            return ml_log(ml_factorial_exact_small((int)x));
        }

        ml_dd_t L = ml_lgamma_positive_dd(x);
        return L.hi + L.lo;
    }

    /* Reflection: log|Gamma(x)| = log(pi/|sin(pi x)|) - lgamma(1-x). */
    double sinpx = ml_sin(ML_PI * x);
    if (sinpx == 0.0) {
        return ml_make_inf(0);
    }

    double absin = ml_fabs(sinpx);
    if (ml_is_half_integer(x)) {
        absin = 1.0;
    }
    if (absin == 0.0) {
        return ml_make_inf(0);
    }

    ml_dd_t logterm = ml_log_dd(ML_PI / absin);
    ml_dd_t Lpos = ml_lgamma_positive_dd(1.0 - x);
    ml_dd_t r = ml_dd_sub(logterm, Lpos);

    return r.hi + r.lo;
}

ML_API double ml_gamma_new(double x) {
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return x > 0.0 ? ml_make_inf(0) : ml_make_nan();

    /* Poles at zero and negative integers. */
    if (x <= 0.0 && x == ml_round(x)) {
        return ml_make_nan();
    }

    if (x > 0.0) {
        if (x > ML_GAMMA_OVERFLOW) {
            return ml_make_inf(0);
        }

        /* Exact small integer cases: gamma(n) = (n-1)!. */
        if (x == ml_round(x) && x <= 23.0) {
            return ml_factorial_exact_small((int)x);
        }

        ml_dd_t L = ml_lgamma_positive_dd(x);
        return ml_exp_dd(L);
    }

    /* Reflection: Gamma(x) = pi / (sin(pi x) * Gamma(1-x)). */
    double sinpx = ml_sin(ML_PI * x);
    if (sinpx == 0.0) {
        return ml_make_nan();
    }

    double absin = ml_fabs(sinpx);
    double sinsign = ml_signbit(sinpx) ? -1.0 : 1.0;

    if (ml_is_half_integer(x)) {
        absin = 1.0;
        sinsign = ml_half_sin_sign(x);
    }

    if (absin == 0.0) {
        return ml_make_nan();
    }

    ml_dd_t Lpos = ml_lgamma_positive_dd(1.0 - x);
    double G = ml_exp_dd(Lpos);

    if (ml_isinf(G)) {
        return ml_copysign(0.0, sinsign);
    }
    if (G == 0.0) {
        return (sinsign < 0.0) ? -ml_make_inf(0) : ml_make_inf(0);
    }

    double denom = sinsign * absin * G;
    return ML_PI / denom;
}
'''


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
        print(f"  [skip] {path}: already patched with {MARKER}")
        return

    idx = text.find(OLD_MARKER)
    if idx < 0:
        fail(
            f"{path}: could not find old gamma marker {OLD_MARKER}. "
            "Source may have drifted."
        )

    new_text = text[:idx] + NEW_GAMMA_SECTION
    write_text(path, new_text)


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
    print("  MATHLIB v12A1: GAMMA/LGAMMA STIRLING DD FIX")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")

    print("\n[1/1] integral.c — replace gamma core with DD Stirling")
    patch_integral(v12, force)

    archive_self(v12, force)

    print("\n---------------------------------------------------------")
    print("  Gamma Stirling DD fix applied.")
    print("")
    print("  What changed:")
    print("    - x >= 8 uses double-double Stirling expansion")
    print("    - 0 < x < 8 recurs upward to x+m >= 8")
    print("    - negative reflection uses double-double log term")
    print("    - half-integer reflection uses |sin(pi x)| = 1")
    print("    - exact integer shortcuts retained for n <= 23")
    print("")
    print("  Rebuild and verify:")
    print("")
    print("    cd v12A1")
    print("    cmake --build build")
    print("    ./build/oracle_check > /tmp/oracle_out3.txt || true")
    print("    grep -n FAIL /tmp/oracle_out3.txt")
    print("    tail -n 30 /tmp/oracle_out3.txt")
    print("")
    print("  Also run edge tests:")
    print("")
    print("    MATHLIB_EDGE_SANITIZERS=1 bash tests/run_edge_tests.sh")
    print("")
    print("  Expected: gamma/lgamma failures collapse toward <= 5 ULP.")
    print("=========================================================")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
