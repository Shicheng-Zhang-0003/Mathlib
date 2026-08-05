#!/usr/bin/env python3
"""
25_a1_gamma_residue.py

Run from the folder that CONTAINS the v12A1 working folder.

Final gamma/lgamma residue fix for the three remaining oracle failures:

    gamma(0.001)  ~10 ULP
    gamma(-0.9)    ~6 ULP
    lgamma(-2.5)  ~11 ULP

Changes:
    - Positive small-x gamma uses product reconstruction:
          gamma(x) = gamma(x+m) / prod_{k=0}^{m-1}(x+k)
      with double-double product. This avoids exponentiating a
      subtracted log sum for small x.

    - lgamma still uses double-double log recurrence.

    - Reflection uses:
          log|gamma(x)| = log(pi) - log|sin(pi x)| - lgamma(1-x)
      with a double-double log(pi) term.

    - sin(pi x) is computed with a compensated pi argument:
          arg = x*pi_hi + x*pi_lo
      before calling ml_sin().

    - Half-integer positive values use exact half-integer structure:
          Gamma(n+1/2) = sqrt(pi) * prod_{j=0}^{n-1}(j+1/2)
          lgamma(n+1/2) = 0.5*log(pi) + sum log(j+1/2)

Targets:
    v12A1/src/integral.c

Usage:
    python3 25_a1_gamma_residue.py
    python3 25_a1_gamma_residue.py --force
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

NEW_MARKER = "MATHLIB_V12A1_GAMMA_STIRLING_V2"

OLD_MARKERS = [
    "MATHLIB_V12A1_GAMMA_STIRLING",
    "MATHLIB_V12A1_GAMMA_DD2",
    "MATHLIB_V12A1_GAMMA_LANCZOS",
]

NEW_GAMMA_SECTION = r'''/* MATHLIB_V12A1_GAMMA_STIRLING_V2 */
/*
* Gamma / log-gamma closure layer.
*
* Strategy:
*   - x >= 8: double-double Stirling expansion.
*   - 0 < x < 8:
*       lgamma(x) = lgamma(x+m) - sum log(x+k)
*       gamma(x)  = gamma(x+m) / prod(x+k)
*     where x+m >= 8.
*   - Half-integer positive values use exact half-integer structure.
*   - Negative reflection uses double-double log(pi) and compensated sin(pi x).
*
* This targets the final residue:
*   gamma(0.001), gamma(-0.9), lgamma(-2.5).
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

/*
* Double-double pi.
*
* ML_PI_HI_D is the nearest double to pi.
* ML_PI_LO_D is the positive residual true_pi - double_pi.
*/
#define ML_PI_HI_D 0x1.921fb54442d18p+1
#define ML_PI_LO_D 0x1.1a62633145c07p-53

/* ------------------------------------------------------------------ */
/* Double-double primitives.                                           */
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

/* Renormalize hi + lo. Uses two-sum so it is safe even after cancellation. */
static inline ml_dd_t ml_dd_renorm(double hi, double lo) {
    return ml_dd_two_sum(hi, lo);
}

static inline ml_dd_t ml_dd_add(ml_dd_t a, ml_dd_t b) {
    ml_dd_t s = ml_dd_two_sum(a.hi, b.hi);
    return ml_dd_renorm(s.hi, s.lo + a.lo + b.lo);
}

static inline ml_dd_t ml_dd_add_d(ml_dd_t a, double b) {
    return ml_dd_add(a, ml_dd_from_d(b));
}

static inline ml_dd_t ml_dd_sub(ml_dd_t a, ml_dd_t b) {
    ml_dd_t nb;
    nb.hi = -b.hi;
    nb.lo = -b.lo;
    return ml_dd_add(a, nb);
}

static inline ml_dd_t ml_dd_mul_d(ml_dd_t a, double b) {
    double p = a.hi * b;
    double e = ML_FMA(a.hi, b, -p) + a.lo * b;
    return ml_dd_renorm(p, e);
}

static inline ml_dd_t ml_dd_mul(ml_dd_t a, ml_dd_t b) {
    double p = a.hi * b.hi;
    double e = ML_FMA(a.hi, b.hi, -p) + (a.hi * b.lo + a.lo * b.hi);
    return ml_dd_renorm(p, e);
}

/* ------------------------------------------------------------------ */
/* Double-double logarithm.                                            */
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

    double num = m - 1.0;
    ml_dd_t den = ml_dd_two_sum(m, 1.0);

    double q = num / den.hi;
    double r = ML_FMA(-q, den.hi, num);
    r = ML_FMA(-q, den.lo, r);
    double q2 = r / den.hi;
    ml_dd_t z = ml_dd_renorm(q, q2);

    ml_dd_t z2 = ml_dd_mul(z, z);

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

    double ed = (double)e;
    double ehi = ed * ML_LN2_HI;
    double elo = ML_FMA(ed, ML_LN2_HI, -ehi) + ed * ML_LN2_LO;
    ml_dd_t eln2 = ml_dd_renorm(ehi, elo);

    return ml_dd_add(lm, eln2);
}

/* Double-double log(pi). */
static ml_dd_t ml_log_pi_dd(void) {
    ml_dd_t lp = ml_log_dd(ML_PI_HI_D);
    return ml_dd_add_d(lp, ML_PI_LO_D / ML_PI_HI_D);
}

/* ------------------------------------------------------------------ */
/* exp(hi + lo).                                                       */
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
/* Stirling expansion for x >= 8.                                      */
/* ------------------------------------------------------------------ */

static ml_dd_t ml_stirling_lgamma_dd(double x) {
    if (!ml_isfinite(x)) {
        return ml_dd_from_d(x);
    }
    if (x <= 0.0) {
        return ml_dd_from_d(ml_make_nan());
    }

    static const double sc[10] = {
        1.0 / 12.0,
        -1.0 / 360.0,
        1.0 / 1260.0,
        -1.0 / 1680.0,
        1.0 / 1188.0,
        -691.0 / 360360.0,
        1.0 / 156.0,
        -3617.0 / 122400.0,
        43867.0 / 244188.0,
        -174611.0 / 125400.0
    };

    double invx = 1.0 / x;
    double invx2 = invx * invx;

    double corr = invx * sc[0];
    double p = invx * invx2;
    for (int i = 1; i < 10; i++) {
        corr += p * sc[i];
        p *= invx2;
    }

    ml_dd_t lx = ml_log_dd(x);
    ml_dd_t w = ml_dd_from_d(x - 0.5);
    ml_dd_t prod = ml_dd_mul(w, lx);

    ml_dd_t L = ml_dd_sub(prod, ml_dd_from_d(x));
    L = ml_dd_add_d(L, ML_HALF_LOG_2PI);
    L = ml_dd_add_d(L, corr);

    return L;
}

/* ------------------------------------------------------------------ */
/* Half-integer helpers.                                               */
/* ------------------------------------------------------------------ */

static int ml_is_half_integer(double x) {
    if (!ml_isfinite(x)) {
        return 0;
    }
    /* Integers are not half-integers. */
    if (x == ml_round(x)) {
        return 0;
    }
    double t = 2.0 * x;
    if (!ml_isfinite(t)) {
        return 0;
    }
    return t == ml_round(t);
}

/*
* For x = n + 0.5, returns n.
* Only safe for small/moderate half-integers; callers use it for x < 8.
*/
static long long ml_half_index(double x) {
    double nd = ml_round(x - 0.5);
    if (!ml_isfinite(nd)) {
        return 0;
    }
    if (ml_fabs(nd) >= 9007199254740992.0) {
        return 0;
    }
    return (long long)nd;
}

/* sin(pi * (n + 0.5)) = (-1)^n. */
static double ml_half_sin_sign(double x) {
    long long n = ml_half_index(x);
    if ((n % 2LL) != 0LL) {
        return -1.0;
    }
    return 1.0;
}

static ml_dd_t ml_lgamma_half_positive_dd(double x) {
    long long n = ml_half_index(x);

    /* 0.5 * log(pi) */
    ml_dd_t L = ml_dd_mul_d(ml_log_pi_dd(), 0.5);

    for (long long j = 0; j < n; j++) {
        double v = (double)j + 0.5;
        L = ml_dd_add(L, ml_log_dd(v));
    }

    return L;
}

static double ml_gamma_half_positive(double x) {
    long long n = ml_half_index(x);

    /* sqrt(pi) = exp(0.5 * log(pi)) */
    ml_dd_t L = ml_dd_mul_d(ml_log_pi_dd(), 0.5);
    double g = ml_exp_dd(L);

    for (long long j = 0; j < n; j++) {
        g *= ((double)j + 0.5);
    }

    return g;
}

/* ------------------------------------------------------------------ */
/* Positive-domain lgamma and gamma.                                   */
/* ------------------------------------------------------------------ */

static ml_dd_t ml_lgamma_positive_dd(double x) {
    if (x >= 8.0) {
        return ml_stirling_lgamma_dd(x);
    }

    if (ml_is_half_integer(x)) {
        return ml_lgamma_half_positive_dd(x);
    }

    double t = 8.0 - x;
    int m = (int)t;
    if ((double)m < t) {
        m++;
    }
    if (m < 1) {
        m = 1;
    }

    double y = x + (double)m;
    ml_dd_t L = ml_stirling_lgamma_dd(y);

    for (int k = 0; k < m; k++) {
        double v = x + (double)k;
        L = ml_dd_sub(L, ml_log_dd(v));
    }

    return L;
}

static double ml_gamma_positive(double x) {
    if (x >= 8.0) {
        ml_dd_t L = ml_stirling_lgamma_dd(x);
        return ml_exp_dd(L);
    }

    if (ml_is_half_integer(x)) {
        return ml_gamma_half_positive(x);
    }

    double t = 8.0 - x;
    int m = (int)t;
    if ((double)m < t) {
        m++;
    }
    if (m < 1) {
        m = 1;
    }

    double y = x + (double)m;

    ml_dd_t Ly = ml_stirling_lgamma_dd(y);
    double gy = ml_exp_dd(Ly);

    /* Product P = x * (x+1) * ... * (x+m-1) in double-double. */
    ml_dd_t P = ml_dd_from_d(1.0);
    for (int k = 0; k < m; k++) {
        double v = x + (double)k;
        P = ml_dd_mul_d(P, v);
    }

    double p = P.hi + P.lo;
    if (p == 0.0) {
        return ml_make_inf(0);
    }

    return gy / p;
}

/* ------------------------------------------------------------------ */
/* Exact small factorials.                                             */
/* ------------------------------------------------------------------ */

static double ml_factorial_exact_small(int n) {
    double f = 1.0;
    for (int k = 2; k <= n - 1; k++) {
        f *= (double)k;
    }
    return f;
}

/* ------------------------------------------------------------------ */
/* Compensated sin(pi x).                                              */
/* ------------------------------------------------------------------ */

static double ml_sinpi(double x) {
    if (ml_isnan(x) || ml_isinf(x)) {
        return ml_make_nan();
    }
    if (x == 0.0) {
        return x;
    }

    double arg_hi = x * ML_PI_HI_D;
    double arg_lo = ML_FMA(x, ML_PI_HI_D, -arg_hi) + x * ML_PI_LO_D;
    double arg = arg_hi + arg_lo;

    if (!ml_isfinite(arg)) {
        return ml_sin(ML_PI * x);
    }

    return ml_sin(arg);
}

/* ------------------------------------------------------------------ */
/* Public APIs.                                                        */
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
        /* lgamma(n) = log((n-1)!) via exact factorial. */
        if (x == ml_round(x) && x <= 23.0) {
            return ml_log(ml_factorial_exact_small((int)x));
        }

        ml_dd_t L = ml_lgamma_positive_dd(x);
        return L.hi + L.lo;
    }

    /* Reflection: log|Gamma(x)| = log(pi) - log|sin(pi x)| - lgamma(1-x). */
    double s = ml_sinpi(x);
    if (s == 0.0) {
        return ml_make_inf(0);
    }

    double absin = ml_fabs(s);
    if (ml_is_half_integer(x)) {
        absin = 1.0;
    }
    if (absin == 0.0) {
        return ml_make_inf(0);
    }

    ml_dd_t logterm = ml_log_pi_dd();
    if (absin != 1.0) {
        logterm = ml_dd_sub(logterm, ml_log_dd(absin));
    }

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

        /* Exact factorial shortcut. */
        if (x == ml_round(x) && x <= 23.0) {
            return ml_factorial_exact_small((int)x);
        }

        return ml_gamma_positive(x);
    }

    /* Reflection: Gamma(x) = pi / (sin(pi x) * Gamma(1-x)). */
    double s = ml_sinpi(x);
    if (s == 0.0) {
        return ml_make_nan();
    }

    double sin_use = s;
    if (ml_is_half_integer(x)) {
        sin_use = ml_half_sin_sign(x);
    }
    if (sin_use == 0.0) {
        return ml_make_nan();
    }

    double G = ml_gamma_positive(1.0 - x);

    if (ml_isinf(G)) {
        return ml_copysign(0.0, sin_use);
    }
    if (G == 0.0) {
        return (sin_use < 0.0) ? -ml_make_inf(0) : ml_make_inf(0);
    }

    return ML_PI_HI_D / (sin_use * G);
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

    if NEW_MARKER in text and not force:
        print(f"  [skip] {path}: already at {NEW_MARKER}")
        return

    search_markers = [NEW_MARKER] + OLD_MARKERS
    start = None

    for marker in search_markers:
        pattern = re.compile(r"/\*\s*" + re.escape(marker) + r"\s*\*/")
        match = pattern.search(text)
        if match:
            start = match.start()
            break

    if start is None:
        fail(
            f"{path}: could not find gamma section marker. "
            "Expected one of: " + ", ".join(search_markers)
        )

    new_text = text[:start] + NEW_GAMMA_SECTION
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
    print("  MATHLIB v12A1: GAMMA/LGAMMA FINAL RESIDUE FIX")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")

    print("\n[1/1] integral.c — gamma residue v2")
    patch_integral(v12, force)

    archive_self(v12, force)

    print("\n---------------------------------------------------------")
    print("  Gamma residue fix applied.")
    print("")
    print("  Key changes:")
    print("    - small positive gamma uses DD product reconstruction")
    print("    - half-integers use exact half-integer structure")
    print("    - reflection uses DD log(pi)")
    print("    - sin(pi x) uses compensated pi argument")
    print("")
    print("  Rebuild and verify:")
    print("")
    print("    cd v12A1")
    print("    cmake --build build")
    print("    ./build/oracle_check > /tmp/oracle_out4.txt || true")
    print("    grep -n FAIL /tmp/oracle_out4.txt")
    print("    tail -n 30 /tmp/oracle_out4.txt")
    print("")
    print("  Edge tests:")
    print("")
    print("    MATHLIB_EDGE_SANITIZERS=1 bash tests/run_edge_tests.sh")
    print("")
    print("  Expected: oracle and edge suites collapse to zero failures.")
    print("=========================================================")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
