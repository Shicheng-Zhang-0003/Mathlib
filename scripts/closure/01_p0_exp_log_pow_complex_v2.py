#!/usr/bin/env python3
"""
01_p0_exp_log_pow_complex_v2.py

Run from the folder that CONTAINS the v11S working folder.

This is the robust replacement for 01_p0_exp_log_pow_complex.py.

It applies the P0 closure fixes using regex-based function replacement
instead of brittle exact-string anchors.

Targets:

    v11S/src/internal/pow_util.h
    v11S/src/exp_log.c
    v11S/src/complex.c
    v11S/docs/CLOSURE_P0_APPLIED.md

Usage:

    python3 01_p0_exp_log_pow_complex_v2.py
    python3 01_p0_exp_log_pow_complex_v2.py --force
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path


POW_UTIL_H = r"""#ifndef ML_INTERNAL_POW_UTIL_H
#define ML_INTERNAL_POW_UTIL_H

#include "ml_core.h"

/*
 * MATHLIB_CLOSURE_P0:
 * Shared exponent-classification helpers for real and complex power paths.
 */

static inline int ml_is_integer_double(double y) {
    if (ml_isnan(y) || ml_isinf(y)) return 0;
    if (y == 0.0) return 1;

    /*
     * At or beyond 2^53, all representable finite doubles are integral.
     */
    if (ml_fabs(y) >= 9007199254740992.0) return 1;

    return ml_round(y) == y;
}

static inline int ml_is_odd_integer_double(double y) {
    long long yi;
    long long ay;

    if (!ml_is_integer_double(y)) return 0;

    /*
     * At or beyond 2^53, representable integers are multiples of two
     * or larger powers of two, so they cannot be odd.
     */
    if (ml_fabs(y) >= 9007199254740992.0) return 0;

    yi = (long long)y;
    ay = (yi < 0) ? -yi : yi;

    return (int)(ay % 2LL);
}

#endif /* ML_INTERNAL_POW_UTIL_H */
"""


NEW_EXP = r"""ML_API double ml_exp(double x) {
    /* MATHLIB_CLOSURE_P0_EXP_GUARD */
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return (x > 0.0) ? ml_make_inf(0) : 0.0;

    if (x == 0.0) return 1.0;
    if (x > 709.78) return ml_make_inf(0);
    if (x < -745.13) return 0.0;

    double n = ml_round(x / ML_LN2);
    double r = x - n * 0.69314718036912381649 - n * 1.90821490974462528503e-10;

    static const double inv_fact[] = {
        1.0, 1.0, 0.5, 0.16666666666666666, 0.041666666666666664,
        0.008333333333333333, 0.001388888888888889, 0.0001984126984126984,
        2.48015873015873e-05, 2.7557319223985893e-06, 2.7557319223985888e-07,
        2.505210838544172e-08, 2.08767569878681e-09, 1.6059043836821613e-10,
        1.1470745597729725e-11, 7.647163731819816e-13, 4.779477332387385e-14,
        2.8114572543455206e-15, 1.5619206968586226e-16, 8.22063524662433e-18
    };

    double result = inv_fact[19];
    for (int i = 18; i >= 1; i--) result = ML_FMA(result, r, inv_fact[i]);
    result = ML_FMA(result, r, 1.0);

    return ml_ldexp_pure(result, (int)n);
}

"""


NEW_LOG = r"""ML_API double ml_log(double x) {
    /* MATHLIB_CLOSURE_P0_LOG_GUARD */
    if (ml_isnan(x)) return x;
    if (x == 0.0) return -ml_make_inf(0);
    if (x < 0.0) return ml_make_nan();
    if (ml_isinf(x)) return x;
    if (x == 1.0) return 0.0;

    int e;
    double m = ml_frexp_pure(x, &e);

    int adjust = (m < 0.7071067811865475);
    m *= (1.0 + adjust);
    e -= adjust;

    double z = (m - 1.0) / (m + 1.0);
    double z2 = z * z;

    double poly = 0.09523809523809523;
    poly = poly * z2 + 0.10526315789473684;
    poly = poly * z2 + 0.11764705882352941;
    poly = poly * z2 + 0.13333333333333333;
    poly = poly * z2 + 0.15384615384615385;
    poly = poly * z2 + 0.18181818181818182;
    poly = poly * z2 + 0.2222222222222222;
    poly = poly * z2 + 0.2857142857142857;
    poly = poly * z2 + 0.4;
    poly = poly * z2 + 0.6666666666666666;
    poly = poly * z2 + 2.0;

    return z * poly + e * ML_LN2;
}

"""


NEW_POW = r"""ML_API double ml_pow(double x, double y) {
    /* MATHLIB_CLOSURE_P0_POW_TREE */
    if (ml_isnan(y)) {
        if (x == 1.0) return 1.0;
        return ml_make_nan();
    }

    if (y == 0.0) return 1.0;
    if (ml_isnan(x)) return ml_make_nan();
    if (x == 1.0) return 1.0;

    if (x == 0.0) {
        if (ml_isinf(y)) {
            return (y > 0.0) ? 0.0 : ml_make_inf(0);
        }

        if (y > 0.0) {
            if (ml_signbit(x) && ml_is_odd_integer_double(y)) {
                return ml_copysign(0.0, -1.0);
            }
            return 0.0;
        }

        if (ml_signbit(x) && ml_is_odd_integer_double(y)) {
            return -ml_make_inf(0);
        }
        return ml_make_inf(0);
    }

    if (ml_isinf(y)) {
        double ax = ml_fabs(x);
        if (ax == 1.0) return 1.0;
        if (y > 0.0) return (ax > 1.0) ? ml_make_inf(0) : 0.0;
        return (ax > 1.0) ? 0.0 : ml_make_inf(0);
    }

    if (ml_isinf(x)) {
        if (x > 0.0) {
            return (y > 0.0) ? ml_make_inf(0) : 0.0;
        }

        if (!ml_is_integer_double(y)) return ml_make_nan();

        if (y > 0.0) {
            return ml_is_odd_integer_double(y) ? -ml_make_inf(0) : ml_make_inf(0);
        }

        return ml_is_odd_integer_double(y) ? ml_copysign(0.0, -1.0) : 0.0;
    }

    if (x < 0.0) {
        if (!ml_is_integer_double(y)) return ml_make_nan();

        double ax = -x;
        double mag = ml_exp(y * ml_log(ax));

        if (ml_is_odd_integer_double(y)) return -mag;
        return mag;
    }

    return ml_exp(y * ml_log(x));
}

"""


NEW_SINH = r"""ML_API double ml_sinh(double x) {
    /* MATHLIB_CLOSURE_P0_SINH_SMALL */
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return x;

    double ax = ml_fabs(x);

    if (ax < 1e-4) return x;

    if (ax > 709.782712893384) {
        return ml_make_inf(x < 0.0);
    }

    double ep = ml_exp(ax);
    double em = ml_exp(-ax);
    double r = 0.5 * (ep - em);

    return (x < 0.0) ? -r : r;
}

"""


NEW_ASINH = r"""ML_API double ml_asinh(double x) {
    /* MATHLIB_CLOSURE_P0_ASINH_LARGE */
    if (ml_isnan(x) || ml_isinf(x)) return x;

    double ax = ml_fabs(x);
    if (ax == 0.0) return x;
    if (ax < 1e-4) return x;

    if (ax > 1e150) {
        double r = ml_log(2.0) + ml_log(ax);
        return (x < 0.0) ? -r : r;
    }

    double r = ml_log(ax + ml_hypot_internal(ax, 1.0));
    return (x < 0.0) ? -r : r;
}

"""


NEW_COMPLEX_ARG = r"""ML_API double ml_cplx_arg(cplx a) {
    /* MATHLIB_CLOSURE_P0_COMPLEX_ARG_ATAN2 */
    return ml_atan2(a.imag, a.real);
}

"""


NEW_COMPLEX_POWER = r"""/* MATHLIB_CLOSURE_P0_COMPLEX_POWER_TREE */
static int ml_cplx_is_nan(cplx z) {
    return ml_isnan(z.real) || ml_isnan(z.imag);
}

static int ml_cplx_is_zero(cplx z) {
    return z.real == 0.0 && z.imag == 0.0;
}

static int ml_cplx_is_one(cplx z) {
    return z.real == 1.0 && z.imag == 0.0;
}

static int ml_cplx_is_pure_real_inf(cplx z) {
    return ml_isinf(z.real) && z.imag == 0.0;
}

static cplx ml_cplx_make_nan(void) {
    double n = ml_make_nan();
    return (cplx){n, n};
}

static cplx ml_cplx_make_zero(void) {
    return (cplx){0.0, 0.0};
}

static cplx ml_cplx_make_one(void) {
    return (cplx){1.0, 0.0};
}

static cplx ml_cplx_make_pos_inf(void) {
    return (cplx){ml_make_inf(0), 0.0};
}

static cplx ml_cplx_make_neg_inf(void) {
    return (cplx){-ml_make_inf(0), 0.0};
}

static cplx ml_cplx_make_neg_zero(void) {
    return (cplx){ml_copysign(0.0, -1.0), 0.0};
}

ML_API cplx ml_cplx_power(cplx a, cplx b) {
    int a_nan = ml_cplx_is_nan(a);
    int b_nan = ml_cplx_is_nan(b);

    int a_zero = ml_cplx_is_zero(a);
    int a_one = ml_cplx_is_one(a);
    int b_zero = ml_cplx_is_zero(b);

    int a_inf = ml_isinf(a.real) || ml_isinf(a.imag);
    int b_inf = ml_isinf(b.real) || ml_isinf(b.imag);

    if (b_zero) {
        return ml_cplx_make_one();
    }

    if (a_one) {
        return ml_cplx_make_one();
    }

    if (a_nan) {
        return ml_cplx_make_nan();
    }

    if (b_nan) {
        return ml_cplx_make_nan();
    }

    if (a_zero) {
        if (b_inf) {
            if (ml_isinf(b.imag)) {
                return ml_cplx_make_nan();
            }

            if (b.real > 0.0) {
                return ml_cplx_make_zero();
            }

            if (b.real < 0.0) {
                return ml_cplx_make_pos_inf();
            }

            return ml_cplx_make_nan();
        }

        if (b.real > 0.0) {
            return ml_cplx_make_zero();
        }

        if (b.real < 0.0) {
            return ml_cplx_make_pos_inf();
        }

        return ml_cplx_make_nan();
    }

    if (a_inf) {
        if (b_inf) {
            return ml_cplx_make_nan();
        }

        if (!ml_cplx_is_pure_real_inf(a)) {
            return ml_cplx_make_nan();
        }

        if (b.imag != 0.0) {
            return ml_cplx_make_nan();
        }

        if (b.real > 0.0) {
            if (a.real < 0.0 && ml_is_odd_integer_double(b.real)) {
                return ml_cplx_make_neg_inf();
            }
            return ml_cplx_make_pos_inf();
        }

        if (b.real < 0.0) {
            if (a.real < 0.0 && ml_is_odd_integer_double(b.real)) {
                return ml_cplx_make_neg_zero();
            }
            return ml_cplx_make_zero();
        }

        return ml_cplx_make_one();
    }

    if (b_inf) {
        return ml_cplx_make_nan();
    }

    cplx log_a = ml_cplx_logarithm(a);

    return ml_cplx_exponential((cplx){
        b.real * log_a.real - b.imag * log_a.imag,
        b.real * log_a.imag + b.imag * log_a.real
    });
}
"""


P0_APPLIED_MD = r"""# MathLib v11S P0 Source Fixes Applied

This file records that the P0 closure source-fix script has been applied.

Script:

    01_p0_exp_log_pow_complex_v2.py

Applied fixes:

- added `src/internal/pow_util.h`
- hardened `ml_exp`
- hardened `ml_log`
- rewrote `ml_pow`
- hardened `ml_sinh`
- hardened `ml_asinh`
- replaced `ml_cplx_arg` with atan2-based logic
- rewrote `ml_cplx_power` with explicit special-case handling

Next validation step:

    cd v11S
    ./closure_gate.sh
"""


PAT_EXP = (
    r"(?ms)^[ \t]*ML_API[ \t]+double[ \t]+ml_exp\(double[ \t]+x\)[ \t]*\{.*?"
    r"(?=^[ \t]*ML_API[ \t]+double[ \t]+ml_log\(|\Z)"
)

PAT_LOG = (
    r"(?ms)^[ \t]*ML_API[ \t]+double[ \t]+ml_log\(double[ \t]+x\)[ \t]*\{.*?"
    r"(?=^[ \t]*ML_API[ \t]+double[ \t]+ml_pow\(|\Z)"
)

PAT_POW = (
    r"(?ms)^[ \t]*ML_API[ \t]+double[ \t]+ml_pow\(double[ \t]+x,[ \t]*double[ \t]+y\)[ \t]*\{.*?"
    r"(?=^[ \t]*ML_API[ \t]+double[ \t]+ml_logb\(|\Z)"
)

PAT_SINH = (
    r"(?ms)^[ \t]*ML_API[ \t]+double[ \t]+ml_sinh\(double[ \t]+x\)[ \t]*\{.*?"
    r"(?=^[ \t]*ML_API[ \t]+double[ \t]+ml_cosh\(|\Z)"
)

PAT_ASINH = (
    r"(?ms)^[ \t]*ML_API[ \t]+double[ \t]+ml_asinh\(double[ \t]+x\)[ \t]*\{.*?"
    r"(?=^[ \t]*ML_API[ \t]+double[ \t]+ml_acosh\(|\Z)"
)

PAT_COMPLEX_ARG = (
    r"(?ms)^[ \t]*ML_API[ \t]+double[ \t]+ml_cplx_arg\(cplx[ \t]+a\)[ \t]*\{.*?"
    r"(?=^[ \t]*/\*|^[ \t]*ML_API[ \t]+cplx[ \t]+ml_cplx_div\(|\Z)"
)

PAT_COMPLEX_POWER_FROM_MARKER = (
    r"(?ms)^/\* MATHLIB_CLOSURE_P0_COMPLEX_POWER_TREE \*/.*\Z"
)

PAT_COMPLEX_POWER_FROM_SIGNATURE = (
    r"(?ms)^[ \t]*ML_API[ \t]+cplx[ \t]+ml_cplx_power\(cplx[ \t]+a,[ \t]*cplx[ \t]+b\)[ \t]*\{.*\Z"
)


def fail(message: str) -> None:
    print("ERROR: " + message)
    sys.exit(1)


def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(content)
    print(f"[write] {path}")


def normalize_text(text: str) -> str:
    return text.replace("\r\n", "\n").replace("\r", "\n")


def ensure_pow_util_include(text: str, filename: str) -> tuple[str, bool]:
    if '#include "internal/pow_util.h"' in text:
        print(f"[skip] {filename} already includes internal/pow_util.h")
        return text, False

    pattern = re.compile(r'(?m)^(#include "internal/hypot.h".*)$')
    new_text, count = pattern.subn(
        lambda m: m.group(1) + '\n#include "internal/pow_util.h"',
        text,
        count=1,
    )

    if count != 1:
        fail(f"{filename}: could not find internal/hypot.h include anchor.")

    print(f"[patch] {filename} include")
    return new_text, True


def patch_function(
    text: str,
    pattern: str,
    replacement: str,
    marker: str,
    description: str,
    force: bool,
) -> tuple[str, bool]:
    if marker in text and not force:
        print(f"[skip] {description}: marker already present")
        return text, False

    regex = re.compile(pattern)
    new_text, count = regex.subn(lambda m: replacement, text, count=1)

    if count != 1:
        print(f"[error] {description}: regex match count = {count}, expected 1.")
        print("This usually means the source file has changed unexpectedly.")
        sys.exit(1)

    print(f"[patch] {description}")
    return new_text, True


def patch_complex_power(text: str, force: bool) -> tuple[str, bool]:
    marker = "MATHLIB_CLOSURE_P0_COMPLEX_POWER_TREE"

    if marker in text and not force:
        print("[skip] ml_cplx_power special-case decision tree: marker already present")
        return text, False

    if marker in text:
        pattern = PAT_COMPLEX_POWER_FROM_MARKER
    else:
        pattern = PAT_COMPLEX_POWER_FROM_SIGNATURE

    regex = re.compile(pattern)
    new_text, count = regex.subn(lambda m: NEW_COMPLEX_POWER, text, count=1)

    if count != 1:
        print("[error] ml_cplx_power special-case decision tree: regex match count = "
              f"{count}, expected 1.")
        sys.exit(1)

    print("[patch] ml_cplx_power special-case decision tree")
    return new_text, True


def patch_exp_log(v11s: Path, force: bool) -> None:
    path = v11s / "src" / "exp_log.c"
    original = path.read_text(encoding="utf-8")
    text = normalize_text(original)
    changed = text != original

    text, ch = ensure_pow_util_include(text, "exp_log.c")
    changed = changed or ch

    text, ch = patch_function(
        text,
        PAT_EXP,
        NEW_EXP,
        "MATHLIB_CLOSURE_P0_EXP_GUARD",
        "ml_exp NaN/Inf guard",
        force,
    )
    changed = changed or ch

    text, ch = patch_function(
        text,
        PAT_LOG,
        NEW_LOG,
        "MATHLIB_CLOSURE_P0_LOG_GUARD",
        "ml_log NaN/Inf guard",
        force,
    )
    changed = changed or ch

    text, ch = patch_function(
        text,
        PAT_POW,
        NEW_POW,
        "MATHLIB_CLOSURE_P0_POW_TREE",
        "ml_pow special-case decision tree",
        force,
    )
    changed = changed or ch

    text, ch = patch_function(
        text,
        PAT_SINH,
        NEW_SINH,
        "MATHLIB_CLOSURE_P0_SINH_SMALL",
        "ml_sinh tiny-input branch",
        force,
    )
    changed = changed or ch

    text, ch = patch_function(
        text,
        PAT_ASINH,
        NEW_ASINH,
        "MATHLIB_CLOSURE_P0_ASINH_LARGE",
        "ml_asinh large-input overflow guard",
        force,
    )
    changed = changed or ch

    if changed:
        write_text(path, text)
    else:
        print(f"[skip] {path} already patched")


def patch_complex(v11s: Path, force: bool) -> None:
    path = v11s / "src" / "complex.c"
    original = path.read_text(encoding="utf-8")
    text = normalize_text(original)
    changed = text != original

    text, ch = ensure_pow_util_include(text, "complex.c")
    changed = changed or ch

    text, ch = patch_function(
        text,
        PAT_COMPLEX_ARG,
        NEW_COMPLEX_ARG,
        "MATHLIB_CLOSURE_P0_COMPLEX_ARG_ATAN2",
        "ml_cplx_arg atan2 replacement",
        force,
    )
    changed = changed or ch

    text, ch = patch_complex_power(text, force)
    changed = changed or ch

    if changed:
        write_text(path, text)
    else:
        print(f"[skip] {path} already patched")


def main() -> int:
    force = "--force" in sys.argv[1:]

    root = Path.cwd()
    v11s = root / "v11S"

    if not v11s.is_dir():
        fail("Run this script from the folder that CONTAINS the v11S directory.")

    print("=========================================================")
    print("  MathLib v11S P0 Source Fixes (Robust v2)")
    print("=========================================================")
    print(f"Root: {root}")
    print(f"v11S: {v11s}")
    print(f"Force: {force}")
    print("---------------------------------------------------------")

    pow_util_path = v11s / "src" / "internal" / "pow_util.h"
    if pow_util_path.exists() and not force:
        print(f"[skip] {pow_util_path} already exists")
    else:
        write_text(pow_util_path, POW_UTIL_H)

    patch_exp_log(v11s, force)
    patch_complex(v11s, force)

    applied_doc = v11s / "docs" / "CLOSURE_P0_APPLIED.md"
    if applied_doc.exists() and not force:
        print(f"[skip] {applied_doc} already exists")
    else:
        write_text(applied_doc, P0_APPLIED_MD)

    try:
        source_script = Path(__file__).resolve()
        archived_script = v11s / "scripts" / "closure" / "01_p0_exp_log_pow_complex_v2.py"

        if source_script != archived_script:
            if archived_script.exists() and not force:
                print(f"[skip] {archived_script} already exists")
            else:
                archived_script.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source_script, archived_script)
                print(f"[archive] {archived_script}")
    except NameError:
        print("[note] Could not archive script because __file__ is unavailable.")

    print("---------------------------------------------------------")
    print("P0 source fixes applied.")
    print("")
    print("Next step:")
    print("    cd v11S && ./closure_gate.sh")
    print("")
    print("If the gate passes, the next script should be:")
    print("    02_p1_simd_index_ci_docs.py")
    print("=========================================================")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
