#!/usr/bin/env python3
"""
07_p0_5_cordic_fuzzer_nonfinite.py

Run from the folder that CONTAINS the v11S working folder.

This fixes Red Team P0-5:

    - double CORDIC must propagate NaN / Inf as NaN
    - fuzz_god_mode must not cast non-finite doubles to integers
    - add regression edge test for CORDIC non-finite behavior

Targets:

    v11S/src/internal/cordic.h
    v11S/tests/fuzz_god_mode.c
    v11S/tests/test_edge_redteam_p0_5.c
    v11S/docs/CLOSURE_PUNCHLIST.md

Usage:

    python3 07_p0_5_cordic_fuzzer_nonfinite.py

Force overwrite:

    python3 07_p0_5_cordic_fuzzer_nonfinite.py --force
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path


CORDIC_MARKER = "MATHLIB_CLOSURE_P2_P0_5_CORDIC_NONFINITE_GUARD"
FUZZ_MARKER = "MATHLIB_CLOSURE_P2_P0_5_FUZZ_CORDIC_NONFINITE_GUARD"
TEST_MARKER = "MATHLIB_CLOSURE_P2_P0_5_TEST"
LOG_MARKER = "<!-- MATHLIB_CLOSURE_P2_P0_5_LOG -->"


NEW_CORDIC_FUNC = r"""static inline void ml_cordic_sincos(double theta, double *sin_out, double *cos_out) {
    /* MATHLIB_CLOSURE_P2_P0_5_CORDIC_NONFINITE_GUARD */
    if (ML_UNLIKELY(sin_out == NULL || cos_out == NULL)) {
        return;
    }

    /*
     * Non-finite inputs must propagate as NaN.
     *
     * The old code allowed NaN to fall through the CORDIC iteration,
     * because comparisons against NaN are false, producing finite-looking
     * garbage instead of NaN.
     */
    if (ml_isnan(theta) || ml_isinf(theta)) {
        *sin_out = ml_make_nan();
        *cos_out = ml_make_nan();
        return;
    }

    /* Range reduction to [-pi, pi] */
    theta = ml_fmod(theta, 2.0 * ML_PI);

    /*
     * Defensive guard:
     * finite input should not produce non-finite reduction, but if it does,
     * fail loudly as NaN instead of casting or iterating on garbage.
     */
    if (ml_isnan(theta) || ml_isinf(theta)) {
        *sin_out = ml_make_nan();
        *cos_out = ml_make_nan();
        return;
    }

    if (theta > ML_PI) theta -= 2.0 * ML_PI;
    if (theta < -ML_PI) theta += 2.0 * ML_PI;

    /* Quadrant mapping: CORDIC only converges in [-pi/2, pi/2] */
    int negate_cos = 0;

    if (theta > ML_PI / 2.0) {
        theta = ML_PI - theta;
        negate_cos = 1;
    } else if (theta < -ML_PI / 2.0) {
        theta = -ML_PI - theta;
        negate_cos = 1;
    }

    double x = CORDIC_GAIN;
    double y = 0.0;
    double z = theta;

    for (int i = 0; i < 24; i++) {
        double x_new, y_new;

        if (z >= 0) {
            x_new = x - (y / (double)(1LL << i));
            y_new = y + (x / (double)(1LL << i));
            z -= cordic_atan[i];
        } else {
            x_new = x + (y / (double)(1LL << i));
            y_new = y - (x / (double)(1LL << i));
            z += cordic_atan[i];
        }

        x = x_new;
        y = y_new;
    }

    if (negate_cos) x = -x;

    *cos_out = x;
    *sin_out = y;
}
"""


NEW_FUZZ_FUNC = r"""void test_fixed_point_cordic() {
    /* MATHLIB_CLOSURE_P2_P0_5_FUZZ_CORDIC_NONFINITE_GUARD */
    printf("--- Fixed-Point CORDIC vs Double CORDIC (1,000 iterations) ---\n");

    for(int i=0; i<1000; i++) {
        double angle = rand_double();

        /*
         * Do NOT cast NaN or Infinity to integer.
         *
         * In C, converting a non-finite floating-point value to an integer
         * type is undefined behavior.
         *
         * Instead, verify that the double CORDIC path correctly returns NaN.
         */
        if (!ml_isfinite(angle)) {
            double d_sin, d_cos;
            ml_cordic_sincos(angle, &d_sin, &d_cos);

            CHECK(ml_isnan(d_sin) && ml_isnan(d_cos),
                  "Double CORDIC non-finite input produces NaN");

            continue;
        }

        /* O(1) range reduction to prevent infinite loop on massive inputs */
        angle = ml_fmod(angle, 2.0 * ML_PI);

        if (!ml_isfinite(angle)) {
            CHECK(0, "ml_fmod(finite, 2*pi) produced non-finite result");
            continue;
        }

        if (angle > ML_PI) angle -= 2.0 * ML_PI;
        if (angle < -ML_PI) angle += 2.0 * ML_PI;

        ml_q16_16_t fixed_angle = (ml_q16_16_t)(angle * 65536.0);

        ml_q16_16_t f_sin, f_cos;
        ml_cordic_sincos_fixed(fixed_angle, &f_sin, &f_cos);

        double d_sin, d_cos;
        ml_cordic_sincos(angle, &d_sin, &d_cos);

        CHECK_NEAR((double)f_sin / 65536.0, d_sin, 1e-3, "Fixed vs Double CORDIC sin");
        CHECK_NEAR((double)f_cos / 65536.0, d_cos, 1e-3, "Fixed vs Double CORDIC cos");
    }
}
"""


TEST_C = r"""/* MATHLIB_CLOSURE_P2_P0_5_TEST */
/* Red Team P2 P0-5: CORDIC non-finite regression tests */
#include "test_harness.h"
#include "ml_core.h"
#include "internal/cordic.h"

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Red Team P0-5");

    double s, c;

    ml_cordic_sincos(ml_make_nan(), &s, &c);
    ASSERT_TRUE(&ctx,
        ml_isnan(s) && ml_isnan(c),
        "CORDIC NaN propagates");

    ml_cordic_sincos(ml_make_inf(0), &s, &c);
    ASSERT_TRUE(&ctx,
        ml_isnan(s) && ml_isnan(c),
        "CORDIC +inf propagates");

    ml_cordic_sincos(-ml_make_inf(0), &s, &c);
    ASSERT_TRUE(&ctx,
        ml_isnan(s) && ml_isnan(c),
        "CORDIC -inf propagates");

    ml_cordic_sincos(0.0, &s, &c);
    ASSERT_NEAR(&ctx, s, 0.0, 1e-6, "CORDIC sin(0)");
    ASSERT_NEAR(&ctx, c, 1.0, 1e-6, "CORDIC cos(0)");

    ml_cordic_sincos(ML_PI / 2.0, &s, &c);
    ASSERT_NEAR(&ctx, s, 1.0, 1e-5, "CORDIC sin(pi/2)");
    ASSERT_NEAR(&ctx, ml_fabs(c), 0.0, 1e-5, "CORDIC cos(pi/2) near zero");

    return ml_test_summary(&ctx);
}
"""


PUNCHLIST_LOG = r"""
<!-- MATHLIB_CLOSURE_P2_P0_5_LOG -->
## Red Team P2 P0-5

- Hardened double CORDIC against NaN / Inf inputs.
- Hardened `fuzz_god_mode` fixed CORDIC test to avoid casting non-finite doubles to integers.
- Added regression suite `tests/test_edge_redteam_p0_5.c`.
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
    print(f"[write] {path}")


def locate_v11s() -> tuple[Path, Path]:
    root = Path.cwd()

    candidate = root / "v11S"
    if candidate.is_dir():
        return root, candidate

    # Convenience: allow running from inside v11S itself.
    if (root / "src" / "internal" / "cordic.h").is_file():
        print("[note] Running from inside v11S; treating current directory as v11S.")
        return root.parent, root

    fail(
        "Run this script from the folder that CONTAINS the v11S directory, "
        "or from inside v11S itself."
    )


def patch_cordic(v11s: Path, force: bool) -> None:
    path = v11s / "src" / "internal" / "cordic.h"

    if not path.is_file():
        fail(f"Missing expected file: {path}")

    text = normalize(path.read_text(encoding="utf-8"))

    if CORDIC_MARKER in text and not force:
        print(f"[skip] {path}: CORDIC non-finite guard already present")
        return

    pattern = re.compile(
        r"(?ms)^static inline void ml_cordic_sincos\("
        r"double theta, double \*sin_out, double \*cos_out\)\s*\{"
        r".*?"
        r"(?=^#endif)"
    )

    patched, count = pattern.subn(
        lambda m: NEW_CORDIC_FUNC + "\n",
        text,
        count=1,
    )

    if count != 1:
        fail(
            f"{path}: expected exactly one ml_cordic_sincos() match, got {count}. "
            "Source may have drifted."
        )

    write_text(path, patched)


def patch_fuzzer(v11s: Path, force: bool) -> None:
    path = v11s / "tests" / "fuzz_god_mode.c"

    if not path.is_file():
        fail(f"Missing expected file: {path}")

    text = normalize(path.read_text(encoding="utf-8"))

    if FUZZ_MARKER in text and not force:
        print(f"[skip] {path}: fuzzer CORDIC non-finite guard already present")
        return

    pattern = re.compile(
        r"(?ms)^void test_fixed_point_cordic\(\)\s*\{"
        r".*?"
        r"(?=^int main\()"
    )

    patched, count = pattern.subn(
        lambda m: NEW_FUZZ_FUNC + "\n\n",
        text,
        count=1,
    )

    if count != 1:
        fail(
            f"{path}: expected exactly one test_fixed_point_cordic() match, got {count}. "
            "Source may have drifted."
        )

    write_text(path, patched)


def write_test(v11s: Path, force: bool) -> None:
    path = v11s / "tests" / "test_edge_redteam_p0_5.c"

    if path.exists() and not force:
        try:
            old = normalize(path.read_text(encoding="utf-8"))
            if TEST_MARKER in old:
                print(f"[skip] {path}: already present")
                return
        except UnicodeDecodeError:
            pass

    write_text(path, TEST_C)


def append_log(v11s: Path) -> None:
    path = v11s / "docs" / "CLOSURE_PUNCHLIST.md"

    if not path.exists():
        write_text(
            path,
            "# MathLib v11S Closure Punchlist\n" + PUNCHLIST_LOG,
        )
        return

    text = normalize(path.read_text(encoding="utf-8"))

    if LOG_MARKER in text:
        print(f"[skip] {path}: P0-5 log already present")
        return

    with open(path, "a", encoding="utf-8", newline="\n") as fh:
        fh.write(PUNCHLIST_LOG)

    print(f"[append] {path}")


def archive_self(v11s: Path, force: bool) -> None:
    try:
        source_script = Path(__file__).resolve()
        archived_script = v11s / "scripts" / "closure" / source_script.name

        if source_script == archived_script:
            return

        if archived_script.exists() and not force:
            print(f"[skip] {archived_script}: already archived")
            return

        archived_script.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_script, archived_script)
        print(f"[archive] {archived_script}")

    except NameError:
        print("[note] Could not archive script because __file__ is unavailable.")


def main() -> int:
    force = "--force" in sys.argv[1:]

    root, v11s = locate_v11s()

    print("=========================================================")
    print("  Red Team P2 P0-5: CORDIC / fuzzer non-finite fix")
    print("=========================================================")
    print(f"Root: {root}")
    print(f"v11S: {v11s}")
    print(f"Force: {force}")
    print("---------------------------------------------------------")

    patch_cordic(v11s, force)
    patch_fuzzer(v11s, force)
    write_test(v11s, force)
    append_log(v11s)
    archive_self(v11s, force)

    print("---------------------------------------------------------")
    print("P0-5 changes applied.")
    print("")
    print("Now rebuild the fuzzer.")
    print("")
    print("If using CMake:")
    print("")
    print("    cd v11S")
    print("    cmake --build build")
    print("    ./build/fuzz_god_mode 123456789")
    print("")
    print("If using Make:")
    print("")
    print("    cd v11S")
    print("    make fuzz_god_mode")
    print("    ./build/fuzz_god_mode 123456789")
    print("")
    print("Also run edge tests:")
    print("")
    print("    MATHLIB_EDGE_SANITIZERS=1 bash tests/run_edge_tests.sh")
    print("=========================================================")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
