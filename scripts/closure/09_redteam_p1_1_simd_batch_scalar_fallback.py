#!/usr/bin/env python3
"""
09_redteam_p1_1_simd_batch_scalar_fallback.py

Run from the folder that CONTAINS the v11S working folder.

Expected layout:

    parent_folder/
        09_redteam_p1_1_simd_batch_scalar_fallback.py
        v11S/
            include/
            tests/
            docs/
            ...

This script implements Red Team P2 item P1-1:

    Align SIMD batch rsqrt fallback semantics with scalar fast-math.

It modifies:

    v11S/include/mathlib/simd_batch.h

It adds:

    v11S/tests/test_edge_redteam_p1_1.c

It appends a closure log entry to:

    v11S/docs/CLOSURE_PUNCHLIST.md

Usage:

    python3 09_redteam_p1_1_simd_batch_scalar_fallback.py

Force overwrite:

    python3 09_redteam_p1_1_simd_batch_scalar_fallback.py --force
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path


NEW_SIMD_BATCH_H = r"""#ifndef MATHLIB_SIMD_BATCH_H
#define MATHLIB_SIMD_BATCH_H

#include "profiles.h"
#include "ml_core.h"
#include "bitwise_fp.h"
#include "fast_math.h"

/* MATHLIB_CLOSURE_P1_SIMD_BATCH_GUARD */
/* MATHLIB_REDP2_P1_1_SIMD_BATCH_SCALAR_FALLBACK */

#if defined(__AVX2__)
#include <immintrin.h>

typedef double ml_vec4d __attribute__((aligned(32), vector_size(32)));

static inline void ml_simd_batch_poly(const double* in, double* out) {
    __m256d v = _mm256_loadu_pd(in);
    __m256d res = _mm256_add_pd(_mm256_mul_pd(v, v), v);
    _mm256_storeu_pd(out, res);
}

static inline void ml_simd_batch_rsqrt(const double* in, double* out) {
    /*
     * MATHLIB_CLOSURE_P1:
     * Exceptional inputs fall back to scalar semantics for the whole batch.
     * The AVX2 fast path is used only for positive finite normal inputs.
     *
     * MATHLIB_REDP2_P1_1:
     * The scalar fallback now uses ml_fast_rsqrt() per lane, so batch
     * exceptional behavior matches the scalar fast-math contract exactly:
     *
     *   NaN             -> NaN
     *   negative finite -> NaN
     *   +inf            -> +0
     *   +0              -> +inf
     *   -0              -> -inf
     *   positive finite -> approximate rsqrt
     */
    for (int i = 0; i < 4; i++) {
        if (ml_isnan(in[i]) || ml_isinf(in[i]) ||
            ml_is_subnormal(in[i]) || !(in[i] > 0.0)) {
            for (int j = 0; j < 4; j++) {
                out[j] = ml_fast_rsqrt(in[j]);
            }
            return;
        }
    }

    __m256d v = _mm256_loadu_pd(in);

    __m256i vi = _mm256_castpd_si256(v);
    __m256i magic = _mm256_set1_epi64x(0x5fe6ec85e7de30daLL);
    __m256i vi_half = _mm256_srli_epi64(vi, 1);
    __m256i yi = _mm256_sub_epi64(magic, vi_half);

    __m256d y = _mm256_castsi256_pd(yi);

    __m256d half = _mm256_set1_pd(0.5);
    __m256d three_half = _mm256_set1_pd(1.5);

    __m256d x_half = _mm256_mul_pd(v, half);
    __m256d y2 = _mm256_mul_pd(y, y);
    __m256d sub = _mm256_sub_pd(three_half, _mm256_mul_pd(x_half, y2));

    y = _mm256_mul_pd(y, sub);

    y2 = _mm256_mul_pd(y, y);
    sub = _mm256_sub_pd(three_half, _mm256_mul_pd(x_half, y2));

    y = _mm256_mul_pd(y, sub);

    _mm256_storeu_pd(out, y);
}

#else

/* Scalar Fallback */

static inline void ml_simd_batch_poly(const double* in, double* out) {
    for (int i = 0; i < 4; i++) {
        out[i] = (in[i] * in[i]) + in[i];
    }
}

static inline void ml_simd_batch_rsqrt(const double* in, double* out) {
    /*
     * MATHLIB_REDP2_P1_1:
     * Use the scalar fast-math contract directly.
     */
    for (int i = 0; i < 4; i++) {
        out[i] = ml_fast_rsqrt(in[i]);
    }
}

#endif

#endif /* MATHLIB_SIMD_BATCH_H */
"""


TEST_C = r"""/* MATHLIB_REDP2_P1_1_TEST */
/* Red Team P2 P1-1: SIMD batch rsqrt scalar-fallback semantics */

#include "test_harness.h"
#include "simd_batch.h"
#include "fast_math.h"
#include "ml_core.h"

#include <stdint.h>
#include <string.h>

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Red Team P1-1");

    double in[4];
    double out[4];

    /*
     * Exceptional batch:
     *   +0, -0, negative finite, NaN
     *
     * This must force scalar fallback for the whole batch.
     */
    in[0] = 0.0;
    in[1] = -0.0;
    in[2] = -1.0;
    in[3] = ml_make_nan();

    ml_simd_batch_rsqrt(in, out);

    ASSERT_TRUE(&ctx,
        ml_isinf(out[0]) && out[0] > 0.0,
        "batch rsqrt(+0) gives +inf");

    ASSERT_TRUE(&ctx,
        ml_isinf(out[1]) && ml_signbit(out[1]) != 0,
        "batch rsqrt(-0) gives -inf");

    ASSERT_TRUE(&ctx,
        ml_isnan(out[2]),
        "batch rsqrt(negative finite) gives NaN");

    ASSERT_TRUE(&ctx,
        ml_isnan(out[3]),
        "batch rsqrt(NaN) gives NaN");

    /*
     * Exceptional batch:
     *   +inf, -inf, smallest subnormal, positive finite
     */
    uint64_t min_bits = 1;
    double min_sub;
    memcpy(&min_sub, &min_bits, sizeof(double));

    in[0] = ml_make_inf(0);
    in[1] = -ml_make_inf(0);
    in[2] = min_sub;
    in[3] = 4.0;

    ml_simd_batch_rsqrt(in, out);

    ASSERT_TRUE(&ctx,
        out[0] == 0.0 && ml_signbit(out[0]) == 0,
        "batch rsqrt(+inf) gives +0");

    ASSERT_TRUE(&ctx,
        ml_isnan(out[1]),
        "batch rsqrt(-inf) gives NaN");

    double scalar_sub = ml_fast_rsqrt(min_sub);
    ASSERT_TRUE(&ctx,
        ml_isfinite(out[2]) && out[2] == scalar_sub,
        "batch rsqrt(subnormal) matches scalar ml_fast_rsqrt");

    double scalar_four = ml_fast_rsqrt(4.0);
    ASSERT_TRUE(&ctx,
        out[3] == scalar_four,
        "batch rsqrt(4) matches scalar ml_fast_rsqrt in fallback");

    /*
     * All-positive normal batch:
     * This may use the AVX2 fast path, so compare approximately.
     */
    double in_pos[4] = {1.0, 4.0, 16.0, 100.0};
    double out_pos[4];

    ml_simd_batch_rsqrt(in_pos, out_pos);

    for (int i = 0; i < 4; i++) {
        double expected = 1.0 / ml_sqrt(in_pos[i]);
        ASSERT_NEAR(&ctx,
            out_pos[i],
            expected,
            1e-12,
            "batch rsqrt positive normal approximate");
    }

    return ml_test_summary(&ctx);
}
"""


PUNCHLIST_LOG = r"""
<!-- MATHLIB_REDP2_P1_1_LOG -->
## Red Team P2 P1-1

- Aligned `ml_simd_batch_rsqrt()` exceptional fallback with scalar `ml_fast_rsqrt()`.
- Batch fallback now preserves scalar fast-math semantics for zero, negative, inf, NaN, and subnormal inputs.
- Added regression suite `tests/test_edge_redteam_p1_1.c`.
"""


def fail(message: str) -> None:
    print("ERROR: " + message)
    sys.exit(1)


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
    if (root / "include" / "mathlib" / "simd_batch.h").is_file():
        print("[note] Running from inside v11S; treating current directory as v11S.")
        return root.parent, root

    fail(
        "Run this script from the folder that CONTAINS the v11S directory, "
        "or from inside v11S itself."
    )


def write_header(v11s: Path, force: bool) -> None:
    path = v11s / "include" / "mathlib" / "simd_batch.h"

    if not path.is_file():
        fail(f"Missing expected file: {path}")

    marker = "MATHLIB_REDP2_P1_1_SIMD_BATCH_SCALAR_FALLBACK"

    try:
        old = path.read_text(encoding="utf-8")
        if marker in old and not force:
            print(f"[skip] {path}: P1-1 marker already present")
            return
    except UnicodeDecodeError:
        pass

    write_text(path, NEW_SIMD_BATCH_H)


def write_test(v11s: Path, force: bool) -> None:
    path = v11s / "tests" / "test_edge_redteam_p1_1.c"

    if path.exists() and not force:
        try:
            old = path.read_text(encoding="utf-8")
            if "MATHLIB_REDP2_P1_1_TEST" in old:
                print(f"[skip] {path}: already present")
                return
        except UnicodeDecodeError:
            pass

    write_text(path, TEST_C)


def append_log(v11s: Path) -> None:
    path = v11s / "docs" / "CLOSURE_PUNCHLIST.md"
    marker = "<!-- MATHLIB_REDP2_P1_1_LOG -->"

    if not path.exists():
        write_text(
            path,
            "# MathLib v11S Closure Punchlist\n" + PUNCHLIST_LOG,
        )
        return

    try:
        old = path.read_text(encoding="utf-8")
        if marker in old:
            print(f"[skip] {path}: P1-1 log already present")
            return
    except UnicodeDecodeError:
        pass

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
    print("  Red Team P2 P1-1: SIMD batch rsqrt scalar fallback")
    print("=========================================================")
    print(f"Root: {root}")
    print(f"v11S: {v11s}")
    print(f"Force: {force}")
    print("---------------------------------------------------------")

    write_header(v11s, force)
    write_test(v11s, force)
    append_log(v11s)
    archive_self(v11s, force)

    print("---------------------------------------------------------")
    print("P1-1 changes applied.")
    print("")
    print("Next verification step:")
    print("")
    print("    cd v11S")
    print("    MATHLIB_EDGE_SANITIZERS=1 bash tests/run_edge_tests.sh")
    print("")
    print("Expected new suite:")
    print("")
    print("    --- test_edge_redteam_p1_1 ---")
    print("    [Edge Red Team P1-1] Passed: ..., Failed: 0")
    print("")
    print("After that passes, the next block is:")
    print("")
    print("    P1-2: complex exp must return NaN when trig components are NaN")
    print("=========================================================")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
