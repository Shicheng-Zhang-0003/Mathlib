#!/usr/bin/env python3
"""
02_p1_simd_index_ci_docs.py

Run from the folder that CONTAINS the v11S working folder.

This script applies the P1 closure hardening layer after P0 has passed.

Targets:

    v11S/include/mathlib/simd_batch.h
    v11S/tests/test_edge_audit_ip1.c
    v11S/src/cpu_dispatch.c
    v11S/include/mathlib/ml_tensor.h
    verify_v11s.sh
    v11S/.github/workflows/ci.yml
    v11S/docs/CLOSURE_P1_APPLIED.md
    v11S/docs/CLOSURE_PUNCHLIST.md

Usage:

    python3 02_p1_simd_index_ci_docs.py
    python3 02_p1_simd_index_ci_docs.py --force
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

/* MATHLIB_CLOSURE_P1_SIMD_BATCH_GUARD */

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
     */
    for (int i = 0; i < 4; i++) {
        if (ml_isnan(in[i]) || ml_isinf(in[i]) ||
            ml_is_subnormal(in[i]) || !(in[i] > 0.0)) {
            for (int j = 0; j < 4; j++) {
                if (ml_isnan(in[j])) {
                    out[j] = ml_make_nan();
                } else if (in[j] > 0.0) {
                    out[j] = 1.0 / ml_sqrt(in[j]);
                } else {
                    out[j] = 0.0;
                }
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
    for(int i=0; i<4; i++) out[i] = (in[i] * in[i]) + in[i];
}

static inline void ml_simd_batch_rsqrt(const double* in, double* out) {
    /* v11S AUDIT IP-1: NaN propagation in scalar fallback */
    for(int i=0; i<4; i++) {
        if (ml_isnan(in[i])) {
            out[i] = ml_make_nan();
        } else if (in[i] > 0.0) {
            out[i] = 1.0 / ml_sqrt(in[i]);
        } else {
            out[i] = 0.0;
        }
    }
}

#endif

#endif
"""


NEW_TEST_AUDIT_IP1_C = r"""/* v11S AUDIT IP-1: regression tests */
#include "test_harness.h"
#include "ml_core.h"
#include "ml_trig.h"
#include "ml_complex.h"
#include "ieee754.h"
#include "fast_math.h"
#include "simd_batch.h"
#include <stdint.h>
#include <string.h>

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Audit IP-1");

    /* atan2 NaN / signed-zero / infinity behavior */
    ASSERT_TRUE(&ctx, ml_isnan(ml_atan2(0.0, ml_make_nan())), "atan2(0,NaN) is NaN");
    ASSERT_TRUE(&ctx, ml_isnan(ml_atan2(ml_make_nan(), 0.0)), "atan2(NaN,0) is NaN");
    ASSERT_TRUE(&ctx, ml_signbit(ml_atan2(-0.0, ml_make_inf(0))) != 0, "atan2(-0,+inf) preserves sign");
    ASSERT_TRUE(&ctx, ml_signbit(ml_atan2(-1.0, ml_make_inf(0))) != 0, "atan2(-finite,+inf) preserves sign");
    ASSERT_NEAR(&ctx, ml_atan2(0.0, -1.0), ML_PI, 1e-15, "atan2(+0,-1)");
    ASSERT_NEAR(&ctx, ml_atan2(-0.0, -1.0), -ML_PI, 1e-15, "atan2(-0,-1)");
    ASSERT_TRUE(&ctx, ml_signbit(ml_atan2(-0.0, 1.0)) != 0, "atan2(-0,+1) preserves sign");

    /* complex power special cases */
    cplx zero = {0.0, 0.0};
    cplx one = {1.0, 0.0};
    cplx nanv = {ml_make_nan(), 0.0};

    cplx p00 = ml_cplx_power(zero, zero);
    ASSERT_TRUE(&ctx, p00.real == 1.0 && p00.imag == 0.0, "cplx pow(0,0) == 1");

    cplx p1n = ml_cplx_power(one, nanv);
    ASSERT_TRUE(&ctx, p1n.real == 1.0 && p1n.imag == 0.0, "cplx pow(1,NaN) == 1");

    /* experimental ieee754 helpers are now safe */
    ASSERT_TRUE(&ctx, ml_isnan(logarithm_ieee754(-1.0)), "ieee log(-1) is NaN");
    ASSERT_TRUE(&ctx, ml_isinf(logarithm_ieee754(0.0)) && logarithm_ieee754(0.0) < 0.0, "ieee log(0) is -inf");
    ASSERT_NEAR(&ctx, logarithm_ieee754(1.0), 0.0, 1e-15, "ieee log(1)");
    ASSERT_NEAR(&ctx, logarithm_ieee754(2.0), ML_LN2, 1e-12, "ieee log(2)");
    ASSERT_NEAR(&ctx, logarithm_ieee754(0.5), -ML_LN2, 1e-12, "ieee log(0.5)");

    ASSERT_NEAR(&ctx, exponential_ieee754(1.0), ML_E, 1e-10, "ieee exp(1)");
    ASSERT_TRUE(&ctx, ml_isinf(exponential_ieee754(1000.0)) && exponential_ieee754(1000.0) > 0.0, "ieee exp(1000) is +inf");
    ASSERT_TRUE(&ctx, exponential_ieee754(-1000.0) == 0.0, "ieee exp(-1000) is 0");

    /* fast math subnormal fallback */
    uint64_t min_bits = 1;
    double min_sub;
    memcpy(&min_sub, &min_bits, sizeof(double));
    ASSERT_TRUE(&ctx, ml_isfinite(ml_fast_rsqrt(min_sub)), "fast rsqrt subnormal finite");

    /* SIMD batch NaN propagation is now unconditional after P1 guard */
    double in[4];
    double out[4];

    in[0] = 1.0;
    in[1] = 4.0;
    in[2] = 0.0;
    in[3] = ml_make_nan();

    ml_simd_batch_rsqrt(in, out);

    /* MATHLIB_CLOSURE_P1_SIMD_BATCH_GUARD_TEST */
    ASSERT_TRUE(&ctx, ml_isnan(out[3]), "batch rsqrt propagates NaN");

    /* fmod exactness basics */
    ASSERT_NEAR(&ctx, ml_fmod(5.0, 3.0), 2.0, 1e-15, "fmod(5,3)");
    ASSERT_NEAR(&ctx, ml_fmod(-5.0, 3.0), -2.0, 1e-15, "fmod(-5,3)");
    ASSERT_TRUE(&ctx, ml_isnan(ml_fmod(1.0, 0.0)), "fmod(x,0) is NaN");
    ASSERT_TRUE(&ctx, ml_signbit(ml_fmod(-0.0, 1.0)) != 0, "fmod(-0,1) preserves sign");

    return ml_test_summary(&ctx);
}
"""


NEW_CPU_DISPATCH_C = r"""#include "cpu_dispatch.h"
#include <stddef.h>

#if defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#endif

/* MATHLIB_CLOSURE_P1_SIZE_T_INDEXING */

static void ml_matmul_scalar(const double* ML_RESTRICT A, const double* ML_RESTRICT B, double* ML_RESTRICT C, int N) {
    size_t n = (size_t)N;

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            double sum = 0.0;
            for (size_t k = 0; k < n; k++) {
                sum += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = sum;
        }
    }
}

#if ML_COMPILE_TIME_AVX2
ML_TARGET_AVX2
static void ml_matmul_avx2(const double* ML_RESTRICT A, const double* ML_RESTRICT B, double* ML_RESTRICT C, int N) {
    size_t n = (size_t)N;

    for (size_t i = 0; i < n; i++) {
        size_t j = 0;

        for (; j + 4 <= n; j += 4) {
            __m256d c_vec = _mm256_setzero_pd();

            for (size_t k = 0; k < n; k++) {
                __m256d a_vec = _mm256_broadcast_sd(&A[i * n + k]);
                __m256d b_vec = _mm256_loadu_pd(&B[k * n + j]);
                c_vec = _mm256_fmadd_pd(a_vec, b_vec, c_vec);
            }

            _mm256_storeu_pd(&C[i * n + j], c_vec);
        }

        for (; j < n; j++) {
            double sum = 0.0;
            for (size_t k = 0; k < n; k++) {
                sum += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = sum;
        }
    }
}
#endif

ML_API void ml_matmul(const double* ML_RESTRICT A, const double* ML_RESTRICT B, double* ML_RESTRICT C, int N) {
    if (ML_UNLIKELY(A == NULL || B == NULL || C == NULL || N <= 0)) return;

    size_t n = (size_t)N;

    /* Defensive guard against index-space overflow. */
    if (ML_UNLIKELY(n > ((size_t)-1) / n)) return;

#if ML_COMPILE_TIME_AVX2
    ml_matmul_avx2(A, B, C, N);
#else
    ml_matmul_scalar(A, B, C, N);
#endif
}

ML_API int ml_cpu_has_avx2(void) {
    return ML_COMPILE_TIME_AVX2;
}

ML_API int ml_cpu_has_fma(void) {
#if defined(__FMA__) && (defined(__x86_64__) || defined(__i386__))
    return 1;
#else
    return 0;
#endif
}

ML_API int ml_cpu_has_sse41(void) {
#if defined(__SSE4_1__) && (defined(__x86_64__) || defined(__i386__))
    return 1;
#else
    return 0;
#endif
}

ML_API int ml_cpu_has_neon(void) {
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    return 1;
#else
    return 0;
#endif
}
"""


NEW_TENSOR_H = r"""#ifndef MATHLIB_ML_TENSOR_H
#define MATHLIB_ML_TENSOR_H

#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include "ml_compiler.h"
#include "ml_types.h"

/* ============================================================================
 * v11S TENSOR & WORKSPACE SYSTEM
 *
 * SAFETY BY DEFAULT: All bounds checks, NULL checks, and dimension validation
 * are UNCONDITIONAL. The MATHLIB_PROFILE_HARDENED flag now ONLY controls
 * memory canaries and NaN poisoning, NOT basic safety.
 *
 * ABI STABILITY: The ml_workspace_t struct has a FIXED size regardless of
 * compile flags. The magic_canary is ALWAYS present to guarantee binary
 * compatibility across translation units.
 * ========================================================================== */

typedef struct {
    uint8_t *base;
    size_t capacity;
    size_t offset;
    uint64_t magic_canary;
} ml_workspace_t;

#define ML_WORKSPACE_CANARY 0xDEADBEEFCAFEBABEULL

ML_INLINE void ml_workspace_init(ml_workspace_t *ws, void *buffer, size_t size) {
    if (ML_UNLIKELY(!ws)) return;

    if (ML_UNLIKELY(!buffer || size == 0)) {
        ws->base = NULL;
        ws->capacity = 0;
        ws->offset = 0;
        ws->magic_canary = 0;
        return;
    }

    uintptr_t raw = (uintptr_t)buffer;
    uintptr_t aligned = (raw + (uintptr_t)31) & ~(uintptr_t)31;
    size_t pad = (size_t)(aligned - raw);

    if (ML_UNLIKELY(pad >= size)) {
        ws->base = NULL;
        ws->capacity = 0;
        ws->offset = 0;
        ws->magic_canary = 0;
        return;
    }

    ws->base = (uint8_t *)aligned;
    ws->capacity = size - pad;
    ws->offset = 0;
    ws->magic_canary = ML_WORKSPACE_CANARY;
}

ML_INLINE void *ml_workspace_alloc(ml_workspace_t *ws, size_t bytes) {
    if (ML_UNLIKELY(!ws)) return NULL;
    if (ML_UNLIKELY(ws->magic_canary != ML_WORKSPACE_CANARY)) return NULL;

    size_t aligned = (bytes + 31) & ~(size_t)31;

    if (ML_UNLIKELY(aligned < bytes)) return NULL;
    if (ML_UNLIKELY(ws->offset > ws->capacity)) return NULL;
    if (ML_UNLIKELY(aligned > ws->capacity - ws->offset)) return NULL;

    void *ptr = ws->base + ws->offset;
    ws->offset += aligned;
    return ptr;
}

ML_INLINE void ml_workspace_reset(ml_workspace_t *ws) {
    if (ML_UNLIKELY(!ws)) return;
    if (ML_UNLIKELY(ws->magic_canary != ML_WORKSPACE_CANARY)) return;
    ws->offset = 0;
}

typedef struct {
    double *data;
    int rows;
    int cols;
} ml_tensor_view_t;

ML_INLINE ml_tensor_view_t ml_tensor_view(double *data, int rows, int cols) {
    return (ml_tensor_view_t){data, rows, cols};
}

/* MATHLIB_CLOSURE_P1_SIZE_T_INDEXING */
ML_INLINE double *ml_tensor_at(ml_tensor_view_t t, int r, int c) {
    if (ML_UNLIKELY(!t.data)) return NULL;
    if (ML_UNLIKELY(r < 0 || r >= t.rows || c < 0 || c >= t.cols)) return NULL;

    return &t.data[(size_t)r * (size_t)t.cols + (size_t)c];
}

#define ML_TENSOR_AT(t, r, c) ((t).data[(size_t)(r) * (size_t)(t).cols + (size_t)(c)])

#endif /* MATHLIB_ML_TENSOR_H */
"""


NEW_VERIFY_SH = r"""#!/bin/bash
set -e
trap 'rm -f /tmp/oracle_check /tmp/ultimate_fuzz' EXIT

# MATHLIB_CLOSURE_P1_VERIFY_GATE

echo "========================================================="
echo "  MATHLIB v11S: THE ULTIMATE VERIFICATION GAUNTLET"
echo "========================================================="

cd v11S

echo "[1/7] Configuring with ASan + UBSan..."
rm -rf build && mkdir -p build && cd build

if ! CMAKE_OUT=$(cmake .. -DMATHLIB_PROFILE=SCIENTIFIC -DCMAKE_BUILD_TYPE=Debug -DMATHLIB_SANITIZERS=ON 2>&1); then
    echo "❌ FAIL: CMake configuration failed."
    echo "$CMAKE_OUT"
    exit 1
fi

echo "[2/7] Building with maximum paranoia..."
if ! BUILD_OUT=$(cmake --build . 2>&1); then
    echo "❌ FAIL: Build failed."
    echo "$BUILD_OUT"
    exit 1
fi

cd ..

echo "[3/7] Running Modular CI Tests..."
./build/test_core > /dev/null
./build/test_trig > /dev/null
./build/test_linalg > /dev/null
./build/test_dsp > /dev/null
echo "✅ PASS: Modular CI Tests passed."

echo "[4/7] Running Edge Tests (ASan + UBSan)..."
if ! MATHLIB_EDGE_SANITIZERS=1 bash tests/run_edge_tests.sh > /dev/null; then
    echo "❌ FAIL: Edge tests failed."
    exit 1
fi
echo "✅ PASS: Edge tests passed."

echo "[5/7] Running Boundary Gauntlet..."
./build/fuzz_boundary > /dev/null
echo "✅ PASS: Boundary Gauntlet passed."

echo "[6/7] Running mpmath Oracle..."
gcc -std=c99 -O3 -fPIE -fsanitize=address,undefined -fno-omit-frame-pointer -Iinclude/mathlib -Isrc -DMATHLIB_HAS_ORACLE_DATA -o /tmp/oracle_check tests/test_oracle.c -Lbuild -lmathc -lm
/tmp/oracle_check

echo "[7/7] Unleashing the Ultimate Fuzzer (ASan/UBSan)..."
gcc -std=c99 -O3 -fsanitize=address,undefined -fno-omit-frame-pointer -Iinclude/mathlib -Isrc -o /tmp/ultimate_fuzz tests/ultimate_fuzzer.c -Lbuild -lmathc -lm
/tmp/ultimate_fuzz "${MATHLIB_ULTIMATE_SEED:-123456789}"

echo "========================================================="
echo "🎉 ALL VERIFICATION CHECKS PASSED. v11S closure candidate verified."
echo "========================================================="
"""


NEW_CI_YML = """name: MathLib A3 CI/CD Matrix

on:
  push:
    branches: [ "main", "v11A3" ]
  pull_request:
    branches: [ "main" ]

jobs:
  build-and-test:
    runs-on: ubuntu-latest

    strategy:
      matrix:
        compiler: [gcc, clang]
        profile: [SCIENTIFIC, GRAPHICS, EMBEDDED]
        sanitizer: [none, asan]

    steps:
      - uses: actions/checkout@v3

      - name: Configure CMake
        run: |
          SANITIZER_FLAG=""
          if [ "${{ matrix.sanitizer }}" == "asan" ]; then
            SANITIZER_FLAG="-DMATHLIB_SANITIZERS=ON"
          fi
          cmake -B build -DCMAKE_C_COMPILER=${{ matrix.compiler }} -DMATHLIB_PROFILE=${{ matrix.profile }} $SANITIZER_FLAG

      - name: Build
        run: cmake --build build

      - name: Run Smoke Tests
        run: ./build/test

      - name: Run Modular CI/CD Tests
        run: |
          ./build/test_core
          ./build/test_trig
          ./build/test_linalg
          ./build/test_dsp

      # MATHLIB_CLOSURE_P1_CI_DETERMINISTIC
      - name: Run Edge Tests
        run: bash tests/run_edge_tests.sh

      - name: Run Boundary Gauntlet
        run: ./build/fuzz_boundary

      - name: Quick Fuzz (deterministic)
        run: ./build/fuzz_god_mode 123456789
"""


P1_APPLIED_MD = r"""# MathLib v11S P1 Safety / Process Fixes Applied

This file records that the P1 closure hardening script has been applied.

Script:

    02_p1_simd_index_ci_docs.py

Applied fixes:

- AVX2 `ml_simd_batch_rsqrt` now falls back to scalar semantics for exceptional inputs
- AVX2 batch NaN propagation test is now unconditional
- matrix multiplication indexing now uses wide size arithmetic internally
- tensor offset arithmetic now uses wide size arithmetic
- official verification script now runs edge tests and boundary gauntlet
- CI now runs edge tests and deterministic quick fuzz

Next validation step:

    cd v11S
    ./closure_gate.sh
"""


PUNCHLIST_LOG = r"""

<!-- MATHLIB_CLOSURE_P1_LOG -->
## Closure Log

- P0 source fixes applied and strict closure gate passed.
- P1 safety/process fixes applied by `02_p1_simd_index_ci_docs.py`.
- AVX2 batch rsqrt now falls back to scalar semantics for exceptional inputs.
- Matrix/tensor indexing now uses wide size arithmetic internally.
- Official verification script now runs edge tests and boundary gauntlet.
- CI now runs edge tests and deterministic quick fuzz.
"""


def fail(message: str) -> None:
    print("ERROR: " + message)
    sys.exit(1)


def write_full(
    path: Path,
    content: str,
    marker: str,
    force: bool,
    executable: bool = False,
) -> bool:
    if path.exists() and not force:
        try:
            old = path.read_text(encoding="utf-8")
            if marker in old:
                print(f"[skip] {path}: marker already present")
                return False
        except UnicodeDecodeError:
            pass

    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(content)

    if executable:
        path.chmod(0o755)

    print(f"[write] {path}")
    return True


def append_punchlist_log(path: Path) -> None:
    marker = "<!-- MATHLIB_CLOSURE_P1_LOG -->"

    if path.exists():
        text = path.read_text(encoding="utf-8")
        if marker in text:
            print(f"[skip] {path}: closure log already present")
            return

        with open(path, "a", encoding="utf-8", newline="\n") as fh:
            fh.write(PUNCHLIST_LOG)

        print(f"[append] {path}")
    else:
        write_full(path, "# MathLib v11S Closure Punchlist\n" + PUNCHLIST_LOG, marker, force=True)


def main() -> int:
    force = "--force" in sys.argv[1:]

    root = Path.cwd()
    v11s = root / "v11S"

    if not v11s.is_dir():
        fail("Run this script from the folder that CONTAINS the v11S directory.")

    print("=========================================================")
    print("  MathLib v11S P1 Safety / Process Fixes")
    print("=========================================================")
    print(f"Root: {root}")
    print(f"v11S: {v11s}")
    print(f"Force: {force}")
    print("---------------------------------------------------------")

    write_full(
        v11s / "include" / "mathlib" / "simd_batch.h",
        NEW_SIMD_BATCH_H,
        "MATHLIB_CLOSURE_P1_SIMD_BATCH_GUARD",
        force,
    )

    write_full(
        v11s / "tests" / "test_edge_audit_ip1.c",
        NEW_TEST_AUDIT_IP1_C,
        "MATHLIB_CLOSURE_P1_SIMD_BATCH_GUARD_TEST",
        force,
    )

    write_full(
        v11s / "src" / "cpu_dispatch.c",
        NEW_CPU_DISPATCH_C,
        "MATHLIB_CLOSURE_P1_SIZE_T_INDEXING",
        force,
    )

    write_full(
        v11s / "include" / "mathlib" / "ml_tensor.h",
        NEW_TENSOR_H,
        "MATHLIB_CLOSURE_P1_SIZE_T_INDEXING",
        force,
    )

    write_full(
        root / "verify_v11s.sh",
        NEW_VERIFY_SH,
        "MATHLIB_CLOSURE_P1_VERIFY_GATE",
        force,
        executable=True,
    )

    write_full(
        v11s / ".github" / "workflows" / "ci.yml",
        NEW_CI_YML,
        "MATHLIB_CLOSURE_P1_CI_DETERMINISTIC",
        force,
    )

    write_full(
        v11s / "docs" / "CLOSURE_P1_APPLIED.md",
        P1_APPLIED_MD,
        "P1 Safety / Process Fixes Applied",
        force,
    )

    append_punchlist_log(v11s / "docs" / "CLOSURE_PUNCHLIST.md")

    try:
        source_script = Path(__file__).resolve()
        archived_script = v11s / "scripts" / "closure" / "02_p1_simd_index_ci_docs.py"

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
    print("P1 safety/process fixes applied.")
    print("")
    print("Next step:")
    print("    cd v11S && ./closure_gate.sh")
    print("")
    print("After that, the tree is ready for Red Team review.")
    print("=========================================================")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
