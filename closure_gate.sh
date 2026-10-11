#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

echo "========================================================="
echo "  MATHLIB v12R2: DEVELOPMENT GATE (R5)"
echo "========================================================="

CC="${CC:-gcc}"
SEED="${MATHLIB_ULTIMATE_SEED:-123456789}"
# R5: honor MATHLIB_BUILD_DIR (default build) so /tmp builds never touch
# the repo tree; never rm -rf an existing tree destructively without opt-in.
BUILD_DIR="${MATHLIB_BUILD_DIR:-build}"
if [ "${MATHLIB_GATE_CLEAN:-0}" = "1" ]; then
    rm -rf "$BUILD_DIR"
fi

echo "[1/10] Configuring with ASan + UBSan ($BUILD_DIR)..."
cmake -B "$BUILD_DIR" \
    -DMATHLIB_PROFILE=SCIENTIFIC \
    -DCMAKE_BUILD_TYPE=Debug \
    -DMATHLIB_SANITIZERS=ON

echo "[2/10] Building..."
cmake --build "$BUILD_DIR"

echo "[3/10] Running modular tests..."
"$BUILD_DIR/test_core"
"$BUILD_DIR/test_trig"
"$BUILD_DIR/test_linalg"
"$BUILD_DIR/test_dsp"

echo "[4/10] Running monolithic smoke test..."
"$BUILD_DIR/test"

echo "[5/10] Running edge tests with sanitizers..."
MATHLIB_EDGE_SANITIZERS=1 MATHLIB_EDGE_BUILD_DIR="$BUILD_DIR" bash tests/run_edge_tests.sh

echo "[6/10] Running boundary gauntlet..."
"$BUILD_DIR/fuzz_boundary"

echo "[7/10] Running god-mode fuzzer (seed $SEED)..."
"$BUILD_DIR/fuzz_god_mode" "$SEED"

echo "[8/10] Running mpmath oracle validation..."
# MATHLIB_V12A1_ORACLE_BUILD_FIX: oracle_check is produced by the CMake
# build above and inherits sanitizer link flags. No manual gcc line.
"$BUILD_DIR/oracle_check"
echo "[9/10] Running ultimate fuzzer..."
# R5: prefer the CMake ultimate_fuzzer target when present (single source
# of truth); fall back to the legacy ad-hoc gcc line otherwise.
if [ -x "$BUILD_DIR/ultimate_fuzzer" ]; then
    "$BUILD_DIR/ultimate_fuzzer" "$SEED"
else
    "$CC" -std=c99 -O3 \
        -fsanitize=address,undefined -fno-omit-frame-pointer \
        -Iinclude/mathlib -Isrc \
        -o "$BUILD_DIR/ultimate_fuzzer" \
        tests/ultimate_fuzzer.c \
        -L"$BUILD_DIR" -lmathc -lm
    "$BUILD_DIR/ultimate_fuzzer" "$SEED"
fi

echo "[10/10] Running modular python runner..."
(cd "$BUILD_DIR" && python3 ../tests/run_tests.py)

echo "========================================================="
echo "  v12R2 DEVELOPMENT GATE PASSED"
echo "========================================================="
