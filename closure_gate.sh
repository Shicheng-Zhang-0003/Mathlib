#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

echo "========================================================="
echo "  MATHLIB v11S: STRICT CLOSURE GATE"
echo "========================================================="

CC="${CC:-gcc}"
SEED="${MATHLIB_ULTIMATE_SEED:-123456789}"

echo "[1/7] Configuring with ASan + UBSan..."
rm -rf build
cmake -B build \
    -DMATHLIB_PROFILE=SCIENTIFIC \
    -DCMAKE_BUILD_TYPE=Debug \
    -DMATHLIB_SANITIZERS=ON

echo "[2/7] Building..."
cmake --build build

echo "[3/7] Running modular tests..."
./build/test_core
./build/test_trig
./build/test_linalg
./build/test_dsp

echo "[4/7] Running edge tests with sanitizers..."
MATHLIB_EDGE_SANITIZERS=1 bash tests/run_edge_tests.sh

echo "[5/7] Running boundary gauntlet..."
./build/fuzz_boundary

echo "[6/7] Running mpmath oracle validation..."
"$CC" -std=c99 -O3 -fPIE \
    -fsanitize=address,undefined -fno-omit-frame-pointer \
    -Iinclude/mathlib -Isrc \
    -DMATHLIB_HAS_ORACLE_DATA \
    -o build/oracle_check \
    tests/test_oracle.c \
    -Lbuild -lmathc -lm

./build/oracle_check

echo "[7/7] Running ultimate fuzzer..."
"$CC" -std=c99 -O3 \
    -fsanitize=address,undefined -fno-omit-frame-pointer \
    -Iinclude/mathlib -Isrc \
    -o build/ultimate_fuzzer \
    tests/ultimate_fuzzer.c \
    -Lbuild -lmathc -lm

./build/ultimate_fuzzer "$SEED"

echo "========================================================="
echo "  STRICT CLOSURE GATE PASSED"
echo "========================================================="
