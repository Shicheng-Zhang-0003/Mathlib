#!/usr/bin/env bash
# v11S CLOSURE IP-19: deterministic soak test (tests/ path-safe version)
set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$ROOT_DIR"

echo "========================================================="
echo "  MATHLIB v11S: DETERMINISTIC SOAK TEST"
echo "========================================================="

BASE_SEED=${MATHLIB_SOAK_SEED:-1000000}
ITERATIONS=${MATHLIB_SOAK_ITERATIONS:-10000}
FAILED=0

mkdir -p tests/regression/artifacts

for i in $(seq 1 "$ITERATIONS"); do
    SEED=$((BASE_SEED + i))
    LOG="tests/regression/artifacts/soak_${SEED}.log"

    if ! ./build/fuzz_god_mode "$SEED" > "$LOG" 2>&1; then
        echo "❌ FAILURE at iteration $i (seed $SEED)"
        python3 tests/run_fuzz_vault.py --logfile "$LOG" --fuzzer ./build/fuzz_god_mode || true
        FAILED=1
        break
    fi

    rm -f "$LOG"

    if [ $((i % 1000)) -eq 0 ]; then
        echo "  [SOAK] $i / $ITERATIONS iterations passed..."
    fi
done

if [ $FAILED -eq 0 ]; then
    echo "🎉 SOAK TEST PASSED: $ITERATIONS iterations, 0 failures."
else
    echo "⚠️  SOAK TEST FAILED."
fi
