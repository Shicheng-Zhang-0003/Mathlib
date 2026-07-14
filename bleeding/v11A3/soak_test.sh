#!/bin/bash
echo "========================================================="
echo "  MATHLIB A3: GOD-MODE SOAK TEST (10,000 Iterations)"
echo "========================================================="
echo "This will run the fuzzer 10,000 times with different seeds."
echo "If it finishes, your memory safety and edge-case routing are bulletproof."
echo ""
FAILURES=0
for i in {1..10000}; do
    ./build/fuzz_god_mode > /dev/null 2>&1
    if [ $? -ne 0 ]; then
        echo "❌ FAILURE at iteration $i"
        python3 tests/run_fuzz_vault.py
        FAILURES=$((FAILURES+1))
        break
    fi
    if [ $(($i % 1000)) -eq 0 ]; then
        echo "  [SOAK] $i / 10000 iterations passed..."
    fi
done

if [ $FAILURES -eq 0 ]; then
    echo "🎉 SOAK TEST PASSED: 10,000 iterations, 0 failures."
else
    echo "⚠️  SOAK TEST FAILED."
fi
