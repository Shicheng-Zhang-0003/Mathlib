#!/usr/bin/env bash
# v11S CLOSURE IP-20: edge test runner
set -u

cd "$(dirname "${BASH_SOURCE[0]}")/.."

CC=${CC:-gcc}

CFLAGS="-std=c99 -O3 -Wall -Wextra -Wconversion -Wshadow -Wpedantic -Werror -fno-fast-math -ffp-contract=off -Iinclude/mathlib -Isrc"

if [ "${MATHLIB_EDGE_SANITIZERS:-0}" = "1" ]; then
    CFLAGS="$CFLAGS -fsanitize=address,undefined -fno-omit-frame-pointer"
fi

SRC="
src/core.c
src/trig.c
src/exp_log.c
src/complex.c
src/linalg.c
src/fft.c
src/cpu_dispatch.c
src/combinatorics.c
src/quadratics.c
src/polynomial.c
src/numerical.c
src/statistics.c
src/integral.c
src/ode.c
src/optimization.c
src/quaternion.c
src/fixed_point.c
"

mkdir -p build

FAILED=0

for t in tests/test_edge_*.c; do
    name=$(basename "$t" .c)
    echo "--- $name ---"

    if ! $CC $CFLAGS -o "build/$name" "$t" $SRC -lm; then
        echo "❌ BUILD FAILED: $name"
        FAILED=1
        continue
    fi

    if ! "./build/$name"; then
        echo "❌ TEST FAILED: $name"
        FAILED=1
    fi
done

if [ "$FAILED" -eq 0 ]; then
    echo "🎉 ALL EDGE TESTS PASSED"
else
    echo "⚠️  SOME EDGE TESTS FAILED"
fi

exit $FAILED
