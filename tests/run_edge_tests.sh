#!/usr/bin/env bash
# v11S CLOSURE IP-20: edge test runner
set -u

cd "$(dirname "${BASH_SOURCE[0]}")/.."

CC=${CC:-gcc}

CFLAGS="-std=c99 -O3 -Wall -Wextra -Wconversion -Wshadow -Wpedantic -Wundef -Werror -fno-fast-math -ffp-contract=off -Iinclude/mathlib -Isrc"

if [ "${MATHLIB_EDGE_SANITIZERS:-0}" = "1" ]; then
    CFLAGS="$CFLAGS -fsanitize=address,undefined -fno-omit-frame-pointer"
fi

EDGE_SANITIZERS="${MATHLIB_EDGE_SANITIZERS:-0}"
LIBRARY_BUILD_DIR="${MATHLIB_BUILD_DIR:-build}"
LIBRARY_ARCHIVE="$LIBRARY_BUILD_DIR/libmathc.a"

# Canonical TU list: must match MATHLIB_SOURCES in CMakeLists.txt (59 TUs).
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
src/orthogonal.c
src/calculus.c
src/numbertheory.c
src/transforms.c
src/optim_n.c
src/ode_sys.c
src/spectral.c
src/stats_inv.c
src/sde.c
src/pde.c
src/harmonic.c
src/mcmc.c
src/manifold.c
src/info.c
src/analytic_nt.c
src/control.c
src/kelvin.c
src/contest.c
src/geometry.c
src/graph.c
src/modeling.c
src/measure.c
src/funcanal.c
src/topology.c
src/diffgeom.c
src/algebraic.c
src/math_extra.c
src/bigint.c
src/exact.c
src/interval.c
src/polyq.c
src/gfp.c
src/special.c
src/linalg_extra.c
src/calculus_extra.c
src/stats_extra.c
src/graph_extra.c
src/optim_extra.c
src/geom_extra.c
src/complexan.c
src/ode_extra.c
src/school.c
"

mkdir -p "${MATHLIB_EDGE_BUILD_DIR:-build}"

FAILED=0

for t in tests/test_edge_*.c; do
    name=$(basename "$t" .c)
    echo "--- $name ---"

    # Reuse the sanitizer-built CMake archive when available. The fallback
    # keeps this script usable on its own and for unsanitized source builds.
    if [ "$EDGE_SANITIZERS" = "1" ] && [ -f "$LIBRARY_ARCHIVE" ]; then
        LINK_INPUTS="$LIBRARY_ARCHIVE"
    else
        LINK_INPUTS="$SRC"
    fi

    if ! $CC $CFLAGS -o "${MATHLIB_EDGE_BUILD_DIR:-build}/$name" "$t" $LINK_INPUTS -lm; then
        echo "❌ BUILD FAILED: $name"
        FAILED=1
        continue
    fi

    if ! "${MATHLIB_EDGE_BUILD_DIR:-build}/$name"; then
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
