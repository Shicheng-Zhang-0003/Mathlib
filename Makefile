CC = gcc
# Strict C99 warnings, dependency tracking, and IEEE-754 determinism
# NOTE: -Wfloat-equal is intentionally omitted (The Float-Equal Paradox).
CFLAGS = -std=c99 -O3 -Wall -Wextra -Wconversion -Wshadow -Wpedantic -Wundef -Werror -fno-fast-math -ffp-contract=off -Iinclude/mathlib -Isrc -MMD -MP
LDLIBS = -lm

SRC = src/core.c src/trig.c src/exp_log.c src/complex.c src/linalg.c src/fft.c src/cpu_dispatch.c src/combinatorics.c src/quadratics.c src/polynomial.c src/numerical.c src/statistics.c src/integral.c src/ode.c src/optimization.c src/quaternion.c src/fixed_point.c src/orthogonal.c src/calculus.c src/numbertheory.c src/transforms.c src/optim_n.c src/ode_sys.c src/spectral.c src/stats_inv.c src/sde.c src/pde.c src/harmonic.c src/mcmc.c src/manifold.c src/info.c src/analytic_nt.c src/control.c src/kelvin.c src/contest.c src/geometry.c src/graph.c src/modeling.c src/measure.c src/funcanal.c src/topology.c src/diffgeom.c src/algebraic.c src/math_extra.c src/bigint.c src/exact.c src/interval.c src/polyq.c src/gfp.c src/special.c src/linalg_extra.c src/calculus_extra.c src/stats_extra.c src/graph_extra.c src/optim_extra.c src/geom_extra.c src/complexan.c src/ode_extra.c src/school.c
DEPS = $(SRC:.c=.d)

TEST_SRC = tests/test.c
BENCH_SRC = benchmarks/bench.c
BENCH_LIBM_SRC = benchmarks/bench_vs_libm.c
FUZZ_GOD_SRC = tests/fuzz_god_mode.c
FUZZ_BOUND_SRC = tests/fuzz_boundary_gauntlet.c
ORACLE_SRC = tests/test_oracle.c

# Align output paths with CMake and soak_test.sh
OUT_DIR = build

all: $(OUT_DIR)/test $(OUT_DIR)/bench $(OUT_DIR)/fuzz_god_mode $(OUT_DIR)/fuzz_boundary $(OUT_DIR)/oracle_check

$(OUT_DIR):
	mkdir -p $(OUT_DIR)

$(OUT_DIR)/test: $(TEST_SRC) $(SRC) | $(OUT_DIR)
	$(CC) $(CFLAGS) -o $@ $(TEST_SRC) $(SRC) $(LDLIBS)

$(OUT_DIR)/bench: $(BENCH_SRC) $(SRC) | $(OUT_DIR)
	$(CC) $(CFLAGS) -o $@ $(BENCH_SRC) $(SRC) $(LDLIBS)

# Head-to-head timing vs system libm (not part of `all`; run: make bench-vs-libm)
$(OUT_DIR)/bench_vs_libm: $(BENCH_LIBM_SRC) $(SRC) | $(OUT_DIR)
	$(CC) $(CFLAGS) -o $@ $(BENCH_LIBM_SRC) $(SRC) $(LDLIBS)

.PHONY: bench-vs-libm
bench-vs-libm: $(OUT_DIR)/bench_vs_libm
	./$(OUT_DIR)/bench_vs_libm

$(OUT_DIR)/fuzz_god_mode: $(FUZZ_GOD_SRC) $(SRC) | $(OUT_DIR)
	$(CC) $(CFLAGS) -o $@ $(FUZZ_GOD_SRC) $(SRC) $(LDLIBS)

$(OUT_DIR)/fuzz_boundary: $(FUZZ_BOUND_SRC) $(SRC) | $(OUT_DIR)
	$(CC) $(CFLAGS) -o $@ $(FUZZ_BOUND_SRC) $(SRC) $(LDLIBS)

$(OUT_DIR)/oracle_check: $(ORACLE_SRC) $(SRC) | $(OUT_DIR)
	$(CC) $(CFLAGS) -DMATHLIB_HAS_ORACLE_DATA -o $@ $(ORACLE_SRC) $(SRC) $(LDLIBS)

clean:
	rm -f $(OUT_DIR)/test $(OUT_DIR)/bench $(OUT_DIR)/bench_vs_libm $(OUT_DIR)/fuzz_god_mode $(OUT_DIR)/fuzz_boundary $(OUT_DIR)/oracle_check
	rm -f $(DEPS)

.PHONY: all clean

# Modular CI/CD Tests
TESTS = $(OUT_DIR)/test_core $(OUT_DIR)/test_trig $(OUT_DIR)/test_linalg $(OUT_DIR)/test_dsp

$(OUT_DIR)/test_core: tests/test_core.c $(SRC) | $(OUT_DIR)
	$(CC) $(CFLAGS) -o $@ $< $(SRC) $(LDLIBS)

$(OUT_DIR)/test_trig: tests/test_trig.c $(SRC) | $(OUT_DIR)
	$(CC) $(CFLAGS) -o $@ $< $(SRC) $(LDLIBS)

$(OUT_DIR)/test_linalg: tests/test_linalg.c $(SRC) | $(OUT_DIR)
	$(CC) $(CFLAGS) -o $@ $< $(SRC) $(LDLIBS)

$(OUT_DIR)/test_dsp: tests/test_dsp.c $(SRC) | $(OUT_DIR)
	$(CC) $(CFLAGS) -o $@ $< $(SRC) $(LDLIBS)

test_modular: $(TESTS)
	@cd $(OUT_DIR) && python3 ../tests/run_tests.py

clean_modular:
	rm -f $(TESTS)

# Include auto-generated dependency files
-include $(DEPS)
