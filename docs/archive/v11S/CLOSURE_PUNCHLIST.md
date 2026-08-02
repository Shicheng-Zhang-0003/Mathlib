# MathLib v11S Closure Punchlist

This file is the authoritative short-term closure queue for v11S.

Current operational state:

> Late v11A3 / v11S release-candidate territory.
> Not yet a true stable v11S release.

---

## P0 Blockers

These must be fixed before v11S can be treated as stable.

### 1. `ml_pow` is not closure-grade

Current implementation is effectively:

    exp(y * log(x))

That is insufficient for:

- negative base with integer exponent,
- signed zero behavior,
- `0^0`,
- `1^NaN`,
- `NaN^0`,
- infinite exponent behavior,
- infinite base behavior.

Required action:

- add integer-exponent classification helpers,
- implement explicit special-case decision tree,
- preserve signed zero where appropriate,
- keep overflow/underflow behavior consistent with `ml_exp`.

---

### 2. `ml_cplx_power` lacks required special cases

Current implementation uses only:

    exp(b * log(a))

That cannot satisfy required special cases such as:

- `0^0 = 1`
- `1^NaN = 1`
- `NaN^0 = 1`
- `0^i = NaN`
- `0^-1 = inf`
- `inf^i = NaN`
- negative real principal-branch powers

Required action:

- add explicit complex special-case decision tree,
- preserve principal-branch behavior for finite nonzero bases,
- document conservative infinite-exponent contract.

---

### 3. `ml_cplx_arg` has a quadrant bug

The third-quadrant path currently divides by the wrong component.

Required action:

- replace hand-rolled logic with `ml_atan2(imag, real)`.

This also aligns implementation with the documented claim that complex arg is atan2-based.

---

### 4. `ml_exp` needs explicit NaN / infinity guards

`ml_exp` must classify NaN and infinity before range reduction.

Required action:

- return NaN for NaN input,
- return positive infinity for positive infinity,
- return positive zero for negative infinity,
- avoid casting NaN to integer.

---

### 5. `ml_log` must handle positive infinity

Current behavior can fall through into the finite path and produce NaN.

Required action:

- return positive infinity for positive infinity,
- preserve negative infinity as domain NaN,
- preserve zero as negative infinity.

---

### 6. Official verification must run edge tests

The edge suite exists but is not part of the main closure gate.

Required action:

- use `closure_gate.sh` as the strict gate,
- later fold the edge suite into the official verification path,
- do not claim closure while edge tests are skipped.

---

## P1 High-Priority Issues

These should be fixed soon after P0.

### 1. `ml_sinh` loses tiny inputs

Add a small-input branch consistent with other hyperbolics.

---

### 2. AVX2 batch `rsqrt` needs semantic guarding

The scalar fallback handles NaN / non-positive inputs.

The AVX2 path should either:

- fall back to scalar semantics for exceptional inputs,
- or explicitly document the restricted fast-path contract.

Preferred:

- guard the AVX2 path and fall back for non-finite / non-positive / subnormal inputs.

---

### 3. Matrix and tensor indexing should use wide size arithmetic internally

Public API may remain `int`, but internal offset arithmetic should avoid signed overflow.

Required action:

- convert matrix loops to size-based indexing internally,
- convert tensor offset arithmetic to wide size arithmetic.

---

### 4. CI should run deterministic fuzzing and edge tests

Required action:

- add edge tests to CI or create a dedicated closure gate job,
- use fixed fuzzer seeds in CI.

---

### 5. Documentation must match implementation

After source fixes:

- update `API_STATUS.md`,
- update `KNOWN_LIMITATIONS.md`,
- update `V11S_CLOSURE_SUMMARY.md`,
- ensure README verification flow matches actual verification scripts.

---

## Script Sequence

### Script 00

`00_closure_bootstrap.py`

Creates:

- closure punchlist,
- team doctrine,
- script-only change policy,
- strict closure gate.

This script does not modify math source.

---

### Script 01

`01_p0_exp_log_pow_complex.py`

Intended to fix:

- `src/internal/pow_util.h`
- `src/exp_log.c`
- `src/complex.c`

This is the first real source-change script.

---

### Script 02

`02_p1_simd_index_ci_docs.py`

Intended to fix:

- SIMD batch guarding,
- internal size arithmetic,
- CI determinism,
- documentation alignment.

---

## Closure Rule

v11S is not stable until:

1. all P0 items are fixed,
2. edge tests pass,
3. sanitizers pass,
4. documentation matches code,
5. strict closure gate passes.


<!-- MATHLIB_CLOSURE_P1_LOG -->
## Closure Log

- P0 source fixes applied and strict closure gate passed.
- P1 safety/process fixes applied by `02_p1_simd_index_ci_docs.py`.
- AVX2 batch rsqrt now falls back to scalar semantics for exceptional inputs.
- Matrix/tensor indexing now uses wide size arithmetic internally.
- Official verification script now runs edge tests and boundary gauntlet.
- CI now runs edge tests and deterministic quick fuzz.

<!-- MATHLIB_CLOSURE_P2_P0_1_LOG -->
## Red Team P2 P0-1

- Fixed `ml_ldexp_pure()` false overflow for subnormal significands.
- Normalization now occurs before overflow judgment.
- Added regression suite `tests/test_edge_redteam_p0_1.c`.

<!-- MATHLIB_CLOSURE_P2_P0_2_LOG -->
## Red Team P2 P0-2

- Replaced `ml_round()` with exact decomposition-based round-half-away-from-zero.
- Fixed misrounding near half-integers caused by the old `x +/- 0.5` trick.
- Added regression suite `tests/test_edge_redteam_p0_2.c`.

<!-- MATHLIB_CLOSURE_P2_P0_3_LOG -->
## Red Team P2 P0-3

- Fixed `ml_exp()` premature overflow threshold.
- Fixed `ml_exp()` premature underflow threshold.
- Added regression suite `tests/test_edge_redteam_p0_3.c`.

<!-- MATHLIB_CLOSURE_P2_P0_3_FIX_LOG -->
## Red Team P2 P0-3 Fix

- Corrected earlier P0-3 marker collision.
- `ml_exp()` now uses `ML_LOG_DBL_MAX` and `ML_LOG_UNDERFLOW`.
- Old coarse thresholds `709.78` and `-745.13` are replaced.
- Regression suite `tests/test_edge_redteam_p0_3.c` remains active.

<!-- MATHLIB_CLOSURE_P2_P0_5_LOG -->
## Red Team P2 P0-5

- Hardened double CORDIC against NaN / Inf inputs.
- Hardened `fuzz_god_mode` fixed CORDIC test to avoid casting non-finite doubles to integers.
- Added regression suite `tests/test_edge_redteam_p0_5.c`.

<!-- MATHLIB_CLOSURE_P2_LOG_EXP_ROUNDTRIP_LOG -->
## Red Team P2 Log/Exp Roundtrip Fix

- Replaced over-tight `log(exp(x))` fuzzer tolerance with ULP-aware tolerance.
- Added `ml_fuzz_ulp_distance()` to `fuzz_god_mode.c`.
- Improved `ml_log()` final reconstruction using `ML_FMA`.
- Added deterministic regression suite `tests/test_edge_redteam_p2_logexp.c`.

<!-- MATHLIB_CLOSURE_P2_P0_4_LOG -->
## Red Team P2 P0-4

- Fixed premature `ml_sinh()` overflow threshold.
- Fixed premature `ml_cosh()` overflow threshold.
- Hyperbolic overflow boundary now uses `log(DBL_MAX) + log(2)`.
- Near-overflow evaluation now uses shifted exponential form `exp(x - ln2)`.
- Added regression suite `tests/test_edge_redteam_p0_4.c`.

<!-- MATHLIB_REDP2_P0_1_LOG -->
## Red Team P2 P0-1

- Fixed `ml_ldexp_pure()` false overflow for subnormal significands.
- Normalization now occurs before overflow judgment.
- Added regression suite `tests/test_edge_redteam_p0_1.c`.

<!-- MATHLIB_REDP2_P0_2_LOG -->
## Red Team P2 P0-2

- Replaced `ml_round()` with exact decomposition-based round-half-away-from-zero.
- Fixed misrounding near half-integers caused by the old `x +/- 0.5` trick.
- Added regression suite `tests/test_edge_redteam_p0_2.c`.

<!-- MATHLIB_REDP2_P1_1_LOG -->
## Red Team P2 P1-1

- Aligned `ml_simd_batch_rsqrt()` exceptional fallback with scalar `ml_fast_rsqrt()`.
- Batch fallback now preserves scalar fast-math semantics for zero, negative, inf, NaN, and subnormal inputs.
- Added regression suite `tests/test_edge_redteam_p1_1.c`.

<!-- MATHLIB_REDP2_P1_1_TEST_TOLERANCE_LOG -->
## Red Team P2 P1-1 Test Tolerance Fix

- Corrected `tests/test_edge_redteam_p1_1.c` tolerance for the approximate fast rsqrt path.
- Positive-normal batch rsqrt now uses a fast-math-appropriate relative tolerance.
- Exceptional-input semantic checks remain strict.
- No SIMD kernel behavior was changed.
