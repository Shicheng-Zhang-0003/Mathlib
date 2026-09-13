# MathLib v12R2 — R2 Refinement Release

**Tag:** v12R2-refinement
**Date:** 2026-09-13
**Oracle:** 212 passed, 0 failed (all functions ≤ 5 ULP vs mpmath 80-digit ground truth)
**Gate:** Full closure gate passed (build, modular, edge, fuzz, boundary, oracle)

---

## What v12R2 Is

v12R2 is the refinement cycle following the v12A1 architectural evolution release.
v12A1 replaced approximations with the real thing. v12R2 fixes critical bugs and improves numerical accuracy.

## Key Fixes in v12R2

### Critical Correctness Fixes

**1. Portable `ml_sqrt`** — Removed x86 inline asm (`sqrtsd`) which violated C99, bypassed MXCSR rounding mode, and could flush subnormals. Now uses `__builtin_sqrt` everywhere.

**2. Fixed `ml_fmod`** — The original word-at-a-time algorithm was mathematically incorrect (integer modulus on significands). Replaced with proper IEEE-754 `fmod(x,y) = x - trunc(x/y)*y` using error-free multiplication via `ML_FMA`. Fixes `sin/cos` range reduction in EMBEDDED profile.

**3. Fixed `ml_ldexp_pure`/`ml_frexp_pure` shift UB** — Added bounds check `sig > (UINT64_MAX >> 1)` before `sig <<= 1` to avoid C99 undefined behavior on shifts ≥64.

**4. Fixed `ml_tanh`/`ml_atanh`/`ml_asinh` small-x threshold** — Changed `1e-4` → `1.5e-8`. At `x=1e-4`, Taylor truncation error was ~3,300,000 ULP. Now correctly returns `x` only when error < 1 ULP.

**5. Fixed `ml_atan` accuracy** — Added argument reduction to `|x| ≤ tan(π/12) ≈ 0.268` with 22 terms (up to x⁴³). Previous 18-term Taylor at `|x|≤0.5` had ~2000 ULP error at x=0.5. Now < 0.001 ULP.

**6. Fixed `ml_asin`/`ml_acos` cancellation near 1** — Replaced unstable formulas:
- `asin(x) = 2*atan(x/(1+sqrt(1-x²)))` → `π/2 - 2*atan(sqrt((1-x)/(1+x)))` for x > 0.9
- `acos(x) = π/2 - asin(x)` → `2*atan(sqrt((1-x)/(1+x)))` for x ≥ 0

**7. Fixed `ml_pow` integer exponent limit** — Raised `|y| ≤ 64` → `|y| ≤ 1023`. Binary exponentiation is exact up to overflow threshold.

**8. Fixed `ml_gamma_new` reflection formula** — Uses full `ML_PI_HI_D + ML_PI_LO_D` instead of `ML_PI_HI_D` only (~1 ULP improvement).

### Numerical Accuracy Improvements

**9. `ML_FMA` software fallback** — Replaced `(a*b)+c` (two roundings) with proper Dekker Two-Product + Two-Sum emulation for single-rounding semantics on non-FMA platforms.

**10. `ml_two_product` Dekker split** — Changed multiplier from `2²⁶+1` (67108865) to `2²⁷+1` (134217729) for exact 53-bit split.

**11. `ml_variance` Welford's algorithm** — Replaced two-pass catastrophic-cancellation algorithm with numerically stable online version.

**12. `ml_binomial_pmf` uses `lgamma`** — Replaced naive `sum log()` with `lgamma(n+1) - lgamma(r+1) - lgamma(n-r+1)`.

**13. `ml_polynomial_eval` uses `ML_FMA`** — Horner evaluation now fused (single rounding per step).

**14. All decimal constants → exact hex floats** — `ML_LN2_HI/LO`, `1/√2`, polynomial coefficients, thresholds (`ML_LOG_DBL_MAX`, `ML_LOG_UNDERFLOW`) now exact.

### Robustness

**15. Complex division** — Added `denom == 0` check in Smith's method to guard against catastrophic cancellation.

**16. Newton-Raphson** — Removed bogus `ml_fabs(dfx) < epsilon` check (epsilon is x-tolerance, not derivative threshold).

**17. Minimax header** — Renamed `maclaurin_*` → `taylor_*` with clarifying comments that these are Taylor, not minimax. True minimax in `minimax_coeffs.h` (DORMANT).

### Test Results

- **Oracle:** 212 passed, 0 failed (sin, cos, exp, log, gamma, lgamma, pow — all ≤ 5 ULP)
- **Edge tests:** 22 suites, all passed (441 assertions)
- **Fuzzers:** fuzz_god_mode, fuzz_boundary — all passed
- **Boundary gauntlet:** 25 passed, 0 failed
- **Soak:** 10,000 iterations available via `--soak`
- **Sanitizers:** ASan + UBSan clean

## Deferred to v12A2

- Minimax coefficient swap (coefficients generated and validated, dormant)
- Interval arithmetic / verified computing
- Automatic differentiation
- Special functions (Bessel, elliptic, hypergeometric)
- Adaptive quadrature with certified error estimates

## What's Not In Scope

- Universal correctly-rounded transcendental functions
- Identical output across every architecture
- Replacement of specialized high-precision libraries

These are design choices, not hidden defects.

---

*v11S shipped 2026-08-02. v12A1 A1 closure completed 2026-08-11. v12R2 refinement completed 2026-09-13.*
*The gamma nightmare is over. The critical bugs are fixed.*
