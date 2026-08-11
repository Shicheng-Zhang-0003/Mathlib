# MathLib v12A1 — A1 Closure Release

**Tag:** v12A1-closure
**Date:** 2026-08-11
**Oracle:** 212 passed, 0 failed (all functions ≤ 5 ULP vs mpmath 80-digit ground truth)
**Gate:** Full closure gate passed (build, modular, edge, fuzz, boundary, oracle, ultimate fuzzer)

---

## What v12A1 Is

v12A1 is the architectural evolution cycle following the v11S stable release.
v11S proved the foundations. v12A1 replaces approximations with the real thing.

## Key Achievements

### Payne-Hanek Range Reduction (V6)
The 1e15 wall is gone. `sin(1e50)`, `cos(1e300)` produce finite, correct results.
Per-term integer/fractional separation with Kahan accumulation and full 66-entry
table sweep. Validated by oracle at 1e10 through 1e300.

### Hybrid Gamma Function
The gamma nightmare is dead. After 16 scripts (30→45), the implementation is:

- *****Lanczos DD** (g=7, n=9) for 0.5 ≤ x < 8
- **Stirling DD** (12-term correction) for x ≥ 8
- **Exact half-integer formulas** (product/sum, no Lanczos approximation error)
- **1-step recurrence** for x < 0.5: `lgamma(x) = lgamma(x+1) − log(x)`
- **Reflection formula** with compensated sinpi for negative arguments

All paths validated ≤ 5 ULP against mpmath ground truth, including the
nightmare inputs gamma(0.1), gamma(0.01), gamma(0.001), and their lgamma
counterparts.

### Extended-Precision pow
Dekker-split log with FMA-compensated product. Integer-exponent fast path
via binary exponentiation. Full IEEE-754 special-case decision tree.

### Error-Free Cody-Waite Reductions
`ml_exp` uses FMA-based 2-term Cody-Waite with canonical ln(2) split.
`ml_log` uses compensated reconstruction with `ML_LN2_HI` / `ML_LN2_LO`.

### Validated Maclaurin Kernels
The active trig/exp/log kernels are Maclaurin (Taylor) series, validated
≤ 5 ULP under FMA Horner evaluation. True minimax coefficients are generated,
validated, and stored in `minimax_coeffs.h` (DORMANT). The swap is deferred
to v12A2.

## The Gamma Saga (Scripts 30–45)

The gamma function went through 16 scripts to reach ≤ 5 ULP:

| Script | Fix |
|--------|-----|
| 30 | Payne-Hanek V6 + Stirling V5 (depth 16, corrected B₂₂/B₂₄) |
| 31 | Lanczos + Stirling hybrid (V6) |
| 32 | Half-integer root cause fix (exact formulas, V7) |
| 33 | Product-then-divide → DD log-subtract-then-exp |
| 37 | gamma_new log-space fix (script 33 missed ml_gamma_new) |
| 38 | exp_dd second-order Taylor (dropped L.lo²/2 = 10 ULP) |
| 39–40 | Recurrence depth 8 → 16 |
| 44 | Direct Lanczos dispatch (kill the recurrence) |
| 45 | 1-step recurrence: lgamma(x) = lgamma(x+1) − log(x) |

Root cause of the final failures: catastrophic cancellation in the multi-step
recurrence formula. `lgamma(x+16) − Σlog(x+k)` subtracts two large DD values
(~28 − ~26) to get a small result (~2.25), amplifying rounding errors by ~12.5×.
The 1-step recurrence has cancellation factor ~1.04.

## Test Results

- **Oracle:** 212 passed, 0 failed (sin, cos, exp, log, gamma, lgamma, pow)
- **Edge tests:** 22 suites, all passed (441 assertions)
- **Fuzzers:** fuzz_god_mode, fuzz_boundary, ultimate_fuzzer — all passed
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

*v11S shipped 2026-08-02. v12A1 A1 closure completed 2026-08-11.*
*The gamma nightmare is over.*
