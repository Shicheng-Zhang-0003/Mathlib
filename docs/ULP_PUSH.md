# ULP Push toward <0.5 (Correct Rounding) — Status 2026-09-30

`<0.5 ULP` means correctly rounded (single rounding error). This is the
Table Maker's Dilemma: proving it for ALL binary64 inputs requires
~100+ bit intermediates + Ziv's onion-peeling fallback + exhaustive
worst-case search (Lefèvre-style, millions of CPU hours per function).
CORE-MATH took 5+ years for ~30 functions. No single session proves it.

## Measured baselines (vs sys + mpmath 80-dps)

- sin/cos/exp/log (benign grid): 0-1 ULP vs sys. Already near-rounding.
- pow: **41 ULP** on `pow(0.1,100)` (integer-path O(n) accumulation +
  DD log 0.38 ULP ×|y|=100 amplification). Fixed this session.
- gamma: 21 ULP at 6.7, 3-5 ULP at 0.1/0.001. lgamma: 4 ULP at 6.7.
- oracle worst: gamma(1e-3) 5 ULP. Gate is ≤5 ULP, not ≤0.5.

## Fixed this session (`d5fa016`)

- atanh `lc`: hex-exact `2/(2k+1)` (auditable, bit-identical).
- pow general: `L=y*logl(x)` in 80-bit LD + `expl` + single round, with
  DD-vs-LD Ziv blend near boundaries. `pow(0.1,100)` 41→~0 ULP.
- pow integer `|y|≤1023`: LD accumulation then round-once
  (`|y|·2^-65` error, <0.05 ULP even at 1023). Exact cases stay exact.
- pow wide (100 cases) vs mpmath80: worst **41→0 ULP**.
- No regressions: oracle 212/0, smoke 30/0, core 5/0.

## Law

`pow_err ≈ |y|·log_err + exp_err`. For `|y|=100`, log needs 0.005 ULP
to hold pow <0.5. DD (106-bit) can do it only if the log itself is
<0.01 ULP; ours was 0.38 ULP → LD path required. Same law governs
`gamma = exp(lgamma)`: `gamma_err ≈ |lgamma|·lgamma_rel + exp_err`;
at 6.7 (`lgamma≈6`) lgamma needs <0.08 ULP for gamma <0.5.

## What remains for PROVEN <0.5

1. Per-function Ziv loop: fast DD/LD path → boundary proximity test
   (`|frac(mantissa)-0.5| < ε`) → MPFR/long-double-128 fallback.
   Mechanism exists in pow (LD-vs-DD blend); needs MPFR backend + proof.
2. Worst-case search: Baart/Lefèvre hunters for sin/cos/exp/log/pow/gamma
   over full binary64 (or documented reduced domains). Millions of core-hours.
3. Gamma/LD Stirling: `x≥8` via LD `logGamma` + `expl` + single round;
   half-integer exact path already 0 ULP; `x<0.5` recurrence + reflection
   need TD (triple-double) to hold 0.08 ULP. Not yet implemented.
4. Range-reduction proof: Payne-Hanek already few-ULP full-range; proving
   <0.5 needs exact `2/π` table + formal error bound (currently empirical).
5. `long double==double` platforms (MSVC/ARM): LD paths collapse to DD;
   proven <0.5 impossible there without software TD/QD. Document as
   platform-scoped guarantee (x86-64 only).

Do not claim `<0.5 ULP` on the release until 1+2 are done per function.
Current honest claim: **pow <0.5 on tested grid (100 mpmath cases 0 ULP),
core 0-1 ULP typical, gate ≤5 ULP proven, gamma documented 21 ULP at 6.7.**
