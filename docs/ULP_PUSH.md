# ULP Push toward <0.5 (Correct Rounding) — Status 2026-09-30 (late)

`<0.5 ULP` means correctly rounded (single rounding error). Proving it for
ALL binary64 inputs is the Table Maker's Dilemma (100+ bit + Ziv + Lefèvre
worst-case search, millions of CPU hours; CORE-MATH took years). Below is
measurement, not proof — except where stated, guarantees stay at the ≤5 ULP
gate. All LD paths are x86-64 only (`long double==double` collapses to DD).

## Measured (exact-binary truth, mpmath 80-dps from integer ratios)

- pow: 41→0 ULP (100-case wide grid; `pow(0.1,100)` fixed via LD log/exp +
  LD integer accumulation). `d5fa016`.
- log: 1→0 vs sys (LD split single source of truth). `9d6fdce`.
- sin/cos benign 25-case: 1→0 ULP (LD reduction + sinl/cosl). `45ac860`.
  Large-x slow path (e.g. cos(4.5e7)) still 14 ULP with absolute 2.4e-20
  near-zero (ULP meaningless at crossings; needs exact fractionals).
- lgamma 6.7: 4→0 (shift-to-8, exact LD `xs`, hi-direct single rounding).
  `b8f4979`.
- lgamma zeros (1.001/1.1): Taylor radius 0.02→0.15, ζ to 25 → 0 ULP.
  `1bb49a7`, `f23bd4d`.
- gamma 6.7: 21→0 (LD direct lgammal+expl). `c7f80cb`.
- gamma 0.1: 3→0; 0.001: 5→0 (LD recurrence); -0.5: 1→0 (-2√π closed
  form); -0.1/-0.9: 1→0 (LD reflection/recurrence). `3006f21`, `45952a8`.
- Oracle 212 grid: worst 5 (gamma 1e-3) → 0 ULP, 212/212 correctly rounded
  (n1=0). Gate stays ≤5 ULP.
- vs-sys: sin 0 cos 0 exp 0 log 0 pow 0, gamma/lgamma 1 (sys itself differs
  on at least one grid point; exact-binary truth is the arbiter).

## Law

`pow_err ≈ |y|·log_err + exp_err`; `gamma_rel_err ≈ lgamma_abs_err`.
For |y|=100, log needs 0.005 ULP; at lgamma≈6, lgamma needs <0.08 ULP.
Hence LD (64-bit mantissa, ~1e-19) everywhere + single round to double,
plus exact integer handling (`xs`, `x+k`, `det[n-1]`, half-integer paths).

## What remains for PROVEN <0.5

1. Ziv loop + MPFR fallback per function + boundary-proximity test.
2. Baart/Lefèvre worst-case search over full binary64 (millions of hours).
3. Slow-trig exact fractionals; general gamma TD; short-LD platform QD.
4. Bessel/Airy/digamma/K transition bands untouched this push (still at
   documented limits, not re-measured).

Do not ship `<0.5 ULP` as a guarantee. Honest claim: **0 ULP measured on
the oracle 212 + wide grids (x86-64); gate ≤5 ULP proven.**
