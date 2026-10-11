# Contest Guide (`ml_contest.h`, 36 functions)

Extravagant contest mathematics, tiered low -> high. Extends the existing
contest core (combinatorics, numbertheory, arith/geom sums, Nim, Karamata,
partition) with the classic missing staples. 35th TU (`src/contest.c`).

## Tier 1 — school / AMC

| Function | Theorem | Sentinel |
|---|---|---|
| `ml_triangular/pentagonal/hexagonal` | closed forms | `UINT64_MAX` overflow |
| `ml_is_triangular` | `8m+1` odd square (integer Newton, no FP) | `0` = no |
| `ml_double_factorial` | `(-1)!!=0!!=1` | `UINT64_MAX` |
| `ml_binom_odd` | Lucas mod 2 (`k&~n==0`) | `0` = even/invalid |
| `ml_stars_bars` | `C(s+b,b)` via `ml_ncr` | `UINT64_MAX` |
| `ml_hanoi_moves` | `2^n-1` (`64` capped: it would equal the sentinel) | `UINT64_MAX` |
| `ml_josephus` | iterative, `O(n)`, cap `1e8` | `UINT64_MAX` |
| `ml_lucas/pell/tribonacci` | checked recurrences | `UINT64_MAX` |

## Tier 2 — AIME / olympiad

Bell (Aitken), Stirling 1st/2nd (DP tables), Eulerian, Lah
(`C(n-1,k-1)n!/k!`), Motzkin (`(n+2)M=(2n+1)M'+3(n-1)M''`), large Schröder,
Narayana (`rows sum to Catalan_n`, pinned), distinct partitions DP
(`n<=200`; Euler: equals odd partitions), Goldbach counts (heap sieve,
`n<=1e6`), Sophie Germain, Lucas-Lehmer (`p<=61`), Tonelli-Shanks (least
root), tau/sigma (`n<=1e12`), Frobenius `ab-a-b` + count `(a-1)(b-1)/2`.

## Tier 3 — Putnam / college

Harmonic via `psi(n+1)+gamma` (reuses staged digamma, ~1 ULP), Basel partial
(Kahan LD, limit `pi^2/6` documented), Faulhaber via Bernoulli `B0..B30`
(`(-1)^j` pairs the `B1=-1/2` convention; LD ~1e-15), Wythoff cold
(`floor(k phi)`, status return), subtraction-game Grundy (mex DP, cap 4096),
Schur degree-1 verifier, AM-GM instance check, cubic Vieta sums.

## Exactness policy

u64 = exact or `UINT64_MAX`; verdicts fail closed to `0`; doubles follow
the `ml_sum_k` exactness-guard pattern where exact, else LD
regression-guard tier (no oracle ULP claimed — `PRECISION_CONTRACT` rule 3).
Pinned by `tests/test_edge_contest.c` (96 assertions).
