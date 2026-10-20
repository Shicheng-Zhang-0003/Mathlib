# Contest Guide (`ml_contest.h`, 42 functions)

Extravagant contest mathematics, tiered low -> high. Extends the existing
contest core (combinatorics, numbertheory, arith/geom sums, Nim, Karamata,
partition) with the classic missing staples. 35th TU (`src/contest.c`).

## Tier 0 — recreational (AMC 8/10)

Digit sums, digital roots (`1+(n-1) mod 9`), palindrome checks, Collatz
total stopping time (`UINT64_MAX` past u64/10^6 steps — no claim on the open
conjecture), normal magic-square verification (`n<=32`).

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
Pinned by `tests/test_edge_contest.c` (96 assertions) and
`tests/test_edge_olympiad_nt.c` (39 assertions: LTE, floor-sum, Legendre,
continued fractions, Pell, recreational).

## Olympiad NT annex (`ml_numbertheory.h`)

Legendre's `v_p(n!)`, AtCoder floor-sum (`__uint128_t`, same-file precedent
as `ml_crt2`), LTE (`v_p(a^n-b^n)` with odd/even-2 cases and precondition
fail-loud `UINT64_MAX`), `sqrt(D)` continued-fraction period (`0` = square,
`-1` = invalid, `-2` = past the 10000 bound), convergents, and the Pell
fundamental solution (first convergent hit wins; `D=61` gives the famous
1766319049/226153980; beyond-u64 solutions report `OVERFLOW`).
