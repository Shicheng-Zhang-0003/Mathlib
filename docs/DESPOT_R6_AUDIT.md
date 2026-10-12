# Despot Round-6 Audit — fuzz upgrade + QR/fft fixes

All builds/tests ran in `/tmp/opencode` (repo tree untouched except the fix).
Seed below is `123456789` unless noted.

## 1. Fuzz coverage gap (the big find)

God-mode covered only the ~15 oldest TUs; 28 newer TUs (all of Batch-1,
contest, geometry, graph, modeling, doctorate) had ZERO randomized
fuzzing — deterministic edge pins only. Fixed by God-mode Part 2: ten new
sections, 51,600 checks with fixed per-iteration budgets:

| Section | Iters x checks | Content |
|---|---|---|
| contest identities | 800x7 = 5600 | Narayana=Catalan rows, hook/class trivials, Lucas parity, Josephus sim, Hanoi, Lucas-Fib |
| olympiad NT | 600x7 = 4200 | floor_sum/LTE/Pell/Tonelli/CF/vp naive differentials |
| graph invariants | 500x9 = 4500 | components, BFS, Dijkstra==BFS, Kruskal==n-1, topo validity, cycle parity |
| modeling forms | 500x9 = 4500 | stationary sum/residual, step sum, Little, SIR residual, conjugate exact |
| geometry predicates | 600x9 = 5400 | R>=2r, area==Heron, Euler/nine-point, hull, power, rect |
| doctorate spot | 400x20 = 8000 | Banach, Neumann, Hilbert, cycle Betti, Bareiss, Sturm/Descartes, Rodrigues, cubic/Vieta |
| stats/info | 600x8 = 4800 | entropy/KL/cross, logistic symmetry, softplus, chi2 roundtrip, sphere |
| spectral/ODE | 300x6 = 1800 | CG/QR residuals, trace, RK4, cubic roots |
| roundtrips | 400x7 = 2800 | Haar/DCT/FFT/spline/Poisson/heat |
| integer blizzard | 2000x5 = 10000 | egcd/modinv/crt2/nCr edge-biased |

Discipline: self-contained xorshift64* stream (Part 1 keeps `rand()`
untouched, so its stream and counts are unchanged); inputs constructed
valid (bounded rejection + deterministic fallback), so Part 2's total is
exactly 51,600 on every host and seed. Gate stays "0 failed", never a count.

Measured: Part 1 = 61,473 (seed 123456789), 61,393 (seed 1), 61,297
(seed 999) — this resolves the old 61393/61473 dispute on the record:
it was seed-dependence all along, not libc drift. New grand totals:
god-mode 113,073/0 (123456789), 112,993/0 (1), 112,897/0 (999).

Boundary gauntlet grew 25 -> 63 (contest/geometry/graph/modeling/
doctorate exact sections). Ultimate fuzzer grew 49,962 -> 89,462
(+19,500 double crash with Inf-aware predicates, +20,000 integer
blizzard). Stale banners fixed (v1.0/v12A1 -> v12R2).

## 2. Real bugs the new fuzz caught

- `ml_qr_iter_eig` (WRONG-ANSWER + non-convergence Kaboom, now fixed):
  the Givens sweep nested similarities instead of forming RQ (~20%
  stall on SPD incl. rotation blocks already in Schur form). Replaced
  with Hessenberg reduction + explicit shifted QR (MGS of p(H),
  full-block similarity). Two further bugs fell out during the fix:
  Householder beta off by v0^2 (trace drift) and MGS orthogonality
  loss near singular p(H) (1e-8 drift; cured by reorthogonalization).
  A final exact period-2 cycle (|p| == 1 on the whole spectrum) is
  broken by a single ad-hoc exceptional shift. Verified: 0
  non-convergences in 2000 general 2x2..6x6, det clean throughout,
  LAPACK-agreeing values (numpy cross-check).
- `ml_fft_is_supported` lied about 2^25+ (returned 1; execute refuses):
  now pow2 AND <= 2^24, pinned in boundary + edge spirit.
- Test-side (library correct): CF k=0 bound, LTE a==b, Neumann
  row-sum contract, Descartes-positive-only, descending cycle edge,
  logistic closed endpoints, CHECK_NEAR strict-eps, magic/LL pins.

## 3. Stale labels swept (all verified by grep)

`BUILD_AND_INSTALL` 34->43 TUs, `MAINTAINERS` 25/25->35/35 suites,
`KNOWN_LIMITATIONS` freeze note 34->43, fuzzer banners v1.0/v12A1->v12R2.
Historical audit reports and release-notes history left intact; only
living docs updated. Header-per-TU check: all 43 TUs have headers
(`cpu_dispatch`/`fft` under their own names, by design).

## 4. Gate on the R6 tree (measured)

- oracle 212/0 worst 0 ULP (gate ≤5 ULP), modular 4/4, smoke 30/0,
  boundary 63/0, god 113073/0 (0 failed), ultimate 89462/0,
  edge 35/35, 3-profile `-Werror` + ASan/UBSan clean.
