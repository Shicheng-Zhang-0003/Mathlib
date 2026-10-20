# MathLib v12R2 — R2 Refinement Release (Public: V1.2-RC2)

**Internal tag:** v12R2
**Public tag:** V1.2-RC2 (`V1.2-RC2 == v12R2`, same commit)
**Date:** 2026-09-30
**Oracle:** 412 passed, 0 failed, worst 0 ULP (26 families: 7 core sin, cos, exp, log, gamma, lgamma, pow + 19 R6 expansion, correctly rounded on the grid vs mpmath 50-dps; gate stays ≤5 ULP, see `docs/PRECISION_CONTRACT.md`)
**Gate:** Full closure gate passed (build, modular, edge 23 suites, fuzz, boundary, oracle, sanitizers)

## Unreleased v12R3 Development Audit (2026-10-19)

The working tree is now version 12.3.0 with 59 translation units. It adds
experimental exact/numerical extensions, a pre-university utility module, and audit fixes for integer conversion,
remainder quotient bits, gamma sign reporting, large negative zeta zeros,
fixed-capacity bigint arithmetic, and special-function boundary behavior.
These additions have directed regression coverage but no oracle-tier accuracy
certification. The v12R2 release measurements below do not apply to them.

The 2026-10-19 core accuracy recheck added 15 atan reference vectors and
midpoint-tail regressions for `expm1` and `log1p`. The current oracle has 406
mpmath vectors across 27 families plus 21 large-argument checks; all 427
checks pass, with a worst reference error of 1 ULP at `gamma(50)`. On a
deterministic 9,508-input atan sample, the extended-precision `atanl` path
reduced 1-ULP results from 69 to 3. The accuracy audit now passes 378 checks.

---

## What v12R2 Is

v12R2 is the refinement cycle following the v12A1 architectural evolution release.
v12A1 replaced approximations with the real thing. v12R2 fixes critical bugs, improves numerical accuracy, achieves full thread-safety, and adds 12 Batch-1 modules + kelvin (34 TUs total).

## Key Fixes in v12R2

### Critical Correctness Fixes

**1. Portable `ml_sqrt`** — Removed x86 inline asm (`sqrtsd`) which violated C99, bypassed MXCSR rounding mode, and could flush subnormals. Now uses `__builtin_sqrt` everywhere.

**2. Fixed `ml_fmod`** — The original word-at-a-time algorithm was mathematically incorrect (integer modulus on significands). Replaced with exact integer-significand long-division (`ax=sx*2^ex`, `rem=((sx*2^d)%sy)*2^ey`, no FP ops, sign-of-x, RNE subnormal handling). Fixes `sin/cos` range reduction in EMBEDDED profile.

**3. Fixed `ml_ldexp_pure`/`ml_frexp_pure` shift UB** — Added bounds check `sig > (UINT64_MAX >> 1)` before `sig <<= 1` to avoid C99 undefined behavior on shifts ≥64.

**4. Fixed `ml_tanh`/`ml_atanh`/`ml_asinh` small-x threshold** — Changed `1e-4` → `1.5e-8`. At `x=1e-4`, Taylor truncation error was ~3,300,000 ULP. Now correctly returns `x` only when error < 1 ULP.

**5. Fixed `ml_atan` accuracy** — Added argument reduction to `|x| ≤ tan(π/12) ≈ 0.268` with 22 terms (up to x⁴³). Previous 18-term Taylor at `|x|≤0.5` had ~2000 ULP error at x=0.5. Now < 0.001 ULP.

**6. Fixed `ml_asin`/`ml_acos` cancellation near 1** — Replaced unstable formulas:
- `asin(x) = 2*atan(x/(1+sqrt(1-x²)))` → `π/2 - 2*atan(sqrt((1-x)/(1+x)))` for x > 0.9
- `acos(x) = π/2 - asin(x)` → `2*atan(sqrt((1-x)/(1+x)))` for x ≥ 0

**7. Fixed `ml_pow` integer exponent limit** — Raised `|y| ≤ 64` → `|y| ≤ 1023`. Binary exponentiation is exact up to overflow threshold.

**8. Fixed `ml_gamma_new` reflection formula** — `ml_lgamma` uses full `ML_PI_HI_D + ML_PI_LO_D` via `ml_log_pi_dd`; `ml_gamma_new` uses HI/LO double-double division `q=HI/(s*G) + (FMA(-q,s*G,HI)+LO)/(s*G)` (~0.5 ULP improvement over HI-only).

### Numerical Accuracy Improvements

**9. `ML_FMA` software fallback** — Replaced `(a*b)+c` (two roundings) with proper Dekker Two-Product + Two-Sum emulation for single-rounding semantics on non-FMA platforms.

**10. `ml_two_product` Dekker split** — Changed multiplier from `2²⁶+1` (67108865) to `2²⁷+1` (134217729) for exact 53-bit split.

**11. `ml_variance` Welford's algorithm** — Replaced two-pass catastrophic-cancellation algorithm with numerically stable online version.

**12. `ml_binomial_pmf` uses `lgamma`** — Replaced naive `sum log()` with `lgamma(n+1) - lgamma(r+1) - lgamma(n-r+1)`.

**13. `ml_polynomial_eval` uses `ML_FMA`** — Horner evaluation now fused (single rounding per step).

**14. `ml_polynomial_newton`** — Removed bogus `fabs(dfx) < epsilon` check (epsilon is x-tolerance, not derivative threshold; only exact `dfx==0` aborts).

**15. All decimal constants → exact hex floats** — `ML_LN2_HI/LO`, `1/√2`, polynomial coefficients, thresholds (`ML_LOG_DBL_MAX`, `ML_LOG_UNDERFLOW`) now exact.

### Robustness & Thread-Safety

**16. Thread-safety: all mutable static scratch buffers removed (2026-09-27 despot audit)** — Verified by `grep "static ...\[" src/*.c` returning no per-call mutable state: `optim_n.c`, `mcmc.c`, `pde.c`, `manifold.c`, `harmonic.c`, `spectral.c`, `calculus.c` (spline Thomas), `linalg.c` (Jacobi `W`), `info.c` (`ml_mi_discrete` marginals), `numbertheory.c` (prime sieve now heap-allocated per call), and `analytic_nt.c` (partition table + zeta Euler table now stack-local) all use stack-local or heap-per-call scratch. Remaining `static` instances are read-only tables (`static const`) or function-linkage helpers, which are thread-safe. Core TUs (trig/exp_log/complex/fft/linalg-solve) remain stateless. The DESIGN_CONTRACT "Stateless & Thread-Safe, No Global State" claim now holds for the full tree.

**17. Complex division** — Added `denom == 0` check in Smith's method to guard against catastrophic cancellation.

**18. Quadratics / cubic** — Discriminant always in `long double`; `ml_cubic` uses trig fallback for casus irreducibilis (3 real roots); saturation cast on `long double`→`double` root.

**19. Optimization / ODE** — Scale-aware convergence (`tol*(1+|x|)`), scale-aware finite-difference step (`sqrt(eps)*(1+|x|)`).

**20. Linear algebra** — `ml_matrix_exp_2x2` per-element overflow fixup: diagonal entries round to signed Inf; off-diagonal `em*b*sh` with `b==0` is exactly 0, not Inf/NaN, even when `em` is +Inf.

**21. SDE** — True Brownian bridge sampler (`ml_brownian_bridge_sample`) with variance `t*(T-t)/T`; OU exact transition density.

**22. SIMD dispatch runtime guard** — `ml_cpu_has_fma/avx2/sse41` use `__builtin_cpu_supports` at runtime; compile-time macro `__FMA__` reflects build flags, not host CPU — old macro test risked SIGILL.

**23. Benchmark rdtsc carve-out** — `benchmarks/bench.c` isolates non-portable `rdtsc` under x86 guard with `clock()` fallback, `BENCH_UNITS` cycles/ticks, excluded from lib build.

### New Batch-1 Modules (v12R2 — no oracle coverage, see `docs/PRECISION_CONTRACT.md`)

| Module | Header | Key Functions |
|--------|--------|---------------|
| **optim_n** | `ml_optim_n.h` | `ml_nelder_mead`, `ml_lbfgs_min`, `ml_adam_min` (n≤16/32) |
| **ode_sys** | `ml_ode_sys.h` | `ml_ode_dp5_sys`, `ml_ode_be2_sys`, `ml_ode_symplectic_verlet` (n≤16) |
| **spectral** | `ml_spectral.h` | `ml_cg_solve`, `ml_gmres_solve`, `ml_power_iter`, `ml_svd_jacobi`, `ml_qr_iter_eig` |
| **stats_inv** | `ml_stats_inv.h` | `ml_gamma_inv`, `ml_beta_inv`, `ml_chi2_inv`, `ml_student_t_inv`, `ml_f_inv`, `ml_normal_logcdf` |
| **sde** | `ml_sde.h` | `ml_sde_euler_maruyama`, `ml_sde_milstein`, `ml_brownian_bridge_sample`, `ml_ou_exact` |
| **pde** | `ml_pde.h` | `ml_heat_explicit`, `ml_heat_implicit`, `ml_wave_leapfrog`, `ml_poisson_1d`, `ml_fem1d_assemble` |
| **harmonic** | `ml_harmonic.h` | `ml_haar_fwt/iwt`, `ml_morlet_cwt`, `ml_fft_real`, `ml_fft2d_pow2` |
| **mcmc** | `ml_mcmc.h` | `ml_mh_sample_ctx`, `ml_kde_gaussian`, `ml_ess`, `ml_gelman_rubin` |
| **manifold** | `ml_manifold.h` | `ml_proj_sphere`, `ml_proj_stiefel`, `ml_exp_sphere`, `ml_sphere_dist` |
| **info** | `ml_info.h` | `ml_entropy`, `ml_kl_div`, `ml_cross_entropy`, `ml_mi_discrete`, `ml_logistic`, `ml_softplus` |
| **analytic_nt** | `ml_analytic_nt.h` | `ml_hurwitz_zeta` (stub for s≤1,a≠1), `ml_dirichlet_eta_cplx`, `ml_theta3`, `ml_partition_p`, `ml_zeta_cplx` |
| **control** | `ml_control.h` | `ml_lqr_gain_2x2`, `ml_kalman_1d`, `ml_lyapunov_2x2_trace` |
| **kelvin** | `ml_kelvin.h` | `ml_kelvin_ber/bei/ker/kei` (series + asym XMAX=20; ber/bei(Inf)=NaN) |

Despot audit 2026-09-30 (see `docs/DESPOT_AUDIT.md`): exp2 subnormal, remainder
long-double quotient, crt2 `%q`, kronecker INT64_MIN, mobius sentinel,
mult_order exact, cross_entropy +Inf, MH rejection, SDE `_ctx` + forward-only,
BE2 fail-loud, CG/GMRES relative tol, SVD/Jacobi scale-invariant stops,
kelvin Inf, fixed-point rounding, DST-I alias, Haar/manifold contracts,
EMBEDDED rsqrt, new status codes, `extern "C"`, 34-TU coherence.

ULP push 2026-09-30 (see `docs/ULP_PUSH.md`): pow 41→0 (LD log/exp + LD
integer accumulation), log 1→0 (LD split single source), sin benign 1→0 (LD
reduction + sinl/cosl), lgamma 6.7 4→0 (shift-to-8 exact LD xs), gamma 6.7
21→0 (LD direct), gamma 0.1 3→0 (Taylor radius ζ to 25), gamma 0.001 5→0 (LD
recurrence), gamma -0.5 1→0 (-2√π closed form), gamma -0.1/-0.9 1→0 (LD
reflection/recurrence). Oracle worst 5→0 ULP on the grid (measured, not a
proof for all inputs — Table Maker's Dilemma, x86-64 only).

---

## Test Results

- **Oracle:** 412 passed, 0 failed, worst 0 ULP (26 families: 7 core + 19 R6 — correctly rounded on the grid vs mpmath 50-dps; gate ≤5 ULP)
- **Wide grids (exact-binary truth):** pow 100 cases 0 ULP (was 41); sin/cos benign 25 cases 0 ULP (was 1); gamma 11-point grid 0 ULP (6.7 was 21, 0.1 was 3, 0.001 was 5, -0.5 was 1); vs-sys sin/cos/exp/log/pow 0, gamma/lgamma 1 (sys itself differs)
- **Edge tests:** 25 suites, all passed (750+ assertions; incl. `test_edge_despot_r4` 21 pins + `test_edge_despot_r5` 42 pins)
- **Fuzzers:** fuzz_god_mode (0 failed `--seed 123456789`: 61473/0 AVX2 host, 61393/0 elsewhere — host-sensitive), fuzz_boundary — all passed
- **Boundary gauntlet:** 25 passed, 0 failed
- **Soak:** 10,000 iterations available via `--soak`
- **Sanitizers:** ASan + UBSan clean (see gate log)
- **Profiles:** SCIENTIFIC / GRAPHICS / EMBEDDED all build `-Werror` clean
- **Thread-safety:** verified — no per-call mutable static state in any TU

## Deferred to v12A2

- Minimax coefficient swap (coefficients generated and validated, dormant)
- Interval arithmetic / verified computing
- Automatic differentiation
- Special functions (Bessel, elliptic, hypergeometric)
- Adaptive quadrature with certified error estimates
- Oracle coverage for Batch-1 modules

## What's Not In Scope

- Universal correctly-rounded transcendental functions
- Identical output across every architecture
- Replacement of specialized high-precision libraries

These are design choices, not hidden defects.

---

*v11S shipped 2026-08-02. v12A1 A1 closure completed 2026-08-11. v12R2 refinement completed 2026-09-30 (public V1.2-RC2, same commit).*
*The gamma nightmare is over. The critical bugs are fixed. The thread-safety audit is complete. 12 Batch-1 modules + kelvin ship (34 TUs).*

## Round-2 despot (2026-10-01)

Re-audit in `/tmp/opencode/mathlib-work/` found and fixed 11 issues: zeta OOB (ASan abort → 212/0 clean), kelvin sign flip (ber -pi/8 split), prime_pi odd miss, acosh 458k ULP → 0, softplus 57k ULP → exact, normal_inv sigma stall → exact, hurwitz B2 10k ULP → 9e-10, catalan C35/C36 exact, log1p coeff, Airy LD Taylor, QR NaN-poison + hurwitz alias, 35→34 TU docs. Oracle stays 212/0 worst 0 ULP (gate ≤5 ULP); god 61393/0 measured; edge 23/23; 3-profile clean. Remaining holes (kelvin 0.5%, Airy/K Temme) documented in KNOWN_LIMITATIONS/DESPOT_AUDIT (both closed in Round-3).

## Round-3 despot (2026-10-01) — closing the remaining limits

Every limit still open after Round-2 is closed with a real method change:

- **Kelvin**: 2-term P/Q replaced by the full DLMF 10.67.3-4 sums over
  `a_k/x^k` with per-function phase rotations (ber/bei `-pi/8 + 3k pi/4`,
  ker/kei `+pi/8 + k pi/4`), least-term truncation. `ker/kei` crossover
  split to 12. ~0.5% -> ~1e-13.
- **Bessel K0/K1**: the 8<x<10 transition hole is closed by an LD
  adaptive-Simpson evaluation of `K_nu(x) = int exp(-x cosh t) cosh(nu t) dt`
  (DLMF 10.32.10) for 4<=x<16. Series below, least-term asymptotic above.
  The integral is positive and non-cancelling, so this needs no Temme
  uniform expansion. ~2e-9 -> 6.8e-15 worst over [1e-3, 300].
- **Airy Ai**: same bridge via `Ai(x) = sqrt(x/3)/pi * K_1/3(zeta)` for
  2.5<x<8.5. ~5e-9 -> 2.1e-13 worst over a 0.1-step sweep of [-14, 200].
- **digamma**: long-double recurrence+Stirling plus an exact half-integer
  shortcut and LD `cot(pi x)`. 10 ULP -> 0 ULP.
- **Y0/Y1**: crossover 14 -> 13 (`Y1(13.9)` 2.56e-12 -> 1.08e-13). What
  remains is absolute error 1.1e-13 at zeros of Y1, which is intrinsic.
- **`long double == double`**: K and Ai measured in plain double at
  2.6e-16 / 2.4e-15 — no degradation. They are now more portable than the
  series/asymptotic hybrids they replaced.
- **Stack**: `ml_fft2d_pow2` 256 KB -> 2 KB; tree max frame now 65,600 B.
- **Guard**: `accuracy_audit` 371 assertions; K/Airy tolerances 5e7/1e8 ->
  65536 ULP with 10 new transition-band pins.

Gate: oracle 212/0 worst 0 ULP, modular 4/4, smoke 30/0, boundary 25/0,
god 61393/0 measured, edge 23/23, 3-profile `-Werror` clean.

## Round-4 despot (2026-10-07) — full system audit + fixes

Mathematical, programming, and operational audit of the entire tree against
external truth (mpmath 80-dps, system libm, DLMF identities, closed-form
physics). Full report: `docs/DESPOT_R4_AUDIT.md`.

- **Math (wrong-answer fix):** `ml_zeta(0)` 0 → -1/2 exact (trivial-zero guard
  wrongly included s==0; reflection needs zeta(1)=Inf so the closed form is
  returned directly). Verified vs DLMF 25.2.2; trivial zeros -2/-4 still 0.
- **Math (continuation closed):** `ml_hurwitz_zeta` 0<s<=1,a!=1 NaN stub →
  shared LD Euler-Maclaurin helper (N=200, B2..B12, Kahan; 0 ULP vs mpmath at
  1.1/1.5/0.9/0.01). Comment fix B6→B12; duplication 2x30 lines → 1 helper.
  s<=0,a!=1 stays fail-loud NaN (documented).
- **Programming:** complex pow shared `ml_cplx_pow_w` helper (was 3x LD
  product) + zero-guard in `ml_cplx_logarithm_ld`; sqrt comment fix
  (well-conditioned sum, not "exact power-of-two"); Stirling-LD TwoSum
  uniformity; `-Wundef` + `float.h` coherence (`LDBL_MANT_DIG` vs
  `__LDBL_MANT_DIG__` documented); `run_edge_tests.sh` honors
  `MATHLIB_EDGE_BUILD_DIR` so `/tmp` builds never touch the repo tree.
- **Operational:** root `*.o`/`*.d` clutter removed (gitignored but violated
  the `/tmp`-only evidence policy); TU lists re-verified 34/34 in all 4
  canonical files; fuzz figure corrected to measured 61393/0 (source
  unchanged since 02/08/26; 61473 was never reproducible).
- **Physics outcomes verified:** ODE RK4 exp(1), SUVAT, projectile range
  v0²/g, SHM cos, PDE FTCS r≤1/2 + CFL≤1 enforcement, Poisson/heat-implicit
  Thomas, Kelvin ber(20) ~1e-13 relative, K/Airy bridges hold.
- **Guard:** new `tests/test_edge_despot_r4.c` 21/21 pins every fix.
  Gate on the current tree: oracle 212/0 worst 0 ULP, modular 4/4, smoke 30/0,
  boundary 25/0, god 0 failed, accuracy_audit 371/0, R4 21/21, edge 24/24,
  3-profile `-Werror` + ASan/UBSan clean.

## Round-5 despot (2026-10-10, Edmonton) — full system audit + fixes

Mathematical, programming, and operational audit of the entire tree; full
report: `docs/DESPOT_R5_AUDIT.md`. R4 history above is left intact; R5
corrects the fuzz record (61473/0 reproduces on AVX2 hosts; 61393/0
elsewhere — `rand()` + conditional CHECKs + AVX make the count
host-sensitive; gate is 0 failed, not an exact count).

- **Math (no wrong-answer on valid inputs; one dead-clause fix):**
  `ml_hurwitz_zeta` dead integer limb removed + |v|>=2^53 guard; s<=0,a!=1
  stays fail-loud NaN (documented). Core + Batch-1 physics re-verified vs
  mpmath 80-dps / DLMF / closed form (see R5 §1).
- **Programming:** non-GNUC scalar fallback + honest `ml_cpu_has_avx2`,
  fixed-point shift/INT64_MIN, LDBL-guard unification, cross3/transpose
  aliasing, iterative overflow-checked egcd, MSVC builtin mapping,
  CONVERGENCE/OVERFLOW precision, DoS caps, RNG warmup, ODE `>=` arrival.
- **Operational:** `ultimate_fuzzer` CMake target; gate 10 steps + smoke/god/
  python-runner + `MATHLIB_BUILD_DIR`; runners honor `MATHLIB_EDGE_BUILD_DIR`
  / sanitizers; soak `/tmp` logs; root tracked ELFs removed; scripts marked
  historical; TU lists 34/34.
- **Guard:** new `tests/test_edge_despot_r5.c` 42/42 pins every fix.
  Gate on the current tree: oracle 212/0 worst 0 ULP, modular 4/4, smoke 30/0,
  boundary 25/0, god 61473/0 AVX2 (0 failed), accuracy_audit 371/0, R4 21/21,
  R5 42/42, edge 25/25, 3-profile `-Werror` + ASan/UBSan clean.

## Contest extravaganza (35th TU: `ml_contest.h` + `src/contest.c`)

36 functions tiered low→high, exactness contracts per house style
(`UINT64_MAX` sentinel, fail-closed verdicts, exactness-guarded doubles;
no oracle ULP claimed). See `docs/CONTEST_GUIDE.md`.

- **Tier 1:** triangular/pentagonal/hexagonal + is-triangular (integer
  Newton), double factorial, Lucas-mod-2 parity, stars-and-bars, Hanoi,
  Josephus, Lucas/Pell/Tribonacci.
- **Tier 2:** Bell (Aitken), Stirling 1st/2nd, Eulerian, Lah, Motzkin,
  large Schröder, Narayana (row sums = Catalan, pinned), distinct
  partitions, Goldbach counts, Sophie Germain, Lucas-Lehmer, Tonelli-Shanks
  (least root), tau/sigma, Frobenius number + count.
- **Tier 3:** harmonic via digamma, Basel partial, Faulhaber via Bernoulli
  `B0..B30` (sign convention pinned by S_1/S_4), Wythoff cold, subtraction
  Grundy, Schur-r1 verifier, AM-GM check, cubic Vieta sums.
- **Guard:** new `tests/test_edge_contest.c` 96/96 (Python big-int +
  mpmath truth; caught a Motzkin-index slip and the Bernoulli `B1` sign
  live). TU lists 35/35 in all 4 canonical files; edge now 26/26.

## Tier coverage: geometry, olympiad NT, graphs, modeling (38 TUs)

- **Geometry (`ml_geometry.h`):** triangle centers via Euler-vector
  orthocenter, `R=abc/4K`/`r=K/s` reusing staged Heron, Euler-line +
  nine-point verifiers, LD predicates, monotone-chain hull, clock, solids.
  Reuses (not duplicates) the hardened `ml_point_plane_dist`.
- **Olympiad NT (`ml_numbertheory.h` annex):** Legendre `v_p(n!)`,
  floor-sum, LTE (odd + 2-adic), CF period/convergents, Pell fundamental
  (`D=61` → 1766319049/226153980); recreational Tier-0 in `ml_contest.h`
  (digits, Collatz, magic squares).
- **Graphs (`ml_graph.h`):** DSU, BFS/DFS, bipartite, components, Kahn DAG,
  O(n²) Dijkstra, Kruskal (SINGULAR when disconnected); edge lists,
  `n<=4096`, nonneg weights.
- **Modeling (`ml_modeling.h`):** Markov stationary (CONVERGENCE-loud;
  imprimitive Leslie pinned), SIR R0/final-size/peak, Leslie step+growth,
  M/M/1 closed forms, Beta-Binomial/Gamma-Poisson + beta credible,
  logistic Lyapunov (ln2 at r=4, -ln2 at r=2.5).
- **Guard:** `test_edge_geometry` 45/45, `test_edge_olympiad_nt` 39/39,
  `test_edge_graph` 29/29, `test_edge_modeling` 40/40. TU lists 38/38;
  edge now 30/30.

## Doctorate slice (43 TUs): measure, funcanal, topology, diffgeom, algebraic

Computable graduate mathematics only (solvers, verifiers, exact integer
algorithms — never proofs). See `docs/DOCTORATE_GUIDE.md`.

- **Measure (`ml_measure.h`):** Dirichlet showcase (exactly 0), simple
  sums, discrete L^p, Hölder/Minkowski/Markov verifiers (fail closed).
- **Funcanal (`ml_funcanal.h`):** Banach solver (certified L + posteriori
  bound), Neumann inverses, Hilbert projection (scale-aware rank test —
  irrational norms never hit exact zero), Rayleigh max via staged power
  iteration.
- **Topology (`ml_topology.h`):** mod-2 Betti from validated closed
  complexes (heap bitsets, never 128KB stack), chain checks + Betti,
  Euler characteristic, exact Bareiss rank/det with 2^64 guards.
- **Diffgeom (`ml_diffgeom.h`):** graph/parametric curvature, Frenet
  (per-order step sizes — h=2e-8 loses ~60% on 2nd derivatives, fixed),
  sphere K=1/9, straight motion SINGULAR.
- **Algebraic (`ml_algebraic.h`) + PNT:** Sylvester resultants (sympy:
  Res(x²-2,x³-3)=1), discriminants, Sturm isolation, Descartes (exact),
  hook/S_n sizes (exact), Rodrigues; Ei/li/psi with li(1e6)≈78627.55.
- **Guard:** `test_edge_measure` 16/16, `test_edge_funcanal` 16/16,
  `test_edge_topology` 21/21, `test_edge_diffgeom` 18/18,
  `test_edge_algebraic` 32/32. TU lists 43/43; edge now 35/35.

## Round-6 ULP push: 26 families at measured 0 ULP (see `docs/ULP_PUSH.md`)

Starting from the easiest: log10 LO-sign fix (2→0), erf 1-erfc route
(6→0/1), log2 hex split, cbrt LD-neighbor Ziv-lite, log1p LD Horner,
atanh/asinh/acosh LD forms, sinh/cosh/tanh/expm1 LD composition,
asin/acos LD + atanl, ml_exp true-DD Cody-Waite reduction (167→3 fails
on [0,1], remainder Table-Maker), atan LD tree + compensated Horner
(200→23), erf LD Taylor (175→0), erfc LD complements, Bessel Y1
harmonic-LD (5→0.85 ULP benign). New invisible-composition law documented:
two 0-ULP doubles still compose to 1 ULP; only single-rounding from
>64-bit evaluation attains 0. Residual Table-Maker boundary catalog
shipped (exact doubles). Oracle grows 212 → 412 vectors over 26 families
(19 new, Bessel J benign-only by generator rule), worst 0 ULP.

## Round-6 fuzz upgrade + QR/fft fixes (see `docs/DESPOT_R6_AUDIT.md`)

- **Coverage:** god-mode Part 2 adds 51,600 fixed-count checks (xorshift
  stream, seed/host-stable total) over Batch-1/contest/geometry/graph/
  modeling/doctorate invariants and differentials (Narayana=Catalan,
  floor/LTE/Pell/Tonelli vs naive, DSU/BFS/Dijkstra/Kruskal agreement,
  Markov/SIR/Little/conjugate, R>=2r, Betti alternating sum, Sturm vs
  Descartes, Rodrigues, roundtrips, integer blizzard). Part-1 stream
  untouched. Totals: god 113,073/0, boundary 63/0 (was 25), ultimate
  89,462/0 (was 49,962; Inf-aware crash predicates).
- **Fixed:** `ml_qr_iter_eig` nested-similarity stall (~20% on SPD;
  Francis Hessenberg + explicit double-shift QR with reorthogonalized
  MGS and single-shift exceptional path; 0 non-convergences in 2000
  general 2x2..6x6, LAPACK-agreeing); `ml_fft_is_supported` now
  enforces the 2^24 cap it documents.
- **Gate:** oracle 412/0 worst 0 ULP, modular 4/4, smoke 30/0, boundary
  63/0, god 113073/0, ultimate 89462/0, edge 35/35, 3-profile `-Werror`
  + ASan/UBSan clean.
