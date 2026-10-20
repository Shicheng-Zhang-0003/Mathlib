# MathLib v12R3 Development Tree

This unreleased working tree is version 12.3.0 and contains 59 C translation
units. The v12R2 release notes below remain historical. The v12R3 tree adds
experimental exact/numerical extensions and a pre-university utility module.

v12R2 followed the v12A1 architectural evolution. Its 43-translation-unit
history is summarized below.

<!-- MATHLIB_V12A1_DOCS_ALIGNMENT -->
## v12R2 refinement summary (historical)

v12R2 key fixes:
- **Removed x86 inline asm** — portable `__builtin_sqrt` with correct rounding
- **Fixed fmod** — proper IEEE-754 `x - trunc(x/y)*y` with error-free multiplication
- **Fixed atan/asin/acos** — argument reduction to |x| ≤ tan(π/12), stable formulas
- **Fixed tanh/atanh/asinh** — small-x threshold 1e-4 → 1.5e-8 (was millions of ULP error)
- **Fixed pow integer path** — limit 64 → 1023 (exact binary exponentiation)
- **Fixed gamma reflection** — uses full π (ML_PI_HI_D + ML_PI_LO_D)
- **Fixed variance** — Welford's online algorithm (no catastrophic cancellation)
- **ML_FMA software fallback** — proper Dekker Two-Product+Two-Sum
- **All decimal constants → exact hex floats** — ln2, π/√2, polynomial coefficients
- **Thread-safety audit** — all mutable static scratch buffers removed (2026-09-27)
- **ml_polynomial_eval / newton** — ML_FMA Horner; bogus dfx<epsilon check removed
- **Quadratics / cubic** — long-double discriminant, trig fallback for casus irreducibilis
- **Optimization / ODE** — scale-aware convergence and finite-difference steps
- **Linear algebra** — per-element overflow fixup in 2x2 exp; QR iterative eigensolver
- **SDE** — true Brownian bridge sampler; OU exact transition
- **SIMD dispatch** — runtime `__builtin_cpu_supports` guard (no SIGILL risk)

New Batch-1 modules (v12R2):
- `ml_optim_n` — Nelder-Mead, L-BFGS, Adam (n≤32, stack workspaces)
- `ml_ode_sys` — DP5, backward-Euler step-doubling, symplectic Verlet (n≤16)
- `ml_spectral` — CG, GMRES, power iteration, Jacobi SVD, QR eig (n≤256/64/16)
- `ml_stats_inv` — inverse CDF for gamma, beta, χ², Student-t, F, normal log-CDF
- `ml_sde` — Euler-Maruyama, Milstein, Brownian bridge (mean + sampler), OU exact
- `ml_pde` — explicit/implicit heat, wave leapfrog, Poisson 1D, FEM P1 element
- `ml_harmonic` — Haar FWT/IWT, Morlet CWT, real FFT, 2D FFT (power-of-two)
- `ml_mcmc` — Metropolis-Hastings (ctx), KDE, ESS, Gelman-Rubin
- `ml_manifold` — sphere/Stiefel projection, sphere exponential, sphere distance
- `ml_info` — entropy, KL, cross-entropy, discrete MI, logistic, softplus
- `ml_analytic_nt` — Hurwitz zeta (stub for s≤1,a≠1), Dirichlet eta, theta3, partition p, complex zeta
- `ml_control` — 2x2 LQR (Newton-Kleinman), 1D Kalman, Lyapunov trace margin
- `ml_kelvin` — ber/bei/ker/kei series + asymptotic (XMAX=20)
- `ml_contest` — contest extravaganza, 42 functions Tier-0→high (recreational, figurate, Bell/Stirling/Eulerian/Lah/Motzkin/Schroder/Narayana, partitions, Goldbach, Tonelli, Faulhaber, Wythoff, Grundy, Schur, Vieta; see `docs/CONTEST_GUIDE.md`)
- `ml_geometry` — triangle centers, R/r, Euler-line + nine-point, shoelace, seg-intersect, point-in-poly, convex hull, clock, solids
- `ml_school` — pre-university unit conversions, rounding modes, percentages, sequences, linear-coordinate helpers, and common measurement formulas
- `ml_graph` — DSU, BFS/DFS, bipartite, components, Kahn topo, Dijkstra, Kruskal
- `ml_modeling` — Markov chains, SIR, Leslie, M/M/1, conjugate Bayes, logistic Lyapunov
- `ml_measure` — Lebesgue Dirichlet, simple sums, Lp, Hölder/Minkowski/Markov verifiers
- `ml_funcanal` — Banach solver, Neumann inverses, Hilbert projection, Rayleigh max
- `ml_topology` — simplicial mod-2 Betti, Euler char, Bareiss rank/det
- `ml_diffgeom` — curvature, Frenet frames, Gaussian/mean curvatures
- `ml_algebraic` — resultants, discriminants, Sturm, hook/S_n sizes, Rodrigues (see `docs/DOCTORATE_GUIDE.md`)

Despot audit (2026-09-30, see `docs/DESPOT_AUDIT.md`): exp2 subnormal, remainder
long-double quotient, crt2 `%q` fix, kronecker INT64_MIN, mobius sentinel,
mult_order exact reduction, cross_entropy +Inf, MH rejection accounting, SDE
ctx variants + forward-only dt, BE2 Newton fail-loud, CG/GMRES relative tol,
SVD/Jacobi scale-invariant stops + convergence signals, kelvin Inf, fixed-point
rounding, DST-I alias, Haar/manifold contracts, EMBEDDED rsqrt, new status
codes, `extern "C"` throughout, 34-TU build coherence.

## Build

```bash
cmake -B build -DMATHLIB_PROFILE=SCIENTIFIC
cmake --build build
```

## Test

```bash
python3 run_all_tests.py
```

## Current status
<!-- MATHLIB_V12A1_README_STATUS_V2 -->
<!-- MATHLIB_V12A1_A1_FREEZE -->

The v12R3 code and audit findings are recorded in [`docs/DESPOT_R7_AUDIT.md`](docs/DESPOT_R7_AUDIT.md).
Subject coverage and the current implementation sequence are in
[`docs/MATH_COVERAGE.md`](docs/MATH_COVERAGE.md); exact-arithmetic and
certificate contracts are in [`docs/EXACT_FOUNDATIONS.md`](docs/EXACT_FOUNDATIONS.md).
The elementary-through-advanced-secondary scope is mapped in
[`docs/PREUNIVERSITY_COVERAGE.md`](docs/PREUNIVERSITY_COVERAGE.md).
The exact-foundations step adds dynamic integers/rationals with finite-decimal
parsing, outward-rounded basic intervals, exact Q[x] operations with a
division-certificate verifier, and prime-field arithmetic. The broader v12R3 extension modules have directed
regression coverage but no oracle-tier accuracy certification. Do not treat
the historical v12R2 accuracy measurements as validation of those additions.

See `docs/V12A1_ROADMAP.md` for the work plan.
See `release_notes.md` for the v12R2 refinement summary.
See `docs/API_STATUS.md` for current module status and `docs/BUILD_AND_INSTALL.md`
for installation details.
