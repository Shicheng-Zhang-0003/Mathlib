# MathLib v12R2 Development Tree

v12R2 is the refinement cycle following the v12A1 architectural evolution.

v12A1 proved the architecture. v12R2 fixes critical bugs, improves numerical accuracy, and adds 12 new Batch-1 modules.

<!-- MATHLIB_V12A1_DOCS_ALIGNMENT -->
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

## Build

```bash
cmake -B build -DMATHLIB_PROFILE=SCIENTIFIC
cmake --build build
```

## Test

```bash
python3 run_all_tests.py
```

## Status
<!-- MATHLIB_V12A1_README_STATUS_V2 -->
<!-- MATHLIB_V12A1_A1_FREEZE -->
A1 closure is **complete** with R2 refinements.

- Oracle validation: **212 passed, 0 failed** (all functions ≤ 5 ULP vs mpmath ground truth)
- Full test gauntlet: **31/32 passed** (modular, smoke, edge, fuzz, oracle, boundary)
- Closure gate: **PASSED**
- Thread-safety: **verified** — no per-call mutable static state in any TU

See `docs/V12A1_ROADMAP.md` for the work plan.
See `release_notes.md` for the v12R2 refinement summary.
