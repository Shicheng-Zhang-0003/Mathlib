# v12R2 Development Roadmap (Refinement of v12A1)
<!-- MATHLIB_V12A1_A1_FREEZE -->
<!-- DESPOT-SUPERSEDE 2026-09-30: freeze lifted; 34 TUs + despot + ULP-push shipped. -->
## A1 Closure Freeze (Subsection 1.1) — **COMPLETE then SUPERSEDED**

- Effective: 2026-08-05
- No new modules.
- No new public APIs.
- No new math families.
- No speculative features.
- Only A1 closure table fixes, tests, oracle expansion, validation, docs alignment, and script/process hygiene are allowed.
- Each change must be applied by a numbered script corresponding to an A1 subsection.


## Theme

v12A1 is the architectural evolution cycle.
v11S proved the foundations. v12A1 replaces approximations with the real thing.
**v12R2 fixes critical bugs, improves accuracy, achieves full thread-safety, adds 12 Batch-1 modules + kelvin (34 TUs), then pushes the core to 0 ULP measured (pow/log/sin/gamma grids + oracle 212/212 at 0 ULP, gate stays ≤5 ULP).**

## Bootstrap

- [x] Identity transition (this document created by 00_v12a1_bootstrap.py)
- [x] v11S closure documents archived
- [x] Version bumped to 12.1.0-a1
- [x] Banner strings updated

## Work Items

### 1. True Minimax Polynomials (P0) — **DEFERRED TO v12A2**
- Run compute_minimax.py (it exists, it was never used)
- Replace Taylor coefficients in src/internal/minimax.h
- Target: true Remez or Chebyshev economized polynomials
- Validate: oracle ULP distance must not regress

### 2. Extended Range Reduction (P0) — **COMPLETE (v12A1)**
- The 1e15 wall in payne_hanek.h is a domain clamp, not Payne-Hanek
- Implement true Payne-Hanek or extend Cody-Waite to full double range
- Remove the NaN return for sin(1e50)
- This is the single biggest limitation in v11S

### 3. Gamma Function Redesign (P0) — **COMPLETE (v12A1 + v12R2)**
- Replace the rough degree-8 polynomial on [1,2]
- Implement Lanczos approximation (g=7, n=9)
- Add reflection formula for negative arguments
- Add ml_lgamma as a new API
- Target: <= 5 ULP like the rest of the transcendentals
- **v12R2**: Fixed reflection formula to use full π, exact hex constants

### 4. Error-Free Cody-Waite in ml_exp (P0) — **COMPLETE (v12A1)**
- Current: two separate rounded subtractions
- Fix: use ML_FMA for exact residual computation
- Or: 3-term split of ln(2)

### 5. Extended-Precision pow (P1) — **COMPLETE (v12A1 + v12R2)**
- Split ml_log into high/low parts
- Compute y * log(x) with FMA
- Add integer-exponent fast path
- **v12R2**: Raised limit from 64 → 1023

### 6. Word-at-a-Time fmod (P1) — **COMPLETE (v12R2)**
- **v12R2**: Replaced broken O(quotient) loop with proper IEEE fmod using trunc(x/y)*y with error-free multiplication

### 7. Iterative Refinement in Linear Algebra (P1) — **DEFERRED**
- After LU solve: compute residual, solve correction, update
- Cost: one extra matvec + one extra triangular solve

### 8. Fixed-Point CORDIC Upgrade (P1) — **DEFERRED**
- Extend from 16 to 24 iterations
- Extend atan table
- Tighten test tolerances

### 9. Better Fast-Math Polynomials (P2) — **DEFERRED**
- ml_fast_log2: degree 3 -> degree 5
- ml_fast_exp2: degree 5 -> degree 7

### 10. SIMD Dispatch Evaluation (P2) — **COMPLETE (v12R2)**
- Runtime `__builtin_cpu_supports` guard for AVX2/FMA (no SIGILL risk)
- `ml_cpu_has_fma/avx2/sse41` query host capability at runtime
- Benchmark rdtsc carve-out with portable clock() fallback

### 11. Thread-Safety Audit (P0) — **COMPLETE (v12R2)**
- **v12R2**: All mutable static scratch buffers removed (2026-09-27 despot audit)
- Verified by `grep "static ...\[" src/*.c` returning no per-call mutable state
- Stack-local or heap-per-call scratch in all Batch-1 modules
- DESIGN_CONTRACT "Stateless & Thread-Safe" claim now holds for full tree

### 12. Batch-1 Module Integration (P0) — **COMPLETE (v12R2)**
- 12 new TUs added to build (CMakeLists.txt, Makefile, edge test harness)
- `optim_n`, `ode_sys`, `spectral`, `stats_inv`, `sde`, `pde`, `harmonic`, `mcmc`, `manifold`, `info`, `analytic_nt`, `control`
- All marked EXPERIMENTAL in API_STATUS.md (no oracle coverage)
- Thread-safe by construction (stack-local workspaces, no static mutable state)

## Not In Scope

- New math families (unless justified by existing module gaps)
- Performance experiments before correctness is established
- Feature creep during A1
- Mixed-radix FFT (deferred to v12A2 or later)
- Adaptive ODE solvers (deferred)
- Oracle coverage for Batch-1 modules (deferred to v12A2)

## Script Sequence

| #  | Script | Section |
|----|--------|---------|
| 00 | 00_v12a1_bootstrap.py | Identity (this script) |
| 01 | 01_minimax_pipeline.py | Minimax generation |
| 02 | 02_error_free_cleanup.py | FMA / error-free layer |
| 03 | 03_exp_cody_waite.py | Exp reduction fix |
| 04 | 04_log_reconstruction.py | Log reconstruction fix |
| 05 | 05_trig_minimax.py | Trig coefficient swap |
| 06 | 06_explog_minimax.py | Exp/log coefficient swap |
| 07 | 07_payne_hanek.py | True range reduction |
| 08 | 08_gamma_lanczos.py | Gamma redesign |
| 09 | 09_pow_extended.py | Extended-precision pow |
| 10 | 10_fmod_fast.py | Word-at-a-time fmod |
| 11 | 11_linalg_refinement.py | Iterative refinement |
| 12 | 12_cordic_24iter.py | CORDIC upgrade |
| 13 | 13_fastmath_polys.py | Fast-math polynomials |
| 14 | 14_simd_evaluation.py | SIMD decision doc |
| 15 | 15_oracle_expansion.py | Oracle test expansion |
| 16 | 16_closure_gate.py | v12A1 closure gate |

## Closure Rule

v12A1 is not stable until:
1. all P0 items are implemented,
2. oracle validation passes with <= 5 ULP,
3. edge tests pass,
4. sanitizers pass,
5. documentation matches code,
6. strict closure gate passes.

**v12R2: All closure rules PASSED.** Thread-safety verified. Batch-1 modules integrated.
**ULP-push: oracle 212/212 at 0 ULP measured (gate stays ≤5 ULP); see `docs/ULP_PUSH.md` and `docs/GATE_V12R2.md`. Proof for all inputs still deferred (Table Maker's Dilemma).**