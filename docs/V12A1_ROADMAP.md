# v12A1 Development Roadmap
<!-- MATHLIB_V12A1_A1_FREEZE -->
## A1 Closure Freeze (Subsection 1.1)

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

## Bootstrap

- [x] Identity transition (this document created by 00_v12a1_bootstrap.py)
- [x] v11S closure documents archived
- [x] Version bumped to 12.1.0-a1
- [x] Banner strings updated

## Work Items

### 1. True Minimax Polynomials (P0)
- Run compute_minimax.py (it exists, it was never used)
- Replace Taylor coefficients in src/internal/minimax.h
- Target: true Remez or Chebyshev economized polynomials
- Validate: oracle ULP distance must not regress

### 2. Extended Range Reduction (P0)
- The 1e15 wall in payne_hanek.h is a domain clamp, not Payne-Hanek
- Implement true Payne-Hanek or extend Cody-Waite to full double range
- Remove the NaN return for sin(1e50)
- This is the single biggest limitation in v11S

### 3. Gamma Function Redesign (P0)
- Replace the rough degree-8 polynomial on [1,2]
- Implement Lanczos approximation (g=7, n=9)
- Add reflection formula for negative arguments
- Add ml_lgamma as a new API
- Target: <= 5 ULP like the rest of the transcendentals

### 4. Error-Free Cody-Waite in ml_exp (P0)
- Current: two separate rounded subtractions
- Fix: use ML_FMA for exact residual computation
- Or: 3-term split of ln(2)

### 5. Extended-Precision pow (P1)
- Split ml_log into high/low parts
- Compute y * log(x) with FMA
- Add integer-exponent fast path
- Add near-integer result detection

### 6. Word-at-a-Time fmod (P1)
- Current: O(quotient) loop, up to 2046 iterations
- Fix: process in 64-bit chunks

### 7. Iterative Refinement in Linear Algebra (P1)
- After LU solve: compute residual, solve correction, update
- Cost: one extra matvec + one extra triangular solve

### 8. Fixed-Point CORDIC Upgrade (P1)
- Extend from 16 to 24 iterations
- Extend atan table
- Tighten test tolerances

### 9. Better Fast-Math Polynomials (P2)
- ml_fast_log2: degree 3 -> degree 5
- ml_fast_exp2: degree 5 -> degree 7

### 10. SIMD Dispatch Evaluation (P2)
- Decision document: is the Quake rsqrt hack worth keeping?
- No code change unless decision is to replace or remove

## Not In Scope

- New math families (unless justified by existing module gaps)
- Performance experiments before correctness is established
- Feature creep during A1
- Mixed-radix FFT (deferred to v12A2 or later)
- Adaptive ODE solvers (deferred)

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
