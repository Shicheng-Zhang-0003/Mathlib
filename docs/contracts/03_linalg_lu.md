# Math Diary: Linear Algebra (Zero-Alloc LU Decomposition)

> **STATUS: OWNED (2026-09-27)**
>
> History: v12R2 DORMANT template (unverified blanks) retired; v11S closure
> retained — LU solve uses partial pivoting, relative singularity
> thresholding, workspace alignment/canary hardening, and invalid-argument
> guards. Content below reverse-engineered from `src/linalg.c`,
> `include/mathlib/ml_linalg.h`, `include/mathlib/ml_tensor.h` and tiered
> against `docs/PRECISION_CONTRACT.md` (smoke/edge tiers; LU has no oracle
> coverage).

---

## 1. The Mathematical Identity
*What fundamental mathematical truth does this module enforce?*

- **Core factorization (Doolittle PA=LU):** `ml_lu_decomp` (`src/linalg.c:8`,
  declared `include/mathlib/ml_linalg.h:9`) computes `P*A = L*U` with row
  partial pivoting. `P` is a permutation vector (`P[i] = i` at entry, pairwise
  row swaps thereafter). `L` is unit-lower-triangular (implicit ones on the
  diagonal, multipliers stored strictly below) and `U` is upper-triangular
  (diagonal plus above); both share one `n x n` `LU` buffer. The defining
  invariant is `A[P[i],j] = sum_k L[i,k]*U[k,j]` up to rounding.
- **Cholesky variant:** `ml_cholesky` (`src/linalg.c:303`) enforces
  `A = L*L^T` for symmetric positive-definite input. Symmetry is checked
  first (`diff <= 1e-12*sc`, asymmetric pair returns `ML_ERR_INVALID_ARG`,
  non-finite returns `ML_ERR_NAN_INPUT`); then `s = A[i,j] - sum_k
  L[i,k]*L[j,k]`, diagonal `L[i,i] = sqrt(s)` with `s > 0` else
  `ML_ERR_SINGULAR`, off-diagonal `L[i,j] = s / L[j,j]`; upper triangle is
  zeroed.
- **QR variant:** `ml_qr_solve` (`src/linalg.c:345`) factors tall-or-square
  `A = Q*R` by Modified Gram-Schmidt. `Q` (orthonormal columns) lives in
  workspace `V`, `R` is upper-triangular; invariant `Q^T*Q = I`,
  `A = Q*R` for full-rank `m >= n`.
- **Determinant identity:** `ml_determinant` (`src/linalg.c:487`) returns
  `det(A) = sign(P) * prod_i U[i,i]`. `sign(P)` comes from inversion-count
  parity of `P` (`inv` pairs `P[i] > P[j]`, `i < j`; odd negates), not from a
  swap counter. Singular factor returns `0.0`; other failures return NaN.
- **2x2 spectral identities:** `ml_eigen2x2` (`src/linalg.c:681`) solves
  `lambda^2 - tr*lambda + det = 0` via the stable q-form: `q = 0.5*(tr +
  sgn(tr)*sqrt(disc))`, large root `l0 = q`, small root `l1 = det/q` in
  `long double`, so no cancellation when `|tr| >> sqrt|det|`. `ml_svd_2x2`
  (`src/linalg.c:837`) returns `s = sqrt(eig(A^T*A))` with `e = a^2+c^2`,
  `f = a*b+c*d`, `g = b^2+d^2`, same stable discriminant path, ordered
  `s0 >= s1`.
- **2x2 exponential:** `ml_matrix_exp_2x2` (`src/linalg.c:869`) is closed-form
  Cayley-Hamilton / Jordan calculus. With `m = (a+d)/2`,
  `ah = (a-d)/2`, `D = ah^2 + b*c`, `em = exp(m)`: for `D >= 0`,
  `exp(A) = em*(cosh(s)*I + sinh(s)/s*(A - m*I))`, `s = sqrt(D)`; for
  `D < 0`, `cosh/sinh` become `cos/sin` with `s = sqrt(-D)`. `s == 0`
  short-circuits `sinh(s)/s` (resp. `sin(s)/s`) to `1`. **v12R2 per-element
  overflow fixup:** diagonal entries round to signed Inf; off-diagonal
  `em*b*sh` with `b==0` (resp. `c==0`) is exactly 0, not Inf/NaN, even when
  `em` is +Inf. Overflow rounded to signed Inf is representable; only NaN
  remains a failure.

---

## 2. The Algorithm (Execution Path)
*Trace the data flow from the raw `double` input to the final `double` output. What are the distinct phases?*

1. **Phase 1: Input Validation & Workspace Setup** (`ml_lu_decomp:8-45`)
   - NULL checks, dimension checks (`A.cols == n`, `LU.rows == n`, `LU.cols == n`)
   - `size_t` overflow guard: rejects `n` where `n*n` would wrap on 32-bit
   - Copies `A` → `LU` with immediate non-finite rejection (`ML_ERR_NAN_INPUT`)
   - Initializes permutation `P[i] = i`

2. **Phase 2: Infinity Norm for Relative Singularity Threshold** (`ml_lu_decomp:51-77`)
   - Scaled row-sum: `rmax = max_j |A[i,j]|`, then `sum_j |A[i,j]|/rmax * rmax`
   - Prevents overflow for `1e308` entries (valid finite matrix → spurious NaN)
   - `matrix_norm == 0` → exact zero matrix → `ML_ERR_SINGULAR`
   - `matrix_norm` infinite → scale too large → `ML_ERR_NAN_INPUT`

3. **Phase 3: Relative Singularity Threshold** (`ml_lu_decomp:89-100`)
   - `threshold = matrix_norm * 2.220446049250313e-16 * n` (machine eps * n)
   - If threshold underflows to 0, use smallest subnormal `4.94e-324`
   - Preserves scale sensitivity for extremely small matrices

4. **Phase 4: LU Decomposition with Partial Pivoting** (`ml_lu_decomp:102-165`)
   - For each column `i`: find `max_row` with max `|LU[k,i]|` for `k >= i`
   - Swap rows `i` and `max_row` in `LU`, update `P`, track parity
   - Pivot `|LU[i,i]| <= threshold` → `ML_ERR_SINGULAR`
   - Compute multipliers `LU[k,i] /= LU[i,i]` for `k > i`
   - Rank-1 update: `LU[k,j] -= LU[k,i] * LU[i,j]` for `k,j > i`

5. **Phase 5: Forward/Back Substitution** (`ml_solve:177-215`)
   - Permute RHS: `b_perm[i] = b[P[i]]`
   - Forward: `y[i] = b_perm[i] - sum_{j<i} LU[i,j]*y[j]`
   - Backward: `x[i] = (y[i] - sum_{j>i} LU[i,j]*x[j]) / LU[i,i]`

---

## 3. The "Dark Arts" (Hardware & IEEE-754 Sympathy)
*Why did the AI write it this way instead of the naive textbook way?*

- **Bitwise / Memory Tricks:**
  - `memcpy` for bitwise IEEE ops (`ml_fabs`, `ml_isnan`, `ml_copysign`) — avoids strict-aliasing UB from unions
  - `ML_TENSOR_AT(A,i,j)` macro with `size_t` indexing — `size_t` overflow guard before any access
  - Workspace `ml_workspace_t` is a bump allocator with canary (`ws->base + ws->size == ws->canary`) — detects buffer overrun without `malloc`

- **Branchless Logic:**
  - Pivot search uses `ml_fabs` + conditional move pattern (compiler lowers to `cmov` on x86)
  - Scaled norm accumulation avoids `if (v > rmax)` in inner loop via `rmax` pre-scan

- **Constants:**
  - `2.220446049250313e-16` = `DBL_EPSILON` (machine epsilon, 2⁻⁵²)
  - `4.9406564584124654e-324` = `DBL_TRUE_MIN` (smallest positive subnormal, 2⁻¹⁰⁷⁴)
  - `ML_LN2_HI/LO` split for Cody-Waite in exp/log (not directly in LU but shared)

---

## 4. Error Budget & Failure Modes
*Where does the physical limit of the silicon break the mathematical ideal?*

- **Maximum ULP Error:** LU has no oracle coverage (smoke/edge tier only). The decomposition is backward stable: `||A - P^T*L*U|| <= O(eps) * ||A||`. No ULP claim may be cited.
- **Catastrophic Cancellation:**
  - Rank-1 update `LU[k,j] -= LU[k,i] * LU[i,j]` can lose precision for ill-conditioned matrices
  - Relative threshold `matrix_norm * eps * n` prevents false singular on scaled matrices
  - Subnormal fallback threshold preserves scale sensitivity
- **Overflow/Underflow Boundaries:**
  - Scaled norm computation prevents `1e308` row-sum overflow
  - `ml_determinant` returns `0.0` for singular, NaN for other failures
  - `ml_matrix_exp_2x2` per-element overflow fixup (v12R2): diagonal→signed Inf, off-diagonal with zero coefficient→exact 0, only NaN is error
- **NaN Propagation:**
  - Non-finite input → `ML_ERR_NAN_INPUT` immediately (before any compute)
  - Singular pivot → `ML_ERR_SINGULAR`
  - Workspace exhaustion → `ML_ERR_WORKSPACE` (canary check)
  - All structural APIs return `ml_status_t`, never silent NaN

---

## 5. Orchestrator's Reverse-Engineering Log
*Paste insights, proofs, and "Aha!" moments from your Red/Green team consultations here.*

- 2026-09-27 - Despot audit thread-safety: Jacobi eigen `W[64][64]` was `static` in `linalg.c:721` (mutable per-call). Converted to stack-local `double W[64][64]` inside function. Verified `grep "static double\[" src/linalg.c` returns zero hits.
- 2026-09-27 - Despot audit `ml_matrix_exp_2x2` overflow: original code returned NaN on any overflow. Fix: per-element fixup distinguishes `0*Inf=0` (off-diagonal with b=0 or c=0) from true Inf (diagonal) and NaN (genuine failure). Matches mathematical limit.
- 2026-09-27 - Relative threshold rationale: absolute `1e-15` fails for tiny matrices (e.g., `1e-20 * I` has norm `1e-20`, threshold `1e-15` → false singular). Relative `norm * eps * n` scales with matrix. Subnormal fallback `DBL_TRUE_MIN` handles norm < `1e-308`.
- 2026-09-27 - `ml_solve_refined` (iterative refinement) stubbed in header but not implemented — deferred to v12A2 per roadmap item 7.
- 2026-09-27 - QR solve uses MGS (not Householder) for workspace simplicity: `Q` in caller workspace, no extra allocation. MGS loses orthogonality for ill-conditioned `A` but acceptable for smoke-tier.
