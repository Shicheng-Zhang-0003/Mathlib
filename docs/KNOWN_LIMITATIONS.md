# Known Limitations

<!-- MATHLIB_V12A1_DOCS_ALIGNMENT -->

v12A1 intentionally does not claim:

- universal correctly rounded transcendental functions
- identical output across every architecture
- replacement of specialized high precision libraries

These limitations are design choices, not hidden defects.

## v12A1 Specific

- **Minimax coefficients are DORMANT.** The generated coefficients in
  `src/internal/minimax_coeffs.h` are validated but not active.
  The running kernels use Maclaurin (Taylor) series, which already
  meet the <=5 ULP oracle gate. The minimax swap is deferred to v12A2.

- **The 1e15 wall is removed.** Payne-Hanek V6 handles the full
  double range up to ~1.8e308. `sin(1e50)` and `cos(1e300)` produce
  finite, correct results.

- **A1 closure freeze is lifted for the despot audit.** v12R2 shipped 12
  Batch-1 modules plus kelvin (34 TUs); the freeze document is historical.
  No further silent scope creep: CMake/Makefile/run_all_tests.py/
  run_edge_tests.sh canonical lists must match (34 TUs).

- **Gamma uses Lanczos g=7 n=9.** This coefficient set has ~1e-15
  intrinsic approximation error. Half-integers bypass Lanczos entirely
  via exact product/sum formulas. The 1-step recurrence for x < 0.5
  avoids the cancellation that plagued the old multi-step recurrence.

- **`ml_integral_traditional` remains experimental.** It is a simple
  Riemann sum integrator and is not part of the validated numerical
  core.

- **Thread-safety: all mutable static scratch buffers removed (2026-09-27
  despot audit).** Verified by `grep "static ...\[" src/*.c` returning no
  per-call mutable state: `optim_n.c`, `mcmc.c`, `pde.c`, `manifold.c`,
  `harmonic.c`, `spectral.c`, `calculus.c` (spline Thomas), `linalg.c`
  (Jacobi `W`), `info.c` (`ml_mi_discrete` marginals), `numbertheory.c`
  (prime sieve is now heap-allocated per call), and `analytic_nt.c`
  (partition table + zeta Euler table are now stack-local) all use
  stack-local or heap-per-call scratch. Remaining `static` instances are
  read-only tables (`static const`) or function-linkage helpers, which are
  thread-safe. Core TUs (trig/exp_log/complex/fft/linalg-solve) remain
  stateless. The DESIGN_CONTRACT "Stateless & Thread-Safe, No Global State"
  claim now holds for the full tree; re-verify with
  `grep -n "static double\|static float\|static int\|static long double\|static cplx\|static uint" src/*.c`
  (expect only `static` helper *functions*, which are stateless).

---

## v12R2 Accuracy Audit

An independent audit (mpmath 50-dps references, ULP-measured) repaired the
following. `tests/test_edge_accuracy_audit.c` pins the results as a
regression guard (282 assertions). Per `PRECISION_CONTRACT.md` this is a
regression guard, not an oracle-tier precision certification.

| Function | Before | After | Note |
|---|---|---|---|
| `ml_erfc` | ~1064 ULP | <=0.7 ULP | Laplace CF + convergent series, long double |
| `ml_bessel_j0/j1/y0/y1` | ~1e-3 rel | <=4e-13 rel | Hankel P/Q, least-term asymptotics |
| `ml_bessel_i0/i1` | exact | exact | |
| `ml_bessel_k0/k1` | ~1e9 rel | <=5e-9 rel | transition band x~8-10; needs Temme for more |
| `ml_airy_ai` | wrong coefficients | <=5e-9 rel | DLMF 9.7.5/9.7.6, recurrence c_k=c_{k-1}(6k-5)(6k-1)/(72k) |
| `ml_digamma` | ~70 ULP | <=2 ULP (pos) | Kahan recurrence + Stirling to x^-12 |
| `ml_jacobi_symbol` | wrong sign | exact | reciprocity test, negative a, n>2^63 |
| `ml_exp10` | ~1e-16 rel | <=1 ULP | corrected residual constant |
| `ml_cosh` | overflow at 710 | <=1 ULP | corrected ML_HALF_EXP_709 |

Known remaining limits (documented, not defects):

- `ml_bessel_k0/k1` in 8 < x < 10: ~1e-9 relative. The ascending series
  cancels like exp(x^2/4) and the asymptotic bottoms out at exp(-2x); the
  elementary method cannot do better. Temme's uniform expansion is required.
- `ml_airy_ai` in 5 < x < 6: ~5e-9 relative. Same reason: the Taylor series
  loses digits to cancellation while the asymptotic is still converging.
- `ml_digamma` for x < 0: ~10 ULP at x=-0.5 via the reflection formula.
- `ml_airy_ai` for |x| > 1e6: degrades as the phase zeta = (2/3)x^{3/2}
  loses digits to rounding; inherent to double precision.

- `ml_lgamma` near its zeros (x=2, x~3.5626): ~700 ULP. The Lanczos sum
  cancels ~3 digits there, and lgamma itself is near zero, so the relative
  error blows up. lgamma is <=1 ULP for x>=8 and <=7 ULP for x in
  [1e-3, 0.5]. A local series expansion around the zeros would be needed
  to do better; this is intrinsic to the Lanczos/Stirling method.
- `ml_gamma` inherits lgamma's error amplified by lgamma(x): ~21 ULP at
  x~6.7 (lgamma~6), consistent with the documented ~1e-15 lgamma error.
