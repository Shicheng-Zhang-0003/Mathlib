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

- **A1 closure freeze is in effect.** No new modules, APIs, or math
  families are permitted. Only closure fixes, tests, validation,
  and documentation alignment are allowed.

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
