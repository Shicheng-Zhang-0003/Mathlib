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
