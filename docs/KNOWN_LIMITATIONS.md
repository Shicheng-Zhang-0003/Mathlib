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

- Large-x trig near zero: absolute error ~2e-20 (e.g. cos(4.5e7)), ULP
  meaningless at crossings. Fractional accumulation is LD but not exact
  integer; full-range proof needs exact fractionals (deferred).
- `long double==double` platforms (MSVC/ARM): LD paths collapse to DD;
  measured 0-ULP numbers are x86-64 only. Guarantees stay at the ≤5 ULP gate.
- Batch-1 still no oracle ULP certification (bisection inverses, Krylov
  without preconditioning, xorshift RNG, Dirichlet-only PDEs, 2x2 control).

Closed by the ULP push (measured exact-binary, x86-64):
`ml_gamma` 6.7/0.1/0.001/-0.5/-0.1/-0.9 (21/3/5/1/1/1 → 0 ULP),
`ml_lgamma` 6.7/0.1/0.001 (4/1/1 → 0 ULP), pow 41→0 (100-case grid),
sin/cos benign 1→0 (25-case grid). Oracle worst 5→0 ULP on the 212 grid.
See `docs/ULP_PUSH.md`. Unchanged (still limited, not re-measured this
push):

- `ml_bessel_k0/k1` in 8 < x < 10: ~1e-9 relative. The ascending series
  cancels like exp(x^2/4) and the asymptotic bottoms out at exp(-2x); the
  elementary method cannot do better. Temme's uniform expansion is required.
- `ml_airy_ai` in 5 < x < 6: ~5e-9 relative. Same reason: the Taylor series
  loses digits to cancellation while the asymptotic is still converging.
- `ml_digamma` for x < 0: ~10 ULP at x=-0.5 via the reflection formula.
- `ml_airy_ai` for |x| > 1e6: degrades as the phase zeta = (2/3)x^{3/2}
  loses digits to rounding; inherent to double precision.

- `ml_lgamma` near its zeros (x=2, x~3.5626): historically ~700 ULP via
  Lanczos cancellation; local Taylor (x=1,2 radius 0.15, ζ to 25) now holds
  the measured grid at 0 ULP, but the general near-zero relative-error
  caution stands (value near zero → relative blows up; absolute is what is
  controlled). A wider certified expansion is still deferred.
- `ml_gamma` historically inherited lgamma's error amplified by lgamma(x)
  (~21 ULP at x~6.7); LD direct paths now hold the measured grid at 0 ULP
  (x86-64). General proof still deferred (Ziv + worst-case search).

## Round-2 despot (2026-10-01) — what was fixed vs what remains

Fixed and pinned: zeta OOB (ASan abort), kelvin sign flip (ber -pi/8 split),
prime_pi odd undercount, acosh 458k ULP, softplus 57k ULP, normal_inv sigma
stall, hurwitz B2 10k ULP, catalan C35/C36 exact via recurrence, log1p coeff,
Airy LD Taylor, QR NaN-poison, hurwitz-margin alias, 35->34 TU docs.

Remains (documented limits, not defects):
- Kelvin 2-term P/Q ~0.5% at 20 (bei 20.1 125562 vs 126161 true). Needs
  DLMF 10.67.5 R/S to 1/x^5 or I0/K0 uniform continuation.
- Airy 5<x<6 Taylor/asym hole ~3e-10..3e-9 (LD halves it); K 8<x<10 ~2e-9.
  Both need Temme uniform expansion (deferred to v12A2).
- Y1 near 14: 25k ULP (crossover cause, same family).
- Batch-1 no oracle; long double==double collapse; heavy stacks — unchanged.
