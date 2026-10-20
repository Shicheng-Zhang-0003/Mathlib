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
  run_edge_tests.sh canonical lists must match (59 TUs in the current tree).

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

## v12R3 Development Limits

The eleven extension modules added in the current 54-TU v12R3 tree have no oracle-
tier certification. `ml_special` uses bounded integer orders/degrees and
fixed-iteration approximations; it is not a general-purpose libm replacement.
`ml_bigint` is fixed at 4096 bits, and its deterministic xorshift generator,
modular arithmetic, and Miller-Rabin helper are not suitable for cryptographic
key generation or secret-dependent computations. `ml_math_extra` documents
library-specific behavior and does not promise full platform-libm equivalence.

An independent audit (mpmath 50-dps references, ULP-measured) repaired the
following. `tests/test_edge_accuracy_audit.c` pins the results as a
regression guard (371 assertions; was 282 before Round-3). Per
`PRECISION_CONTRACT.md` this is a
regression guard, not an oracle-tier precision certification.

| Function | Before | After | Note |
|---|---|---|---|
| `ml_erfc` | ~1064 ULP | <=0.7 ULP | Laplace CF + convergent series, long double |
| `ml_bessel_j0/j1/y0/y1` | ~1e-3 rel | <=4e-13 rel | Hankel P/Q, least-term asymptotics |
| `ml_bessel_i0/i1` | exact | exact | |
| `ml_bessel_k0/k1` | ~1e9 rel | <=6.8e-15 rel | Round-3: K-quadrature bridge for 4<=x<16 closes the 8-10 hole (was 2e-9) |
| `ml_airy_ai` | wrong coefficients | <=2.1e-13 rel | Round-3: Ai=sqrt(x/3)/pi*K_1/3 quadrature for 2.5<x<8.5 closes the 5-6 hole (was 5e-9) |
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
- `ml_hurwitz_zeta` scope (Round-4 update): s>0 fully implemented via LD
  Euler-Maclaurin N=200 B2..B12 (0 ULP vs mpmath 80-dps, including
  0<s<=1,a!=1 continuation). Only s<=0,a!=1 stays a fail-loud NaN stub
  (needs the Hurwitz functional equation). `ml_zeta(0)=-1/2` exact.
- Batch-1 still no oracle ULP certification (bisection inverses, Krylov
  without preconditioning, xorshift RNG, Dirichlet-only PDEs, 2x2 control).
- `long double==double` platforms (MSVC/ARM): core LD paths collapse to DD;
  K0/K1 and Airy-Ai are exempt (quadrature is LD-independent, measured
  2.6e-16/2.4e-15 in emulated double).

## Round-4 despot (2026-10-07)

Full mathematical/programming/operational re-audit; see
`docs/DESPOT_R4_AUDIT.md`. Fixed and pinned by
`tests/test_edge_despot_r4.c` (21 assertions): `ml_zeta(0)=-1/2` (was 0),
Hurwitz 0<s<=1 continuation (was NaN stub) via shared EM helper, B2..B12
comment fix, complex sqrt Hull-Beebe wording, complex pow shared w helper +
zero-guard, Stirling-LD TwoSum, edge-runner `/tmp` override, root artifact
cleanup. Measured: oracle 212/0 worst 0 ULP, god 0 failed, boundary 25/0,
accuracy_audit 371/0, R4 21/21, 3-profile `-Werror` + ASan/UBSan clean.

Closed by the ULP push (measured exact-binary, x86-64):
`ml_gamma` 6.7/0.1/0.001/-0.5/-0.1/-0.9 (21/3/5/1/1/1 → 0 ULP),
`ml_lgamma` 6.7/0.1/0.001 (4/1/1 → 0 ULP), pow 41→0 (100-case grid),
sin/cos benign 1→0 (25-case grid). Oracle worst 5→0 ULP on the 212 grid.
See `docs/ULP_PUSH.md`. Unchanged (still limited, not re-measured this
push):

- ~~`ml_bessel_k0/k1` in 8 < x < 10: ~1e-9 relative.~~ **CLOSED (Round-3)**
  by the K-quadrature bridge; worst over [1e-3, 300] is 6.8e-15.
- ~~`ml_airy_ai` in 5 < x < 6: ~5e-9 relative.~~ **CLOSED (Round-3)** by the
  `K_1/3` bridge; worst over a 0.1-step sweep of [-14, 200] is 2.1e-13.
- ~~`ml_digamma` for x < 0: ~10 ULP at x=-0.5.~~ **CLOSED (Round-3)** —
  exact half-integer shortcut + long-double `cot(pi x)`; now 0 ULP.
- `ml_bessel_y0/y1`: worst **absolute** error 1.1e-13 (Y0 @13, Y1 @12.75).
  Relative error reaches 1e-10 at zeros of Y1 (`Y1(11.75) = -1.96e-4`).
  Intrinsic: the ascending series cancels 4.7 digits; near a zero, no
  elementary method does better. Same regime as large-x trig.
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

Remains after Round-2 (all closed in Round-3 except Batch-1/LD-collapse/stacks):
- ~~Kelvin 2-term P/Q ~0.5% at 20 (bei 20.1 125562 vs 126161 true). Needs
  DLMF 10.67.5 R/S to 1/x^5 or I0/K0 uniform continuation.~~ **CLOSED
  (Round-3)** by DLMF 10.67.3-4 full least-term sums (~1e-13 relative).
- ~~Airy 5<x<6 Taylor/asym hole ~3e-10..3e-9 (LD halves it); K 8<x<10 ~2e-9.
  Both need Temme uniform expansion (deferred to v12A2).~~ **CLOSED
  (Round-3)** by K-quadrature bridge (K 6.8e-15, Ai 2.1e-13).
- ~~Y1 near 14: 25k ULP (crossover cause, same family).~~ **CLOSED (Round-3)**
  crossover 14->13 (1.08e-13); residual at zeros is intrinsic absolute error.
- Batch-1 no oracle; long double==double collapse (core paths; K/Ai now
  exempt); heavy stacks — unchanged (see Round-3/4 notes below).
- `ml_hurwitz_zeta` scope (Round-4): s>0 fully implemented via LD
  Euler-Maclaurin N=200 B2..B12 (0 ULP vs mpmath 80-dps, including
  0<s<=1,a!=1 continuation). Only s<=0,a!=1 stays a fail-loud NaN stub
  (needs the Hurwitz functional equation). `ml_zeta(0)=-1/2` exact.
- `ml_hurwitz_zeta` scope (Round-5): integer test guarded for |v|>=2^53;
  dead a<=0 limb removed; s<=0,a!=1 stub stays (documented).

## Round-4 despot (2026-10-07)

Full mathematical/programming/operational re-audit; see
`docs/DESPOT_R4_AUDIT.md`. Fixed and pinned by
`tests/test_edge_despot_r4.c` (21 assertions): `ml_zeta(0)=-1/2` (was 0),
Hurwitz 0<s<=1 continuation (was NaN stub) via shared EM helper, B2..B12
comment fix, complex sqrt Hull-Beebe wording, complex pow shared w helper +
zero-guard, Stirling-LD TwoSum, edge-runner `/tmp` override, root artifact
cleanup. Measured: oracle 212/0 worst 0 ULP, god 0 failed, boundary 25/0,
accuracy_audit 371/0, R4 21/21, 3-profile `-Werror` + ASan/UBSan clean.

## Round-5 despot (2026-10-10, Edmonton)

Full math/programming/operational re-audit; see `docs/DESPOT_R5_AUDIT.md`.
Fixed and pinned by `tests/test_edge_despot_r5.c` (42 assertions): CPU
non-GNUC scalar fallback + honest `ml_cpu_has_avx2`, fixed-point shift/INT64_MIN
hardening, LDBL-guard unification, cross3/transpose aliasing, iterative
overflow-checked egcd, MSVC builtin mapping, CONVERGENCE/OVERFLOW precision,
DoS caps (majorizes/ess/DCT/Pollard), RNG warmup unification, ODE `>=`
arrival, Hurwitz large-guard, `ultimate_fuzzer` CMake target, gate/runner
`/tmp` honors, soak `/tmp` logs, root-binary cleanup. Measured: oracle 212/0
worst 0 ULP, god 61473/0 AVX2 (0 failed is the gate), boundary 25/0,
accuracy_audit 371/0, R4 21/21, R5 42/42, edge 25/25, 3-profile `-Werror` +
ASan/UBSan clean.

## Round-7 despot (2026-10-19)

See `docs/DESPOT_R7_AUDIT.md`. Fixed general-order Bessel sign/stability,
Jacobi inversion near m=1, Airy Bi/derivative branches, Dawson's least-term
handling, Struve H0 oscillatory range, high-order polygamma, and related
domain/parity cases. `tests/test_edge_despot_r7.c` pins 973 regression
assertions. Full gate: 48/48 passed; SCIENTIFIC, GRAPHICS, and EMBEDDED builds
passed, as did strict Clang and the installed CMake package consumer.

These extensions still lack comprehensive oracle-tier certification. General
Bessel Y/K near-integer and unsafe transition regions, Struve H0 outside
`|x|<=10000`, and unconverged bounded series return NaN by contract. Actual
MSVC/Windows and non-x86 runtime behavior were not tested. The GCC analyzer
has one advisory on a size-checked allocation-size multiplication in
`src/analytic_nt.c`.

## Exact-foundations increment (2026-10-19)

The current tree adds memory-bounded dynamic integers/rationals, exact
univariate Q[x], finite prime fields up to 32-bit characteristic, basic
outward-rounded binary64 intervals, and an exact polynomial-division witness
verifier. This is the first reusable foundation layer, not a general CAS or
proof assistant. Polynomial arithmetic is dense and univariate; field
extensions, multivariate Gröbner bases, symbolic expression DAGs, rigorous
transcendental intervals, and general proof certificates remain future work.
