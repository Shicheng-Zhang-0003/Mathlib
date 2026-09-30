# API Guide (v12R2 / V1.2-RC2, 34 TUs)

Public headers are grouped by subsystem (all carry `extern "C"` guards):

core
- fundamental types (`ml_types.h`: `ml_status_t` now includes
  `OVERFLOW/CONVERGENCE/UNSUPPORTED`, ABI-additive)
- IEEE helpers (`ml_core.h`, `bitwise_fp.h`)
- constants, compiler/SIMD/profile routing (`ml_compiler.h`, `profiles.h`,
  `simd*.h`, `fast_math.h`, `cpu_dispatch.h`, `version.h`)

transcendental (measured 0 ULP on oracle 212 + wide grids, x86-64; gate ≤5 ULP)
- trigonometric functions (`ml_trig.h`: LD reduction + sinl/cosl)
- exponential/logarithmic + pow (`ml_exp_log.h`: LD split/pow, integer LD)
- complex (`ml_complex.h`: Annex-G sqrt)
- gamma/lgamma/digamma (`ml_integral.h`: LD shift/Taylor/direct/reflection)
- kelvin (`ml_kelvin.h`: ber/bei NaN at Inf, XMAX)

numerical
- integration, statistics (+covariance NaN, relative inv stops),
  combinatorics, numerical solvers (NaN fail-loud)
- optimization/optim_n, ode/ode_sys (BE2 fail-loud), spectral (relative tol),
  stats_inv, sde (`_ctx` variants, forward-only), pde (Dirichlet-only),
  harmonic (Haar contract), mcmc (rejection accounting), manifold (unit gates),
  info (cross-entropy +Inf), analytic_nt (hurwitz stub), control (Hurwitz margin),
  orthogonal/polynomial/quadratics, fixed-point (rounded constants),
  numbertheory (crt2 `%q`, mobius 2=UNRESOLVED, exact mult_order)

signal processing
- FFT (pow2 Cooley-Tukey), DSP primitives, DCT/DST (`ml_dst1` alias for DST-I),
  Laplace table

linear algebra
- matrix and vector operations (LU/Cholesky/QR, Jacobi SINGULAR on
  non-convergence, determinant finite guard)

Consumers should include only required subsystem headers.
Status per header: `docs/API_STATUS.md`. Accuracy tiers: `docs/PRECISION_CONTRACT.md`.
Measured ULP push: `docs/ULP_PUSH.md`. Gate evidence: `docs/GATE_V12R2.md`.
