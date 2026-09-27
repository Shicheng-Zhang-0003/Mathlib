# MathLib Design Contract (v11S Closure Candidate — v12R2 Verified)
<!-- v11S CLOSURE IP-22: docs alignment -->

This document defines the architectural and operational boundaries of MathLib v11S.

Any future contributions must adhere to these policies.

---

## Changelog Note (v11S)

* The `MATHLIB_PROFILE_HARDENED` CMake option has been removed.
* Safety checks (NULL bounds, dimension validation, workspace canary checks) are now **unconditional**.
* Core IEEE primitives were hardened with an exact decomposition/recomposition layer.
* `ml_fmod` now uses exact integer-significand modulo for finite nonzero inputs.
* Directed edge-case suites and deterministic fuzz seeds are now part of the closure process.
* Documentation and version metadata were aligned with implemented behavior.

**v12R2 Verification (2026-09-27):**
* Thread-safety audit complete: all mutable static scratch buffers removed.
* Runtime SIMD dispatch guards (`__builtin_cpu_supports`) prevent SIGILL on non-AVX2 hosts.
* 12 Batch-1 modules integrated with thread-safe stack-local workspaces.
* DESIGN_CONTRACT "Stateless & Thread-Safe, No Global State" claim now verified for full tree.

---

## 1. Memory Policy

* **Zero Internal Allocation:** Core APIs (`ml_solve`, `ml_fft_execute`, etc.) must never call `malloc`, `calloc`, or `free`.
* **Client-Provided Scratchpads:** Any operation requiring temporary memory must accept a `ml_workspace_t` bump allocator from the caller.
* **Legacy Isolation:** Heap-heavy legacy modules are quarantined and not part of the core static library.
* **Batch-1 Modules:** Use stack-local fixed-size workspaces (e.g., `double S[17][16]` in `ml_nelder_mead`, `double V[32][32]` in `ml_proj_stiefel`) — no heap allocation in hot paths.

---

## 2. Threading Policy

* **Stateless & Thread-Safe:** Core math functions are pure and stateless.
* **No Global State:** There are no global variables, hidden caches, or thread-local storage.
* MathLib is inherently thread-safe by design.

**v12R2 Verification:** Confirmed by `grep -n "static double\|static float\|static int\|static long double\|static cplx\|static uint" src/*.c` — only static helper *functions* (stateless) and read-only `static const` tables remain. All per-call mutable buffers (`harmonic.c` buf/tmp, `spectral.c` CG/GMRES/Lanczos/SVD scratch, `calculus.c` spline Thomas, `linalg.c` Jacobi `W`, `info.c` MI marginals, `numbertheory.c` prime sieve, `analytic_nt.c` partition/zeta tables) converted to stack-local or heap-per-call.

---

## 3. Determinism Policy

* **Deterministic Within Build/Configuration:** The `SCIENTIFIC` profile guarantees deterministic output for a given build/configuration.
* Different SIMD/FMA paths may differ by bounded ULP unless explicitly validated.
* The build system enforces:
  - `-fno-fast-math`
  - `-ffp-contract=off`
* **Runtime SIMD dispatch** (`ml_cpu_has_fma/avx2/sse41`) uses `__builtin_cpu_supports` — the code path is selected at program start based on host capability, then remains fixed.

---

## 4. Error Handling Policy

* **Core Math:** Returns standard IEEE-754 `double` (`NaN` / `Inf` for domain errors).
* **Structural APIs:** Return `ml_status_t` for explicit failure signaling:
  - `ML_SUCCESS`
  - `ML_ERR_SINGULAR`
  - `ML_ERR_WORKSPACE`
  - `ML_ERR_INVALID_ARG`
  - `ML_ERR_NAN_INPUT`
  - `ML_ERR_INTERNAL`
* **Convergence failure** (optim, ODE, spectral) returns `ML_ERR_SINGULAR` (no dedicated non-convergence code; never returns success un-converged).

---

## 5. Precision Policy

* Validated core transcendentals target **≤ 5 ULP** deviation from ground-truth `mpmath` (80-digit) under the documented domain.
* Exact bitwise operations (`ml_isnan`, `ml_fabs`, etc.) are 100% IEEE-754 exact.
* Fast math functions are approximate and must not be treated as correctly-rounded libm replacements.
* **Oracle tier only** certifies precision (see `docs/PRECISION_CONTRACT.md`). Smoke/edge/modular tolerances are loose and not ULP claims.

---

## 6. Known Limitations

See [`docs/KNOWN_LIMITATIONS.md`](KNOWN_LIMITATIONS.md) for the explicit limitation set.
