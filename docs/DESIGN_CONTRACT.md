# MathLib Design Contract (v11S Closure Candidate)
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

---

## 1. Memory Policy

* **Zero Internal Allocation:** Core APIs (`ml_solve`, `ml_fft_execute`, etc.) must never call `malloc`, `calloc`, or `free`.
* **Client-Provided Scratchpads:** Any operation requiring temporary memory must accept a `ml_workspace_t` bump allocator from the caller.
* **Legacy Isolation:** Heap-heavy legacy modules are quarantined and not part of the core static library.

---

## 2. Threading Policy

* **Stateless & Thread-Safe:** Core math functions are pure and stateless.
* **No Global State:** There are no global variables, hidden caches, or thread-local storage.
* MathLib is inherently thread-safe by design.

---

## 3. Determinism Policy

* **Deterministic Within Build/Configuration:** The `SCIENTIFIC` profile guarantees deterministic output for a given build/configuration.
* Different SIMD/FMA paths may differ by bounded ULP unless explicitly validated.
* The build system enforces:
  - `-fno-fast-math`
  - `-ffp-contract=off`

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

---

## 5. Precision Policy

* Validated core transcendentals target **≤ 5 ULP** deviation from ground-truth `glibc libm` under the documented domain.
* Exact bitwise operations (`ml_isnan`, `ml_fabs`, etc.) are 100% IEEE-754 exact.
* Fast math functions are approximate and must not be treated as correctly-rounded libm replacements.

---

## 6. Known Limitations

See [`docs/KNOWN_LIMITATIONS.md`](KNOWN_LIMITATIONS.md) for the explicit v11S limitation set.
