# Despot Round-7 Audit — extension math and package integration

This audit follows the user-requested mathematical, programming, and
operational review of the v12R3 working tree. It is an engineering audit
with regression and differential evidence, not a formal proof of every
routine over every binary64 input.

## Scope

- Reviewed the current C99 library, its 54 translation units, public headers,
  CMake/Makefile and Python/shell test entry points, package export, and the
  new extension APIs.
- Ran targeted mpmath comparisons for the new special-function paths and
  added directed regression cases in `tests/test_edge_despot_r7.c`.
- Built SCIENTIFIC, GRAPHICS, and EMBEDDED configurations, plus GCC and Clang
  builds. Detailed final gate results are recorded after the complete suite
  finishes.

## Findings repaired

### Numerical behavior

- General-order Bessel J and I series mishandled negative noninteger orders;
  the gamma sign is now included and summation uses compensated long-double
  accumulation. Integer J uses stable recurrence, while I/Y/K use recurrence
  or transition-region quadrature where the old asymptotic/series choices
  were unstable.
- Jacobi inverse incomplete-F used unconstrained Newton steps and failed near
  `m=1`; it now brackets and bisects the amplitude.
- Airy Bi and derivative paths used inaccurate local series; Bessel
  identities now provide the positive and negative argument branches.
- Dawson's asymptotic expansion could continue after its least term. It now
  stops at the least term, and its small-argument series uses compensated
  long-double summation.
- Struve H0's former approximation lost accuracy for larger oscillatory
  arguments. It now uses its integral representation over the documented
  range.
- Si/Ci asymptotic expansions stop at the least term; invalid negative-real
  Ci/Chi inputs return NaN; spherical functions now preserve their negative
  argument parity and infinite limits.
- High-order polygamma uses a convergent tail instead of a divergent
  Bernoulli asymptotic expansion.

Unsupported numerical regions fail explicitly: near-integer noninteger Y/K,
unsafe asymptotic transition zones, Struve H0 outside `|x| <= 10000`, and
series that do not converge within their iteration budget return NaN. These
are documented in `include/mathlib/ml_special.h`; they remain scope limits,
not evidence of full-domain support.

### Programming and package behavior

- The installed static CMake target now propagates the sanitizer runtime link
  flags when the library archive was built with ASan/UBSan. An installed
  consumer fixture compiles, links, and calls a public extension API through
  `mathlib::mathc`.
- The compiler portability header maps the used builtin math functions for
  MSVC and provides portable int64 overflow helpers. This was reviewed and
  built on GCC/Clang; it has not been compiled with MSVC itself.
- Regression coverage was extended for the repaired edge cases, including
  negative fractional Bessel orders, large-order recurrence, Jacobi near
  `m=1`, Airy branches, Dawson, Struve range, and spherical parity.

## Evidence and limits

- mpmath differential sweeps found maximum sampled relative errors of about
  `3e-11` for the general-order Bessel J grid, `2e-14` for I, `1e-9` for Y,
  and `6e-15` for K. Sampled Airy Bi/derivatives, Dawson, Fresnel, and other
  special families also matched their references within the regression
  tolerances. These are finite sample grids, not worst-case proofs or an
  oracle-tier certification.
- Existing core oracle guarantees are unchanged. The new extension modules
  have no comprehensive oracle tier yet; directed tests demonstrate only
  their pinned cases. See `docs/PRECISION_CONTRACT.md`.
- `ml_bigint` remains fixed at 4096 bits. Its xorshift generator and
  Miller-Rabin helper are not cryptographic. Several numerical extension
  APIs intentionally have bounded domains or iteration limits.
- The GCC analyzer reports an advisory imprecise-floating-point warning on a
  bounded allocation-size multiplication in `src/analytic_nt.c`; the size is
  range-checked before allocation. It is not a sanitizer or runtime failure.
- Platform evidence covers Linux GCC 13.3 and Clang 18.1.3. The three
  profiles were built, but actual Windows/MSVC and non-x86 runtime behavior
  remain unverified.

## Sandbox and temporary-file diagnosis

The full runner configures and builds under `build/`, and edge binaries also
go there by default. GCC itself uses `TMPDIR` (or `/tmp` if unset) for
intermediate assembly files. Thus a filesystem sandbox that permits writes
to the repository but blocks `/tmp` can produce compiler failures even
though the project build paths are in the repository. Setting `TMPDIR`,
`TMP`, and `TEMP` to `build/.tmp` keeps compiler temporaries in the checkout.
The optional soak scripts also default their evidence logs to `/tmp`; set
`MATHLIB_SOAK_DIR` to a directory under the checkout when enabling soak.

There was one earlier audit build invocation before these overrides were
applied; GCC created transient assembly files in `/tmp`. I noticed and
stopped that run, did not inspect or clean anything outside the repository,
and redirected all later build/test temporaries into `build/.tmp`.

## Validation record

- `run_all_tests.py --jobs 4`: **48 passed, 0 failed**. This includes 37
  edge executables, the 973-assertion R7 regression executable, fuzzers,
  oracle check, Python runner, and installed CMake package consumer. The
  optional soak loop was not run.
- SCIENTIFIC sanitizer build and installed consumer passed. GRAPHICS and
  EMBEDDED builds passed; strict Clang 18.1.3 build passed. Public-header
  C++17 syntax checks passed with and without AVX2.
- GCC analyzer build completed. Its only reported advisory was the bounded
  allocation-size multiplication in `src/analytic_nt.c:217`.
- CMake, `run_all_tests.py`, and `tests/run_edge_tests.sh` each enumerate the
  same 54 translation units. `git diff --check` passed.
