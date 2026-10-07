# Despot Round-4 Audit — 2026-10-07

Full mathematical, programming, and operational audit of the entire
MathLib v12R2 tree (34 TUs, 45 public headers, 24 edge suites). Every claim
below is implemented in code and pinned by `tests/test_edge_despot_r4.c`
(21 assertions), the existing gauntlet, or a named probe in
`/tmp/opencode/despot-r4/probes/`.

Policy: all builds and test objects lived in `/tmp/opencode/despot-r4/` only
(`build-sci/graph/emb/asan`, `edge/`, `stack/`, `probes/`). The repo tree
holds sources alone. Ground truth: mpmath 80-dps (integer-ratio exact-binary),
system libm (`glibc`), DLMF identities, and closed-form physics.

## 1. Mathematical audit (vs external truth)

| # | Location | Lie found | Truth | Fix | Measured |
|---|---|---|---|---|---|
| M1 | `integral.c ml_zeta` | `zeta(0)` returned `0` (trivial-zero guard `s<=0` included `s==0`) | DLMF 25.2.2: `zeta(0)=-1/2`; reflection needs `zeta(1)=Inf` | Explicit `s==0 → -0.5`; guard narrowed to `s<0` | 0 ULP; `-2/-4` still `0` |
| M2 | `analytic_nt.c ml_hurwitz_zeta` | `0<s<=1,a!=1` returned NaN (e.g. `hurwitz(0.9,2)`) | Same LD Euler-Maclaurin is the analytic continuation there (pole term `Na^{1-s}/(s-1)`; ≤130x cancellation at s=0.01, harmless in LD) | Shared `ml_hurwitz_em_ld` helper (N=200, B2..B12, Kahan, `(long double)n+(long double)a`) | 0 ULP at (1.1,1)/(1.5,2)/(0.9,2)/(0.01,3.7) vs mpmath 80-dps |
| M3 | `analytic_nt.c` comment | "through B6" / "B2..B6" | `BB[1..6]` = B2,B4,B6,B8,B10,B12; loop `k<=6` = B2..B12 | Comment fixed | — |
| M4 | `complex.c ml_cplx_logarithm` | `log(hypot)` 0.35-ULP amplified 88x into `log(1+0.1i).real` | `0.5*logl(x²+y²)` in LD; `m-1` Sterbenz-exact on [0.5,2] + `log1pl` ≈0.6 ULP | Already in tree (Round-4 base); R4 adds zero-guard + pins | 0 ULP real/imag |
| M5 | `complex.c ml_cplx_sqrt` | `im=sqrt((r-x)/2)` cancelled to 0 for `sqrt(1e300+1e200i)` | Hull-Beebe division form `t=sqrt(2(r±x)); re=t/2 or |y|/t` | Already in tree; R4 fixes comment ("`*2` exact, `r+x` well-conditioned", not "power-of-two scaling") | `im` ~5e49 (1 ULP) vs `cmath` |
| M6 | `complex.c ml_cplx_power` | `w=b*log(a)` rounded in double (3.5e-13 at |w|~2500; 885 ULP case) | `w` in LD + `expl/sinl/cosl` (|wi|≤1e15 window; past ~1e19 direction is unrepresentable — glibc agrees) | Already in tree; R4 dedups 3x product into `ml_cplx_pow_w` + zero-guard | `(1+0.1i)^37` 0 ULP |
| M7 | `integral.c Stirling-LD` | `w=x-0.5L` single rounding (DD path used exact TwoSum) | TwoSum in LD: `w_hi=x-0.5L; w_lo=(x-w_hi)-0.5L` | Fixed; `L=(w_hi*lx+w_lo*lx)-x+…` | Uniform ~1e-19 (below 1e-18 truncation) |
| M8 | `integral.c ml_gamma_p/q` | `P(1/2,Inf)=NaN` via `erf(1e300)` | DLMF 8.2.4: `P(a,Inf)=1, Q(a,Inf)=0` | Already in tree; R4 pins | Exact |
| M9 | `core.c ml_remainder` | Rounded quotient lost ~1000 bits (`remainder(1e300,1000)=3.5e283`) | Integer `2F vs Y` doubling loop; `|y|/2` exact; final `f∓ay` Sterbenz-exact | Already in tree; R4 pins bit-exactness | Bit-exact vs `remainder(3)` |
| M10 | `integral.c zeta s>1` vs `hurwitz(s,1)` | Suspected N=50-double vs N=200-LD inconsistency | Measured: agree to 0-1 ULP (zeta tail B2..B12 + Kahan suffices for s>1) | No change; documented | `s=1.1..2.9` diff ≤4.4e-16 |

### Physics outcomes verified (closed forms + stability contracts)

- ODE: RK4 `dy/dt=y` → `e` to 1e-12 at dt=0.01; Euler converges at dt=0.001.
- SUVAT `s=ut+½at²`, projectile `R=v₀²sin2θ/g`, SHM `x=Acos(ωt+φ)` — exact.
- PDE: `ml_heat_explicit` rejects `r>1/2` (`ML_ERR_SINGULAR`, FTCS Von Neumann);
  `ml_heat_implicit` Thomas unconditionally stable; `ml_wave_leapfrog` rejects
  `cfl>1` (CFL); `ml_poisson_1d` Thomas `-u''=f`; `ml_fem1d_assemble` element-only
  (documented, caller assembles).
- Kelvin `ber(20)=47489.37026504` vs mpmath `…06176`: 5.3e-13 relative (~1e-13
  tier; ULP is meaningless at magnitude 5e4 — docs correctly use relative).
- Bessel K bridge holds 6.8e-15; Airy Ai holds 2.1e-13; digamma/erfc 0 ULP.
- Linalg/spectral/control: LU/Cholesky/QR/Jacobi/SVD contracts + CG/GMRES
  relative `tol*‖b‖` + Kalman/LQR 2x2 scope verified by existing suites
  (despot_check 17/17; modular 4/4). No new discrepancies found.

## 2. Programming audit

- **Static mutable state:** `grep "static ...[" src/*.c` (minus `static const`
  and helpers) → zero hits. `ml_compiler.h` `ML_INLINE static…` + `static inline`
  helpers only. 8-thread smoke (`sin/exp/gamma/mean`) deterministic, no races.
- **Warnings:** 3-profile `-Werror` clean incl. `-Wundef` (added Round-4 base).
  `float.h` for `LDBL_MANT_DIG` present; `trig.c`/`internal/*` use the
  predefined `__LDBL_MANT_DIG__` (no header needed) — both spellings documented
  in `ml_compiler.h`; no silent-fallback risk remains.
- **Sanitizers:** ASan+UBSan oracle 212/0 + R4 21/21 clean. `fanalyzer` clean on
  `core.c`/`integral.c`.
- **Complexity/perf:** Hurwitz 4x `powl` (N=200) accepted for the 0-ULP tier;
  `rising/fact` recomputed per call (12x12 mults — negligible vs `powl`);
  `remainder` ≤2100-iteration doubling loop (µs); `zeta` Borwein LD+Kahan n=30.
  Duplication removed: Hurwitz EM 2x30 lines → 1 helper; pow LD product 3x → 1
  helper (`ml_cplx_pow_w`); Stirling coeffs intentionally duplicated (double vs
  LD precision — documented, not a defect).
- **Edge:** `t==0` guards in `csqrt` effectively dead (post `r==0` return) —
  kept as defense-in-depth; `Gd==0` branches in negative-gamma dead
  (`pos_arg≥2 ⇒ Γ≥1`) — kept, harmless.

## 3. Operational audit

- **TU coherence:** 34/34 identical order in `CMakeLists.txt` / `Makefile` /
  `run_all_tests.py LIB_SOURCES` / `tests/run_edge_tests.sh` / `src/*.c` on disk.
  (`rg src/` 36-hit false positive from 2 guard strings already documented.)
- **Artifacts:** root `*.o`/`*.d` removed (gitignored yet present; violated the
  `/tmp`-only evidence policy in `GATE_V12R2.md`). `run_edge_tests.sh` now honors
  `MATHLIB_EDGE_BUILD_DIR` (was hardcoded `build/`); audit edge objects went to
  `/tmp/opencode/despot-r4/edge/`.
- **Fuzz figure:** god-mode measured **61393/0** (`--seed 123456789`) on the
  committed + R4 tree. The 61473 figure in Round-2 drafts is not reproducible
  (`tests/fuzz_god_mode.c` unchanged since 02/08/26); corrected everywhere to
  the measured 61393. Boundary 25/0 holds.
- **Stack:** `-fstack-usage` max `ml_fft_real` 65,600 B; 32-33 KB frames
  (`haar`, `spline`, `mi_discrete`, `jacobi_eigen`) + 27 KB GMRES + 16 KB PDE
  remain thread-safe but heavy for embedded — heap-workspace API still deferred
  (documented, not hidden).
- **Docs:** 5 stale spots fixed (fuzz 61473→61393; Kelvin 0.5%→1e-13;
  Hurwitz stub scope; B6→B12; edge 23→24 suites). This file + `release_notes.md`
  Round-4 + `KNOWN_LIMITATIONS.md` Round-4 + `PRECISION_CONTRACT.md` R4 guard
  note the new state. `GATE_V12R2.md` (2026-09-30 evidence) left historical.

## 4. Validation (this tree, `/tmp` objects only)

- Builds: SCIENTIFIC / GRAPHICS / EMBEDDED `-Werror` clean; ASan+UBSan clean.
- Oracle: 212/0 worst 0 ULP (gate ≤5 ULP). Modular 4/4. Smoke 30/0.
- God-mode 61393/0 (`--seed 123456789`). Boundary 25/0.
- Edge: 24/24 (23 existing + `test_edge_despot_r4` 21/21), incl. ASan R4 clean.
- `accuracy_audit` 371/0 (K/Airy 65536-ULP guard tier, measured 0-8 ULP).

## 5. Remaining (honest, not hidden)

- `ml_hurwitz_zeta(s<=0,a!=1)` fail-loud NaN (needs Hurwitz functional equation).
- Batch-1 no oracle ULP tier (structural).
- Core LD paths collapse where `long double==double` (K/Ai exempt by quadrature).
- Y0/Y1 absolute ~1e-13 near zeros; large-x trig crossings absolute regime.
- Heavy stack frames (heap-workspace API deferred).
- Proven <0.5 ULP for all binary64 still needs Ziv + worst-case search (Table
  Maker's Dilemma; see `ULP_PUSH.md`).
