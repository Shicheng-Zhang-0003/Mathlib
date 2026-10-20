# MathLib subject coverage contract

## Purpose

MathLib aims to collect mathematical algorithms and contest mathematics in
one C library. “Coverage” means a callable implementation with stated input
domain, output semantics, complexity/bounds, and a validation tier. A header,
example, identity checker, or named showcase alone does not establish coverage
of a whole subject.

This map uses the major subject areas of MSC2020 as a navigation aid. MSC is a
research-literature taxonomy maintained jointly by Mathematical Reviews and
zbMATH; it is not a finite checklist of all mathematical knowledge. See the
[MSC2020 classification](https://msc2020.org/MSC_2020.pdf) and the current
[MathSciNet MSC page](https://mathscinet.ams.org/msc/msc.html).

## Coverage labels

- **Implemented** — a useful family of algorithms exists, with public API and
  regression coverage. This label does not itself claim formal proof or full
  subject coverage.
- **Partial** — selected algorithms/examples exist; broad use of the subject
  is not supported.
- **Absent** — no meaningful public computational API exists yet.
- **Certified** — a result carries a checkable proof/certificate or a
  mathematically justified enclosure. No entire subject is certified merely
  because some routines return status codes.

Numerical quality is tracked separately in `PRECISION_CONTRACT.md` and
`API_STATUS.md`: a regression test is not an oracle-tier accuracy guarantee.
The level-by-level school scope is tracked separately in
`PREUNIVERSITY_COVERAGE.md`.

## Current map

| MSC area | Current status | Representative current surface | Main missing coverage |
|---|---|---|---|
| 00–01 General topics, history, exposition | **Partial** | Recreational/contest routines and engineering docs | Searchable mathematical knowledge base, history, glossary, proofs, and curricula |
| 03 Logic and foundations | **Absent** | None beyond C-level validation | Propositional/first-order logic, set/type theory, proof objects, proof checking, computability, model theory |
| 05–08 Combinatorics, order, algebraic systems | **Partial** | Contest combinatorics, graph algorithms, basic ordered/data structures | Extremal/probabilistic/additive combinatorics, designs/matroids, lattices, semigroups, general algebraic structures |
| 11–20 Number theory and algebra | **Partial** | Integer arithmetic, selected modular/contest NT, prime fields F_p, exact Q[x] arithmetic and division certificates, resultants, discriminants, polynomial roots, dynamic exact integers/rationals | Extension fields, algebraic number theory, abstract algebra, Galois/representation theory, multivariate Gröbner bases and broader symbolic algebra |
| 14 Algebraic geometry | **Absent** | Resultants and polynomial tools are adjacent building blocks | Varieties/schemes, ideals, Gröbner/elimination workflows, curves, divisors, arithmetic geometry |
| 15 Linear and multilinear algebra | **Implemented, bounded** | Dense factorizations, SVD, iterative solvers, matrix functions, tensor views | Sparse/distributed methods, exact and structured matrices, multilinear/tensor decompositions, stronger conditioning/error contracts |
| 18–19 Category theory, homological algebra, K-theory | **Absent** | None | Categories/functors, chain complexes, derived functors, spectral sequences, K-theory |
| 20–22 Groups and topological/Lie groups | **Partial** | Small symmetric-group hook/class formulas, quaternions, 2x2 control examples | General groups, presentations, subgroup/representation algorithms, Lie groups/algebras and harmonic analysis on groups |
| 26–33 Real/complex analysis, measure, potential theory, special functions | **Partial** | Numerical calculus/integration, complex arithmetic, selected contour routines and special functions | Rigorous limits/convergence, Lebesgue measure/integration, distributions, several complex variables, potential theory, broad analytic continuation and certified special functions |
| 34–39 ODE/PDE, dynamics, functional/difference equations | **Partial** | Scalar/system ODE solvers, sample BVP/PDE routines, SDE, Lorenz/Rössler, selected functional equations | Stiff/DAE/event systems, broad nonlinear PDEs, global FEM assembly, weak forms, bifurcation/ergodic theory, functional/difference equation solvers |
| 40–45 Series, approximation, harmonic analysis, transforms, integral equations | **Partial** | Series acceleration/Pade, FFT, DCT/DST, Haar, selected integral routines | Wavelet libraries, distributions, generalized Fourier/spectral methods, integral-equation solvers, approximation theory and certified error estimates |
| 46–49 Functional/operator analysis, variational methods, optimization | **Partial** | A few Banach/Hilbert examples, gradient methods, low-dimensional control | Banach/Hilbert/operator theory as a computational system, calculus of variations, convex analysis, LP/QP/conic/integer optimization and certificates |
| 51–58 Geometry, topology, global analysis | **Partial** | Euclidean/computational geometry, curvature examples, mod-2 simplicial Betti numbers | Convex/projective geometry, point-set topology, homotopy/cohomology, persistent homology, bundles/characteristic classes, knots, global differential geometry |
| 60–62 Probability and statistics | **Partial** | Common distributions, descriptive statistics, inverse CDFs, MCMC and simple models | General stochastic-process framework, multivariate inference, regression/GLMs, time series, survival/causal/nonparametric statistics and validated RNGs |
| 65 Numerical analysis | **Partial** | Floating-point kernels, quadrature, root finding, ODE/PDE and matrix solvers, outward-rounded basic intervals | Arbitrary-precision interval/ball transcendental functions, proof-producing solver bounds, reproducible cross-platform rounding, adaptive/sparse scalable methods |
| 68 Computer-science mathematics | **Partial** | Graph algorithms, FFT, exact integer algorithms, basic information routines | Automata/formal languages, complexity, cryptographic math (with safe primitives), coding theory, SAT/SMT, symbolic verification |
| 70–86 Mechanics and mathematical physics | **Partial** | Quaternion mechanics, basic kinematics, control, SDE, example dynamical systems | Rigid-body/continuum/fluid mechanics, elasticity, electromagnetism, thermodynamics, quantum/statistical mechanics, relativity and validated simulation |
| 90–94 Operations research, economics, biology, systems, information | **Partial** | Graph flow/matching, MCMC, Markov/SIR/queueing, small control and info-theory routines | General optimization/decision/game theory, econometrics/finance, population/ecological models, information/coding theory and communications |
| 97 Mathematics education | **Absent** | API/test documentation only | Curriculum map, course sequences, authored exercises, hint/proof pathways, assessment and misconception checks |

## Contest coverage contract

Contest support is a cross-cutting layer, not a synonym for `ml_contest`.
Current coverage includes many elementary identities, integer sequences,
selected olympiad number theory, graph tasks, Euclidean constructions and
verification routines. `ml_contest.h` currently exports 42 functions across
recreational, school, olympiad, and Putnam-style tiers. It does not yet include
a broad problem corpus, proof
strategy library, proof checker, or verified solution certificates. The
largest missing contest families are synthetic geometry, inequalities,
functional equations, extremal/combinatorial arguments, probability, and
proof-writing. A numeric answer checker cannot certify a proof.

The pre-university implementation scope and remaining school-curriculum gaps
are listed in `PREUNIVERSITY_COVERAGE.md`.

## Shared foundation milestones

These foundations unlock multiple subject areas and should be built before
adding isolated named APIs:

1. **Exact arithmetic — first usable layer implemented:** dynamically sized
   integers and canonical rationals now live in `ml_exact.h`; finite fields
   F_p arithmetic now lives in `ml_gfp.h`; extension fields and
   algebraic-number representations follow.
2. **Certified numerics — first usable layer implemented:** `ml_interval.h`
   provides outward-rounded basic arithmetic and exact rational-to-binary64
   enclosures. Transcendental interval functions and checkable solver
   residual/error certificates follow.
3. **Symbolic representation — first exact polynomial layer implemented:**
   `ml_polyq.h` supplies univariate Q[x] arithmetic, evaluation, division,
   gcd, and exact identity checking. Immutable expression DAGs, multivariate
   polynomials, and rewrite systems follow.
4. **Proof/certificate layer — first algebraic witness verifier implemented:**
   quotient/remainder certificates are checked exactly in Q[x]. Verifiers
   for optimization witnesses and contest solution objects follow.
5. **Subject modules:** abstract algebra, topology, analysis, optimization,
   and mathematical physics built on those shared primitives.

Each milestone must add a public API contract, independent regression cases,
and at least one independent cross-check strategy. New areas remain marked
**Partial** until those conditions hold; broad “all mathematics” claims must
not be inferred from module count.

## Current bounded examples

- `ml_bigint_t` is fixed at 4096 bits; it is not dynamically arbitrary
  precision. The new `ml_mpz_t` / `ml_mpq_t` APIs are dynamic and limited by
  available memory instead.
- `ml_topology` computes selected finite simplicial invariants over F2; it is
  not a general topology or homology package.
- `ml_measure`, `ml_funcanal`, and `ml_diffgeom` contain selected instances
  and verifiers, not comprehensive subject frameworks.
- The documented FEM routine returns one element's matrices rather than
  assembling/solving a global finite-element problem.
- New numerical extensions are regression-tested on pinned cases but do not
  inherit oracle-tier accuracy guarantees.
