# Doctorate Guide (43 TUs): the computable slice of graduate mathematics

Five modules plus a PNT extension. Rule throughout: implement what is
*computable* (solvers, verifiers, exact integer algorithms); never claim
proofs. No oracle ULP certification (PRECISION_CONTRACT rule 4).

## Measure (`ml_measure.h`, 7 fns)

Lebesgue Dirichlet showcase (integral exactly 0), simple-function sums
(LD), discrete L^p geometry, Hölder / Minkowski / Markov instance
verifiers on uniform samples (fail closed, tol-gated).

## Functional analysis (`ml_funcanal.h`, 4 fns)

Banach contraction solver with caller-certified L and an a posteriori
`L^n/(1-L)` bound (CONVERGENCE-loud), Neumann series inverses
(`||I-A||_inf < 1`, n<=16), Hilbert best approximation (MGS with
scale-aware rank test + back substitution), Rayleigh maxima wrapping the
staged power iteration.

## Topology (`ml_topology.h`, 6 fns)

Simplicial mod-2 Betti numbers from closed complexes (validated closure,
heap bitset elimination — a 1024² bitset is 128KB, never stack),
chain-complex checks + Betti, Euler characteristic, and exact Bareiss
rank/determinant over int64 (LD-exact below 2^64 with 2^64 guards;
OVERFLOW past it). Pinned: filled (1,0,0), hollow (1,1,0), tet surface
(1,0,1).

## Differential geometry (`ml_diffgeom.h`, 4 fns)

Graph/parametric curvature, Frenet frames (per-derivative step sizes —
h=2e-8 loses ~60% to rounding on second derivatives in double precision;
fixed during review), sphere K=1/9 and signed H pinned, plane K=H=0.
Straight motion is SINGULAR with kappa 0 (normal plane undefined).

## Computational algebra (`ml_algebraic.h`, 7 fns) + PNT (`ml_analytic_nt.h`)

Sylvester resultants + discriminants via staged Bareiss (exact),
Sturm isolation (LD tier, deg<=8), Descartes (exact), Specht dimensions
and S_n class sizes (hook formula, exact), Rodrigues so(3). Plus Ei (LD
series), li, Chebyshev psi (heap sieve) — PNT illustration pinned:
li(1e6)≈78627.55 vs pi(1e6)=78498.

## Explicitly deferred (still doctorate, not yet computable here)

Gröbner bases, spectral PDE methods, EKF/MPC control, Gaussian processes,
Bellman-Ford, Floyd-Warshall, 2D Poisson, category-theoretic diagram
checkers, Lie brackets beyond Rodrigues.
