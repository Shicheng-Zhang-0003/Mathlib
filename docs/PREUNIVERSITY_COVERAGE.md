# Pre-university mathematics coverage

## Scope and standard

This inventory covers the computational topics commonly taught from primary
school through advanced secondary school, including precalculus, introductory
calculus where offered, and school-level mathematics contests. National
curricula differ in order and depth, so these are topic families rather than a
single official syllabus. “Covered” means that a useful callable tool exists;
it does not mean the library teaches the topic, supplies exercises, or checks a
written proof.

The current API is strongest as a numerical and exact-arithmetic toolkit.
Every level below is still **partial** as a curriculum: mathematical language,
worked examples, visual exploration, adaptive practice, and proof feedback are
not part of the C library.

## Level map

| Level | Existing and newly available tools | Status | Largest gaps |
|---|---|---|---|
| Primary / elementary | Dynamic exact integers and rationals with exact finite-decimal parsing (`ml_exact.h`), base-10 digit and figurate contest routines (`ml_contest.h`), length/area/mass/time/liquid/solid-volume/temperature conversions, explicit rounding, percent-of, rectangle/triangle/circle measurement (`ml_school.h`), descriptive statistics (`ml_statistics.h`, `ml_stats_extra.h`) | **Partial** | Place-value and written-operation workflows, fraction/decimal/percent teaching sequences, exact money rounding, clocks/calendars, charts, and age-appropriate problem sets |
| Middle school | GCD/LCM, modular arithmetic and primes, ratios/percent change, linear equations/inequalities, arithmetic/geometric sequences, slope and line intersection, probability distributions, polygon/triangle geometry, graph and counting tools | **Partial** | Absolute-value/rational inequalities, proportional and piecewise models, transformations, angle-chasing lessons, constructions, data plots, and explanatory feedback |
| High school | Quadratic/cubic roots and inequalities, exact Q[x], vectors/matrices and systems, coordinate and Euclidean geometry, trigonometry, sequences, compound growth, descriptive/inferential statistics, common distributions, numerical differentiation/integration | **Partial** | Conic sections, function composition/inverses/transformations, rational/exponential/log equation solvers, symbolic limits/derivatives/integrals, and proof-oriented geometry/algebra |
| Advanced secondary / contests | AMC-style recreational arithmetic, selected AIME/olympiad number theory and combinatorics, contest geometry, exact divisibility and modular routines, selected inequality and identity verifiers | **Partial** | Broad theorem/lemma library, synthetic geometry proof objects, inequalities beyond a few verifiers, functional-equation solving, extremal combinatorics, probability proofs, and a calibrated problem corpus with hint and proof feedback |

## Added pre-university API

`ml_school.h` fills several frequently used formula gaps:

- Metric and common US customary length, area, mass, time, liquid-capacity,
  and solid-volume conversions; Celsius/Fahrenheit/Kelvin conversions with
  below-absolute-zero rejection.
- Percentage calculations and percentage change.
- Exact finite-decimal-to-rational parsing (`ml_mpq_set_decimal`), decimal-place
  rounding with explicit tie rules, plus significant-digit
  rounding when the equivalent decimal-place count is within [-308,308].
- Arithmetic and geometric sequence terms and finite sums; simple and
  compound amount models.
- One-variable linear equations and linear/quadratic inequality solution sets,
  slope, and line intersection.
- Triangle, rectangle, trapezoid, circle-sector, rectangular-prism, pyramid,
  and cylinder measurement formulas.

The interface is finite-input checked and leaves outputs unchanged on failure.
Calculations use binary64 and are not exact financial accounting; exact
fraction work should use `ml_mpq_t`. US customary volume units are liquid
measures. A year in the time converter is fixed at 365 days.

## Remaining work, ordered by student value

1. Add exact money rounding policies and reusable fraction/decimal/percent
   teaching conversions; current display rounding operates on binary64 values.
2. Extend inequality sets to absolute-value, rational, and simple polynomial
   forms; improve near-degenerate quadratic classification with certified
   discriminant bounds.
3. Add conics and function transformations, inverses/composition, piecewise
   functions, and graph sampling with discontinuity handling.
4. Expand school probability and data displays, including exact finite sample
   spaces and reproducible classroom-friendly simulations.
5. Build a separate education layer: topic-indexed problems, worked examples,
   hints, misconception tags, and proof objects that can be checked.

The proof-oriented contest items require mathematical certificates, not just
numeric answer matching. They should remain explicitly partial until that
verification layer exists.
