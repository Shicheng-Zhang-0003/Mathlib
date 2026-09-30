# v12A1 A1 Feature Freeze — SUPERSEDED by v12R2 / V1.2-RC2

<!-- MATHLIB_V12A1_A1_FREEZE -->
<!-- DESPOT-SUPERSEDE 2026-09-30: A1 freeze lifted. v12R2 shipped 12 Batch-1
     modules + kelvin (34 TUs) plus despot audit fixes. This document is
     historical. Canonical scope is now CMake/Makefile/run_all_tests.py/
     run_edge_tests.sh (34 TUs) + docs/DESPOT_AUDIT.md. -->

Subsection: 1.1
Effective: 2026-08-05

v12A1 is now in **A1 closure freeze**.

## Allowed during A1 closure

Only the following are permitted:

1. Correctness fixes required by the A1 closure table.
2. Test and oracle expansion needed to prove those fixes.
3. Validation and gate work.
4. Documentation alignment.
5. Script/process hygiene directly tied to A1 closure.

## Not allowed during A1 closure

The following are frozen:

1. New modules.
2. New public APIs.
3. New math families.
4. New performance experiments.
5. Speculative refactors.
6. Feature creep of any kind.
7. Anything whose primary purpose is to make A2 easier.

## Script rule

Every A1 closure change must be applied through a numbered script in:

```text
v12A1/scripts/v12a1/
```

and that script must correspond to a specific A1 closure subsection.

## Closure mindset

The goal is no longer expansion.

The goal is:

> make the current tree true, tested, documented, and stable enough to close A1.
