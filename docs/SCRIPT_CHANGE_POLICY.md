# MathLib Script-Only Change Policy

All future codebase changes conducted by Blue Team should be applied through a single script whenever practical.

The script is run from the folder that CONTAINS the `v11S` working folder.

Expected layout:

    parent_folder/
        some_change_script.py
        v11S/
            src/
            include/
            tests/
            docs/
            ...

---

## Goals

1. Atomicity
   A change set is applied as one logical operation.

2. Reproducibility
   The same script can be rerun on a clean tree.

3. Auditability
   The script itself documents what changed.

4. Safety
   Scripts should avoid silent manual edits.

5. Orchestrator control
   The human orchestrator remains in control of every change batch.

---

## Script Naming Convention

Scripts should be numbered:

    00_closure_bootstrap.py
    01_p0_exp_log_pow_complex.py
    02_p1_simd_index_ci_docs.py
    03_...

The number defines intended order.

---

## Script Requirements

Every script should:

- verify that it is run from the correct parent directory,
- verify that `v11S/` exists,
- print what it writes or modifies,
- avoid interactive prompts,
- be idempotent where possible,
- use `--force` if overwrite is desired,
- preserve LF line endings,
- avoid introducing fast-math, heap allocation, or global state.

---

## Source Modification Rules

Scripts that modify source code should:

- prefer small, explicit replacements,
- fail loudly if expected text is missing,
- not apply partial changes if a required anchor is absent,
- update tests and docs in the same change set when behavior changes,
- preserve public API signatures unless explicitly authorized.

---

## Forbidden During A3

Unless explicitly authorized:

- no new modules,
- no new public APIs,
- no new math families,
- no performance experiments,
- no speculative refactors,
- no feature creep.

A3 is for closure.
