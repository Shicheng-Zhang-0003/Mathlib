# MathLib Closure Scripts

This directory contains scripts that apply controlled change sets to the v11S tree.

Scripts are run from the folder that CONTAINS `v11S`.

Example:

    python3 00_closure_bootstrap.py

---

## Current Script Sequence

### 00_closure_bootstrap.py

Creates closure control documents and strict closure gate.

Does not modify math source.

Creates:

    v11S/docs/CLOSURE_PUNCHLIST.md
    v11S/docs/TEAM_DOCTRINE.md
    v11S/docs/SCRIPT_CHANGE_POLICY.md
    v11S/scripts/closure/README.md
    v11S/closure_gate.sh

---

### 01_p0_exp_log_pow_complex.py

Intended next script.

Purpose:

- fix P0 exp/log/pow/complex blockers.

Likely targets:

    v11S/src/internal/pow_util.h
    v11S/src/exp_log.c
    v11S/src/complex.c

---

### 02_p1_simd_index_ci_docs.py

Intended follow-up script.

Purpose:

- SIMD semantic guarding,
- internal size arithmetic,
- CI determinism,
- documentation alignment.

---

## Rule

Do not manually drift the tree when a script can express the change atomically.
