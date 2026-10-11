# MathLib v12A1 Development Scripts (HISTORICAL)

This directory contains the surviving controlled change sets once applied
to the v12A1 tree (now v12R2).

R5 despot: `00_v12a1_bootstrap.py` is referenced below but missing on disk,
and numbering has gaps (no 34-36, 42-43, 51). The surviving `30.py`..`52.py`
are kept for provenance; do not treat the sequence as a runnable pipeline.
Current gates are `closure_gate.sh` + `run_all_tests.py`.

Scripts are run from the folder that CONTAINS the tree (now repo root, not `v12A1`).

Example (historical):

    python3 30.py

---

## Script Sequence

See `docs/V12A1_ROADMAP.md` for the full table.

Scripts are numbered in dependency order.
Scripts 07, 08, 10, 11, 12, 14 can run in parallel after bootstrap.

---

## Rule

Do not manually drift the tree when a script can express the change atomically.
