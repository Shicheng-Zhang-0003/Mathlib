# MathLib v12A1 Development Scripts

This directory contains scripts that apply controlled change sets to the v12A1 tree.

Scripts are run from the folder that CONTAINS `v12A1`.

Example:

    python3 00_v12a1_bootstrap.py

---

## Script Sequence

See `docs/V12A1_ROADMAP.md` for the full table.

Scripts are numbered in dependency order.
Scripts 07, 08, 10, 11, 12, 14 can run in parallel after bootstrap.

---

## Rule

Do not manually drift the tree when a script can express the change atomically.
