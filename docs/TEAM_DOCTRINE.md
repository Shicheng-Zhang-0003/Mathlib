# MathLib Team Doctrine

This document defines how the multi-model development workflow operates.

The goal is not merely to generate code.

The goal is to generate, verify, understand, and own the system.

---

## Release Philosophy

All releases follow:

    S -> A1 -> A2 -> A3 -> S

Meaning:

- `S` = stable release
- `A1` = first major development wave
- `A2` = second major development wave
- `A3` = freeze, audit, bug discovery, hardening
- next `S` = promoted stable release

Current state:

> v11S is operationally between late A3 and true S.

Therefore:

- no new features,
- no new math families,
- no new public APIs,
- only correctness, tests, documentation, and closure hygiene.

---

## Blue Team

Owner: Qwen

Purpose:

- large-scale implementation,
- refactoring,
- test scaffolding,
- documentation generation,
- script-based codebase changes,
- turning contracts into code.

Rules:

- Blue implements against contracts.
- Blue does not silently invent new contracts.
- Blue preserves zero allocation, thread safety, and C99 compatibility.
- Blue changes are applied through single scripts whenever possible.

---

## Red Team

Members:

- ChatGPT
- Gemini
- GLM
- Kimi, when available

Purpose:

- mathematical correctness review,
- IEEE-754 special-case review,
- standards review,
- test-gap analysis,
- documentation/code mismatch detection.

Rules:

- Red does not generate large implementation patches unless explicitly asked.
- Red produces defect lists, severity rankings, and missing tests.
- Red disagreements become test cases.
- Red is invoked after pseudo-blue-red and before final release promotion.

---

## Yellow Team

Owner: Microsoft Copilot

Purpose:

- Windows portability,
- MSVC / clang-cl compatibility,
- CMake presets,
- install profiles,
- packaging,
- Windows CI,
- DLL / static library validation.

Rules:

- Yellow does not change numerical behavior.
- Yellow preserves the deterministic scientific contract.
- Yellow must validate behavior, not merely compilation.

---

## Violet Team

Owner: Grok

Purpose:

- extreme adversarial review,
- hostile input discovery,
- UB traps,
- security and robustness edge cases,
- fuzzer blind-spot discovery.

Rules:

- Violet does not propose feature creep during A3.
- Violet produces attack vectors, not roadmaps.
- Violet findings become tests or documented limitations.

---

## GitHub Copilot

Optional role:

- holistic reviewer,
- documentation drift checker,
- small regression-test writer,
- learning companion.

It is not the final mathematical authority.

---

## Decision Rules

1. Tests decide correctness.
2. If Red models disagree, convert the disagreement into a test.
3. If behavior cannot be fixed now, document it as a known limitation.
4. Documentation must never overclaim implementation.
5. During A3, correctness beats novelty.
6. A module is mathematically owned only when the orchestrator can explain its contract, failure modes, and tests.

---

## Handoff Packet

Before real Red Team review, prepare:

1. changed files,
2. contract decisions,
3. test evidence,
4. sanitizer output,
5. known limitations,
6. deferred items.

No tree should be handed to Red Team without these.
