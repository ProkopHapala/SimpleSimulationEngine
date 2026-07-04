---
name: doc-read-navigate
description: Before writing new code — search existing implementations, topical audits, READMEs to avoid duplication
---

## Before You Write — Check What Exists

1. **Search topical audits**: `doc/TopicalAudit/*.md` — one file per scientific topic, lists all implementations across languages with status and parity
2. **Search CODEMAP.md** (if exists): file locations and module relationships
3. **Grep codebase**: search for function/class names matching your planned implementation
4. **Read README.md** in the target folder — may describe what's already there

## Decision

- **Found exact match**: reuse it. Don't duplicate.
- **Found similar**: generalize it instead of writing new. If generalization risks breaking existing code: **stop and report for approval**.
- **Found nothing**: proceed, but remember to document the new topic post-implementation (see `doc-task-summary` skill).

## Separation of Concerns Check

Before writing, verify your plan doesn't mix:
- Compute logic with plotting/diagnostics
- Backend with GUI/CLI
- Test scripts reimplementing shared module functions

If you need plotting/debugging: check for shared utilities first (e.g., `pyBall/plot_utils.py`, `TestUtils.py`).

## Where to Find Things

- `doc/TopicalAudit/` — cross-language implementation maps per scientific topic
- `doc/AGENTs/skills/` — task-specific skills (debugging, OpenCL, parity, etc.)
- `doc/AGENTs/protocols/` — domain-specific protocols (forcefields, topology, QM)
- `README.md` in any folder — local index of contents
