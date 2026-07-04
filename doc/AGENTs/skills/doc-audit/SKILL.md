---
name: doc-audit
description: Dedicated documentation work — OKF format, topical audits, extract-overview and inline-doc workflows
---

## Documentation System

Two layers, built bottom-up:

### README.md per folder
- No frontmatter — it's an index, not a concept
- 1-3 sentences: what this folder contains
- Bullet list of key files/subfolders with one-line descriptions
- Build from folders you're working in, not top-down

### Topical Audit (`doc/TopicalAudit/`)
Cross-language maps connecting implementations of the same scientific topic. One file per topic.

```yaml
---
type: TopicalAudit
title: <Topic Name>
tags: [topic, cross-language]
---
```

Body:
- **Summary**: 1 paragraph
- **Implementations**: table `Language | Location | Status | Notes`
- **Parity Status**: verified pairs, tolerance, test reference
- **Open Issues**: TODOs, unported features

Status values: `active`, `experimental`, `deprecated`, `unfinished`.

## OKF Principles
- Markdown + YAML frontmatter, no tooling required
- Frontmatter: `type` (required), `title`, `description`, `tags`, `timestamp` (optional)
- Cross-link with markdown links (absolute `/path/to.md` preferred)
- Structural markdown (headings, lists, tables, code blocks) over freeform prose

## Extract Documentation Workflow
Generate concise human-readable overview from source files:
1. **Scope**: one cohesive unit (single file or file+header pair)
2. **Read**: skim imports, structure, function bodies for algorithmic essence
3. **Compose**: Background (3-8 sentences), File/Class Summary, Functions (one-liner per function — purpose only, no args/returns), Tutorial (CLI flags / GUI controls / library snippet), Notes & Pitfalls
4. **Discover usage**: detect argparse/main for CLI, check `*.sh` for run commands, find event bindings for GUI
5. **Validate**: run with small inputs if feasible; avoid full API dumps

## Inline Documentation Workflow
Insert concise standardized comments into code:
- **C/C++/OpenCL**: `///` one-liners (Doxygen); avoid `/* */`
- **Python**: single-line docstrings `"""One-line summary"""`
- File header: 1-3 sentences background + AUTO-DOC markers
- Functions: one-line `@brief` above definition (C++) or docstring inside (Python)
- **AUTO-DOC markers**: `/// === AUTO-DOC BEGIN ===` / `/// === AUTO-DOC END ===` (C++), `# === AUTO-DOC BEGIN ===` / `# === AUTO-DOC END ===` (Python) — allows safe regeneration
- Never delete existing comments; keep both if conflicting, mark new as AUTO-DOC
- No parameter/return enumeration unless critical; focus on purpose and essence
- Minimal edits — don't refactor while documenting

## Related Skills
- `doc-read-navigate` — check existing implementations before writing
- `doc-task-summary` — quick doc updates after implementation
