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

## Writing Style: Essence Over Obvious

The goal of documentation is to convey what a reader **cannot easily see** from the code.
Stating the obvious is lazy — documentation should contain high-context background,
motivation, design rationale, and highlights that capture the **essence** of the file.

### What to Write

- **Essence**: What is this file *for*, at the highest level? One sentence that a new
  contributor can understand without reading the code. Example: "Fast function
  approximations for hot loops in physics simulations" — not "Math utility functions."
- **Motivation**: Why does this file exist? What problem does it solve that the
  alternatives (std::unordered_map, naive approach, library X) don't? Example:
  "Open addressing with linear probing — unlike std::unordered_map's linked lists,
  this keeps probe sequences cache-friendly for spatial hash tables."
- **Design decisions and trade-offs**: What were the choices, and why this one?
  Example: "Power-of-2 batch size so division becomes bit-shift — faster on
  architectures without hardware divide."
- **Non-obvious algorith, desing decision and principle**: Things that would surprise someone reading
  the code. Example: "hits[] array tracks collision count per slot so getAllInBox()
  can stop early — this is the key optimization over naive open addressing."
- **Use cases**: Where is this actually used in the engine? Example: "Used as
  edgesOfVerts in Builder2 for O(1) 'which edges touch vertex i?' queries."
- **Background/context**: Algorithmic or mathematical context that informs the design.
  Example: "Same principle as NEB/string methods in molecular dynamics, applied here
  to geometric implicit surfaces."

### What NOT to Write

- **Obvious restatements of code**: "This file defines a HashMap class" — the reader
  can see that. Instead: "Open-addressing hash map optimized for integer box-index keys
  in spatial bucketing."
- **Parameter/return type listings**: Unless critical for understanding, skip them.
  The code already shows types. Focus on *purpose*, not *signature*.
- **Empty formalities**: "This file provides various utility functions for..." — say
  what the utilities are *for* and *why they're here*, not that they exist.
- **Fluff and filler**: "This powerful and flexible module..." — adjectives are noise.
  Describe what it does, not how impressive it is.
- **Tutorial-style walkthroughs**: The header doc is not a tutorial. It's context for
  someone who is about to read or modify the code.

### Format Rules

- Use `///` Doxygen style for C/C++ headers (not `/* */`)
- Start with `/// @file filename.h` then `/// @brief <one-sentence essence>`
- Follow with `///` paragraphs for: motivation, design decisions, key functions/structs,
  use cases, background — only the ones relevant to this file
- Use `**bold**` for key terms and struct/function names within the doc
- Use `-` bullet lists for enumerating design choices or key functions
- Keep total header doc to 10-25 lines — dense, not verbose
- Every sentence must carry information that is NOT visible in the code

### README.md Style

- One bullet per file in the folder
- Each bullet: **filename.h** — one-sentence essence, then comma-separated key details
- No frontmatter, no headings beyond the folder name
- Order: most important/foundational files first, implementation files last
- Include .cpp files — they get a one-liner saying what they implement

### Examples

Good `@brief`:
> Open-addressing hash map for spatial bucketing — optimized for box-indexed lookups.

Bad `@brief`:
> A hash map data structure.

Good body excerpt:
> Unlike std::unordered_map (separate chaining with linked lists), this uses open
> addressing with linear probing. The key insight: this map is designed for spatial
> hash tables where the key is a cell/box index (integer), not an arbitrary object.

Bad body excerpt:
> This class provides methods for inserting, removing, and finding elements in a
> hash table. It supports generic types through templates.

## Related Skills
- `doc-read-navigate` — check existing implementations before writing
- `doc-task-summary` — quick doc updates after implementation
