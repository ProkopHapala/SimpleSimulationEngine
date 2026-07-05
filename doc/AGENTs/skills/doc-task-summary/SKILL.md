---
name: doc-task-summary
description: After implementing — update file headers, README, topical audits; note duplication
---

## After You Implement — Update Docs

Same session as the code, while context is fresh. Style rules: `doc-audit` skill.

### 1. File / module headers (changed files only)

At top of each new or materially changed source file:

- **Essence** — one sentence: what problem this file solves (not a file list).
- **Design** — only non-obvious choices (SSOT, data flow, why this approach).
- **Open issues / caveats** — bullets for landmines: partial/deferred state, ordering deps, duplicate parallel implementations, strict/fail-loud APIs, import/export limits, out-of-scope TODOs.
- **Functions** — one-line purpose comment only where the name is not enough (language style per `doc-audit`: `///` etc.).

Do not delete existing comments. Do not document every parameter. Do not refactor while documenting.

### 2. README.md in the folder you touched

- If missing: create one (1–3 sentences + bullet list of key files, see `doc-audit` skill for format)
- If exists: add/update entries for new or changed files
- Keep it a quick index, not a manual

### 3. Topical audit (`doc/TopicalAudit/`)

- New take on an existing topic: add a row to the implementations table
- Topic file doesn't exist yet: create one (see `doc-audit` skill for format)
- Mark old implementations as `deprecated` if superseded
- Cross-cutting caveats from headers → audit **Open Issues** when relevant

### 4. Duplication check

- Search for similar logic across languages — if found, note in topical audit or file caveats
- Don't consolidate now unless trivial — just record it

## Don't Over-Do It

- **Small edit**: caveat in header if behavior changed; README bullet if user-facing
- **New feature**: headers + caveats on new files; README update; audit row if new topic
- **Full documentation**: only when explicitly asked — use `doc-audit`

## Related Skills

- `doc-audit` — format and inline-comment conventions per language
- `doc-read-navigate` — read existing headers (caveats) before extending code
