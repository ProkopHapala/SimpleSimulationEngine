---
name: doc-task-summary
description: After implementing — update README.md and topical audits, check for duplication
---

## After You Implement — Update Docs

1. **README.md** in the folder you touched:
   - If missing: create one (1-3 sentences + bullet list of key files, see `doc-audit` skill for format)
   - If exists: add/update entries for new or changed files
   - Keep it a quick index, not a manual

2. **Topical audit** in `doc/TopicalAudit/`:
   - If your implementation is a new take on an existing topic: add a row to the implementations table
   - If topic file doesn't exist yet: create one (see `doc-audit` skill for format)
   - Mark old implementations as `deprecated` if superseded

3. **Duplication check**:
   - Search for similar logic across languages — if found, note in topical audit
   - Don't consolidate now unless trivial — just record it

## Don't Over-Do It

- Small edits: just update the README.md bullet for the file you changed
- New feature: add topical audit entry + README.md update
- Don't write full documentation unless explicitly asked — that's the `doc-audit` skill's job
