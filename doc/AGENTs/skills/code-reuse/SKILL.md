---
name: code-reuse
description: Writing new code — inventory-first, module-vs-script placement, no-new-files, STOP triggers for duplication
trigger:
  glob:
    - "**/*.py"
    - "**/*.js"
    - "**/*.mjs"
    - "**/*.cpp"
    - "**/*.h"
    - "**/*.cl"
    - "**/scripts/**"
    - "**/test*/**"
---

## Core Rule: No Duplication, No Ad-Hoc Functions, No New Files Unless Explicitly Asked

Before writing ANY new function, search the codebase for existing implementations. This is a blocking step — do not skip it. See also `doc-read-navigate` skill for where to search.

## Procedure

### Step 1: Inventory Before Writing (BLOCKING)
1. Search the codebase for existing functions that do the same or similar thing (grep_search, code_search, find_by_name)
2. If you find a match: import and use it. Do not write a new one.
3. If you find a near-match: generalize the existing function (with approval if risky)
4. If NO match exists: write the function in the appropriate **existing** shared module, not in a script

### Step 2: Module vs Script Placement
- **Shared modules** (e.g. `pyBall/`, `cpp/common/`, `web/common_js/`): Export reusable functions. This is the ONLY place reusable logic lives.
- **Scripts** (e.g. `tests/*.py`, `scripts/*.mjs`, `*.sh`): Thin wrappers that import from modules and orchestrate. May contain only test-specific glue (logging, file I/O, test case wiring).
- **Hard rule**: If a function in a script is not test-specific glue, it MUST be moved to a shared module and imported. No exceptions.

### Step 3: No New Files Unless Explicitly Requested
- **Do NOT create new modules or scripts** unless the user explicitly asks for them.
- **Do NOT create new helper files**, utility files, or wrapper scripts without explicit user approval.
- When a new function is needed, place it in the **most suitable existing module** — do not create a new file for it.
- If you are unsure which module to use, search for related functionality and choose the module where similar functions already live.
- **Consolidation principle**: When refactoring, prefer merging multiple scripts into one unified entry point rather than creating new per-feature scripts. Mark old scripts as deprecated — do not delete them.

### Step 4: Pre-Commit Checklist
- [ ] Searched for existing implementations before writing each function?
- [ ] All reusable functions are in shared modules (not in script files)?
- [ ] Test/script files contain only glue + imports?
- [ ] No copy-pasted logic from elsewhere?
- [ ] No new files created unless explicitly requested?

### Step 5: STOP Triggers
- About to write a function and think "this might exist already" → STOP and search
- About to export something from a script file → STOP, move it to a module
- About to copy-paste code → STOP, refactor into a shared function
- About to create a new file (module, script, helper) → STOP, ask: "Can this go in an existing module?" If yes, use the existing module. If no, ask the user for permission.
- About to create multiple scripts for related functionality → STOP, consider a single unified CLI with subcommands instead
