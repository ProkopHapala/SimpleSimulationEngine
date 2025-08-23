# Workflow: Extract Documentation (Language‑agnostic)

## Purpose
Generate a concise, didactic, human‑readable document from source files (`.py`, `.h/.cpp`, `.cl`, …). Cover scientific background, summarize API (one‑liners per function), and provide a short tutorial (CLI/GUI/library usage).

## Inheritance
- Inherits `doc/CodingRules/global_rules.md`
- Follows `doc/CodingRules/workflow_template.md`

## Workflow Steps
1. Scoping (target selection)
   - Identify the target file(s) and, if applicable, the main class/module/kernel to document.
   - Prefer documenting one cohesive unit at a time (single file or file+header pair).

2. Read & understand code
   - Skim top comments, imports/includes, and file names for domain context (physics/chemistry/geometry/graphics).
   - Parse structure:
     - Python: module‑level functions/classes; look for `if __name__ == "__main__":` and `argparse`.
     - C/C++: declarations in `.h`, definitions in `.cpp`; scan `main()` for CLI; note build/run scripts in `tests_bash/`.
     - OpenCL: `__kernel` entry points and argument roles; locate host code bindings (e.g., Python `OpenCLBase`).
   - For each function/class/kernel: read the body to infer the purpose and the essence of the algorithm; ignore signature details.

3. Compose the document (output skeleton)
   - Title: `<File or Module Name>: Overview`
   - Overview & Background
     - 3–8 concise sentences explaining what the code implements and why (physics/mathematics/chemistry context, key equations/ideas in words; keep symbols brief).
   - File/Class Summary
     - Short paragraph describing the file’s role and any key classes/structs.
   - Functions / Kernels (API summary)
     - Bullet list, one line per item: `name` — 1–2 sentences on purpose and algorithmic essence (no args/returns).
       - ommit argumnet list or return values
   - Tutorial / How to Use
     - If CLI/script: how to run; list CLI flags (bullets) with short explanations.
     - If GUI: list main widgets/buttons and keyboard/mouse controls (bullets) and what they do.
     - If library: minimal code snippet to initialize and call core entry points; how to extract results.
   - Notes & Pitfalls
     - Edge cases, assumptions, data ranges, performance hints.

4. Discovering usage details
   - CLI: detect `argparse` (Python) or `int main(int argc, char** argv)` parsing (C++); list flags and defaults.
   - Run scripts: check `tests_bash/` or local `*.sh` to infer build/run commands; include the exact command if stable.
   - GUI controls: search for event bindings (e.g., key/mouse callbacks) and summarize.

5. Validate & finalize
   - Optional: run the code with small inputs to confirm CLI defaults and produce a minimal example in the tutorial.
   - Avoid long tables or full API dumps.

## Output format (example skeleton)
```markdown
# <Module/File>: Overview

## Background
<Short didactic explanation of model/physics/math>

## File/Class Summary
<Brief role of the file and key classes>

## Functions
- `foo` — Computes X using Y heuristic; used to update Z per step.
- `bar` — Builds neighbor list via grid hashing; O(N) expected.
- `step` — Advances state by one time step using symplectic update.

## Tutorial
### Run (CLI)
<command>
- --dt: time step
- --n: number of particles

### GUI Controls
- LMB: select
- RMB drag: rotate view

### Library Example (Python)
```python
from pkg import mod
obj = mod.Solver(params)
obj.init(state)
obj.run(n=100)
res = obj.result()
```

## Notes
- Stable for dt < 0.1; uses float64.
- Input assumes SI units.
```

## Language‑specific notes
- Python
  - Prefer module doc extraction from `__main__` and `argparse` for CLI. Note `np.set_printoptions` and plotting flags.
- C/C++
  - Reflect build/run via repo bash scripts (e.g., `tests_bash/.../*.sh` or `C/build.sh`). Document `-g/-O*` flags only if relevant.
- OpenCL
  - List kernels by name; summarize what each kernel computes and buffers it reads/writes. Mention host entry point(s).

## Pitfalls
- Do not dump signatures or all parameters; focus on purpose and essence.
- Keep bullets one‑liners; avoid multi‑line prose per function.
- Prefer domain explanation over implementation minutiae.
