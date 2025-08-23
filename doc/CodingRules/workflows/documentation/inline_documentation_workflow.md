# Workflow: Inline Documentation (Language‑agnostic)

## Purpose
Insert concise, standardized inline documentation into code using per‑language conventions. Use the extracted doc as source. Keep comments brief (one‑liners), focusing on purpose and algorithmic essence, not exhaustive param/return listings.

## Inheritance
- Inherits `doc/CodingRules/global_rules.md`
- Follows `doc/CodingRules/workflow_template.md`

## Workflow Steps
1. Inputs & scope
   - Start from the freshly created documentation (see `extract_documentation_workflow.md`).
   - Choose a small, cohesive target (one module/file/class) per pass to keep diffs minimal.

2. Comment style by language (standard)
   - C/C++/OpenCL: Doxygen one‑liners using `///` per global/file/class/function; avoid `/* ... */`.
   - Python: single‑line docstrings (`"""One‑line summary"""`) at module, class, and function start; keep to one line.

3. What to add (and where)
   - File/module header (top of file)
     - Brief background (1–3 sentences) explaining purpose and domain context.
     - Mark with update tag for future regeneration.
   - Class/struct/kernels
     - One‑line summary above the declaration/definition.
   - Functions
     - One‑line `@brief` style summary above the definition (C/C++/OpenCL) or as a one‑line docstring inside the function (Python).
   - Optional anchors
     - Surround auto‑generated headers with markers for safe updates:
       - C/C++/OpenCL: `/// === AUTO‑DOC BEGIN ===` / `/// === AUTO‑DOC END ===`
       - Python: `# === AUTO‑DOC BEGIN ===` / `# === AUTO‑DOC END ===`

4. Application rules
   - Never delete existing comments; if conflicting, keep both, but mark the new one as AUTO‑DOC.
   - prefer one‑liners as per repo style.
   - Do not enumerate parameters/returns unless critical; focus on purpose and essence of the algorithm.
   - Minimal edits: avoid reflowing code

5. Validate
   - Build/run quickly to ensure no syntax errors:
     - C/C++: use repo bash scripts (e.g., `C/build.sh`, `tests_bash/.../*.sh`).
     - Python: import module and run a tiny call path if feasible.
     - OpenCL: ensure comments don’t break kernel pragmas; rebuild host once.

6. Record & update strategy
   - Include timestamp and tool tag in the file/module header, e.g., `/// AUTO‑DOC (2025‑08‑23) by agent`.
   - Future runs should replace only within the AUTO‑DOC markers; outside content remains intact.

## Examples
- C++ function
```cpp
/// Computes pairwise forces using cell lists; O(N) expected per step.
void computeForces(const double* x, const double* y, double* fx, double* fy, int n);
```

- C++ file header
```cpp
/// === AUTO-DOC BEGIN ===
/// @file md.cpp — Minimal MD integrator (Verlet); SI units; educational reference.
/// Implements pair forces, integration, and simple thermostat.
/// === AUTO-DOC END ===
```

- Python function
```python
def step(x, v, dt):
    """Advance state by one step using symplectic update."""
    ...
```

- Python module header
```python
# === AUTO-DOC BEGIN ===
"""Minimal MD demo: integrates N particles with Lennard-Jones; plots energy and trajectory."""
# === AUTO-DOC END ===
```

- OpenCL kernel
```c
/// Updates velocities using accumulated forces; dt in seconds; float32.
__kernel void update_vel(__global float2* v, __global const float2* f, const float dt) { ... }
```

## Language‑specific notes
- Python
  - Prefer module/class/function one‑liners; avoid long doctrings; no parameter matrices unless crucial.
- C/C++
  - Use `///` everywhere for Doxygen; avoid block comments.
  - Place `@file` header at top; per‑function `///` brief above each definition.
- OpenCL
  - Use `///` safe comments before `__kernel` definitions; avoid interfering with pragmas.

## Pitfalls
- Do not refactor while documenting; keep diffs tiny.
- Avoid repeating the same long background in many files—file header should be succinct.
- Ensure markers are consistent to allow safe regeneration.
