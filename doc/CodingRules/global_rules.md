# Global Coding Rules

Applies to all languages (Python, C/C++, OpenCL) and all workflows. These are the non‑negotiable baseline rules for this repository.

## Core principles
- Minimal dependencies: stdlib + NumPy/Matplotlib unless explicitly requested
- Concise, modular code: small, testable, step‑by‑step changes; divide‑and‑conquer
- Reuse existing functions; avoid duplication (check first, report if adaptation needed)
- Pure, data‑oriented functions with explicit inputs/outputs; default named args to avoid long call strings
- Fail loudly: no silent handling; assertions for invariants; crashes with stack trace preferred to masking
- Comment out deprecated/experimental code instead of deleting; mark with TODO/DEBUG

## NEVER DO THIS
- Never delete, override, or rearrange existing code without explicit permission in the user prompt
- Never perform random aesthetic/style edits unrelated to the task
- Never apply "quick‑fixes" that hide root causes (e.g., hard‑coded outputs)

## Debugging first
- Debuggability > UX; do not hide issues
- Initially add debug prints for key variables and flow; remove later after debugging is done
- Avoid broad try/except as they mask bugs; prefer loud crashes which write stack trace
- Make small, testable changes and run after every change
- Log function entries and key conditions when helpful to track the flow
- Mark unfinished/experimental code clearly (e.g., # TODO, # DEBUG)

## Performance (general)
- Prefer data‑oriented code that is cache‑friendly and avoids overheads
- Preallocate and reuse buffers; avoid repeated allocation in hot paths
- Be explicit about dtypes/shapes; prefer contiguous memory where possible

## Style (cross‑language)
- Prefer concise/compact code; avoid bloated structures and unnecessary empty lines
- Short variable names OK (math/physics symbols like E, T, m) when locally clear
- Prefer one‑liner expressions, assume unlimited line-width
- Avoid line wrapping that hurts readability of expressions (assume wide editor)
- Inline comments behind the code line for rationale; reserve header-comments for modules/classes/functions
- Doxygen: use `///`; avoid `/* ... */`
- C++: use `printf` for debugging over `std::cout`; prefer plain C arrays (e.g., double*) in hot paths
- Vector math in C/C++: use `Vec3.h`, `Vec2.h` with `Vec3d`, `Vec2d`, and helpers like dot(), norm(), cross()

## Visualization
- Separate compute vs plotting; no plotting in core algorithms
- Plotting optional via flags (e.g., --noPlot, --saveFig)
- `plt.show()` only in CLI/main, never in library code
- Prefer shared plotting helpers (e.g., plot_utils.py) to avoid duplication
