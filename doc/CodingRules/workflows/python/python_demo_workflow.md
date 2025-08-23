# Workflow: Python Demo

## Purpose
Quick sketch of a scientific/numerical algorithm (physics, chemistry, geometry, graphics) using NumPy/Matplotlib. Optimized for clarity, minimal code, and simple visualization.

## Inheritance
- Inherits `doc/CodingRules/global_rules.md`
- Complements `doc/CodingRules/Python_Numerics_Rules.md`

## Workflow Steps
1. Setup
   - Imports: `import numpy as np`, `import matplotlib.pyplot as plt`
   - Printing: `np.set_printoptions(linewidth=np.inf)` (wide output)
   - Dependencies: keep to NumPy/Matplotlib unless explicitly needed
   - Imports at module top; avoid importing inside functions
   - CLI entry guarded by `if __name__ == "__main__":`; expose `--verbosity` to control diagnostics
   - CLI: add `argparse` only if parameters exceed a few toggles
2. Coding strategy
   - Small, pure functions; default named args to shorten call sites
   - Keep core computation separate from plotting
   - Prefer vectorization and preallocation; avoid Python loops in hot paths
   - Reuse existing helpers/utilities; avoid duplication (cite or wrap if extension is needed)
   - Prefer tuples/data-classes over dictionaries when feasible (lighter, clearer APIs); never use a Python dict as a replacement for a C struct (use dataclasses or NumPy structured dtypes when fixed layout is needed)
3. Debugging strategy
   - Add end-of-run sanity checks (e.g., conservation laws)
   - Gate prints by `verbosity`; print shapes/min/max to catch NaN/Inf
   - Fail loudly: use assertions for invariants; avoid broad try/except
   - Avoid debug prints inside hot loops; keep them in setup/teardown or gated blocks
   - Do not store debug-only data in core objects by default
   - Prefer simple prints over heavy logging frameworks
   - Use helper functions for multi-line diagnostics (keep hot path clean)
4. Testing & validation
   - Validate against analytical cases or downscaled inputs
   - Add quick asserts for array shapes and invariants
   - Check analytical vs numerical derivatives where applicable (see `python/pyMolecular/plotUtils.py`: `numDeriv`, `plot1d(bNumDeriv=True)`, `plot_funcs`)
5. Performance considerations
   - Avoid temporary copies; prefer views and contiguous arrays
   - Compute only what is needed for the taken branch
   - Prefer NumPy arrays over Python lists in hot paths; preallocate shapes
   - Use boolean masks and advanced indexing instead of branchy loops
   - For growing data structures, use list/dict/set comprehensions
   - Be explicit about dtype and memory layout; default `float64`; for PyOpenCL prefer `float32`
6. Visualization/reporting
   - Minimal plotting; optional via flags (`--noPlot`, `--saveFig`)
   - Reuse `python/pyMolecular/plotUtils.py` helpers: `plot1d`, `plot_func`, `plot_funcs` (and molecule plotting) instead of re‑implementing
   - Call `plt.show()` only in `__main__`

### Minimal example: derivative check

```python
import numpy as np
from pyMolecular.plotUtils import plot1d  # repo helper

x = np.linspace(-2.0, 2.0, 2001)
k = 3.0
E = 0.5*k*x**2            # energy
F = -k*x                  # analytical force = -dE/dx

fig, _ = plot1d(x, [E], derivs=[F], labels=['harmonic'], bNumDeriv=True)
# Numerical derivative (force) overlay is drawn by plot1d via numDeriv()
# Defer plt.show() to CLI entrypoint
```

### Minimal example: optional plotting guard

```python
plot, save = True, False   # or from argparse
if plot:
    fig, _ = plot1d(x, [E], derivs=[F], labels=['harmonic'], bNumDeriv=True)
    if save: fig.savefig('out.png', dpi=150)
    # call plt.show() only in __main__
```

## Notes
- Keep demo code compact (~50–100 LOC where possible)
- Prefer inline comments after code; avoid verbose headers
- Short, clear variable names (physics/math symbols welcome, e.g., `E`, `F`, `T`) and compact expressions
- See also `doc/CodingRules/Python_Numerics_Rules.md` for extended guidance
- For instrumentation patterns and regression snippets see `doc/CodingRules/workflows/python/python_debugging_workflow.md`
