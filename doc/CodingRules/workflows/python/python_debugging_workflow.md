# Workflow: Python Debugging

## Purpose
Diagnose and fix incorrect behavior in NumPy/Matplotlib numerical code. Emphasizes visibility and reproducibility.

## Inheritance
- Inherits `doc/CodingRules/global_rules.md`
- Complements `doc/CodingRules/Python_Numerics_Rules.md`

## Workflow Steps
1. Setup
   - Add `--verbosity` flag; default quiet
   - Enable wide prints: `np.set_printoptions(linewidth=np.inf)`
   - Keep plotting optional via flags; do not mix with core logic
2. Coding strategy
   - Instrument function entries/exits; print key shapes/min/max/nan counts
   - Keep debug prints one‑liners; if multi‑line, wrap in a small helper
   - Comment out old code rather than delete; tag with `# TODO` / `# DEBUG`
3. Debugging strategy
   - Reproduce with smallest input; add assertions for invariants
   - Bisect the pipeline to localize fault; print intermediate arrays (first few elems)
   - Prefer deterministic seeds for reproducibility
4. Testing & validation
   - Add quick unit‑like checks; compare to analytical/simple cases
   - After fix, preserve a regression snippet (tiny input + expected metric)
5. Performance considerations
   - Favor clarity while debugging; loops acceptable temporarily
   - Avoid excessive prints in hot loops; sample or summarize
6. Visualization/reporting
   - Minimal diagnostic plots behind flags; mark NaN/Inf clearly
   - `plt.show()` only in CLI

## Common guidance (applies to any Python program)
- Imports at module top; avoid importing inside functions
- CLI main guard: `if __name__ == "__main__":`
- Minimal dependencies (NumPy/Matplotlib); compute separated from plotting
- Plotting optional via flags (`--noPlot`, `--saveFig`); call `plt.show()` only in CLI
- Prefer tuples/dataclasses over dictionaries for small records; never use a Python dict as a replacement for a C struct
- Prefer vectorization/preallocation; be explicit about dtype and shapes (CPU default `float64`, PyOpenCL prefers `float32`)
- Avoid Python loops and if/then/else branches in hot paths; use NumPy array ops and boolean masks. If container growth is needed, use list/dict/set comprehensions outside hot loops

## Debug-specific techniques
- Instrumentation one-liners: print shapes/min/max/NaN counts; wrap multi-line in a tiny helper
- Bisect pipeline: validate outputs after each stage on tiny inputs; assert invariants
- Reproduce deterministically: fixed seeds; store a tiny regression case after fixing
- Avoid broad try/except; let stack traces surface; comment out old code (keep context)
- Gate prints using a `verbosity` integer; avoid prints inside hot loops or sample sparsely

### Minimal helpers (concise)

```python
import numpy as np

def stats(name, a):
    print(name, a.shape, np.nanmin(a), np.nanmax(a), np.count_nonzero(~np.isfinite(a)))

def vprint(v, *args, level=1):
    if v >= level: print(*args)
```

### Derivative sanity checks
- Use `python/pyMolecular/plotUtils.py` to compare analytical vs numerical derivatives.
- Note: `numDeriv(x,y)` returns `-dy/dx` (force from energy); keep sign conventions explicit.

```python
import numpy as np
from pyMolecular.plotUtils import plot1d

x = np.linspace(-2, 2, 2001)
E = 0.5*x**2
F = -x
fig, _ = plot1d(x, [E], derivs=[F], labels=['test'], bNumDeriv=True)
# inspect mismatch between analytical and numerical overlays
```

### GPU/CPU dtype mismatch quick check
- Common bug when debugging PyOpenCL vs NumPy baselines:
  - Assert host arrays: `assert a.dtype == np.float32` before upload (if kernels use float)
  - Cast explicitly: `a32 = np.asarray(a, dtype=np.float32, order='C')`
  - Compare with tolerance: `np.allclose(gpu, cpu.astype(np.float32), rtol=1e-5, atol=1e-6)`

## Cross-links
- See `doc/CodingRules/Python_Numerics_Rules.md` for broader numerics/plotting rules
- See `doc/CodingRules/workflows/python/python_demo_workflow.md` for initial sketching patterns and compact examples

## Notes
- Avoid broad try/except; fail loudly with stack traces
- Keep core function signatures stable while instrumenting
