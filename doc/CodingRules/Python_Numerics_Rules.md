# Python numerical coding rules (NumPy/Matplotlib)

Concise rules tailored for developing numerical algorithms in python.

## Core principles
- Minimal dependencies: standard library + NumPy/Matplotlib, unless user ask otherwise
- Keep code concise and minimal; prefer small, testable changes, step-by-step, divide-and-conquer
- Reuse existing functions; avoid duplication; 
   - check if usefull functions exist and report if they need modification
- Prefer pure, data-oriented functions with explicit inputs/outputs
   - use default named arguments to avoit exesively long call-strings
- Fail loudly: avoid broad try/except; use assertions for invariants
- Comment out deprecated code when experimenting (do not delete it); mark with TODO/DEBUG

## Performance (NumPy-first)
- Avoid Python loops and branches in hot paths, they are slow
   - prefer pre-allocated numpy arrays over lists
   - prefer using numpy boolen array masks and advanced indexing over if branches 
   - if you need grow datastructute use list/dict/set comprehensions
- prefer tuples and data-classes over dictionaries
- never use Python dict as a replacement for a C struct; for fixed-layout records use `dataclasses` or NumPy structured dtypes
- be mindful of memory layout and dtype (default float64, for pyOpenCL use float32)
- prefer importain library at the start of script/module rather then in functions

## Debugging/logging, test-driven development, (debuggability > UX)
- When you desing algorithm thing about how to debug it and test it from the start
  - for LLM we need information in terminal output in text-form
    - make sure debug prints are easy to read not exesively long, like a table
  - for user the plots and visializations are preferable
  - Add quick assertions and end-of-run sanity checks (area/mass/energy)
- at initial stages of development add lot of debug-prints which help follow the algorithm and spot errors
   - be carefull adding too many debug prints inside hot loops  
- Gate all debug prints by a verbosity level; keep prints as one-liners
- If multi-line diagnostics are needed, wrap them in a small helper (to avoid bloating the code with debug-prints)
- Do not store debug-only data in core objects by default
- Prefer simple prints to heavy logging frameworks
- for key numpy arrays print overall debug information (e.g. min,max,shape) which help use spot strange values (nan, inf, etc.) 

## Plotting (Matplotlib)
- Separate compute and visualization; no plotting in core algorithms
- Plotting should be optional: add flags like `--noPlot`, `--saveFig`
- Call `plt.show()` only in CLI/main, never in library code
- Prefer to reuse helpers from `python/pyMolecular/plotUtils.py` (avoid re‑implementing): `plot1d`, `plot_func`, `plot_funcs`, molecule plotting helpers
- make plotting code concise
   - unless asked otrwise, fancy titples, labels and formating is not necessary, do not spend many lines of code on that 

### Derivative checks (analytical vs numerical)
- When a function has an analytical derivative, always cross‑check with numerical.
  - Use `plotUtils.numDeriv(xs, y)` to get a numerical derivative (note: returns `-dy/dx`, i.e. force from energy).
  - Overlay in plots via `plot1d(xs, [y], derivs=[dy], labels=[...], bNumDeriv=True)` or use `plot_funcs` when the callable returns `(E, F, [S])`.
  - Keep sign conventions explicit in labels (e.g., `F = -dE/dx`).

## API/CLI
- Add `if __name__ == "__main__":` guard to CLI code instead of `main()` function
- Use argparse for CLI arguments with most global parameters unless too long 
- Provide `--verbosity` to control diagnostics; default to quiet

## Style
- Short, clear variable names, idealy representing using physical or mathematicall symbols (`E` for Energy, `T` for temparature etc.) 
- Compact code, avoiding empty lines, prefering one-liners
- Avoid line wrapping when it harms readability of expressions (assume wide editor)
   - at the start of all files set numpy print options allow infinite line-length `np.set_printoptions(linewidth=np.inf)` 
- Inline comments behind code line 
   - reserve header comments only for high-level structures (functions, classes, modules)
