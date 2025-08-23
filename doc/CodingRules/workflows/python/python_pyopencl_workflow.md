# Workflow: Python + PyOpenCL

## Purpose
Develop GPU kernels in `.cl` with a Python host. Manage buffers, macro substitutions, and validation against CPU NumPy.

## Inheritance
- Inherits `doc/CodingRules/global_rules.md`

## Workflow Steps
1. Setup
   - Use project `OpenCLBase` utilities (context, queue, buffers)
   - Keep kernels under `kernels/`; prepare macro markers for metaprogramming
2. Coding strategy
   - Define kernel args clearly; prefer `float32` on GPU
   - Auto‑parse kernel headers where possible; generate buffer maps from dicts
   - Host data modeling: prefer tuples/dataclasses over dictionaries; never use a Python dict as a replacement for a C struct. For fixed-layout data, use NumPy structured dtypes or small dataclasses
   - Host vectorization: avoid Python loops and if/then/else branches in hot paths; use NumPy array ops and boolean masks. If container growth is needed, use list/dict/set comprehensions outside hot loops
3. Debugging strategy
   - Start with tiny sizes; compare GPU vs CPU reference arrays
   - Always print kernel build logs on compile errors; dump expanded macros
   - Use `cl.enqueue_copy` to inspect device buffers
4. Testing & validation
   - Add deterministic tests with fixed seeds; tolerance‑based comparisons
   - Record timing and basic throughput for sanity
5. Performance considerations
   - Align work‑group sizes; minimize host↔device transfers
   - Reuse device buffers; avoid reallocation in loops
6. Visualization/reporting
   - Plot correctness deltas and timing curves in Python

## OpenCLBase helpers (project patterns)
- __Buffer management__: allocate from a dict of byte sizes via `OpenCLBase.try_make_buffers({'name': nbytes, ...})`; creates `<name>_buff` and mirrors into `buffer_dict`. Use `create_buffer()` for explicit names and `check_buf()` to grow/shrink (0 releases).
- __Transfers__: name-based `toGPU('name_buff', host)` / `fromGPU('name_buff', host)` for convenience; low-level `toGPU_(buf, host, offset)` / `fromGPU_(buf, host, offset)` when holding `cl.Buffer` and using byte offsets.
- __Headers → args__: `load_program()` auto-fills `kernelheaders`; put scalars in `self.kernel_params` (match header names), keep device-arg names in `buffer_dict`, then call `generate_kernel_args('kernel')` to get the argument tuple in correct order.
- __Source templating__: specialize kernels by injecting macros/functions at exact-line sentinels using `preprocess_opencl_source()` (e.g., `//<<<MARKER`). Extract reusable blocks from `.cl` libraries marked by `//>>>function NAME` / `//>>>macro NAME` via `parse_cl_lib()`.
- __Work sizes__: choose `nloc` (preferred local size) and pad with `roundUpGlobalSize(N)`; minimize transfers inside loops.

Examples in this repo follow these patterns (see `python/pyMolecular/OCL/MolecularDynamics.py` for persistent buffers + header‑driven args; `python/pyMolecular/OCL/run_scanNonBond.py` for macro injection and quick compile‑run scans). Keep usage generic; avoid hard‑coding shapes/types where not required.

## Cross-links
- Align dtype guidance with the demo workflow: default `float64` on CPU (NumPy), prefer `float32` on GPU buffers (PyOpenCL).
- For plotting and derivative sanity checks, reuse `python/pyMolecular/plotUtils.py` helpers (`plot1d`, `plot_func`, `plot_funcs`, `numDeriv`) per `doc/CodingRules/Python_Numerics_Rules.md` and `doc/CodingRules/workflows/python/python_demo_workflow.md`.

## Notes
- Keep kernel code simple first; optimize after correctness
- For full API and concise examples see `docs/OpenCLBase.md`
