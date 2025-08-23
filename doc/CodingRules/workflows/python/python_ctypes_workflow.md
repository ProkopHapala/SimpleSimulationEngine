# Workflow: Python + C/C++ (ctypes)

## Purpose
Maintain and synchronize a C/C++ shared library and its Python `ctypes` wrapper with tests and plots on the Python side.

## Inheritance
- Inherits `doc/CodingRules/global_rules.md`

## Workflow Steps
1. Setup
   - Prefer auto-build + load via `python/pyMeta/cpp_utils.py`: set `cpp_utils.BUILD_PATH` to lib build dir; `lib = cpp_utils.loadLib('<name>')` (auto-recompile if needed); capture logs
   - Alternatively use project bash scripts (`tests_bash/.../*.sh`, `C/build.sh`) when appropriate
   - Load functions and set `argtypes`/`restype` explicitly
2. Coding strategy
   - Mirror C/C++ headers in Python types (`ctypes.c_double`, pointers); ensure NumPy dtypes/contiguity match
   - Keep wrappers thin; convert inputs once at boundary
   - Data modeling (Python side): prefer tuples/dataclasses over dictionaries; never use a Python dict as a replacement for a C struct. For interop, use `ctypes.Structure` or NumPy structured dtypes (AoS) with explicit field layout
   - Vectorization: avoid Python loops and if/then/else branches in hot paths; use NumPy array ops and boolean masks. If container growth is needed, use list/dict/set comprehensions outside hot loops
3. Debugging strategy
   - On crash, verify pointer alignment, shapes, strides, dtypes
   - Print addresses and small array slices; add C-side `printf` with guards
4. Testing & validation
   - Maintain a small Python harness comparing CPU NumPy vs C/C++ outputs
   - Update tests on signature change; keep wrappers and headers in sync
5. Performance considerations
   - Zero‑copy where possible; pass `.ctypes.data_as(...)` for buffers
   - Batch calls to reduce Python↔C overhead
6. Visualization/reporting
   - Plot results in Python; parse C/C++ logs if needed

## Repository pattern (generic) + examples
- Python wrapper module: `python/<package>/<Name>.py`
  - Prefer `python/pyMeta/cpp_utils.py` to auto-rebuild and load shared libs: set `cpp_utils.BUILD_PATH = <.../cpp/Build/libs/<Domain>>`; `lib = cpp_utils.loadLib('<libname>')`
  - Define `np.ctypeslib.ndpointer` and set `argtypes`/`restype`; use a helper like `_np_as` for pointer casts when useful
  - Variant: direct `ctypes.CDLL` with a local `make` step is acceptable (see examples), but using `cpp_utils` is recommended for consistency
- C ABI bridge source: `cpp/libs/<Domain>/<Name>_lib.cpp` or `cpp/libs/<Domain>/<Name>.cpp`
  - Export `extern "C"` functions (no C++ name mangling) wrapping C++ engine classes
  - Expose plain C signatures; return raw pointers for zero‑copy buffers when needed
- C++ implementation: `cpp/common/<domain>/**/*.h/.cpp` (engine classes, heavy logic)

Examples in this repo (non‑exhaustive):
- eFF: `python/pyMolecular/eFF.py` ↔ `cpp/libs/Molecular/eFF_lib.cpp` ↔ `cpp/common/molecular/eFF.h`
- RigidMol: `python/pyMolecular/RigidMol.py` ↔ `cpp/libs/Molecular/RigidMol.cpp` ↔ engine headers under `cpp/common/...`
- CLCFGO: `python/pyMolecular/CLCFGO.py` ↔ `cpp/libs/Molecular/CLCFGO_lib.cpp` ↔ engine headers under `cpp/common/...`

Example (Python setup):
```python
from pyMeta import cpp_utils
cpp_utils.BUILD_PATH = cpp_utils.PACKAGE_PATH + '../../../cpp/Build/libs/<Domain>'
lib = cpp_utils.loadLib('<libname>')  # e.g., 'eFF_lib', 'RigidMol', 'CLCFGO_lib'
```
Auto-build toggles (in `python/pyMeta/cpp_utils.py`):
- `clean_build` (default True): if True, runs a clean rebuild; set False to avoid full rebuilds during iteration/tests
- `recompile_glob` (default True): if True, `loadLib(..., recompile=True)` triggers `make(<name>)`; set False to skip auto-make

## Mapped buffers and switches (general pattern)
- __Buffer registries (C/C++)__: store named pointers and expose via simple getters; populate in `makeDefaultBuffers()`.
```cpp
std::unordered_map<std::string,double*> buffers;
std::unordered_map<std::string,int*>    ibuffers;
// filled in makeDefaultBuffers(); expose extern "C" getBuff(name)/getIBuff(name)
```

- __Python views (ctypes + NumPy)__: create zero‑copy arrays from names and shapes.
```python
apos = getBuff("apos", (n,3)); atype = getIBuff("atype", (n,))
epos = getBuff("epos", (norb, perOrb, 3))
```

- __Tri‑state switches (-1/0/+1)__: -1 disable, +1 enable, 0 keep.
```cpp
#define SET_BOOL(b,i) do{ if(i>0) b=true; else if(i<0) b=false; }while(0)
// extern "C" void setSwitches_(int normalize, int kinetic, int coulomb, ...){
//     SET_BOOL(solver.bNormalize, normalize); /* ... */
// }
```
```python
lib.setSwitches_.argtypes = [c_int]*N
lib.setSwitches_(+1, 0, -1, ...)  # enable, keep, disable
```
Notes: arrays must be contiguous with matching dtype/shape; buffer ownership remains on C/C++ side.

## ABI / interop rules
- Use `extern "C"` in `*_lib.cpp` for all exported functions; keep signatures C‑friendly
- Pass contiguous NumPy arrays with matching `dtype` and shape; avoid implicit copies
- Expose buffers as raw pointers when zero‑copy is desired; document ownership semantics
- Mind 64‑bit pointers and platform ABI; keep struct layout stable if exposed

## Notes
- Sync triad: C/C++ header ↔ ctypes wrapper ↔ test script
- Document ABI expectations (calling convention, struct layout) if non‑trivial

## Sync / change workflow
- When changing C++ engine (`cpp/common/.../*.h`/`.cpp`):
  - Update the `extern "C"` bridge (`cpp/libs/.../*_lib.cpp`) accordingly
  - Rebuild via `cpp_utils.loadLib(...)` (auto) or project bash script
  - Update Python `argtypes`/`restype` in the wrapper; run the small harness
- Keep diffs minimal; verify with a tiny problem size; plot/print summaries on Python side if helpful
