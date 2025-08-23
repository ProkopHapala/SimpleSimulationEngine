# OpenCLBase: Minimal OpenCL helper for FireCore

File: `pyBall/OCL/OpenCLBase.py`

OpenCLBase is a small, pragmatic layer over PyOpenCL tailored for physics simulations. It removes boilerplate around device/context setup, buffer management, kernel header parsing, kernel-argument assembly, and simple source templating/macros.

- __Focus__: correctness, simplicity, small API; favors lists/tuples and straight-line code.
- __Used by__: `HubbardSolver.py`, `NonBondFitting.py`, `MolecularDynamics.py`, `run_scanNonBond.py`.

---

## Quick start (minimal)

```python
from pyBall.OCL.OpenCLBase import OpenCLBase
import numpy as np

ocl = OpenCLBase(nloc=64)
ocl.load_program(kernel_path="path/to/kernels.cl")

N = 100000
sizes = { 'pos': N*16, 'vel': N*16, 'force': N*16 }
ocl.try_make_buffers(sizes)                   # creates pos_buff, vel_buff, force_buff (+ buffer_dict entries)

pos = np.zeros((N,4), np.float32)
vel = np.zeros((N,4), np.float32)
force = np.zeros((N,4), np.float32)
ocl.toGPU('pos_buff', pos); ocl.toGPU('vel_buff', vel)

ocl.kernel_params = {'N': N, 'dt': 1e-3}
args = ocl.generate_kernel_args('integrate')  # [pos_buff, vel_buff, force_buff, N, dt]

g = (ocl.roundUpGlobalSize(N),)
ocl.prg.integrate(ocl.queue, g, (ocl.nloc,), *args)

ocl.fromGPU('force_buff', force)
```

---

## Architecture and conventions

- __Context/queue__: created in `__init__()` via `select_device()`; prints device info (`clu.get_cl_info`).
- __Attributes__:
  - `self.nloc`: preferred local size (work-group size).
  - `self.ctx`, `self.queue`: PyOpenCL context/queue.
  - `self.prg`: built `cl.Program`.
  - `self.buffer_dict`: name→`cl.Buffer` map (also mirrors buffers created by `try_make_buffers`).
  - `self.kernelheaders`: kernel name→header string, auto-filled by `load_program()`.
- __Naming__: make buffer names match kernel parameter names for device arguments; scalars are provided in `self.kernel_params`.

---

## Function reference

Top-level utilities:
- `print_devices(platforms=None)`: Quickly list available platforms/devices. Use it to confirm device names and indices before selecting one.
- `select_device(platforms=None, preferred_vendor='nvidia', bPrint=False, device_index=0)`: Create a context on a device whose name contains `preferred_vendor` (case-insensitive); falls back to default selection. Optionally prints the chosen device for transparency.

Class `OpenCLBase`:
- `__init__(self, nloc=32, device_index=0)`: Initialize context and command queue, print basic device info, and set the preferred local size `nloc`. Also prepares registries (`buffer_dict`, `kernelheaders`, `prg`) used by other helpers.
- `load_program`: Read a `.cl` file and build a `cl.Program`. Optionally extracts kernel headers immediately so later you can assemble arguments automatically, and can print them for debugging.
- `extract_kernel_headers`: Parse the source text and capture each `__kernel` signature (robust to comments and multi-line formatting). Produces `{kernel_name: header_string}` used by the header→args pipeline.
- `create_buffer`: Allocate a device buffer of a given byte size and store it under `buffer_dict[name]`. Prefer this when you want explicit dictionary-managed buffers without attribute helpers.
- `check_buf`: Ensure `buffer_dict[name]` exists and has at least `required_size` bytes; releases the buffer if `required_size==0`. Ideal for dynamic data that grows/shrinks between runs without leaking memory.
- `try_make_buff`: Ensure `self.<buff_name>` exists and has exactly `sz` bytes; re-allocates only on size change. Returns `(buf, created)` so you can perform one-time initialization when a buffer was newly created.
- `try_buff`: Convenience wrapper that only creates `name+suffix` if `name` is present in a provided list. Useful when enabling optional features without branching elsewhere.
- `try_make_buffers`: Batch-ensure multiple buffers given `{name: size_in_bytes}`. Creates attributes `<name>_buff` and mirrors new allocations into `buffer_dict` for easy I/O.
- `toGPU_`: Low-level host→device copy into a given `cl.Buffer` (no name lookup), with optional `byte_offset`. Use inside tight loops when you already hold the buffer object.
- `fromGPU_`: Low-level device→host copy; if `host_data` is `None` it allocates one with given `shape`/`dtype` and returns it. Handy for quick reads without pre-allocating.
- `toGPU`: Name-based host→device copy via `buffer_dict[buf_name]`. The simplest path when using the naming conventions established by `try_make_buffers`.
- `fromGPU`: Name-based device→host copy via `buffer_dict[buf_name]`; mirror of `toGPU` for downloads.
- `bufflist`: Return a list of buffers by name; useful for compactly passing many buffers to kernels or validation utilities.
- `roundUpGlobalSize`: Pad a global work-size to a multiple of `nloc` so kernels can assume full groups and avoid guard branches.
- `parse_kernel_header`: Extract an ordered list of `(arg_name, kind)` from a kernel header, robust to comments and line breaks; `kind=0` for device objects (`__global` buffers/images), `kind=1` for scalars. This is the bridge from source code to runtime argument wiring.
- `generate_kernel_args`: Build the final argument tuple in header order by resolving device objects from `buffer_dict` and scalars from `self.kernel_params`. Reduces boilerplate and catches missing names early with optional verbose printing.
- `preprocess_opencl_source`: Lightweight metaprogramming: replace exact-line sentinels with file/function/macro content to specialize kernels (e.g., different pair potentials) without duplicating code. Exact matching avoids accidental replacements inside comments/docs.
- `parse_cl_lib`: Read a `.cl` "library" and extract blocks marked by `//>>>function NAME` and `//>>>macro NAME` into dictionaries. Pair with `preprocess_opencl_source` to inject model-specific code where needed.
- `parse_forces_cl`: Backward-compatible alias of `parse_cl_lib` kept for existing callers.

---

## Buffer management patterns

- __Minimal dynamic sizing__ (Hubbard/MolecularDynamics style):
```python
sz = N*16
ocl.try_make_buffers({'pos': sz, 'vel': sz, 'force': sz})  # ensures pos_buff, vel_buff, force_buff
```

- __Explicit dictionary-managed buffer__:
```python
o.cl.create_buffer('grid', nx*ny*4)
o.cl.toGPU('grid', grid_np)
```

- __Resize-or-release__:
```python
o.cl.check_buf('tmp', required_size=M*4)   # grows if needed; size 0 releases
```

- __Batch set kernel buffers__:
```python
bufs = ocl.bufflist(['pos_buff','vel_buff','force_buff'])
```

Notes:
- Sizes are bytes; ensure dtype-stride matches your kernel.
- Keep names aligned with kernel parameter names for device arguments.

---

## Data transfer

```python
ocl.toGPU('pos_buff', pos_np)
ocl.fromGPU('force_buff', force_np)
# low-level
ocl.toGPU_(ocl.buffer_dict['pos_buff'], pos_np)
```

- Use `_` variants when you already have the `cl.Buffer` object.
- Avoid transfers inside hot loops; upload once and keep on device.

---

## Kernel headers → arguments without boilerplate

1) Build program and auto-extract headers:
```python
o.cl.load_program(kernel_path='kernels.cl')
```

2) Provide scalars in `self.kernel_params` using exact header names:
```python
o.cl.kernel_params = {'N': N, 'dt': dt}
```

3) Ensure device-arg names exist in `buffer_dict` (or via `try_make_buffers`).

4) Generate args in the correct order:
```python
args = ocl.generate_kernel_args('integrate', bPrint=True)
```

5) Launch:
```python
g = (ocl.roundUpGlobalSize(N),)
ocl.prg.integrate(ocl.queue, g, (ocl.nloc,), *args)
```

Header parsing handles multi-line and comments; device-args include `__global` buffers and images; scalars are everything else.

---

## OpenCL code-generation (Metaprogramming)

Often it happens that we want to specialized general OpenCL kernel into many variants (e.g. evaluation of non-bodning interactions can use different potential-energy function like Morse, Coulomb, Lenard-Jones etc.). OpenCLBase class provide convenient tools to do that using string-substitution.

Use exact-line sentinels for macros to avoid accidental replacements (safer and verified in practice).

- __Parse a CL library for blocks__:
```python
lib = ocl.parse_cl_lib('cpp/common_resources/cl/Forces.cl')
macro_body = lib['macros']['MODEL_MorseQ_PAIR']
```

- __Inject blocks into a template__:
```python
subs = {
  'files':     {},
  'functions': {},
  'macros':    { 'MODEL_PAIR_ACCUMULATION': macro_body }  # matches line: //<<<MODEL_PAIR_ACCUMULATION
}
code = ocl.preprocess_opencl_source('cpp/common_resources/cl/FitREQ.cl', subs, output_path='build/FitREQ_templated.cl', bPrint=True)
```

- __Build the preprocessed source__:
```python
o.cl.prg = cl.Program(o.cl.ctx, code).build()
o.cl.kernelheaders = ocl.extract_kernel_headers(code)
```

Guidelines:
- Macro sentinels must be the entire line: `//<<<MARKER` (no trailing text).
- Prefer small, composable macros/functions. Keep injection points stable.

---

## Work sizes

```python
global_size = (ocl.roundUpGlobalSize(N),)
local_size  = (ocl.nloc,)
```

Pick `nloc` that matches device’s preferred multiple; `roundUpGlobalSize` pads to the next multiple.

---

## Usage in FireCore modules

- __HubbardSolver (`pyBall/OCL/HubbardSolver.py`)__
  - Reallocation via `try_make_buffers` in `realloc_*` flows.
  - Kernels: `solve_MC_2phase`, `solve_MC_neigh` with setup helpers mapping buffers/scalars.
  - Pattern: upload spin/state arrays once; regenerate arg tuples via headers.

- __NonBondFitting (`pyBall/OCL/NonBondFitting.py`)__
  - Model templating: `preprocess_opencl_source` injects `MODEL_*` macro blocks into `evalSampleDerivatives_template`.
  - Switch between template and baseline kernel; manages parameter/sample buffers.

- __MolecularDynamics (`pyBall/OCL/MolecularDynamics.py`)__
  - Uses `generate_kernel_args` to assemble calls; keeps positions/velocities/forces as persistent buffers.

- __run_scanNonBond (`pyBall/OCL/run_scanNonBond.py`)__
  - Macro substitution + compile + run scans; quick experiments with nonbond models.

---

## Best practices and pitfalls

- __Name alignment__: device-arg names in headers must match keys in `buffer_dict` (e.g., `pos_buff` vs `pos` → keep consistent).
- __Scalars in `kernel_params`__: ensure present before `generate_kernel_args`.
- __Byte sizes__: always compute sizes in bytes (`N * dtype.itemsize * vec_len`).
- __Debug mapping__: use `generate_kernel_args(..., bPrint=True)` to verify order and values.
- __Macro safety__: only replace exact-line sentinels; do not use raw string replace inside comments/file scope.
- __Resize discipline__: `check_buf(name, required_size)` to grow/shrink; 0 releases.

---

## Minimal end-to-end example (concise)

OpenCL:
```c
__kernel void saxpy(__global float* y, __global const float* x, float a, int n){
  int i = get_global_id(0);
  if(i<n) y[i] = a*x[i] + y[i];
}
```

Python:
```python
ocl = OpenCLBase(nloc=128)
ocl.load_program(kernel_path='kernels.cl')

n = 1_000_000
x = np.ones(n, np.float32)
y = np.zeros(n, np.float32)
sz = n*4
ocl.try_make_buffers({'x': sz, 'y': sz})
ocl.toGPU('x_buff', x); ocl.toGPU('y_buff', y)

ocl.kernel_params = {'a': 2.0, 'n': n}
args = ocl.generate_kernel_args('saxpy')

g = (ocl.roundUpGlobalSize(n),)
ocl.prg.saxpy(ocl.queue, g, (ocl.nloc,), *args)
ocl.fromGPU('y_buff', y)
```

---

## Appendix: full API

- Top-level: `print_devices()`, `select_device()`.
- Class:
  - `__init__`, `load_program`, `extract_kernel_headers`.
  - Buffers: `create_buffer`, `check_buf`, `try_make_buff`, `try_buff`, `try_make_buffers`, `bufflist`.
  - Transfer: `toGPU_`, `fromGPU_`, `toGPU`, `fromGPU`.
  - Kernel args: `parse_kernel_header`, `generate_kernel_args`, `roundUpGlobalSize`.
  - Source: `preprocess_opencl_source`, `parse_cl_lib`, `parse_forces_cl`.

If you want cross-links or diagrams, we can add them later. The goal here is to give you fast, correct patterns with minimal code.
