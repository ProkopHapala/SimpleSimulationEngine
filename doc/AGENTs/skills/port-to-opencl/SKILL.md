---
name: port-to-opencl
description: Porting Fortran/C++→PyOpenCL — OpenCLBase, kernel caching, persistent buffers, GPU performance
trigger:
  glob:
    - "**/*.cl"
    - "**/fortran/**/*"
    - "**/cpp/**/*.cpp"
    - "**/cpp/**/*.h"
    - "**/pyBall/OCL/**/*"
---

## Why PyOpenCL First

Always start with PyOpenCL before CUDA. Benefits:
- On-the-fly compilation: edit `.cl` files, reload instantly
- NumPy integration: zero-copy buffers, seamless data transfer
- Rich diagnostics: build logs, kernel introspection, printf debugging
- Visualization: matplotlib plots directly from GPU results
- Fast iteration: no rebuild cycles, test changes immediately

**GPU is always single-precision (float32).** Use `%f` instead of `%g`, avoid double for numerical speed.

## Base Class: OpenCLBase

Inherit from `pyBall/OCL/OpenCLBase.py` for efficient GPU resource management:

**Kernel caching:** Compile once during `__init__`, cache in `self.prg`. Skip recompilation on subsequent calls.
```python
if not self.load_program(rel_path="../../cpp/common_resources/cl/FitREQ.cl", base_path=base_path):
    exit(1)
```

**Build option caching:** If your GUI toggles compile-time flags (collision on/off, debug prints on/off, dynamics vs relaxation), cache the last-used build flags tuple and only recompile when flags actually change. Do not recreate the program object on every checkbox click — multi-second freezes mask real bugs and make interactive debugging impossible.

**Persistent buffer management:** Allocate once, reuse across calls. Use `try_make_buffers()` which checks size and only reallocates if needed.
```python
buffs = {"input": sz, "output": sz}
self.try_make_buffers(buffs, suffix="_buff")  # stored in self.buffer_dict
```

**bTryAllocate guards:** Guard buffer allocation with `if bTryAllocate:` to skip dict creation and allocation in hot paths.
```python
def my_kernel(self, data, bTryAllocate=True):
    if bTryAllocate:
        buffs = {"input": sz, "output": sz}
        self.try_make_buffers(buffs, suffix="_buff")
    self.toGPU_(self.input_buff, data)
    self.prg.my_kernel(self.queue, gs, ls, self.input_buff, self.output_buff)
    return self.fromGPU_(self.output_buff)
```

First call: allocates buffers. Subsequent calls with same sizes: skips allocation.

## Porting Principles (Fortran→OpenCL)

- Reference is truth: If outputs differ, OpenCL is wrong.
- Preserve identity: Use unique keys (e.g., (iatom, ineigh, mbeta)). Don't collapse by (iatom, jatom). Respect periodic shifts.
- Match inputs exactly: Copy rotation frames, epsilon, twister inputs per term. Don't reuse generic frames.
- Preserve structure: Replicate overwrite patterns (e.g., transpose/conjugate for reverse slots), not just accumulation.
- Separate checks: Keep structural (integer) checks separate from float parity. Expect exact zeros only in structural checks.

**Workflow:** Setup references → Scan microtests → Structural integrity → Block-level parity → Dense reassembly → End-to-end term → Debug drill (bisect).

## GPU Kernel Perforamnce Guidelines

- Design around memory latency and cache efficiency.
- Prefer "gather" operations over "scatter" designs where possible.
- Avoid branching, atomics, and unnecessary synchronization in GPU kernels.
- Maximize the efficiency of shared/local memory and workgroups.
- Avoid unnecessary host-device data transfers.

