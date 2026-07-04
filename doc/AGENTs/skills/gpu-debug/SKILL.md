---
name: gpu-debug
description: Debugging GPU/OpenCL/CUDA kernels — gated prints, CPU↔GPU tracing, barrier pitfalls
trigger:
  glob:
    - "**/*.cl"
    - "**/*.cu"
    - "**/apps_OCL/**/*"
    - "**/apps_CUDA/**/*"
    - "**/pyBall/OCL/**/*"
---

## Gated Debug Macros

Define at kernel top:
```c
#define DBG_UFF       1    // Enable/disable all prints
#define IDBG_SYS      0    // Which system to print
#define IDBG_ATOM     0    // Which atom to print
```

Gate prints in kernel:
```c
#if DBG_UFF
if(get_global_id(0)==IDBG_ATOM && get_global_id(1)==IDBG_SYS){
    printf("GPU_BOND %d: L %f K %f F %f\n", i_bond, L, K, F);
}
#endif
```

## Verbosity & Targeting

Gate OpenCL `printf` with compile-time flags and runtime verbosity levels:

| Level | Purpose | Cost |
|-------|---------|------|
| 0 | Silent (CI) | None |
| 1 | Events (start/stop/toggle) | Minimal |
| 2 | Per-workgroup summary (topology, bbox, COG) | Low |
| 3 | Per-atom dumps | High — always target |

**Component bitmasks** to enable only specific logs:
```c
#define DEBUG_COMPONENTS 0x01  // collision
#define DEBUG_COMPONENTS 0x02  // port constraints
#define DEBUG_COMPONENTS 0x04  // rotation
#define DEBUG_COMPONENTS 0x08  // momentum
```

**Workgroup targeting** — only print atoms you care about:
```c
-DENABLE_DEBUG_PRINTS
-DDEBUG_VERBOSITY=3
-DDEBUG_TARGET_WG=5
-DDEBUG_GID_START=320 -DDEBUG_GID_END=384
```

**Build option caching** — if the GUI toggles debug prints on/off, cache the last build-flags tuple and only recompile when flags actually change. Multi-second recompilation freezes mask real bugs.

## Synchronized Tracing

When comparing CPU vs GPU, inject **identical** printf format in both implementations for automated diff checks.

**GPU is always single-precision (float32).** Use `%f` instead of `%g`, avoid double for numerical speed.

**CPU side:**
```cpp
if(iatom==IDBG_ATOM && isys==IDBG_SYS){
    printf("CPU_BOND %d: L %g K %g F %g\n", i_bond, L, K, F);
}
```

**GPU side:** Match format exactly (same variable order, same precision specifiers).

## Barrier Deadlock

**Critical:** Never return early before a barrier. Work-group size is rounded up to local size multiples, creating "dummy" work-items that must reach all barriers.

**Wrong:**
```c
if(iG >= natoms) return;  // DEADLOCK
barrier(CLK_LOCAL_MEM_FENCE);
```

**Right:**
```c
const bool bValid = (iG < natoms);
if(bValid){ /* compute */ }
barrier(CLK_LOCAL_MEM_FENCE);  // All work-items reach here
```

## GPU-Specific Pitfalls

- **Struct alignment:** OpenCL `float4` (16 bytes) vs C++ `Vec3d` (24 bytes) or `Vec3f`. Verify strides.
- **System replicas:** OpenCL is Multi-System (`nSystems`). C++ Reference often Single-System. Compare System 0.
