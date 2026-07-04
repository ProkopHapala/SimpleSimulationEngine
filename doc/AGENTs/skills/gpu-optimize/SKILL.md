---
name: gpu-optimize
description: Optimizing GPU/OpenCL kernel performance — memory hierarchy, gather vs scatter, tiling, auxiliary arrays, local memory reductions, kernel launch overhead
trigger:
  glob:
    - "**/*.cl"
    - "**/*.cu"
    - "**/apps_OCL/**/*"
    - "**/apps_CUDA/**/*"
    - "**/pyBall/OCL/**/*"
---

## Memory Latency & Cache

- Design around memory latency — coalesced global reads, cache-friendly access patterns
- Prefer contiguous access patterns; avoid random-access into global memory
- Use `__attribute__((aligned(16)))` for structs to match cache line boundaries

## Data Layout & Aligned Types

- Prefer `float4` / `float2` / `int4` — hardware loads 16-byte aligned vectors in single transaction
- **Avoid structs in kernels** — OpenCL allows them but they cause poor coalescing and padding waste
- **Structure-of-Arrays (SoA) over Array-of-Structs (AoS)**: separate `float4* pos`, `float4* vel`, `float* charge` — not `struct Atom{ float4 pos; float4 vel; float q; }* atoms`
- SoA enables coalesced reads: all threads read `pos[i]` from contiguous memory. AoS scatters reads across struct stride.
- **Pack heterogeneous data into single `float4*` or `int4*`**: e.g. `{x, y, z, type}` in one float4, `{charge, mass, flags, _pad}` in another — one load fetches multiple fields
- We often use flat `float4*` arrays for positions, velocities, forces — packed as `{x, y, z, w}` where w carries extra scalar data

```
// Bad: AoS — strided access, poor coalescing
struct Atom { float4 pos; float4 vel; float q; };
__global struct Atom* atoms;  // each thread reads atoms[i].pos = 48-byte stride

// Good: SoA — coalesced, each field contiguous
__global float4* pos;   // atoms[i].pos = pos[i], 16-byte stride
__global float4* vel;   // atoms[i].vel = vel[i]
__global float*  q;     // atoms[i].q   = q[i]

// Also good: pack into float4 arrays
__global float4* pos;   // {x, y, z, type}
__global float4* props; // {charge, mass, flags, _pad}
```

## Gather Over Scatter

- Prefer "gather" (each thread reads from fixed indices) over "scatter" (each thread writes to arbitrary indices)
- Scatter requires atomics or serialization; gather is naturally parallel
- If scatter is unavoidable, use atomic operations sparingly or restructure as gather + reduction

## Local Memory & Workgroups

- NVIDIA: 48KB local memory per compute unit (shared across workgroups on same CU)
- Workgroup size: 32 (warp) minimum, 64/128/256 typical, 1024 max — must be multiple of 32
- Tune workgroup size to occupancy — too small wastes SMs, too large causes register spilling
- Maximize shared/local memory usage for data reused across work-items
- Use `barrier(CLK_LOCAL_MEM_FENCE)` to synchronize local memory access within workgroup
- Never return early before a barrier (see skill:`gpu-debug` for barrier deadlock patterns)
- Consider if algorithm needs global reductions — e.g. Jacobi avoids local memory reductions that CG requires

## Minimize Host-Device Transfers & Kernel Launches

- Avoid unnecessary `clEnqueueReadBuffer`/`clEnqueueWriteBuffer` in hot paths
- Batch multiple operations on GPU before reading results back
- Use persistent buffers (allocate once, reuse) — see `OpenCLBase.try_make_buffers()`
- Guard allocation with `bTryAllocate` to skip dict creation in hot paths
- **PyOpenCL kernel launch overhead is large** — minimize number of kernel executions in iterative methods (MD timesteps, Jacobi iterations)
- Parallelize as much work as possible in a single kernel (broad parallelism) rather than calling many small kernels in sequence
- Fuse multiple kernel operations into one where possible

## Branching & Atomics

- Avoid divergent branching (warp/workgroup divergence) — all threads should take same path
- **Avoid atomics especially in OpenCL** — they serialize execution and performance is unpredictable across vendors
- Prefer reduction trees (see Local Memory Reduction below) or separate kernels over atomic counters

## Precision

- GPU is always single-precision (float32). Use `%f` not `%g` in printf
- Avoid double precision unless absolutely required — <10x slower on my GPU
- Be aware of float32 precision limits in accumulation — use Kahan summation or compensated reduction if needed

## Thread Saturation

- Typical GPU has ~10k threads — ensure problem size saturates this
- If problem is too small, GPU is underutilized — consider running multiple systems/replicas in parallel
- If workgroups exceed thread count, they queue — balance workgroup count vs size

## Profiling

- Use `clGetEventProfilingInfo` to measure kernel execution time
- Identify bottlenecks: memory-bound vs compute-bound kernels
- Compare actual GFLOPs/s to theoretical peak to assess optimization headroom

## Related Skills

- skill:`port-to-opencl` — porting workflow, OpenCLBase class, kernel caching
- skill:`gpu-debug` — gated debug macros, CPU↔GPU tracing, barrier pitfalls

---

## Memory Hierarchy & Register Pressure

GPU memory hierarchy (fast → slow): **registers** (per-thread) → **local/shared memory** (per-workgroup) → **global memory** (device-wide). Each level is ~10-100x slower than the previous. Algorithm is usually limited by the smallest level.

| Level | Size | Latency | Scope |
|-------|------|---------|-------|
| Registers | ~255 per thread | 1 cycle | per-thread |
| Local memory | 48KB per CU | ~20 cycles | per-workgroup |
| Global memory | GBs | ~400 cycles | device-wide |

- **Register spilling**: too many live variables → spill to local memory (slow). Too much local memory → spill to global (very slow)
- **Scoping**: enclose independent computation blocks in `{ }` to release registers between blocks — variables from one block don't occupy registers in the next
- **Memory-conscious design**: prefer algorithms that use less memory at each level. Reading from small 1MB local memory is much faster than same number of random reads from 1GB global memory
- **Don't materialize large grids**: e.g. don't build full 3D basis function grid — materialize only 1D radial part per atom type, reconstruct 3D on the fly
- Workgroup occupancy is limited by register count per thread: more registers/thread → fewer concurrent workgroups → lower latency hiding

```
// Bad: registers i,j,k live simultaneously
float sum = 0;
for(int i=0; i<N; i++) sum += A[i];
float prod = 1;
for(int j=0; j<M; j++) prod *= B[j];
float norm = 0;
for(int k=0; k<K; k++) norm += C[k];

// Good: scopes release registers between blocks
float sum, prod, norm;
{
    float s = 0;
    for(int i=0; i<N; i++) s += A[i];
    sum = s;
}
{
    float p = 1;
    for(int j=0; j<M; j++) p *= B[j];
    prod = p;
}
{
    float n = 0;
    for(int k=0; k<K; k++) n += C[k];
    norm = n;
}
```

# GPU Design Patterns / Templates

Reusable kernel patterns for common GPU computation patterns.

## Tiled Design

Most N-body / matrix operations can use tiling instead of naive O(N²) global reads. Divide work into tiles that fit in local memory, fill collaboratively, then compute from local — amortizes global reads by factor of tile size.

```
// Tiled matrix multiply: C[i,j] = sum_k A[i,k]*B[k,j]
// Tile size TS (e.g. 32), workgroup = (TS, TS)
__local float A_tile[TS][TS], B_tile[TS][TS];

int tx = get_local_id(0), ty = get_local_id(1);
int i  = get_global_id(0), j  = get_global_id(1);
float acc = 0.0f;

for(int kt = 0; kt < N; kt += TS){
    // Collaborative load: each work-item loads one element of each tile
    A_tile[ty][tx] = A[i * N + kt + tx];
    B_tile[ty][tx] = B[(kt + ty) * N + j];
    barrier(CLK_LOCAL_MEM_FENCE);

    // Compute from local memory — TS reads from global amortized into 1
    for(int k = 0; k < TS; k++){
        acc += A_tile[ty][k] * B_tile[k][tx];
    }
    barrier(CLK_LOCAL_MEM_FENCE);  // before overwriting tile
}
C[i * N + j] = acc;
```

Same pattern applies to: N-body force calculation (tile of neighbor positions), pairwise distance matrices, convolution with shared tiles.

## Auxiliary Arrays & Ping-Pong (Avoid Scatter Without Atomics)

When each thread needs to write results to multiple neighbors (scatter), use auxiliary arrays written by one kernel and assembled by another — no atomics needed.

**Example: Newton's 3rd law in force calculation** — each bond/interaction produces forces on both atoms i and j. Instead of atomic-add to `f[atom]`, write recoils to auxiliary array, then sum in a second kernel.

```
// Kernel 1: compute forces, write to auxiliary arrays (no atomics)
// Each thread handles one bond (i,j), writes force on i and recoil on j
int ibond = get_global_id(0);
int i = bond_i[ibond], j = bond_j[ibond];
float3 f = computeForce(pos[i], pos[j], ...);

f_i[ibond] = f;        // force on atom i from this bond
f_j[ibond] = -f;       // recoil on atom j (Newton's 3rd law)

// Kernel 2: assemble — gather all contributions per atom
// Each thread handles one atom, sums all bond contributions
int iatom = get_global_id(0);
float3 total = (float3)(0,0,0);
for(int b : bonds_of[iatom]){   // gather, not scatter
    total += f_i[b] + f_j[b];   // or separate assembly kernels
}
force[iatom] = total;
```

**Ping-pong variant**: for iterative methods, alternate between two buffers (read A → write B, then read B → write A) to avoid read-write conflicts without synchronization.

```
// Iterative solver: alternate src/dst each iteration
// Even iterations: read bufA, write bufB
// Odd iterations:  read bufB, write bufA
__kernel void iterate(__global float* src, __global float* dst, ...){
    int i = get_global_id(0);
    dst[i] = update(src, i);  // no race: src and dst are different buffers
}
// Host: swap src/dst pointers between iterations
```

Use when: force assembly, neighbor contributions, iterative solvers (Jacobi, Gauss-Seidel), any scatter pattern that can be decomposed into write-then-gather.

## Local Memory Reduction

When you need a sum/max/min across work-items, use tree reduction in local memory — much faster than atomics.

```
// Reduce sum of N values within a workgroup (N = workgroup size, power of 2)
__local float sdata[WG_SIZE];
int tx = get_local_id(0);

sdata[tx] = input[get_global_id(0)];  // load
barrier(CLK_LOCAL_MEM_FENCE);

// Tree reduction: halve active threads each step
for(int stride = WG_SIZE / 2; stride > 0; stride >>= 1){
    if(tx < stride){
        sdata[tx] += sdata[tx + stride];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
}

// Thread 0 has the result
if(tx == 0) output[get_group_id(0)] = sdata[0];
```

For non-power-of-2: pad with zeros or handle remainder in first stride. For multi-group reductions: write per-group results to global, launch a second small kernel to reduce those.