---
name: cpu-perf
description: Optimizing CPU compiled code (C/C++/Fortran/Rust) — cache hierarchy, SIMD, data-oriented design, loop optimization, memory layout
trigger:
  glob:
    - "**/*.cpp"
    - "**/*.h"
    - "**/*.c"
    - "**/*.f"
    - "**/*.f90"
    - "**/*.rs"
---

## Memory Hierarchy & Cache

CPU memory hierarchy (fast → slow): **L1** (~32KB, ~1ns) → **L2** (~256KB, ~4ns) → **L3** (~8MB, ~12ns) → **DRAM** (GBs, ~100ns). Each level ~10x slower.

| Level | Size | Latency | Scope |
|-------|------|---------|-------|
| Registers | ~16 per thread | <1ns | per-thread |
| L1 cache | 32KB | ~1ns | per-core |
| L2 cache | 256KB | ~4ns | per-core |
| L3 cache | 8MB | ~12ns | shared |
| DRAM | GBs | ~100ns | system |

- **Cache lines are 64 bytes** — always access data in 64-byte aligned chunks. One cache miss fetches entire line.
- **Spatial locality**: access memory sequentially — adjacent elements ride same cache line for free
- **Temporal locality**: reuse data while it's still in cache — re-access within same loop iteration
- **Cache blocking/tiling**: for O(N²) or O(N³) loops, block into chunks that fit in L1/L2 — same principle as GPU tiled design

```
// Naive matrix multiply: thrashes cache for large N
for(int i=0; i<N; i++)
  for(int j=0; j<N; j++){
    float acc = 0;
    for(int k=0; k<N; k++) acc += A[i*N+k] * B[k*N+j];
    C[i*N+j] = acc;
  }

// Cache-blocked: BS chosen so A[i:i+BS, k:k+BS] fits in L1
for(int ii=0; ii<N; ii+=BS)
  for(int jj=0; jj<N; jj+=BS)
    for(int kk=0; kk<N; kk+=BS)
      for(int i=ii; i<ii+BS; i++)
        for(int j=jj; j<jj+BS; j++){
          float acc = C[i*N+j];
          for(int k=kk; k<kk+BS; k++) acc += A[i*N+k] * B[k*N+j];
          C[i*N+j] = acc;
        }
```

## Data Layout & Memory Access

- **Flat arrays over pointer-chasing**: `float* arr` with index arithmetic beats `std::vector<std::vector<float>>` or linked structures
- **Structure-of-Arrays (SoA) over Array-of-Structs (AoS)**: separate `float* x, *y, *z` — not `struct Point{float x,y,z;}* pts` — enables SIMD and coalesced cache access
- **Align to 16/32/64 bytes**: `alignas(32) float arr[N]` — enables SIMD load instructions
- **Avoid false sharing**: when multiple threads write to adjacent fields in same struct, cache lines bounce between cores. Pad to cache line boundary.

```
// Bad: false sharing — threads writing adjacent elements invalidate each other's cache lines
struct Counter { int n; };  // 4 bytes, multiple per cache line
Counter counters[NTHREADS]; // threads[i].n++ causes cache line ping-pong

// Good: pad to cache line
struct alignas(64) Counter { int n; };  // each counter gets own cache line
```

## SIMD / Vectorization

- **Prefer auto-vectorization over manual SIMD** — compilers (-O3, -march=native) are usually better at register scheduling than humans. Focus on making data layout vectorization-friendly, not writing intrinsics.
- **Vectorization-friendly data layout**: contiguous arrays, unit stride, aligned to 16/32 bytes, SoA not AoS (see Data Layout)
- **Vectorization-friendly code**: no branches in inner loop, no function calls (or `__attribute__((always_inline))`), `__restrict__` pointers to tell compiler no aliasing
- **`#pragma omp simd`**: portable hint to vectorize loop even when compiler is conservative
- **Manual SIMD (SSE/AVX/AVX-512 intrinsics) only as last resort**: when compiler refuses to vectorize despite good layout, and you've verified with compiler reports (`-fopt-info-vec`)

```
// Auto-vectorizable: contiguous, unit stride, no branches, restrict pointers
void add(float* __restrict__ a, float* __restrict__ b, float* __restrict__ c, int n){
    #pragma omp simd
    for(int i=0; i<n; i++) c[i] = a[i] + b[i];
}
// Compiler generates AVX2/AVX-512 automatically with -O3 -march=native
```

## Loop Optimization

- **Hoist invariant computations** out of loops
- **Unroll small loops**: compiler does this at -O3, but manual unroll helps when you know trip count
- **Fuse loops** with same iteration space to improve cache locality
- **Split loops** if body is too large (register pressure)
- **Avoid branches in inner loops**: replace `if` with arithmetic masks or separate loops

```
// Bad: branch in inner loop prevents vectorization
for(int i=0; i<n; i++){
    if(mask[i]) c[i] = a[i] + b[i];
    else        c[i] = a[i] - b[i];
}

// Good: branchless — arithmetic select
for(int i=0; i<n; i++){
    float sign = mask[i] ? 1.0f : -1.0f;
    c[i] = a[i] + sign * b[i];
}

// Also good: split into two loops
for(int i=0; i<n; i++) if(mask[i]) c[i] = a[i] + b[i];
for(int i=0; i<n; i++) if(!mask[i]) c[i] = a[i] - b[i];
```

## Memory Allocation

- **No heap allocation in hot paths** — `malloc`/`new` are ~100ns+ each
- **Preallocate**: allocate buffers once, reuse across iterations
- **Stack allocation**: `float buf[N]` or `alloca()` for temporary arrays — effectively free
- **Arena/pool allocators**: for many small objects of same type
- **Avoid `std::vector` resize in loops**: reserve capacity upfront, `v.reserve(n)` once

## Branch Prediction

- **Sort data by branch outcome**: if 90% of data takes same branch, sort so all true cases are contiguous — branch predictor succeeds
- **Likely/unlikely hints**: `__builtin_expect(!!(x), 1)` or `[[likely]]`/`[[unlikely]]` (C++20)
- **Eliminate branches**: arithmetic masks, lookup tables, bit tricks

## Precomputation & Sub-expression Reuse

- **Precompute small reusable quantities**: if a value is used multiple times, compute once and reuse — especially invariants across loop iterations
- **Hoist sub-expressions**: `exp(-b*x)` appears in both value and derivative — compute once, reuse

```
// Bad: exp called twice
double morse(double x, double D, double a, double r0){
    double dx = x - r0;
    return D * (exp(-2*a*dx) - 2*exp(-a*dx));
}
double morse_deriv(double x, double D, double a, double r0){
    double dx = x - r0;
    return 2*D*a * (exp(-2*a*dx) - exp(-a*dx));
}

// Good: precompute e=exp(-a*dx), reuse for value AND derivative
struct Morse { double D, a, r0; };
inline void morse_eval(const Morse& m, double x, double& E, double& F){
    double dx = x - m.r0;
    double e  = exp(-m.a * dx);   // single exp call
    E = m.D * (e*e - 2*e);
    F = 2*m.D*m.a * (e*e - e);    // reuse e, no second exp
}
```

- **Precompute lookup tables** for repeated expensive evaluations (e.g. grid of spline values)
- **Manual preload** from large arrays into local variables: similar to GPU tiled design, though CPU cache does this automatically — only helps when access pattern is non-obvious to prefetcher

## Avoid Transcendental Functions

`sin`, `cos`, `exp`, `log` cost ~20-100 cycles each. Avoid when possible:

- **Polynomial approximations**: Taylor expansion, Chebyshev, spline fits — see `cpp/common/math/fastmath.h` for existing fast approximations
- **Replace sin/cos series with complex multiplication**: rotating by angle θ repeatedly → multiply by `(cosθ + i*sinθ)` each step. No trig per iteration.

```
// Bad: recompute sin/cos each step
for(int i=0; i<N; i++){
    double c = cos(theta * i);
    double s = sin(theta * i);
    // use c, s...
}

// Good: complex multiplication — no trig after first step
double c = cos(theta), s = sin(theta);  // compute once
double re = 1.0, im = 0.0;
for(int i=0; i<N; i++){
    // use re, im...
    double nr = re*c - im*s;  // complex multiply: rotate by theta
    double ni = re*s + im*c;
    re = nr; im = ni;
}
```

- **Prefer algebraic representations over angles**: complex numbers, quaternions, vectors, rotation matrices — avoid trig entirely. Store rotation as quaternion, not Euler angles. Compose with quaternion multiplication (16 mul) instead of matrix from angles (2 sin + 2 cos + mul).
- **When transcendental is unavoidable**: precompute outside loop, or use fast polynomial approximation from `fastmath.h`

## Multi-threading (OpenMP)

- **`#pragma omp parallel for`**: trivial parallelization of counting loops
- **False sharing**: ensure each thread writes to separate cache lines (see Data Layout)
- **NUMA awareness**: first-touch policy — initialize data on thread that will use it
- **Avoid oversubscription**: don't nest OpenMP parallel regions

### OpenMP Fork/Join Latency — Restructure Algorithms

**Opening/closing `#pragma omp parallel` has significant latency** (~microseconds). Do NOT fork-join per iteration. Instead, open `#pragma omp parallel` once around the outer loop, use `#pragma omp for` (no `parallel`) inside.

This requires **restructuring serial algorithms** from per-iteration function calls to per-atom work units inside a single parallel region.

```
// BAD: fork-join every iteration (latency dominates for small N)
for(int itr=0; itr<niter; itr++){
    #pragma omp parallel for reduction(+:E)
    for(int i=0; i<n; i++) E += evalForce(i);
    #pragma omp parallel for
    for(int i=0; i<n; i++) moveAtom(i);
}

// GOOD: fork once, join once — restructure to per-atom work units
#pragma omp parallel
for(int itr=0; itr<niter; itr++){
    #pragma omp for reduction(+:E)
    for(int i=0; i<n; i++) E += evalSingleAtom(i);
    // implicit barrier after for — all atoms evaluated before move

    #pragma omp for
    for(int i=0; i<n; i++) moveSingleAtom(i);
    // implicit barrier — all atoms moved before next iteration
}
```

### Assembly Step Pattern (Avoid Atomics)

Same as GPU: when forces from multiple interactions write to the same atom, use auxiliary buffers + assembly kernel instead of atomics.

```
// Phase 1: each thread evaluates forces for its atoms, writes to per-atom buffers
#pragma omp for reduction(+:E)
for(int ia=0; ia<natoms; ia++){
    E += evalSingleAtom(ia);  // writes to fapos[ia] (thread-local, no race)
}
// implicit barrier

// Phase 2: assemble — gather contributions, apply constraints, move
#pragma omp for
for(int ia=0; ia<natoms; ia++){
    assemble_atom(ia);  // collect neighbor recoils, apply constraints
}
// implicit barrier
```

Real example: `MolWorld_sp3::run_omp` — serial `evalAllForces(); moveAllAtoms()` restructured to `#pragma omp parallel` with per-atom `evalSingleAtom` + `assemble_atom` + `moveSingleAtom` phases.

```
// Simple parallel loop (one-shot, not iterative)
#pragma omp parallel for schedule(static)
for(int i=0; i<n; i++) c[i] = a[i] + b[i];

// Reduction
#pragma omp parallel for reduction(+:sum)
for(int i=0; i<n; i++) sum += a[i] * b[i];
```


## Related Skills

- skill:`gpu-optimize` — GPU-specific optimization (tiling, local memory, atomics avoidance)
- skill:`python-perf` — Python performance (vectorization, NumPy, when Python loops are OK)
- skill:`cpp-build` — build system, ASAN vs opt modes
