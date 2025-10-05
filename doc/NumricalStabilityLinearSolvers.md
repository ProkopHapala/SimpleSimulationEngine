# Iterative Solvers for Truss Dynamics: Jacobi and Gauss–Seidel on displacements

## 1. Motivation

When simulating truss or spring-mass systems with implicit integration (e.g. implicit Euler), we must solve a large linear system at every timestep:

$$
A p = b
$$

where:

* $p$ is the vector of unknown vertex displacements,
* $A$ is the **system matrix** containing both inertia and stiffness:
  $$
  A = \frac{M}{\Delta t^2} + K,
  $$
  with $M$ the diagonal mass matrix (scaled by timestep $\Delta t$) and $K$ the stiffness matrix assembled from springs,
* $b$ is the **right-hand side (RHS)** containing inertial prediction and external forces:
  $$
  b = \frac{M}{\Delta t^2} \big(p_0 + v_0 \Delta t + f_{ext} \Delta t^2 \big) + f_{int},
  $$
  where $p_0$ and $v_0$ are reference position and velocity, $f_{ext}$ are external forces (e.g. gravity), and $f_{int}$ are internal spring contributions.

### 1.1 Update loop in C++

Implicit Euler integration loop with linear solver. Velocity is taken as $v = \frac{p_{k+1}-p_{k}}{\Delta t}$.

```C++
    void run_LinSolve(int niter) {
        memcpy(ps_cor, points, nPoint * sizeof(Vec3d)); // backup current positions
        double dt2    = dt * dt;
        double inv_dt = 1/dt;
        double cdamp  = fmax( 1-damping*dt, 0. );
        for (int iter = 0; iter < niter; iter++) {
            getForce( forces, points, vel, dt );
            //  Predictor step
            for (int i=0;i<nPoint;i++){ 
                ps_pred[i] = points[i].f + vel[i].f*dt + forces[i].f*(dt2/points[i].w);  // inertial prediction
            }
            // Apply fixed constraints
            if(kFix) for (int i = 0; i < nPoint; i++) { if (kFix[i] > 1e-8 ) { ps_pred[i] = points[i].f; } }
            // solve linear system for implicit euler step
            linsolve( ps_pred, ps_cor );            
            //   Corrector step
            double l2sum = 0.0;
            for (int i=0;i<nPoint;i++) {
                Vec3d v   = ps_cor[i] - points[i].f;
                v.mul(inv_dt);
                v.mul( cdamp );  // damp velocity
                if(kFix){ if( kFix[i] > 1e-8){ v=Vec3dZero; } }
                vel[i].f    = v;
                points[i].f = ps_cor[i];  // update positions
            }            
            time += dt;
        } // for iter ( time stepping )
    }
```

### 1.2 Iterative Linear solver

```C++
    void TrussDynamics_d::updateIterativeMomentum( Vec3d* psa, Vec3d* psb ){
        updatePD_RHS(psa, bvec );
        for (int i = 0; i < nSolverIters; i++) {  
            switch( (LinSolveMethod)linSolveMethod ){
                case LinSolveMethod::JacobiMomentum:    { updateJacobi_lin( psa, psb, bvec ); } break;
                case LinSolveMethod::JacobiFlyMomentum: { updateJacobi_fly( psa, psb );       } break;
                case LinSolveMethod::GSMomentum:        { 
                    for (int j=0; j<nPoint; j++){ psb[j]=psa[j]; }
                    updateGaussSeidel_lin( psb, bvec ); 
                } break; 
                case LinSolveMethod::GSFlyMomentum:     { 
                    for (int j=0; j<nPoint; j++){ psb[j]=psa[j]; }
                    updateGaussSeidel_fly( psb ); 
                } break; 
            }
            double bmix = mixer.get_bmix( i );                 // bmix_k = bmix(k)
            if( (i==0) || (i>=(nSolverIters-1) ) ) bmix = 0.0; // make sure we stop momentum mixing for the last iteration      
            for(int j=0; j<nPoint; j++){ 
                Vec3d p = psb[j]; 
                p.add_mul( interia[j], bmix );   //    p_{k+1} = p'_k + bmix d_k
                interia[j] = p - psa[j];         //    d_{k+1} = p_{k+1} - p_k
                psa[j]         = p;
            }
        }
        for (int i=0; i<nPoint; i++){ psb[i]=psa[i]; } // final update
    }
```

Direct factorization of $A$ is too expensive for GPU-based, interactive, or large-scale problems. Instead, we use **iterative relaxation methods**: Jacobi or Gauss–Seidel. These methods update vertex displacements locally and converge progressively.

## 2. Working with Displacements

To avoid numerical cancellation, we always operate on **relative displacements** rather than absolute positions:

$$
dp = p - p_0
$$

* Computing $dp$ is **cheap**, but it is **precision-critical** because it isolates small differences from potentially huge world coordinates.
* This subtraction and other geometry-sensitive steps (like edge lengths, inertial prediction, or time integration with Verlet/Leapfrog) are best done in **double precision**.
* The expensive, highly parallel iterative solver can then safely run in **single precision** on GPU using $dp$.

## 3. Jacobi vs. Gauss–Seidel

### Jacobi

* Each vertex update uses **only the values from the previous iteration**.
* Updates can be computed in parallel for all vertices.
* Convergence is usually slow.

### Gauss–Seidel

* Each vertex update uses the **most recently updated neighbors** (in-place update).
* Converges faster than Jacobi.
* To parallelize safely, we need **graph coloring**: partition vertices into independent sets ("colors") such that no two neighbors share the same color. Vertices of one color can be updated in parallel, then we synchronize before moving to the next color.

## 4. Diff vs. Fly Variants

Both Jacobi and Gauss–Seidel can be implemented in two flavors:

### 4.1. Diff (Linear Precomputed RHS)

* Before iterations, compute per-vertex coefficients:

  * $b_i$ = RHS vector
  * $A_{ii}$ = diagonal entry of system matrix
* These are typically computed in **double precision** on CPU or in a GPU pre-pass, then cast to float.
* Iterative update rule (Jacobi form):
  $$
  dp_i^{new} = \frac{b_i + \sum_{j \in \mathcal{N}(i)} k_{ij} dp_j}{A_{ii}}
  $$
* Pros:

  * Cheaper per iteration (just sums and divisions).
  * Stable if precomputation is accurate.
* Cons:

  * Only correct for the linearized system around the reference $p_0$.
  * If deformation is large, may converge slowly.

### 4.2. Fly (Nonlinear Recomputed RHS)

* RHS $b_i$ and diagonal $A_{ii}$ are recomputed **at every iteration** from the current displacements.
* Includes nonlinear terms like:
  $$
  d_{ij} = \frac{p_i - p_j}{|p_i - p_j|}, \quad b_i \sim k_{ij} \left(\frac{l_0}{|d_{ij}|} - 1\right)(p_i - p_j)
  $$
* Pros:

  * Adapts to nonlinearity automatically.
  * No need for a separate precomputation kernel.
  * Often more robust for large deformations or chaotic motion.
* Cons:

  * More expensive per iteration (norms, divisions).
  * Sensitive to numerical precision (requires clamping and sometimes double precision for lengths).
  * No longer a strict linear solver; behaves more like nonlinear fixed-point iteration.

## 5. Numerical Stability and Precision Considerations

* **Edge length normalization:** avoid division by zero. Use
  $$
  \text{len} = \max(|p_i - p_j|, \epsilon)
  $$
  with $\epsilon \approx 10^{-9}$–$10^{-12}$.
* **Large coordinates:** always work in **relative displacements** to keep magnitudes small.
* **Double precision use cases:**

  * Computing edge lengths and reciprocals (for Fly variant).
  * Precomputing RHS and diagonal (for Diff variant).
* **Regularization:** add a small diagonal shift
  $$
  A_{ii} \leftarrow A_{ii} + \lambda
  $$
  to prevent instabilities when $A_{ii}$ is very small.
* **Jacobi vs GS:** Jacobi is fully parallel but converges slowly; GS requires coloring but converges faster.

## 6. Mapping to GPU Kernels

* **Local memory:** load displacements into local memory once, iterate there, then write back.
* **Loop inside kernel:** avoid kernel launch overhead by doing all iterations within a single kernel.
* **Coloring for GS:** precompute color partitions on host; pass `color_list` and `color_offsets` arrays to kernel.
* **CSR neighbor storage:** use compressed sparse row format (`nbr_off`, `nbr_len`, `nbrs`) for neighbor traversal.
* **Two versions per method:**

  * **Diff kernels:** simple arithmetic, lighter, rely on precomputed $b_i, A_{ii}$.
  * **Fly kernels:** heavier, recompute per iteration, but better adapt to nonlinear deformation.

## 7. Summary

* **Jacobi-Diff**: simple, fully parallel, converges slowly, but cheap.
* **Jacobi-Fly**: parallel, more robust for nonlinear motion, heavier per iteration.
* **GS-Diff**: converges faster, needs coloring, still cheap per iteration.
* **GS-Fly**: fastest convergence in practice for nonlinear deformation, but most expensive per iteration.

These four solvers form a useful spectrum of trade-offs between speed per iteration and robustness to large deformations. In practice:

* Use **GS-Diff** for small, well-conditioned problems.
* Use **GS-Fly** when deformations are large or unpredictable.
* Compare against Jacobi versions as baselines for performance and accuracy.

# Jacobi solvers in OpenCL

## Key translation points

1. **Jacobi vs Gauss–Seidel**

   * Jacobi is naturally parallel: each vertex update uses **only the previous iterate** `ps_in` and writes to `ps_out`.
   * So we don’t need graph coloring, no data hazards. That’s why it’s the baseline for comparison.

2. **Displacements, not absolute positions**

   * As you said: operate on displacements `dp = p - p0`. That way we don’t lose precision when positions are large.
   * In code: you keep a double-precision buffer `p0` on CPU/GPU. Jacobi runs on `dp` in single-precision (`float3`/`float4`).

3. **Right-hand side preparation (`bvec`)**

   * As in your C++: the expensive parts (like `updatePD_dRHS`) should be computed in **double precision** before the Jacobi iterations, then passed in as `bvec`.
   * `bvec[i] = (f_i, Aii)` stores RHS vector and diagonal coefficient. That allows GPU solver to stay in float.

4. **Neighbor sums**

   * Each thread accumulates $\sum_j k_{ij} dp_j$ for its vertex. That’s the Jacobi off-diagonal contribution.

5. **Loop structure**

   * Outer loop over iterations stays inside kernel (to avoid kernel enqueue overhead).
   * At each iter:

     * All threads read `ps_in` (local copy).
     * Compute new `ps_out[i]`.
     * Swap buffers.

6. **Local memory**

   * Since the whole system is small (128–512 verts), we can keep both `ps_in` and `ps_out` in local memory.
   * Only global write at end.

## OpenCL kernel (Jacobi iterations in local memory)

```c
// Jacobi iterative solver for truss dynamics (on displacements, not absolute positions).
// Single workgroup handles the entire system (up to ~512 verts).
// Inputs:
//   - dp_global : initial displacements (float4, w unused) [nverts]
//   - bvec      : RHS array (float4, xyz = b_i, w = Aii) [nverts]
//   - neigh_off, neigh_len, neighs, neigh_k : CSR-like neighbor storage
// Config:
//   - nIters    : number of Jacobi iterations

inline float3 to3(float4 v) { return (float3)(v.x,v.y,v.z); }
inline float4 to4(float3 v) { return (float4)(v.x,v.y,v.z,0.0f); }

__kernel void jacobi_truss(
    __global float4* dp_global,      // in/out displacements
    __global const float4* bvec,     // RHS: xyz = b_i, w = Aii
    __global const int* neigh_off,
    __global const int* neigh_len,
    __global const int* neighs,
    __global const float* neigh_k,   // stiffness per neighbor entry
    const int nverts,
    const int nIters
){
    const int lid = get_local_id(0);
    const int lsize = get_local_size(0);

    __local float4 dp_in[1024];
    __local float4 dp_out[1024];

    // Load initial displacement into local buffer
    for (int i = lid; i < nverts; i += lsize) {
        dp_in[i] = dp_global[i];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    for (int iter=0; iter<nIters; ++iter) {
        for (int i = lid; i < nverts; i += lsize) {
            float3 sum_j = (float3)(0.0f);
            int no = neigh_off[i];
            int nl = neigh_len[i];
            for (int jj=0; jj<nl; ++jj) {
                int idx = no + jj;
                int j   = neighs[idx];
                float k = neigh_k[idx];
                float3 dpj = to3(dp_in[j]);
                sum_j += k * dpj;
            }
            float4 bi = bvec[i];
            float3 rhs = to3(bi) + sum_j;
            float Aii = bi.w;
            float3 new_dp = rhs * (1.0f / Aii);
            dp_out[i] = to4(new_dp);
        }
        barrier(CLK_LOCAL_MEM_FENCE);

        // Swap dp_in/out
        for (int i = lid; i < nverts; i += lsize) {
            dp_in[i] = dp_out[i];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // Write back final displacement
    for (int i = lid; i < nverts; i += lsize) {
        dp_global[i] = dp_in[i];
    }
}
```

### Commentary

* **RHS (`bvec`)**: This must be computed beforehand (CPU or GPU, in double precision). Each entry packs $(b_i, A_{ii})$. That makes the Jacobi loop lightweight.
* **Local memory size**: I declared `[1024]` — adjust for your max system size. Host must ensure `nverts <= local_size_limit`.
* **Single WG assumption**: For your target sizes (128–512 vertices, WG=32), a single workgroup can handle all vertices easily.
* **Numerical note**: Since we iterate on displacements, accuracy is good even with `float`. The double-precision work (RHS assembly) ensures correct setup.

### What the “fly” variant does

* Normally, in PD/VBD, you **precompute** the RHS (`bvec`) and diagonal (`Aii`) in double precision, then run iterative solver in float.
* The **fly variant** skips that precomputation:

  * At each iteration, each vertex thread recomputes its own `b_i` and `Aii` from the current displacement/position.
  * This includes nonlinear terms like $d_{ij} = \frac{p_i - p_j}{|p_i - p_j|}$ (normalized edge direction).

This means:

* It’s **closer to a nonlinear Gauss–Newton scheme**, because $b$ and $A$ are recomputed each iteration instead of being frozen.
* The solver does not strictly solve the linearized system $A dp = b$ anymore; instead, it’s “self-updating” the RHS each iteration.


### Recommendation

* If positions are stored as absolute coordinates (large values, like 1e+6) → do the `|p_i - p_j|` and normalization in **double**.
* If you instead work with **displacements relative to reference** (values ~1e-1 .. 1e2) → single precision is fine.
* Either way: add a small epsilon in denominator (`1e-8f` or so).

## OpenCL kernel for Jacobi “fly” variant

This kernel recomputes $b_i$ and $Aii$ on the fly each iteration.

```c
// Jacobi "fly" solver: recompute RHS (b_i) and diagonal (Aii) each iteration.
// Operates directly on displacements dp[] in local memory.
// Each iteration: 
//   dp_new[i] = ( M_i/dt^2 * p'_i + sum_j [ spring terms ]) / Aii
//
// Inputs:
//   - dp_global : initial displacements (float4, w unused) [nverts]
//   - rest_len  : rest lengths per edge (float) [nedges]
//   - stiff     : stiffness per edge (float) [nedges]
//   - neigh_off, neigh_len, neighs : neighbor CSR
//   - points_w  : mass per vertex (points[i].w)
// Config:
//   - dt        : timestep
//   - nIters    : iterations
//
// Assumes: system fits in one workgroup (<= 512 verts).

inline float3 to3(float4 v) { return (float3)(v.x,v.y,v.z); }
inline float4 to4(float3 v) { return (float4)(v.x,v.y,v.z,0.0f); }

__kernel void jacobi_truss_fly(
    __global float4* dp_global,      // in/out displacements
    __global const float* points_w,  // vertex masses
    __global const int* neigh_off,
    __global const int* neigh_len,
    __global const int* neighs,
    __global const float* stiff,     // stiffness per neighbor entry
    __global const float* rest_len,  // rest length per neighbor entry
    const int nverts,
    const float dt,
    const int nIters
){
    const int lid = get_local_id(0);
    const int lsize = get_local_size(0);

    __local float4 dp_in[1024];
    __local float4 dp_out[1024];

    // Load initial displacement into local buffer
    for (int i = lid; i < nverts; i += lsize) {
        dp_in[i] = dp_global[i];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    float inv_dt2 = 1.0f / (dt*dt);

    for (int iter=0; iter<nIters; ++iter) {
        for (int i = lid; i < nverts; i += lsize) {
            float3 pi = to3(dp_in[i]);

            // inertial contribution: M_i/dt^2 * p_i
            float Ii = points_w[i]*inv_dt2;
            float3 bi = pi * Ii;
            float Aii = Ii;

            int no = neigh_off[i];
            int nl = neigh_len[i];
            for (int jj=0; jj<nl; ++jj) {
                int idx = no + jj;
                int j   = neighs[idx];
                float3 pj = to3(dp_in[j]);
                float k   = stiff[idx];
                float l0  = rest_len[idx];

                float3 dij = pi - pj;
                // compute length in double for stability
                double dx = (double)dij.x;
                double dy = (double)dij.y;
                double dz = (double)dij.z;
                double len2 = dx*dx + dy*dy + dz*dz;
                double len  = sqrt(len2);
                if(len < 1e-12) len = 1e-12;

                float invlen = (float)(1.0/len);

                // Add RHS contributions
                // bi += k*( l0/len - 1 ) * (pi - pj)
                float coeff = k * ( (float)(l0*invlen) - 1.0f );
                bi += coeff * dij;

                // bi += k * pj
                bi += k * pj;

                Aii += k;
            }

            // Jacobi update
            float invA = 1.0f / Aii;
            float3 new_dp = bi * invA;
            dp_out[i] = to4(new_dp);
        }
        barrier(CLK_LOCAL_MEM_FENCE);

        // Swap
        for (int i = lid; i < nverts; i += lsize) {
            dp_in[i] = dp_out[i];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // Write back final
    for (int i = lid; i < nverts; i += lsize) {
        dp_global[i] = dp_in[i];
    }
}
```

### Behavior summary

* **Validity**: Yes, the “fly” variant is a valid substitute. It’s no longer a strict linear solver, but rather an iterative nonlinear scheme. In practice it often behaves like a damped Newton method.
* **Numerical stability**:

  * Sensitive to very short edges (`len → 0`). Must clamp denominator.
  * If positions are large, length computation should be done in double (as above).
* **Performance**: More arithmetic per iter than precomputed Jacobi, but saves memory and handles nonlinearity better.

# Gauss–Seidel in OpenCL

## Goal

Present the design choices and motivations for implementing iterative solvers for truss (spring/edge) dynamics on GPUs. We assume implicit time stepping (implicit Euler) so we solve a linear(ized) system at each time step. The solver works on **displacements** (differences from a reference $p_0$) to avoid catastrophic cancellation when world coordinates are large.

We present two GPU-ready Gauss–Seidel (GS) variants (coloring-based):

1. **GS (linear / precomputed)** — the standard linear Jacobi/Gauss–Seidel where the RHS $b$ and diagonal $A_{ii}$ are precomputed (preferably in double precision), then iterations run in single precision on displacements.
2. **GS (fly / nonlinear)** — a variant that recomputes local RHS and diagonal every iteration from current iterates (nonlinear / “on-the-fly”). This often better handles large deformations or avoids precomputation cost, but it has different stability/accuracy trade-offs.

These are compared to Jacobi kernels (baseline) and a vertex-block Newton (VBD) solver (more work per vertex but faster convergence).

## Why Gauss–Seidel + coloring?

* GS uses the most recent updates (in-place), so it usually converges faster than Jacobi for the same number of local solves.
* Coloring allows parallel GS: all vertices of a color are mutually independent (no edge between them), so they can be updated simultaneously.
* For chain-like structures, choose coloring or ordering that matches topology (BFS from anchor) to propagate influence quickly.

## High-level algorithm (Gauss–Seidel color-by-color)

* Precompute a vertex coloring such that vertices of the same color are not neighbors (no direct springs between same-color vertices). All vertices of a color can be updated in parallel without data hazard.
* Copy displacements ($dp$) into local memory for fast access.
* Repeat for `nIters`:

  * For each color `c`:

    * Threads in the workgroup process distinct indices from color `c` (work assigned by local id with stride).
    * For each vertex `i` processed:

      * Accumulate gradient/contributions from **all** neighbors using the *current* (local) $dp$ values.
      * For linear variant: use precomputed `b_i` and `Aii`; compute `dp_i = (b_i + \sum_j k_{ij} dp_j)/Aii`.
      * For fly variant: recompute `b_i` and `Aii` from geometry using current $dp$, then compute update.
      * Write the new `dp_i` into local memory (in-place), so subsequent threads (other colors) see the updated value for later colors.
    * Barrier to ensure all updates for color `c` are visible.
* Write local `dp` back to global memory.

## Design goals for GPU kernels

* **Use local memory** for the working displacement arrays — reduces global memory traffic inside inner loops.
* **Minimize branching** in the hot loop; use contiguous color lists and a compact CSR neighbor list.
* **Compute the full per-vertex Hessian/diagonal (or Aii) and gradient before solving** — do not solve using partial information.
* **Keep per-vertex work balanced**: each thread should handle several vertices if `color_size > local_size`.
* **Add small regularizer** to diagonals when necessary (prefer `Aii += reg` over skipping updates).
* **Precision policy**:

  * Precompute `b_i` and `Aii` in **double** if input positions are large or if you compute geometry terms relative to absolute coordinates.
  * Iterations on `dp` can run in single precision (float) for throughput.
  * For fly variants, compute norms/reciprocals in double if coordinates have large magnitude or if you want extra stability.

## Numerical stability and operations in double precision

* Danger zones:

  * Division by very small edge length during normalization. Use `safe_len = max(len, eps)`; consider computing `len` in double.
  * Subtraction of large similar numbers when using absolute positions — that's why we iterate on displacements.
* Recommendation:

  * If host can compute `b_i` and `Aii` in double (RHS assembly), do it and pass floats to the kernel.
  * If using fly variant, compute length and reciprocal in double inside kernel (`double dx = (double)...`) then cast result to float for accumulation.
* Regularization:

  * Use `Aii += reg_factor * max(1.0f, trace(H))` in VBD or simply `Aii += reg_const` before dividing to avoid blow-ups.

## When to use which variant

* **Linear precomputed GS**: prefer when you can assemble `b`/`Aii` cheaply in double or on CPU. Faster per iteration and predictable.
* **Fly GS**: prefer if assembling `b` is expensive or you need to adapt to nonlinearity dynamically. More robust for large instantaneous deformations, but slower per iteration and heavier arithmetic.
* **VBD (local Newton)**: prefer when convergence speed per iteration is primary; more complex but fewer iterations overall.

## Data-layout reminders (host preprocess)

* Precompute and send to GPU:

  * `color_list[]` (concatenated vertex indices per color).
  * `color_offsets[]` (length `num_colors+1`).
  * CSR-like neighbor arrays: `nbr_off[i]`, `nbr_len[i]`, `nbrs[idx]`.
  * `nbr_k[idx]` and `nbr_l0[idx]` aligned with `nbrs`.
  * For linear variant: `bvec[i]` as `float4` where `.xyz = b_i` and `.w = Aii`.
  * Masses scaled by `1/dt^2` optionally pre-multiplied.
* Ensure per-color vertex lists are contiguous for branchless indexing.

## kernel Gauss-Seidel (linear diff)

```c
// Gauss-Seidel (linear diff) kernel using coloring, operates on displacements dp (dp = p - p0).
// Precondition: host computed bvec[i] (float4) where .xyz = b_i, .w = Aii (preferably in double, then cast to float)
// Data layout / inputs:
//   pos_dp_global : float4 dp (in/out)   [nverts]
//   bvec_global   : float4 (b,Aii)       [nverts]
//   color_list, color_offsets           // concatenated list per color
//   nbr_off, nbr_len, nbrs, nbr_k       // CSR neighbor storage aligned
//   num_colors, nIters, det_eps, reg_eps
//
// Assumptions:
//   - Entire vertex set assigned to single workgroup (nverts <= LOCAL_MAX).
//   - color_list holds contiguous blocks for each color: indices = color_list[color_offsets[c] .. color_offsets[c+1]-1]
//   - mass/dt^2 contributions already folded into bvec by host (if desired).

#pragma OPENCL EXTENSION cl_khr_fp64 : enable

inline float3 to3(float4 v) { return (float3)(v.x, v.y, v.z); }
inline float4 to4(float3 v) { return (float4)(v.x, v.y, v.z, 0.0f); }

__kernel void gs_linear_colors_dp(
    __global float4* dp_global,       // in/out displacements (dp = p - p0)
    __global const float4* bvec_global, // .xyz = b_i, .w = Aii
    __global const int* color_list,
    __global const int* color_offsets, // length num_colors+1
    const int num_colors,
    __global const int* nbr_off,
    __global const int* nbr_len,
    __global const int* nbrs,
    __global const float* nbr_k,       // stiffness aligned with nbrs
    const int nverts,
    const int nIters,
    const float reg_eps)               // regularizer scale
{
    const int lid = get_local_id(0);
    const int lsize = get_local_size(0);

    // local buffers (adjust size if needed; host must ensure nverts <= LOCAL_MAX)
    __local float4 loc_dp[1024];

    // load initial dp into local memory
    for (int i = lid; i < nverts; i += lsize) {
        loc_dp[i] = dp_global[i];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // main GS loop
    for (int iter = 0; iter < nIters; ++iter) {
        // sweep colors
        for (int c = 0; c < num_colors; ++c) {
            int off = color_offsets[c];
            int off_next = color_offsets[c+1];
            int color_count = off_next - off;

            // each thread handles multiple vertices in this color
            for (int local_idx = lid; local_idx < color_count; local_idx += lsize) {
                int vid = color_list[off + local_idx];

                // load current dp and bvec
                float3 dpi = to3(loc_dp[vid]);
                float4 bi4 = bvec_global[vid];
                float3 b_i = to3(bi4);
                float Aii = bi4.w;

                // accumulate sum_j k_ij * dp_j (current dp values in loc_dp allow GS)
                float3 sum_j = (float3)(0.0f, 0.0f, 0.0f);
                int no = nbr_off[vid];
                int nl = nbr_len[vid];
                for (int e = 0; e < nl; ++e) {
                    int idx = no + e;
                    int j = nbrs[idx];
                    float k = nbr_k[idx];
                    float3 dpj = to3(loc_dp[j]); // uses the latest dp (GS behavior)
                    sum_j.x += k * dpj.x;
                    sum_j.y += k * dpj.y;
                    sum_j.z += k * dpj.z;
                }

                // compute new dp: (b_i + sum_j) / Aii
                // regularize Aii slightly if needed
                float reg = reg_eps * fmax(1.0f, Aii);
                float invA = 1.0f / (Aii + reg);
                float3 new_dp = (b_i + sum_j) * invA;

                // store in local buffer (in-place GS)
                loc_dp[vid] = to4(new_dp);
            } // per-color per-thread loop

            // synchronize all threads before next color
            barrier(CLK_LOCAL_MEM_FENCE);
        } // colors
    } // iters

    // write back final dp
    for (int i = lid; i < nverts; i += lsize) {
        dp_global[i] = loc_dp[i];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
}
```

## kernel Gauss-Seidel (fly diff)

```c
// Gauss-Seidel (fly) kernel using coloring, operates on displacements dp (dp = p - p0).
// This recomputes local RHS and diagonal Aii each iteration from current dp (non-linear).
// Inputs:
//   dp_global : float4 dp (in/out)
//   points_w  : mass per vertex (float) OR mass*inv_h2 if pre-scaled on host
//   color_list, color_offsets
//   nbr_off, nbr_len, nbrs, nbr_k, nbr_l0
//   num_colors, nIters, dt
// Notes:
//   - Some geometry ops (lengths, inv len) are done in double for stability.
//   - For small edge lengths, we clamp length to eps to avoid div-by-zero.
//   - Regularization of diagonal is applied similarly.

#pragma OPENCL EXTENSION cl_khr_fp64 : enable

inline float3 to3(float4 v) { return (float3)(v.x, v.y, v.z); }
inline float4 to4(float3 v) { return (float4)(v.x, v.y, v.z, 0.0f); }

__kernel void gs_fly_colors_dp(
    __global float4* dp_global,       // in/out displacements
    __global const float* points_w,   // mass per vertex (or mass/dt^2 pre-scaled)
    __global const int* color_list,
    __global const int* color_offsets,
    const int num_colors,
    __global const int* nbr_off,
    __global const int* nbr_len,
    __global const int* nbrs,
    __global const float* nbr_k,
    __global const float* nbr_l0,
    const int nverts,
    const float dt,
    const int nIters,
    const float reg_eps)
{
    const int lid = get_local_id(0);
    const int lsize = get_local_size(0);

    __local float4 loc_dp[1024];

    // load initial dp into local
    for (int i = lid; i < nverts; i += lsize) {
        loc_dp[i] = dp_global[i];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    float inv_dt2 = 1.0f / (dt*dt);
    const double EPS_LEN = 1e-12;

    for (int iter = 0; iter < nIters; ++iter) {
        for (int c = 0; c < num_colors; ++c) {
            int off = color_offsets[c];
            int off_next = color_offsets[c+1];
            int color_count = off_next - off;

            for (int local_idx = lid; local_idx < color_count; local_idx += lsize) {
                int vid = color_list[off + local_idx];

                // read local dp
                float3 pi = to3(loc_dp[vid]);

                // inertial contribution
                float Ii = points_w[vid] * inv_dt2; // if host provided mass, else pre-scale
                float3 bi = (float3)(pi.x * Ii, pi.y * Ii, pi.z * Ii);
                float Aii = Ii;

                int no = nbr_off[vid];
                int nl = nbr_len[vid];
                for (int e = 0; e < nl; ++e) {
                    int idx = no + e;
                    int j = nbrs[idx];
                    float k = nbr_k[idx];
                    float l0 = nbr_l0[idx];

                    float3 pj = to3(loc_dp[j]);
                    float3 dij = (float3)(pi.x - pj.x, pi.y - pj.y, pi.z - pj.z);

                    // compute length in double
                    double dx = (double)dij.x;
                    double dy = (double)dij.y;
                    double dz = (double)dij.z;
                    double len2 = dx*dx + dy*dy + dz*dz;
                    double len = sqrt(len2);
                    if (len < EPS_LEN) len = EPS_LEN;
                    double invlen = 1.0 / len;

                    // coefficient and contributions
                    float coeff = (float)(k * ((double)l0 * invlen - 1.0));
                    bi.x += coeff * dij.x;
                    bi.y += coeff * dij.y;
                    bi.z += coeff * dij.z;

                    // also the k * pj term used in the original code's "fly" update
                    bi.x += k * pj.x;
                    bi.y += k * pj.y;
                    bi.z += k * pj.z;

                    Aii += k;
                } // neighbors

                // regularize diagonal slightly
                float reg = reg_eps * fmax(1.0f, Aii);
                float invA = 1.0f / (Aii + reg);

                float3 new_dp = (float3)(bi.x * invA, bi.y * invA, bi.z * invA);

                // write back updated dp (in-place GS)
                loc_dp[vid] = to4(new_dp);
            } // local_idx

            // ensure all updates for this color are visible before next color
            barrier(CLK_LOCAL_MEM_FENCE);
        } // colors
    } // iters

    // write back
    for (int i = lid; i < nverts; i += lsize) {
        dp_global[i] = loc_dp[i];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
}
```

