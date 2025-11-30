# Design Document: GPU Particle-In-Cell Accumulation Algorithms
**Application:** Steady-state Fusion/Plasma Simulation (Cylindrical Symmetry)
**Constraint:** High-Resolution Grid, Limited Particle Count ($N \approx N_{threads}$), Bandwidth-Bound.

---

## 1. Problem Statement & Core Philosophy
The objective is to compute a noise-free, time-averaged steady-state density or potential map from a stochastic particle simulation.
*   **Bottleneck:** Global Memory Bandwidth. Writing to a high-res global grid every time-step saturates the bus and causes race conditions.
*   **Goal:** Maximize arithmetic intensity by keeping data in **Registers** (L0) and **Shared Memory** (L1) as long as possible.
*   **Approach:** We reject "Scatter" (random global writes) in favor of two distinct strategies based on particle dynamics.

---

## 2. Mathematical Representation (The Basis Choice)

Before defining the kernel architecture, we must define how a "cloud" or "group" of charge is represented in memory.

### The Verdict: Local Cartesian Grid (Bilinear/Bicubic)
We compare three options for representing local charge distribution:
1.  **Multipole Expansion ($Y_{lm}$ / Taylor):** Rejected. While efficient for Far-Field, it is numerically unstable in the Near-Field (inside the cloud). It cannot robustly represent arbitrary shapes (e.g., dumbbells) without high-order terms that diverge at $r \to 0$.
2.  **Complex Meshes (FCC/BCC):** Rejected. While isotropically superior, the $O(1)$ indexing logic of Cartesian grids outweighs the minor isotropy gains on a GPU.
3.  **Local Cartesian Patch (Selected):** A small dense grid (e.g., $16 \times 16$) stored in Shared Memory.
    *   **Method:** Bilinear interpolation of particles onto this local patch.
    *   **Math:** Equivalent to a grid of overlapping Gaussian/Triangular basis functions.
    *   **Benefit:** Unconditionally stable, trivial to implement, handles any particle distribution shape, and allows fast $O(1)$ memory access.

---

## 3. Strategy A: The "Coherent Cloud" (Local Patch)
**Use Case:** Particles move slowly or remain spatially grouped (laminar flow, tight beams).

### Architecture
We treat a CUDA Workgroup (e.g., 32-256 threads) as a coherent cloud.
1.  **Initialization:** The Workgroup claims a fixed rectangular region (Patch) in **Shared Memory** that covers the expected trajectory of the group for `nSteps`.
2.  **Integration Loop:** Particles evolve in **Registers**. At each sub-step, they `atomicAdd` their charge to the **Shared Memory Patch**.
3.  **Projection:** Once the loop finishes, the aggregated Shared Memory Patch is added to the Global Grid *once*.

### Trade-offs
*   **Pros:** Minimal global bandwidth (write 1 patch per 100 steps). High arithmetic intensity.
*   **Cons:** Fails if particles diverge (ergodicity). If a particle leaves the patch bounds, its contribution is lost or requires a fallback slow-path.
*   **Optimization:** The patch is **Eulerian** (fixed in space for the duration of the kernel) to ensure correct time-integration without smearing history.

### Pseudocode (Kernel 1)
```cpp
// Shared Memory: A small local grid (e.g., 16x16 float)
__shared__ float sm_patch[16][16]; 

void kernel_coherent_cloud(...) {
    // 1. Initialize Patch
    // Calculate bounding box of this workgroup's particles
    float2 patch_origin = ...; 
    // Zero out shared memory
    for(int i=tid; i < 256; i+=blockDim.x) sm_patch[i] = 0.0f;
    __syncthreads();

    // 2. Physics Loop (In Registers)
    float2 pos = load_particle();
    for(int step = 0; step < n_steps; step++) {
        update_verlet(pos, ...);

        // 3. Accumulate to Local Patch (Bilinear Splat)
        float u = (pos.x - patch_origin.x) * inv_dx;
        float v = (pos.y - patch_origin.y) * inv_dy;
        
        if (in_bounds(u, v)) {
            // Fast Shared Memory Atomics
            atomicAdd_shared(&sm_patch[int(v)][int(u)],   weight_tl * q);
            atomicAdd_shared(&sm_patch[int(v)][int(u)+1], weight_tr * q);
            // ... (bl, br)
        }
    }

    __syncthreads();

    // 4. Flush Patch to Global Grid
    // Each thread takes one pixel of the patch
    for(int i=tid; i < 256; i+=blockDim.x) {
        if(sm_patch[i] > 0) {
           int global_idx = map_patch_to_global(i, patch_origin);
           atomicAdd_global(&global_grid[global_idx], sm_patch[i]);
        }
    }
}
```

---

## 4. Strategy B: The "Trajectory Sorter" (Deferred Splatting)
**Use Case:** Particles are Ergodic (High temperature, fast mixing, chaotic).

### Architecture
We decouple simulation from accumulation using a **Scatter-Gather** approach.
1.  **Emitter Kernel:** Simulates particles. Instead of accumulating, it simply **streams** position samples (every $k$ steps) to a massive linear Global Buffer (`pos`, `q`).
2.  **Sort Phase:** We use `cub::DeviceRadixSort` to sort these millions of anonymous samples based on their **Grid Cell Index**.
3.  **Rasterizer Kernel:** We launch thread blocks mapped to **Grid Tiles** (e.g., $16 \times 16$ global pixels). Each tile loads its relevant range of *sorted* samples from the buffer, accumulates them in Shared Memory, and writes to the grid.

### Trade-offs
*   **Pros:** Handles infinite velocity/divergence perfectly. Zero race conditions. High memory coherency during rasterization.
*   **Cons:** Higher memory footprint (storing history buffer). Two-pass overhead (Write + Sort + Read).
*   **Note:** Emitting "Points" is preferred over "Line Segments" for performance, assuming the sample rate is sufficient to avoid aliasing.

### Pseudocode

**Phase 1: Emitter (Simulation)**
```cpp
void kernel_emit_samples(float4* sample_buffer, ...) {
    for(int step=0; step < n_steps; step++) {
        update_physics(pos, ...);
        
        // Emit sample every k steps
        if (step % stride == 0) {
            int idx = atomicAdd(&global_counter, 1);
            sample_buffer[idx] = make_float4(pos.x, pos.y, q, cell_index(pos));
        }
    }
}
```

**Phase 2: Sort (Host)**
```cpp
// Use CUB or Thrust to sort sample_buffer by cell_index
cub::DeviceRadixSort::SortPairs(..., sample_buffer_keys, sample_buffer_values, ...);
```

**Phase 3: Rasterizer (Tiled Gather)**
```cpp
// One Block = One 16x16 Tile of the Global Grid
__shared__ float sm_tile[16][16];

void kernel_rasterize_sorted(float4* sorted_samples, int* tile_offsets) {
    // 1. Identify range of samples belonging to this tile
    int start_idx = tile_offsets[blockIdx.x];
    int end_idx   = tile_offsets[blockIdx.x + 1];

    // 2. Cooperative Load & Accumulate
    // Threads loop over the sorted samples for this tile
    for (int i = start_idx + threadIdx.x; i < end_idx; i += blockDim.x) {
        float4 sample = sorted_samples[i];
        
        // Project sample to Shared Memory Tile
        float u = (sample.x - tile_origin.x) * inv_dx;
        float v = (sample.y - tile_origin.y) * inv_dy;
        
        atomicAdd_shared(&sm_tile[v][u], sample.z * weight); 
        // ... bilinear neighbors
    }
    
    __syncthreads();

    // 3. Write Tile to Global Memory (Coalesced)
    if (valid_pixel(threadIdx)) {
        global_grid[global_pixel_idx] += sm_tile[threadIdx.y][threadIdx.x];
    }
}
```

---

## 5. Summary of Caveats & Optimization

1.  **Shared Memory Bank Conflicts:** In Strategy A & B, `atomicAdd` to Shared Memory is fast *unless* all threads write to the same address. Ensure the patch size is larger than a Warp size or use padding to minimize conflicts.
2.  **Global Atomics:** Strategy A still uses global atomics at the end. This is safe because `N_Workgroups` (e.g., 100) is much smaller than `N_GridPoints` (65k), making collision probability zero.
3.  **Sort Overhead:** For Strategy B, sorting is the most expensive step. Ensure you pack data (e.g., `float2 pos` + `half2 q_idx`) to minimize bandwidth during the sort.
4.  **Sampling Theorem:** For Strategy B, the `sample_stride` must be small enough that a particle does not "skip" over important features of the potential, but large enough to avoid filling VRAM with correlated data.
5.  **Steady State Convergence:** Both methods mathematically converge to the same result: $\int \rho(t) dt$. Strategy A does it "in-flight"; Strategy B does it "post-mortem".

## 6. Final Recommendation

*   **Start with Strategy B (Trajectory Sorting).**
    *   It is more robust.
    *   It decouples physics bugs from grid bugs.
    *   It handles the "ergodic" case natively, which is likely for long-running steady-state convergence.
    *   The "overhead" of sorting is often less than the penalty of uncoalesced memory writes in a naive Scatter approach.
*   **Switch to Strategy A** only if you are memory-capacity bound (cannot store the trajectory buffer) or if the particles are known to be extremely stiff/slow-moving.