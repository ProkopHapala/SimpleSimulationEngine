---
type: TopicalAudit
title: Soft Body & Truss Dynamics
tags: [topic, cross-language, physics, soft-body, truss, projective-dynamics]
---

## Summary

Soft body dynamics via truss/mass-spring systems. C++ implementation in `SoftBody.h` and `TrussDynamics_d.h/.cpp` with projective dynamics formulation and 14 linear solver variants. Python reimplementation in `pyTruss/` with NumPy. Truss generation/mesh building in `docs/TrussGeneration/`. Vertex Block Descent algorithm documented separately. Application: `BlockHouseTactics` game (bridge builder style).

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/common/dynamics/SoftBody.h` | active | Soft body dynamics base: bonds (springs) and kinks (angular springs). `SoftBodyLinearized` for small-deformation linear elasticity. |
| C++ | `cpp/common/dynamics/TrussDynamics_d.h/.cpp` | active | Projective dynamics truss solver: CG, Cholesky, Jacobi, GaussSeidel, momentum variants, fly/diff variants. `LinSolveMethod` enum with 14 methods. |
| C++ | `cpp/common/dynamics/TrussDynamics_f.h/.cpp` | active | Single-precision (float) mirror of TrussDynamics_d for GPU upload. Adds `impuls[]` buffer for correction impulses. |
| C++ | `cpp/common/dynamics/SoftPolyLine2D.h` | active | 2D deformable polyline: radial (spring) and angular (bending) forces. Complex-number rotation reference. |
| C++ | `cpp/common/dynamics/SoftPolyLine3D.h` | active | 3D deformable polyline: radial, angular (cross-product bending), torsion forces. Per-segment forward/up vectors. |
| C++ | `cpp/common/dynamics/Chain2D.h` | active | 2D chain with fixed-length constraints via iterative relaxation from pivot. |
| C++ | `cpp/common/dynamics/LinearElasticity.h` | active | Assembly of 3D global stiffness matrix from bond connectivity (6x6 blocks per bond). For static solve or modal analysis. |
| C++ | `cpp/apps/BlockHouseTactics/` | active | Truss simulation game (BridgeBuilder-like) |
| Python | `python/pyTruss/` | active | Full Python reimplementation: `truss_solver.py` (1249 lines), `IterativeLinearSolvers.py` |
| Python | `python/pyTruss/ARCHITECTURE.md` | active | Architecture documentation |
| Python | `python/pyTruss/IMPLEMENTATION_SUMMARY.md` | active | Implementation report |
| Python | `python/pyTruss/README_new_solver.md` | active | New solver documentation |
| Python | `python/pyTruss/Vivace_graphColoring_GaussSeidel.md` | active | Graph coloring for parallel GS |
| Doc | `doc/VertexBlockDescent.md` | active | Vertex Block Descent algorithm (19KB) |
| Doc | `doc/Markdown/cpp/common/dynamics/TrussDynamics_d.h.md` | active | Auto-generated API doc |
| Doc | `doc/Markdown/cpp/common/dynamics/ProjectiveDynamics.md` | active | Projective dynamics notes |
| Doc | `doc/Markdown/cpp/common/OCL/ProjectiveDynamicsOCL.md` | active | OpenCL projective dynamics |
| Doc | `docs/TrussGeneration/` | active | Truss mesh generation: `ConstructionBlock.md`, `MeshBuilder2.md`, `truss_low_to_hi.md`, `MeshFaceGenerationDiscrepancy.md` |
| Test | `tests_bash/Orbital/` | active | Truss regression tests with LDLT reference data |

## Parity Status

- **C++ `TrussDynamics_d::updateIterativeMomentum` ↔ Python `solve_iterative_momentum`**: Architecture matches. Both use momentum mixing schedule (`bmix`), sub-method selection (jacobi_diff, gs_diff, jacobi_fly, gs_fly). C++ uses `Vec3d*` arrays, Python uses NumPy.
- **C++ `updateJacobi_lin` ↔ Python `_update_jacobi_lin`**: Direct port. Same formula: `ps_out[i] = (bi.f + sum_j) / bi.w`.
- **C++ `updateGaussSeidel_lin` ↔ Python `_update_gs_lin`**: Direct port. In-place update.
- **C++ `updatePD_dRHS` ↔ Python `_build_linear_diff_system`**: Both compute RHS for displacement formulation.
- **Test reference**: `tests_bash/Orbital/` has LDLT decomposition reference data (`LDLT_L_dense.txt`, `LDLT_LT_reconstructed.txt`).

## Open Issues

- C++ `linSolveMethod` default is `2` (Cholesky) — iterative methods need more testing
- `MeshFaceGenerationDiscrepancy.md` documents a known mesh generation bug
- Graph coloring for parallel Gauss-Seidel: Python-only (`Vivace_graphColoring_GaussSeidel.md`), not in C++
- `ProjectiveDynamicsOCL.md` suggests OpenCL port exists or planned — status unclear
- Linear elasticity notes: `cpp/common/dynamics/Notes/LinarElasticity.md` (low quality)
- See `doc/NumricalStabilityLinearSolvers.md` for stability analysis of these solvers

---

## Algorithm Review: Soft Body & Truss Dynamics

This section reviews algorithms for soft body / truss / mass-spring dynamics — implicit time integration, linear and nonlinear solvers, acceleration schemes, and mesh generation. For each algorithm we describe: merit, use cases, strengths, weaknesses, GPU parallelization potential, and implementation status across platforms (C++, Python, OpenCL, Julia).

### A. Time Integration Formulations

#### A1. Implicit Euler (Variational Form / Projective Dynamics)

- **Merit**: Unconditionally stable time integration. Each step minimizes global energy $G(x) = \frac{1}{2h^2}\|x-y\|_M^2 + E(x)$ where $y = x^t + hv^t + h^2 a^{ext}$. Velocity recovered as $v^{t+1} = (x^{t+1}-x^t)/h$. Allows large timesteps without instability.
- **Use cases**: All truss/soft-body simulations in the project. Bridge builder game. Cloth simulation. Molecular mechanics.
- **Strengths**: Unconditionally stable. Energy-decreasing (for reasonable dt). Same framework supports springs, collisions, friction. Local-global splitting enables iterative solvers.
- **Weaknesses**: First-order accuracy — energy dissipation (numerical damping). Nonlinear when constraints are re-projected ("fly" variant). Requires solving a linear (or nonlinear) system each step.
- **GPU parallelization**: The system solve is the bottleneck — parallelization depends on the solver choice (see B1–B6). The predictor and corrector steps are embarrassingly parallel.
- **Implemented**:
  - C++: `TrussDynamics_d::run_LinSolve()` in `cpp/common/dynamics/TrussDynamics_d.cpp`. Predictor → `linsolve()` → corrector loop. Documented in `doc/NumricalStabilityLinearSolvers.md:30-57`.
  - Python: `TrussSolver.run_solver()` in `python/pyTruss/truss_solver.py`. Same structure with pluggable solver callbacks.
  - OpenCL: `OCL_Orb::run_projective_dynamics()` — predictor/corrector kernels in `cpp/common_resources/cl/spacecraft.cl`. Documented in `doc/Markdown/cpp/common/OCL/ProjectiveDynamicsOCL.md`.
  - Julia: `Julia/try_julia_plot/ProjectiveDynamics.jl` — mirrors same structure.
  - See also: `doc/TopicalAudit/projective-dynamics.md` (dedicated audit).

#### A2. Explicit Euler / Leap-Frog (Not Used for Truss)

- **Merit**: Simplest integrator. No system solve needed.
- **Use cases**: Not used for truss dynamics (too stiff). Used in molecular dynamics and orbital mechanics.
- **Strengths**: O(N) per step. No linear solve.
- **Weaknesses**: Conditionally stable — dt limited by stiffest spring (dt < 2√(m/k)). Useless for stiff trusses.
- **GPU parallelization**: Excellent — fully parallel force evaluation.
- **Implemented**: Not used in truss. `TrussDynamics_d` has `move_MD()` and `move_LeapFrog()` but these are for relaxation/optimization, not production dynamics. See `doc/TopicalAudit/ode-integration.md` for ODE integrator inventory.

#### A3. Position-Based Dynamics (PBD)

- **Merit**: Directly manipulates positions to satisfy constraints (not forces). Iteratively project each constraint. Very stable, allows large dt. Popular for cloth and deformable bodies in games.
- **Use cases**: Cloth, real-time interactive simulation. Molecular cluster relaxation.
- **Strengths**: Unconditionally stable. No stiffness matrix. Intuitive constraint formulation. Fast convergence for visual quality.
- **Weaknesses**: Stiffness depends on iteration count (not physical). Not energy-conserving. Linear convergence only. No natural support for momentum/inertia.
- **GPU parallelization**: Good — constraint projections are parallel within a color group. Same coloring approach as VBD.
- **Implemented**:
  - C++: `cpp/sketches_SDL/Molecular/test_PBD_LJ_cluster.cpp` — PBD for Lennard-Jones cluster relaxation. Not integrated into `TrussDynamics_d`.
  - Python/OpenCL: Not implemented for truss.
  - Gap: PBD is not used in the main truss pipeline. Could be an alternative for real-time game scenarios.

#### A4. XPBD (Extended Position-Based Dynamics)

- **Merit**: Extension of PBD that decouples stiffness from iteration count by introducing compliance parameters (inverse stiffness C = 1/k). Physically accurate — matches implicit Euler in the limit. Supports constitutive models (elasticity, plasticity). From Macklin & Müller (2013).
- **Use cases**: Production cloth and soft body simulation. Game physics. FEM-like deformation with PBD stability.
- **Strengths**: Stiffness is physically meaningful (not iteration-dependent). Unconditionally stable. Supports arbitrary constitutive models via strain energy. GPU-friendly (same parallel structure as PBD).
- **Weaknesses**: More complex than PBD (requires compliance matrix per constraint). Still first-order. Convergence depends on iteration count for accuracy (though not for stability).
- **GPU parallelization**: Excellent — same as PBD/VBD, constraint projections are parallel within color groups.
- **Implemented**: **Not implemented**. No XPBD code exists in the repo. This is a gap — XPBD would be a natural complement to the existing PD framework, offering simpler implementation and better game-realism tradeoffs.

---

### B. Linear & Nonlinear Solvers for the PD System

#### B1. Linearized Jacobi (Precomputed Stiffness)

- **Merit**: Solve $A p = b$ where $A = M/h^2 + K$ (precomputed). Jacobi update: $p_i^{new} = (b_i + \sum_j k_{ij} p_j) / A_{ii}$. Simplest parallel solver.
- **Use cases**: GPU default solver. Fast approximate solve.
- **Strengths**: Embarrassingly parallel. No synchronization needed. O(N) per iteration. Simple implementation.
- **Weaknesses**: Slow convergence (linear rate). Spectral radius ρ close to 1 for stiff systems → many iterations needed. No error correction.
- **GPU parallelization**: Excellent — each vertex updated independently. Already implemented in OpenCL.
- **Implemented**:
  - C++: `TrussDynamics_d::updateJacobi_lin()` in `TrussDynamics_d.cpp`. Uses precomputed `bvec` and `kngs`.
  - Python: `solve_jacobi_diff()` in `truss_solver.py:807`. Displacement-based variant.
  - OpenCL: `updateJacobi_lin` kernel in `spacecraft.cl`. GPU `ialg=3`.
  - Julia: Available in `ProjectiveDynamics.jl`.

#### B2. Jacobi "Fly" (Nonlinear Constraint Projection)

- **Merit**: Re-evaluates constraint projection each iteration: $p'_{ij} = p_j + d_{ij} \cdot (l_0/l)$. Goes beyond linear $Ap=b$ — updates RHS every iteration with nonlinear projection. Handles asymmetric stiffness (kPress ≠ kPull).
- **Use cases**: When constraints are nonlinear (large deformation). When precomputed stiffness is stale.
- **Strengths**: More accurate than linearized Jacobi for large deformations. No matrix assembly needed. Handles tension/compression asymmetry.
- **Weaknesses**: More computation per iteration (norm, division per neighbor). Still linear convergence. No guarantee of energy decrease.
- **GPU parallelization**: Excellent — same parallel structure as B1, slightly more arithmetic per neighbor.
- **Implemented**:
  - C++: `TrussDynamics_d::updateJacobi_fly()` in `TrussDynamics_d.cpp`.
  - Python: `solve_jacobi_fly()` in `truss_solver.py:871`.
  - OpenCL: `updateJacobi_neighs` kernel in `spacecraft.cl`. GPU `ialg=0`.

#### B3. Gauss-Seidel (In-Place Update)

- **Merit**: Uses updated values immediately: $p_i^{new} = (b_i + \sum_{j<i} k_{ij} p_j^{new} + \sum_{j>i} k_{ij} p_j^{old}) / A_{ii}$. Faster convergence than Jacobi (roughly 2×).
- **Use cases**: CPU high-accuracy solve. When fewer iterations are needed.
- **Strengths**: 2× faster convergence than Jacobi. Simple implementation on CPU.
- **Weaknesses**: Sequential dependency — not trivially parallel. Requires graph coloring for GPU. More memory traffic (in-place updates).
- **GPU parallelization**: Requires graph coloring — vertices in same color are independent. Python has coloring; OpenCL does not implement GS. VBD uses coloring for its parallel coordinate descent (see B5).
- **Implemented**:
  - C++: `TrussDynamics_d::updateGaussSeidel_lin()` and `updateGaussSeidel_fly()`.
  - Python: `solve_gs_diff()` in `truss_solver.py:914` and `solve_gs_fly()` at `:973`.
  - OpenCL: **Not implemented** — gap. Graph coloring exists in Python (`Vivace_graphColoring_GaussSeidel.md`) but no OpenCL GS kernel.
  - Doc: `doc/Parallel_GaussSeidel.md` (303 lines) — detailed OpenCL kernel design for colored GS with local memory, analytic 3×3 solves, and synchronization. Ready for implementation.

#### B4. Momentum Mixing (SmartMixer / Heavy-Ball)

- **Merit**: Accelerates Jacobi/GS by mixing current iterate with previous delta: $p^{new} = p^{jacobi} \cdot (1-b) + p^{prev\_delta} \cdot b$. Heavy-ball / Polyak momentum method. `SmartMixer` ramps `bmix` from 0 → 0.75 over iterations 3→10.
- **Use cases**: Default acceleration method in C++ and GPU. Reduces iteration count by 2–5×.
- **Strengths**: Simple — one extra vector blend. No inner products (GPU-friendly). Adaptive schedule. Composable with any base solver.
- **Weaknesses**: `bmix` schedule is hand-tuned (not adaptive to problem). `SmartMixer::get_bmix()` simplified to step function (linear ramp commented out). May overshoot for some problems.
- **GPU parallelization**: Excellent — one vector blend per iteration, fully parallel.
- **Implemented**:
  - C++: `TrussDynamics_d::updateIterativeMomentum()` with `SmartMixer` struct in `TrussDynamics_d.h:76-92`. `bmix` schedule: pure Jacobi for itr<3, then b_end=0.75.
  - Python: `solve_iterative_momentum()` in `truss_solver.py`. Matches C++ schedule.
  - OpenCL: `updateJacobi_mix` kernel (ialg=1) and `updateJacobi_smart` (ialg=2, adaptive). Documented in `ProjectiveDynamicsOCL.md:40-42`.

#### B5. Chebyshev Semi-Iterative Acceleration

- **Merit**: 3-term recurrence: $x^{(n)} = \omega_n (\tilde{x}^{(n)} - x^{(n-2)}) + x^{(n-2)}$ with $\omega_{n} = 4/(4 - \rho^2 \omega_{n-1})$. Optimal polynomial for reducing error over $[-\rho, \rho]$. Quadratic convergence improvement over plain iteration.
- **Use cases**: Accelerating any stationary iterative solver (Jacobi, GS, VBD). Particularly effective for PD where spectral radius is stable.
- **Strengths**: No inner products needed (unlike CG) — fully GPU-parallel. Only needs two previous iterates. Theoretical optimality for linear systems. Well-documented in `doc/ChebyshevAndInterialAccelerationOfLinearSolvers.md` (427 lines).
- **Weaknesses**: Requires spectral radius estimate ρ (via power iteration or manual tuning). If ρ too high → divergence. If too low → no acceleration. Nonlinear systems (fly Jacobi) violate linear assumptions but work in practice (heuristic).
- **GPU parallelization**: Excellent — no global reductions, just vector operations. Ideal for GPU.
- **Implemented**:
  - Python: `IterativeLinearSolvers.py` — `chebyshev_accelerator()`, `chebyshev_omega_update()`, `solve_jacobi_chebyshev()`, `solve_gauss_seidel_chebyshev()`. `test_Chebyshev_accel.py` for convergence testing. `projective_dynamics_iterative.py` — `solve_pd_chebyshev()`.
  - C++: **Not implemented** as standalone. `SmartMixer` (B4) serves a similar role but is simpler. The VBD document describes Chebyshev for VBD acceleration (Sec. 3.8) but it's not in C++ `TrussDynamics_d`.
  - OpenCL: **Not implemented** — gap. Documented as part of VBD algorithm but no kernel exists.
  - VBD: Chebyshev acceleration is part of the VBD algorithm specification (`doc/VertexBlockDescent.md:228-238`) — skip for colliding vertices. Not yet implemented in Python VBD or OpenCL VBD kernels.

#### B6. Cholesky / LDLT Direct Factorization

- **Merit**: Factor $A = L D L^T$ once, then solve by forward/back substitution each step. Exact (to machine precision). O(N³) factorization, O(N²) solve.
- **Use cases**: Small systems (< few thousand vertices). Reference solution for testing iterative solvers. When matrix doesn't change between steps (constant stiffness).
- **Strengths**: Exact. No iteration count to tune. Best accuracy. Reusable factorization across timesteps if K is constant.
- **Weaknesses**: O(N³) factorization — doesn't scale. Dense matrix storage. No GPU implementation. Fill-in for sparse systems.
- **GPU parallelization**: Poor — dense factorization is hard to parallelize efficiently on GPU. cuSOLVER could be used but not integrated.
- **Implemented**:
  - C++: `TrussDynamics_d::prepare_LinearSystem()` builds `PDmat`, `LDLT_L`, `LDLT_D`. `LinSolveMethod::Cholesky` (default, `linSolveMethod=2`). Sparse variant: `CholeskySparse` using `PDsparse`.
  - Python: `solve_cholesky()` in `truss_solver.py:560`. Uses NumPy `np.linalg.solve(L, rhs)`.
  - OpenCL: **Not available** — CPU-only. Experimental `fwd_subs` kernel exists but not integrated.
  - Test: `tests_bash/Orbital/` has LDLT reference data (`LDLT_L_dense.txt`, `LDLT_LT_reconstructed.txt`).

#### B7. Conjugate Gradient (CG)

- **Merit**: Krylov subspace method for SPD systems. Optimal for symmetric positive definite A. Converges in at most N iterations (exact arithmetic). Superlinear convergence typical.
- **Use cases**: Medium-large systems where Cholesky is too expensive but Jacobi is too slow. Good preconditioner candidate.
- **Strengths**: Optimal for SPD. Better convergence than Jacobi/GS. Memory: only 4 vectors (x, r, p, Ap).
- **Weaknesses**: Requires inner products (global reductions on GPU — synchronization barrier). Matrix must be SPD. Sensitive to conditioning (needs preconditioner).
- **GPU parallelization**: Moderate — sparse mat-vec is parallel, but dot products require global reduction (slow on GPU). Preconditioned CG (PCG) with Jacobi preconditioner is common.
- **Implemented**:
  - C++: `CGsolver` in `TrussDynamics_d.h:159`. `LinSolveMethod::CG` (0) and `CGsparse` (1). `cpp/common/math/CG.h`. Hardcoded `maxIters=10`, `maxErr2=1e-5` — not configurable.
  - Python: Uses `np.linalg.solve` as direct reference. No iterative CG in truss solver (relies on `IterativeLinearSolvers.py` framework).
  - OpenCL: **Experimental only** — `dot_mat_vec_loc` kernel exists but not integrated into PD driver.

#### B8. Vertex Block Descent (VBD)

- **Merit**: Coordinate descent over vertices — minimize local energy $G_i(x)$ per vertex using a single Newton step: $\Delta x_i = -H_i^{-1} f_i$ where $H_i$ is 3×3 Hessian (analytic inverse). Gauss-Seidel-like but with per-vertex Newton instead of linear update. From A. Chen (2024).
- **Use cases**: Cloth, soft body, contact-rich simulation. GPU-friendly alternative to global solvers. Handles nonlinear constraints naturally (spring projection, collisions).
- **Strengths**: No global system assembly — purely local 3×3 solves. Naturally parallel via graph coloring. Handles nonlinear energies (springs, collisions, friction) in unified framework. Chebyshev acceleration applicable. Adaptive initialization (scaled acceleration prediction). Collision handling via penalty energy + friction.
- **Weaknesses**: Sequential over color groups (not fully parallel). Convergence depends on coloring quality. Per-vertex Newton may overshoot for highly coupled systems. No global convergence guarantee for nonlinear case. Requires Hessian computation per vertex per iteration (more arithmetic than Jacobi).
- **GPU parallelization**: Good — vertices within a color are independent. One workgroup per vertex, threads sum force contributions. Chunking needed to fit neighbor lists in local memory. `build_vbd_chunks()` in `truss_ocl.py` handles this.
- **Implemented**:
  - Python (CPU): `solve_vbd()` in `truss_solver.py:438-547`. Serial loop over all vertices. Computes 3×3 Hessian + gradient per vertex, solves via `np.linalg.solve(H, grad)`. No coloring (serial). No Chebyshev. No adaptive init. No collisions.
  - OpenCL: `vbd_vertex_serial` kernel in `python/pyTruss/truss.cl:351` — single-vertex serial kernel. `vbd_vertex_chunk` kernel at `:270` — workgroup-parallel kernel with local memory for chunked vertices. `solve3x3()` inline function at `:67`. `build_vbd_chunks()` in `truss_ocl.py:32` splits colors into workgroup-sized batches.
  - C++: **Not implemented** — gap. VBD is documented (`doc/VertexBlockDescent.md`, 505 lines) but no C++ code exists. The document includes OpenCL kernel code and design notes for modular integration.
  - Doc: `doc/VertexBlockDescent.md` — full algorithm (global optimization, local solver, collisions, friction, initialization, Chebyshev acceleration, parallelization, Algorithm 1 pseudocode). Two versions: original rewrite + Google Gemini version.

#### B9. Augmented VBD (AVBD)

- **Merit**: Extension of VBD that adds augmented Lagrangian terms to improve convergence for stiff constraints and inextensibility. Multiplier updates enforce constraints more strongly than pure penalty methods.
- **Use cases**: Stiff cloth (near-inextensible). Volume preservation. Rigid body coupling.
- **Strengths**: Better constraint enforcement than penalty-only VBD. Maintains VBD's parallel structure.
- **Weaknesses**: More complex (multiplier fields, update schedule). Limited literature.
- **GPU parallelization**: Good — same structure as VBD with additional multiplier updates per vertex.
- **Implemented**: **Not implemented**. No AVBD code or documentation exists. Gap — would be relevant if VBD proves insufficient for stiff constraints in the bridge builder game.

#### B10. Multigrid Solver

- **Merit**: Solve on hierarchy of coarser grids — smooth on fine grid, restrict residual to coarse grid, solve, prolongate correction back. Optimal O(N) complexity for elliptic problems.
- **Use cases**: Large soft body systems (10K+ vertices). When O(N) scaling is needed. Cloth with high resolution meshes.
- **Strengths**: Optimal O(N) convergence. Geometry-independent convergence rate (for elliptic problems). Well-established theory.
- **Weaknesses**: Requires hierarchy of meshes (coarsening). Complex implementation (restriction/prolongation operators, V-cycles/W-cycles). Hard to parallelize on GPU (coarse grids are small). Not natural for unstructured truss meshes (need algebraic multigrid).
- **GPU parallelization**: Moderate — fine grid smoothing is parallel, but coarse levels have fewer vertices. AMG (algebraic multigrid) can help but is complex.
- **Implemented**: **Not implemented**. No multigrid code exists in the repo. Mentioned only in passing in unrelated docs. Gap — would be relevant for large-scale truss simulations but requires significant implementation effort.

---

### C. Acceleration & Optimization Methods

#### C1. FIRE (Fast Inertial Relaxation Engine)

- **Merit**: Adaptive timestep damped dynamics optimizer. If force-velocity alignment is positive, grow dt and reduce damping (accelerate). If negative, reset velocity, shrink dt, increase damping (brake). From Bitzek et al. (2014).
- **Use cases**: Static relaxation of truss structures. Molecular mechanics geometry optimization. Trajectory optimization.
- **Strengths**: Superlinear convergence for smooth landscapes. No line search needed. Adaptive dt. Simple implementation.
- **Weaknesses**: Not a time integrator (no energy conservation). Heuristic parameters (finc, fdec, falpha). Can oscillate near saddle points. Not suitable for dynamic simulation (only relaxation).
- **GPU parallelization**: Force evaluation is parallel. The FIRE control logic (dt adjustment, velocity reset) is scalar — done on CPU.
- **Implemented**:
  - C++: `DynamicOpt::move_FIRE()` in `cpp/common/dynamics/DynamicOpt.h`. `TrussDynamics_d::FIRE_update()` and `run_omp()` with OpenMP parallelization. Uses cosine similarity threshold (cs<0.02) instead of vf<0.
  - Python: Not in pyTruss. Used in `python/pyMolecular/RigidMol.py` via ctypes binding to C++ `DynamicOpt`.
  - OpenCL: Not implemented.
  - Doc: `doc/TopicalAudit/optimization.md` — dedicated audit of FIRE and other optimizers.

#### C2. Gradient Descent & Damped MD (Relaxation)

- **Merit**: Simple force-based position update: $x \leftarrow x + dt \cdot v$, $v \leftarrow (1-damp) \cdot v + dt \cdot f/m$. Velocity reset when $f \cdot v < 0$.
- **Use cases**: Quick relaxation. Fallback when more sophisticated optimizers are unavailable.
- **Strengths**: Trivial implementation. Stable with enough damping.
- **Weaknesses**: Very slow convergence. No adaptive timestep. Sensitive to damping parameter.
- **GPU parallelization**: Excellent — fully parallel.
- **Implemented**:
  - C++: `TrussDynamics_d::move_GD()` and `move_MD()` in `TrussDynamics_d.h:447-449`.
  - Python: `gradient_descent_relax()` and `dynamical_relaxation()` in `python/subPixelContour/subPixelContour.py`.

---

### D. Differential Formulation (Displacement-Based Solve)

#### D1. Displacement-Based PD (Diff Solvers)

- **Merit**: Solve for displacements $dp = p - p_0$ instead of absolute positions. RHS: $b = K \cdot d - A \cdot p'$. Better conditioned for single-precision (avoids large coordinate cancellation).
- **Use cases**: GPU single-precision solvers. Large world coordinate offsets. When float32 precision is insufficient for absolute position solve.
- **Strengths**: Better numerical stability in single precision. Smaller dynamic range. Compatible with Jacobi, GS, and momentum variants.
- **Weaknesses**: Requires storing neutral positions $p_0$. Slightly more complex RHS assembly. Still O(N) per iteration.
- **GPU parallelization**: Same as base solver (Jacobi/GS/momentum) — no additional parallelization challenges.
- **Implemented**:
  - C++: `TrussDynamics_d::updatePD_dRHS()` assembles displacement RHS. `updateIterativeJacobiDiff()`, `updateIterativeMomentumDiff()`, `updateIterativeExternDiff()`. `LinSolveMethod::JacobiDiff` (11), `MomentumDiff` (12), `ExternDiff` (13).
  - Python: `solve_jacobi_diff()` in `truss_solver.py:807`. `solve_gs_diff()` at `:914`. `_build_linear_diff_system()` helper.
  - OpenCL: **Not implemented** — major gap. `ProjectiveDynamicsOCL.md:54-58` documents the TODO: needs `updatePD_dRHS` kernel and displacement buffers (`ibuff_ps0`, `ibuff_dpos`).
  - Doc: `doc/NumricalStabilityLinearSolvers.md` (770 lines) — full analysis of displacement-based iterative solvers.

---

### E. Collision & Contact Handling

#### E1. Edge-Vertex Collision (TrussDynamics)

- **Merit**: Detect and resolve collisions between truss vertices and edges. Edge-vertex bond constraint with interpolation parameter `c` and stiffness `K`. Damping via relative velocity projection.
- **Use cases**: Self-collision in truss structures. Vertex-edge contact in bridge builder game.
- **Strengths**: Integrated into force evaluation. Handles damping. Interpolated contact point on edge.
- **Weaknesses**: Penalty-based (not hard constraint). Collision damping is a free parameter. No CCD (continuous collision detection) — discrete only.
- **GPU parallelization**: Moderate — collision detection is parallel, but pair resolution may have conflicts.
- **Implemented**:
  - C++: `TrussDynamics_d::evalEdgeVert()` in `TrussDynamics_d.h:318-357`. `EdgeVertBond` struct at `:69-74`. Bounding boxes via `Buckets` for broad-phase.
  - Python/OpenCL: Not implemented for truss dynamics.
  - Doc: `doc/TopicalAudit/collision-detection.md` — broad-phase collision audit (hash grid, sweep-and-prune, AABB tree).

#### E2. VBD Collision & Friction Model

- **Merit**: Quadratic penalty energy $E_c = \frac{1}{2} k_c d^2$ for contacts, with friction force from tangential relative motion. Gradient and Hessian added to per-vertex local system. Smooth static-to-dynamic friction transition via $f_1(u)$.
- **Use cases**: Cloth self-collision. Contact-rich soft body simulation. VBD framework.
- **Strengths**: Unified with VBD energy framework — collisions are just another energy term. Friction handled in same Newton step. Hessian contributes to local 3×3 solve.
- **Weaknesses**: Penalty-based (stiffness k_c is a parameter). Normal assumed constant during differentiation. Requires CCD/DCD for detection (not part of solver).
- **GPU parallelization**: Good — collision energy terms are per-vertex, added to local Hessian/gradient. No global synchronization.
- **Implemented**:
  - Doc only: `doc/VertexBlockDescent.md:153-196` — full equations for collision energy, gradient, Hessian, and friction. Algorithm 1 includes DCD/CCD steps.
  - Python VBD: `solve_vbd()` in `truss_solver.py` does **not** implement collisions — only spring energy.
  - OpenCL VBD: `vbd_vertex_serial` and `vbd_vertex_chunk` kernels do **not** include collision terms.
  - C++: `TrussDynamics_d` has its own collision system (E1), not VBD-style.
  - Gap: VBD collision model is documented but not implemented in any code.

---

### F. Parallelization & Graph Coloring

#### F1. Graph Coloring for Parallel Gauss-Seidel

- **Merit**: Color vertices so no two adjacent vertices share a color. Process each color in parallel (vertices within a color are independent). Enables parallel GS.
- **Use cases**: Parallel GS on GPU. VBD parallelization.
- **Strengths**: Enables GS-level convergence with Jacobi-level parallelism. Well-studied (graph theory).
- **Weaknesses**: Number of colors = sequential steps (typically 3–10 for mesh graphs). Coloring is a preprocessing step (NP-hard in general, but greedy works for meshes). Load balancing across colors.
- **GPU parallelization**: Essential — this IS the parallelization strategy. One kernel launch per color.
- **Implemented**:
  - Python: `python/pyTruss/Vivace_graphColoring_GaussSeidel.md` (119 lines) — graph coloring design. `test_graph_coloring.py` — coloring tests. `truss_ocl.py:build_vbd_chunks()` — chunk packing for VBD.
  - C++: **Not implemented** — gap. GS is sequential on CPU.
  - OpenCL: **Not implemented for GS** — but VBD kernels use coloring (`vbd_vertex_chunk` processes per-color chunks).
  - Doc: `doc/Parallel_GaussSeidel.md` (303 lines) — detailed OpenCL kernel design with local memory, coloring, and synchronization. Ready for implementation.

#### F2. VBD Chunk Packing

- **Merit**: Split each color group into workgroup-sized chunks that fit in local memory. Each chunk includes vertex indices and their complete neighbor lists. Respects workgroup size (32) and neighbor limit (128).
- **Use cases**: GPU VBD with local memory optimization.
- **Strengths**: Eliminates global memory traffic during inner loops. Compact neighbor lists. Workgroup-sized batches.
- **Weaknesses**: Preprocessing overhead (chunk building). Irregular chunk sizes. May waste threads if chunks are small.
- **GPU parallelization**: This IS the GPU strategy — chunks map to workgroups.
- **Implemented**:
  - Python: `build_vbd_chunks()` in `python/pyTruss/truss_ocl.py:32-50`. `VBD_WORKGROUP_SIZE=32`, `VBD_NEIGHBOR_LIMIT=128`.
  - OpenCL: `vbd_vertex_chunk` kernel in `truss.cl:270` uses chunk buffers.

---

### G. Mesh Generation & Truss Construction

#### G1. ConstructionBlock Mesh Builder

- **Merit**: Generate truss meshes from construction blocks — structured grids, bridges, towers. Low-to-high resolution refinement.
- **Use cases**: Bridge builder game (`BlockHouseTactics`). Test meshes for solver comparison.
- **Strengths**: Structured, parameterized. Supports refinement.
- **Weaknesses**: Known mesh face generation discrepancy (bug documented).
- **GPU parallelization**: N/A — preprocessing.
- **Implemented**:
  - C++: `Mesh::Builder2` (referenced in `docs/TrussGeneration/MeshBuilder2.md`).
  - Doc: `docs/TrussGeneration/ConstructionBlock.md`, `MeshBuilder2.md`, `truss_low_to_hi.md`, `MeshFaceGenerationDiscrepancy.md`, `DrawUV.md`.
  - Python: `Truss` class in `python/pyTruss/truss.py` with builders `make_grid`, `make_rope`.

---

### H. Implementation Gap Summary

| Algorithm | C++ | Python | OpenCL | Julia | Gap |
|-----------|-----|--------|--------|-------|-----|
| Implicit Euler (PD) | ✅ | ✅ | ✅ | ✅ | — |
| Explicit Euler / Leap-Frog | ✅ (relax) | ❌ | ❌ | ❌ | Not used for truss |
| PBD | ✅ (molecular only) | ❌ | ❌ | ❌ | Not in truss pipeline |
| XPBD | ❌ | ❌ | ❌ | ❌ | Not implemented |
| Jacobi (linearized) | ✅ | ✅ | ✅ | ✅ | — |
| Jacobi "fly" (nonlinear) | ✅ | ✅ | ✅ | ❌ | — |
| Gauss-Seidel | ✅ | ✅ | ❌ | ✅ | OpenCL gap (design doc exists) |
| Momentum (SmartMixer) | ✅ | ✅ | ✅ | ❌ | — |
| Chebyshev acceleration | ❌ | ✅ | ❌ | ❌ | C++/OCL gap |
| Cholesky / LDLT | ✅ | ✅ | ❌ | ❌ | OCL gap (experimental only) |
| Conjugate Gradient | ✅ | ❌ (ref only) | ❌ (exp.) | ❌ | Not integrated on GPU |
| Vertex Block Descent (VBD) | ❌ | ✅ (CPU) | ✅ | ❌ | C++ gap |
| Augmented VBD (AVBD) | ❌ | ❌ | ❌ | ❌ | Not implemented |
| Multigrid | ❌ | ❌ | ❌ | ❌ | Not implemented |
| FIRE optimizer | ✅ | ❌ (pyTruss) | ❌ | ❌ | Not in Python truss |
| Gradient Descent / Damped MD | ✅ | ✅ | ❌ | ❌ | — |
| Diff (displacement) solvers | ✅ | ✅ | ❌ | ❌ | OCL gap (documented TODO) |
| Edge-vertex collision | ✅ | ❌ | ❌ | ❌ | Python/OCL gap |
| VBD collision + friction | ❌ (doc) | ❌ | ❌ | ❌ | Doc only, not implemented |
| Graph coloring (parallel GS) | ❌ | ✅ | ❌ (VBD only) | ❌ | C++/OCL gap for GS |
| VBD chunk packing | ❌ | ✅ | ✅ | ❌ | — |
| Mesh generation | ✅ | ✅ | N/A | ❌ | — |

### I. Related Documentation

| Document | Path | Content |
|----------|------|---------|
| VBD algorithm | `doc/VertexBlockDescent.md` | Full VBD specification (505 lines): global optimization, local Newton solver, collisions, friction, adaptive init, Chebyshev acceleration, parallelization, Algorithm 1 pseudocode, OpenCL 3×3 solve code |
| PD audit | `doc/TopicalAudit/projective-dynamics.md` | Dedicated PD topical audit: system matrix, fly Jacobi, SmartMixer, diff formulation, CPU↔GPU kernel table |
| PD OpenCL mapping | `doc/Markdown/cpp/common/OCL/ProjectiveDynamicsOCL.md` | CPU↔GPU kernel correspondence, absolute vs differential, momentum details, TODO for GPU diff solvers |
| Chebyshev design | `doc/ChebyshevAndInterialAccelerationOfLinearSolvers.md` | 427-line design doc: theory, update formula, spectral radius estimation, nonlinear adaptation, parameter tuning |
| Numerical stability | `doc/NumricalStabilityLinearSolvers.md` | 770-line analysis: Jacobi/GS on displacements, C++ code, stability of single-precision |
| Parallel GS design | `doc/Parallel_GaussSeidel.md` | 303-line OpenCL kernel design: local memory, coloring, 3×3 analytic solves, synchronization |
| Graph coloring | `python/pyTruss/Vivace_graphColoring_GaussSeidel.md` | Graph coloring for parallel GS on GPU |
| Linear algebra audit | `doc/TopicalAudit/linear-algebra-solvers.md` | Cross-language solver inventory: CG, BiCG, Gauss elimination, Jacobi eigenvalue |
| Optimization audit | `doc/TopicalAudit/optimization.md` | FIRE, gradient descent, random optimization, line search |
| Collision audit | `doc/TopicalAudit/collision-detection.md` | Broad-phase: hash grid, sweep-and-prune, AABB tree |
| TrussDynamics API | `doc/Markdown/cpp/common/dynamics/TrussDynamics_d.h.md` | Auto-generated C++ API documentation |
| Truss generation | `docs/TrussGeneration/` | `ConstructionBlock.md`, `MeshBuilder2.md`, `truss_low_to_hi.md`, `MeshFaceGenerationDiscrepancy.md`, `DrawUV.md` |
| pyTruss architecture | `python/pyTruss/ARCHITECTURE.md` | Python solver architecture |
| pyTruss implementation | `python/pyTruss/IMPLEMENTATION_SUMMARY.md` | Implementation report |
| pyTruss new solver | `python/pyTruss/README_new_solver.md` | New solver documentation |

### J. Key Discrepancies & Open Opportunities

1. **VBD not in C++**: VBD is fully documented and implemented in Python + OpenCL, but not in C++ `TrussDynamics_d`. The design notes in `doc/VertexBlockDescent.md:17-54` describe a modular solver interface that would allow VBD to plug in alongside existing PD solvers.
2. **Chebyshev not in C++/OpenCL**: Chebyshev acceleration is well-documented and implemented in Python `IterativeLinearSolvers.py`, but the C++ and GPU paths use only `SmartMixer` (heavy-ball momentum). Chebyshev could offer better convergence for the GPU Jacobi/VBD path.
3. **Diff solvers not on GPU**: The displacement-based formulation (better for float32) is implemented on C++ and Python but not OpenCL. This is the most impactful GPU gap — single-precision drift on large coordinates.
4. **Gauss-Seidel not on GPU**: Graph coloring exists in Python, and a detailed OpenCL kernel design exists (`doc/Parallel_GaussSeidel.md`), but no implementation. GS would roughly halve iteration count vs Jacobi.
5. **VBD collisions not implemented**: The VBD collision/friction model is fully specified in the doc but not implemented in any code. The Python `solve_vbd()` only handles spring energy.
6. **XPBD not considered**: XPBD would be a natural complement — simpler than PD for game scenarios, physically accurate stiffness, GPU-friendly. No code or discussion exists.
7. **Multigrid not considered**: For large-scale truss systems (10K+ vertices), multigrid would offer O(N) scaling. No code or discussion exists. Would require algebraic multigrid for unstructured truss meshes.
8. **AVBD not considered**: Augmented VBD could improve convergence for stiff constraints (near-inextensible cloth, volume preservation). No code or discussion exists.
9. **CG not integrated on GPU**: Experimental CG kernels (`dot_mat_vec_loc`) exist but are not wired into the PD driver. CG with Jacobi preconditioner could outperform momentum-accelerated Jacobi for medium-sized problems.
10. **SmartMixer simplified**: `SmartMixer::get_bmix()` was simplified from a linear ramp to a step function (commented out code at `TrussDynamics_d.h:82-86`). This may affect convergence quality.
