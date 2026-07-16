---
type: TopicalAudit
title: Projective Dynamics
tags: [topic, cpp, opencl, projective-dynamics, jacobi, truss, implicit-euler, constraint-solving]
---

## Summary

Projective Dynamics (PD) for truss/spring-mass systems: implicit Euler time stepping with iterative linear solvers for the stiffness system `A·p = b` where `A = M/dt² + K`. CPU implementation in `TrussDynamics_d` with multiple solver variants (Jacobi, Gauss-Seidel, Cholesky, CG, momentum-accelerated). GPU/OpenCL implementation in `OCL_Orb` with Jacobi family kernels. Local-global optimization: each Jacobi step re-projects constraint positions (`p'_{ij}`) and solves weighted average.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/common/dynamics/TrussDynamics_d.h` | active | Main class: `points`, `vel`, `forces`, `neighs`, `params` (l0,kPress,kPull,damping). `SmartMixer` for adaptive momentum. `LinSolveMethod` enum: Jacobi, JacobiFly, JacobiMomentum, GS, Cholesky, CG, Diff variants. `kFix` for fixed points. |
| C++ | `cpp/common/dynamics/TrussDynamics_f.h` | active | Float-precision mirror of TrussDynamics_d for GPU upload. Same solver structure with float arithmetic. Adds `impuls[]` buffer. |
| C++ | `cpp/common/dynamics/TrussDynamics_d.cpp` | active | `run_LinSolve()` — main implicit Euler loop (predictor → solve → corrector). `updatePD_RHS()` — assembles `b_i` and `A_ii`. `updatePD_dRHS()` — displacement-based RHS for single-precision stability. `updateJacobi_lin()` — linearized Jacobi using precomputed `bvec`/`kngs`. `updateJacobi_fly()` — on-the-fly constraint projection (nonlinear). `updateGaussSeidel_lin()` — in-place GS. `dotPD()` — matrix-vector product `f = A·p`. `evalTrussForce()` — nonlinear spring force. `FIRE_update()` — FIRE optimizer variant. |
| C++ | `cpp/common/OCL/OCL_Orb.h` | active | GPU PD: `run_projective_dynamics(nSolverIters, ialg)`. Predictor/corrector kernels. Jacobi variants: `run_updateJacobi_neighs()` (ialg=0, fly), `run_updateJacobi_mix()` (ialg=1, momentum), `run_updateJacobi_smart()` (ialg=2, adaptive), `run_updateJacobi_lin()` (ialg=3, linearized). `run_SolverConvergence()` — per-iteration residual download for benchmarking. |
| C++ | `cpp/common/OCL/OCL_Orb_cpp.h` | active | C++ OpenCL wrapper version (cl::Buffer API). Same kernel structure as `OCL_Orb.h` (C API). |
| OpenCL | `cpp/common_resources/cl/spacecraft.cl` | active | Kernels: `PD_perdictor` (leap-frog velocity + position predict), `PD_corrector` (velocity from position change), `updateJacobi_neighs` (fly Jacobi with constraint projection), `updateJacobi_mix` (momentum blending), `updateJacobi_lin` (linearized using `bvec`/`kngs`), `updatePD_RHS` (assemble diagonal + RHS), `getTrussForces`, `evalTrussForce2`. |
| Doc | `doc/Markdown/cpp/common/OCL/ProjectiveDynamicsOCL.md` | doc | Detailed CPU↔GPU mapping with kernel correspondence table, absolute vs differential formulations, momentum mixing details, TODO for GPU diff solvers |
| Doc | `doc/NumricalStabilityLinearSolvers.md` | doc | Iterative solver convergence analysis (see numerical-stability.md) |

## Sub-topics

### PD System Matrix

`A = M/dt² + K` where:
- `M` = diagonal mass matrix
- `K` = stiffness matrix from spring/bond constraints
- `A_ii = M_i/dt² + Σ_j k_ij` (diagonal)
- `b_i = M_i/dt² · p'_i + Σ_j k_ij · d_ij` (RHS with inertial prediction + constraint offsets)

Jacobi update: `p_i^{new} = (b_i + Σ_j k_ij · p_j) / A_ii`

### "Fly" Jacobi (Nonlinear)

`updateJacobi_fly()` / `updateJacobi_neighs` kernel: re-evaluates constraint projection each iteration:
- `d_ij = p_i - p_j`, `l = |d_ij|`
- `p'_{ij} = p_j + d_ij · (l0/l)` (ideal position satisfying constraint)
- `k = (l < l0) ? kPress : kPull` (asymmetric stiffness)
- `p_i^{new} = (M_i/dt² · p_i + Σ_j k_ij · p'_{ij}) / (M_i/dt² + Σ_j k_ij)`

This goes beyond linear `Ap=b` — updates `b` every iteration with nonlinear projection.

### Momentum Mixing (SmartMixer)

`updateJacobi_mix()` / `updateIterativeMomentum()`:
- `p_i^{new} = p_i^{jacobi} · bmix.x + dps_i · bmix.y`
- `dps_i = p_i^{new} - p_i` (store delta for next iteration)
- `SmartMixer`: `bmix.y` ramps from 0.55→0.75 over iterations 3→10, pure Jacobi before iteration 3
- GPU `run_updateJacobi_smart()` mirrors this schedule

### Differential Formulation (CPU-only)

`updatePD_dRHS()` assembles `db = b - A·p₀` for displacement-based solve:
- `bi = Σ_j k_ij · d_ij · (l0/l - 1)` (displacement RHS)
- Better conditioned for single-precision (avoids large coordinate cancellation)
- **Not ported to GPU** — major gap

### CPU↔GPU Kernel Correspondence

| CPU LinSolveMethod | GPU ialg | GPU Kernel | Notes |
|---|---|---|---|
| Jacobi (linear) | 3 | `updatePD_RHS` + `updateJacobi_lin` | Uses precomputed `bvec`/`kngs` |
| JacobiFly | 0 | `updateJacobi_neighs` | On-the-fly constraint projection |
| JacobiMomentum | 1 | `updateJacobi_mix` | Momentum blending via `dps` |
| JacobiMomentum+Smart | 2 | `updateJacobi_smart` | Adaptive `bmix` schedule |
| JacobiDiff/MomentumDiff | — | Not implemented | GPU gap: needs `updatePD_dRHS` |
| GaussSeidel | — | Not available | CPU-only (sequential dependency) |
| Cholesky | — | Not available | CPU-only (LDLT factorization) |
| CG | — | Experimental only | `dot_mat_vec_loc` kernel exists but not integrated |

## Parity Status

- **CPU `updateJacobi_lin()` ↔ GPU `updateJacobi_lin` kernel**: Same formula `(bi + Σ_j k·p_j) / Aii`. Documented in `ProjectiveDynamicsOCL.md`. No automated parity test but `run_SolverConvergence()` enables comparison.
- **CPU `updateJacobi_fly()` ↔ GPU `updateJacobi_neighs`**: Same nonlinear constraint projection. No automated test.
- **CPU `updateIterativeMomentum()` ↔ GPU `updateJacobi_mix`**: Same momentum formula. `SmartMixer` constants match. No automated test.
- **CPU `run_LinSolve()` ↔ GPU `run_projective_dynamics()`**: Same predictor-solve-corrector structure. No end-to-end parity test.
- **Differential solvers**: CPU-only, no GPU counterpart — documented gap.

## Open Issues

- GPU lacks differential (`Diff`) solver variants — single-precision drift on large coordinates
- GPU lacks Gauss-Seidel, Cholesky, CG — only Jacobi family available
- `PD_corrector` kernel has hardcoded 2D debug: `{ pi.z = 0.0f; v_new.z = 0.0f; }` — should be conditional
- `SmartMixer::get_bmix()` simplified to step function (commented out linear ramp) — may affect convergence
- `run_PDcl()` hardcodes `ialg=3` (linear Jacobi) — other variants not tested in PDcl path
- No automated CPU↔GPU convergence comparison test
- `OCL_Orb.h` (C API) and `OCL_Orb_cpp.h` (C++ API) duplicate logic — maintenance burden
- Experimental CG kernels (`dot_mat_vec_*`, `fwd_subs`) not integrated into PD driver
