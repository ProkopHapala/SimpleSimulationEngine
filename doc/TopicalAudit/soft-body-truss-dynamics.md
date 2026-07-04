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
| C++ | `cpp/common/dynamics/SoftBody.h` | active | Soft body dynamics base |
| C++ | `cpp/common/dynamics/TrussDynamics_d.h/.cpp` | active | Projective dynamics truss solver: CG, Cholesky, Jacobi, GaussSeidel, momentum variants, fly/diff variants. `LinSolveMethod` enum with 14 methods. |
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
