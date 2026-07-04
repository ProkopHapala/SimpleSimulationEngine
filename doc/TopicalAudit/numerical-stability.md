---
type: TopicalAudit
title: Numerical Stability & Precision
tags: [topic, cross-language, numerical-stability, precision, gpu, float, double]
---

## Summary

Numerical stability concerns across the repo: iterative linear solver convergence (Jacobi/Gauss-Seidel/CG), catastrophic cancellation in single-precision GPU physics, and compensated arithmetic (Knuth's TwoSum) for split-float high-precision on GPU. The main documentation covers stability of iterative solvers for truss dynamics (implicit Euler) and a design document for high-precision physics on single-precision GPU hardware.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `doc/NumricalStabilityLinearSolvers.md` | doc | 770-line treatise on iterative solver stability for truss dynamics. Covers implicit Euler system matrix `A = M/dt² + K`, Jacobi/Gauss-Seidel iteration matrices, spectral radius convergence, Chebyshev acceleration, momentum mixing. Includes C++ code for `run_LinSolve()` and `updateIterativeMomentum()`. |
| C++ | `doc/High_Precision_Physics_on_singlepoint_GPU_hardware.md` | doc | Design document for split-float (compensated arithmetic via Knuth's TwoSum) and split-int (fixed-point hybrid) on single-precision GPUs. Covers catastrophic cancellation in large-world coordinates, `P = P_coarse + P_fine` decomposition, integration into MD loop, comparison table. |
| C++ | `doc/ChebyshevAndInterialAccelerationOfLinearSolvers.md` | doc | Chebyshev acceleration and Heavy Ball momentum for iterative linear solvers |
| C++ | `doc/Parallel_GaussSeidel.md` | doc | Parallel Gauss-Seidel variants for GPU/multi-core |
| Python | `python/IterativeLinearSolvers.py` | active | `get_iteration_matrix()` — computes spectral radius of Jacobi/GS iteration matrices. `jacobi_step()`, `gauss_seidel_step()`, `chebyshev_accelerator()`, `momentum_accelerator()`, `solve_iterative()` — general framework with convergence tracking (residual + error norms). See also linear-algebra-solvers.md. |
| C++ | `cpp/common/dynamics/DynamicOpt.h` | active | `ff_safety = 1e-32` — prevents division by zero in FIRE velocity mixing. `scale_dt` limits prevent force/velocity blow-up. |
| C++ | `python/EulerianImpacFluid/EulerianImpacFluid.cl` | active | Positivity floors in Rusanov flux: `p = max(p, p_floor)`, `rho = max(rho, rho_floor)` — prevents negative density/pressure in compressible flow. |

## Sub-topics

### Iterative Solver Convergence

- System matrix: `A = M/dt² + K` (inertia + stiffness)
- Jacobi iteration matrix: `T_J = -D⁻¹(L+U)`, converges if `ρ(T_J) < 1`
- Gauss-Seidel iteration matrix: `T_GS = -(D+L)⁻¹U`
- Chebyshev acceleration: `x_{k+1} = x_{k-1} + ω(x̃_{k+1} - x_{k-1})` with optimal `ω` schedule
- Heavy Ball momentum: `x_{k+1} = x̃_{k+1} + β(x_k - x_{k-1})`
- Convergence criterion: residual norm `< tol` or force norm `< ftol`

### GPU Single-Precision Problem

- `float` has ~7 decimal digits → catastrophic cancellation for large coordinates
- Example: `p_i = 10000.123`, `p_j = 10000.000` → `dij = 0.12` (lost 0.003)
- Split-float: `P = P_coarse + P_fine` using Knuth's TwoSum for error tracking
- Split-int: fixed-point hybrid with integer subtraction (zero error for `dij`)
- Recommendation: Split-float for most physics (simpler, robust TwoSum integration)

### Positivity Preservation (Eulerian Fluid)

- Pressure floor: `p = max(p, p_floor)` in Rusanov flux
- Density floor: `rho = max(rho, rho_floor)`
- These prevent negative values that would destabilize the compressible Euler solver

## Parity Status

- **Python `IterativeLinearSolvers.py` ↔ C++ truss dynamics solvers**: Python implements Jacobi/GS steps with NumPy; C++ implements same algorithms in `TrussDynamics_d`. No formal parity test.
- **Split-float design ↔ implementation**: Design document only — no C++ or OpenCL implementation found in repo.

## Open Issues

- Split-float / split-int not implemented — only design document exists
- No automated convergence rate comparison between solvers
- `ff_safety = 1e-32` in `DynamicOpt` may be too small for some problems
- No energy drift monitoring for long-time integration (orbital mechanics, MD)
- Python `IterativeLinearSolvers.py` uses `np.linalg.solve(A,b)` as reference — O(n³) not scalable for benchmarking
- No condition number estimation for system matrices
