---
type: TopicalAudit
title: Linear Algebra & Iterative Solvers
tags: [topic, cross-language, linear-algebra, solvers]
---

## Summary

Dense and sparse linear algebra routines (matrix multiply, transpose, Gauss elimination, CG, BiCG, Jacobi eigenvalue) in C++ `Lingebra.h`. Sparse matrix support in `SparseMatrix.h`/`SparseMatrix2.h`. Python implementations in `pyTruss/IterativeLinearSolvers.py` (Jacobi, Gauss-Seidel, Chebyshev/momentum acceleration) and `pyTruss/truss_solver.py` (momentum-accelerated iterative solvers matching C++ architecture). Key application: truss dynamics linear system solve in `TrussDynamics_d.cpp` with 14 solver variants (CG, Cholesky, Jacobi, Gauss-Seidel, momentum variants, fly/diff variants).

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/common/math/Lingebra.h` / `.cpp` | active | Dense: `linSolve_gauss`, `linSolve_CG`, `linSolve_BCG`, `eig_Jacobi`, `leastSquareFit_Gauss`. Template `genLinSolve_CG` for custom dot funcs. |
| C++ | `cpp/common/math/SparseMatrix.h` | active | Sparse matrix CG solver |
| C++ | `cpp/common/math/SparseMatrix2.h` | active | Alternate sparse format |
| C++ | `cpp/common/dynamics/TrussDynamics_d.h/.cpp` | active | 14 solver methods via `LinSolveMethod` enum: CG, CGsparse, Cholesky, CholeskySparse, Jacobi, GaussSeidel, JacobiMomentum, GSMomentum, JacobiFlyMomentum, GSFlyMomentum, Force, JacobiDiff, MomentumDiff, ExternDiff |
| Python | `python/pyTruss/IterativeLinearSolvers.py` | active | `jacobi_step`, `gauss_seidel_step`, `chebyshev_accelerator`, `momentum_accelerator`, `solve_iterative` framework |
| Python | `python/pyTruss/truss_solver.py` | active | `solve_iterative_momentum` matching C++ `updateIterativeMomentum`, spectral radius estimation |
| Python | `python/constrain_solver.py` | experimental | Constraint solver (PGS-style?) |
| Python | `python/PGS.py` | experimental | Projected Gauss-Seidel |

## Parity Status

- **C++ `updateIterativeMomentum` ↔ Python `solve_iterative_momentum`**: Architecture matches (momentum mixing schedule, sub-method selection). Python uses NumPy arrays, C++ uses raw `Vec3d*`. Tolerance and test reference: `tests_bash/Orbital/` has truss regression data.
- **Chebyshev acceleration**: Python-only in `IterativeLinearSolvers.py`. C++ uses momentum (heavy-ball) instead. No direct parity.
- **CG/BiCG**: C++ `Lingebra.cpp` has both. Python relies on `numpy.linalg.solve` as reference in `solve_iterative`.

## Open Issues

- C++ `linSolve_CG`/`linSolve_BCG` have hardcoded `maxIters=10` and `maxErr2=1e-5` — not configurable
- C++ uses `delete` (not `delete[]`) for arrays allocated with `new[]` — potential UB
- Chebyshev acceleration exists only in Python, not ported to C++
- `SparseMatrix2.h` relationship to `SparseMatrix.h` unclear — possible duplication
- `constrain_solver.py` and `PGS.py` not audited — need to verify if they use same algorithmic framework
- See `doc/NumricalStabilityLinearSolvers.md` and `doc/ChebyshevAndInterialAccelerationOfLinearSolvers.md` for theoretical analysis
