# pyTruss

Truss/cloth simulation using Projective Dynamics with iterative linear solvers (Jacobi, Gauss-Seidel, VBD, Chebyshev, momentum) on CPU and GPU (OpenCL).

## Core modules

- **truss.py** — Truss data structure: points, bonds, masses, stiffness, fixed-point constraints. Single source of truth for geometry and physical parameters.
- **sparse.py** — Sparse linear algebra: Jacobi/Gauss-Seidel (dense + sparse), graph coloring, neighbor list building, `build_rest_length_dense`, general `linsolve_iterative` with momentum.
- **IterativeLinearSolvers.py** — Modular solver framework: `jacobi_step`, `gauss_seidel_step`, `chebyshev_omega_update`, `momentum_bmix_update`, `estimate_spectral_radius`, composable `solve_iterative` wrapper.
- **projective_dynamics.py** — PD system assembly: `make_pd_matrix`, `make_pd_rhs`, `solve_pd` (direct), `makeSparseSystem`. The mathematical foundation for all iterative solvers.
- **projective_dynamics_iterative.py** — PD with iterative solvers: `solve_pd_jacobi`, `solve_pd_chebyshev`. Uses shared components from `sparse.py` and `IterativeLinearSolvers.py`.
- **truss_solver.py** — CPU solver with pluggable backends: VBD, Jacobi-diff/fly, GS-diff/fly, momentum, Chebyshev. `TrussSolver` class, `SOLVERS` dict, `run_solver_suite` for multi-solver comparison.
- **plot_utils.py** — Matplotlib visualization: `plot_truss` (LineCollection-based), `plot_graph_coloring`.

## GPU solvers

- **truss_ocl.py** — Original monolithic OpenCL solver (Jacobi, GS, VBD, diff). Legacy, kept for reference.
- **truss_solver_ocl.py** — Old GPU wrapper (VBD-serial only). Uses standalone `OCLRuntime`.
- **truss_solver_ocl_new.py** — Refactored GPU solver (VBD, Jacobi-diff, Jacobi-fly). Inherits `OpenCLBase` from pyMolecular. Pre-allocated buffers, kernel source in `truss.cl`.

## Run scripts

- **run_vbd_cloth.py** — CLI runner: CPU + old GPU backend. Grid truss, selectable solver, before/after plot.
- **run_vbd_cloth_new.py** — CLI runner: CPU + new GPU backend. Adds `--gpu`, `--nloc`, `--device`, `--sub-method` options.
- **run_vbd_cloth_old.py** — Legacy runner using monolithic `truss_ocl.py` directly.
- **run_solver_debug.py** — CPU-only solver debugging: controlled perturbation, per-iteration logging, `--solver-suite` for multi-solver convergence comparison.
- **example_ocl_new.py** — Minimal usage example for the new GPU solver API.

## Test scripts

- **test_graph_coloring.py** — Graph coloring visualization on grid trusses.
- **test_Chebyshev_accel.py** — Convergence comparison: plain vs Chebyshev vs momentum on random SPD matrices.
- **test_Jacobi_Chebyshev_convergence.py** — Convergence comparison: Jacobi vs Chebyshev/inertial on a PD system.

## Other files

- **truss.cl** — OpenCL kernel source for GPU solvers.
- **ARCHITECTURE.md** — Old vs new GPU architecture comparison.
- **IMPLEMENTATION_SUMMARY.md** — Implementation notes for the new GPU backend.
- **backup_png/** — Pre-refactoring output images.
- **refactored_png/** — Post-refactoring output images for review.
