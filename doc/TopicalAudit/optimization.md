---
type: TopicalAudit
title: Optimization Methods
tags: [topic, cross-language, optimization, fire, gradient-descent, line-search, random]
---

## Summary

Optimization methods across the repo: FIRE (Fast Inertial Relaxation Engine) as the primary workhorse for molecular mechanics, truss dynamics, and trajectory optimization; gradient descent and damped dynamics for curve fitting; random/stochastic optimizers for global search; line search (Regula Falsi with Illinois correction) for root finding. C++ `DynamicOpt` is the central FIRE implementation, reused by multiple subsystems via `bindArrays()`. Python implementations mirror the same algorithms for NumPy-based workflows.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/common/dynamics/DynamicOpt.h` / `.cpp` | active | FIRE optimizer: `move_FIRE()`, `move_GD()`, `move_MD()`, `move_LeapFrog()`. Adaptive dt (finc/fdec/falpha), damping schedule, force/velocity limits. `bindArrays()` for zero-copy. `optimize()` convergence loop. |
| C++ | `cpp/common/dynamics/TrussDynamics_d.cpp` | active | `FIRE_update()` — inline FIRE variant for OpenMP parallel truss relaxation. Uses cosine similarity threshold (cs<0.02) instead of vf<0. `run_omp()` with `#pragma omp` parallelization. |
| C++ | `cpp/common/math/lineSearch.h` | active | `LineSearch` class: Regula Falsi with Illinois algorithm correction. Bounded/unbounded modes. Linear extrapolation for ray-marching style extension. Used for root finding (ray-shape intersection). |
| C++ | `cpp/common/math/optimizer_random.h` | active | `OptimizerRandom` — random perturbation + accept/reject. `OptimizerRandom_2` — with momentum term and decay. Stochastic global optimization. |
| C++ | `cpp/common/Orbital/TrajectoryVariation.h` | active | `getTrajectoryVariation()` — objective function for trajectory optimization (minimize thrust integral along B-spline) |
| C++ | `cpp/apps/OrbitalWar/test_OptContinuousThrust.cpp` | active | FIRE + B-spline trajectory optimization demo (20 control points) |
| C++ | `cpp/apps/MolecularEditor/doc/GlobalOptimization.md` | doc | Global optimization strategy for molecular conformations |
| Python | `python/pyMolecular/RigidMol.py` | active | `setOptFIRE()` — ctypes binding to C++ `DynamicOpt::setOptFIRE()`. Parameters: dt_max, dt_min, damp_max, minLastNeg, finc, fdec, falpha, kickStart. |
| Python | `python/pySimE/space/exp/OrbitalTransferOpt/Random_optimization.py` | active | `MC_Run()` — Monte Carlo random search with adaptive step size. `MCBias_Run()` — MC with bias momentum direction. Used for orbital transfer optimization. |
| Python | `python/subPixelContour/subPixelContour.py` | active | `solve_with_gradient_descent()` — GD with clamped error for implicit contour fitting. `solve_with_damped_dynamics()` — damped dynamics with inertia, velocity reset on f·v<0. `gradient_descent_relax()` / `dynamical_relaxation()` — model-agnostic relaxators. Rational model optimization (coefficients + non-uniform weights). |
| Python | `python/IterativeLinearSolvers.py` | active | Chebyshev acceleration and Heavy Ball momentum for iterative linear solvers (see linear-algebra-solvers.md) |

## Sub-topics

### FIRE (Fast Inertial Relaxation Engine)

Core algorithm in `DynamicOpt::move_FIRE()`:
1. Compute `ff = Σ f²`, `vv = Σ v²`, `vf = Σ v·f`
2. If `vf < 0`: reset velocity, increase damping, shrink dt
3. Else: `v = (1-damp)·v + damp·sqrt(vv/ff)·f`; if enough steps since last negative, grow dt and shrink damping
4. Leap-frog position update with scaled dt

Parameters: `finc=1.1`, `fdec=0.5`, `falpha=0.98`, `damp_max=0.2`, `minLastNeg=5`

Variant in `TrussDynamics_d::FIRE_update()`: uses cosine similarity `cs = vf/sqrt(vv·ff)` with threshold 0.02 instead of simple `vf<0` check. Also adds velocity magnitude guard (`vv>1600` → reset).

### Gradient Descent & Damped Dynamics (Python)

- `solve_with_gradient_descent()`: standard GD with clamped error, optional early stop on gradient norm
- `solve_with_damped_dynamics()`: inertia-based with `v = μ·v + lr·f`, velocity reset when `f·v < 0` (same principle as FIRE but simpler)
- `gradient_descent_relax()` / `dynamical_relaxation()`: model-agnostic versions accepting `force_fn(u) -> (f, E)` callback
- Rational model: `F = (A@(w·c))/(A@w)` with separate gradients for coefficients `c` and weights `w`

### Random/Stochastic Optimization

- `OptimizerRandom`: uniform random perturbation within `dxMax` bounds, accept if fitness improves
- `OptimizerRandom_2`: adds momentum direction with decay on failure
- `MC_Run()` / `MCBias_Run()`: Python MC with adaptive step size (grow on hit, shrink on miss ratio), optional bias direction

## Parity Status

- **C++ `DynamicOpt::move_FIRE()` ↔ Python `setOptFIRE()`**: Python calls C++ via ctypes. Same parameters. Direct parity.
- **C++ `DynamicOpt::move_FIRE()` ↔ `TrussDynamics_d::FIRE_update()`**: Variant implementations. TrussDynamics uses cosine similarity threshold and OpenMP parallelization. No formal parity test.
- **C++ `DynamicOpt` ↔ Python `dynamical_relaxation()`**: Same physics (damped dynamics with velocity reset on f·v<0). Python is model-agnostic, C++ is array-based. No formal parity test.
- **C++ `OptimizerRandom` ↔ Python `MC_Run()`**: Same concept (random perturbation + accept/reject). Python adds adaptive step size. No formal parity test.

## Open Issues

- `DynamicOpt::move_FIRE()` has an older commented-out version (lines 84-133) with force limit handling — current active version (line 136+) removed force limit check
- `TrussDynamics_d::FIRE_update()` uses hardcoded `vv>1600` velocity guard — should be parameterized
- `LineSearch.insertPoint()` not implemented (placeholder `return -1`)
- No L-BFGS or conjugate gradient implementation in C++ (CG exists only for linear systems in `Lingebra.h`)
- No automated convergence comparison between FIRE and GD
- `optimizer_random.h` uses C `randf()` — not thread-safe, quality limited
- Python `MC_Run()` uses Python 2 print statements — needs porting to Python 3
