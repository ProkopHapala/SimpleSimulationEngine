---
type: TopicalAudit
title: ODE Integration
tags: [topic, cross-language, ode, runge-kutta, rkf45, euler, integrator]
---

## Summary

ODE integration for time-evolution of physical systems. C++ `ODEintegrator.h` provides base Euler and RKF45 (Runge-Kutta-Fehlberg 4(5)) with adaptive timestep. Used by orbital mechanics (`SpaceWorld`, `SpaceBodyIntegrator`), molecular dynamics, and other dynamic systems. Python has scattered ODE tests in `pySimE`. The RKF45 implementation supports both function-pointer and virtual-object derivative callbacks.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/common/dynamics/ODEintegrator.h` | active | `ODEintegrator` (Euler base) + `ODEintegrator_RKF45` (RKF45 with adaptive step). 6-stage RKF45 with error estimator. Adaptive dt: SAFETY=0.2, PGROW=1.2, PSHRINK=0.7, MAX_ADAPT=10. `dt_min=0.001`, `dt_max=0.1`. Supports `ODEderivFunc` (C function pointer) and `ODEderivObject` (virtual method). |
| C++ | `cpp/common/Orbital/SpaceBodies.h` | active | `SpaceBodyIntegrator : ODEderivObject, ODEintegrator_RKF45` — single-body gravitational integration with multiple centers. `evalTrj()` populates trajectory arrays. |
| C++ | `cpp/common/Orbital/SpaceWorld.h` | active | `SpaceWorld : ODEderivObject` — N-body gravitational + thrust. Uses `ODEintegrator_RKF45 ode` member. `predictTrjs()` integrates all bodies forward. |
| Python | `python/pySimE/` | partial | Various ODE test scripts (not audited in detail) |

## Sub-topics

### RKF45 (Runge-Kutta-Fehlberg 4(5))

6-stage method with embedded 4th/5th order error estimator:
- Stages: `a = {0, 1/4, 3/8, 12/13, 1, 1/2}`
- 5th order solution: `cs1 = {16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55}`
- 4th order solution: `cs2 = {25/216, 0, 1408/2565, 2197/4104, -1/5}`
- Error: `Yerr = dt·(dYdt2 - dYdt1)`

Adaptive step: if `max|Yerr·invMaxYerr| < 1` → accept, grow dt if error < SAFETY. Else → reject, shrink dt by PSHRINK.

### Interface Design

Dual callback mechanism:
- `ODEderivFunc` — C-style function pointer: `void (*)(double t, int n, double* Y, double* dY)`
- `ODEderivObject` — C++ virtual interface: `getDerivODE(double t, int n, double* Y, double* dY)`

Both checked in `step_RKF45()`: if `getDerivs` non-null use it, else if `derivObj` non-null use virtual method.

## Parity Status

- **C++ `ODEintegrator_RKF45` ↔ Python `pySimE`**: Python ODE tests not audited in detail. No formal parity test found.
- **`SpaceBodyIntegrator` ↔ `SpaceWorld`**: Both use RKF45 but different force models (single-body multi-center vs. full N-body). No consistency check.

## Open Issues

- No Verlet/leap-frog integrator in `ODEintegrator.h` (exists only in `DynamicOpt` as `move_LeapFrog`)
- No higher-order methods (RK8, symplectic integrator) for long-term orbital stability
- `invMaxYerr` array must be set by user — not initialized automatically
- `MAX_STEPS=10000` hardcoded — may be too low for some problems
- No energy conservation monitoring for Hamiltonian systems
- Python ODE tests in `pySimE` not audited — may be stale or incomplete
