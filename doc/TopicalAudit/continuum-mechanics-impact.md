---
type: TopicalAudit
title: Continuum Mechanics & Impact Fluid Dynamics
tags: [topic, python, opencl, eulerian, compressible, impact, nuclear, level-set, eos]
---

## Summary

Compressible Eulerian fluid dynamics for hypervelocity impact simulation (e.g. uranium projectile hitting liquid hydrogen — nuclear bomb implosion scenario). 5-equation model (Allaire et al.) with level-set interface tracking, stiffened gas EOS for high density ratios (100x). OpenCL GPU implementation with 8x8 tiled local memory. Rusanov flux. Redistancing for level-set maintenance. Python driver with PyOpenCL. Design discussion covers Level Set Method, Ghost Fluid Method, Mie-Grüneisen EOS.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| Python | `python/EulerianImpacFluid/EulerianImpacFluid.py` | active | PyOpenCL driver: buffer management, ping-pong swap, `step()`, `redistance()`, `get_data()`. NVIDIA GPU preference. |
| OpenCL | `python/EulerianImpacFluid/EulerianImpacFluid.cl` | active | 5-equation model kernel: `update_fluid` (Rusanov flux, level-set advection, positivity floors), `redistance_phi` (re-initialization). 8x8 tiles, 32 threads, local memory halo. Stiffened gas EOS mixture. |
| Python | `python/EulerianImpacFluid/test_eulerian_fluid.py` | active | Test driver: uranium projectile (ρ=19050) into liquid hydrogen (ρ=71), v=5000 m/s, p₀=1 GPa. Matplotlib animation. Mass/energy diagnostics. NaN detection. |
| Python | `python/EulerianImpacFluid/EulerianImpacFluid copy.py` | deprecated | Backup copy — exclude |
| Python | `python/EulerianImpacFluid/EulerianImpacFluid copy.cl` | deprecated | Backup copy — exclude |
| Doc | `python/EulerianImpacFluid/EulerianImpacFluid.md` | active | LLM chat: Level Set Method theory, 5-equation model derivation, Ghost Fluid Method, step-by-step update loop, OpenCL kernel design (2758 lines) |
| Doc | `doc/Continum_Mechanics_Solver.md` | active | Continuum mechanics solver design notes |

## Physics Model

### 5-Equation Model (Allaire et al.)

1. **Mass conservation (fluid 1)**: `∂(α₁ρ₁)/∂t + ∇·(α₁ρ₁u) = 0`
2. **Mass conservation (fluid 2)**: `∂(α₂ρ₂)/∂t + ∇·(α₂ρ₂u) = 0`
3. **Global momentum**: `∂(ρu)/∂t + ∇·(ρu⊗u + pI) = 0`
4. **Global total energy**: `∂E/∂t + ∇·((E+p)u) = 0`
5. **Interface advection (level set)**: `∂φ/∂t + u·∇φ = 0`

### Equation of State (Stiffened Gas Mixture)

- Material 1 (H₂): γ₁=1.4, P∞₁=0
- Material 2 (Uranium): γ₂=3.0, P∞₂=40 GPa
- Mixture pressure: `p = (e_int - Σ αᵢγᵢP∞ᵢ/(γᵢ-1)) / (Σ αᵢ/(γᵢ-1))`
- Volume fraction from level-set: `α₂ = H_ε(φ)`, `α₁ = 1 - α₂`

### Numerical Method

- **Flux**: Rusanov (Local Lax-Friedrichs) — `F = 0.5*(F_L + F_R) - 0.5*S_max*(U_R - U_L)`
- **Level-set advection**: Non-conservative upwind: `φ_new = φ - dt*(u*∂φ/∂x + v*∂φ/∂y)`
- **Redistancing**: Solve `∂φ/∂τ + S(φ₀)(|∇φ| - 1) = 0` for 5 iterations
- **Positivity floors**: Minimum density, minimum internal energy, vacuum momentum kill
- **Tiling**: 8x8 cells per workgroup, 32 threads, 10x10 local memory (with halo)

## Parity Status

- **OpenCL kernel ↔ Python test driver**: Driver initializes state, launches kernel, reads back diagnostics. CPU-side EOS diagnostic in `test_eulerian_fluid.py` verifies pressure initialization.
- **No CPU reference solver** — OpenCL is the only implementation. No parity test against analytical solution or external code.
- **Design doc (`EulerianImpacFluid.md`) ↔ implementation**: Doc describes full 5-equation model with MUSCL/WENO reconstruction; implementation uses first-order Rusanov (simpler, more dissipative but stable).

## Open Issues

- First-order reconstruction (no MUSCL/WENO slope limiter) — diffusive but stable
- No HLLC or Roe solver — Rusanov only (more dissipative)
- No Ghost Fluid Method — uses mixture EOS instead (simpler but less accurate for extreme density ratios)
- `EulerianImpacFluid copy.py` and `copy.cl` are backups — should be excluded/cleaned
- No automated regression test or analytical comparison
- No C++ counterpart — Python/OpenCL only
- Single-precision GPU — concerns for high-pressure regimes (see `doc/High_Precision_Physics_on_singlepoint_GPU_hardware.md`)
- Application: nuclear bomb implosion (uranium compression), impact fusion, hypervelocity impact
- `doc/Continum_Mechanics_Solver.md` — separate design doc, relationship to implementation unclear
