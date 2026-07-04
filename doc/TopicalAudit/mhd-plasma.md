---
type: TopicalAudit
title: MHD & Plasma Physics
tags: [topic, javascript, mhd, plasma, nozzle, coil, flux-conservation, elliptic-integrals]
---

## Summary

Axisymmetric MHD plasma nozzle simulation in JavaScript. Models plasma as a ring of current-carrying nodes in r-z plane, with superconducting (SC) coil, parabolic cage coils, and plasma ring. Flux conservation solver computes induced currents via mutual inductance matrix. Dynamic mode integrates plasma motion under Lorentz + gas pressure forces. Three.js rendering with GLSL shader for B-field visualization.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| JS | `js/mhd_demo/physics.js` | active | Core physics: `coilField()` (B-field from circular loop via elliptic K,E), `mutualInductance()` (Neumann formula), `selfInductance()`, `solveFluxConservation()` (Gaussian elimination for M·I = Φ₀), `initParabolicNozzle()`, `initTwoStates()`, `stepSimulation()`. Plasma as spring-connected nodes with gas pressure (adiabatic P·V^γ = const). |
| JS | `js/mhd_demo/render.js` | active | Three.js renderer: orthographic camera, coil ring visualization, B-field vector grid, GLSL fragment shader for HSV-colored B-field magnitude/direction background. Elliptic integrals in GLSL. |
| JS | `js/mhd_demo/ui.js` | active | UI: trajectory computation (interpolated + dynamic modes), energy plot, parameter sliders, coil editor. Dynamic mode: single plasma coil Euler integration with Lorentz + gas pressure forces. |
| JS | `js/mhd_demo/main.js` | active | Entry point: init state, renderer, UI, animation loop |
| JS | `js/mhd_demo/test_solver.js` | active | Offline Node.js test: init state, step simulation, check for NaN in currents |
| Doc | `docs/BiotSavart/MHD_Plasma_Nozzle_webgl.md` | doc | Design document for MHD plasma nozzle web visualization |
| Doc | `docs/BiotSavart/AeroPotentialFlow_simulation_javascript.md` | doc | Related: potential flow simulation in JavaScript |

## Sub-topics

### Coil Field Calculation

B-field from circular current loop at (r,z) using elliptic integrals:
- `k² = 4·a·r / ((a+r)² + z²)` (elliptic modulus squared)
- `ellipticK(m)` via AGM (arithmetic-geometric mean) iteration — 20 iterations
- `ellipticE(m)` via Carlson symmetric form
- `Br = common · (z/r) · (-K + (a²+r²+z²)/denom2 · E)`
- `Bz = common · (K + (a²-r²-z²)/denom2 · E)`
- Guard: `r < eps → Br = 0`, `k²` clamped to `[0, 1-1e-9]`

### Flux Conservation Solver

At t=0: `Φ₀ = M₀ · I₀` where only SC has current.
At time t: solve `M_t · I_t = Φ₀` for induced currents.
- Mutual inductance: `M_ij = μ₀·√(r_i·r_j)·[(2/k - k)·K - (2/k)·E]`
- Self inductance: `L ≈ μ₀·a·(ln(8a/r_w) - 2)`
- Gaussian elimination with partial pivoting (`gaussianSolve()`)
- RHS: `Φ₀_i - Σ_sc M_i_sc · I_sc` (SC contribution subtracted)

### Dynamic Simulation

Single plasma coil as point mass:
- Forces: Lorentz (`F = I·L·B` where `L = 2π·r`) + gas pressure (`P = P₀·(V₀/V)^γ`)
- Euler integration: `v += F/m·dt`, `r += v·dt`
- Guard: `r < r_min → r = r_min, v = -0.5·v` (axis collision)
- Parameters: `dt=5e-6s`, `mass=0.1kg`, `P₀=1e6 Pa`, `γ=5/3`

### GLSL B-field Shader

Fragment shader computes B-field on a quad mapped to world coordinates:
- Up to 64 coils with `[r, z, I, 0]` uniforms
- Elliptic integrals in GLSL (10-iteration AGM)
- HSV coloring: hue = field direction, value = field magnitude
- Alternative grayscale magnitude mode

## Parity Status

- **JS `solveFluxConservation()` ↔ Python `demo_plasma_dynamics_flux.py`**: UI code references Python implementation. Same algorithm (mutual inductance + Gaussian elimination). No formal parity test.
- **JS `coilField()` ↔ GLSL `coilField()`**: Same formula, different precision (JS double vs GLSL float). No formal parity test.
- **Dynamic mode ↔ Static solver mode**: Two trajectory modes (interpolated vs dynamic). No consistency check.

## Open Issues

- `ellipticK()` and `ellipticE()` in `physics.js` have unused `c` variable and `sum` not used in `ellipticK`
- GLSL shader uses `float uCoils[256]` (64×4) but WebGL1 uniform array size limits may be exceeded on some hardware
- Dynamic mode uses simplified single-coil plasma model — not multi-node
- `gaussianSolve()` skips singular columns (`|pivot| < 1e-12`) — may produce wrong results for ill-conditioned matrices
- No energy conservation check in dynamic mode
- `test_solver.js` only checks for NaN — no accuracy verification
- No Python counterpart for the full multi-node plasma simulation
- `stepSimulation()` not fully audited — spring forces and gas pressure model need verification
