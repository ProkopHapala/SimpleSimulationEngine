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

---

## Algorithm Review: MHD & Plasma Simulation

This section reviews algorithms relevant to axisymmetric MHD plasma simulations — magnetic field computation, flux conservation, current solving, plasma dynamics, and visualization. For each algorithm we describe: merit, use cases, strengths, weaknesses, GPU parallelization potential, and implementation status across platforms (Python, JavaScript, GLSL).

### A. Magnetic Field Computation (Biot–Savart for Circular Loops)

#### A1. Elliptic Integral B-field from Circular Current Loop

- **Merit**: Exact analytic expression for the magnetic field of a circular current loop at arbitrary (r,z) using complete elliptic integrals K(m) and E(m). The fundamental building block of all axisymmetric B-field calculations.
- **Use cases**: Computing B-field from SC coils, cage coils, plasma current rings. Field visualization. Force computation. Flux computation.
- **Strengths**: Exact (to machine precision given accurate K,E). O(1) per evaluation point per coil. Well-known standard formula. Handles off-axis points (not just on-axis).
- **Weaknesses**: Elliptic integrals are expensive (iterative AGM or special functions). Singularity at r=0 (Br component) and at coil location (r=a, z=z₀). Requires regularization (ε guards, clamping k² to [0, 1-ε]).
- **GPU parallelization**: Good — each (pixel, coil) pair is independent. The AGM iteration for K,E is a fixed-loop kernel (10–20 iterations). Already implemented in GLSL fragment shader for visualization. For many coils × many evaluation points, this is embarrassingly parallel. Cost = O(N_coils × N_pixels × N_AGM_iters).
- **Implemented**:
  - JS: `coilField(r, z, I, a)` in `js/mhd_demo/physics.js:54-76`. Uses AGM (20 iterations) for K, Carlson symmetric form for E.
  - GLSL: `coilField()` in `js/mhd_demo/render.js:148-163` (fragment shader). AGM with 10 iterations. Float precision.
  - Python: `field_loop_rz(a, z0, I, r_grid, z_grid)` in `doc/python/MHD/inductance_core.py:283-326`. Uses `scipy.special.ellipk/ellipe` (machine precision). Vectorized over grid.
  - Python (on-axis): `coil_axis_Bz(a, z0, I, z_eval)` in `inductance_core.py:442-464`. Simple closed-form for r=0.
  - C++/OpenCL: Not implemented.

#### A2. Elliptic Integrals via AGM (Arithmetic-Geometric Mean)

- **Merit**: Iterative method to compute K(m) and E(m). Each iteration roughly doubles the number of correct digits. 10 iterations → ~1000 digits theoretically; 20 iterations → float64 precision.
- **Use cases**: Computing K,E when `scipy.special` is not available (JS, GLSL, embedded).
- **Strengths**: Simple loop, no lookup tables. Converges quadratically. Works in any language.
- **Weaknesses**: Slower than hardware-optimized `scipy.special` on CPU. Float precision in GLSL limits accuracy (10 iterations is marginal for float32). The `sum` variable in `ellipticK` is computed but unused (known bug in `physics.js`).
- **GPU parallelization**: Excellent — fixed-iteration loop, no branching, thread-independent. Already running in GLSL.
- **Implemented**:
  - JS: `ellipticK(m)` and `ellipticE(m)` in `js/mhd_demo/physics.js:11-52`. 20 iterations.
  - GLSL: `ellipticK(m)` and `ellipticE(m)` in `js/mhd_demo/render.js:123-146`. 10 iterations.
  - Python: `scipy.special.ellipk/ellipe` (not AGM). Also `fast_eliptke_Abramowitz.py` has polynomial approximation for speed.
  - C++/OpenCL: Not implemented.

#### A3. Point Dipole Field (Axial)

- **Merit**: Closed-form B-field for an axial magnetic dipole. Simple rational function, no elliptic integrals.
- **Use cases**: Approximating distant coils as dipoles. Diamagnetic bubble simulations. Dipole-source boundary conditions.
- **Strengths**: O(1), no special functions. Simple formula.
- **Weaknesses**: Only valid for far-field (r >> a). Loses near-field structure.
- **GPU parallelization**: Excellent — trivially parallel, simple arithmetic.
- **Implemented**:
  - Python: `field_dipole_rz(m, z_m, r_grid, z_grid)` in `inductance_core.py:329-366`. Vectorized.
  - JS/GLSL/C++: Not implemented.

---

### B. Inductance Matrix & Flux Conservation

#### B1. Mutual Inductance (Neumann Formula)

- **Merit**: Exact mutual inductance between two coaxial circular loops via elliptic integrals. M_ij = μ₀√(r_i·r_j)·[(2/k - k)·K - (2/k)·E].
- **Use cases**: Building the inductance matrix K for flux conservation solver. Computing coupling between coils.
- **Strengths**: Exact analytic formula. O(1) per pair. Symmetric: M_ij = M_ji.
- **Weaknesses**: Same elliptic integral cost as A1. Singularity when loops coincide (i=j requires self-inductance instead). k→0 (coincident, same z) and k→1 (r→0) need guards.
- **GPU parallelization**: Good — each (i,j) pair is independent. O(N²) parallelism for N coils. For typical N < 100, not worth GPU overhead. But for large N (many plasma nodes), could benefit.
- **Implemented**:
  - JS: `mutualInductance(r1, z1, r2, z2)` in `js/mhd_demo/physics.js:695-718`. Uses AGM-based K,E.
  - Python: `calc_mutual_inductance(r1, z1, r2, z2)` in `inductance_core.py:26-52`. Uses `scipy.special`.
  - C++/OpenCL: Not implemented.

#### B2. Self-Inductance (Thin-Wire Approximation)

- **Merit**: L ≈ μ₀·a·(ln(8a/r_wire) - c) where c depends on the approximation (2 for simple, 1.75 for Rosa-Kirchhoff). Regularizes the i=j diagonal of the inductance matrix.
- **Use cases**: Diagonal entries of K matrix. Needed for flux conservation solver.
- **Strengths**: O(1), simple logarithm. No elliptic integrals.
- **Weaknesses**: Only valid for thin wire (r_wire << a). Constant varies by source (JS uses -2, Python uses -1.75). Inaccurate for thick conductors or finite cross-section.
- **GPU parallelization**: Trivial — O(N), simple arithmetic.
- **Implemented**:
  - JS: `selfInductance(r, wireRadius)` in `js/mhd_demo/physics.js:722-728`. Uses ln(8a/r_w) - 2.
  - Python: `calc_self_inductance(R, r_wire)` in `inductance_core.py:55-70`. Uses ln(8R/r_wire) - 1.75.
  - C++/OpenCL: Not implemented.

#### B3. Flux Conservation Solver (Direct Matrix Solve)

- **Merit**: In a perfectly conducting plasma/cage, magnetic flux through each loop is conserved. At t=0: Φ₀ = K₀·I₀. At time t: solve K_t·I_t = Φ₀ for induced currents. This is the physically correct approach for "fast" plasma expansion (flux frozen into conductors).
- **Use cases**: Computing induced currents in cage and plasma coils when plasma moves. Static equilibrium solver. Dynamic trajectory computation.
- **Strengths**: Physically grounded. Exact (to matrix solve precision). Fast for small N (< 100 coils). Handles fixed-current sources (SC) via partitioned solve: K_uu·I_u = Φ₀ - K_uf·I_fixed.
- **Weaknesses**: O(N³) for direct solve (Gaussian elimination). Requires non-singular K (overlapping loops or r→0 cause singularity). No iterative refinement. JS `gaussianSolve()` skips singular columns (|pivot| < 1e-12) — may silently produce wrong results.
- **GPU parallelization**: Moderate — the matrix assembly (B1, B2) is parallel, but the dense linear solve is not trivially parallel. For small N (< 100), CPU is faster. For large N, could use cuSOLVER/MAGMA. The matrix is dense and symmetric positive definite (could use Cholesky instead of Gaussian elimination).
- **Implemented**:
  - JS: `solveFluxConservation(state)` in `js/mhd_demo/physics.js:779-872`. Uses `gaussianSolve()` (partial pivoting) at `:731-774`. Default solver.
  - Python: `solve_flux_conserving(K1, Phi0, fixed_mask, I_fixed)` in `inductance_core.py:133-186`. Uses `np.linalg.solve` (LAPACK). Supports fixed-current partitioning.
  - Python: `build_inductance_matrix(rs, zs)` in `inductance_core.py:73-99`. Assembles K matrix.
  - C++/OpenCL: Not implemented.

#### B4. Control-Point Iterative Solver (Gradient Descent)

- **Merit**: Instead of flux conservation, specify target B-field values at control points (plasma interior: B=0 diamagnetic; cage surface: B=B_initial flux conservation). Minimize E = Σ(B_target - B_actual)² + λ|I|² via gradient descent.
- **Use cases**: Flexible field shaping where flux conservation is too restrictive. Allows explicit control of field profile at specific locations.
- **Strengths**: Flexible — can impose different boundary conditions at different points. Regularization (λ term) prevents unphysical currents. Iterative — can monitor convergence.
- **Weaknesses**: Only 8 iterations in current implementation — may not converge. Gradient descent is slow (first-order). No line search or trust region. Sensitive to control point placement (singularities near coils). No guarantee of flux conservation.
- **GPU parallelization**: Good — error computation and gradient accumulation are parallel over control points. Each iteration is a matrix-vector product. But N is small (control points + coils < 200), so GPU overhead may dominate.
- **Implemented**:
  - JS: `stepSimulation(state)` in `js/mhd_demo/physics.js:476-686`. The iterative solver is at `:543-573`. Uses `alphaI` (step size) and `lambdaI` (regularization). Precomputes unit fields via `precomputeUnitFields()` at `:470-474`.
  - Python: Not implemented (Python uses direct flux conservation only).
  - C++/OpenCL: Not implemented.

---

### C. Plasma Dynamics (Time Integration)

#### C1. Explicit Euler Integration (Single-Loop Model)

- **Merit**: Simplest time integrator: v += F/m·dt, x += v·dt. At each step: solve flux-conserving currents, compute Lorentz + gas pressure forces, integrate.
- **Use cases**: Dynamic plasma trajectory simulation. Single plasma coil as point mass in (r,z).
- **Strengths**: Trivial to implement. One force evaluation per step. Matches Python reference exactly.
- **Weaknesses**: First-order — energy drifts. No conservation of energy/momentum. Requires small dt for stability (dt=5e-6s, 400 steps = 2ms total). No adaptive timestep. Axis collision handled by crude bounce (v = -0.5·v).
- **GPU parallelization**: Poor for single loop (N=1, no parallelism). For multi-node plasma (N nodes), force computation per node is parallel, but the current solve (B3) is sequential.
- **Implemented**:
  - JS: Dynamic trajectory in `js/mhd_demo/ui.js:277-481`. Single-loop integrator matching Python. Parameters: dt=5e-6, mass=0.1, P₀=1e6, γ=5/3, L_eff=1.0.
  - Python: `run_plasma_dynamics()` in `doc/python/MHD/demo_plasma_dynamics_flux.py:116-235`. Same algorithm, uses `np.linalg.solve` for currents. Also includes magnetic self-pressure (hoop stress) term `Fr_self = mag_self_coeff * I² / r` — **not in JS**.
  - C++/OpenCL: Not implemented.

#### C2. Multi-Node Surface Dynamics (Spring + Pressure + Lorentz)

- **Merit**: Plasma modeled as N nodes connected by springs, with gas pressure normal forces and Lorentz forces per node. More physically realistic than single-loop model.
- **Use cases**: Full plasma surface evolution. Deformation, compression, expansion of plasma boundary.
- **Strengths**: Captures surface deformation. Spring forces maintain topology. Gas pressure drives expansion. Lorentz forces provide magnetic confinement.
- **Weaknesses**: Not validated — `stepSimulation()` multi-node mode is not used for dynamic demo (caused "teleport bug"). Spring constant and damping are ad-hoc. Volume computation via `revolveVolume()` assumes closed surface. No energy conservation check. Explicit Euler with damping factor (0.995) — energy not conserved.
- **GPU parallelization**: Good for force computation (per-node parallel). Current solve is sequential (B3/B4). Spring forces are per-edge parallel. Integration is per-node parallel.
- **Implemented**:
  - JS: `stepSimulation(state)` in `js/mhd_demo/physics.js:617-685`. Multi-node mode (springs + pressure + Lorentz). Currently disabled for dynamic demo (`staticSolver = true` returns before dynamics).
  - Python: Not implemented (Python uses single-loop only).
  - C++/OpenCL: Not implemented.

#### C3. Higher-Order Time Integration (Not Implemented)

- **Merit**: Verlet, Leapfrog, or RK4 integrators that conserve energy (symplectic) or have higher accuracy.
- **Use cases**: Long-time plasma dynamics where energy conservation matters. Production simulations.
- **Strengths**: Energy conservation (symplectic integrators). Higher accuracy (RK4). Larger stable timesteps.
- **Weaknesses**: More complex. More force evaluations per step (RK4: 4 evaluations).
- **GPU parallelization**: Same as C1/C2 — force evaluation is the parallel part.
- **Implemented**: Not implemented. Python docstring suggests "consider Verlet/RK4 for production". Design doc lists "Verlet/Leapfrog" as future work (`MHD_Plasma_Nozzle_webgl.md:182`).

---

### D. Force Computation

#### D1. Lorentz Force on Current Loop

- **Merit**: F = I·L×B where L = 2πr (circumference). For axisymmetric loop: F_r = I·(2πr)·Bz, F_z = -I·(2πr)·Br. Field from all other coils (not self).
- **Use cases**: Plasma dynamics. Force on plasma ring from SC + cage fields.
- **Strengths**: Simple, exact (given accurate B-field). O(N) per target loop (sum over source coils).
- **Weaknesses**: Excludes self-force (hoop stress) in JS version. Python includes separate hoop stress term. Near-coil singularity when source and target are close.
- **GPU parallelization**: Good — each (target, source) pair is independent. O(N²) parallelism.
- **Implemented**:
  - JS: In `ui.js:448-460` (dynamic trajectory) and `physics.js:656-662` (multi-node). Excludes self-force.
  - Python: `compute_lorentz_force_loop(idx_target, rs, zs, I)` in `inductance_core.py:508-537`. Excludes self-force. Python dynamic demo adds `Fr_self = mag_self_coeff * I² / r` separately.
  - C++/OpenCL: Not implemented.

#### D2. Gas Pressure Force (Adiabatic)

- **Merit**: P = P₀·(V₀/V)^γ where V = π·r²·L_eff (cylindrical volume). Force = P·(2πr·L_eff) (outward radial). Models adiabatic expansion/compression of plasma.
- **Use cases**: Plasma expansion driving force. Prevents collapse to r=0.
- **Strengths**: Simple. Physically motivated (adiabatic compression). One scalar computation.
- **Weaknesses**: Simplified geometry (cylinder, not actual plasma shape). L_eff is a free parameter. No axial pressure gradient (Fz=0). γ=5/3 assumes ideal gas.
- **GPU parallelization**: Trivial — O(1).
- **Implemented**:
  - JS: In `ui.js:462-465` (dynamic) and `physics.js:627-653` (multi-node, uses surface area not cylinder).
  - Python: In `demo_plasma_dynamics_flux.py:190-194`. Same cylindrical model.
  - C++/OpenCL: Not implemented.

#### D3. Magnetic Self-Pressure (Hoop Stress)

- **Merit**: F_self ~ I²/r (pinch force). Current loop attracts itself radially inward.
- **Use cases**: Plasma collapse prevention/acceleration. Important for high-current plasmas.
- **Strengths**: Simple. Physically correct scaling.
- **Weaknesses**: Coefficient `mag_self_coeff` is a free parameter (not derived from first principles). Only radial component.
- **GPU parallelization**: Trivial — O(1).
- **Implemented**:
  - Python: `Fr_self = mag_self_coeff * I_p² / max(r_p, r_min)` in `demo_plasma_dynamics_flux.py:197-198`. Default coeff = 1e-7.
  - JS: **Not implemented** — gap. The JS dynamic trajectory does not include hoop stress.
  - C++/OpenCL: Not implemented.

---

### E. Vector Potential & Flux Function

#### E1. Azimuthal Vector Potential A_φ

- **Merit**: A_φ for circular loop and axial dipole. B = curl(A) — can compute B from A via finite differences. Flux Ψ = r·A_φ (field lines are contours of Ψ).
- **Use cases**: Flux surface visualization. Streamline computation. Alternative B-field computation via numerical curl. Diamagnetic bubble solver (find dipole moment such that Ψ_total = 0 on boundary).
- **Strengths**: Exact analytic formula (elliptic integrals). Ψ contours directly show field line topology.
- **Weaknesses**: Same elliptic integral cost. Singularity at r=0 (A_φ → 0 but formula has 1/√r).
- **GPU parallelization**: Good — same as A1, embarrassingly parallel over grid.
- **Implemented**:
  - Python: `Aphi_loop_rz()` in `inductance_core.py:399-439`. `Aphi_dipole_rz()` at `:369-396`. Vectorized.
  - Python: `psi_loop_rz()` at `:483-492`. `psi_dipole_rz()` at `:495-503`. Flux function Ψ = r·A_φ.
  - Python: `compute_numerical_B_from_A()` at `:542-556`. Verifies curl(A) = B via finite differences.
  - JS/GLSL/C++: Not implemented. Gap — no flux function visualization in JS.

#### E2. Dipole Strength Solver (Diamagnetic Bubble)

- **Merit**: Find dipole moment m such that Ψ_total = Ψ_seed + Ψ_dipole = 0 on a boundary. Models a diamagnetic plasma bubble that expels magnetic flux.
- **Use cases**: Diamagnetic plasma equilibrium. Simple alternative to full flux conservation solver.
- **Strengths**: Single scalar unknown (m). Simple root-finding. Physically intuitive (dipole cancels external flux).
- **Weaknesses**: Only models plasma as a single dipole — no surface structure. Assumes spherical boundary. Not a full MHD solution.
- **GPU parallelization**: Poor — single scalar root-find. But evaluation of Ψ on boundary points is parallel.
- **Implemented**:
  - Python: `solve_dipole_strength()` in `demo_dipole_gemini.py:20`. Also used in `demo_diamagnetic_dipole_bubble.py`.
  - JS/C++: Not implemented.

---

### F. Visualization

#### F1. GLSL B-field Shader (HSV Coloring)

- **Merit**: Fragment shader computes B-field per pixel on GPU. Up to 64 coils. HSV color: hue = field direction, value = magnitude. Real-time visualization.
- **Use cases**: Interactive B-field visualization. Background field map.
- **Strengths**: GPU-parallel — O(1) per pixel, millions of pixels in real-time. No CPU bottleneck. Compact code.
- **Weaknesses**: Float precision (10-iteration AGM is marginal). Uniform array limit (`float[256]` may exceed WebGL1 limits on some hardware). No log-scale option. No streamline rendering.
- **GPU parallelization**: Excellent — this IS the GPU implementation. Each pixel is independent.
- **Implemented**:
  - GLSL: Fragment shader in `js/mhd_demo/render.js:110-212`. HSV + grayscale modes. Up to 64 coils.
  - Python: Uses matplotlib contour/quiver (CPU). No GPU visualization.
  - C++/OpenCL: Not implemented.

#### F2. CPU Vector Field Grid

- **Merit**: Sample B-field at N×N grid points on CPU, draw as colored line segments. Fallback when shader is disabled.
- **Use cases**: Debugging visualization. Vector field arrows.
- **Strengths**: Simple. No shader complexity. Works on all hardware.
- **Weaknesses**: O(N² × N_coils) on CPU — slow for large grids. 32×32 = 1024 points only.
- **GPU parallelization**: Could be done with instanced rendering or compute shader, but shader (F1) is better.
- **Implemented**:
  - JS: `updateField(sim)` in `js/mhd_demo/render.js:405-482`. 32×32 grid, LineSegments.
  - Python: matplotlib quiver (in demo scripts).
  - C++/OpenCL: Not implemented.

---

### G. Geometry Generators

#### G1. Parabolic Cage Generator

- **Merit**: Generate cage coil positions along a parabolic wall r² = A·z + B. Equidistant in R (denser near throat).
- **Use cases**: Magnetic nozzle cage geometry. Equidistant-R ensures high field precision near throat.
- **Strengths**: Simple. Physically motivated (parabolic nozzle shape). Dense sampling where field gradients are highest.
- **Weaknesses**: Fixed parabolic shape. No spline or arbitrary wall profile.
- **GPU parallelization**: N/A — too few points (8–12 coils).
- **Implemented**:
  - JS: `makeParabolicCage(rMin, rMax, zStart, zLen, count)` in `js/mhd_demo/physics.js:94-117`.
  - Python: `generate_parabolic_rings()` in `inductance_core.py:202-215`. Also `generate_parabolic_cage_rings()` wrapper.
  - C++/OpenCL: Not implemented.

#### G2. Spherical Plasma Generator

- **Merit**: Generate plasma nodes on a spherical shell in (r,z) via θ parameterization.
- **Use cases**: Initial plasma shape for multi-node solver.
- **Strengths**: Simple. Produces closed surface (arc from pole to pole).
- **Weaknesses**: Only spherical. No ellipsoid, toroid, or arbitrary shape.
- **GPU parallelization**: N/A — too few points.
- **Implemented**:
  - JS: `makeSphericalPlasma(centerZ, radius, count)` in `js/mhd_demo/physics.js:119-147`.
  - Python: `generate_spherical_rings()` in `inductance_core.py:235-255`. Also `generate_spherical_plasma_loops()` wrapper.
  - C++/OpenCL: Not implemented.

#### G3. Other Geometry Generators (Python Only)

- **Merit**: Disk, tube, driver, projectile geometry generators for various MHD configurations (Gauss gun, etc.).
- **Use cases**: Gauss gun simulation, disk armature, tube driver coils.
- **Strengths**: Flexible. Parameterized.
- **Weaknesses**: Python only — no JS counterpart.
- **GPU parallelization**: N/A.
- **Implemented**:
  - Python: `generate_disk_rings()` in `inductance_core.py:263-270`. `generate_tube_rings()` at `:273-280`. `build_driver_tube()` at `:561-572`. `build_projectile_disk()` at `:575-582`. `build_geometry_vectors()` at `:585-595`.
  - JS/C++: Not implemented.

---

### H. Verification & Testing

#### H1. Flux Consistency Check

- **Merit**: Verify that numerical integration of A_φ over a disk gives the same flux as K·I (mutual inductance formula). Cross-checks two independent computations.
- **Use cases**: Code verification. Ensuring B-field and inductance formulas are consistent.
- **Strengths**: Catches sign errors, normalization bugs, elliptic integral parameter confusion (m vs k).
- **Weaknesses**: Not automated in CI. Must be run manually.
- **GPU parallelization**: N/A — verification code.
- **Implemented**:
  - Python: `check_mutual_consistency()` in `check_B_consistency.py:50`. `check_on_axis()` at `:75`. `check_flux_self_inductance()` in `check_B_kernels.py:22`. `check_on_axis_B()` at `:46`.
  - JS: `test_solver.js` only checks for NaN — no accuracy verification. Gap.

#### H2. Numerical curl(A) = B Verification

- **Merit**: Compute B via finite-difference curl of A_φ and compare to direct B-field formula. Verifies that A_φ and B formulas are mutually consistent.
- **Use cases**: Code verification for vector potential and field formulas.
- **Strengths**: Independent computation path. Catches formula errors.
- **Weaknesses**: Finite-difference accuracy limited by step size. Not automated.
- **GPU parallelization**: N/A.
- **Implemented**:
  - Python: `compute_numerical_B_from_A()` in `inductance_core.py:542-556`. `check_Aphi_vs_B()` in `demo_diamagnetic_dipole_bubble.py:42`. `run_Aphi_profile_tests()` at `:73`.
  - JS: Not implemented. Gap.

---

### I. Implementation Gap Summary

| Algorithm | Python | JS | GLSL | C++/OpenCL | Gap |
|-----------|--------|-----|------|------------|-----|
| Elliptic integral B-field (loop) | ✅ (scipy) | ✅ (AGM) | ✅ (AGM, float) | ❌ | C++/OCL missing |
| Point dipole B-field | ✅ | ❌ | ❌ | ❌ | JS/GLSL gap |
| Mutual inductance (Neumann) | ✅ | ✅ | ❌ | ❌ | GLSL/C++ gap |
| Self-inductance (thin-wire) | ✅ (-1.75) | ✅ (-2.0) | ❌ | ❌ | Constant mismatch |
| Flux conservation (direct solve) | ✅ (LAPACK) | ✅ (Gaussian elim.) | ❌ | ❌ | C++/OCL missing |
| Control-point iterative solver | ❌ | ✅ (gradient descent) | ❌ | ❌ | Python gap |
| Explicit Euler (single-loop) | ✅ | ✅ | ❌ | ❌ | C++/OCL missing |
| Multi-node surface dynamics | ❌ | ✅ (untested) | ❌ | ❌ | Python gap, untested |
| Higher-order integration (Verlet/RK4) | ❌ | ❌ | ❌ | ❌ | Not implemented |
| Lorentz force | ✅ | ✅ | ❌ | ❌ | C++/OCL missing |
| Gas pressure (adiabatic) | ✅ | ✅ | ❌ | ❌ | C++/OCL missing |
| Hoop stress (self-pressure) | ✅ | ❌ | ❌ | ❌ | JS missing |
| Vector potential A_φ | ✅ | ❌ | ❌ | ❌ | JS/GLSL/C++ gap |
| Flux function Ψ = r·A_φ | ✅ | ❌ | ❌ | ❌ | JS gap (no flux surface viz) |
| Dipole strength solver | ✅ | ❌ | ❌ | ❌ | JS gap |
| GLSL B-field shader (HSV) | ❌ | N/A | ✅ | ❌ | Python has no GPU viz |
| CPU vector field grid | ✅ (matplotlib) | ✅ | N/A | ❌ | — |
| Parabolic cage generator | ✅ | ✅ | ❌ | ❌ | — |
| Spherical plasma generator | ✅ | ✅ | ❌ | ❌ | — |
| Disk/tube/driver generators | ✅ | ❌ | ❌ | ❌ | JS gap |
| Flux consistency check | ✅ | ❌ (NaN only) | ❌ | ❌ | JS lacks accuracy tests |
| curl(A)=B verification | ✅ | ❌ | ❌ | ❌ | JS gap |

### J. Key Discrepancies Between Python and JS

1. **Self-inductance constant**: Python uses ln(8R/r_w) - 1.75 (Rosa-Kirchhoff), JS uses ln(8a/r_w) - 2 (simpler). This causes systematic differences in flux conservation solutions.
2. **Hoop stress**: Python dynamic includes `Fr_self = coeff * I² / r`, JS does not. JS trajectory will diverge from Python for high currents.
3. **Elliptic integral precision**: Python uses `scipy.special` (machine precision), JS uses 20-iteration AGM (~1e-15), GLSL uses 10-iteration AGM (~1e-6 float). GLSL may show visible artifacts for near-singular field points.
4. **Gaussian elimination**: JS `gaussianSolve()` skips singular columns silently. Python `np.linalg.solve` raises an exception. JS may produce wrong results without warning.
5. **Multi-node model**: JS has a multi-node surface model (springs + pressure + Lorentz) that is not validated and not used for dynamic demo. Python has no multi-node model.
