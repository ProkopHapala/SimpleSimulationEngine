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

## Related Audits

- **`fluid-dynamics.md`** — Overview of all fluid methods (Eulerian, potential flow, Poisson). This audit is the detailed deep-dive into the Eulerian compressible solver.
- **`parallel-particle-cell.md`** — Grid infrastructure (`CubeGridRuler`), cell lists, neighbor search. The Eulerian solver uses fixed Cartesian grids; PIC methods use the same grid infrastructure.
- **`aerodynamics-hydrodynamics.md`** — Potential flow / vortex methods for aerodynamics. Complementary regime (incompressible, irrotational) to the compressible Eulerian solver here.
- **`soft-body-truss-dynamics.md`** — Lagrangian elastic dynamics (truss/mass-spring). Different paradigm: explicit vs. implicit time integration, structural vs. fluid mechanics.
- **`collision-detection.md`** — Broad-phase collision detection (hash grid, sweep-and-prune, AABB tree). Relevant for PIC particle-to-grid interactions.
- **`radiosity-raytracing-scattering.md`** — Radiation transport. Coupled to hydrodynamics in nuclear implosion scenarios (see `doc/python/Burn1D/SphericalImplosion.md`).

---

## Algorithm Review: Heterogeneous Material Simulation Methods

This section reviews algorithms for simulating multiple immiscible substances (heterogeneous materials) — their interfaces, mixing, and phase separation. For each method: merit, use cases, strengths, weaknesses, GPU parallelization, and implementation status in this codebase.

### A. Interface Tracking Methods (Eulerian)

#### A1. Level Set Method (LSM)

- **Merit**: Track material interfaces using a signed-distance function $\phi(x,t)$ where $\phi=0$ is the interface. Volume fraction $\alpha = H_\epsilon(\phi)$ via smoothed Heaviside. Interface normal $n = \nabla\phi/|\nabla\phi|$ is always available. Topological changes (merging/splitting of droplets) handled automatically.
- **Use cases**: Two-material compressible flow with high density ratios (e.g. uranium/hydrogen = 100:1). Hypervelocity impact. Free-surface flows. Any problem where interface geometry matters more than exact mass conservation.
- **Strengths**: Sharp interface representation (zero-crossing is topologically robust). Automatic topology changes. Smooth normals for EOS and surface tension. Simple advection: $\partial\phi/\partial t + u \cdot \nabla\phi = 0$. Well-established theory (Osher & Sethian, 1988).
- **Weaknesses**: **Not mass-conserving** — $\phi$ advection drifts mass. Requires redistancing (re-initialization to signed-distance) every few steps. Redistancing is iterative (5-10 sub-steps). Interface smearing over $\epsilon \approx 1.5\Delta x$. Mass loss/gain over long simulations. Needs separate partial density tracking for mass conservation.
- **GPU parallelization**: Excellent — advection is local stencil, redistancing is local stencil. Both map to tiled kernels with local memory halos. Already implemented.
- **Implemented**:
  - Python+OpenCL: `python/EulerianImpacFluid/EulerianImpacFluid.cl` — `heaviside()`, `smoothed_sign()`, `redistance_phi` kernel (Godunov-type HJ solver). `EulerianImpacFluid.py` — `redistance()` method (5 iterations, $d\tau = 0.5dx$). 8x8 tiles, 2-cell halo, local memory.
  - C++: **Not implemented** — no C++ level set code exists.
  - Doc: `EulerianImpacFluid.md` (2758 lines) — full derivation of LSM, 5-equation model, Ghost Fluid Method, redistancing theory. `doc/Continum_Mechanics_Solver.md` — discusses LSM as interface tracking option.
  - Status: **Active, GPU-only.** The only interface tracking method implemented in the codebase.

#### A2. Volume of Fluid (VOF)

- **Merit**: Track volume fraction $\alpha \in [0,1]$ per cell per material. Advect $\alpha$ conservatively. Interface reconstructed from $\alpha$ contour (e.g. PLIC — Piecewise Linear Interface Calculation). Mass-conserving by construction.
- **Use cases**: Incompressible multiphase flow. Free-surface flows where mass conservation is critical. Breaking waves, droplet dynamics.
- **Strengths**: **Exact mass conservation** (conservative advection). Simple storage (one scalar per material per cell). Well-established (Hirt & Nichols, 1981). No redistancing needed.
- **Weaknesses**: Interface reconstruction is complex (PLIC geometry). Interface normals from $\nabla\alpha$ are noisy (stair-step artifacts). Topology changes require special treatment. $\alpha$ smearing over several cells without compression. No natural signed-distance property.
- **GPU parallelization**: Moderate — advection is parallel, but PLIC interface reconstruction is sequential and geometrically complex. VOF compression fluxes need careful synchronization.
- **Implemented**: **Not implemented**. Discussed in `EulerianImpacFluid.md` as an alternative to level set. The `apply_sharpening_conservative` kernel (VOF compression via $\alpha(1-\alpha)$ flux) is designed in the doc but not implemented in `.cl` file. Gap — VOF would complement level set for mass conservation.

#### A3. Coupled Level Set / VOF (CLSVOF)

- **Merit**: Combine strengths of both: VOF for mass conservation, level set for sharp interface geometry and normals. Two-way feedback: $\phi$ guides interface reconstruction, VOF corrects $\phi$ to match mass distribution. State-of-the-art for multiphase flow (Sussman & Puckett, 2000).
- **Use cases**: High-accuracy multiphase simulation. Problems requiring both mass conservation and smooth interface normals (surface tension, curvature).
- **Strengths**: Mass-conserving (VOF) + sharp interface (LS). Best of both worlds. Accurate curvature from $\phi$. Redistancing can use VOF mass as constraint.
- **Weaknesses**: Complex implementation (two interface representations to maintain). Computational cost of both methods. Coupling logic is non-trivial.
- **GPU parallelization**: Moderate — both components are stencil-based but coupling adds synchronization.
- **Implemented**: **Not implemented**. Discussed extensively in `EulerianImpacFluid.md` — the `correct_phi_from_mass` kernel (feedback: $\phi \leftarrow \phi + \lambda(\alpha - H(\phi))$) is designed in the doc but not implemented. The user's intuition (described in the LLM chat) matches CLSVOF exactly: use $\phi$ for topology, use $\rho$ for mass, couple them with a correction step. Gap — this is the natural next step for the Eulerian solver.

#### A4. Ghost Fluid Method (GFM)

- **Merit**: Treat material interface as a discontinuity. Create "ghost" states on each side of the interface (e.g. ghost uranium in hydrogen region with hydrogen's pressure/velocity but uranium's entropy). Solve Riemann problems separately for each material. Eliminates pressure oscillations at interfaces with extreme density ratios.
- **Use cases**: Extreme density ratios (>100:1). Shock-interface interaction. Hypervelocity impact where mixture EOS fails.
- **Strengths**: No pressure oscillations at interfaces. Each material uses its own EOS. Sharp interface (no smearing). Well-established (Fedkiw et al., 1999).
- **Weaknesses**: Complex implementation (interface detection, ghost state construction). Requires level set for interface location. Multi-material (>2) is difficult. No mass mixing at interface (pure materials only).
- **GPU parallelization**: Moderate — ghost state construction is local but requires interface detection.
- **Implemented**: **Not implemented**. Discussed in `EulerianImpacFluid.md` and listed as an open issue in this audit: "No Ghost Fluid Method — uses mixture EOS instead (simpler but less accurate for extreme density ratios)." The current solver uses a mixture stiffened-gas EOS which works for moderate density ratios but may produce pressure errors at extreme ratios. Gap — GFM would improve accuracy for the uranium/hydrogen impact scenario.

#### A5. Phase Field Method

- **Merit**: Track order parameter $c(x,t) \in [0,1]$ governed by Cahn-Hilliard equation: $\partial c/\partial t + u \cdot \nabla c = \nabla \cdot (M \nabla \mu)$ where $\mu = \delta F/\delta c$ is chemical potential. Interface has finite thickness controlled by $\epsilon$. Based on free energy functional.
- **Use cases**: Micro-scale two-phase flow with surface tension. Phase separation. Cahn-Hilliard dynamics. Materials science (spinodal decomposition).
- **Strengths**: Physically rigorous (thermodynamic free energy). Natural surface tension from chemical potential. Topology changes automatic. Smooth interface.
- **Weaknesses**: **Too diffusive for hypervelocity impact** (explicitly stated in `EulerianImpacFluid.md`). Requires very fine grid to resolve interface thickness $\epsilon$. Fourth-order PDE (Cahn-Hilliard) is computationally expensive. Mass conservation depends on discretization.
- **GPU parallelization**: Moderate — Cahn-Hilliard is stencil-based but 4th-order derivatives need wider stencils.
- **Implemented**: **Not implemented**. Explicitly rejected in `EulerianImpacFluid.md`: "Phase Field Method: This uses variational derivatives of a free energy functional (Cahn-Hilliard equation). It is physically rigorous for micro-scale surface tension but generally too diffusive for hypervelocity impacts." Correct assessment — phase field is for micro-scale, not impact physics.

#### A6. Mass Scavenging / Donor-Acceptor (PBD-style)

- **Merit**: Detect "leaked" material (small $\alpha$ in wrong phase) and conservatively transfer it to nearest "master" cell with high $\alpha$. Inspired by Position-Based Dynamics: project mass back to constraint manifold ($\phi > 0$). Exact conservation of mass, momentum, energy.
- **Use cases**: Cleaning up numerical diffusion of material interfaces. Preventing "fog" of small amounts of material spreading into vacuum. Hypervelocity impact where material fragments must stay coherent.
- **Strengths**: Exact conservation (mass, momentum, energy all transferred). Simple logic (threshold + nearest neighbor). Topological guarantee (all mass within $\Delta x$ of $\phi=0$). Works as post-processing step after any solver.
- **Weaknesses**: Grid bias (tends to align interfaces with grid axes). Kills small physical droplets (atomization). No curvature information. Heuristic — not a PDE, so no convergence theory. Threshold parameter is tunable.
- **GPU parallelization**: Moderate — gather pattern with 3x3 stencil, but requires two passes (leak detection + master accumulation). Race conditions possible if multiple leaks target same master.
- **Implemented**: **Designed but not implemented**. The `scavenge_matter` kernel is fully specified in `EulerianImpacFluid.md` with OpenCL code. Uses `is_leak()` threshold (alpha < 0.05), 3x3 neighbor search for master, conservative mass/momentum/energy transfer. Not in `EulerianImpacFluid.cl`. Gap — would improve interface sharpness for the existing solver.

---

### B. Lagrangian & Hybrid Methods

#### B1. Particle-In-Cell (PIC)

- **Merit**: Material mass, momentum, and identity stored on Lagrangian particles. Pressure/temperature solved on Eulerian grid. Particles deposit (P2G) to grid, grid solves pressure, forces interpolated (G2P) back to particles. O(N) scaling (each particle interacts with $2^D$ grid cells).
- **Use cases**: Impact/armor penetration. Shaped charges. Explosions. Plasma simulation. Any problem with large deformation where pure Eulerian methods suffer from numerical diffusion.
- **Strengths**: Natural multi-material tracking (each particle knows its material). No interface smearing (particles carry material identity). Handles large deformation without remeshing. Mass-conserving (particles carry mass). Well-established (since 1950s).
- **Weaknesses**: Particle noise (statistical fluctuation in P2G). Needs many particles per cell for smooth fields. Particle sorting/deposition is bandwidth-limited. No direct particle-particle interactions (pressure only via grid). Momentum not exactly conserved (interpolation errors).
- **GPU parallelization**: Moderate — P2G deposition requires atomic writes or sorting. G2P interpolation is parallel. Grid solve is parallel. Particle sorting is the bottleneck. See `doc/Parallel_Particle_To_Cell_accumulation.md` (175 lines) for GPU PIC accumulation strategies.
- **Implemented**:
  - C++: `cpp/common/dynamics/MechPIC2D.h` — 2D PIC with bilinear interpolation, ideal gas EOS ($p \propto \rho^{5/3}$), `particlesToCells()`, `moveMD()`. `CompressibleMaterial` struct per material. Self-correction variant (`moveMD_selfCorrect`) subtracts particle's own contribution. **Early prototype** — EOS marked "probably wrong", no shock viscosity, no energy tracking.
  - C++: `cpp/common/dynamics/MechPIC2D_Temperature.h` — Extension with temperature tracking, particle volumes, moles_new buffer. Also early prototype.
  - C++: `cpp/common/dynamics/TEMP/PIC.h` — 3D PIC for visco-elastic material (impact/armor/shaped charges). Uses `Grid3D` + `SoftBody` (bonds). Very early — mostly declarations, `SoftBody` class reused from truss dynamics.
  - C++: `cpp/common/dynamics/CompressiveParticles.h` — Compressible particle method: spherical particles with variable radius, adiabatic compression ($pV^\kappa = const$), pairwise overlap forces. 3D. Distinct from PIC — pure Lagrangian, no grid. Application: impact, explosions, nuclear.
  - Python: `doc/python/Burn1D/implosion_solver.py` — 1D spherical Lagrangian implosion with Numba. Multi-material (Gas/Pusher/Tamper), staggered grid, artificial viscosity, DT fusion, neutron diffusion. Not PIC — pure Lagrangian.
  - OpenCL: **Not implemented** — no GPU PIC kernel exists. `doc/Parallel_Particle_To_Cell_accumulation.md` designs GPU strategies (coherent cloud, scatter-sort-gather).
  - Tests: `tests_bash/Particle_In_Cell/README.md` — documents PIC algorithm (projection, grid-field solver, interpolation) but no test scripts.
  - Status: **Early C++ prototypes, no GPU, no production use.**

#### B2. FLIP (Fluid-Implicit-Particle)

- **Merit**: Extension of PIC where particle velocities are updated by the **change** in grid velocity (not the absolute grid velocity). This preserves particle momentum history and reduces numerical diffusion. $v_p^{new} = v_p^{old} + (v_g^{new} - v_g^{old})$ interpolated to particle.
- **Use cases**: Fluid simulation with complex interfaces. Liquid animation (graphics). Incompressible flow with free surfaces. Problems where PIC is too diffusive.
- **Strengths**: Much less numerical diffusion than PIC. Preserves particle momentum. Good for free-surface flows. Standard in graphics fluid simulation (Zhu & Bridson, 2005).
- **Weaknesses**: Particle noise can accumulate (needs occasional PIC blending). Energy conservation is poor (no implicit dissipation). Not ideal for shocks (designed for incompressible). Requires velocity grid in addition to pressure grid.
- **GPU parallelization**: Same as PIC — P2G and G2P are the bottleneck.
- **Implemented**: **Not implemented**. No FLIP code exists. Would be a natural extension of the existing PIC prototypes if incompressible fluid simulation is needed.

#### B3. Material Point Method (MPM)

- **Merit**: Generalization of PIC for solid mechanics. Particles carry deformation gradient, stress, and material state. Grid solves momentum equation with constitutive model. Handles large deformation, contact, fracture. Modern variants (MLS-MPM, G2P2G) improve accuracy.
- **Use cases**: Solid mechanics with large deformation. Impact/penetration. Fracture and fragmentation. Snow, sand, visco-elasto-plastic materials. Graphics animation.
- **Strengths**: Handles solids AND fluids in unified framework. Natural multi-material (particle material identity). Large deformation without remeshing. Contact handled automatically (particles from different bodies share grid). Modern variants are competitive with FEM for large deformation.
- **Weaknesses**: Memory-heavy (particles carry many state variables). Particle-grid transfer introduces noise. CFL-limited by grid + particle motion. Stress computation per particle is expensive. Not as accurate as FEM for small deformation.
- **GPU parallelization**: Moderate — same P2G/G2P bottleneck as PIC. Modern MLS-MPM maps well to GPU. Stress computation is per-particle (parallel).
- **Implemented**: **Not implemented**. No MPM code exists. The PIC prototypes (`MechPIC2D.h`, `PIC.h`) are the closest but lack deformation gradients and constitutive models. Gap — MPM would be the natural choice for impact/armor/shaped charge simulation with material strength.

#### B4. Smoothed Particle Hydrodynamics (SPH)

- **Merit**: Pure Lagrangian meshless method. Particles carry mass, momentum, energy. Fields approximated by kernel interpolation: $A(r) = \sum_j m_j A_j / \rho_j W(|r-r_j|, h)$. Derivatives via kernel gradients. No grid needed at all.
- **Use cases**: Free-surface flows. Astrophysics (star formation, galaxy collisions). Fluid-structure interaction. Problems with very large deformation and topology changes.
- **Strengths**: No grid — no mesh tangling. Natural multi-material (particle identity). Mass-conserving by construction. Handles free surfaces naturally. Lagrangian — no advection errors.
- **Weaknesses**: **O(N²) without neighbor search** (or O(N) with cell lists). Particle clustering/voids. Tensile instability (artificial stress needed). Boundary conditions difficult. No exact conservation (kernel truncation). Sound speed stiffness limits dt. Many variants (WCSPH, ISPH, DFSPH) — no consensus.
- **GPU parallelization**: Good — particle-particle interactions are parallel with cell-list neighbor search. See `parallel-particle-cell.md` for grid infrastructure. But neighbor list construction is the bottleneck.
- **Implemented**: **Not implemented**. Listed as an open issue in `fluid-dynamics.md`: "No SPH (Smoothed Particle Hydrodynamics) implementation." The grid infrastructure (`CubeGridRuler`, `HashMap2D`) and neighbor search exist in `parallel-particle-cell.md` but no SPH kernel or force computation. Gap — SPH would be relevant for free-surface and astrophysical flows.

#### B5. Lagrangian Finite Volume (Hydrocode)

- **Merit**: Mesh moves with material. Each cell/triangle has fixed mass. Volume changes with deformation. Pressure forces from $\mathbf{F} = -P \nabla V$. Material interfaces are mesh edges — perfectly sharp, no tracking needed. Artificial viscosity for shocks.
- **Use cases**: Hypervelocity impact. Shaped charges. Nuclear implosion. Any problem with sharp material interfaces and strong shocks. Standard in production hydrocodes (LS-DYNA, CTH, ALE3D).
- **Strengths**: **Perfect interface tracking** (material boundaries are mesh edges). Exact mass conservation (fixed mass per cell). Natural multi-material (each cell has one material). Handles shocks with artificial viscosity. Well-established (Von Neumann & Richtmyer, 1950s).
- **Weaknesses**: **Mesh tangling** — cells invert under large shear, requiring remeshing. Remeshing is complex and diffusive (edge flipping mixes materials). 3D remeshing is very hard. Not suitable for fluid mixing or turbulence. Topological changes (fragmentation) require special treatment.
- **GPU parallelization**: Poor — mesh connectivity is irregular, force computation is sparse matrix. Remeshing is inherently sequential. Not commonly GPU-accelerated.
- **Implemented**:
  - C++: `doc/Continum_Mechanics_Solver.md` (367 lines) — full design for 2D cylindrical Lagrangian hydrocode: triangular mesh, Pappus theorem for volume, $\mathbf{F} = P\nabla V$, artificial viscosity (Von Neumann-Richtmyer), edge flipping with conservative mixing, Stiffened Gas / Mie-Grüneisen / Thomas-Fermi EOS. **Design only, no code.**
  - Python: `doc/python/Burn1D/implosion_solver.py` — 1D spherical Lagrangian implosion. Multi-material (Gas/Pusher/Tamper), staggered grid, artificial viscosity (linear + quadratic), DT fusion, neutron diffusion. **Active, 1D only.**
  - Python: `doc/python/Burn1D/lagrangian_tube_solver.py` — 1D Lagrangian tube solver with variable cross-section, 5-species combustion, wall forces. **Active, 1D only.**
  - Python: `doc/python/Burn1D/wave_tube_solver.py` — 1D Lagrangian acoustic wave solver. **Active, 1D only.**
  - C++: `CompressiveParticles.h` — Spherical compressible particles (variable radius, pairwise overlap). Pure Lagrangian, no mesh. **Early prototype.**
  - Status: **1D Python implementations active. 2D/3D C++ is design-only.**

---

### C. Eulerian Multi-Material Formulations

#### C1. 5-Equation Model (Allaire et al.)

- **Merit**: Quasi-conservative multi-material model. Separate mass conservation for each material ($\alpha_k \rho_k$), shared momentum and energy. Volume fraction from level set. Mixture EOS. Avoids pressure oscillations at material interfaces.
- **Use cases**: Two-material compressible flow with high density ratios. Hypervelocity impact. Multi-phase flow.
- **Strengths**: Simple extension of single-material Euler equations. Handles density ratios up to ~100:1 with mixture EOS. Conservative for mass and momentum. Level set provides volume fractions.
- **Weaknesses**: Volume fraction advection is non-conservative (level set drift). Mixture EOS smears interface properties. Not as accurate as Ghost Fluid for extreme ratios. Energy equation can produce oscillations at interfaces.
- **GPU parallelization**: Excellent — all updates are local stencils. Already implemented with tiled local memory.
- **Implemented**:
  - Python+OpenCL: `python/EulerianImpacFluid/` — full implementation. See `continuum-mechanics-impact.md` for details.
  - Status: **Active, GPU-only.** The only multi-material Eulerian model implemented.

#### C2. 7-Equation Model (Baer-Nunziato)

- **Merit**: Full two-phase model with separate velocity, pressure, and energy for each material. 7 equations: 2× mass, 2× momentum, 2× energy, 1× volume fraction. Each material has its own velocity field — handles velocity non-equilibrium at interfaces.
- **Use cases**: Deflagration-to-detonation. Reactive flow with multiple velocities. Granular flow. Problems where materials have different velocities at the interface.
- **Strengths**: Most general multi-material model. Handles velocity and pressure non-equilibrium. No mixture EOS needed (each material has own state). Well-posed hyperbolic system.
- **Weaknesses**: Complex — 7 equations vs 5. Many source terms (relaxation: velocity, pressure, temperature). Stiff relaxation terms need implicit treatment. More memory per cell. Harder to stabilize.
- **GPU parallelization**: Good — stencil-based, but relaxation terms may need iteration.
- **Implemented**: **Not implemented**. Not discussed in any doc. The 5-equation model is simpler and sufficient for the current uranium/hydrogen impact scenario where velocity equilibrium is a reasonable assumption. Gap — would be needed for problems with velocity non-equilibrium (e.g. granular flow, detonation).

#### C3. Diffuse Interface / Mixture EOS

- **Merit**: Simplest multi-material approach. Single velocity, single pressure. Materials share cell via mixture EOS. Volume fraction from level set or mass fraction. No separate Riemann solvers per material.
- **Use cases**: Moderate density ratios (<100:1). Quick prototyping. When interface physics is not critical.
- **Strengths**: Simplest implementation. Single set of conservative variables. No ghost states. Works with any single-material solver.
- **Weaknesses**: Pressure errors at extreme density ratios. Interface smearing. No material-specific shock treatment. Mixture sound speed can be wrong.
- **GPU parallelization**: Excellent — same as single-material solver.
- **Implemented**: **This is what the current solver does.** `EulerianImpacFluid.cl` uses mixture stiffened-gas EOS with volume fractions from level set. Listed as a weakness in the audit: "No Ghost Fluid Method — uses mixture EOS instead."

---

### D. Lattice Methods

#### D1. Lattice Boltzmann Method (LBM)

- **Merit**: Mesoscopic method — solve Boltzmann equation on lattice. Distribution functions $f_i(x,t)$ stream and collide. Recover Navier-Stokes in continuum limit. Natural multi-phase via Shan-Chen pseudo-potential or color-gradient models. Local operations — excellent for GPU.
- **Use cases**: Incompressible / weakly compressible flow. Porous media flow. Micro-channel flow. Multi-phase flow with complex geometry. Problems where boundary geometry is complex but Mach numbers are low.
- **Strengths**: **Excellent GPU parallelization** — all operations are local (stream + collide). Natural multi-phase (Shan-Chen, color-gradient, free-energy). Handles complex boundaries (bounce-back). No Poisson solver needed (pressure from EOS). Simple implementation.
- **Weaknesses**: **Not for compressible flow / shocks** — LBM is limited to low Mach numbers (Ma < 0.3 typically). Lattice artifacts (discrete velocities). Memory-heavy (19 or 27 distribution functions per cell in 3D). Stability issues at high Reynolds. Not suitable for hypervelocity impact or shock physics. Compressible LBM variants exist but are less mature.
- **GPU parallelization**: **Best of all methods** — perfectly local, no global reductions, no neighbor lists. Ideal for GPU. This is why LBM is popular for GPU fluid simulation.
- **Implemented**: **Not implemented.** No LBM code or discussion exists in the codebase. Gap — LBM would be relevant for low-speed multi-phase flow (e.g. spacecraft coolant, biological fluids) but **not for the hypervelocity impact applications** that are the primary focus. The compressible Eulerian solver is the correct choice for impact physics.

---

### E. Interface Sharpening & Mass Correction

#### E1. VOF Compression (Anti-Diffusion)

- **Merit**: Add artificial compression velocity $u_c \propto \nabla\alpha$ to counteract numerical diffusion of volume fraction. Flux $F = u_c \cdot \alpha(1-\alpha)$ — only active at interface (vanishes in pure cells). Conservative.
- **Use cases**: Sharpening smeared interfaces in VOF or coupled LS-VOF. Preventing material "fog" in vacuum.
- **Strengths**: Conservative. Simple flux formulation. Self-limiting ($\alpha(1-\alpha) \to 0$ in pure cells). Compatible with any Eulerian solver.
- **Weaknesses**: Compression strength is tunable (not physical). Can over-compress (distort interface). CFL condition affected by compression velocity. Not a complete solution (only slows diffusion).
- **GPU parallelization**: Good — local stencil flux computation.
- **Implemented**: **Designed but not implemented.** `apply_sharpening_conservative` kernel fully specified in `EulerianImpacFluid.md`. Not in `.cl` file.

#### E2. Level Set Mass Correction (Phi Feedback)

- **Merit**: After advection, correct $\phi$ to match mass distribution: $\phi^{new} = \phi^{pred} + \lambda(\alpha - H(\phi))$. If mass moved in ($\alpha > H(\phi)$), expand $\phi$. If mass moved out ($\alpha < H(\phi)$), shrink $\phi$. Two-way coupling between mass and geometry.
- **Use cases**: Coupling level set (topology) with conservative mass tracking. CLSVOF-style correction.
- **Strengths**: Simple correction. Maintains $\phi$ as mass-consistent. No extra advection (just nudging). Preserves zero-crossing topology.
- **Weaknesses**: Relaxation factor $\lambda$ is tunable. Not exact mass conservation (nudging, not flux). May need several iterations for large mismatches.
- **GPU parallelization**: Excellent — point-wise operation.
- **Implemented**: **Designed but not implemented.** `correct_phi_from_mass` kernel fully specified in `EulerianImpacFluid.md`. Not in `.cl` file. This is the key missing piece for CLSVOF coupling.

---

### F. Implementation Gap Summary

| Method | C++ | Python | OpenCL | Doc | Gap |
|--------|-----|--------|--------|------|-----|
| Level Set Method | ❌ | ❌ | ✅ | ✅ (2758 lines) | C++ gap |
| Volume of Fluid (VOF) | ❌ | ❌ | ❌ | ✅ (designed) | Not implemented |
| CLSVOF (coupled LS+VOF) | ❌ | ❌ | ❌ | ✅ (designed) | Not implemented — natural next step |
| Ghost Fluid Method | ❌ | ❌ | ❌ | ✅ (discussed) | Not implemented — needed for extreme density ratios |
| Phase Field | ❌ | ❌ | ❌ | ✅ (rejected) | Correctly rejected for hypervelocity |
| Mass Scavenging (PBD-style) | ❌ | ❌ | ❌ | ✅ (designed) | Not implemented — would improve interface sharpness |
| VOF Compression | ❌ | ❌ | ❌ | ✅ (designed) | Not implemented |
| LS Mass Correction | ❌ | ❌ | ❌ | ✅ (designed) | Not implemented — key for CLSVOF |
| 5-Equation Model (Allaire) | ❌ | ✅ (driver) | ✅ (kernels) | ✅ | C++ gap |
| 7-Equation Model (Baer-Nunziato) | ❌ | ❌ | ❌ | ❌ | Not implemented or discussed |
| Diffuse Interface / Mixture EOS | ❌ | ✅ | ✅ | ✅ | Current approach — works but limited |
| PIC (Particle-In-Cell) | ✅ (prototype) | ❌ | ❌ | ✅ | Early prototype, no GPU, no production |
| FLIP | ❌ | ❌ | ❌ | ❌ | Not implemented |
| MPM (Material Point Method) | ❌ | ❌ | ❌ | ❌ | Not implemented — natural for impact/armor |
| SPH | ❌ | ❌ | ❌ | ❌ | Not implemented — grid infra exists but no SPH |
| Lagrangian FV (Hydrocode) | ❌ (design) | ✅ (1D) | ❌ | ✅ (367 lines) | 2D/3D C++ is design-only |
| Compressive Particles | ✅ (prototype) | ❌ | ❌ | ❌ | Early prototype, pure Lagrangian |
| LBM | ❌ | ❌ | ❌ | ❌ | Not implemented — not suitable for compressible impact |

### G. Key Discrepancies & Open Opportunities

1. **CLSVOF is designed but not implemented**: The `EulerianImpacFluid.md` chat (2758 lines) designs a complete CLSVOF pipeline: level set advection + `correct_phi_from_mass` (feedback) + `redistance_phi` (topology) + `apply_sharpening_conservative` (VOF compression) + `scavenge_matter` (PBD-style cleanup). Only the level set advection and redistancing are implemented. The mass correction, VOF compression, and scavenging kernels are fully specified but not in the `.cl` file. This is the most impactful gap for the Eulerian solver.
2. **No Ghost Fluid Method**: The current mixture EOS works for moderate density ratios but produces pressure errors at extreme ratios (>100:1). GFM would solve this but requires separate Riemann solves per material at the interface. The 5-equation model is a simpler alternative that works for the current uranium/hydrogen scenario.
3. **PIC is early prototype**: `MechPIC2D.h` has the basic structure (P2G, grid EOS, G2P) but the EOS is marked "probably wrong", there's no shock viscosity, no energy tracking, and no GPU implementation. The `CompressiveParticles.h` is a different approach (pure Lagrangian, no grid) that may be more suitable for impact physics.
4. **No MPM**: MPM would be the natural choice for impact/armor/shaped charge simulation with material strength. It combines PIC's multi-material tracking with solid mechanics constitutive models. The existing PIC prototypes could be extended to MPM by adding deformation gradients and stress computation.
5. **No SPH**: The grid infrastructure (`CubeGridRuler`, `HashMap2D`) and neighbor search exist, but no SPH kernel or force computation. SPH would be relevant for free-surface flows and astrophysical simulations but is less suitable for shock physics.
6. **Lagrangian hydrocode is 1D-only**: The 1D Python implementations (implosion, tube, wave) are active and useful. The 2D/3D C++ design (`Continum_Mechanics_Solver.md`) is detailed (triangular mesh, edge flipping, artificial viscosity, EOS) but not implemented. This is a significant effort but would provide the most accurate multi-material simulation for impact physics.
7. **No LBM**: Correctly absent — LBM is not suitable for compressible flow with shocks. Would only be relevant for low-speed multi-phase applications (coolant, biological).
8. **No 7-equation model**: The 5-equation model assumes velocity equilibrium. For problems with velocity non-equilibrium (granular flow, detonation), the 7-equation Baer-Nunziato model would be needed. Not currently relevant for the uranium/hydrogen impact scenario.
9. **GPU PIC accumulation designed but not implemented**: `doc/Parallel_Particle_To_Cell_accumulation.md` (175 lines) designs two GPU strategies (coherent cloud with shared memory patch, scatter-sort-gather with radix sort) for PIC particle-to-grid deposition. Neither is implemented. This is the key infrastructure needed for GPU PIC/MPM.
10. **Burn1D is separate from EulerianImpacFluid**: The 1D Lagrangian solvers (`implosion_solver.py`, `lagrangian_tube_solver.py`) and the 2D Eulerian solver (`EulerianImpacFluid`) are completely separate codebases with no shared infrastructure. They solve different dimensionalities and use different formulations (Lagrangian vs Eulerian). A unified multi-material framework would allow comparison and cross-validation.
