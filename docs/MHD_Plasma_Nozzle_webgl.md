

I would like to make simple demo solving magnetohydrodynamics in web browser writen in javascript. It should also solver potantial flow aerodynamics as it is mathematically closely related. I will give you two references which can be very usefull 

1) for reference how to implement the web page and use WebGL look here:
/home/prokop/git/SimpleSimulationEngine/js/spacecraft_editor
it is rather complex project for 3D modeling meshes. Our project is simpler but should reuse main components. Notice that many shared components are here
/home/prokop/git/SimpleSimulationEngine/js/common_js
in particular
@Vec3.js @MeshRenderer.js @GUIutils.js 

2) We have rather comprehesive desing document about how to implement this from numerical and physical perspective, read it carefully
@MHD_Plasma_Nozzle_Simulation.md 

we also have as a reference some function for differnt basis function in potential flow, mainly the eliptical integral to compute field of circular curent-loop (coil) in any point in space
@PotentialFlow.md 
@PotentialFlow.h 
@elliptic_integral.c 

please before you start write markdown with comprehesive plan which I can review

---

## Goals & Scope
- Axisymmetric (r–z) magnetohydrodynamics demo in WebGL (compute on CPU JS; render via WebGL).
- Potential flow analogies for verification (vortex/coil fields, elliptic-integral reference).
- Reuse shared components from `js/common_js` and patterns from `js/spacecraft_editor`.

## Architecture Overview
1. **Data & Math Core (JS)**
   - Vector math via `Vec3.js` (extend with minimal helpers if needed).
   - Field evaluators: magnetic field of circular loop (use elliptic integral reference in `elliptic_integral.c`/`PotentialFlow.h`), potential-flow basis (sources/vortices).
   - Particle/segment structs for plasma shell rings.
2. **Simulation Layer**
   - Axisymmetric ring model: plasma surface as ring nodes with springs (no remeshing), cage rings, superconducting rings.
   - Control-point system (Method 2: overdetermined fit with regularization) to solve induced currents.
   - Time stepping: explicit/semi-implicit Verlet or simple leapfrog; per-step iterative current relaxation (gradient descent) to avoid matrix inversion.
   - Forces: Lorentz (I×B), gas pressure (uniform ideal gas, volume from ring polygon of revolution), springs for node spacing.
3. **Rendering & UI**
   - Use `MeshRenderer.js` for simple line/mesh rendering of rings, field samples.
   - UI with `GUIutils.js`: sliders for coil currents, learning rate, regularization, pressure, time step; toggles for showing field lines/force vectors.
   - 2D axisymmetric view (r–z) plus optional 3D revolve preview (low-poly torus bands).
4. **Integration with `spacecraft_editor` patterns**
   - Reuse shader setup, camera controls, render loop structure.
   - Lightweight scene graph: arrays of drawable primitives (lines/points), updated each frame.

## Physics Kernels (MVP)
- **Magnetic field of a loop**: port compact JS version of coil field using elliptic integrals (K,E) parameterized by m=k²; validate on-axis case vs analytic.
- **Superposition**: sum over SC coils (fixed I), cage coils (solved I), plasma rings (solved I).
- **Control-point solve** (overdetermined):
  - Build A (B at control points from unit currents).
  - Targets: zero B inside plasma, flux-preserving or zero-penetration near cage.
  - Iterative update: I ← I - α * (Aᵀ·error) - λ·I (regularization).
- **Volume & pressure**: compute toroidal volume from ring polygon; P = P0*(V0/V)^γ.
- **Forces**:
  - Magnetic: F_r = I·2πr·B_z, F_z = I·2πr·(-B_r).
  - Gas pressure: normal from adjacent segments; force ∝ pressure * ring area.
  - Springs: keep arc-length near rest length.
- **Integration**: per node v,p update; clamp dt, optional damping.

## Rendering Plan
- Field preview: sample grid in r–z; map |B| to color; optional streamline seeds.
- Geometry: draw rings (plasma, cage, SC), control points, force vectors (scaled).
- UI panels: presets (Helmholtz pair, nozzle), play/pause, step, reset.

## MVP Milestones
1. **Setup & reuse**
   - Bootstrap page with render loop using `MeshRenderer` + camera controls.
   - UI scaffolding with `GUIutils`.
2. **Field kernel**
   - Implement coil B-field (elliptic integral) in JS; unit tests vs on-axis formula.
   - Grid sampling + visualization.
3. **Control-point solver**
   - Hardcode simple cage + plasma ring positions; iterative current fit with regularization; visualize resulting B on grid.
4. **Dynamics**
   - Add plasma ring nodes with springs + pressure; Lorentz forces from solved currents; integrate motion.
   - Basic stability tuning: α, λ, dt, spring k.
5. **Potential-flow mode**
   - Reuse basis functions to show analogy (vortex/sink sources), toggle mode for verification.
6. **Polish**
   - Presets, parameter UI, tooltips, simple logging overlay for convergence/error.

## File/Code Organization
- `js/mhd_demo/`
  - `main.js` (init, loop, UI hooks)
  - `physics.js` (fields, solver, forces, integration)
  - `render.js` (meshes/lines, field textures)
  - `ui.js` (GUI setup)
  - Reuse/import from `js/common_js/Vec3.js`, `MeshRenderer.js`, `GUIutils.js`.

## Validation Checklist
- On-axis B matches analytic coil formula.
- Field symmetry about axis holds numerically.
- Currents remain bounded with regularization when plasma nodes get close.
- Volume/pressure evolution consistent with γ-law; energy not obviously exploding for stable dt.

## TODO (update during work; check only after your visual confirmation)
- [ ] Create `js/mhd_demo/` scaffold (main loop, renderer hookup, UI skeleton reusing Vec3/MeshRenderer/GUIutils).
- [ ] Implement coil magnetic field kernel in JS using elliptic integrals; validate on-axis case (initial implementation, basic).
- [ ] Add field sampling grid and visualization in r–z view.
- [ ] Build control-point iterative solver (overdetermined fit + regularization) for cage/plasma currents (initial gradient-descent fit).
- [ ] Add plasma ring model with springs + pressure; compute toroidal volume and gas law update.
- [ ] Integrate forces (Lorentz, pressure, springs) with stable timestep; basic damping.
- [ ] Hook UI controls (currents, learning rate, regularization, pressure, dt, toggles) and presets.
- [ ] Add potential-flow mode toggle sharing basis functions for verification.
- [ ] Polish visuals (force vectors, field magnitude colormap, simple logging overlay).
