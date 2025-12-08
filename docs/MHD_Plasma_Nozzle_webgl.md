# MHD Plasma Nozzle WebGL Demo – Design & Implementation Document

## Goals & Scope
- **Axisymmetric (r–z) Magnetohydrodynamics**: Compute plasma equilibrium in a magnetic nozzle configuration using JS on CPU, visualized via WebGL.
- **Verification Analogies**: Compare results with potential flow analogies (vortex/coil fields).
- **Interactive Web Demo**: Real-time control of coil currents, geometry, and solver parameters.

## Architecture Overview

### 1. Data & Math Core (JS)
- **Vector Math**: `Vec3.js` from `common_js`.
- **Field Evaluators**:
  - `coilField(r, z, I, a)`: Computes magnetic field of a circular loop using elliptic integrals (Carlson symmetric form via AGM).
  - Validated against analytic on-axis formulas.
- **Simulation State**:
  - Managed in `physics.js`.
  - Supports two distinct states: `t=0` (Initial small plasma) and `t=1` (Expanded plasma).

### 2. Simulation Layer (Solvers)
Two solver strategies were implemented to determine induced currents in the plasma and cage coils:

#### **Method 1: Flux Conservation (Direct Matrix Solve)**
- **Principle**: In a perfect conductor (or high-conductivity plasma), magnetic flux through any closed loop is conserved.
- **Implementation**:
  - `mutualInductance(r1, z1, r2, z2)`: Computes generic mutual inductance between loops.
  - `selfInductance(r, a)`: Computes self-inductance of a loop.
  - **System**: Moves flux conservation equation $\Phi_{total} = \Phi_{initial}$ to linear system $M \cdot I = \Phi_{initial} - \Phi_{external}$.
    - $M$: Mutual/Self inductance matrix of unknown coils (Cage + Plasma).
    - $\Phi_{initial}$: Flux at $t=0$ computed from SC coil only (cage/plasma currents = 0).
    - $\Phi_{external}$: Flux from fixed Superconductor coils at current positions.
  - **Solver**: Gaussian elimination with partial pivoting.
  - **Physics**: At $t=0$, only SC has current. As plasma expands, induced currents in cage/plasma conserve the initial flux.
- **Pros**: Fast, robust, physically grounded for "fast" expansions. **Now the default solver.**

#### **Method 2: Control Points (Iterative Regularized Fit)**
- **Principle**: Find currents that satisfy specific magnetic field constraints at chosen spatial locations.
- **Implementation**:
  - **Control Points**: Generated strategically:
    - **Axis**: $B_r = 0$ (symmetry).
    - **Vertex**: $B_r = 0$, $B_z \approx 0$.
    - **Plasma Focus**: $B = 0$ (diamagnetism).
    - **Cage Surface**: $B \approx B_{initial}$.
    - **Plasma Interior**: $B = 0$ (diamagnetism).
  - **Algorithm**: Gradient descent to minimize error $E = \sum (B_{target} - B_{actual})^2 + \lambda |I|^2$.
- **Pros**: Flexible, allows shaping field profile explicitly.

### 3. Rendering & UI
- **Renderer**: `DemoRenderer` class in `render.js`.
  - **WebGL**: `THREE.WebGLRenderer` for 3D lines (coils) and field visualization.
  - **CSS2D**: `THREE.CSS2DRenderer` for text labels (`I=50.0kA`) overlay.
  - **GLSL Shader**: GPU-accelerated B-field visualization with HSV coloring (hue = direction, value = magnitude).
  - **Visualization**:
    - **Asymmetric Vector Field**: Grid points (dots) with vectors extending along field lines. Length $\propto |B|$.
    - **Coil Geometry**: Renders both 3D rings (for perspective) and Profile lines (for r-z shape).
    - **Symmetric View**: Covers full $r \in [-2.5, 2.5]$ range.
- **User Interface**: `ui.js`.
  - **Tabbed Layout**: "Settings" and "Editor" tabs.
  - **Settings Tab**: State toggle ($t=0/1$), currents, radii, solver method selection, verbosity, trajectory slider, auto-update checkbox.
  - **Editor Tab**: Two textareas for manual coil editing in `TYPE R0 Z0 R1 Z1 I0` format.
  - **Trajectory Control**: Slider interpolates between $t=0$ and $t=1$, with pre-computed trajectory data.
  - **Energy Plot**: Real-time canvas plot showing magnetic energy and total current vs. time.

## Implementation Details & Insights

### Key Takeaways
1.  **Coordinate Systems in Physics vs. Rendering**:
    -   **Physics**: Axisymmetric $(r, z)$ where $z$ is the symmetry axis.
    -   **Three.js**: $y$ is typically "up". We mapped Physics $z \to$ World $x$, Physics $r \to$ World $y$.
    -   **Insight**: Consistent naming (`r`, `z`) in physics code is crucial. Converting only at the last mile (rendering) prevents confusion.
    -   **Correction**: Initial field visualization was asymmetric/shifted because `coilField` assumed coils at $z=0$. FIX: Always pass relative $z$ ($z_{point} - z_{coil}$).

2.  **Library Assumptions (`Draw3D`)**:
    -   `Draw3D.drawCircleAxis` expects a **unit vector** for the starting direction. Passing a vector of length $R$ resulted in $R^2$ scaling, making rings visible as "cones".
    -   **Fix**: Normalized inputs to helper functions.

3.  **Solver Stability**:
    -   **Method 1 (Flux)** requires a non-singular M matrix. Overlapping nodes or $r \to 0$ can cause singularity. Added $\epsilon$ checks.
    -   **Method 2 (Control Points)** requires careful placement of points. Placing points *exactly* on current loops causes singularities ($1/d$ distance).
    -   **Strategy**: Place control points slightly offset from wire locations (e.g., `r - epsilon`).

4.  **Visualization as Debugging**:
    -   Visualizing the **Grid Points** as dots and vectors starting from them immediately revealed symmetry issues that arrow helpers masked.
    -   **Text Labels** on coils ($I=...$) were essential to verify if the solver was actually doing anything (e.g., distinguishing $0$ from $10^{-6}$).

5.  **Flux Conservation Implementation**:
    -   Critical insight: $\Phi_0$ must be computed at **initial positions** ($t=0$), not current positions.
    -   At $t=0$: Only SC coil has current, cage/plasma have zero current.
    -   At time $t$: Solve $M_t \cdot I_t = \Phi_0 - \Phi_{SC,t}$ where $M_t$ is inductance matrix at current positions.
    -   This matches the Python reference implementation exactly.

6.  **Cage Geometry - Equidistant R**:
    -   Switched from equidistant-Z to equidistant-R sampling.
    -   Formula: $z = z_{start} + A(r^2 - r_{min}^2)$ where $A = (z_{end} - z_{start})/(r_{max}^2 - r_{min}^2)$.
    -   Result: Dense coil spacing near the throat (vertex), sparse at exit.
    -   **Why**: Throat region needs high field precision for accurate flux shielding.

## MVP Status Checklist

### Core Features (Completed)
- [x] **Project Structure**: Scaffolded `js/mhd_demo/` with separation of concerns (Physics, Render, Main, UI).
- [x] **B-Field Kernel**: Implemented elliptic integral solver for circular loops.
- [x] **Visualization**:
    - [x] Symmetric r-z grid sampling.
    - [x] Magnitude-scaled vector field.
    - [x] CSS2D Text labels for coil currents.
    - [x] Separate side-panel UI with tabs.
    - [x] GLSL shader for GPU-accelerated field rendering (HSV mode).
- [x] **Solvers**:
    - [x] **Flux Conservation**: Matrix inverse approach (Method 1) - **DEFAULT**.
    - [x] **Control Points**: Iterative regularized fit (Method 2).
- [x] **Geometry Modeling**:
    - [x] Parabolic Cage (equidistant-R sampling for vertex density).
    - [x] Spherical Plasma (nodes for $t=0$ and $t=1$).
- [x] **Interaction**:
    - [x] State switching ($t=0 \leftrightarrow t=1$).
    - [x] Parametric input for geometry (radius, focus, etc.).
    - [x] Trajectory slider with pre-computed interpolation.
    - [x] Auto-update checkbox for automatic re-solving.
    - [x] Mouse wheel support on all numeric inputs.
    - [x] Manual coil editing via textareas.
    - [x] Energy/current plot canvas.

### Next Phase (TODO)

#### High Priority
- [ ] **Trajectory Enhancements**:
    - [ ] Add user-configurable number of steps (currently hardcoded to 20).
    - [ ] Display current/energy values on hover over energy plot.
    - [ ] Export trajectory data as CSV.
- [ ] **GLSL Shader Refinements**:
    - [ ] Add toggle for HSV vs. magnitude-only coloring.
    - [ ] Implement log-scale option for better dynamic range.
    - [ ] Add colorbar/legend overlay.
- [ ] **UI/UX Polish**:
    - [ ] Syntax highlighting in coil editor textareas.
    - [ ] Validation feedback for invalid coil definitions.
    - [ ] Responsive layout for mobile/tablet.
    - [ ] Add tooltips explaining each parameter.
- [ ] **Performance**:
    - [ ] Profile trajectory computation, optimize matrix solver.
    - [ ] Consider WebWorker for background solving.
    - [ ] Add progress indicator for long computations.

#### Medium Priority
- [ ] **Advanced Visualization**:
    - [ ] Streamlines/field lines overlay.
    - [ ] Flux surface contours.
    - [ ] Pressure contour visualization (requires pressure model).
- [ ] **Export/Import**:
    - [ ] Save/load full simulation state as JSON.
    - [ ] Export screenshots/animations.
    - [ ] Generate shareable URLs with encoded parameters.
- [ ] **Testing**:
    - [ ] Unit tests for `solveFluxConservation`.
    - [ ] Unit tests for `makeParabolicCage`.
    - [ ] Integration tests for UI parameter bindings.
    - [ ] Regression tests comparing to Python reference.

#### Low Priority / Future
- [ ] **Full MHD Dynamics**:
    - [ ] Implement Lorentz force ($\vec{J} \times \vec{B}$).
    - [ ] Add gas pressure model ($P = P_0 (V_0/V)^\gamma$).
    - [ ] Time integration (Verlet/Leapfrog) for node motion.
- [ ] **Potential Flow Mode**: Toggle to switch math kernel to fluid sources/sinks.
- [ ] **Advanced Interaction**: Drag-to-move coils/nodes with mouse.
- [ ] **Multi-Physics**: Couple to thermal model, radiation losses.

---

## Detailed Implementation Notes

### Flux Conservation Solver
The flux conservation solver now correctly matches the Python reference (`demo_coil_motion_flux.py`):

1. **Initial State ($t=0$)**:
   - Only the SC (seed) coil carries current $I_{SC}$.
   - Cage and plasma currents are zero.
   - Compute initial flux through each loop: $\Phi_{0,i} = \sum_{SC} M_{i,SC} \cdot I_{SC}$.

2. **At Time $t$**:
   - Plasma nodes interpolate: $\vec{r}(t) = (1-t)\vec{r}_0 + t\vec{r}_1$.
   - Build inductance matrix $M_t$ at current positions.
   - Compute SC contribution at current positions: $\Phi_{SC,t,i} = \sum_{SC} M_{i,SC}(t) \cdot I_{SC}$.
   - Solve: $M_t \cdot I_t = \Phi_0 - \Phi_{SC,t}$.

3. **Result**:
   - At $t=0$: Cage/plasma currents are zero (as expected).
   - As plasma expands: Induced currents appear to conserve flux.
   - Cage currents shield the SC from changing plasma field.

### Cage Generator (Equidistant-R)
The parabolic cage is now generated with **equidistant radius** sampling:

```javascript
const A = (zEnd - zStart) / (rMax*rMax - rMin*rMin);
for (let i = 0; i < nCage; i++) {
    const t = i / (nCage - 1);
    const r = rMin + (rMax - rMin) * t;  // Equidistant in R
    const z = zStart + A * (r*r - rMin*rMin);
    // ...
}
```

This ensures dense coil spacing near the throat ($r = r_{min}$), which is critical for accurate field shaping.

### GLSL Shader
The B-field shader computes the magnetic field on the GPU for every pixel:

- **Uniforms**: Coil positions/currents packed as `float[256]` (64 coils × 4 floats).
- **Field Calculation**: Elliptic integrals approximated using AGM method (10 iterations).
- **Coloring**: HSV mode maps field direction to hue, magnitude to brightness.
- **Saturation**: `b-max` parameter controls the reference field strength for normalization.

### Auto-Update System
When the "Auto-update" checkbox is enabled:

1. Any parameter change triggers `autoSolveIfEnabled()`.
2. This calls `computeTrajectory(20)` to pre-compute 21 states from $t=0$ to $t=1$.
3. Each state stores: positions, currents, energy, total current.
4. The trajectory slider looks up the nearest pre-computed state.
5. Energy plot and shader update automatically.

### Mouse Wheel Input
All numeric inputs support mouse wheel:

```javascript
el.addEventListener('wheel', (e) => {
    e.preventDefault();
    const step = parseFloat(el.step) || 0.1;
    const delta = e.deltaY < 0 ? step : -step;
    el.value = (parseFloat(el.value) + delta).toFixed(4);
    el.dispatchEvent(new Event('change'));
}, { passive: false });
```

### Energy Plot
The canvas plot shows:
- **Green line**: Magnetic energy $W = \frac{1}{2} \sum L_{ij} I_i I_j$.
- **Orange line**: Total plasma current $I_{total} = \sum I_{plasma}$.
- **White vertical line**: Current $t$ position marker.

---

## Reference Links
- [Python Reference Implementation](../doc/python/MHD/demo_coil_motion_flux.py)
- [Spacecraft Editor](file:///home/prokop/git/SimpleSimulationEngine/js/spacecraft_editor/js/SpaceCraftRenderer.js) (Reference for rendering patterns)
- [Vector Math](file:///home/prokop/git/SimpleSimulationEngine/js/common_js/Vec3.js)
