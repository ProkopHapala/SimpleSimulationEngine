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
    - $\Phi_{initial}$: Flux at $t=0$ (usually 0 for plasma, or SC-only field for cage).
    - $\Phi_{external}$: Flux from fixed Superconductor coils.
  - **Solver**: Gaussian elimination with partial pivoting.
- **Pros**: Fast, robust, physically grounded for "fast" expansions.

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
  - **Visualization**:
    - **Asymmetric Vector Field**: Grid points (dots) with vectors extending along field lines. Length $\propto |B|$.
    - **Coil Geometry**: Renders both 3D rings (for perspective) and Profile lines (for r-z shape).
    - **Symmetric View**: Covers full $r \in [-2.5, 2.5]$ range.
- **User Interface**: `ui.js`.
  - Side-panel layout (non-overlapping).
  - Controls: State toggle ($t=0/1$), currents, radii, solver method selection, verbosity.

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

## MVP Status Checklist

### Core Features (Completed)
- [x] **Project Structure**: Scaffolded `js/mhd_demo/` with separation of concerns (Physics, Render, Main, UI).
- [x] **B-Field Kernel**: Implemented elliptic integral solver for circular loops.
- [x] **Visualization**:
    - [x] Symmetric r-z grid sampling.
    - [x] Magnitude-scaled vector field.
    - [x] CSS2D Text labels for coil currents.
    - [x] Separate side-panel UI.
- [x] **Solvers**:
    - [x] **Flux Conservation**: Matrix inverse approach (Method 1).
    - [x] **Control Points**: Iterative regularized fit (Method 2).
- [x] **Geometry Modeling**:
    - [x] Parabolic Cage (radially sampled for vertex density).
    - [x] Spherical Plasma (nodes for $t=0$ and $t=1$).
- [x] **Interaction**:
    - [x] State switching ($t=0 \leftrightarrow t=1$).
    - [x] Parametric input for geometry (radius, focus, etc.).

### Next Phase (Planned)
- [ ] **Full Dynamics**:
    -   Implement forces: Lorentz ($J \times B$) and Gas Pressure ($P = P_0 (V_0/V)^\gamma$).
    -   Time integration (Verlet/Leapfrog) for node motion.
- [ ] **Potential Flow Mode**: Toggle to switch math kernel to fluid sources/sinks.
- [ ] **Advanced Interaction**: Drag-to-move coils/nodes.

## Reference Links
- [Spacecraft Editor](file:///home/prokop/git/SimpleSimulationEngine/js/spacecraft_editor/js/SpaceCraftRenderer.js) (Reference for rendering patterns)
- [Vector Math](file:///home/prokop/git/SimpleSimulationEngine/js/common_js/Vec3.js)
