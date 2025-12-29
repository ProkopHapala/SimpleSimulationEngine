# Aerodynamics & Potential Flow – JavaScript Plan (WIP)

This document outlines a plan to **reimplement selected aerodynamics / potential-flow demos and apps in JavaScript**, primarily inside `js/spacecraft_editor`.

The focus is on:

- Potential-flow around wings / lifting surfaces.
- Building aircraft from panels (low-res and UV-generated surfaces).
- Visualizing velocity fields and streamlines.

---

## 1. Reference C++ Demos and Headers

These are the main C++ sources we mirror conceptually:

- **Potential flow / vortex lattice:**
  - `test_VortexLattice` (C++ app)
  - `test_VortexLattice.cpp`
  - `PotentialFlow.h`

- **Aircraft geometry & design:**
  - `AeroCraftDesign.h` – higher-level aircraft layout, wings, tails, etc.
  - `AeroSurf.h` – 3D aerodynamic surfaces built from a small number of panels per wing.
  - `AeroSurf2D.h` – 2D polars (lift/drag vs angle of attack), used by 2D sail/airfoil sims.

- **Apps using these pieces:**
  - `AeroCombat_main.cpp` – aircraft / combat simulator using `AeroSurf`.
  - `AeroCombatOGL3.cpp` – SDL/OpenGL front-end for AeroCombat, good for understanding rendering & controls.
  - `SailWar_main.cpp` – 2D sail/airfoil game using `AeroSurf2D.h` and polar curves.

- **Mesh generation helpers (airfoil geometry):**
  - `DrawOGL3.h` – especially:
    - `Teardrop2Mesh(...)`
    - `NACASegment2Mesh(...)`
  - These show how airfoil / teardrop cross-sections are turned into meshes (segments, UV surfaces).

All JS work should conceptually sit on top of these ideas, but implemented in `MeshBuilder` / JS math.

---

## 2. Geometry: Building Aircraft from Panels

Goal: build wings / aircraft as **collections of surface panels** that can be used by a potential-flow solver.

### 2.1 Low-res panel wings (AeroSurf-style)

- Mirror `AeroSurf.h` approach:
  - Each wing is represented by **a few coarse panels** (e.g. root/mid/tip segments), each panel carrying:
    - Corner points in 3D.
    - Local normal, area, and possibly control-point positions.
  - Good for fast prototypes and visualizing basic lift distribution.

- JS implementation idea:
  - Add a simple `AeroSurfJS` structure in `spacecraft_editor`:
    - Store a list of wing panels (quads/triangles) with geometric data.
    - Provide helpers to add wings via a few parameters (span, chord, sweep, dihedral, twist).

### 2.2 High-res wings via UV surfaces (DrawOGL3-style)

- Use ideas from `DrawOGL3.h`:
  - `Teardrop2Mesh` and `NACASegment2Mesh` turn airfoil profiles into **meshes along a parametric span**.
  - In JS, we generate:
    - A UV grid over a wing surface (chord × span).
    - Vertices from an airfoil profile (e.g. NACA) swept along the span with twist/dihedral.

- JS implementation idea:
  - In `MeshesUV` / `MeshBuilder`:
    - Add a `WingUV` generator that takes:
      - Airfoil function (e.g. NACA 4-digit in JS).
      - Span, root/tip chords, sweep, dihedral, twist distributions.
      - Subdivision counts (nChord, nSpan).
    - Emit:
      - A triangle/quad mesh for visualization.
      - A derived list of **panels** (one per quad) with centroids/normals for potential-flow.

### 2.3 Emitting panels from surfaces

Common step for both low-res and UV wings:

- Given a mesh (quads/triangles) on the wing:
  - Choose a **panel definition** (usually quads):
    - For each quad: store vertices, centroid, normal, area.
  - Optionally merge small quads into larger effective panels for solver stability.

- JS side:
  - Implement a helper that converts selected faces in `MeshBuilder` into a `PanelSet` structure used by the solver.

---

## 3. Potential Flow System: Formulation and Solution

Objective: solve for potential / circulation strengths such that **normal velocity at sample points near the surface is zero** (no flow through the wing).

### 3.1 Basis functions

- Choose a simple basis compatible with JS performance:
  - **Panel-based sources/vortices:** one unknown strength per panel (or per panel edge/line).
  - Possibly separate **source** (divergence) and **vortex** (circulation) components for thickness vs lift.

- JS structure:
  - `BasisFunction` objects associated with panels/edges, each with:
    - Position/line definition.
    - Influence function to evaluate induced velocity at a point.

### 3.2 Sample points and boundary condition

- For each panel we generate one or more **sample points** slightly offset from the surface along the normal.
- Boundary condition: at each sample point `p_i`:
  - `(v_free(p_i) + v_induced_from_basis(p_i)) · n_i = 0`
  - Where `v_free` is freestream + any prescribed background flow.

- There should be **more sample points than basis functions** for over‑constrained, more stable systems.

### 3.3 Linear system assembly

- Write the system as `A * x = b`:
  - `x`: unknown strengths of basis functions.
  - `A[i,j]`: normal component at sample `i` due to basis `j`.
  - `b[i]`: `- (v_free(p_i) · n_i)`.

- JS implementation:
  - Small dense matrices at first (hundreds of unknowns) to keep things simple.
  - Provide two solution paths:
    - **Direct:** Gauss elimination / LU (for small test cases).  
    - **Iterative:** simple Jacobi / Gauss–Seidel / conjugate-gradient if we want to push resolution.

- Numerical goals:
  - Keep everything simple and visual; precision is less important than intuitive behavior and interactivity.

---

## 4. Velocity Field Evaluation and Rendering

Once we have basis strengths, we can evaluate the induced velocity field anywhere in space.

### 4.1 CPU-side field sampling and streamlines

- Precompute velocity samples on a **3D grid or set of probe lines**:
  - For visualization of cross-sections, downwash behind wings, etc.

- Streamlines (pure JS, CPU):
  - Emit polyline trajectories:
    - Start seeds upstream or near interesting regions.
    - Integrate `dx/dt = v(x)` using fixed or adaptive step.
  - Store as vertex buffers in `MeshBuilder` and render as lines.

### 4.2 Shader-based particle advection

- Alternative, more dynamic visualization:
  - Use **particles** moving in the flow field in a vertex shader or simple GPU step.
  - Each particle emits a short trail (polyline) that approximates a streamline.

- Implementation sketch:
  - Represent the flow field as:
    - Either direct evaluation of basis functions in the shader (small number of basis elements).  
    - Or a 3D texture / grid of velocities precomputed on the CPU.
  - In each frame, in the vertex shader:
    - Read previous particle position.
    - Sample velocity field.
    - Advance position by small time step.

### 4.3 UI hooks in `js/spacecraft_editor`

- Add a dedicated **PotentialFlow / Aero panel**:
  - Parameters: freestream velocity, angle of attack, wing geometry presets.
  - Buttons to:
    - Rebuild wing geometry (low-res/UV).
    - Rebuild panel set and sample points.
    - Assemble and solve potential-flow system.
    - Toggle visualization modes (streamlines, particles, pressure/velocity color maps).

---

## 5. Extensions and Related Sims

After the basic potential-flow demo works, we can consider:

- **AeroCraft / AeroCombat-style simulators** (from `AeroCraftDesign.h`, `AeroCombat_main.cpp`, `AeroCombatOGL3.cpp`):
  - Use the JS potential-flow results (lift, drag, moments) to drive a simple 6-DOF rigid-body aircraft model.
  - Provide a minimal “aerocraft sandbox” mode in `spacecraft_editor`.

- **2D sail / airfoil games** (from `SailWar_main.cpp`, `AeroSurf2D.h`):
  - Start with tabulated polars (lift/drag vs AoA) in JS.
  - Simple 2D rigid-body plus wind field, later upgraded to 2D potential flow.

The core priority remains: **get a clean, interactive potential-flow demo around wings in JS**, with geometry clearly derived from panels and UV surfaces, and a transparent link back to the C++ reference files listed above.
