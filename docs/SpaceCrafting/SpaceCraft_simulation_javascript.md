# SpaceCraft Simulation – JavaScript Plan (WIP)

This document outlines a **lightweight JavaScript reimplementation** of selected SimpleSimulationEngine demos for easy web publishing and debugging. It is a companion to:

- `SpaceCrafting.md` – legacy spacecraft pipeline overview.
- `SpaceCrafting_new.md` – updated high‑level design and mesh pipeline.
- `SpaceCraftConstructionProblems.md` – detailed geometry/physics problems (telescopic truss, plates, welds, JS prototyping tasks).

The goal is to build small, focused JS/Three.js demos that exercise the same ideas as the C++/SDL tests, but **run directly in the browser**.

---

## 1. Context: JS Spacecraft Editor

All of this should live in **`js/spacecraft_editor`**:

- Interactive mesh/truss playground and prototype for spacecraft‑like structures.
- `MeshesUV` / `MeshBuilder` for tubes, slabs, plates, ropes, stick materials.
- Visualization of vertex/edge selections, SDF tools, and index‑based rails.

Each C++ test we port becomes a small JS “mode” or scene inside this editor, reusing the same mesh and rendering infrastructure.

---

## 2. JS Reimplementation of 3D Physics / Engineering Demos

We use the existing `sketches_SDL/3D` tests as a catalogue (`docs/sketches_SDL_3D.md`, `docs/test_3D_Notes.md`, `tests_bash/sketches_SDL/test_3D.sh`). For web use we focus on **clear, self‑contained cases** that demonstrate either:

- Spatial data structures and mesh topology.
- Radiation / scattering / radiosity ideas.
- Aerodynamics / potential flow relevant to spacecraft.

### 2.1 High‑priority physics / engineering demos to port

From `docs/test_3D_Notes.md` and `docs/sketches_SDL_3D.md`:

- **`test_Radiosity`**  
  - Full radiosity simulation in an orthogonal corridor map.  
  - JS version: start with a **small rectangular room / corridor network** and simple form‑factor approximation or matrix solver.  
  - Target: visualize convergence and color bleeding; later reuse for spacecraft interiors and radiator/shield heating.

- **`test_RayScattererMMC`**  
  - Monte‑Carlo volumetric scattering in tetrahedral/volumetric cells.  
  - JS version: simplified 2D/3D volume with stochastic ray paths, histogram of path length / angular distribution.  
  - Target: intuition for volumetric exhaust plumes, gas clouds, radiation fog around ships.

- **`test_Scatterer`**  
  - Surface‑based scattering over a network of channels connecting surface elements.  
  - JS version: discrete graph of patches with flux transmission matrix; visualize multiple bounces and equilibrium.  
  - Target: approximate **surface‑to‑surface energy/particle transport** for radiators, shields, and complex spacecraft hulls.

- **`test_SphereSampling`**  
  - Procedural planet/asteroid generation, and more importantly **sampling spherical functions** (directional radiation intensity, environment maps).  
  - JS version: focus on **spherical sampling schemes** (icosahedral/hex meshes, equal‑area sampling, directional histograms) and map them to textures / point clouds.  
  - Target: representation of directional radiation fields, sky brightness, threat directions around a spacecraft.

- **`test_VortexLattice`**  
  - Potential‑flow aerodynamics around a wing, using vortex lattice method.  
  - JS version: 2D/3D minimal panel method over a discrete wing or fin, visualizing circulation and induced velocity.  
  - Target: intuitive tool for fin/wing design, reentry control surfaces, maybe radiator aerodynamics.

These five form the **core physics/radiation/aero JS demos**.

### 2.2 Supporting geometry / mesh demos

We also want small JS demos mirroring the **mesh / spatial‑structure tests**:

- **`test_Mesh`**  
  - Mesh editing operations: edge finding, edge collapse, convex hull / from‑planes generators.  
  - JS version: operate on `MeshBuilder` meshes in `spacecraft_editor`, with interactive selection and edge operations.  
  - Target: debug and prototype **mesh editing workflows** for spacecraft hulls, plates, telescopic trusses, and welds.

- **`test_AABBTree`**  
  - AABB tree visualization for many objects.  
  - JS version: build an AABB tree over mesh chunks/triangles and visualize bounding volumes and ray queries.  
  - Target: foundation for **raytracing, picking, and radiation occlusion** in web demos.

- **`test_SphereSampling`** (geometry link)  
  - In addition to planetary noise, explicitly use it as **basis for spherical discretization** of radiation fields – link with `test_Scatterer` and `test_RayScattererMMC`.

These geometry demos concentrate on **topology and data structures** that later support spacecraft radiation and damage models.

### 2.3 Radiation / scattering variants – what each is good for

For radiation and scattering we have several distinct C++ tests; the JS versions should keep their roles clearly separated:

- **`test_Radiosity` – diffuse surface–surface exchange**  
  - **Model:** classic radiosity between diffuse patches using form factors (matrix or iterative solver).  
  - **Good for:** closed or semi‑closed environments (corridors, internal bays, radiator farms) where light/heat bounces many times between walls.  
  - **JS focus:** blocky corridor / room meshes from `MeshBuilder`, convergence visualization, color bleeding, equilibrium temperatures.

- **`test_Scatterer` – network of surface channels**  
  - **Model:** flux moving along a predefined network of channels connecting surface elements (graph‑like transport instead of full matrix over all pairs).  
  - **Good for:** structured surfaces (radiator stacks, shield layers, truss‑mounted panels) where we care about **how power flows through a designed network** rather than arbitrary view‑factors.  
  - **JS focus:** patch graph + flux on edges, compare with radiosity on the same geometry.

- **`test_RayScattererMMC` – volumetric Monte‑Carlo transport**  
  - **Model:** individual rays/particles undergoing random scattering/absorption events inside a volume (Monte‑Carlo random walks).  
  - **Good for:** plumes, gas clouds, reactor exhaust, dust/nebulae – any **volumetric medium** around the spacecraft.  
  - **JS focus:** simple voxel/tetra volume, ray path visualization, path‑length and exit‑angle histograms.

- **`test_SphereSampling` – representation of directional fields**  
  - **Model:** sampling scalar functions on a sphere using structured meshes (icosahedral/hex, octahedral, etc.).  
  - **Good for:** storing and visualizing **directional radiation intensity**, sky brightness, threat directions, sensor coverage.  
  - **JS focus:** sphere discretizations + mapping between direction samples and radiation/scattering outputs from the other demos.

In `js/spacecraft_editor` these tools should be **composable**:

- Geometry from tubes/slabs/plates → used by AABBTree, radiosity, and surface scatterer.
- Directional results from radiosity/scattering → projected onto sphere samplings from `test_SphereSampling`.
---

## 4. Relation to SpaceCrafting Pipeline Docs

### 4.1 Link to `SpaceCrafting.md` and `SpaceCrafting_new.md`

- Those docs describe the **C++ pipeline**:
  - Lua → `SpaceCraft` components → `BuildCraft_blocks` / `Mesh::Builder2` → truss / radiosity / scattering solvers.
- The JS demos are **front‑end playgrounds** that:
  - Use analogous mesh structures (`MeshBuilder`, UV‑based generators, stick materials).  
  - Prototype **topology patterns** (tubes, slabs, plates, sliders, telescopic trusses) before hardening them in C++.
- Many ideas in `SpaceCrafting_new.md` §5–§9 (sliders, tubes, plates, welds, radiation) already have JS counterparts in `spacecraft_editor`; this document makes explicit which **C++ tests** we want JS versions of.

### 4.2 Link to `SpaceCraftConstructionProblems.md`

`SpaceCraftConstructionProblems.md` defines detailed geometry/physics problems and **JS‑first tasks** (J0–J6) that are mostly implemented in `js/spacecraft_editor`:

- Telescopic truss recoil damper (`SlabTube` / `TubeSheet` tubes, sliders on polygonal rails).  
- Plates on girders and ropes (radiators, sails, shields) via parametric quad/triangulation helpers.  
- Welds / automatic connections using SDF‑based selections and proximity bonds.  
- JS visualizations for SlabTube, TubeSheet, rail selection strategies (index‑based, angle‑based, SDF‑based).

The **new JS demos listed here** (radiosity, scatterers, sphere sampling, AABB tree, mesh editing) should be seen as the **simulation‑side companions** to those geometry problems:

- JS geometry tools create the meshes (tubes, plates, hulls, telescopic trusses).  
- JS physics/radiation demos operate on those meshes (flux transport, visibility, occlusion, heating).

---

## 5. Next Steps (Implementation Sketch)

Short‑term JS tasks:

- Implement minimal **radiosity demo** on a blocky corridor mesh generated in JS.  
- Implement **surface scattering** and **volumetric scattering** prototypes with simple UIs for source placement and material parameters.  
- Add **sphere‑sampling utilities** in JS (icosahedral grids, HEALPix‑like approximations) and link them to scattering/radiosity as input/output spaces.  
- Add **AABB tree + ray query visualizer** attached to MeshBuilder meshes.  
- Extend the **mesh editing panel** in `spacecraft_editor` to explicitly mirror selected `test_Mesh` operations.

Once stable, these patterns can be ported or mirrored in C++ (`TriangleRayTracer`, `Radiosity`, `Scatterer`, `Scatterer2`) so that both **desktop** and **web** paths use conceptually the same building blocks.
