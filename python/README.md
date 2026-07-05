Python simulation modules — PyOpenCL GPU solvers, ctypes C++ bindings, scientific computing, and visualization tools for physics, chemistry, aerospace, and graphics.

## Subfolders with README

- **EulerianImpacFluid/** — 2D Eulerian multi-material compressible flow solver with level-set interface tracking (PyOpenCL). Stiffened-gas EOS for high-velocity impact simulations (uranium into liquid hydrogen). Has [README.md](EulerianImpacFluid/README.md).
- **pyTruss/** — Truss/cloth simulation using Projective Dynamics with iterative linear solvers (Jacobi, Gauss-Seidel, VBD, Chebyshev) on CPU and GPU (OpenCL). Has [README.md](pyTruss/README.md).

## Subfolders without README

- **GLCL/** — Monolithic PyQt5+OpenGL+OpenCL N-body demo. Superseded by GLCL2.
- **GLCL2/** — Modular OpenGL+OpenCL simulation framework: separate `OGLsystem.py` (GL), `OCLsystem.py` (CL), `GLCLGUI.py` (widget), `GLCLBrowser.py` (script loader). Successor to GLCL.
- **GUI/** — `GUITemplate.py` — reusable PyQt5 GUI template for Python apps.
- **LandCraft/** — `landcraft.py` — 2D terrain/landscape generation and simulation.
- **OrbitalWar/** — OpenCL-based orbital trajectory optimizer + space combat tactics model. Includes `.cl` kernel and detailed `.md` derivation.
- **Radiosity/** — 2D/3D radiative heat transfer solver: surface discretization, view factors, radiosity equation solver, polygon rasterization, octahedral sphere mapping, triangle occlusion raytracer. Has [README.md](Radiosity/README.md).
- **benchPython/** — Python performance benchmarking scripts.
- **particleCollisions2D/** — 2D particle system with Morse+Coulomb forces, Projective Dynamics constraints, sparse iterative solvers. Three versions (class-based, O(1), functional). Known bug in `particle_system_new.py:11`.
- **patterns/** — `patterns.py` — design pattern implementations and examples.
- **pyApprox/** — Numerical approximation methods: integration, radial basis, rational, Taylor series.
- **pyCombatModels/** — Combat models: air, land, space engagement math (Lanchester-type equations).
- **pyFlight/** — ctypes bindings to C++ flight simulation library (`libFlight.so`).
- **pyMeta/** — Python metaprogramming utilities.
- **pyMetaC/** — C code generation via Python metaprogramming.
- **pyMolecular/** — Molecular simulation toolkit: Gaussian orbitals, Boys function, eFF (electron force field), CLCFGO, MMFF, reactive force fields, rigid molecules, OpenCL MD. Largest Python module (39 items).
- **pyPapers/** — Bibtex/Mendeley bibliography utilities.
- **pyRay/** — Modular SDF ray tracer: scene, image, OpenCL acceleration, common utilities. Has [README.md](pyRay/README.md).
- **pyScatter/** — X-ray / neutron radiation scattering simulation: spacecraft geometry, attenuation, Monte Carlo design doc, OpenCL kernels, radiosity pipeline. Has [README.md](pyScatter/README.md).
- **pyShaderToy/** — ShaderToy-like GPU shader explorer with formula nodes, expression parser, GLSL code generation, 43+ shader files.
- **pyShock/** — ctypes bindings for C++ shock physics solver (incomplete — no `__init__.py`).
- **pySimE/** — Simulation engine package: constants, physics, chemistry, and `space/` (92 items: rockets, railguns, orbital transfer optimization, space elevators, sky hooks, asteroid analysis).
- **pySymGLSL/** — Symbolic GLSL simulation framework: importable package with pipelines, shaders, GUI widget. More structured than pyShaderToy.
- **pyVis3D/** — ctypes bindings for 3D visualization C++ library.
- **pyplot3D/** — Simple 3D plotting utilities.
- **subPixelContour/** — Sub-pixel contour extraction algorithm (42K).
- **terrain_ocl/** — Terrain erosion simulation with OpenCL (680-line `terrain.py`, 40K kernel). Has docs on watershed transform and accumulative SDF erosion.
- **terrain_ocl_2/** — Refactored terrain erosion (416-line `terrain.py`, 3K kernel). Successor to terrain_ocl.
- **utils/** — C++ header extraction utilities.

## Cross-Module: Radiosity ↔ pyScatter ↔ pyRay

These three modules share the geometric problem of **ray-triangle intersection for occlusion/visibility**, but differ in physics and output:

| | **Radiosity** | **pyScatter** | **pyRay** |
|---|---|---|---|
| **Physics** | Thermal radiosity (diffuse IR, Stefan-Boltzmann) | X-ray/neutron attenuation + scattering (Beer-Lambert) | Rendering (SDF sphere tracing) |
| **Geometry** | Explicit triangles (mesh) | Explicit triangles + tubes | Implicit SDF primitives |
| **Core query** | "Can element i see j?" (binary occlusion) | "What's the attenuation along this ray?" (continuous μ·L) | "What's the first hit?" (closest surface) |
| **Occlusion** | N² ray-triangle, binary mask (legacy) / sparse top-K channels (new) | M·K ray-primitive, accumulated μ·L | Implicit in ray marching |
| **GPU kernel** | `radiosity.cl::occlusion_matrix` | `SpacecraftSimulation.cl::compute_sparse_channels` + `Xray.cl::xray_simulation` | `rayScene_basic.cl::rayTrace_basic` |
| **Acceleration** | None (brute force) / sparse top-K + cluster culling (new) | Spatial clustering + local memory + cone culling | SDF sphere tracing |
| **Output** | Temperature distribution | Detector image / attenuation map | Rendered image |

**Shared bottleneck**: ray-triangle intersection between points in 3D space. Radiosity needs binary N² visibility; pyScatter needs accumulated attenuation along rays. The intersection primitive (Möller-Trumbore) is duplicated across `radiosity.cl` and `SpacecraftSimulation.cl` — should be unified into a single shared kernel.

**Already shared**: `SimulationPipeline.py` in pyScatter imports and reuses Radiosity's `rasterize_face_to_elems`, `assemble_elements`. Sparse radiosity solve (`solve_radiosity_sparse`, `calculate_temperatures_sparse`) is implemented locally in `SimulationPipeline.py`.

**pyRay is fundamentally different** — it uses SDF ray marching, not explicit triangle intersection. Not directly applicable to the radiosity/scattering occlusion problem.

**Needed**: unified ray-triangle intersection primitive. BVH acceleration still identified as future work for candidate scan phase.

## End-to-End Status (tested on `ship_ICF_marksman_2.obj`)

Real spacecraft geometry: 1271 vertices, 4268 edges, **zero faces** (wireframe/truss, ~200×400×1100 m).

| Capability | Status | Full ship? | Notes |
|---|---|---|---|
| X-ray attenuation (tubes) | ✅ Complete | ✅ Yes | 4268 tubes, 67 clusters, 0.006s on Intel GPU. Scales as M×K. |
| X-ray attenuation (triangles) | ✅ Complete | N/A | No faces in OBJ; tested on `high_res_mesh.obj` (84 tris) |
| Radiosity (view factors + solve) | ✅ Complete | ✅ Yes | NumPy vectorized; scales as N² in memory but not GPU-limited |
| Radiosity occlusion (GPU) | ✅ Sparse top-K channels | ⚠️ Subset tested | `compute_sparse_channels`: O(N·K·C), row-chunked. Tested: 120 elems (T⁴=68.59 vs 68.82 dense), 2000 elems (T⁴=63.13, 16.7s). Full-ship not yet tested. |
| Edge→surface conversion | ⚠️ Ad hoc | — | Inline in test scripts; need reusable utility |
| Morton/Z-order sorting | ❌ Missing | — | Current cluster sort is simple z+x; Morton code would improve locality |
| GPU sparse solver | ❌ Missing | — | Current `solve_radiosity_sparse` is CPU gather; should move to GPU for large N |
| Monte Carlo scattering | ❌ Skeleton | — | `Scatter.cl` incomplete; design doc only |
| pyRay integration | ❌ N/A | — | SDF ray marching, different paradigm |

**Bottom line**: X-ray/attenuation path is complete and works end-to-end on the full ship. Radiosity now uses sparse top-K channels (no dense N×N matrix) — tested on subsets up to 2000 elements with close parity to dense method. Full-ship run pending. Monte Carlo scattering is not implemented.

## Root-level files

- **BaseGUI.py** — Base PyQt5 GUI class for simulation apps.
- **PGS.py** — Projected Gauss-Seidel constraint solver.
- **constrain_solver.py** — General constraint solver implementation.
- **try_CG.py** — Conjugate gradient method experiments.
- **try_Linearized_Move.py** — Linearized movement solver experiment.
- **test_*.py** (17 files) — Test scripts for modules across pyMolecular, pySimE, pyFlight, pyVis3D, chemistry, physics.
