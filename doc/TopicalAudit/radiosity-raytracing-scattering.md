---
type: TopicalAudit
title: Radiosity, Ray Tracing & Radiation Scattering
tags: [topic, cross-language, radiosity, ray-tracing, scattering, radiation, gpu]
---

## Summary

Multiple approaches to radiative heat transport, radiosity, and particle scattering (X-ray, neutron). C++ implementations: `TriangleRayTracer.h` (occlusion), `Radiosity.h` (matrix-based radiosity), `Scatterer.h` (matrix-based scattering), `Scatterer2.h` (channel-network scattering). Python: 2D/3D radiosity solvers in `python/Radiosity/`, OpenCL X-ray attenuation in `python/pyScatter/`. Shared OpenCL kernel `radiosity.cl` for GPU occlusion. Primary application: spacecraft design — radiation from engines, weapons, stars; structural thermal loading.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/common/dynamics/TriangleRayTracer.h` | active | Triangle sampling to `SurfElement`s, sequential occlusion test via `rayInTriangle()`. No spatial acceleration. Base class for Radiosity and Scatterer. |
| C++ | `cpp/common/dynamics/Radiosity.h` | active | Dense coupling matrix `M[i][j]` with Lambertian cosine/r² + occlusion. `step_Direct` iterative update. `processTriangles` TODO. Toy/demo quality. |
| C++ | `cpp/common/dynamics/Scatterer.h` | unfinished | Matrix-based scattering with direction bins, `ScatterElem`. Many TODOs, unclear normalization. Largely abandoned. |
| C++ | `cpp/common/dynamics/RayScatter.h` | active | Monte Carlo ray scattering: `ScatterMaterial` (mean free path, angular distribution), `RayPath` (path probability), `RayScatterWorld` (multi-mesh tracing with stochastic scattering). For neutron transport, participating media. |
| C++ | `cpp/common/dynamics/Scatterer2.h` | unfinished | Channel-network approach: `Channel`/`HalfChannel` between `ScatterElem2`s. `step_Direct` has scattering + transfer phases. **Bug**: uses undefined `n` (should be `scattelems.size()`), `backChans` not populated. |
| C++ | `cpp/sketches_SDL/3D/test_Radiosity.cpp` | active | Demo: load mesh, sample elements, build coupling, iterate, visualize panels + coupling lines |
| C++ | `cpp/sketches_SDL/3D/test_Scatterer.cpp` | active | Demo for Scattering2: two triangles, flux sources, channels, visualization |
| C++ | `cpp/sketches_SDL/3D/test_RayScattererMMC.cpp` | experimental | Ray-based Monte Carlo scatterer test |
| OpenCL | `cpp/common_resources/cl/radiosity.cl` | active | GPU occlusion kernel: brute-force + broad-phase (group-level) + narrow-phase design. Local memory tiling, barrier sync. Used by Python `TriangleOcclusionRaytracer.py`. |
| Python | `python/Radiosity/Radiosity.py` | active | 2D radiosity: line elements, view factors, solve (I-W)B=P, T⁴ formulation. No occlusion. |
| Python | `python/Radiosity/Radiosity3D.py` | active | 3D radiosity: hex-rasterize mesh faces to elements, 3D view factors, solve, visualize. No occlusion. |
| Python | `python/Radiosity/TriangleOcclusionRaytracer.py` | active | PyOpenCL occlusion matrix between element COGs using `radiosity.cl`. Separate from Radiosity3D — integration TODO. |
| Python | `python/Radiosity/PolygonRasterization.py` | active | Hex rasterization of 2D polygons: clipping, centroid/area, boundary merging. Used by Radiosity3D. |
| Python | `python/Radiosity/OctahedralSphereMaping.py` | active | Octahedral sphere parameterization for direction sampling |
| Python | `python/pyScatter/xray_sim.py` | active | OpenCL X-ray attenuation: Beer-Lambert through clustered tubes/triangles. Deterministic, no MC scatter. |
| Python | `python/pyScatter/cl/SpacecraftSimulation.cl` | active | OpenCL kernel for spacecraft radiation simulation |
| Python | `python/pyScatter/cl/cl/Scatter.cl` | unfinished | Placeholder kernel — no real scattering logic, needs full rewrite |
| Python | `python/pyScatter/AttenuationTest.py` | active | Attenuation test driver |
| Python | `python/pyScatter/SimulationPipeline.py` | active | Pipeline orchestration for pyScatter |
| Python | `python/pyScatter/SpacecraftGeometry.py` | active | Spacecraft geometry generation for scattering tests |
| Doc | `doc/RayTracing.md` | active | CPU vs OpenCL implementation comparison for radiosity coupling |
| Doc | `docs/RadiosityRaytracing/ScatteringSolverComparison.md` | active | Detailed comparison: `Scatterer.h` (matrix) vs `Scatterer2.h` (channel) |
| Doc | `docs/RadiosityRaytracing/RadiosityAndScatteringOverviewDiscussion.md` | active | LLM chat: full review of all scattering/radiosity files, status, gaps |
| Doc | `docs/RadiosityRaytracing/TriangleOcclusionRaytracer.md` | active | Design plan: GPU occlusion, broad-phase/narrow-phase hierarchy |
| Doc | `docs/RadiosityRaytracing/SphereSampling.md` | active | Icosahedral diamond mapping for S² sampling (used for direction discretization) |
| Doc | `python/Radiosity/Radiosity.md` | active | 2D radiosity tutorial/design |
| Doc | `python/Radiosity/Radiosity3D.md` | active | 3D radiosity design note |
| Doc | `python/Radiosity/OctahedralSphereMaping.md` | active | Octahedral mapping design |
| Doc | `python/Radiosity/PolygonRasterization.md` | active | Hex rasterization design |
| Doc | `python/pyScatter/MonteCarloScattere.md` | active | MC scattering design dialogue (ideas only, no code) |

## Parity Status

- **C++ `Radiosity.h` ↔ Python `Radiosity3D.py`**: Same algorithm (coupling matrix + iterative solve). C++ has occlusion via `TriangleRayTracer`, Python does not (separate `TriangleOcclusionRaytracer.py` not integrated).
- **C++ `Scatterer.h` ↔ `Scatterer2.h`**: Different approaches to same problem. `Scatterer2.h` is newer but buggy. `Scatterer.h` is abandoned.
- **OpenCL `radiosity.cl` ↔ C++ `TriangleRayTracer.h`**: GPU vs CPU occlusion. GPU uses local memory tiling + broad/narrow phase. CPU is brute-force sequential. No automated parity test.
- **Python `xray_sim.py` ↔ C++**: Python-only, no C++ counterpart for X-ray attenuation.

## Open Issues

- `Scatterer2.h` `step_Direct` uses undefined `n` — should be `scattelems.size()`
- `Scatterer2.h` `backChans` allocated but never populated — backward channel mapping missing
- `Scatterer.h` largely unfinished, many TODOs — candidate for deprecation
- `Scatter.cl` is a non-functional placeholder
- Python `Radiosity3D.py` lacks occlusion — `TriangleOcclusionRaytracer.py` exists but is separate
- `Radiosity.h` `processTriangles` is commented out (TODO)
- `Xray.cl` assumes `cluster_size ≤ 64` (CACHE_CAPACITY), no caps on tubes
- No cross-language parity tests for any radiosity/scattering result
- Spacecraft application: radiation from engines/stars/weapons uses radiosity+scattering tools with different source terms (see `encyclopedia/space_warfare/07_Physical_model_and_equations.md`)

## Algorithm Review: Radiation Transport Methods

This section reviews algorithms for solving radiation transport — both thermal (visible/IR) and ionizing (neutron/gamma/X-ray) — in two regimes relevant to spacecraft simulation: **dense environments** (reactor interiors, shielding layers) and **sparse environments** (girders, thin sheets, mostly empty space between components). For each method we note merits, use cases, strengths, weaknesses, GPU suitability, and implementation status within this codebase.

### A. Radiosity (Matrix-Based View Factor Method)

**Merit**: Solves equilibrium thermal radiation exchange between diffuse (Lambertian) surfaces. Builds a coupling matrix $F_{ij}$ (view factors) and solves $(I - \rho F) B = P$ for outgoing radiosity $B$ given emissive power $P$ and reflectivity $\rho$.

**Use cases**: Thermal equilibrium of spacecraft surfaces (radiators, hull, engine heat), laser heating distribution, IR signature estimation. Best when surfaces are diffuse and angular distribution is smooth (cosine law suffices — no per-element angular storage needed).

**Strengths**:
- Exact for diffuse surfaces; no angular discretization needed
- Linear system is sparse-ish (only coupled by visibility); iterative solvers (Jacobi, CG) converge fast
- Deterministic — no noise, reproducible
- GPU-friendly occlusion pass (brute-force or tiled)

**Weaknesses**:
- O(N²) coupling matrix for N elements — memory bottleneck for large scenes
- Cannot handle specular reflection, transmission, or anisotropic scattering
- No angular resolution — assumes diffuse emission/reflection
- Requires occlusion/visibility computation as a separate pass

**GPU parallelization**: Occlusion matrix can be computed on GPU (brute-force or tiled local memory). Matrix assembly and solve are CPU-side currently. Could parallelize matrix-vector products on GPU.

**Implementation status**:
- **Python 2D**: `python/Radiosity/Radiosity.py` — active, complete (line elements, view factors, T⁴ formulation, solve). No occlusion.
- **Python 3D**: `python/Radiosity/Radiosity3D.py` — active, complete (hex-rasterized mesh faces, 3D view factors, solve, visualize). No occlusion integrated.
- **C++**: `cpp/common/dynamics/Radiosity.h` — active but toy/demo quality. Dense coupling matrix with optional occlusion via `TriangleRayTracer`. `processTriangles` is TODO.
- **GPU occlusion**: `cpp/common_resources/cl/radiosity.cl` — active. Brute-force + broad-phase design. Used by `TriangleOcclusionRaytracer.py` but not integrated into `Radiosity3D.py`.
- **Gap**: Occlusion not wired into Python 3D radiosity pipeline. No GPU-accelerated solve. No spectral/energy-dependent view factors.

### B. Beer-Lambert Attenuation (Deterministic Ray Marching)

**Merit**: Computes transmission $T = e^{-\sum \mu_i \cdot d_i}$ along a straight ray through thin materials. Each intersected primitive (triangle sheet or cylinder) contributes $\mu \cdot \text{thickness}$ to optical depth.

**Use cases**: X-ray/gamma shadow imaging, first-order shielding estimation, detector flux maps. Ideal for sparse spacecraft geometry where most rays pass through few thin objects.

**Strengths**:
- Extremely simple and fast — just intersection + accumulation
- Deterministic, no noise
- GPU-friendly with cluster tiling and cone culling
- Handles projected thickness correctly (tilted sheet → longer path → more attenuation)

**Weaknesses**:
- No scattering — only straight-line attenuation
- No energy degradation (monoenergetic only)
- No secondary radiation (fluorescence, bremsstrahlung)
- Cluster size limited by local memory (CACHE_CAPACITY=64 in current kernel)

**GPU parallelization**: Excellent. Workgroup = detector tile. Cluster tiling into local memory with barrier sync. Cone-sphere culling to skip irrelevant clusters. All implemented in `xray_sim.py` + `SpacecraftSimulation.cl`.

**Implementation status**:
- **Python+OpenCL**: `python/pyScatter/xray_sim.py` + `cl/SpacecraftSimulation.cl` — active. Beer-Lambert through clustered tubes/triangles with cone culling. Deterministic, no MC scatter.
- **C++**: No direct counterpart. `RayScatter.h` has Beer-Lambert as part of ray-path evaluation but integrated with scattering.
- **Gap**: No multi-energy attenuation (spectrum hardening). No scattered flux contribution. Cluster size hardcoded to 64.

### C. Monte Carlo Path Tracing with Next Event Estimation (NEE)

**Merit**: Stochastic sampling of scattering paths. For single scattering: pick a scatterer object, pick a point on it, compute weight = $L \cdot T_{in} \cdot \sigma_s \cdot G \cdot T_{out} / \text{PDF}$. NEE directly connects source → scatter point → detector, avoiding wasted rays into empty space.

**Use cases**: Ionizing radiation scattering (neutron, gamma) in sparse spacecraft geometry. Single-scatter NEE is efficient when scattering is rare (thin materials, mostly empty space). Multi-scatter extensions needed for dense shielding.

**Strengths**:
- Unbiased — converges to exact answer with enough samples
- Handles arbitrary geometry (triangles, cylinders, any primitive with area sampling)
- Naturally parallel: 1 thread = 1 sample path, no inter-thread dependencies
- Importance sampling can dramatically reduce variance (sample by area, by expected contribution)
- Extensible to multi-bounce, energy degradation, anisotropic scattering

**Weaknesses**:
- Noisy at low sample counts — need many samples for smooth results
- Transmission computation (attenuation along source→scatter and scatter→detector) is the bottleneck — requires ray-scene intersection against all objects
- Variance increases with scene complexity (many blockers → most paths have T≈0)
- Multi-scatter paths exponentially harder to sample efficiently

**GPU parallelization**: Excellent for single scattering. Each thread = 1 MC sample. Attenuation via tile-based cone-walking (cooperative cluster loading into local memory). Workgroup restricted to small spatial cluster → all rays form narrow cone → only relevant blocker clusters loaded. Designed in detail in `MonteCarloScattere.md`.

**Implementation status**:
- **Design only**: `python/pyScatter/MonteCarloScattere.md` — extensive design dialogue (1742 lines). Covers NEE math, cone-sphere culling, tile-based local memory loading, GPU kernel pseudocode, data structures (AoS Primitive, Cluster), multi-scatter extension ideas. **No kernel implemented.**
- **C++ CPU MC**: `cpp/common/dynamics/RayScatter.h` + `test_RayScattererMMC.cpp` — experimental. CPU ray-based MC with `ScatterMatrial` (mean free path, angular CDF with `nScatAng=16` bins, Lorenzian distribution). `RayScatterWorld` traces rays through triangle meshes with stochastic scattering. `RayPath::eval` computes path weight with Beer-Lambert + angular scattering probability. **Works but CPU-only, no GPU.**
- **Gap**: No GPU MC scattering kernel. `Scatter.cl` is a non-functional placeholder. The NEE design in `MonteCarloScattere.md` is the most promising path forward.

### D. Multi-Scatter / Bidirectional Path Tracing / Metropolis Light Transport

**Merit**: Extends single-scatter MC to multiple scattering events. Bidirectional path tracing builds paths from both source and detector, connecting in the middle. Metropolis Light Transport (MLT) uses Markov Chain MC to efficiently sample hard-to-reach paths (e.g., detector hidden behind many layers).

**Use cases**: Deep shielding penetration, radiation through complex labyrinths, multi-bounce neutron/gamma scattering in dense materials. Not needed for sparse spacecraft single-scatter but essential for reactor interiors.

**Strengths**:
- Bidirectional: much more efficient than unidirectional for hidden detectors
- Metropolis: excels when most paths have zero contribution (caustics, deep penetration)
- Physically exact — handles any scattering order

**Weaknesses**:
- High implementation complexity
- Metropolis has initial bias (burn-in period), chain correlation
- Hard to parallelize efficiently (MLT paths are correlated)
- Overkill for sparse geometry where single-scatter NEE suffices

**GPU parallelization**: Bidirectional path tracing is parallelizable (independent path pairs). Metropolis is harder — typically uses multiple independent chains. For spacecraft application, the design doc explicitly recommends **against** Metropolis for single scattering and suggests standard importance sampling instead.

**Implementation status**:
- **Not implemented**. Mentioned as future direction in `MonteCarloScattere.md`: "Later we will consider multiple-scattering, something like bidirectional path integration used e.g. for Metropolis light transport."
- **Gap**: No multi-scatter GPU kernel. No bidirectional path tracing. No Metropolis. These are long-term goals for dense shielding simulation.

### E. Neutron Diffusion (Multi-Group, Finite Volume)

**Merit**: Solves the diffusion approximation of the neutron transport equation: $-\nabla \cdot (D \nabla \phi) + \Sigma_a \phi = S$, where $D = 1/(3\Sigma_{tr})$ and $\Sigma_{tr} = \Sigma_a + \Sigma_s(1-g)$. Multi-group energy discretization with downscattering. Solved as a tridiagonal system per energy group.

**Use cases**: Reactor physics — neutron flux distribution in fissile assemblies, criticality ($k_{eff}$), burnup. Ideal for **dense, optically thick** materials where diffusion approximation is valid (many scatterings isotropize the flux). 1D spherical geometry for implosion devices.

**Strengths**:
- Extremely efficient: tridiagonal solve is O(N) per energy group
- Naturally conserves particles (finite volume)
- Handles multi-energy groups with downscattering
- No angular discretization needed (diffusion smears angles)
- Couples naturally with hydrodynamics (moving mesh, compression)

**Weaknesses**:
- **Invalid in optically thin / sparse media** — diffusion assumes many scatterings per mean free path
- **Invalid near boundaries** — requires Marshak/vacuum boundary corrections
- No angular information — cannot model directional beams or streaming
- Cannot capture strong anisotropy (forward-peaked scattering)
- 1D spherical only in current implementation

**GPU parallelization**: Tridiagonal solve is sequential per group but embarrassingly parallel across energy groups. Could use GPU batched tridiagonal solvers. Currently CPU-only (SciPy `solve_banded`).

**Implementation status**:
- **Python**: `doc/python/Burn1D/neutron_diffusion_solver.py` — active. Multi-group (G=8), 1D spherical, moving mesh with implosion, fission + burnup + poison products. Uses `scipy.linalg.solve_banded`. Coupled to `implosion_solver.py` for hydrodynamics.
- **Design doc**: `doc/python/Burn1D/RadiationTransferSpherical.md` — extensive design dialogue covering diffusion theory, finite volume discretization, boundary conditions, and Python implementation.
- **Gap**: 1D only. No 2D/3D diffusion. No GPU. No coupling with spacecraft scattering code. Diffusion approximation invalid for sparse spacecraft geometry.

### F. Discrete Ordinates (S_N) / Channel-Based Propagation

**Merit**: Discretizes angular space into discrete directions (ordinates) and solves the transport equation along each direction. The `Scatterer2.h` channel-network approach is a variant: explicit channels between scattering elements carry flux in specific directions, with scattering events redistributing flux between channels.

**Use cases**: Radiation transport in intermediate optical thickness where diffusion is invalid but MC is too noisy. Can handle anisotropic scattering. Suitable for shielding calculations with known angular distributions.

**Strengths**:
- Deterministic — no noise
- Captures angular dependence (unlike diffusion)
- Can handle forward-peaked scattering
- More efficient than MC for optically thick, smooth problems

**Weaknesses**:
- Angular discretization introduces ray effects (artifacts along discrete directions)
- O(N²) channel/matrix storage for N elements × M directions
- Convergence requires many ordinates for high angular resolution
- Complex to implement correctly (channel bookkeeping, backward mapping)

**GPU parallelization**: Channel-based propagation is parallelizable across elements. Matrix-based approach can use GPU sparse matrix operations. Not currently GPU-implemented.

**Implementation status**:
- **C++ matrix-based**: `cpp/common/dynamics/Scatterer.h` — unfinished, largely abandoned. Dense coupling matrix with direction bins (`nScatAng=16`), `ScatterElem` with ray arrays. Many TODOs, unclear normalization.
- **C++ channel-based**: `cpp/common/dynamics/Scatterer2.h` — unfinished. `Channel`/`HalfChannel` between `ScatterElem2`s. `step_Direct` has scatter + transfer phases. **Bug**: uses undefined `n` (should be `scattelems.size()`), `backChans` not populated. `ScatterElem2` has anisotropic kernel (`rot`, `thicks`, `areas`, `beta`).
- **Comparison doc**: `docs/RadiosityRaytracing/ScatteringSolverComparison.md` — detailed comparison of the two approaches.
- **Gap**: Both C++ implementations are broken/unfinished. No GPU version. No Python port. Channel bookkeeping (back-channels, indexing) needs completion.

### G. P_N (Spherical Harmonics) Method

**Merit**: Expands angular dependence in spherical harmonics $Y_l^m(\Omega)$, truncating at order N. Converts the transport equation into a coupled system of spatial equations for expansion coefficients. Eliminates ray effects of S_N methods.

**Use cases**: Dense media with smooth angular distributions. Reactor core calculations. Good where angular flux is smooth (many scatterings).

**Strengths**:
- No ray effects (unlike S_N)
- Galerkin convergence — spectral accuracy in angle
- Natural for problems with near-isotropic flux

**Weaknesses**:
- High order needed for forward-peaked or streaming problems
- Coupled system is larger and denser than diffusion
- Complex implementation (recurrence relations, boundary conditions)
- Not suitable for sparse media with directional beams

**GPU parallelization**: Matrix assembly and solve can be parallelized. Not currently implemented.

**Implementation status**:
- **Not implemented**. Mentioned in `OctahedralSphereMaping.md` as a comparison: "Spherical Harmonics can be computationally expensive for high-frequency data" and "prone to ringing artifacts."
- **Gap**: No P_N implementation. Not a priority for sparse spacecraft geometry where MC/NEE is more appropriate.

### H. Hybrid / Clumped Attenuation Heuristics

**Merit**: Approximate attenuation through clustered geometry using precomputed "optical tensors" per cluster — mass density + projected-area tensor (shadow ellipsoid). Avoids per-primitive ray intersection. Wigner-style blend between geometric shadow and mass attenuation.

**Use cases**: Fast first-order shielding estimates for spacecraft. When exact per-element attenuation is too expensive but pure Beer-Lambert through individual primitives is too slow for large scenes.

**Strengths**:
- Very fast — precompute per-cluster optical properties, then just evaluate tensor along ray
- Handles anisotropic clustering (shadow ellipsoid captures orientation)
- GPU-friendly — just a few dot products per cluster per ray

**Weaknesses**:
- Approximate — loses fine geometric detail
- Requires precomputation of cluster optical tensors
- Validation needed against exact Beer-Lambert

**Implementation status**:
- **Design only**: Discussed in `MonteCarloScattere.md` as "clumped-attenuation" mode. Precompute per-cluster: `mu_iso` (mass attenuation), `sigma_geo` (projected-area tensor), mean chunk size. Mode switch between exact geometric and clumped.
- **Gap**: Not implemented. Needs precomputation pipeline and kernel mode switch.

## Angular Sampling Methods

Representing functions on the sphere (directional distributions, scattering phase functions, environment maps) is essential for radiation transport. Several approaches are implemented or designed:

### 1. Icosahedral Diamond Mapping

**What**: Seamless parameterization of S² via 10 rhombi (diamonds) derived from the icosahedron. Each rhombus = 2 triangles. Barycentric coordinates within spherical triangles, with edge-flip to ensure continuity across seams.

**Where**: `cpp/common/maps/SphereSampling.h` — `sampleIcosa2quads()`, `diTri2cartes()`, `icosa2cartes()`, `getIcosaFace()`. Rendering: `cpp/common_SDL/SDL2OGL/DrawSphereMap.h` — `drawIcosaMap()`, `drawDiTri_seam()`. Demo: `cpp/sketches_SDL/3D/test_SphereSampling.cpp`. Doc: `docs/RadiosityRaytracing/SphereSampling.md`.

**Status**: Complete for forward mapping and rendering with continuous seams. `evalIcosaGrid()` and `accumIcosaGrid()` (inverse projection for accumulation) **not yet implemented** — needed for storing/querying directional distributions in scattering solvers.

**Pros**: Seamless (no cracks), consistent edge orientation, good area uniformity (20 triangles ≈ equal area), bijectable.
**Cons**: Complex implementation (plane tests, barycentric, diamond flips). 10 rhombi × N² grid = 10N² total samples.

### 2. Octahedral Mapping

**What**: Projects S² onto octahedron surface (L1 normalization: $p' = p / \|p\|_1$), then unfolds 8 triangular faces into a square UV map. Upper hemisphere → central diamond, lower hemisphere → 4 corners. Barycentric interpolation with alternating diagonal split for continuity.

**Where**: `python/Radiosity/OctahedralSphereMaping.py` — complete Python implementation with `map_3d_to_uv()`, `find_interpolation_data()`, barycentric interpolation, visualization. Doc: `python/Radiosity/OctahedralSphereMaping.md` — detailed mathematical derivation. C++: `SphereSampling.h::sampleOctahedron()` — minimal helper (octant ID + |x|,|y| only, no full mapping).

**Status**: Python implementation complete with interpolation and visualization. C++ minimal/incomplete. Octahedral mode in `test_SphereSampling.cpp` uses `octa_decode_dir()` and `drawSphere_oct()`.

**Pros**: Simple projection (L1 norm), single square texture (easy GPU texture lookup), good for environment maps and direction encoding. Alternating diagonal ensures continuity.
**Cons**: Non-uniform sampling (area distortion near corners), more samples wasted near octahedron vertices. Not seamless in C++ yet.

### 3. Fibonacci / Golden Ratio Sphere Sampling

**What**: Generates quasi-uniform points on S² using the golden angle: $\phi_i = i \cdot \pi(3 - \sqrt{5})$, $z_i = 1 - 2i/n$, $r_i = \sqrt{1 - z_i^2}$. Produces a sunflower/spiral pattern on the sphere with excellent area uniformity.

**Where**: `cpp/sketches_GLView/invSphereMap.cpp` — `SphereSamplesFibonachi(int n, Vec3d* ps)`. Also includes inverse spherical Fibonacci mapping (commented out, from Keinert et al. 2015). References: extremelearning.com.au, arxiv.org/abs/0912.4540.

**Status**: Forward sampling implemented (generate N points). Inverse mapping (direction → nearest sample index) is commented out / referenced but not active.

**Pros**: Extremely simple (3 lines of code), excellent uniformity, no poles, deterministic. Good for MC direction generation and quasi-Monte Carlo.
**Cons**: No hierarchical structure (can't easily refine locally). Inverse lookup is non-trivial (requires Fibonacci lattice search). Not a grid — can't do barycentric interpolation between neighbors easily.

### 4. Random Sphere Sampling

**What**: Uniform random sampling: $z \sim U(-1,1)$, $\phi \sim U(0,2\pi)$, then $(r\cos\phi, r\sin\phi, z)$ where $r = \sqrt{1-z^2}$.

**Where**: `cpp/common/math/Vec3.h:446-460` — `setHomogenousSphericalSample(u, v)` and `fromRandomSphereSample()`.

**Status**: Complete, used throughout for random direction generation (crater positions in `test_SphereSampling.cpp`, etc.).

**Pros**: Simplest possible, unbiased, uniform area distribution.
**Cons**: No structure, noisy (MC variance), no interpolation possible.

### 5. Angular CDF Bins (Tabulated Phase Function)

**What**: Discretize scattering angle $\theta$ into N bins. Store cumulative distribution function (CDF) and density (PDF) per bin. Sample by inverting CDF: pick random $r \in [0,1]$, find bin via linear/binary search, interpolate within bin.

**Where**: `cpp/common/dynamics/RayScatter.h` — `ScatterMatrial` struct with `nScatAng=16` bins, `Sang[]` (CDF), `Dang[]` (density). `getScatterCos()` inverts CDF. `setLorenz(width)` initializes Lorenzian angular distribution. `interpolateDens(x)` interpolates density at arbitrary cosine.

**Status**: Complete and functional. Used in `test_RayScattererMMC.cpp` for CPU MC scattering.

**Pros**: Compact (16 floats per material), fast sampling (linear search, could be binary), supports arbitrary phase functions (tabulated). Lorenzian model good for forward-peaked scattering.
**Cons**: Only 1D angular (scattering angle cosine, not full 3D direction). Azimuthal angle sampled uniformly (assumes azimuthal symmetry). Limited resolution with 16 bins.

### 6. Spherical Harmonics (Not Implemented)

**What**: Expand angular function as $\sum_{l=0}^{N} \sum_{m=-l}^{l} a_l^m Y_l^m(\theta, \phi)$. Rotational invariance, spectral convergence.

**Status**: Not implemented. Discussed in `OctahedralSphereMaping.md` as comparison: "computationally expensive for high-frequency data", "prone to ringing artifacts". Not suitable for spacecraft scattering with sharp shadows.

## Summary: Method Selection Guide

| Method | Regime | Best For | Implemented? | GPU? |
|--------|--------|----------|--------------|------|
| **Radiosity** | Thermal, diffuse, any density | Spacecraft thermal equilibrium, laser heating | Yes (Python 2D/3D, C++) | Occlusion only |
| **Beer-Lambert** | Ionizing, sparse, single pass | X-ray shadow, first-order shielding | Yes (Python+OpenCL) | Yes |
| **MC NEE (single scatter)** | Ionizing, sparse, 1 bounce | Spacecraft neutron/gamma scatter | Design only (MonteCarloScattere.md) | Designed, not coded |
| **MC multi-scatter / bidirectional** | Ionizing, dense, multi-bounce | Deep shielding, reactor penetration | Not implemented | Not implemented |
| **Neutron diffusion** | Dense, optically thick, reactor | Reactor criticality, burnup | Yes (Python 1D spherical) | No |
| **Discrete ordinates (S_N)** | Intermediate, anisotropic | Shielding with angular detail | Unfinished (C++ Scatterer2.h) | No |
| **P_N (spherical harmonics)** | Dense, smooth angular | Reactor core, near-isotropic | Not implemented | No |
| **Clumped attenuation** | Sparse, fast approx | Quick shielding estimates | Design only | Not implemented |

## Key Gaps and Priorities

1. **MC NEE GPU kernel** — The single most impactful gap. Design is complete in `MonteCarloScattere.md` (cone-walking, cluster tiling, local memory). Needs OpenCL implementation. Would enable ionizing radiation scattering on spacecraft.
2. **Occlusion integration into Python 3D radiosity** — `TriangleOcclusionRaytracer.py` exists but is separate from `Radiosity3D.py`. Wiring them together would complete the thermal pipeline.
3. **Scatterer2.h channel bookkeeping** — Fix `backChans` population and `n` variable. Would enable deterministic scattering in C++.
4. **`Scatter.cl` rewrite** — Current placeholder is non-functional. Port the NEE design or the Xray.cl attenuation pattern with scatter extension.
5. **`evalIcosaGrid()` / `accumIcosaGrid()`** — Missing inverse projection for icosahedral sphere map. Needed to store/query directional distributions in scattering solvers.
6. **C++ octahedral mapping completion** — Only `sampleOctahedron()` (minimal) exists. Python version is complete; C++ needs full seamless mapping with reconstruction.
7. **Multi-energy attenuation** — Current Beer-Lambert is monoenergetic. Need spectrum hardening (energy-dependent $\mu$) for realistic gamma shielding.
8. **2D/3D neutron diffusion** — Current solver is 1D spherical only. 2D/3D would enable reactor core simulations beyond spherical implosion.
9. **Multi-scatter GPU kernel** — Extension of NEE to 2+ bounces. Needed for dense shielding. Bidirectional path tracing is the long-term goal.

## Related Audits

- **`continuum-mechanics-impact.md`** — Eulerian compressible multi-material solver. Contains the **Algorithm Review: Heterogeneous Material Simulation Methods** covering all multi-material fluid/solid algorithms. The neutron diffusion solver (`neutron_diffusion_solver.py`) couples with the implosion hydrocode there.
- **`fluid-dynamics.md`** — Fluid dynamics overview. The 1D Lagrangian solvers (`implosion_solver.py`, `lagrangian_tube_solver.py`) couple radiation transport with hydrodynamics.
- **`parallel-particle-cell.md`** — Grid infrastructure and neighbor search. Relevant for spatial clustering of radiation scattering primitives and GPU cell-list acceleration.
- **`spacecraft-design-combat.md`** — Spacecraft design and combat models. The radiation transport tools here serve the spacecraft damage/thermal models described there.
- **`noise-procedural.md`** — Procedural noise generation. Fibonacci sphere sampling and blue noise are related to quasi-MC direction sampling for radiation transport.
