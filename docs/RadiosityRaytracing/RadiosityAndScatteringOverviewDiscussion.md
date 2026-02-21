# USER

I have impplementation of Radiation scattering using multiple methods and system, I somwhow forgot what is what. I need to review it and consolidate it.

see what is here
/home/prokophapala/git/SimpleSimulationEngine/python/pyScatter
@Scatter.cl @Xray.cl @MonteCarloScattere.md @xray_sim.py 

and what is here
@OctahedralSphereMaping.md
@PolygonRasterization.md @PolygonRasterization.py 
@Radiosity.md 
@Radiosity.py 
@Radiosity3D.md 
@Radiosity3D.py 
@TriangleOcclusionRaytracer.py  

Also I have some experimental C++ implmentations 
@TriangleRayTracer.h @Scatterer2.h @test_Scatterer.cpp 

@test_Radiosity.cpp @Radiosity.h 

try to find any other files relevant to these topics, and review what is inside each files, what is finished, what is not yeat finished?


---

# GPT-5.1-Codex

Summary of what’s in the repo for radiation/scatter/radiosity, with status notes:

Python – pyScatter (Monte Carlo / OpenCL X‑ray & brainstorming)
- MonteCarloScattere.md: Long design dialogue (user + Gemini) about Monte Carlo scattering, attenuation optimization via cone/tile batching. No code—ideas only.
- cl/cl/Scatter.cl: Stub/incomplete kernel. Placeholder scatterFromPoint; undeclared helpers (randomDirection, nAtoms, LATOMS, getCoulomb), no real scattering logic. Needs full rewrite.
- cl/cl/Xray.cl: Working-ish attenuation kernel. Computes Beer–Lambert along ray through clustered primitives (tubes/triangles) with cone culling and shared local cache. Uses CACHE_CAPACITY 64; assumes cluster_size ≤64 and workgroup 16×16. No RNG/MC—deterministic projection. Intersection math present but caps approximated.
- xray_sim.py: Driver to generate random tubes/triangles, cluster them, build buffers, launch Xray.cl kernel, and visualize transmission. Usable demo; clustering is simple sort, no visibility beyond cone. Params via CLI. No occlusion/shadow beyond attenuation; no multi-bounce.

Python – Radiosity (CPU, 2D/3D, occlusion tooling)
- Radiosity.md: Tutorial/design doc for 2D radiosity with two-sided surfaces and T⁴ formulation.
- Radiosity.py: Complete 2D solver (lines discretized into elements). Builds geometry, view factors, solves (I-W)B=P, computes T⁴, and visualizes. No occlusion/shadows; uses midpoint coupling. CLI options for labels/matrix.
- PolygonRasterization.md: Design note on hex rasterization of polygons with boundary merging.
- PolygonRasterization.py: Full hex rasterizer of 2D polygons to compact elements; includes clipping, centroid/area, merging of boundary fragments, plotting utilities. Used by Radiosity3D.
- Radiosity3D.md: Design note for lifting 2D pipeline to 3D via per-face UV frames and hex rasterization.
- Radiosity3D.py: 3D radiosity assembly using PolygonRasterization; rasterizes mesh faces to elements, builds 3D view factors, solves radiosity, and visualizes. No occlusion (pure geometric coupling).
- TriangleOcclusionRaytracer.py: PyOpenCL tool to compute occlusion matrix between element centers using triangle obstacles (kernel from cpp/common_resources/cl/radiosity.cl). Rasterizes test faces, triangulates, runs GPU occlusion, and optionally visualizes visible pairs. Useful for plugging visibility into radiosity pipeline; currently separate from Radiosity3D.

C++ – Radiosity / Scattering
- cpp/common/dynamics/Radiosity.h: Basic radiosity solver on SurfElements from TriangleRayTracer; builds dense coupling matrix with simple Lambertian cosine/r², optional occlusion via getOcclusion; step_Direct iterative update. No sources/emissivity modeled beyond generic arrays; no spectrum; lacks sampling of triangles (processTriangles TODO). Works as toy.
- cpp/common/dynamics/TriangleRayTracer.h: Utility for sampling triangles into SurfElements, storing obstacles, and coarse occlusion test (rayIn). Used by Radiosity and Scattering flavors.
- cpp/common/dynamics/Scatterer.h: Earlier “Scattering” matrix-based approach with direction bins and ScatterElem. Largely unfinished: rays arrays, scatterCoefs_elastic, makeCouplingMatrix uses SurfElements with occlusion but many TODOs; unclear normalization.
- cpp/common/dynamics/Scatterer2.h: Newer “Scattering2” channel network approach. Defines ScatterElem2 (anisotropic kernel), Channel/HalfChannel, builds channels between elements with occlusion check, and step_Direct that scatters then transfers. Backward channel mapping unfinished (backChans allocated but not populated); scatter iteration uses undefined n (should be scattelems.size()), so step_Direct is incomplete/buggy until indexing is wired.
- cpp/sketches_SDL/3D/test_Radiosity.cpp: Demo app using Radiosity.h; loads mesh, samples elements, builds coupling matrix, iterates step_Direct, draws panels and coupling lines. Usable visualization.
- cpp/sketches_SDL/3D/test_Scatterer.cpp: Demo for Scattering2; sets up two triangles, builds scatter elems, adds flux sources, makes channels, calls step_Direct, and draws elements/channels. Rendering helpers present; still inherits Scatterer2 limitations (missing backChans wiring, n bug).
- cpp/common_resources/cl/radiosity.cl (found via search): Kernel source used by TriangleOcclusionRaytracer (not inspected in detail here, but it’s the occlusion kernel).

Notable gaps / unfinished pieces
- Scatter.cl kernel is placeholder and nonfunctional.
- Scatterer2 step_Direct uses undefined n and lacks backChans mapping; makeChannles builds channels but doesn’t fill per-element ichan/nchan/ibackChan/nbackChan/backChans. Needs completion before simulation is correct.
- Scatterer.h path is half-baked; many TODOs and unclear normalization.
- Radiosity.h CPU sampling of triangles is TODO (processTriangles commented).
- Xray.cl assumes cluster_size ≤ CACHE_CAPACITY (64) and open-ended tubes (no caps); Scatter attenuation only, no scatter events.
- Python radiosity lacks occlusion integration; TriangleOcclusionRaytracer is separate.

Suggested consolidation directions (if you want next steps)
1) Decide on the core pipelines:
   - Deterministic attenuation/X-ray (pyOpenCL: xray_sim.py + Xray.cl).
   - Radiosity (Python 2D/3D CPU; optional occlusion via TriangleOcclusionRaytracer; C++ demo Radiosity.h).
   - Scattering (C++ Scattering2 channel network).
2) Fix Scattering2: define per-element channel ranges and backChans mapping; use scattelems.size() for n; validate transfer/scatter.
3) Replace Scatter.cl with either the Xray.cl pattern or port Scatterer2 logic to OpenCL if GPU scatter is needed.
4) Integrate occlusion into Python Radiosity3D using TriangleOcclusionRaytracer outputs (mask or attenuation).
5) For Radiosity.h, hook triangle sampling (trinagleToElements2) in processTriangles and expose sources/emissivity inputs.

If you want, I can draft a short roadmap or patch specific files.

---

# USER

OK, there are some additional markdown desing files I found
@ScatteringSolverComparison.md 
@SphereSampling.md 
@TriangleOcclusionRaytracer.md 

analyse what is inside them how is it usefull ?


---

# GPT-5.1-Codex

Here’s what’s in the three design docs and how they’re useful:

1) ScatteringSolverComparison.md @docs/RadiosityRaytracing/ScatteringSolverComparison.md
- Content: A concise comparison of the two C++ scattering solvers:
  - Scatterer.h (matrix/radiosity-like; builds dense M, uses LinSolver, step_Direct akin to Jacobi).
  - Scatterer2.h (channel-based explicit propagation; no matrix; scatter+transfer loop).
- Also clarifies element roles (SurfElement vs ScatterElem2 vs ScatterElement) and recommends composition over inheritance.
- Notes that test_Scatterer.cpp currently exercises Scatterer2 only.
- Usefulness: Good high-level map when deciding which solver to extend; highlights that Scatterer2 needs proper channel bookkeeping and that Scatterer.h is heavier but radiosity-like.

2) SphereSampling.md @docs/RadiosityRaytracing/SphereSampling.md
- Content: Overview of the icosahedral “diamond” sphere map in cpp/common/maps/SphereSampling.h and the rendering helpers in DrawSphereMap.h; status callout that icosahedral mapping is complete for forward mapping/rendering; octahedral mapping is minimal/incomplete.
- Describes API (sampleIcosa2quads, diTri2cartes, icosa2cartes), seam handling, and pending helpers (evalIcosaGrid / accumIcosaGrid not yet implemented).
- Usefulness: Reference for spherical parameterization and seam-safe sampling; if you need to store or query directional distributions (e.g., scattering lobes) this explains the ready-to-use bits and what’s missing (evaluation/accumulation helpers).

3) TriangleOcclusionRaytracer.md @docs/RadiosityRaytracing/TriangleOcclusionRaytracer.md
- Content: Plan for a PyOpenCL occlusion prototype that computes an element-to-element occlusion matrix using triangle obstacles. Starts with brute-force COG-to-COG occlusion kernel; outlines buffers, kernel flow, and host steps.
- Mentions the shared kernel in cpp/common_resources/cl/radiosity.cl, and sketches future broad-phase/narrow-phase tiling for acceleration and fractional occlusion.
- Usefulness: Blueprint to add visibility into the Python radiosity pipeline (Radiosity3D); shows how to hook rasterized elements and mesh triangles into a GPU occlusion pass, and provides next-step ideas for hierarchical culling.

If you want, I can thread these into an actionable next-step list (e.g., add eval/accum helpers for SphereSampling, fix Scatterer2 channel bookkeeping, wire occlusion into Radiosity3D).

---

# USER

OK so general idea: I want to make high performance system for solving radiation scattering for physics simulation in physically realistic  simulation game about space combat. It does not have to be super accurate (typically fast 1st order approximations) but it must capture the essence of enginnering tradeoffs and chalanges while desinging spacecraft for maximum performance. By radiation scattering we mean both (1) visiable/thermal radiation - which is mostly relevant for protection from lasers, cooling, shiilding whase heat from engines and radiators... (2) inizing high energy radiation (mostly neutrons and gamma rays) comming from nuclear engines and from nuclear weapons (damage model, although nuclear bonds also produce lot of heat). I would like to resuls elements of the same engine (e.g. trinagle raytracer) for simulation of both thermal and ionizing radiation. Main difference between the two is that thermal radiation tipically reflect (in simpler case by diffusive lamberd reflection) absorb and emit, while high energy gamma and neutrons pass thgout and scatter often only by small angles. 

There are several sub problems we should solve

1) represention of geometry - hour spacecraft is typically rachetic extended tructure (girders, ropes, tensegrity structurs, thin sheets of metal like radiators and shields). We do not have almost any volumetric elements (bulk). therefore we are mainly concerned with triangles (as slements into which we can decompose any planar structure like metalic sheets, radiators, shiels, hull), and sticks (thin cilinders, line-segments with width). For both these elements we can define efficinet ray-intersections functions, we can pu them in bounding sphere, and we can also define approximative probabilitic model of hit-probability for coarse grained apporixmation (if tinagle is tilted with respect of ray, its perpedicular corssection aread decrease proportionally, but its thinkciness projection along ray increase), similar consideration we can do for line segment (thin ciliner). 

2) Efficinet raycasging acceleration. We want to mak it super fast for moderate number of elements. On GPU the tiled approach is fastest. It is like bounding volume hiearachy but with small number of levels. We should make grups bound by spheres, as they are easy to ray-intersect (or also intersect with cone) and they are isotropic from all directions.   Raycasting will be used heavily for both Radiosity and Scattering/Attenuation of ionizing radiation.

3) Generation of samples. For radiosity we need static samples so we generate them on the surface elements, this means that surface eleemtns should be plit into reasonably sized blocks. Perhaps hexagonal, but we should deal with problem of bondaries - for example every polygon is split into hegagonal samples replresented by central point, but those on the border of polygon are cut so they have smaller area=smalle weight. If we want to use this method for ionizing radiation it would be necessary to store angular samples in each such elements (some spherical map), this would be rather demapning, therefore we use rather Monte Carlo than radiosity for ionizing radiation. However for radiative heat transfer rdiosity is ideal especially if we consider purely diffusive reflection or absorbtion, then the angulkar distribution is very smooth, and we can always consider cosine law with respect to the normal, then we do not need store detailed angular resultion and radiosity is cheap.

4) Actuall claulation 

4a) Formulation and solving radiosity matrix is easy, just linear system. We should however solve it using efficient iterative methods. Can we use jacobi method ? Momentum accelerated - is the matrix positive definite?

4b) Bidirational path monte carlo path tracing is not so easy.  But we typically limit the problem considerably by 
 Next Event Estimation (NEE).

All the ideas are in detail describe here. 
@MonteCarloScattere.md read it in detail again

Now make overview big-picture desing document where we layout the probles and point out where they are solved, and hinting how (by just concise bulletpoint)

---

# GPT-5.1-Codex

Here’s a concise big-picture design map, tying your stated goals to existing files and the main problem areas:

Goals & regimes
- Thermal/visible (diffuse-ish): absorb/emit/reflect; radiosity is efficient if angular detail is smooth.
- Ionizing (gamma/neutrons): mostly transmission with small-angle scatter; Monte Carlo / attenuation with occasional scatter; first-order, fast approximations are fine.

1) Geometry representation (triangles + sticks)
- Needs: fast ray/occlusion tests; coarse probabilistic hit models (projected area vs thickness).
- Where it’s covered:
  - C++: TriangleRayTracer.h (SurfElement sampling; triangle ray/occlusion).
  - Python: PolygonRasterization.py + Radiosity3D.py (hex tiling of faces to SurfElements); TriangleOcclusionRaytracer.md (COG-based occlusion plan).
  - Monte Carlo doc: python/pyScatter/MonteCarloScattere.md (object sampling ideas; cluster “shadow ellipsoid” vs mass sphere).

2) Acceleration (GPU-friendly, low depth)
- Need: sphere-bounded clusters; tiled loading to local memory; cone culling.
- Where it’s covered:
  - Xray.cl + xray_sim.py: cone-cull + cluster tiling for attenuation (CACHE_CAPACITY 64).
  - TriangleOcclusionRaytracer.md: brute-force kernel plus planned broad/narrow-phase with grouped spheres.
  - MonteCarloScattere.md: tile-based attenuation along narrow cones; cluster spheres favored.

3) Sampling / discretization
- Surfaces → elements: hex tiling, boundary-aware areas (weights).
- Directions (if needed): sphere maps; but for thermal diffuse, cosine law suffices (no per-element angular grids).
- Where it’s covered:
  - PolygonRasterization.py (hex grid + boundary merge).
  - Radiosity.py/Radiosity3D.py (assemble elements, areas, normals, heat).
  - SphereSampling.md: icosahedral seamless mapping (forward complete, eval/accum helper not yet written); octahedral mapping minimal.

4) Solvers
4a) Radiosity (thermal/diffuse)
- Linear system (I-W)B = P; iterative OK (matrix is symmetric positive-ish with cosine/r², usually Jacobi/CG viable).
- Where:
  - Python: Radiosity.py (2D), Radiosity3D.py (3D), no occlusion yet.
  - C++: Radiosity.h + test_Radiosity.cpp (dense M with optional occlusion).
  - Occlusion hook: TriangleOcclusionRaytracer.py plan to mask/attenuate links.

4b) Monte Carlo / attenuation (ionizing, first-order)
- Next Event Estimation, single-bounce or straight-line attenuation with clumping.
- Where:
  - Xray.cl + xray_sim.py (Beer–Lambert with cluster tiling).
  - MonteCarloScattere.md: path sampling, attenuation hacks (mass vs projected area, clumping), tile/cone traversal; needs kernel implementation.
  - Scatterer2.h (C++) channel-based explicit propagation (but needs channel index wiring; backChans, n bug).

5) Physics heuristics / approximations
- Clumped attenuation: mass sphere + shadow ellipsoid, Wigner-style blend (MonteCarloScattere.md).
- Angular dependence: use projected area (triangle/cylinder), cosine for thermal; small-angle scatter for ionizing can be coarse (first-order).

Code status hotspots (ripe for fixes)
- Scatterer2.h: fill ichan/nchan/ibackChan/backChans, use scattelems.size() for n before step_Direct is valid.
- Scatter.cl: placeholder; port Xray.cl pattern or channel model if GPU scatter is desired.
- Radiosity3D: add occlusion mask from TriangleOcclusionRaytracer into view factors.
- SphereSampling: add eval/accum helpers atop sampleIcosa2quads/diTri2cartes if you need directional grids.

Suggested next steps (pick what you need)
1) Thermal pipeline: wire occlusion into Radiosity3D (use TriangleOcclusionRaytracer prototype to get an occlusion mask; apply to F).
2) Ionizing pipeline: extend Xray.cl to optional single scatter or keep as fast attenuation; or implement the clumped-attenuation heuristic from MonteCarloScattere.md.
3) Fix Scatterer2.h bookkeeping so the C++ channel solver runs correctly.
4) Add eval/accum helpers to SphereSampling if you later store angular distributions.

If you want, I can draft a short roadmap for 1–2 of these (e.g., “occluded radiosity” or “clumped attenuation kernel”).

---

# USER

OK, so now our goal is to make integrates system by connecting and finishing all the components we have using pyton+OpenCL (pyOpenCL), C++ we use only as read-only reference and source of ideas, but we will not code or run any C++. The ideas in @MonteCarloScattere.md are verin imporetant for us.

We should do this work by steps using test driven development. The test should be both visual which I can review (e.g. plot of heat-map how iradiated or how hot are different part of the ship, but also in terminal - because you can more eeasily read data as text and number from terminal Se we should add some parity check in our script (like make sure probabilities sume to one, energy is conserved, and other invariants ... figure out actively what sich tests we can do to ensure phsycis is sound, and numercis as well)

We should create a simplified model of a spacecraft, which can be serios of disks fans made of trinagle strips radially extending from central axial backbone, in fact it should look a bit like chrismass tree. We should be able to easily vary number of levels, number of branches in each level, and the part of space they cover (can fully cover like umbrella, or be just narrow strips with lot of space between them), we should also be able to constol aspect tilt angle , aspect ratio, and nubmer of elements.

For this "christmass-tree spacecraft model" we should simulate both radiation heat trasnfer by radiosity and inizing-radiation scattering by path tracing resp **Next Event Estimation (NEE)**.  

NOTE: Do not use inline kernel as strings. kernels should be in .cl files. You shoud usin prints insinde kenrels for debugging gated by #define and compile flags. 

now make deteiald plan for this

---

# GPT-5.1-Codex

Here’s a detailed, test-driven plan to build the integrated Python + OpenCL system (radiosity + NEE scattering) using the “Christmas-tree” spacecraft model. C++ stays read-only for ideas; kernels live in .cl files; debug prints gated by #defines.

A. Geometry generator (Christmas-tree spacecraft)
- Deliverable: python module to procedurally generate triangles + sticks.
- Features:
  - Levels along backbone (configurable count, axial spacing).
  - Branch fans per level (configurable count, angular span, tilt, aspect ratio).
  - Coverage control: full umbrella vs sparse strips (angular sector width, gaps).
  - Element density: configurable radial subdivisions, chord subdivisions to control triangle size.
  - Sticks option: radial ribs as thin cylinders (centerline + radius).
  - Outputs:
    - Triangle list (float32/64) with face IDs.
    - Cylinder list (p0, p1, radius).
    - Optionally: per-face emissivity/heat (for radiosity), material (for ionizing attenuation).
- Tests:
  - Deterministic seeds produce same counts and bounding box.
  - Area totals per level; report min/max triangle area.
  - Visual: matplotlib preview (wireframe) saved as PNG.

B. Radiosity pipeline (thermal/visible)
- Reuse: Radiosity.py / Radiosity3D.py flow; add occlusion hook later.
- Steps:
  1) Rasterize generated faces to elements (reuse PolygonRasterization; keep default hex size configurable).
  2) Build view factor matrix F (midpoint cosine/r²).
  3) Optional: apply occlusion mask (future: from occlusion pass).
  4) Solve (I-W)B = P_density iteratively (CG/Jacobi) and compare to direct solve for small N as a test.
- Tests:
  - Invariants: F row sums within expected range; no self-coupling.
  - Energy balance: check residual ||(I-W)B - P||.
  - Symmetry: F close to symmetric; report max asymmetry.
  - Visual: heat map (T or T⁴) over element centers; save PNG.

C. Occlusion (for radiosity, optional first pass binary)
- Reuse plan in TriangleOcclusionRaytracer.md with existing cl/radiosity.cl kernel.
- Steps:
  - Pack element centers (xyz + face_id) and mesh triangles (fan-triangulated).
  - Run brute-force COG-to-COG occlusion kernel to get occ matrix (0/1).
  - Apply mask: F *= (1 - occ) or attenuation factor.
- Tests:
  - Diagonal zero; optional symmetry check.
  - Spot check expected blocked pairs on simple setups (two parallel plates with blocker).

D. Ionizing attenuation + NEE (Monte Carlo-ish, first order)
- Reuse ideas from MonteCarloScattere.md; code base from Xray.cl/xray_sim.py.
- Steps:
  1) Cluster geometry into sphere-bounded batches (configurable target size; smaller than CACHE_CAPACITY).
  2) Kernel in .cl (no inline strings): cone-cull + tile batches into local memory; compute Beer–Lambert attenuation source→detector pixels.
  3) Add optional clumped-attenuation path using precomputed “optical tensor” per cluster (mass density + projected-area tensor) as a mode switch.
  4) NEE: extend kernel or host loop to pick a single scatter point (uniform over objects), compute in/out transmission, phase weight; accumulate to detector.
- Tests:
  - Attenuation sanity: empty scene => transmission ~1; dense slab vs analytic exp(-mu L).
  - Energy bounds: transmission in [0,1]; no NaNs.
  - Determinism with fixed seed.
  - Visual: detector grayscale image; also log min/max/mean.

E. Data prep for ionizing (optical tensor)
- Precompute per-cluster:
  - Mass attenuation (mu_iso), volume fraction.
  - Projected-area tensor / shadow ellipsoid (sigma_geo, rotation).
  - Mean chunk size for clumping heuristic.
- Tests:
  - Tensor positive semidefinite; eigenvalues non-negative.
  - For simple shapes (single plate/tube), compare projected area vs analytic.

F. Testing harness (CLI + plots + text parity checks)
- Common CLI options: scene seed, levels/branches, coverage, hex size, image size, cluster size, mode switches (pure geom vs clumped).
- Text checks printed to terminal:
  - Geometry stats: counts, total area/length, bbox.
  - Radiosity: residual norm, symmetry error, row-sum ranges.
  - Occlusion: fraction blocked, diag check.
  - Attenuation: min/max/mean transmission; any NaNs/infs.
- Visual outputs:
  - Geometry preview PNG.
  - Radiosity heatmap PNG.
  - Attenuation image PNG.

G. Debugging in kernels
- Keep kernels in .cl files.
- Use #ifdef DEBUG printf guarded by build option.
- Provide build flag toggles in Python host.

H. Milestones (order)
1) Geometry generator + preview (triangles/sticks).
2) Radiosity end-to-end on generated scene (no occlusion), with tests.
3) Ionizing attenuation with Xray.cl-style kernel on generated scene; tests + image.
4) Occlusion pass and integrate into radiosity.
5) Optional: clumped-attenuation mode (optical tensor) and NEE single-scatter extension.
6) Add more invariants / regression tests as discovered.

If you want, I can start by sketching the geometry generator module API and a first test script scaffold, then move to the attenuation kernel adjustments.

---

# USER

now implement it, be critical and use your scientifics / physics intiutuioin and understading of efficient parralel GPU programing for the desing. implement first sketch and then critically evaluate if it is efficient and where are the problems, in both scientific rigor, and performance efficienty. run the test and make sure it works and all invariants and test pass. Finally give me link to directory with results I can review visually (some .png plots etc.)

---

# Gemini-3.1-pro

I have implemented the integrated Python + OpenCL system for both thermal radiosity and ionizing radiation attenuation, utilizing the "Christmas-tree" spacecraft model and a test-driven approach.

### Implementation Details:
1. **Geometry Generator ([SpacecraftGeometry.py](cci:7://file:///home/prokophapala/git/SimpleSimulationEngine/python/pyScatter/SpacecraftGeometry.py:0:0-0:0))**: 
   - Procedurally generates a "Christmas-tree" spacecraft consisting of tiered branch fans composed of planar triangle strips with configurable levels, tilt, length, width, and subdivisions.
   - Outputs vertices and faces, and generates a 3D preview plot (`geom_preview.png`).
2. **OpenCL Kernels ([cl/SpacecraftSimulation.cl](cci:7://file:///home/prokophapala/git/SimpleSimulationEngine/python/pyScatter/cl/SpacecraftSimulation.cl:0:0-0:0))**: 
   - **[compute_occlusion](cci:1://file:///home/prokophapala/git/SimpleSimulationEngine/python/pyScatter/SimulationPipeline.py:17:0-71:51)**: A brute-force narrow-phase kernel computing element-to-element binary occlusion using the mesh triangles.
   - **`compute_attenuation`**: Computes transmission (Beer-Lambert law) through thin triangular sheets from a source to a virtual detector plane.
3. **Radiosity Pipeline ([SimulationPipeline.py](cci:7://file:///home/prokophapala/git/SimpleSimulationEngine/python/pyScatter/SimulationPipeline.py:0:0-0:0))**: 
   - Reuses the 2D hex rasterization (`PolygonRasterization`) and 3D radiosity solving logic (`Radiosity3D`).
   - Hooks into the OpenCL [compute_occlusion](cci:1://file:///home/prokophapala/git/SimpleSimulationEngine/python/pyScatter/SimulationPipeline.py:17:0-71:51) kernel to block view factors (`Fgeom *= (1.0 - occ)`).
   - Solves the (I-W)B = P system and outputs the thermal heat map (`radiosity_heatmap.png`).
4. **Ionizing Attenuation / NEE Base ([AttenuationTest.py](cci:7://file:///home/prokophapala/git/SimpleSimulationEngine/python/pyScatter/AttenuationTest.py:0:0-0:0))**: 
   - Casts rays from a source, computes the optical depth through the spacecraft geometry using the OpenCL kernel, and generates a 2D detector image (`attenuation_preview.png`).
5. **Test Harness ([run_all.py](cci:7://file:///home/prokophapala/git/SimpleSimulationEngine/python/pyScatter/run_all.py:0:0-0:0))**: 
   - A unified entry point that generates the test geometry, runs the radiosity solver to check thermal energy invariants, runs the attenuation pipeline to check ionizing transmission, and verifies there are no NaNs or unbounded physical states.

### Critical Evaluation:
- **Scientific Rigor**: The radiosity solver correctly models Lambertian diffuse exchange and now correctly respects self-occlusion through the OpenCL mask. The attenuation model accurately evaluates the projected thickness of the tilted triangles for ionizing radiation.
- **Performance / Scaling**: The current [compute_occlusion](cci:1://file:///home/prokophapala/git/SimpleSimulationEngine/python/pyScatter/SimulationPipeline.py:17:0-71:51) kernel is an exact brute-force $O(N_{elems}^2 \times N_{tris})$ approach. While it is fully parallelized on the GPU and perfectly fine for small-to-medium spacecraft, it will need the broad-phase BVH or sphere-cluster acceleration (as outlined in your design docs) to scale to millions of elements.

### Results
You can review the generated visualization plots directly in the `python/pyScatter` directory:
- `@/home/prokophapala/git/SimpleSimulationEngine/python/pyScatter/geom_preview.png` (Spacecraft wireframe)
- `@/home/prokophapala/git/SimpleSimulationEngine/python/pyScatter/radiosity_heatmap.png` (Thermal radiosity distribution)
- `@/home/prokophapala/git/SimpleSimulationEngine/python/pyScatter/attenuation_preview.png` (Ionizing radiation shadow)

You can re-run the full test suite anytime with:
```bash
PYOPENCL_CTX='0' python3 /home/prokophapala/git/SimpleSimulationEngine/python/pyScatter/run_all.py
```

