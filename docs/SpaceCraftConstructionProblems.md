# SpaceCraft Construction Problems (WIP)

This document refines some of the higher‑level challenges listed in `SpaceCrafting_new.md` into more concrete geometric / algorithmic problems that should be solved in the code. Each problem is meant to translate into:

- **Ship‑design side** (Lua / `SpaceCraft` components + parameters).
- **Mesh‑generation side** (`BuildCraft_blocks`, `Mesh::Builder2`, UV‑based generators).
- **Simulation side** (truss stiffness, sliders, welds, damping).

For now the main focus is a **recoil damper for pulsed propulsion nozzles / pusher plates** implemented as a telescopic truss.

---

## 0. Overview and TODO checklist

This section is a compact checklist / roadmap. Details are in the referenced sections.

### 0.1 Geometry and mechanics problems

- [ ] **[P1] Telescopic truss recoil damper design**  (see §1.1–1.5)
  - Hexagonal outer tube + triangular inner tube, slider rails, collision avoidance.
- [ ] **[P2] Plates on girders and ropes**  (see §2)
  - Triangular / quad plates between polylines, matching tessellation to truss.
- [ ] **[P3] Welds and automatic connections between components**  (see §3)
  - Weld manager operating on selections, rigid vs soft connections.

### 0.2 JS editor / prototyping tasks

Implemented items (JS / editor only):

- [x] **[J0a] `dirMask` parsing and binary input**  (see §1.9)
- [x] **[J0b] Slab offsets vs thickness split**  (see §1.9)
- [x] **[J0c] Shared GUI parameter dictionaries and `buildCommonArgs`**  (see §1.9)
- [x] **[J0d] Length unification and torus radii (`L`, `R + thick`)**  (see §1.9)
- [x] **[J0e] Grid / axis visibility wiring in JS editor**  (see §1.9)
- [x] **[J0f] `SlabTube` twist control wired from GUI**  (see §1.9)

Pending JS‑first tasks:

- [ ] **[J1] SlabTube / TubeSheet parameter playground panel**  (see §1.8–1.9)
- [ ] **[J2] Visualizers for index / angle / SDF rail selection strategies**  (see §1.7–1.8)
- [ ] **[J3] TubeSheet post‑processor for second‑neighbor longitudinal rails**  (see §4.1)
- [ ] **[J4] JS selection manager for named vertex/edge/face sets**  (see §4.2)
- [ ] **[J5] JS stick‑material support and edge coloring by material**  (see §4.3)

---

## 1. Telescopic Truss Recoil Damper for a Nozzle / Pusher Plate

### 1.1 Physical scenario and requirements

We want a **recoil damper** for a large plate or nozzle (e.g. Orion‑like nuclear pulse pusher plate, or a magnetic nozzle) using a **telescopic mechanism**:

- **Two coaxial truss tubes sliding inside each other**, like a telescopic shock absorber.
- Inner tube carries the **pusher plate / nozzle**, outer tube is attached to the **main ship truss**.
- Purpose: absorb high‑frequency impulses while keeping the plate aligned.

Key **requirements**:

- **Low probability of collision** between sliding trusses.
- **Clear central corridor** inside the outer tube (no sticks crossing the axis) for shock‑absorber internals or cabling.
- **Straight, perfectly axial sliding edges** where the two tubes touch.
- **Large clearance everywhere else**: non‑sliding nodes should stay as far apart as practical.
- Sliding surfaces should be **geometrically simple** (preferably straight axial rails) so we can:
  - Parameterize them as `Slider` paths.
  - Generate them robustly from UV / slab primitives.

### 1.2 Proposed geometry: hexagonal outer truss, triangular inner truss

A simple and symmetric configuration:

- **Outer tube**: hollow **hexagonal truss tube**.
- **Inner tube**: slender **triangular truss tube**.

Rationale:

- Hexagon provides **six straight axial rails** (at the 6 vertices) with a natural **triangular subset** (every second vertex) for contact with the inner triangle.
- Triangle is **statically determinate and stiff** while leaving lots of empty space inside the hexagon.
- We can choose **3 contact rails out of 6** to get:
  - Three **contact lines** where the inner triangle touches the outer hex tube (via sliders / bushings).
  - Three **spare rails** for other purposes, or simply left free for better clearance.

#### 1.2.1 Inner triangular tube (slender girder)

Concept:

- A **triangular prism** extruded axially.
- Each of the **3 rectangular faces** is braced by **diagonal X‑bracing** to make a stiff lattice.

Mesh / generator mapping (desired):

- Logically: a `Girder`‑like component with **3 corner rails** + internal cross‑bracing.
- Could be approximated or prototyped via:
  - A **triangular `TubeSheet` / `SlabTube`** variant with `nSides = 3`.
  - Or by building 3 **offset ropes / girders** and then adding diagonal sticks between them.

Constraints:

- Inner triangle must stay **slender** and well inside the inscribed circle of the outer hexagon.
- **No diagonals** should protrude into the clearance corridor needed for sliding against the outer hex.

#### 1.2.2 Outer hexagonal tube (double‑layer truss tube)

Concept:

- Build a **double‑shell hexagonal tube**:
  - **Inner hexagon** of vertices.
  - **Outer hexagon** (larger radius), rotated **by 30°** (half of 60° period) relative to inner.
  - Axially, the inner and outer shells are also **offset by half a segment**.
- Connect inner ↔ outer hexagon by **short diagonal struts**, 4 per inner vertex, forming local pyramids.
- Equivalent to taking a **2‑layer slab structure** and wrapping it into a tube.

This resembles an FCC/BCC‑like pattern folded into a cylinder:

- Each node has several **short, roughly equal‑length struts**, giving isotropic stiffness.
- The **contact rails** (sliding edges) are specific chosen vertex columns that are kept straight and smooth.

Mesh / generator mapping:

- Conceptually: `SlabTube` or `TubeSheet` with extra **radial & diagonal connections** between an inner and an outer layer.
- Implementation ideas using existing code:
  - Use `SlabTube` (in `DrawUV.h`) with a suitable `dirMask` to get:
    - Axial edges (tube length).
    - Azimuthal edges (around circumference).
    - Diagonals and radial links between inner / outer rings.
  - Alternatively, create a `TubeSheet` for a single layer, then a second, slightly larger tube, and connect via `slabEdges`‑like logic.

We must ensure **specific azimuthal tracks** (e.g. every second vertex) are:

- Perfectly **axial** (no zig‑zags in radius or angle).
- Have **uniform spacing** suitable for smooth slider motion.

### 1.3 Sliding contact design and collision‑avoidance constraints

We want to define **three pairs of matching axial rails**:

- On the outer hex tube: pick 3 out of 6 vertex tracks.
- On the inner triangular tube: its 3 vertices.

For each sliding pair we require:

- **Collinearity / coaxiality**: inner‑tube rail vertex path should lie inside or very near the projection of outer‑tube rail.
- **Monotonic axial coordinate**: vertex indices increase monotonically along the axis.
- **Constant (or controlled) radial offset** between inner and outer rails to define the bushing gap.

Collision avoidance:

- All **non‑slider edges / nodes** must stay outside some **safety radius** around each contact rail.
- Ideally, the UV parametrization of both tubes allows us to explicitly tag **rail vertices** and keep diagonal bracing shy of the rail corridor.

For use with existing `Slider` machinery (`SpaceCraftComponents.h`):

- Sliding rails become **`Path`s** over mesh vertices.
- Each `Slider` node has:
  - An anchor vertex on the inner or outer tube (`Slider::ivert`).
  - A `Path` built from a list of vertex indices on the counter‑part tube.
- Dynamics glue (`sliders2edgeverts`, `applySliders2sim`) then keeps anchors constrained to those rails while allowing motion along the axis.

### 1.4 How existing generators can help

Relevant generators and where they live:

- **`TubeSheet`** – `DrawUV.h`
  - Generates a **single‑layer tube** from UV coordinates using `ConeUVfunc`.
  - DirMask controls which edges exist: axial, azimuthal, diagonals.
  - Already used in `constructionBlockApp` via `-TubeSheet` CLI option.
  - Good building block for **one shell** of a truss tube.

- **`SlabTube`** – `DrawUV.h`
  - Builds a **double‑layer tube (slab)**: inner and outer surfaces with connections.
  - Uses `UV_slab` + `slabEdges` to create axial, circumferential, and diagonal edges across the two layers.
  - Very close to what we need for the **double‑shell hex tube**, but currently tuned for a circular / conical section.

- **`QuadSheet` / `QuadSlab`** – `DrawUV.h`
  - `QuadSheet`: single‑layer quad mesh with configurable edge directions and optional clipping.
  - `QuadSlab`: double‑layer quad slab with a thickness offset and flexible connectivity.
  - These are the **flat equivalents** of `TubeSheet` / `SlabTube` and can be used to prototype the lattice in flat space before folding into a tube.

- **`TorusSheet` / Parabola / etc.** – `DrawUV.h`
  - Dishes / toroidal structures for nozzles and pusher plates themselves.
  - Out of immediate scope here, but relevant for attaching the telescopic damper to a plate.

In `constructionBlockApp.cpp` there are ready‑made **CLI test harnesses**:

- `-TubeSheet`, `-SlabTube`, `-QuadSheet`, `-QuadSlab`, `-TorusSheet`.
- Useful for interactively iterating on the **UV parameters** and **dirMask** for:
  - Getting clean **axial rails**.
  - Verifying that cross‑bracing does not intrude into the sliding corridors.

### 1.5 Concrete sub‑tasks for implementation

This section is meant as a checklist that can later become code TODOs.

- **[T1] Hexagonal double‑shell tube prototype**
  - Use `SlabTube` or a custom variant to build a **hexagon‑like** tube:
    - Either by a dedicated `ConeUVfunc` variant for 6‑sided cross‑sections.
    - Or by starting from a circular tube and sampling only 6 azimuthal directions.
  - Ensure **6 straight axial vertex tracks** exist and can be addressed by simple index formulas.

- **[T2] Triangular inner tube prototype**
  - Create a `TubeSheet`‑like function specialized to **3 sides**.
  - Add **X‑bracing** on each of the 3 rectangular faces (can reuse `QuadSlab` / `UV_panel` patterns).

- **[T3] Rail tagging and slider path extraction**
  - For both tubes, define a **convention** (e.g. vertex index modulo `nSides`) to identify rail vertices.
  - Implement helpers in `Mesh::Builder2` or higher level that:
    - Given `nSegments` and `nSides`, return a list of vertex indices for a specific rail.
  - Wire these into `Slider::updatePath` or a new utility so sliders can be automatically attached along rails.

- **[T4] Clearance and collision checks**
  - Define a **clearance radius** around each rail.
  - During generation, avoid placing diagonal edges/vertices inside this radius.
  - Optionally, add a debug mode that uses `Mesh::Builder2` selections (e.g. SDF cylinders) to verify no vertices violate the clearance.

- **[T5] Integration with `SpaceCraft` components**
  - Decide on a **component representation**:
    - One option: new `TelescopicGirder` (inner + outer as sub‑components).
    - Simpler option for now: build separate `Girder`‑like components for inner and outer tubes, plus `Slider`s linking them.
  - Provide **Lua builders** exposing key geometric parameters: lengths, radii, `nSegments`, `nSides` (3/6), clearances, etc.

### 1.6 Using `SlabTube` for polygonal tubes (n = 3, 4, 6)

`SlabTube` in both C++ (`DrawUV.h`) and JS (`MeshesUV.js`) is more general than a “circular tube”:

- **Cross‑section is defined by a finite azimuthal subdivision count** `n.y`.
- For small `n.y`, the circular section becomes effectively **polygonal**:
  - `n.y = 3` → effectively a **triangular** cross‑section.
  - `n.y = 4` → effectively a **square** cross‑section.
  - `n.y = 6` → effectively a **hexagonal** cross‑section.
- The radial profile is still defined by `ConeUVfunc` (linear interpolation between `R1`, `R2`), but with constant radii (`R1 ≈ R2`) it behaves like a cylinder sampled at `n.y` angular directions.

In C++ (`DrawUV.h`):

- `SlabTube` uses a twisted UV mapping over a 2D grid and builds two layers:
  - Inner layer at radii `Rs.a, Rs.b`.
  - Outer layer at `Rs.a+up.z, Rs.b+up.z`.
  - The twist `dudv = 0.5*(n.x-1)/(n.y-1)` skews azimuthal coordinate slightly with axial position.
- Geometry is built via `UV_slab`:
  - `UV_slab_verts` creates vertices on both inner and outer shells.
  - `slabEdges` uses `dirMask` and `stickTypes` to add axial, azimuthal, radial and diagonal edges between the two layers.

In JS (`MeshesUV.js`):

- `SlabTube` mirrors the C++ logic:
  - Same `dudv` and `ConeUVfunc` mapping.
  - Uses `UV_slab` with the JS versions of `UV_slab_verts` and `slabEdges`.
- With small `n.y` the resulting mesh is already a **polygonal double‑shell tube** suitable as a base for:
  - Inner triangular tube: e.g. `n.y = 3`.
  - Outer hex tube: e.g. `n.y = 6`.

For the telescopic recoil damper this suggests a **“parameter‑specialized”** approach instead of new generators:

- Inner **triangular tube**:
  - Use `SlabTube`/`TubeSheet` with `n.y = 3` and appropriate `dirMask`.
  - Add extra X‑bracing using `UV_panel` / `QuadSlab` patterns if needed.
- Outer **hexagonal double shell**:
  - Use `SlabTube` with `n.y = 6` and `up.z > 0` for shell separation.
  - Tune `dirMask` and `stickTypes` to get a dense, roughly isotropic truss between inner and outer shells while keeping some rails clean for sliding.

**JS playground idea:** add a small GUI panel in `js/spacecraft_editor` to play with `SlabTube` parameters (`n`, `Rs`, `L`, `up`, `dirMask`) and regenerate the tube interactively, to build intuition before committing to fixed patterns in C++.

### 1.7 Strategies for selecting straight edge / vertex paths

For sliders and telescopic motion we need **straight, smooth paths** along the edges of tubes/girders. There are three complementary strategies already supported or partially supported in the C++ infrastructure.

#### 1.7.1 Strategy A – Index/stride based paths (generator‑aware)

**Idea:**

- Use knowledge of the **exact vertex layout** produced by a specific generator to identify paths by index pattern.
- Example in `SpaceCraftComponents.h`:
  - `Girder::sideToPath` assumes 4 vertices per segment (`mseg = 4`) and contiguous packing in `pointRange`.
  - It builds a rail for `side ∈ {0,1,2,3}` by `inds[i] = i0 + 4*i + side`.
  - `Ring::sideToPath` follows the same pattern.

**Pros:**

- Very cheap and deterministic.
- Does not need `edgesOfVerts` or any graph traversal.

**Cons:**

- **Brittle w.r.t. generator changes**: any change in vertex ordering, `mvert`, or how segments are laid out breaks the convention.
- Requires a clearly documented contract per generator: `nSegments`, `mvert`, and mapping from `(segmentIndex, sideIndex)` to a global vertex index.

**Implication for tubes:**

- For `SlabTube` / `TubeSheet` we should design and document similar index formulas so that:
  - Given `(nSegAxial, nSides, layer)` we can compute `iv(ia, iside, ilayer)`.
  - A rail is then defined by fixing `iside, ilayer` and running `ia = 0..nSegAxial-1`.

#### 1.7.2 Strategy B – Angle‑based edge walks on the mesh graph

**Idea:**

- Treat the mesh as a graph and **walk along edges** that deviate as little as possible from a target direction.
- Start from an edge (or vertex + preferred direction) and extend the path by repeatedly choosing the next edge whose direction keeps the cosine with the current direction above some `cosMin`.

Relevant helpers in `Mesh::Builder2` (`MeshBuilder2.h`):

- `selectVertEdgeAngle(iv, hdir, cosMin, ie_ignore)` – pick a neighbor edge of vertex `iv` whose direction has `dot(hdir, eDir) ≥ cosMin`.
- `selectEdgeStrip(ie, iv, cosMin, nmax, bUpdateDir)` – follow a single edge strip starting from edge `ie`.
- `selectEdgeStrip2(ie, cosMin, nmaxs, bUpdateDir)` – more flexible variant for longer or branched strips.

**Pros:**

- Works for **arbitrary meshes**, even if generated geometry changes.
- Naturally finds “straightest” paths following actual geometry instead of assumed indexing.

**Cons:**

- Requires good `edgesOfVerts` connectivity info.
- On symmetric meshes (e.g. hex tubes) it may pick a different rail than intended unless guided carefully.

**Use for telescopic tubes:**

- Once the tube axis is known, choose a starting vertex near the desired rail (e.g. via SDF selection), then:
  - Use `selectEdgeStrip2` with a high `cosMin` to trace a nearly straight axial edge path.
  - Record the resulting ordered vertex list as a candidate slider path.

#### 1.7.3 Strategy C – Geometric SDF‑based selection (lines / cylinders)

**Idea:**

- Define an ideal **geometric line** (axis or polyline) and select all vertices lying within a certain distance from that line.
- Implement this via signed‑distance functions (SDFs): e.g. distance to an infinite line or to a capsule.

Existing tools in `Mesh::Builder2` and apps:

- `selectVertsAlongLine(p0, p1, r, bSort)` – select vertices inside a cylinder of radius `r` around the line from `p0` to `p1`, optionally sorted.
- `selectVertsAlongPolyline(r, bSort, n, edges)` – generalization to polylines.
- `selectVertsBySDF(sdf, threshold)` – fully generic selection using an arbitrary SDF (can be used to implement line/cylinder selection in a more flexible way).

**Pros:**

- Purely geometric and **independent of generator internals**.
- Directly expresses “vertices near this ideal axis”.

**Cons:**

- Produces an unordered set of vertices, which then needs to be **sorted** (e.g. by projection onto the axis) to become a path.
- Choice of radius is delicate: too small → gaps; too large → multiple rails get merged.

**Use for telescopic tubes:**

- Define the mechanical axis (e.g. as a line between two design points).
- Use `selectVertsAlongLine` (or equivalent SDF‑based selection) to get candidate rail vertices.
- Sort them by projection onto the axis direction and optionally intersect with the connectivity graph to ensure they form a valid edge path.

### 1.8 Debugging and visualization plan (JS spacecraft editor)

Before encoding anything into `SpaceCraft` / C++, we can prototype and compare all three strategies in the JS editor (`js/spacecraft_editor`).

- **[D1] SlabTube / TubeSheet parameter playground**
  - Add a small GUI panel to control `MeshesUV.SlabTube` / `TubeSheet` parameters:
    - `n.x`, `n.y` (especially `n.y = 3,4,6`).
    - `Rs.x`, `Rs.y`, `L`, `up.z`.
    - `dirMask`, `stickTypes`.
  - Regenerate the tube on each change and render:
    - All edges.
    - Optional labels for vertex/edge indices (similar to `MeshBuilder2Draw` in C++).

- **[D2] Visualize index‑based rails (Strategy A)**
  - Implement a JS helper that, given `(nSegAxial, nSides, layer)`, computes `iv(ia, iside, ilayer)` according to the chosen convention for `SlabTube`.
  - Draw those vertices/edges in a highlight color to verify that the stride logic matches the actual mesh.

- **[D3] Visualize angle‑based edge walks (Strategy B)**
  - Build `edgesOfVerts` in JS for the current mesh.
  - Provide a tool to pick a starting edge or vertex, then:
    - Run a `selectEdgeStrip2`‑like algorithm with adjustable `cosMin` and `nmax`.
    - Draw the resulting strip as a highlighted path with ordered labels.
  - Compare visually with the stride‑based rail from [D2] on the same tube.

- **[D4] Visualize SDF line selections (Strategy C)**
  - Allow definition of a line/axis in the scene (e.g. by clicking two points or tying it to the tube axis).
  - Implement a JS version of `selectVertsAlongLine` to:
    - Select all vertices within radius `r` around the axis.
    - Sort them by axial coordinate and draw them with index labels.
  - Optionally, overlay this SDF‑selected path with the paths from [D2] and [D3] to see where each method agrees or diverges.

The result should be a **practical toolkit** in JS to debug and tune:

- Which `SlabTube` parameter regimes produce good polygonal tubes.
- How robust each path‑selection strategy is on those meshes.
- What conventions we can safely bake into C++ (`SpaceCraftComponents` and mesh builders) for long‑term use.

### 1.9 JS prototyping checklist (current status)

**Implemented in JS (MeshGenTestGUI, MeshesUV, editor):**

- [x] **[J0a] `dirMask` parsing and binary input**
  - GUI accepts binary strings (e.g. `1011`) and JS now parses them with an explicit radix (base 2) before passing them into `MeshesUV` generators.
- [x] **[J0b] Slab offsets and thickness split**
  - GUI parameters separate **UV offsets** (`offsetX`, `offsetY`) from **radial separation** (`thick`) for slab/tube generators and use them consistently.
- [x] **[J0c] Shared GUI parameter dictionaries and `buildCommonArgs`**
  - `MeshGenTestGUI` uses shared `defaultParams` / `params_tube` / `params_slab` and a common argument builder to reduce boilerplate and keep Tube/Slab/Quad/Torus shapes consistent.
- [x] **[J0d] Length unification and torus radii**
  - Quad generators use `L` instead of ad‑hoc `scale`, and `TorusSheet` uses a base radius `R` plus `thick` (`R_major = R + thick`) instead of separate `R_major`/`R_minor` controls.
- [x] **[J0e] Grid / axis visibility wiring in JS editor**
  - `chkShowGrid` / `chkShowAxis` checkboxes now correctly control the Three.js helpers and apply their initial state on startup.
- [x] **[J0f] `SlabTube` twist control**
  - JS `SlabTube` accepts a `twist` parameter (defaulting to the old hardcoded value) and the GUI passes its `twist` slider into the generator, making twist behavior consistent with `TubeSheet`.

**Pending JS prototyping tasks (before C++ hardening):**

- [ ] **[J1] SlabTube / TubeSheet playground panel**
  - Implement the parameter playground from [D1] directly in the editor, with live regeneration and optional index/edge labels for tubes.
- [ ] **[J2] Visualizers for Strategies A/B/C**
  - Implement JS versions of [D2]–[D4] (index‑based rails, angle‑based edge walks, SDF line selections) and integrate them into the editor UI for side‑by‑side comparison.

---

## 2. Plates on Girders and Ropes (radiators, sails, shields) – outline only

**Problem statement (sketch):**

- Attach **plates** (radiators, sails, shields, radiators‑as‑slabs) onto existing **girder / rope networks**.
- Two main geometric cases:
  - Polylines sharing a **common vertex** (`L` or `V` shape) → **triangular plate**.
  - Two disjoint polylines (`| |`) → **quad / trapezoid**.
- Need robust methods to:
  - Match segment counts and parameterizations along each polyline.
  - Generate triangular / quad grids (thin sheets or slabs) between them using `QuadSheet`, `QuadSlab`, etc.
  - Keep tessellation compatible with underlying truss for load transfer and damage.

**Relevant tools:**

- `QuadSheet`, `QuadSlab`, `UV_panel`, `UV_slab`, `slabEdges` in `DrawUV.h`.
- `plateBetweenEdges` GUI function in `constructionBlockApp.cpp`.

(*Detailed design to be filled in later.*)

---

## 3. Welds and Automatic Connections Between Components – outline only

**Problem statement (sketch):**

- Need a **selection / weld manager** to connect existing components:
  - Pick vertex sets (by SDF / selection tools) on two components.
  - Create a `Weld` that introduces bonds between them (within some radius `Rmax`).
- Welds should support both:
  - **Rigid connections** (short, stiff sticks).
  - **Soft / damped connections** (using `StickMaterial` with high damping, or explicit ropes).

**Relevant tools:**

- `Weld` component and `weld_to_mesh` in `SpaceCraft2Mesh2.h`.
- Selection and SDF tools in `Mesh::Builder2` (`selectVertsBySDF`, `bridgeTriPatch`, etc.).

(*To be elaborated once the telescopic truss problem is clearer and integrated.*)

---

## 4. JS editor infrastructure for truss tubes and mesh debugging

This section lists JS‑first features that are useful for prototyping telescopic trusses, sliders, and welds. They can later be ported into C++ once patterns stabilize.

### 4.1 Post‑processor for TubeSheet second‑neighbor rails

**Motivation:**

- `TubeSheet` with a small number of azimuthal subdivisions (e.g. `n.y = 2, 3`) and `twist = 0.5` naturally produces **polygonal tubes**:
  - `n.y = 2` → square‑like cross‑section.
  - `n.y = 3` → hex‑like cross‑section.
- With `twist = 0.5`, every **second longitudinal segment** in UV maps to vertices that are rotated by an integer multiple of `2π / n.y`, i.e. coincide angularly.
- To turn these into **clean, straight rails for sliders** (girders, recoil dampers), we want explicit **second‑neighbor longitudinal edges** connecting vertices separated by a configurable segment stride (e.g. `dSeg = 2`).

**Design (JS prototype):**

- Keep `TubeSheet` itself simple: no new `dirMask` bits for second‑neighbor edges.
- Add a **post‑processor** that runs immediately after tube generation:
  - Knows `n.x`, `n.y` and vertex packing (strides) for the current tube.
  - For each vertex on the base layer, optionally add **exactly one extra edge** to a vertex in segment `i + dSeg`.
  - `dSeg` should be configurable (default `2`), and the mapping can cyclically wrap indices in azimuth to follow the twisted polygons correctly (e.g. connect `1 → 2`, `3 → 1`, etc.).
- This post‑processor should operate purely on the mesh graph (vertices + edges) so it can be reused for both **circular** and **polygonal** tubes as long as the indexing convention is known.
- Long term, these second‑neighbor edges form the **primary slider rails** (stiffer material) while the original `TubeSheet` edges remain lighter supporting structure.

**TODO:**

- [ ] **[J3] Implement a JS TubeSheet post‑processor that adds second‑neighbor longitudinal edges for polygonal tubes** (parameterized by `dSeg` and azimuthal offset, operating on the editor’s mesh representation).

### 4.2 Selection manager in the JS editor

**Motivation:**

- Many algorithms (welds, rail extraction, sliders, plate generation) need to operate on **named sets of vertices / edges / faces**.
- A flexible selection manager in the editor makes it easier to:
  - Prototype these algorithms.
  - Debug selections visually.
  - Reuse the same selections across multiple operations.

**Design (JS prototype):**

- Introduce a **selection registry** in JS:
  - Stores multiple selections, each with an **ID / name** and a **type** (vertices, edges, later faces).
  - Each selection is a set or ordered list of indices into the current mesh.
- GUI support:
  - A small **selection panel** with dropdowns or list boxes to choose the current **active selection** for each type.
  - Buttons to **save**, **replace**, **append to**, and **clear** selections from the current editor state or algorithm output.
- Algorithms (e.g. future **weld manager**) should:
  - Take one or two selection IDs as input (e.g. `Selection A`, `Selection B`).
  - Operate on their indices (e.g. create springs/welds between nearby vertices of A and B).
- Rendering:
  - Draw active selections with highlight colors or thicker lines to make debugging intuitive.

**TODO:**

- [ ] **[J4] Implement a JS selection manager that can store, name, and visualize multiple vertex/edge selections and feed them into algorithms (e.g. welds, rail extraction).**

### 4.3 Edge coloring by stick material (JS `SpaceCraftWorkshop`)

**Motivation:**

- In C++ (`SpaceCraftComponents.h`), `SpaceCraftWorkshop` manages **stick materials** used for different edge types (longitudinal, radial, diagonals, etc.).
- Bringing a similar concept into JS and **coloring edges by material** would:
  - Make it much easier to debug mesh generators.
  - Visually distinguish important rails (e.g. second‑neighbor slider paths) from lightweight bracing.

**Design (JS prototype):**

- Implement a lightweight **`SpaceCraftWorkshop`‑like module in JS** that:
  - Defines a set of **stick material types** (IDs) with mechanical meaning (stiff rail, light brace, crosslink, weld, etc.).
  - Assigns each material a **distinct color** and basic rendering style.
- Extend JS mesh data structures so each edge can carry a **material ID**.
- Update mesh generators and post‑processors in JS to:
  - Tag edges according to their role (e.g. rails vs support, radial vs azimuthal).
  - In particular, edges added by the TubeSheet post‑processor in §4.1 should get a **stronger rail material**.
- Rendering:
  - Modify the editor’s renderer to draw edges using the color associated with their material.
  - Optionally add a **legend** in the GUI that shows material names and colors.

**TODO:**

- [ ] **[J5] Implement JS stick‑material support and color edges by material/type in the editor, mirroring the `SpaceCraftWorkshop` concept.**
