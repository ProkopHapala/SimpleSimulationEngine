






## Status (JS editor, 2025-12-29)

[x] Node blocks (addCube)  
[x] Girder bridges (bridgeFacingPolygons)  
[x] Rope generation (`mesh.rope`)  
[x] Radiator panels between rails (ParametricQuadPatch) with welded edges  
[x] Weld/auto-connections between components (edge→girder)  
[~] Ring (Wheel) generation  
[~] Slider anchors/paths (placeholder anchor + path viz; not exported)  
[ ] Wheel/ring slider binding  
[ ] Shield/pusher-plate/dish (Parabola/Slab)  
[ ] Nozzle/dish attachments (parabolic patch + connections)  
[ ] Material/stick catalogs in JS craft  
[ ] GUI hooks: run BuildCraft_blocks_js, load sample script  
[ ] Worker commands: Node/Girder/Rope/Ring/Slider/Plate/Material  
[ ] Connectivity check + OBJ export in GUI  

Status from docs (all in `docs/SpaceCrafting`):
- `SpaceCrafting_new.md` is the current, detailed overview of the whole pipeline (Lua → SpaceCraft → BuildCraft_blocks → Mesh::Builder2 → OBJ/truss) and open issues. It already captures legacy vs new paths, sliders, and export gaps.
- `SpaceCraftingWithBlocks_new.md` is a concise write‑up of the block-based construction implementation and wheel/slider attachment workflow.
- [SpaceCrafting.md](cci:7://file:///home/prokophapala/git/SimpleSimulationEngine/docs/SpaceCrafting/SpaceCrafting.md:0:0-0:0) is the older architectural overview; good background but partly superseded by the two “_new” docs.
- `SpaceCraftingWithBlocks.md` and `SpaceCraftConstructionProblems.md` were reviewed; the “_new” versions are the authoritative ones. The problems doc is mostly folded into section 9 of `SpaceCrafting_new.md`.

Can we already generate a spacecraft?
- Yes, via C++: Lua scripts → `BuildCraft_blocks` → OBJ/truss (headless `spaceCraftMeshExport`, or interactive apps). Core structural components (nodes as blocks, girders as bridges, rings, sliders, ropes) are implemented. Plates/tanks/welds are not yet wired in for the block path.
- In JS editor (`js/spacecraft_editor`): now uses `BuildCraft_blocks_js` to generate blocks, bridged girders, ropes, rings, radiators (ParametricQuadPatch), welds, and triangle faces rendered in MeshRenderer. Sliders are placeholder anchors/paths only. Materials/catalog, OBJ export, and GUI hooks to invoke the new builder are still TODO.


## Main subproblems (from docs and code):

1) Migration to block-based generator:
   - C++ block path covers nodes/girders/rings/sliders/ropes; plates, welds, tanks, thrusters are TODO in `BuildCraft_blocks` @cpp/common/Orbital/SpaceCraft2Mesh2.h#555-594.
   - Standalone truss-view path (`SpaceCraft2Mesh_blocks.h`) lacks sliders/plates/welds.
2) Sliders & export gaps:
   - Sliders’ EdgeVertBond constraints are not serialized in the text truss export (only points/sticks) @cpp/common/Orbital/SpaceCraft2Mesh_blocks.h#130-201; reloading a `.truss` loses sliders.
3) Plates/shields/radiators:
   - Implemented in legacy truss path, but commented out in block path; need attachment strategy on block/girder faces.
4) Thrusters/dishes/nozzles:
   - Geometry generators exist elsewhere; not wired into block pipeline; need attachment/damping pattern.
5) JS editor feature gap:
   - Uses legacy truss mesh builder and minimal commands (Node, Girder). No blocks/bridges, no sliders, no materials/catalog, no Lua ingestion.
6) Integrity/duplication:
   - Mesh checks exist (`checkSpaceCraftMesh`) but block path still allows TODO components; duplicate vertices/edges avoidance remains important per project norms.

What’s missing / suggested next steps:
C++ pipeline:
- Wire plates (radiators/shields) into `BuildCraft_blocks`, using `plateOnGriders`/slab generators; record pointRange/stickRange.
- Add welds/tanks/thrusters generation in block path, or explicitly mark unsupported to fail fast.
- Extend `exportSimToFile` to optionally write slider EdgeVertBond data (and read counterpart if needed).
- Add regression script using `spaceCraftMeshExport` on a sample Lua ship to validate block path + connectivity check.

JS editor:
- Port editor pipeline to block-based builder: expose block node creation and face-bridging girders (match C++ `nodeBlock_to_mesh` + `girderBlocks_to_mesh`).
- Expand worker command set: materials, rings, sliders, ropes; map shadow IDs robustly; carry material/stick data into mesh gen.
- Add OBJ export from MeshBuilder so browser can inspect results; add connectivity diagnostics akin to `checkSpaceCraftMesh`.
- Integrate Lua loading (if desired) or mirror the C++ build order logic for bound nodes/sliders.

Overall readiness:
- Headless/desktop C++: can already generate spacecraft meshes (block-based core is in place) and export OBJ/truss; missing components limit full ships.
- Web editor: only simple truss ships with nodes/girders; block-based generation and richer components remain to be implemented.

## Actionable plan for JS block-based spacecraft generation (web editor)

1) Replace legacy JS builder with block pipeline
- Implement `BuildCraft_blocks_js(mesh, craft)` in `js/spacecraft_editor/js/SpaceCraft2Mesh.js`:
  - Nodes → blocks: for each abstract node, call `mesh.addCube(pos, size, true)` (or octahedron if available) and store chunk/vert ranges for face picking.
  - Girders → bridges: for each girder, pick most-facing faces between its two node chunks (use `mesh.bridgeFacingPolygons(pA, pB, chRangeA, chRangeB, nseg, stickTypes, stickTypesFaces)`) to generate a bridged truss instead of single edge. Reuse nseg/up parameters from girder data.
  - Ropes: reuse `mesh.rope()` between node anchor verts (from node ranges).
  - Panels (radiators/shields): between selected edges/rails, use `mesh.ParametricQuadPatch` or `mesh.triPlateBetweenEdges` to span quads/strips; map to component ranges.
  - Rings/sliders/wheels: add ring generator (reuse girder/wheel from MeshBuilder if present) and placeholder slider anchors (can use `mesh.make_anchor_point` analog if exposed); at minimum reserve data structures so later constraints can bind paths.
  - Return/update component → mesh range maps so GUI can select and export.

2) Extend worker command set and script ingestion
- Expand `SpaceCraftWorker` to accept Node/Girder/Rope/Ring/Slider/Plate commands (similar to Lua schema). Mirror `ship_ICF_marksman_2.lua` structure: materials, nodes, girders, rings + sliders, radiators, dish/nozzle panels.
- Add material/stick catalogs to JS craft (even minimal) so edges get type IDs for coloring/export.

3) Hook GUI to new builder
- Add dropdown/actions in `GUI.js` to run `BuildCraft_blocks_js` instead of `BuildCraft_truss` and to load a sample script (JS or Lua-parsed JSON) into the worker.
- Keep MeshGenTest panel for diagnostics; add “Build Craft (blocks)” button next to existing run script.
- Add “Show triangles” toggle if desired; faces already render via MeshRenderer.faceMesh.

4) Reuse existing mesh generators (already tested)
- From `MeshGenTestGUI` (`js/spacecraft_editor/js/MeshGenTestGUI.js`):
  - Panels/slabs/dishes: `QuadSheet`, `QuadSlab`, `ParametricQuadPatch`, `ParabolaSheet`, `ParabolaSlab_wrap`, `Parabola_Wire_new`, `ParametricParabolaPatch`.
  - Tubes/slabs: `TubeSheet`, `TubeSheet_swapped`, `TubeSheetBond`, `SlabTube`, `SlabTube_wrap`.
  - Ropes + patch mix: `QuadParametricRopes`, `addSkipEdges`, `addEdgesFromPairs`.
  - Bridging: `bridge_quads` via `bridgeFacingPolygons`.
  - Blocks: `addCube` (returns chunk range for faces).
- From `constructionBlockTests.js`:
  - `mesh.addCube`, `mesh.bridgeFacingPolygons` (auto face selection), `mesh.girder1`, `mesh.rope`, `mesh.triPlateBetweenEdges`, `mesh.ParametricQuadPatch`, `mesh.ParametricParabolaPatch`.

5) Connectivity and export
- Add OBJ export from MeshBuilder in GUI to inspect builds.
- Keep notes that slider/path constraints are not yet exported; plan a JSON sidecar to serialize EdgeVertBond equivalents once sliders are implemented.

6) Test cases to wire first
- Recreate `ship_ICF_marksman_2.lua` in JS commands: nodes as cubes, bridged girders, ropes, a ring + placeholders for sliders, one radiator panel between two rails, one dish/nozzle using `ParametricParabolaPatch`.
- Regression buttons in GUI: (a) minimal 2-node bridge, (b) cube bridge test (from ConstructionBlockTests), (c) radiator quad patch, (d) rope V + plate, (e) nozzle dish.

## Design philosophy (higher-level mesh view)
- The spacecraft is a tensegrity-like higher-level mesh:
  - Nodes = blocks = mesh “vertices.”
  - Girders/ropes = mesh “edges.”
  - Plates (radiators/shields) = mesh “faces” spanning between two rails/edges.
  - Dishes/nozzles/tanks = added a posteriori as larger surface/volume elements attached to the frame.
- Sliders enable dynamic reconfiguration: they let one edge endpoint slide along another edge (girder/rope), providing controlled motion; paths must be explicit and visually distinct.
- Wheels are a special girder in circular form: sliders along the wheel enable rotational motion (momentum wheels / landing gear analogs).
- This hierarchy should guide both data model (abstract modules) and mesh builders (block/bridge/face generators).
- Think of `SpaceCraft` → mesh as **mesh refinement**: the abstract graph is the low-res mesh; `BuildCraft_blocks`/`BuildCraft_truss` perform refinement by replacing each high-level edge (girder/rope) with many detailed edges/verts and each node with a block; see also `docs/TrussGeneration/truss_low_to_hi.md` for the low→high truss pattern.

## Testing strategy (JS, aligned with C++ pipeline)
- Ship definition script (JS), mirroring `SpaceCraft` / `SpaceCraftComponents` (@js/spacecraft_editor/js/SpaceCraft.js ↔ @cpp/common/Orbital/SpaceCraft.h, @cpp/common/Orbital/SpaceCraftComponents.h):
  - Create abstract modules: Node, Girder, Rope, Ring/Wheel, Slider (with path side), Radiator/Shield panel, Nozzle/Dish.
  - Use as sample for the GUI “MakeShip” button so it can be run repeatedly.
- Mesh build step: call `BuildCraft_blocks_js` from @js/spacecraft_editor/js/SpaceCraft2Mesh.js to convert the abstract ship into concrete mesh; keep behavior close to C++ `BuildCraft_blocks` in @cpp/common/Orbital/SpaceCraft2Mesh2.h.
   - The mesh should be cleared before building.
- Slider visualization: when building mesh, emit:
  - Path segments as colored lines (distinct color from structure) along the chosen girder side/ring side.
  - Slider anchor vertex as a highlighted point; optional small glyph/line to show the slider body.
- Hook “MakeShip” in GUI:
  - Button executes the JS ship script, then runs `BuildCraft_blocks_js`, then triggers renderer update.

## Slider (rail + moving joint) design notes — JS editor

What is a slider (concept)
- A slider is a joint that lets one structural component move along a defined path on another structural component. Think “carriage on a rail”: the rail component provides the vertex path; the sliding component contributes a single vertex that travels along that path while the rest of the component stays rigid. The slider holds the relative transform along the rail via a floating parameter (`cur`) and enforces a constraint/weld between the sliding vertex and the current edge segment on the rail.

What the C++ reference does (current state)
- Slider is Node-derived in C++, with `Path path` (ps[], n, cur, closed), `updatePath(rail, side)` using rail->sideToPath, `move` integrating `cur` with drive params, and `StructuralComponent::updateSlidersPaths` picking nearest side, sharing paths, and snapping to nearest path point.
- Path semantics: ordered vertex indices along a rail; `cur` = edge index + interpolation; `closed` for rings; nearest-point search used to initialize `cur`.
- Rail candidates: girder, ring, rope (any StructuralComponent implementing sideToPath/findNearestPoint). Rings set `closed=true`.

Design correction for JS (per new requirement)
- Slider is **not** a Node anymore (nodes are cubes). Slider is a movable joint between two structural components:
  - Rail component: provides the path (edge loop).
  - Sliding component: owns the sliding vertex that is constrained/welded to a point on the path.
- Slider stores references: rail {id,type}, sliding {id,type}, sliding vertex index (in sliding component’s pointRange), path {ps[], closed, cur}, side, method flag (stride vs SDF), plus drive params (force/power/speed/springK/Kdv/maxDist).

Path construction methods (JS)
1) Stride-based (bPickSideByStrides=true)
   - Use rail pointRange and known vertex ordering to derive edge loops by side.
   - Implement `sideToPath`/`getEdgePathVert` for girders (open loop) and rings (closed).
2) SDF-based (bPickSideByStrides=false)
   - Use cylindrical SDF along target edge(s) via `selectVertsBySDF` (SDfuncs cylinder).
   - Sort selected verts along edge direction to form ordered ps[]; mark closed for ring.

Initialization flow
- Choose rail + side (explicit or inferred from nearest rail vertex modulo layout).
- Build path via chosen method; set `closed` if rail is ring.
- Find nearest point on path to the sliding vertex position → set `cur`; optionally snap sliding vertex to interpolated position.

Constraint/welding behavior
- Sliding vertex belongs to sliding component; at runtime/mesh build, bind it to the path segment defined by `cur` (interpolate between path vertices i,i+1; wrap if closed).
- For debugging/logging: store/report path length, ps list, side, method used, current segment indices, and welded vertex IDs.

Data model fields to add in JS
- Slider: { railId, railType, sliderCompId, sliderCompType, sliderVertId, side, methodFlag(bPickSideByStrides), path:{ps:[], closed:false, cur:0}, drive:{maxDist, forceMax, powerMax, maxSpeed, springK, Kdv, mass}, shareId? for shared paths }.
- StructuralComponent helpers: `sideToPath(side)` (stride) and `pathFromSDF(side)` to populate ps[]; `findNearestPoint` for initialization.

Planned JS build step (BuildCraft_blocks_js)
- After both rail and sliding components are built (so pointRange known), build Slider:
  - Build path via stride/SDF per global flag.
  - Initialize `cur` by nearest point; snap sliding vertex position if desired.
  - Record mapping: slider → {rail pointRange, sliding vertex id}.
  - Emit weld/constraint: connect sliding vertex to path segment endpoints used by `cur` (MeshBuilder weldListToRange or dedicated constraint hook).

Test scenario (to validate both methods)
- Two nearby girders; pick one girder edge as rail.
- Sliding component = second girder; pick one of its verts as sliding vertex.
- Run twice: (a) stride path, (b) SDF cylinder path. Log ps[], cur, chosen segment vertices, weld result.
## Slider implementation progress (JS, 2025-12-29)
- **Data model**: Slider now stores rail, sliding component, sliding vertex id, side, methodFlag (stride/SDF), path (ps/cur/closed), weldDist, pathRadius; worker/engine/API carry these params end-to-end.
- **Path selection**:
  - **Stride**: Approximates `edgePathVert` by sampling along the chosen girder edge with tiny radius/threshold so only that edge is picked.
  - **SDF**: Uses a corner-offset cylinder with very small radius/threshold; ring remains SDF (closed loop). Fallbacks/logs if empty.
- **Visualization**:
  - Rail path drawn **green** (stick type 5).
  - Sliding vertex linked to interpolated path point in **cyan** (stick type 2).
  - Test script offsets sliding girder in Z for visual clarity.
- **GUI**:
  - New **Components** panel with slider selection dropdown and `c_along` control.
  - Live updates via `engine.rebuildMesh()` (no worker round-trip).
  - Slider list refreshed after each script run.
- **Diagnostics**:
  - Node test `js/spacecraft_editor/tests/diag_slider.js` validates SDF radius/threshold and offsets.
- **Issues solved**:
  - **NaNs** from missing `rail.up` -> default up + robust normalization.
  - **SDF picking all verts** -> reduced radius (0.05) + edge offset.
  - **GUI syntax error** from stray brace fixed.
  - **Worker/engine API** expanded to support joint-style slider parameters.
