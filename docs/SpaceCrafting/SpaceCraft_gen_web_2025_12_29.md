






## Status

[ ] Node blocks (addCube/addCMesh octa)  
[ ] Girder bridges (bridgeFacingPolygons)  
[ ] Rope generation (`mesh.rope`)  
[ ] Slider anchors/paths (make_anchor_point equivalent + path export)  
[ ] Wheel/ring generation and slider binding  
[ ] Radiator panel generation between rails (Quad/Parametric)  
[ ] Shield/pusher-plate/dish (Parabola/Slab)  
[ ] Nozzle/dish attachments (parabolic patch + connections)  
[ ] Weld/auto-connections between components  
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
- In JS editor (`js/spacecraft_editor`): the engine still uses the legacy `BuildCraft_truss` pipeline (`SpaceCraftEngine.processCommands` calls `BuildCraft_truss(this.mesh, this.craft)`) and only supports Nodes + Girders (materials TODO, no blocks/bridges/sliders/plates). So the browser editor can render a simple truss, but not the new block-based craft yet.

Main subproblems (from docs and code):
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