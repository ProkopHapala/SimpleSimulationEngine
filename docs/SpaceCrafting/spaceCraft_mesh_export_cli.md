# SpaceCraft Mesh Export CLI
Standalone command‑line tool for testing spacecraft generation without SDL/OpenGL rendering or truss simulation.

## Goal

- **Input**: Lua script that builds a spacecraft (e.g. `tests_bash/Orbital/data/ship_ICF_marksman_2.lua`).
- **Processing pipeline**:
  - Lua → C++ builders in `EditSpaceCraft.h`.
  - Populate `SpaceCraft` from `SpaceCraft.h` / `SpaceCraftComponents.h`.
  - Convert `SpaceCraft` into a `Mesh::Builder2` truss/blocks mesh via `SpaceCraft2Mesh2.h` (e.g. `BuildCraft_blocks`).
  - Optionally run `checkSpaceCraftMesh` for integrity diagnostics.
- **Output**: `.obj` file written by `Mesh::Builder2::write_obj` (from `MeshBuilder2.{h,cpp}`) and stdout debug log controlled by global `verbosity` from `globals.h`.

Designed for **fast, headless regression tests** of spacecraft construction and mesh conversion.

## Reused components (IN scope)

- **Lua→SpaceCraft construction**
  - `cpp/common/Orbital/EditSpaceCraft.h`
    - Owns global Lua state `theLua`, `SpaceCraft *theSpaceCraft`, and `SpaceCraftWorkshop workshop` in namespace `SpaceCrafting`.
    - Provides `initSpaceCraftingLua()` that:
      - creates `theLua`, opens Lua libs, registers C API functions:
        - `Material`, `StickMaterial`, `Node`, `BoundNode`, `Rope`, `Rope2`, `Girder`, `Ring`, `Ring2`, `Slider`, `Weld`, `Gun`, `Thruster`, `Tank`, `Radiator`, `Shield`, `Balloon`, `Rock`.
      - sets `theSpaceCraft->workshop = &workshop`.
    - Uses helpers from `LuaHelpers.h`.
  - `cpp/common/Orbital/SpaceCraft.h`
    - `SpaceCraft` container with `nodes`, `girders`, `ropes`, `rings`, `sliders`, `radiators`, `shields`, `tanks`, `pipes`, `thrusters`, `balloons`, `rocks`, `welds`, `components`, etc.
    - Helper methods used by Lua wrappers: `make_Rope`, `make_Girder`, `make_Ring`, `make_Ring2`, `add_Node`, `add_Rope`, `add_Girder`, `add_Ring`, `add_Radiator`, `add_Shield`, `add_Tank`, `add_Thruster`, `add_Gun`.
    - Diagnostics: `checkIntegrity()`, `find_mesh_element()`, `printAll_*()`.
  - `cpp/common/Orbital/SpaceCraftComponents.h`
    - Definitions of all component types used in construction and mesh conversion (`Node`, `Girder`, `Rope`, `Ring`, `Slider`, `Radiator`, `Shield`, `Tank`, `Thruster`, `Weld`, etc.).

- **SpaceCraft → Mesh conversion**
  - `cpp/common/Orbital/SpaceCraft2Mesh2.h`
    - **Used**:
      - `BuildCraft_blocks(Builder2& mesh, SpaceCraft& craft, double max_size=-1, double node_scale=1.0)`
        - Block‑based mesh generation using node meshes (`nodeBlock_to_mesh`), `girderBlocks_to_mesh`, `ring_to_mesh`, `slider_to_mesh`, `rope_to_mesh`.
      - `checkSpaceCraftMesh(const Builder2& mesh, const SpaceCraft& craft, ...)` for connectivity diagnostics.
    - **Not used in this CLI**:
      - `exportSim(TrussDynamics_d&, ...)`, `exportSim(TrussDynamics_f&, ...)`.
      - `applySliders2sim`, `sliders2edgeverts`, `makeBBoxes`, `makePointCunks` (simulation‑related).

- **Mesh representation & OBJ write‑out**
  - `cpp/common/geometry/MeshBuilder2.h` + `MeshBuilder2.cpp`
    - `Mesh::Builder2` with:
      - `verts`, `edges`, `tris`, `blocks`, `chunks`, `strips`.
      - Topology helpers used by SpaceCraft conversion (`vert`, `edge`, `rope`, `wheel`, `girder1`, `plateOnGriders`, `make_anchor_point`, `addCMesh`, etc.).
      - `void write_obj(const char* fname, uint8_t mask=0xFF) const;`
    - Diagnostics: `printSizes()`, `checkAllPointsConnected()`, selection helpers, etc.

- **Global utilities and logging**
  - `cpp/common/globals.h`
    - `static int verbosity = 1;`
    - `_assert` macro with `exit_on_error` depending on `DEBUGBUILD`.
    - Used everywhere for conditional logging and asserts.
  - `LuaHelpers.h`, `IO_utils.h`, `fastmath.h`, `Vec2.h`, `Vec3.h`, `Mat3.h`, `quaternion.h`, etc. (header‑only math/utility infrastructure).

- **CMake integration**
  - Top‑level: `cpp/CMakeLists.txt`
    - Already defines options `WITH_LUA`, `WITH_SDL`, etc.
    - Adds `cpp/apps` subdirectory only when `WITH_SDL` is ON (we will keep the new target under that umbrella for now so Lua is already configured in `apps/OrbitalWar/CMakeLists.txt`).
  - `cpp/apps/OrbitalWar/CMakeLists.txt`
    - Uses `find_package(Lua52 REQUIRED)`, adds Lua include dirs and `-DLUA`.
    - Provides mesh/spacecraft‑related object libraries via:
      - `$<TARGET_OBJECTS:MeshBuilder2>`
      - plus others (Truss, SDL2OGL, TrussDynamics) for GUI apps.

## Explicitly excluded components (OUT of scope for this CLI)

- **Rendering / GUI / SDL/OpenGL**
  - `AppSDL2OGL_3D`, `Draw2D.h`, `Draw3D.h`, `SDL_utils.h`, `GUI.h`, `SpaceCraftGUI.h`, etc.
  - Visual apps:
    - `spaceCraftEditor`, `SpaceCraftEditorNew`, `spaceCraftDynamics`, `spaceCraftDynamicsOCL`, `constructionBlockApp`, and all other SDL‑based apps.
  - Any OpenGL display lists, event handling, or camera control (`SpaceCraftDynamicsApp`, `SpaceCraftEditorApp`).

- **Truss simulation and dynamics**
  - `TrussDynamics_d`, `TrussDynamics_f`, `Truss.h`, `SoftBody`, `DynamicOpt`, and OCL code.
  - Any `runSim`, `exportSim`, `makeBBoxes`, `sliders2edgeverts`, `applySliders2sim`, solver setup (`LinSolverTrussDynamics`, `linSolveMethod`, etc.).
  - Spacecraft dynamics wrapper `SpaceCraftDynamicsApp.h` and all control/slider motion.

- **Unrelated physics/visualization utilities**
  - Raytracing / radiosity visualization used by `spaceCraftEditor.cpp` (`TriangleRayTracer`, `Radiosity`, `visualizeRayTracing`, etc.).
  - Any code that depends on real‑time interaction or animation.

## New CLI application design

### Target name and location

- **Source file**: `cpp/apps/OrbitalWar/spaceCraftMeshExport.cpp`.
- **Executable target**: `spaceCraftMeshExport` (built into `cpp/Build/apps/OrbitalWar/`).
- **Build conditions**:
  - Only when `WITH_SDL` and `WITH_LUA` are both `ON` (same as existing OrbitalWar apps, to reuse Lua configuration).
  - Does **not** link SDL/OpenGL libraries; only common core + Lua.

### Command‑line interface

- Program name: `spaceCraftMeshExport`.
- Basic usage:

  ```bash
  spaceCraftMeshExport -s data/ship_ICF_marksman_2.lua -o out/ship_ICF_marksman_2.obj [-v 0|1|2]
  ```

- Options (implemented via existing `LambdaDict`/`process_args` in `IO_utils.h` / `argparse.h`):
  - **`-s <lua_file>`** (required)
    - Path to Lua spacecraft script, e.g. `data/ship_ICF_marksman_2.lua`.
  - **`-o <obj_file>`** (optional, default: `ship.obj`)
    - Path to OBJ file to write (`write_obj`).
  - **`-v <int>`** (optional)
    - Sets global `verbosity` (from `globals.h`).
    - `0` – quiet (only fatal errors).
    - `1` – normal (key steps and summary; default).
    - `2+` – verbose diagnostics (e.g. integrity checks, mesh size, some component stats).

- Output channels:
  - **OBJ file**: geometry only; exactly what `Mesh::Builder2::write_obj` writes.
  - **stdout**: textual log suitable for `tee` in bash tests (e.g. `| tee OUT-spaceCraftMeshExport`).

### Main steps executed by the CLI

1. **Initialize stdout**
   - Disable buffering (`setbuf(stdout, NULL)`) so debug output appears in real time.

2. **Parse arguments**
   - Use `LambdaDict` / `process_args(argc, argv, funcs)` as in `spaceCraftDynamics.cpp` and `spaceCraftEditor.cpp`.
   - Validate that `-s` (Lua path) is provided; if missing, print usage and exit non‑zero.

3. **Initialize Lua + SpaceCraft**
   - Allocate spacecraft:
     - `theSpaceCraft = new SpaceCraft();` (global from `EditSpaceCraft.h`).
   - Call `initSpaceCraftingLua();` once:
     - Creates `theLua`.
     - Registers Lua C functions (Material, Node, Girder, Rope, Ring, Slider, etc.).
     - Connects `theSpaceCraft` with `workshop`.

4. **Load Lua ship and build SpaceCraft**
   - Clear previous content:
     - `theSpaceCraft->clear();`
   - Execute Lua file:
     - `if (Lua::dofile(theLua, shipPath)) { printf("ERROR ..."); return 1; }`.
   - Sanity checks / logging (depending on `verbosity`):
     - `theSpaceCraft->checkIntegrity();`
     - Optionally: counts of nodes, girders, ropes, rings, sliders, radiators, shields, tanks, thrusters, welds.

5. **Build mesh from spacecraft**
   - Create builder:
     - `Mesh::Builder2 mesh; mesh.clear();`
   - Call block‑based construction:
     - `BuildCraft_blocks(mesh, *theSpaceCraft, 30.0);` (or configurable `max_size` if needed later).
   - Optionally dump sizes / connectivity:
     - `mesh.printSizes();`
     - `checkSpaceCraftMesh(mesh, *theSpaceCraft);` (non‑fatal; prints warnings about unconnected vertices).

6. **Write OBJ**
   - `mesh.write_obj(objPath);` with default mask (verts + normals + uvs + tris + edges as defined in `ObjMask`).
   - On success, print a one‑line summary if `verbosity >= 1`:
     - e.g. `Wrote OBJ '...': nverts=..., nedges=..., ntris=...`.

7. **Exit code**
   - `0` on success.
   - Non‑zero on fatal issues (missing args, Lua file error, hard `_assert` failures in debug builds, failed write, etc.).

## CMake integration plan

- **File**: `cpp/apps/OrbitalWar/CMakeLists.txt`.
- Inside existing `if( WITH_LUA )` block, after other executables, add:

  ```cmake
  add_executable( spaceCraftMeshExport
      spaceCraftMeshExport.cpp
      $<TARGET_OBJECTS:MeshBuilder2>
  )
  target_link_libraries( spaceCraftMeshExport ${LUA_LIBRARY} )
  ```

- Reuse `include_directories` and Lua configuration already set at the top of the `WITH_LUA` block:
  - `${LUA_INCLUDE_DIR}`
  - `${COMMON_SDL_SRCS}/Lua` (for `LuaHelpers.h`).

## Bash test script: OBJ + optional truss export

- **File**: `tests_bash/Orbital/spaceCraftMeshExport.sh`.
- Responsibilities:
  - Set up symlinks to `cpp/apps/OrbitalWar/data` and `cpp/common_resources` (same pattern as `spaceCraftDynamics.sh`).
  - Rebuild `spaceCraftMeshExport` target using all CPUs‑1 via `make -j` in `cpp/Build/apps/OrbitalWar`.
  - Create a local symlink `spaceCraftMeshExport.x` to the built binary.
  - Run the tool on a chosen Lua ship file and capture stdout.

- Current default run mode (OBJ + `.truss`):

  ```bash
  ./$name.x -s data/ship_ICF_marksman_2.lua \
      -o ship_ICF_marksman_2.obj \
      -t ship_ICF_marksman_2.truss \
      -v 1 | tee OUT-spaceCraftMeshExport
  ```

- This mirrors the workflow of `tests_bash/Orbital/spaceCraftDynamics.sh` but **without** any SDL window, OpenGL context, or truss dynamics.

## Testing & debugging guidelines

- **Rebuild before each run**
  - Always use the bash script so the app is recompiled before execution, ensuring tests run with the latest C++ changes.

- **Use `verbosity` for diagnostics**
  - `verbosity=0`: use when comparing OBJ files or logging noise is undesirable.
  - `verbosity>=1`: show pipeline stages, mesh sizes, and basic stats.
  - `verbosity>=2`: enable additional integrity prints from Lua builders and SpaceCraft components (many Lua/C++ wrappers already gate prints on `verbosity`).

- **Compare OBJ and logs across changes**
  - Use `diff` or scripted comparisons on generated `.obj` files for regression testing of spacecraft geometry.
  - Use `OUT-spaceCraftMeshExport` logs (via `tee`) as additional regression signal.

-- **Future extensions (not required now)**
  - Allow multiple `-s` inputs and batch export.
  - Add flag to choose between `BuildCraft_truss` and `BuildCraft_blocks`.
  - Optionally export additional debug formats (e.g. JSON summaries of component counts).

## Truss export/import: final decoupled design

### On the export side (mesh / spacecraft pipeline only)

- **Header**: `cpp/common/Orbital/SpaceCraft2Mesh_blocks.h`.
- Function:

  ```cpp
  inline void exportSimToFile(
      const char* fname,
      const Mesh::Builder2& mesh,
      const SpaceCrafting::SpaceCraftWorkshop& shop
  );
  ```

- Key properties:
  - Lives next to `BuildCraft_blocks` and related mesh-building helpers.
  - **Does not** depend on `TrussDynamics_*` or any simulation state.
  - Works directly from:
    - `mesh.verts[i].pos` (geometry)
    - `mesh.edges[ib]` (topology + material index `e.w`)
    - `shop.stickMaterials.vec[e.w]` (`preStrain`, `linearDensity`, `Kpush`, `Kpull`, `damping`).
  - Derives per-point masses and per-stick spring parameters and writes a text file in the **`TrussSim v1`** format:
    - Header:
      - `# TrussSim v1`
      - `meta <npoint> <nbond>`
    - Points:
      - `p <id>  <x> <y> <z>  <mass> <nneigh>`
    - Sticks:
      - `s <id>  <i> <j>  <kpull> <kpush> <l0> <damping>`

- This exporter is used by `spaceCraftMeshExport` when `-t <out.truss>` is provided.

### On the import side (simulation only)

- **Header**: `cpp/common/dynamics/TrussDynamics_d.h`.
- **Implementation**: `cpp/common/dynamics/TrussDynamics_d.cpp`.
- Free function:

  ```cpp
  void loadSimFromFile(const char* fname, TrussDynamics_d& sim);
  ```

- Key properties:
  - Depends only on `TrussDynamics_d` and the standard library.
  - Reads the same `TrussSim v1` format produced by `exportSimToFile`.
  - Workflow:
    - First pass: read `meta` to get `npoint`, `nbond` and compute `nneighmax` from all `s` lines.
    - Allocate storage via `sim.recalloc(np, nneighmax, nb)`. 
    - Second pass: fill `sim.points`, `sim.bonds`, `sim.bparams` and neighbor lists; reset forces/velocities.
  - Has **no knowledge of spacecraft or meshes**.

This split keeps the exporter close to the mesh/material information, while the importer lives entirely in the dynamics module, matching the preference that **mesh and simulation stay decoupled** and I/O helpers do not create new cross‑dependencies.

## New headless truss simulation CLI: `trussSimBatch`

### Purpose

- Minimal, non‑interactive tool that:
  - Loads a `.truss` file (in `TrussSim v1` format).
  - Initializes a `TrussDynamics_d` simulation via `loadSimFromFile`.
  - Advances the simulation for a user‑specified number of time steps.
  - Writes a trajectory in XYZ format for offline visualization/analysis.
  - Has **no dependency** on Lua, `SpaceCraft`, `Mesh::Builder2`, or any rendering.

### Implementation

- **Source**: `cpp/apps/OrbitalWar/trussSimBatch.cpp`.
- **CMake target** (in `cpp/apps/OrbitalWar/CMakeLists.txt`):

  ```cmake
  add_executable( trussSimBatch
      trussSimBatch.cpp
      $<TARGET_OBJECTS:TrussDynamics_d>
  )
  target_link_libraries( trussSimBatch )
  ```

- CLI usage:

  ```bash
  trussSimBatch \
    -i ship_ICF_marksman_2.truss \
    -o ship_ICF_marksman_2_traj.xyz \
    [-n steps] [-dt dt] [-psave k] [-m method] [-v verbosity]
  ```

- Options:
  - **`-i <truss.txt>`** (required)
    - Path to truss file written by `exportSimToFile`.
  - **`-o <traj.xyz>`** (optional, default: `traj.xyz`)
    - Output XYZ file: one frame every `-psave` steps.
  - **`-n <steps>`** (optional, default: `1000`)
    - Number of integration steps.
  - **`-dt <dt>`** (optional, default: `5e-4`)
    - Time step size.
  - **`-psave <k>`** (optional, default: `10`)
    - Save one XYZ frame every `k` simulation steps.
  - **`-m <method>`** (optional, reserved)
    - Integration/solver method selector (hook for future extensions; currently always uses `sim.run(1, dt, damp)`).
  - **`-v <int>`** (optional)
    - Sets global `verbosity` for diagnostics.

- Output format (per frame): standard XYZ

  ```text
  N
  step <istep> time <time>
  x y z
  ...
  ```

### Bash wrapper for batch simulation

- **File**: `tests_bash/Orbital/trussSimBatch.sh`.
- Responsibilities:
  - Rebuild `trussSimBatch` with `make -j` in `cpp/Build/apps/OrbitalWar` (ASan build in the current setup).
  - Create local symlink `trussSimBatch.x` to the built binary.
  - Run the batch simulator on a specified `.truss` file and capture stdout:

  ```bash
  #!/bin/bash

  name=trussSimBatch
  dir=../../cpp/Build/apps/OrbitalWar

  # build
  ncpu=`nproc`
  ncpu=$(($ncpu - 1))
  if [ $ncpu -lt 1 ]; then ncpu=1; fi

  make -C ../../cpp/Build/apps/OrbitalWar -j$ncpu $name || exit 1

  ln -sf $dir/$name ./$name.x

  # default input/output (user can edit as needed)
  IN_TRUSS=${1:-ship_ICF_marksman_2.truss}
  OUT_TRAJ=${2:-ship_ICF_marksman_2_traj.xyz}

  ./$name.x -i "$IN_TRUSS" -o "$OUT_TRAJ" -n 1000 -dt 5e-4 -psave 10 -v 1 | tee OUT-trussSimBatch
  ```

This mirrors the `spaceCraftMeshExport.sh` pattern: always rebuild before running, run one canonical test case by default, and tee logs for regression checks.

## End‑to‑end workflow: Lua → OBJ + TRUSS → headless dynamics → XYZ

1. **Generate mesh and truss from Lua**
   - Use `spaceCraftMeshExport` via the bash script:

     ```bash
     cd tests_bash/Orbital
     ./spaceCraftMeshExport.sh
     ```

   - This will (by default):
     - Build the app.
     - Run:

       ```bash
       ./spaceCraftMeshExport.x -s data/ship_ICF_marksman_2.lua \
           -o ship_ICF_marksman_2.obj \
           -t ship_ICF_marksman_2.truss \
           -v 1 | tee OUT-spaceCraftMeshExport
       ```

     - Produce:
       - `ship_ICF_marksman_2.obj` – mesh geometry.
       - `ship_ICF_marksman_2.truss` – truss parameters for simulation.

2. **Run truss dynamics in batch mode**
   - Use the batch simulator script:

     ```bash
     ./trussSimBatch.sh ship_ICF_marksman_2.truss ship_ICF_marksman_2_traj.xyz
     ```

   - This rebuilds `trussSimBatch` if needed, then runs headless dynamics and writes `ship_ICF_marksman_2_traj.xyz`.

3. **Post‑processing / visualization**
   - The `.xyz` file can be visualized with external tools (e.g. VMD, Ovito, custom Python scripts) without any dependency on the original Lua or mesh code.

## Design preferences and takeaways from this setup

- **Strict decoupling of responsibilities**
  - Export of truss parameters lives with the mesh/spacecraft code (`SpaceCraft2Mesh_blocks.h`) and **does not** pull in `TrussDynamics_*`.
  - Import of truss parameters lives with the dynamics (`TrussDynamics_d.{h,cpp}`) and **does not** know about meshes or spacecraft.
  - The batch simulator (`trussSimBatch`) depends only on `TrussDynamics_d` + the import helper.

- **No monolithic `TrussIO` module**
  - Export and import helpers are intentionally not grouped into one header that would create new cross‑dependencies.
  - Each side owns only the pieces it actually needs.

- **Headless, script‑friendly tools**
  - Both CLIs avoid SDL/OpenGL and GUI dependencies, making them suitable for:
    - Automated tests.
    - Batch experiments.
    - Running on remote machines or CI without graphics.

- **Consistent CLI style**
  - Uses existing `LambdaDict`/`process_args` pattern for options.
  - Global `verbosity` from `globals.h` drives all diagnostics.
  - Bash wrappers always **rebuild before run** and use `tee` for log capture.

- **File formats chosen for simplicity and debuggability**
  - Human‑readable `.truss` and `.xyz` files make it easy to:
    - Diff results between runs.
    - Inspect and debug simulations and truss generation.

These conventions should be followed for future tools in this ecosystem: keep modules loosely coupled, favor small, composable CLIs, and prioritize formats and logging that are easy to inspect and regress‑test.

