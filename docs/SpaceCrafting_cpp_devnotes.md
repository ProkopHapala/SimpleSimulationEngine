# SpaceCrafting C++ Dev Notes

These notes summarize recent issues in the C++ SpaceCrafting tools (editors + dynamics app), how they were fixed, and what to watch for in future work.

---

## 1. Slider paths not initialized ("Slider path not initialized! Skipping.")

**Symptom**
- When loading a Lua spacecraft and running the dynamics/editor apps, `sliders2edgeverts()` printed:
  - `Slider path not initialized! Skipping.`
- Slider wheels were not bound into the truss simulation, so wheel controls had no effect.

**Root cause**
- `Slider::path` (the path along which the slider moves) was never populated before calling `sliders2edgeverts()`.
- The required call `updateSlidersPaths()` on structural components (e.g. `Ring`) was missing during reload.

**Fix**
- In the editor and dynamics code paths, after constructing the truss from the spacecraft:
  - Call `updateSlidersPaths(true, true, sim.points)` on each ring/component holding sliders.
  - Only then call `sliders2edgeverts(*theSpaceCraft, sim)`.

**Takeaways**
- Whenever we create or reload a `SpaceCraft` and then build a `TrussDynamics_*` from it, the sequence must be:
  1. Build mesh/truss from `SpaceCraft`.
  2. Call `updateSlidersPaths()` on all relevant components.
  3. Call `sliders2edgeverts()` to bind sliders to edges/points.
- Document this as an invariant near `sliders2edgeverts()` and in any reload/init helpers.

---

## 2. Lua ship loading overwritten by default shape

**Symptom**
- Running `SpaceCraftEditorNew` with `-s data/ship_ICF_marksman_2.lua` produced no visible ship or overrode it with a default.

**Root cause**
- `SpaceCraftEditorNew::reloadShip()` used `simulator->initSimDefault()` which creates a generic shape, discarding the loaded Lua spacecraft.
- Command-line parsing for `-s` was not wired in `SpaceCraftEditorNew` the same way as in other apps.

**Fix**
- Wire `process_args` / `LambdaDict` in `SpaceCraftEditorNew::main()` to handle `-s` and pass the Lua file to `reloadShip`.
- Change `reloadShip()` to:
  - Call `simulator->reloadShip(fname)` to load Lua.
  - Call `simulator->initSimulators()` (which exports the loaded spacecraft into `TrussDynamics_d`/`_f`).
  - Then run the slider-path/sliders2edgeverts sequence.

**Takeaways**
- Avoid mixing "default-test" initialization (`initSimDefault`) with real content loading.
- After Lua reloads, always re-run the full init chain that builds the simulation from the loaded ship.

---

## 3. Picker crash in `PickerUI::pick()` (null/invalid picker & ray)

**Symptom**
- ASan SEGV when clicking in `SpaceCraftEditorNew`.
- Crash inside `PickerUI::pick()` / `picker->pick_nearest()`.

**Root cause**
- `picker.picker` backend not guaranteed to be initialized before use in the new editor.
- Ray origin/direction (`ray0`, `hray`) were not set from the camera before calling `pick()`.

**Fix**
- In `SpaceCraftEditorNew` constructor:
  - After `bindSimulators(simulator)`, assign `picker.picker = _sim` (once `_sim` is valid).
- In mouse event handling:
  - Set `picker.ray0` and `picker.hray` from `cam.pos` and camera direction before picking.
  - Only call `picker.pick()` if `picker.picker` is non-null.

**Takeaways**
- For any GUI picking code, treat the backend picker pointer and ray parameters as required preconditions.
- Guard calls with null-checks and make the camera → ray mapping explicit.

---

## 4. Heap-buffer-overflow in `plotSliders()`

**Symptom**
- ASan reported a heap-buffer-overflow on read inside `SpaceCraftDynamicsApp::plotSliders()`.
- Backtrace pointed to indexing `sim.points[...]` while drawing slider bonds.

**Root cause**
- `plotSliders()` assumed:
  - `craft.sliders.size()` matches `sim.nEdgeVertBonds`.
  - All slider `pointRange.x` and `EdgeVertBond.verts` indices are valid in `[0, sim.nPoint)`.
- When this was not true (e.g. some sliders had invalid or out-of-range indices), we read past `sim.points` or `sim.edgeVertBonds`.

**Fix**
- Harden `plotSliders()`:
  - Early return if `craft.sliders.empty()`, `sim.nEdgeVertBonds <= 0`, or `sim.edgeVertBonds == nullptr`.
  - Use `n = min(craft.sliders.size(), sim.nEdgeVertBonds)` for the main loop.
  - For each slider, validate:
    - `ip = o->pointRange.x` in `[0, sim.nPoint)`.
    - `iv0 = ev.verts.x`, `iv1 = ev.verts.y` both in `[0, sim.nPoint)`.
  - If any index is invalid:
    - If `verbosity > 0`, print a descriptive error line.
    - If `exit_on_error` is true, exit immediately to make the bug obvious in debugging builds.
    - Otherwise `continue` and skip drawing that slider.

**Takeaways**
- Drawing/visualization functions must not assume perfect consistency from complex upstream geometry/simulation builders.
- Always bound-check any index coming from data structures that can be out-of-sync (e.g. when reloading a ship or partially rebuilding a simulation).
- Use `verbosity` and `exit_on_error` from `globals.h` to turn silent data bugs into clear diagnostics during development.

---

## 5. Wheel dynamics control (numpad) not affecting wheels in `SpaceCraftEditorNew`

**Symptom**
- In the legacy editor and `spaceCraftDynamics` app, numpad keys controlled wheels/sliders.
- In `SpaceCraftEditorNew`, the same keys appeared to do nothing.

**Root cause**
- `SpaceCraftDynamicsApp::keyStateHandling()` was already updating `simulator->wheel_speed` from numpad state.
- However, the simulation in `SpaceCraftEditorNew` never applied `wheel_speed` to sliders each step:
  - No `sim.user_update` callback was configured.
  - `applySliders2sim()` was not called from within the simulation step.

**Fix**
- Introduce a global pointer to the editor simulator and a control callback, e.g.:

  ```cpp
  static SpaceCraftSimulator* gEditorSimulator = nullptr;

  void SpaceCraftControl(double dt){
      if(gEditorSimulator){
          applySliders2sim( *theSpaceCraft, gEditorSimulator->sim,
                            (double*)&gEditorSimulator->wheel_speed );
      }
  }
  ```

- In `SpaceCraftEditorNew` constructor:
  - After creating the simulator and binding it, set `gEditorSimulator = simulator`.
- In `SpaceCraftEditorNew::reloadShip()` after `initSimulators()` and slider setup:
  - Set `simulator->sim.user_update = SpaceCraftControl;`.

**Takeaways**
- For interactive control, it is not enough to update input state (`wheel_speed`); the simulation loop must call a user update hook (`sim.user_update`) every step.
- Make the control pipeline explicit:
  - Input (keyboard) → high-level control variables (`wheel_speed`) → user update callback → apply to `TrussDynamics` / sliders.

---

## 6. Using `globals.h` diagnostics (`verbosity`, `exit_on_error`, `_assert`)

**Pattern**
- `globals.h` provides:
  - `verbosity` (int) to control how much gets printed.
  - `exit_on_error` (bool) to decide whether to terminate on errors.
  - Macros and helpers like `_assert` and `DBG(...)` for richer debugging.

**How we used it**
- In `plotSliders()`, invalid indices now:
  - Print a clear error message when `verbosity > 0`.
  - Optionally `exit(0)` when `exit_on_error` is true.

**Takeaways**
- Prefer using these globals instead of ad-hoc `printf`s spread across the code.
- For new debug checks:
  - Keep messages short but specific (which structure, which index, what range).
  - Respect `verbosity` so that release/normal runs stay quiet.
  - In debug builds, consider enabling `exit_on_error` so mistakes fail fast.

---

## General Lessons

- **Define and respect initialization order**
  - For complex objects (SpaceCraft → Mesh → Truss → Sliders), write down the required init sequence and reuse it.
  - Avoid having multiple half-overlapping init functions that may or may not run (e.g. `initSimDefault` vs `initSimulators`).

- **Trust, but verify indices**
  - Any time we convert higher-level objects into indexed arrays (vertices, edges, sliders), follow up with bound checks where data is consumed (rendering, control).

- **Centralize control hooks**
  - Use `sim.user_update` (and similar hooks) as the single place where input state influences the simulation each step.
  - This avoids scattered calls to `applySliders2sim()` and makes behavior easier to understand.

- **Use diagnostics aggressively during development**
  - Leverage `verbosity`, `exit_on_error`, and `_assert` to make data inconsistencies fail loudly while iterating.
  - Once stable, lower `verbosity` and disable fatal exits as needed for production-like runs.

---

## 7. Slider anchor indices: `ivert` vs `pointRange`

**Observation**
- `Slider` inherits `pointRange` from `ShipComponent`, but in the current mesh/export pipeline sliders do **not** own a contiguous vertex range in the truss.
- Instead, `Slider::ivert` is the meaningful anchor into `sim.points`.
- As a result, `pointRange.x` for sliders is often left at the default `-1`.

**Bug**
- The first version of `SpaceCraftDynamicsApp::plotSliders()` used `slider->pointRange.x` to look up the slider point, which for sliders was `-1`, leading to invalid indices and ASan errors.

**Fix**
- Change `plotSliders()` to use `slider->ivert` as the slider anchor vertex:
  - Validate `ivert` in `[0, sim.nPoint)` and log/abort via `globals.h` if out of range.
  - Keep using `EdgeVertBond.verts.{x,y}` for the two end vertices of the edge.

**Takeaway**
- For sliders, treat `ivert` as the canonical index into the truss; `pointRange` is for components that actually map to a contiguous vertex segment (girders, rings, ropes, etc.).
- When adding new drawing or control code, always confirm which index field is semantically valid for a given component type before using it.
