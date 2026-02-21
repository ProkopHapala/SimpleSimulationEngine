## AeroCombat code map

Summary of `.cpp`/`.h` files under `cpp/apps/AeroCombat`, with notes on how they use core dynamics classes (`AeroSurf`, `AeroCraft`, `AeroCraftControl`, `Body`).

### Top-level runtime / GUI
- **AeroCombat_main.cpp**  
  SDL2/OpenGL 1.x flight sandbox. Own `AeroCraftGUI` (AppSDL2OGL_3D) that loads an `AeroCraftWarrior`, registers it with a `Shooter` world, and wires both manual input (`AeroControler1`) and autopilot (`AeroCraftControler`). Uses `AeroSurf` panel debug records for HUD plots and logging. Physics comes from `Body` via `AeroCraft` inheritance; control surfaces updated through `AeroCraftControl` helpers.
- **AeroCombatOGL3.cpp**  
  Modern OpenGL 3.x variant of the GUI (`AeroCraftGUI` : AppSDL2OGL3 + SceneOGL3). Creates terrain, shaders, and meshes, spawns `AeroCraftWarrior` plus backup `AeroCraft`, attaches `AeroCraftControler` autopilot and `AeroControler1` manual controller, and steps the world (`Shooter`). Uses `AeroSurf` debug records for left/right ailerons. Relies on `Body` integration through `AeroCraft`.
- **AeroCraftGUI.h / AeroCraftGUI.cpp**  
  Older GUI wrapper around `AeroCraftWorld`/`Shooter`. Provides camera, skybox, and HUD drawing. Renders the active `AeroCraft` using its `AeroSurf` panels; applies control via mouse/keyboard, optionally autopilot (`autoPilot` flag). Physics update is delegated to `world->update`, which ultimately uses `Body` motion on the craft.

### Editors / design tools
- **AeroCraft_editor.cpp**  
  Legacy OpenGL editor for visualizing aerodynamic forces. Instantiates `AeroCraftWarrior`, manipulates AoA/side-slip, calls `AeroCraft::applyAeroForces`, and draws per-panel `AeroSurf` forces via `AeroSurfaceDebugRecord`. `Body` motion is not stepped; focus is on force visualization.
- **AeroCraft_editor2.cpp**  
  Extended editor adding potential-flow (`PotentialFlowSystem`) evaluation via `AeroCraftDesign`. Still builds an `AeroCraftWarrior`, applies aero forces, and visualizes flow fields/forces. Uses `AeroSurf` panels and debug records; relies on `AeroCraft`/`Body` for state but runs mainly static evaluations.
- **AeroCraftDesign.h**  
  Geometry/weight aggregation for aircraft designs (fuselages, wings, gun slots). Converts wing sections into a vortex-lattice `PotentialFlowSystem` (panels/vortices). Uses `Mat3d` transforms; independent of runtime physics but meant to feed `AeroCraft`-like models and their `AeroSurf` discretization.

### Control / flight logic
- **AeroControler1.h**  
  Simple PID-like controller operating on an `AeroCraft` and a backup baseline craft. Computes roll/pitch/yaw corrections and applies them to control surfaces (left/right aileron, rudder, elevator) by rotating their `AeroSurf` local frames. Uses current craft orientation (from `Body`/`AeroCraft`) to chase desired up/forward vectors.
- **AeroCombatHelpers.h**  
  Utility functions to evaluate aerodynamic forces at different orientations (`AeroSurf` panels) and build terrain display lists. `evalAeroFoceRotated/AtRotations` call `AeroCraft::applyAeroForces` and map forces/torques back to world frame. Terrain/skybox helpers are graphics-only. Depends on `AeroCraft` (`Body`-based motion) and `AeroSurf` force computation.
- **AeroTest.h**  
  `AeroTester` collects trajectories for an `AeroCraft` driven by an `AeroCraftControler` autopilot. Reallocates histories, advances the craft by repeatedly calling `applyAeroForces` and `move` (from `Body`), and logs positions/forces/attitudes.

### Rendering helpers
- **AeroDraw.h**  
  Immediate-mode drawing helpers for an `AeroCraft`: renders `AeroSurf` panels, orientations, and per-panel force/velocity debug vectors. Relies on the craft’s `rotMat`/`pos` from `Body` state.
- **AeroCraftWarrior.h**  
  Concrete aircraft used in combat sims; inherits `AeroCraft` and `Warrior3D`. Implements `move_warrior` by accumulating gravity, applying aerodynamic forces (`AeroSurf` panels), then stepping `Body` integration.

### Notes
- Temp files under `temp/` (AeroMath, AeroStaticTest, AeroSurf, ProkopMath) look experimental and are not wired into the main builds; they follow similar patterns (aero math, surface tests).
- All runtime simulators (`AeroCombat_main`, `AeroCombatOGL3`, `AeroCraftGUI`, editors) depend on `AeroSurf` for surface discretization, `AeroCraft` for vehicle state (inherits `Body`), and `AeroCraftControl`/`AeroControler1` for control surface actuation.
