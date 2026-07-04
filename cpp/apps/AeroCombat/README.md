# AeroCombat

Flight simulator with realistic aerodynamics and air combat. Player controls an aircraft with mouse and keyboard, flying over terrain while shooting at targets (icosahedra on sticks). Includes an aircraft model editor for designing the aerodynamic surfaces (wings, fuselage, control surfaces).

## Physics

- Aerodynamic forces computed from wing surface panels using `AeroSurf.h` and `AeroCraft.h`
- Lift, drag, and moments calculated per-surface, aggregated to aircraft-level dynamics
- Rigid body dynamics via `Body.h`
- Terrain collision detection

## Files

- **AeroCombat_main.cpp** — main application: flight sim with shooting, terrain, targets
- **AeroCombatOGL3.cpp** — OpenGL 3+ renderer version of the combat sim
- **AeroCraft_editor.cpp** — interactive editor for aircraft aerodynamic models
- **AeroCraft_editor2.cpp** — newer version of the aircraft editor
- **AeroCraftGUI.cpp / .h** — GUI panel for aircraft parameters (lift/drag visualization, control surfaces)
- **AeroCraftDesign.h** — aircraft design data: wing layout, surface parameters, control surfaces
- **AeroControler1.h** — AI autopilot controller for aircraft
- **AeroCraftWarrior.h** — wrapper combining aircraft + combat behavior
- **AeroCombatHelpers.h** — shared helper functions for combat simulation
- **AeroDraw.h** — rendering helpers for aircraft (surfaces, control surfaces, vectors)
- **AeroTest.h** — test/configuration helpers
- **CMakeLists.txt** — build targets: `AeroCombat_main`, `AeroCraft_editor`, `AeroCraft_editor2`, `AeroCombatOGL3`
- **Notes/** — design notes
- **data/** — aircraft model data
- **python/** — Python analysis scripts

## Related Documentation

- [docs/AeroCombat/AeroCombat.md](../../docs/AeroCombat/AeroCombat.md) — aerodynamic combat model overview
- [docs/AeroCombat/AeroSurf.md](../../docs/AeroCombat/AeroSurf.md) — aerodynamic surface (wing panel) force calculation
