# MinimalTactics

Minimal tactical combat simulation for testing simplified combat models. A lightweight version of LandTactics/FormationTactics focused on core combat mechanics without complex terrain or formation systems.

## Files

- **MinimalTactics_main.cpp** — main application: unit placement, combat, rendering
- **TacWorld.cpp / .h** — tactical world: units, factions, combat resolution
- **Unit.cpp / .h** — unit: position, health, attack, movement
- **UnitType.cpp / .h** — unit type: speed, attack power, defense
- **Faction.cpp / .h** — faction: unit grouping, basic AI
- **MinimalTacticsCommon.h** — shared types
- **CMakeLists.txt** — build target: `MinimalTactics_main`
- **data/** — unit definitions
