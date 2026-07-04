# LandTactics

Tactical simulation of modern land battle (WWII-era). Units use terrain cover, simulate mobility across different terrain types, and compute line-of-sight in 2D top-down view. Squads of soldiers with different unit types (infantry, vehicles) engage in combat with realistic damage models.

## Simulation Features

- 2D top-down tactical view with terrain rendering
- Line-of-sight computation considering terrain elevation
- Unit mobility varies by terrain type (forest, road, open ground)
- Squad-level organization with unit types and stats
- Combat resolution with armor, penetration, and damage models
- Faction-level AI for unit orders and tactical behavior
- Shelters and cover mechanics

## Files

- **LandTactics_main.cpp** — main application: tactical view, unit control, combat
- **LTWorld.cpp / .h** — world model: terrain, line-of-sight, combat resolution, factions
- **LTSquad.cpp / .h** — squad: group of units, formation, orders, movement
- **LTUnit.cpp / .h** — individual unit: position, health, combat actions
- **LTUnitType.cpp / .h** — unit type definitions: speed, armor, weapons, terrain modifiers
- **LTFaction.cpp / .h** — faction: squads, AI strategy, victory conditions
- **LTSurroundings.cpp / .h** — terrain analysis: cover, line-of-sight, mobility
- **LTShelter.h** — shelter/cover objects (buildings, trenches, vegetation)
- **LTdraw.h** — rendering: units, terrain, line-of-sight overlays, UI
- **LTrender.h** — additional rendering helpers
- **LTcommon.h** — shared types and constants
- **CMakeLists.txt** — build target: `LandTactics_main`
- **data/** — unit type definitions, terrain data, scenarios
- **notes/** — design notes

## Related Documentation

- [docs/LandTactics/LandTactics.md](../../docs/LandTactics/LandTactics.md) — LandTactics overview and architecture
- [docs/LandTactics/LTWorld.md](../../docs/LandTactics/LTWorld.md) — world model: terrain, line-of-sight, combat
- [docs/LandTactics/LTUnit.md](../../docs/LandTactics/LTUnit.md) — unit model: position, health, combat actions
- [docs/LandTactics/LTUnitType.md](../../docs/LandTactics/LTUnitType.md) — unit type definitions: speed, armor, weapons
- [docs/LandTactics/LTFaction.md](../../docs/LandTactics/LTFaction.md) — faction model: squads, AI strategy
- [docs/LandTactics/LTShelter.md](../../docs/LandTactics/LTShelter.md) — shelter and cover mechanics
- [docs/LandTactics/LTcommon.md](../../docs/LandTactics/LTcommon.md) — shared types and constants
- [docs/LandTactics/AirCombatModel.md](../../docs/LandTactics/AirCombatModel.md) — air combat model for combined-arms scenarios
