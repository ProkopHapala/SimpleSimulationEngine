# FormationTactics

Sophisticated TotalWar-like simulation of historical land battles. Units are organized into formations (battle lines, factions) with individual soldiers that have morale, stamina, and combat stats. Soldiers fight with melee weapons, formations maintain cohesion, and factions respond to tactical orders.

## Simulation Features

- Formation-based unit control: battle lines, columns, skirmish lines
- Individual soldier AI: melee combat, morale, fatigue, routing
- Soldier types with different weapons, armor, and stats
- Faction-level tactical AI
- 2D top-down rendering with unit icons and formation indicators

## Files

- **FormationTactics_main.cpp** — main application: battle setup, rendering, interaction
- **FormationWorld.cpp / .h** — world model: factions, formations, collision, combat resolution
- **Formation.cpp / .h** — formation logic: movement, cohesion, formation shapes
- **Faction.cpp / .h** — faction management: units, AI orders, victory conditions
- **Soldier.cpp / .h** — individual soldier: health, morale, combat, movement
- **SoldierType.cpp / .h** — soldier type definitions: weapons, armor, stats, abilities
- **BattleLine.cpp / .h** — battle line formation (stub)
- **PolyLineFormation.h** — polyline-based formation shape
- **FormationTacticsCommon.h** — shared types and constants
- **CMakeLists.txt** — build target: `FormationTactics_main`
- **Notes/** — design notes
- **data/** — unit type definitions, battle scenarios

## Related Documentation

- [docs/FormationTactics/DamagePhysicsModel.md](../../docs/FormationTactics/DamagePhysicsModel.md) — damage physics model for melee combat (with Python implementation)
- [docs/FormationTactics/PolearmsPhysicsModel.md](../../docs/FormationTactics/PolearmsPhysicsModel.md) — polearms physics model: reach, swing dynamics, hit detection (with Python implementation)
