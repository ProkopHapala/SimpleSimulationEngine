# CombatModels

Shared libraries wrapping combat simulation models for Python/ctypes access. Provides C API (`extern "C"`) wrappers around C++ combat, land tactics, land craft, and space combat models. Each library is built as a `.so` and loaded via `ctypes` from Python.

## Libraries

- **libCombatModels.so** — general combat models: weapon types, unit types, unit states, division-level combat (`MinimalDivisionLevel.h`)
- **libSpaceCombatLib.so** — space combat: `CombatAssembly` with projectile types, space guns, Whipple shields, projected targets (`spaceCombat.h`)
- **libLandCombatLib.so** — land combat: LT gun types, unit types, units from `apps/LandTactics/` (links `Body2D` objects)
- **libLandCraftLib.so** — LandCraft world C API: terrain, hydraulics, roads, economy, pathfinding, vehicles (wraps `apps/LandCraft/LandCraftWorld`, links `TerrainGrid2D`, `Noise`)

## Files

- **CombatModels.cpp** — `CombatModels` library: weapon/unit type registries, combat resolution using `MinimalDivisionLevel.h`
- **SpaceCombatLib.cpp** — `SpaceCombatLib` library: `CombatAssembly` with space guns, projectiles, Whipple shields, damage models
- **LandCombatLib.cpp** — `LandCombatLib` library: C API for LT gun/unit types, unit creation, combat (wraps `apps/LandTactics/LTUnit*.h`)
- **LandCraftLib.cpp** — `LandCraftLib` library: C API for LandCraft world (terrain gen, hydraulics, roads, economy, pathfinding), buffer registry for numpy/ctypes
- **CMakeLists.txt** — build targets: `CombatModels`, `SpaceCombatLib`, `LandCombatLib`, `LandCraftLib` (all SHARED)

## Related Documentation

- [docs/SpaceTactics/spaceCombat.md](../../docs/SpaceTactics/spaceCombat.md) — space combat mechanics: weapons, shields, damage
- [docs/FormationTactics/DamagePhysicsModel.md](../../docs/FormationTactics/DamagePhysicsModel.md) — damage physics model with armor penetration
- [docs/LandTactics/LTUnitType.md](../../docs/LandTactics/LTUnitType.md) — land tactics unit type definitions
- [docs/LandCraft/LandCraft.md](../../docs/LandCraft/LandCraft.md) — LandCraft architecture overview
- [doc/TopicalAudit/python-cpp-bindings.md](../../doc/TopicalAudit/python-cpp-bindings.md) — ctypes binding patterns used by these libraries
