# Tanks

Tank simulation game with turret control, vehicle movement on terrain, and armor penetration calculation. Terrain is a rectangular grid using cubic splines. Each tank has armor plates with specific thickness, and guns can be equipped in fixed ports.

## Features

- Tank movement on cubic-spline terrain
- Turret rotation and gun aiming
- Armor penetration model: projectile vs. armor thickness per polygon
- Multiple gun ports with customizable loadout
- 2D side-view rendering

## TODO (from source)

- Terrain → rectangular grid using cubic splines
- Armor → each mesh polygon has particular thickness
- Guns → fixed number of ports, which can be equipped with different guns and repositioned

## Files

- **Tanks_main.cpp** — main application: tank control, terrain, combat, rendering
- **Tank.cpp / .h** — tank model: hull, turret, guns, armor, movement, aiming
- **TankHelpers.h** — helper functions: projectile-ballistics, armor penetration calculation
- **TODO.md** — development TODO list
- **CMakeLists.txt** — build target: `Tanks_main`
- **data/** — tank configuration data
- **notes/** — design notes (2 files, including ArmorPenetration.md)

## Related Documentation

- [docs/FormationTactics/DamagePhysicsModel.md](../../docs/FormationTactics/DamagePhysicsModel.md) — damage physics model with armor penetration (shared combat model)
