# StrategicBattleModels

Experimental strategic-level battle models for air and ground combat. Defines data structures for unit types, weapons, armor penetration, hit boxes, and damage resolution. Work-in-progress — contains header files with type definitions and combat formulas, not a runnable application.

## Models

- **Unit types**: speed, hit boxes, guns, armor
- **Aircraft units**: speed, climb rate, guns
- **Damage model**: armor penetration vs. hit box armor, kill probability, repair cost/time
- **Vehicle function**: product of component health (all vital components must function)
- **Rival factions**: support, spearhead, bombers, fighters with order states

## Files

- **src/Rival.h** — rival faction model: unit brands, order enums, gun/armor/damage types
- **src/Combat.h** — combat resolution interface (stub)
- **src/Battlefield.h** — battlefield conditions interface (stub)
- **CMakeLists.txt** — build target: `MinimalTactics_main` (reuses MinimalTactics code)
