# BlockHouseTactics

Tactical building game where the player constructs blockhouses from boxes stored in a 3D spatial hashmap. Each box has faces that can be selected and manipulated. Experimental prototype for a panel-based housing/construction system.

## Files

- **BlockHouseTactics_main.cpp** — main application: interactive 3D box placement and manipulation
- **BlockHouseWorld.cpp / .h** — world model: 3D hashmap of boxes, face selection, block operations
- **FormationWorld.h** — shared formation/world interface (stub)
- **CMakeLists.txt** — build target: `BlockHouseTactics_main`
