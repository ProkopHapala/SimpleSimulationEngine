# libSail

Shared library for sailing ship hydrostatics and stability. Provides C API for computing 2D hull volume, moment, and static stability from triangulated hull cross-sections. Uses Brent line search to find waterline for a target displacement.

## Files

- **libSail.cpp** — `libSail` library: hull volume calculation, hull moment, static stability analysis, waterline search via `lineSearch.h`
- **CMakeLists.txt** — build target: `libSail` (SHARED)

## Related Documentation

- [cpp/apps/SailWar/README.md](../../apps/SailWar/README.md) — SailWar application using sailing physics
