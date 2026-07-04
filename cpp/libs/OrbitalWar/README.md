# OrbitalWar

Shared library wrapping the space tactics simulation world for Python/ctypes access. Provides C API (`extern "C"`) for initializing `SpaceWorld`, registering body trajectory buffers, running orbital dynamics, and retrieving simulation state. Uses singleton pattern with buffer registry mirroring the LandCraftLib approach.

## Files

- **SpaceTacticsLib.cpp** — `SpaceTacticsLib` library: singleton `SpaceWorld` instance, buffer registry for numpy/ctypes, body buffer registration, trajectory export, ASan-aware debugging
- **CMakeLists.txt** — build target: `SpaceTacticsLib` (SHARED), links `Body` object library

## Related Documentation

- [docs/SpaceTactics/SpaceTactics.md](../../docs/SpaceTactics/SpaceTactics.md) — space tactics game design overview
- [docs/SpaceTactics/SpaceWorld.md](../../docs/SpaceTactics/SpaceWorld.md) — space world model
- [docs/SpaceTactics/SpaceBodies.md](../../docs/SpaceTactics/SpaceBodies.md) — celestial body definitions
- [docs/SpaceTactics/SplineManager.md](../../docs/SpaceTactics/SplineManager.md) — trajectory spline management
- [doc/TopicalAudit/python-cpp-bindings.md](../../doc/TopicalAudit/python-cpp-bindings.md) — ctypes binding patterns
