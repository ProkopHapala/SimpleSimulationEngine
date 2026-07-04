# Shock

Shared library for 1D shock wave simulation. Provides C API for Python/ctypes access to `ShockSystem1D` — simulates shock propagation through layered materials with velocity, mass, and volume tracking per cell.

## Files

- **Shock.cpp** — `Shock` library: `ShockSystem1D` singleton, material initialization, Python buffer interface (`init_python`, `getArray`), time stepping
- **CMakeLists.txt** — build target: `Shock` (SHARED, links `ShockWaves` objects)
