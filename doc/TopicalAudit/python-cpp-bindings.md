---
type: TopicalAudit
title: Python-C/C++ Bindings
tags: [topic, python, cpp, ctypes, bindings, c-api, numpy, shared-library]
---

## Summary

Python ↔ C/C++ interop via ctypes with automatic interface generation. C++ libraries expose flat C API functions; Python loads shared libraries via ctypes, defines function signatures, and wraps them with numpy buffer views. A code generation utility (`cpp_utils.py`) parses C function headers and emits Python ctypes boilerplate. The primary use case is LandCraftLib (terrain/world simulation), with patterns reusable for other C++ libraries.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/libs/CombatModels/LandCraftLib.cpp` | active | C API wrapper around `LandCraftWorld` class. Exposes init, buffer access, terrain, hydraulics, roads, vehicles, economy, pathfinding. Global buffers for numpy/ctypes interop. |
| Python | `python/LandCraft/landcraft.py` | active | ctypes bindings to LandCraftLib. Loads `.so`, defines function signatures, wraps with Python functions. numpy buffer views for ground/water arrays. |
| Python | `python/pyMeta/cpp_utils.py` | active | Automatic ctypes interface generator. Parses C function headers, translates types, emits `argtypes`/`restype`/wrapper code. Also handles library compilation and loading. |
| Python | `python/pySimE/` | active | Python simulation engine wrapping C++ core via ctypes. |
| Doc | `doc/CodingRules/workflows/python/python_ctypes_workflow.md` | doc | Workflow documentation for ctypes bindings. |
| Doc | `doc/AGENTs/skills/ctypes-bindings/SKILL.md` | doc | Skill definition: auto-build, zero-copy buffers, memory layout. |
| Test | `tests_bash/LandCraft/test_landcraft.py` | active | Integration test: init buffers, retrieve numpy views, save PGM images, plot hydraulics. |

## Sub-topics

### LandCraftLib C API

`LandCraftLib.cpp` exposes a flat C interface:
- `landcraft_init()` — initialize world with grid dimensions
- `landcraft_getGround()`, `landcraft_getWater()` — return raw pointers to global buffers
- `landcraft_genTerrain()`, `landcraft_stepHydraulics()` — simulation steps
- `landcraft_buildRoad()`, `landcraft_addVehicle()` — world modification
- `landcraft_stepEconomy()`, `landcraft_findPath()` — economy and pathfinding
- Global static buffers managed in C++, Python obtains numpy views via `ctypes.data_as()`

### Python ctypes Bindings

`landcraft.py` pattern:
1. Load shared library: `ctypes.CDLL(BUILD_PATH + "/libLandCraftLib.so")`
2. Define function signatures: `lib.func.argtypes = [c_int, c_float_p, ...]`, `lib.func.restype = c_void_p`
3. Wrap each function with Python convenience method
4. Use `numpy.ctypeslib.as_array()` or `arr.ctypes.data_as()` for zero-copy buffer access

### Automatic Interface Generation

`cpp_utils.py` provides:
- `parseFuncHeader(s)`: Parses C function signature string → (name, ret_type, arg_types, arg_names)
- `translateTypeName(t)`: Maps C types to ctypes types (`int*` → `c_int_p`, `float` → `c_float`)
- `writeFuncInterface(parsed)`: Emits Python code: `argtypes` assignment, `restype` assignment, wrapper function
- `writeFuncInterfaces(headers)`: Batch code generation from list of C headers
- `compile_lib(name)`: Compiles C++ source to `.so` via `g++` with `-fPIC -shared`
- `loadLib(name, recompile)`: Compiles (optional) and loads via ctypes
- `make(what)`: Runs `make` in build directory

### Build Integration

- `cpp_utils.py` expects shared libraries at `cpp/Build/libs/CombatModels/`
- `BUILD_PATH` computed relative to `cpp_utils.py` location: `../../../cpp/Build/libs/CombatModels`
- Compilation flags: `-std=c++11 -O3 -ftree-vectorize -unroll-loops -ffast-math`
- Optional `CPP_BUILD_PATH` environment variable override (per `cpp-build` skill)
- `clean_build` and `recompile_glob` flags control automatic recompilation on import

## Parity Status

- **LandCraftLib.cpp ↔ landcraft.py**: Python bindings mirror C API 1:1. Function signatures manually defined. No automatic sync — if C API changes, Python must be updated manually.
- **cpp_utils.py code generator ↔ actual bindings**: `landcraft.py` uses manually defined signatures, not auto-generated. Code generator exists but is not actively used for LandCraft.
- **pySimE ↔ C++ core**: Separate binding, not audited in detail.

## Open Issues

- No automatic synchronization between C API and Python bindings — manual maintenance required
- `cpp_utils.py` code generator exists but `landcraft.py` doesn't use it — duplication of effort
- Global static buffers in `LandCraftLib.cpp` are not thread-safe
- `compile_lib()` uses `os.system()` — no error handling for compilation failures
- `try/except: pass` in `compile_lib()` cleanup silently swallows errors (violates fail-fast principle)
- No memory ownership documentation — unclear who owns buffers returned to Python
- `BUILD_PATH` is hardcoded relative to `cpp_utils.py` — fragile if package structure changes
- No ctypes validation tests (verifying function signatures match actual C functions)
