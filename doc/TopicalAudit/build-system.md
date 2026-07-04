---
type: TopicalAudit
title: Build System
tags: [topic, cpp, cmake, make, gcc, build, asan, opencl, sdl2, opengl, python-ctypes]
---

## Summary

Multi-language build system centered on CMake for C++ (primary), shell scripts for C, and Python ctypes for runtime compilation. C++ build uses CMake with feature toggles (SDL, OpenCL, OpenMP, ASAN, Lua, networking, music), conditional subdirectory inclusion, and OBJECT libraries for code sharing. Build output goes to `cpp/Build/` (symlink to `Build-opt/` or `Build-asan/`). Python bindings use `cpp_utils.py` for runtime `g++` compilation or `make` invocation. `compile_commands.json` exported for clangd.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| CMake | `cpp/CMakeLists.txt` | active | Root C++ build. 179 lines. Options: RELEASE, WITH_SDL, WITH_OMP, WITH_OPENCL, WITH_LUA, WITH_ASAN, WITH_NET, WITH_MUSIC. C++20, `-fPIC`, `-fno-strict-aliasing`. |
| CMake | `cpp/common/CMakeLists.txt` | active | Common library subdirectories: algorithms, math, dynamics, molecular, maps, dataStructures, engine, CombatModels, utils, geometry. |
| CMake | `cpp/common_SDL/CMakeLists.txt` | active | SDL2/OpenGL subdirectories: SDL2, SDL2OGL, SDL2OGL3, network, Lua. Requires OpenGL, GLU, SDL2. |
| CMake | `cpp/libs/CMakeLists.txt` | active | Shared libraries: libFlight, Molecular, Shock, KosmoSuite, CombatModels, OrbitalWar. |
| CMake | `cpp/libs/CombatModels/CMakeLists.txt` | active | `CombatModels`, `SpaceCombatLib`, `LandCombatLib`, `LandCraftLib` shared libraries. Uses OBJECT libraries for cross-library code sharing. |
| CMake | `cpp/apps/CMakeLists.txt` | active | Application executables: 17 app subdirectories (SailWar, AeroCombat, FormationTactics, LandTactics, LandCraft, OrbitalWar, Tanks, MolecularEditor, etc.). |
| CMake | `cpp/apps_OCL/CMakeLists.txt` | active | OpenCL apps: MolecularEditorOCL. Requires OpenCL + SDL + OpenGL. |
| CMake | `cpp/libs/Molecular/CMakeLists.txt` | active | `Molecular`, `RigidMol`, `ReactiveFF`, `eFF_lib` shared libraries. Uses `MolecularEngine` and `DynamicOpt` OBJECT libraries. |
| Shell | `C/build.sh` | active | C build script. `gcc -std=c99 -Og -Wall`. Links GL, SDL2, tcc, dl. Simple single-file compilation. |
| Python | `python/pyMeta/cpp_utils.py` | active | Runtime compilation via `g++` or `make`. Loads `.so` via ctypes. Auto-generates ctypes interface from C headers. |
| Doc | `doc/AGENTs/skills/cpp-build/SKILL.md` | doc | Build skill: Build symlink (opt/asan), cmake targets, ASAN setup, common failures. |

## Sub-topics

### CMake Configuration

`cpp/CMakeLists.txt` root file:

**Build options** (all OFF by default):
- `RELEASE`: `-O3 -Ofast -march=native` vs debug `-g -Og`
- `WITH_SDL`: GUI/3D graphics (SDL2+OpenGL)
- `WITH_OMP`: OpenMP parallelization
- `WITH_OPENCL`: OpenCL GPU acceleration
- `WITH_LUA`: Lua scripting
- `WITH_ASAN`: AddressSanitizer memory debugging
- `WITH_NET`: Networking (multiplayer)
- `WITH_MUSIC`: Music visualizer

**Compiler flags**:
- Global: `-std=c++20 -fPIC -fno-strict-aliasing`
- Debug: `-g -Og` (+ ASAN flags if enabled)
- Release: `-w -Ofast -march=native`
- Warning suppression: extensive `-Wno-*` list
- Warnings-as-errors: return-type, init-self, uninitialized, implicit-fallthrough, tautological-compare, delete-non-virtual-dtor, overloaded-virtual

**Conditional build tree**:
- `common/` always built
- `libs/` always built
- `libs_OCL/` built if `WITH_OPENCL`
- `common_SDL/`, `libs_SDL/`, `sketches_SDL/`, `apps/` built if `WITH_SDL`
- `apps_OCL/`, `sketches_OCL/` built if `WITH_OPENCL AND WITH_SDL`

### OBJECT Libraries

CMake OBJECT libraries are used to share source between executables and shared libraries:
- `MolecularEngine`, `DynamicOpt` — used by MolecularEditor, MolecularEditorOCL, Molecular/RigidMol libs
- `SDL2OGL` — used by all SDL2/OpenGL apps
- `Body2D`, `TerrainGrid2D`, `Noise` — used by CombatModels shared libraries

### ASAN Support

Two mechanisms:
1. **Global** (CMake 3.13+): `add_compile_options(-fsanitize=address)` + `add_link_options(-fsanitize=address)`
2. **Legacy**: `CMAKE_EXE_LINKER_FLAGS`, `CMAKE_SHARED_LINKER_FLAGS`, `CMAKE_MODULE_LINKER_FLAGS`

Build symlink approach (per `cpp-build` skill):
- `cpp/Build` → `cpp/Build-opt/` (production) or `cpp/Build-asan/` (debugging)
- Must rebuild all targets when switching — mixed `.o` files cause ASAN symbol errors
- `LD_PRELOAD=$(g++ -print-file-name=libasan.so)` for ASAN runs

### clangd Integration

- `CMAKE_EXPORT_COMPILE_COMMANDS ON`
- `compile_commands.json` symlinked from build dir to source dir and parent dir
- Enables clangd LSP for code navigation and completion

### Python Runtime Compilation

`cpp_utils.py`:
- `compile_lib(name)`: `g++ -std=c++11 -O3 -ftree-vectorize -unroll-loops -ffast-math -c -fPIC` then `g++ -shared`
- `make(what)`: Runs `make` in `cpp/Build/libs/CombatModels/`
- `loadLib(name, recompile)`: Compiles (optional) + loads via `ctypes.CDLL`
- `BUILD_PATH = ../../../cpp/Build/libs/CombatModels` (relative to `cpp_utils.py`)

### C Build (Legacy)

`C/build.sh`:
- Simple single-file: `gcc $CFLAGS $IFLAGS $target.c -o $target.x $LFLAGS`
- Links: GL, SDL2, m, tcc, dl
- No build system — just a convenience script

## Parity Status

- **CMake build ↔ Python runtime compilation**: Different toolchains. CMake uses system-configured compiler flags; `cpp_utils.py` hardcodes `-std=c++11 -O3 -ffast-math`. Potential mismatch if CMake uses C++20 and Python uses C++11.
- **CMake ASAN ↔ cpp-build skill**: Skill documents symlink approach; CMake implements global ASAN flags. Both approaches valid but must not be mixed.

## Open Issues

- All build options OFF by default — bare `cmake ..` produces minimal build with no apps
- `MolecularEditor2` commented out in `cpp/apps/CMakeLists.txt` — dead code
- `cpp_utils.py` uses `-std=c++11` while CMake uses `-std=c++20` — version mismatch
- `cpp_utils.py` uses `os.system()` for compilation — no error checking, silent failures
- `try/except: pass` in `cpp_utils.py` cleanup violates fail-fast principle
- No `install()` targets — no way to install built libraries system-wide
- No CI/CD configuration found
- `C/build.sh` uses `mkdir $BUILDDIR` without `-p` — fails if directory exists
- `C/build.sh` uses `rm $BINAME` without `-f` — fails if binary doesn't exist
- No cross-compilation support
- No package manager integration (vcpkg, conan) — dependencies found via `find_package`
- 51 `CMakeLists.txt` files — complex build tree, no build documentation beyond skill doc
