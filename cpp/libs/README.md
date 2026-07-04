# cpp/libs — Shared Libraries (Python/ctypes)

Shared libraries (`.so`) that wrap C++ simulation models and expose a C API (`extern "C"`) for Python/ctypes access. Each library is built with `-fPIC` for position-independent code. Built via CMake, loaded from Python using `ctypes.CDLL`.

## Library Index

### Combat Models

| Library | Description | Key Dependencies |
|---------|-------------|------------------|
| [CombatModels](CombatModels/README.md) | General combat: weapon/unit types, division-level combat | `MinimalDivisionLevel.h` |
| [SpaceCombatLib](CombatModels/README.md) | Space combat: guns, projectiles, Whipple shields, damage | `spaceCombat.h` |
| [LandCombatLib](CombatModels/README.md) | Land combat: LT gun/unit types, unit creation | `apps/LandTactics/`, `Body2D` |
| [LandCraftLib](CombatModels/README.md) | LandCraft world: terrain, hydraulics, roads, economy, pathfinding | `apps/LandCraft/`, `TerrainGrid2D`, `Noise` |

### Molecular

| Library | Description | Key Dependencies |
|---------|-------------|------------------|
| [Molecular](Molecular/README.md) | MolecularWorld with MMFF, force field fitting, optimization | `MolecularEngine`, `DynamicOpt` |
| [RigidMol](Molecular/README.md) | Rigid molecule dynamics: MMFF, distance hierarchy | `MolecularEngine`, `DynamicOpt` |
| [ReactiveFF](Molecular/README.md) | Rigid Atom Reactive Force Field (RARFFarr) | — |
| [eFF_lib](Molecular/README.md) | Electron Force Field with Gaussian interactions | `eFF.h`, `InteractionsGauss.h` |

### Space & Nuclear Physics

| Library | Description | Key Dependencies |
|---------|-------------|------------------|
| [KosmoSuite](KosmoSuite/README.md) | N-body, orbital mechanics, launch trajectories, EM fields, fission, neutrons, shock waves | `ODEintegrator.h`, `spline_hermite.h` |
| [SpaceTacticsLib](OrbitalWar/README.md) | Space tactics world: orbital dynamics, trajectory buffers | `SpaceWorld.h`, `Body` |
| [Shock](Shock/README.md) | 1D shock wave simulation through layered materials | `Shock1D.h`, `ShockWaves` |

### Flight & Naval

| Library | Description | Key Dependencies |
|---------|-------------|------------------|
| [Flight](libFlight/README.md) | Aircraft flight sim: aerodynamics, control surfaces, projectiles | `AeroSurf`, `AeroCraft`, `Body` |
| [libSail](libSail/README.md) | Sailing ship hydrostatics: hull volume, moment, static stability | `lineSearch.h` |

## Build

All libraries are built as SHARED libraries via CMake:
```bash
cd cpp/Build
cmake .. -DWITH_SDL=ON
make -j$(nproc) LandCraftLib    # build specific library
```

## Python Usage

Libraries are loaded from Python via ctypes:
```python
import ctypes
lib = ctypes.CDLL("cpp/Build/libs/CombatModels/libLandCraftLib.so")
lib.landcraft_init(512, 512)
```

See `doc/TopicalAudit/python-cpp-bindings.md` for the full ctypes binding pattern.

## Related Documentation

- [doc/TopicalAudit/python-cpp-bindings.md](../../doc/TopicalAudit/python-cpp-bindings.md) — ctypes binding patterns and C API design
- [doc/TopicalAudit/build-system.md](../../doc/TopicalAudit/build-system.md) — CMake build system, `-fPIC` requirements
- [docs/SpaceTactics/](../../docs/SpaceTactics/) — space tactics and combat design
- [docs/LandCraft/](../../docs/LandCraft/) — LandCraft architecture and subsystems
- [docs/AeroCombat/](../../docs/AeroCombat/) — aerodynamic model documentation
