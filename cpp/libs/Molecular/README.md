# Molecular

Shared libraries wrapping molecular force fields and dynamics for Python/ctypes access. Provides C API (`extern "C"`) wrappers around various molecular mechanics models: MMFF, rigid molecules, reactive force fields, electron force field (eFF), and multi-center Gaussian orbital force field (CLCFGO).

## Libraries

- **libMolecular.so** — `MolecularWorld` with MMFF force field, molecular optimization, point cloud comparison (links `MolecularEngine`, `DynamicOpt`)
- **libRigidMol.so** — rigid molecule dynamics: MMFF + `MMFFBuilder`, `DistanceHierarchy`, `AtomicConfiguration` (links `MolecularEngine`, `DynamicOpt`)
- **libReactiveFF.so** — Rigid Atom Reactive Force Field using `RARFFarr.h`
- **libeFF_lib.so** — Electron Force Field (eFF) using `InteractionsGauss.h` and `eFF.h`
- **libCLCFGO_lib.so** — Compact Linear Combination of Floating Gaussian Orbitals (disabled by default, requires SDL2 for `Plot2D.h`)

## Files

- **Molecular.cpp** — `Molecular` library: `MolecularWorld`, force field fitting (`FitFF.h`), point cloud comparison, optimization
- **RigidMol.cpp** — `RigidMol` library: rigid molecule dynamics with MMFF, distance hierarchy, atomic configuration I/O
- **ReactiveFF.cpp** — `ReactiveFF` library: Rigid Atom Reactive Force Field (`RARFFarr.h`)
- **eFF_lib.cpp** — `eFF_lib` library: Electron Force Field with Gaussian interactions, buffer registry for ctypes
- **CLCFGO_lib.cpp** — `CLCFGO_lib` library: multi-center floating Gaussian orbital force field, AO integrals, grid basis (disabled — requires SDL2)
- **CMakeLists.txt** — build targets: `Molecular`, `RigidMol`, `ReactiveFF`, `eFF_lib` (all SHARED). `CLCFGO_lib` commented out.

## Related Documentation

- [docs/MolGUI_web.md](../../docs/MolGUI_web.md) — molecular GUI web version
- [doc/TopicalAudit/molecular-editor.md](../../doc/TopicalAudit/molecular-editor.md) — molecular editor audit
- [doc/TopicalAudit/python-cpp-bindings.md](../../doc/TopicalAudit/python-cpp-bindings.md) — ctypes binding patterns
