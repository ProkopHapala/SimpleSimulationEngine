# MolecularEditor

Interactive molecular editor and molecular dynamics simulator. Build molecules atom-by-atom, assign force field parameters, run MD simulations, and visualize results in 3D. Supports MMFF force field, conformational search, and various force field demos.

## Features

- Interactive 3D molecular editing (add/remove atoms, bonds)
- MMFF force field parameter assignment and energy evaluation
- Molecular dynamics simulation with real-time visualization
- Conformational search via global optimization (FIRE algorithm)
- Multiple force field demos (eFF, CLCFGO, MMFF)
- 3D rendering of atoms, bonds, electron orbitals

## Files

- **MolecularEditor_main.cpp** — main application: editor GUI, 3D rendering, MD simulation
- **MolecularEditor_old_main.cpp** — older version of the editor
- **ConfSearch.cpp** — conformational search: global optimization of molecular geometry
- **MolecularDraw.h** — 3D rendering of molecules: atoms, bonds, orbitals, forces
- **test_MMFF.cpp** — MMFF force field test/demo
- **test_CLCFGO_OGL3.cpp** — CLCFGO force field demo with OpenGL 3
- **CMakeLists.txt** — build targets: `MolecularEditor_main`, `test_MMFF`, `test_CLCFGO_OGL3`
- **CMakeLists-old.txt** — older CMake configuration
- **data/** — molecular data files
- **doc/** — documentation (global optimization, force field notes)
- **inputs/** — molecular input files (48 files: .xyz, .mol, etc.)
- **py/** — Python helper scripts
- **notes/** — development notes

## Related Documentation

- [docs/MolGUI_web.md](../../docs/MolGUI_web.md) — molecular GUI web version documentation
