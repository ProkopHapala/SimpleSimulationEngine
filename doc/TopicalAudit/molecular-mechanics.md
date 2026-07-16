---
type: TopicalAudit
title: Molecular Mechanics & Molecular Dynamics
tags: [topic, cross-language, physics, molecular, forcefield]
---

## Summary

Molecular mechanics forcefield (`MMFF.h`) for soft-body molecular simulation, and rigid-body molecular dynamics (`MolecularWorld.h`) for coarse-grained faster simulation. Multiple molecular editor apps (C++ and OCL). Python wrapper in `pyMolecular/`. Domain protocols in `doc/AGENTs/protocols/domain/` define forcefield validation, noncovalent interactions, intramolecular forcefields, quantum mechanics, and topology building.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/common/dynamics/Forces.h` | active | Pairwise force kernels: LJ, Morse, Coulomb, spring, repulsion, angular terms. Smoothstep cutoffs (r^2-based). Foundation for molecular dynamics force evaluation |
| C++ | `cpp/common/dynamics/MMFF.h` | active | MMFF forcefield implementation |
| C++ | `cpp/common/dynamics/MolecularWorld.h` | active | Rigid-body molecular dynamics (coarse-grained) |
| C++ | `cpp/apps/MolecularEditor/` | active | Rigid-body MD editor with Lua scripting |
| C++ | `cpp/apps/MolecularEditor2/` | active | Soft-body molecular mechanics editor |
| C++ | `cpp/apps_OCL/MolecularEditorOCL/` | experimental | GPU-accelerated molecular editor |
| C++ | `cpp/apps/MolecularEditor/doc/GlobalOptimization.md` | active | Global optimization for molecular simulation |
| Python | `python/pyMolecular/` | active | Python wrapper for rigid-body MD |
| Doc | `doc/AGENTs/protocols/domain/intramolecular_forcefields.md` | active | Protocol: intramolecular forcefield validation |
| Doc | `doc/AGENTs/protocols/domain/noncovalent_interactions.md` | active | Protocol: noncovalent interaction validation |
| Doc | `doc/AGENTs/protocols/domain/quantum_mechanics.md` | active | Protocol: QM reference calculations |
| Doc | `doc/AGENTs/protocols/domain/molecule_surface.md` | active | Protocol: molecular surface representation |
| Doc | `doc/AGENTs/protocols/domain/topology_building.md` | active | Protocol: molecular topology construction |
| Doc | `doc/AGENTs/skills/forcefield-validation/SKILL.md` | active | Skill: forcefield validation workflow |
| Doc | `docs/MolGUI_web.md` | active | Web-based molecular GUI documentation |

## Parity Status

- **C++ `MMFF.h` ↔ Python `pyMolecular/`**: Python wraps C++ rigid-body MD via ctypes. Forcefield params shared.
- **Rigid-body MD ↔ Soft-body MD**: Two different approaches (coarse-grained rigid bodies vs. full MMFF forcefield). Not directly comparable.
- **Domain protocols** define validation criteria but no automated parity test found.

## Open Issues

- Global optimization algorithms planned but incomplete (see README "Planned" section)
- Multipole expansion and FMM for long-range interactions: planned, not implemented
- `MolecularEditorOCL` status unclear — GPU acceleration experimental
- Forcefield validation skill exists but no regression test suite found
- `cpp/apps/MolecularEditor/notes/GlobalOpt.md` — notes on global optimization
