---
name: forcefield-validation
description: Implementing/debugging forcefields (UFF, SPFF) — invariants, switch semantics, buffer parity
trigger:
  glob:
    - "**/SPFF*"
    - "**/UFF*"
    - "**/common_resources/**/*.dat"
    - "**/common_resources/**/*.xyz"
    - "**/common_resources/**/*.mol"
---

## General Forcefield Invariants

Apply these before SPFF/UFF-specific checks:

### NaN Padding — Invalid-Access Alarm Bell
- Fill all padding/ghost slots with `NaN` (or `inf`, or a sentinel like `1e20`).
- Any accidental read poisons downstream math and makes the bug **visible immediately**.
- Never overload `NaN` with semantics like "pinned" — use a separate flag (`fixmask`).
- **Check**: after every GPU→CPU download, `assert np.isfinite(pos[real]).all()` and crash with a structured dump if it fails.

### Index Mapping — Bidirectional Consistency
For every forward mapping, verify the backward mapping exists and is consistent:
- `neighs[i,k] = j` → `bkSlots[j,slot] = i` → `revSlot[i,k] = slot_in_j`
- Build a Python truth table from the original topology (bond list, angles, exclusions).
- Parse kernel logs and assert every expected pair appears with the correct action (`EXCLUDE_BOND`, `COLLIDE`, `IGNORE_FAR`).

### Conservation Laws — The Physics Floor
When no reference exists, verify physical reasonableness:
- **Linear momentum**: `|sum(m_i * v_i)| < tol` after one step with no external forces.
- **Angular momentum**: `|sum(r_i × p_i) + sum(I_i * ω_i)| < tol` for rotation-capable nodes.
- **Energy drift** (MD only): track total energy; monotonic drift = unstable integrator.
- **Technique**: start with *nodes-only* (no capping atoms) to test core rigid-body logic, then add caps and verify cap recoils balance node rotation.

### Momentum Reset on Constraint Discontinuities
- Iterative solvers store momentum buffers (`dpos_mom`, `dquat_mom`, `vel`, `omega`).
- After pin/unpin, drag start/end, teleport, or mode switch: **zero these buffers** before the next step.
- Without this, stale momentum causes sudden jumps or divergence.
- **Critical for interactive GUIs** where dragging is a sequence of discontinuous constraint target changes.

### Isolation Before Combination
Test subsystems independently before enabling both:
1. **Collision-only**: `ENABLE_COLL=1, ENABLE_PORT=0`
2. **Port-only**: `ENABLE_COLL=0, ENABLE_PORT=1`
3. **Combined**: both enabled
If combined fails but isolation passes, the bug is in interaction logic, not either subsystem.

## Architecture

- Reference: `pyBall/SPFF_multi.py` (wrapper for C++ `SPFFmulti_lib.cpp`)
- Backends: `iParalel=0` (C++ CPU), `iParalel=2` (C++ OpenCL), `pyBall/OCL/UFF.py` (Python OpenCL)
- Data files: `AtomTypes.dat`, `BondTypes.dat`, `AngleTypes.dat` in `common_resources/`

## Switch Semantics

0=keep, >0=force-on, <0=force-off. Passing 0 does nothing.
- UFF: `setSwitchesUFF(DoBond, DoAngle, DoDihedral, DoInversion, DoAssemble, ...)`. `bSubtractBondNonBond` must align with `setSwitchesUFF_NB`.
- SPFF: `setSwitches2(..., SPFF, Angles, PiSigma, PiPiI)`. `bSPFF` is base class flag (True even for UFF). Check `bUFF` to distinguish.

## Buffer Parity

- Call `init_buffers(bUFF=...)` to populate C++ `buffers`/`ibuffers` maps. Read `ndims` first for shapes.
- UFF critical: Check `bonParams`, `angParams`, `a2f`. CPU `angParams` layout is `[K, c0, c1, c2, c3]`; kernels expect split `angParams1=[c0..c3]` and `angParams2_w=[K]`.
- SPFF critical: Check `apos` (Pi-orbitals at `natoms:natoms+nnode`), `bkNeighs` (Back Neighbors), `Ksp`/`Kpp`.

## Pitfalls

- Pi-orbitals: In SPFF, `apos` contains atoms AND pi-nodes. Loop bounds: `natoms` vs `natoms+nnode`.
- Node/cap layout: Current C++ builder sets `nnode=natoms`, allocates one pi slot per atom. Caps have `Ksp/Kpp=0` but occupy `nvecs`. `bkNeighs` sized `nSystems*nvecs`.

## Bonded-Only Parity Flow

`cleanF → getSPFFf4 → updateAtomsSPFFf4(dt=0)` to assemble recoil from `fneigh` via `bkNeighs` without moving atoms. Set `bSubtractVdW=0` when NonBonded off. Propagate switches to `ffls[isys]`. Copy `fapos` back to shared buffers.