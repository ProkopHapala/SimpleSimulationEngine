---
description: Domain-specific test catalog for bonded force fields: bonds, angles, dihedrals, and orbital-alignment surrogates.
---

# Intramolecular Forcefield Validation Skill

## Scope

Force fields that describe covalent geometry: bond stretching, angle bending, torsional barriers, and conjugation effects. This skill provides the oracle for what "correct" looks like.

---

## Canonical Test Progression

Test systems are ordered by complexity. Do not advance to the next system until the current one passes.

| System | DOFs Tested | Key Checks |
|--------|-------------|------------|
| **Diatomic** (e.g., H-H, N-N) | Bond stretch | Equilibrium distance, binding energy, force at dissociation limit |
| **Triatomic bent** (e.g., H-O-H) | Bond + angle | Angle ~104.5°, O-H ~0.96 Å, non-linear |
| **Triatomic linear** (e.g., O-C-O) | Bond + linear angle | 180° angle, symmetric forces |
| **Tetrahedral** (e.g., N-H₃) | Bond + angle | Pyramidal, H-X-H ~107°, inversion barrier |
| **Planar sp²** (e.g., H₂C=O) | Bond + angle + planarity | C=O double bond, all atoms in one plane |
| **Aromatic ring** (e.g., benzene) | Bond + angle + conjugation | Planar, equal bond lengths, delocalization energy |
| **Heteroaromatic** (e.g., pyridine) | As above + polarization | C-N bond shorter than C-C, dipole present |
| **Flexible chain** (e.g., butane) | Dihedral / torsion | Rotational barrier, gauche/trans minima |
| **Cyclic dimer** (e.g., carboxylic acid dimer) | H-bond geometry | O-H···O distance 1.6–2.0 Å, cyclic symmetry |

---

## Term-Specific Parity

### Bond Stretching

| Property | Expected Behavior | Failure Mode |
|----------|-------------------|--------------|
| Equilibrium length | Matches parameter table or reference | Wrong l0 → shifted minimum |
| Stiffness | Force scales linearly near minimum for harmonic | Wrong k → wrong vibrational frequency |
| Dissociation | Energy → finite limit or +∞ (no unphysical well) | Missing repulsion → artificial minimum at 0 |
| Force symmetry | `F(r0 + δ) = -F(r0 - δ)` | Sign error in derivative |

### Angle Bending

Two formulations exist in the wild. The agent must know which the reference uses:

| Formulation | Range | Behavior > 90° | When to Use |
|-------------|-------|-----------------|-------------|
| `cos θ` | Fast, few FLOPs | Not quasi-harmonic | Speed-critical, small angles only |
| `cos(θ/2)` | Requires `sqrt` | Quasi-harmonic everywhere | General purpose, all angles |

**Check**: Run angle scan 0°–180°. If energy curve folds back or force reverses sign near 90°, the wrong formulation is used.

### Dihedrals / Torsions

| Approach | Physics | GPU Cost | Parity Check |
|----------|---------|----------|--------------|
| Explicit 4-body | `E(φ) = Σ V_n cos(nφ)` | High (divergent threads) | Barrier heights, periodicity |
| π–π alignment surrogate | `E ~ k_pp (1 - cos²θ_pp)` | Low (2-body) | Planarity, conjugation preserved |
| π–σ orthogonality | `E ~ k_sp (cos²θ_sp)` | Low (2-body) | σ-bond plane ⊥ π-orbital |

**Critical**: If the reference uses explicit dihedrals and the target uses π-alignment surrogates, exact parity is impossible. Switch to qualitative parity: check planarity, barrier order-of-magnitude, and correct minimum positions.

---

## Recoil & Force Assembly

For GPU implementations that avoid atomic writes:

1. **Verify Newton's third law**: For every bonded pair (A,B), `F_A_from_B = -F_B_from_A`.
2. **Capping atom handling**: Terminal atoms (H, lone pairs) may receive forces only via recoil from their host. Verify they feel the correct equal-and-opposite force.
3. **No self-forces**: An atom must never feel a force from itself through a topology loop.

---

## Qualitative Checks

| System | Expected Geometry | Tolerance |
|--------|-------------------|-----------|
| Water | Bent, O-H ~0.96 Å, H-O-H ~104.5° | ±0.05 Å, ±5° |
| Ammonia | Pyramidal, N-H ~1.01 Å, H-N-H ~107° | ±0.05 Å, ±5° |
| Formaldehyde | Planar, C=O ~1.21 Å, H-C-H ~116° | ±0.05 Å, ±5° |
| Benzene | Planar, C-C ~1.40 Å, all angles 120° | ±0.05 Å, ±3° |
| Carboxylic dimer | Cyclic, O-H···O ~1.7 Å, symmetric | ±0.2 Å |

---

## Rules

1. **Start with diatomic, end with dimer.** Do not skip steps.
2. **Know your angle formulation.** `cos θ` and `cos(θ/2)` are not interchangeable.
3. **Explicit dihedrals and π-alignment surrogates are different methods.** Do not demand exact parity between them.
4. **Recoil forces must balance exactly.** Even tiny asymmetries accumulate in MD and break conservation laws.
