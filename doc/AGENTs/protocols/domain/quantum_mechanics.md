---
description: Domain-specific test catalog for quantum mechanical methods: basis sets, orbital ordering, SCF convergence, and electronic structure.
---

# Quantum Mechanics Validation Skill

## Scope

Electronic structure methods: tight-binding, DFT, Hartree-Fock, and approximate quantum methods. This skill defines what "correct" means for wavefunctions, energies, and derived properties.

---

## Canonical Test Progression

| System | Basis | What It Tests |
|--------|-------|---------------|
| **Diatomic (e.g., H₂)** | s-only (minimal) | Bond formation, σ orbital, binding energy |
| **Diatomic (e.g., N₂)** | sp | Multiple bonds, σ/π separation, triple bond |
| **Heteronuclear (e.g., CO)** | sp | Polarization, dipole moment, charge asymmetry |
| **Triatomic bent (e.g., H₂O)** | sp | Angular geometry, lone pairs, dipole |
| **Tetrahedral (e.g., NH₃)** | sp | Pyramidal, inversion barrier |
| **Planar sp² (e.g., benzene)** | sp + π | Delocalization, aromaticity, ring current |
| **Charged system (e.g., Na⁺)** | sp | Ionization energy, charge conservation |

---

## Basis Set Ordering

Different platforms use different conventions for p-orbital ordering. The agent must know and explicitly verify the mapping:

| Convention | Order | Notes |
|------------|-------|-------|
| Chemistry (Gaussian) | s, px, py, pz | Standard Cartesian |
| Physics (spherical) | s, py, pz, px | Ortega/Fireball convention |
| Grid projection | px, py, pz, s | Used in some GPU kernels |

**Parity check**: After porting, compare orbital-by-orbital overlaps or Hamiltonian matrix elements. A silent permutation breaks all physics without changing total energy.

---

## SCF Convergence

| Check | Criterion | Failure Mode |
|-------|-----------|--------------|
| Energy change | ΔE < 1e-6 energy units | Oscillation, divergence |
| Density matrix change | max|ΔP| < 1e-5 | Charge sloshing |
| Orbital occupancy | Integer within 1e-3 | Fractional occupancy in closed-shell |
| Number of iterations | < 100 for simple systems | Bad mixing, no damping |

**Physical fallback**: If SCF fails to converge, check:
1. Initial guess quality (atomic vs. random)
2. Charge neutrality (total electrons = sum of nuclear charges)
3. Level shifting or damping parameters

---

## Orbital Energy Ordering

For a closed-shell system, orbital energies must satisfy:

```
E(occupied) < E(virtual)
E(homo) < E(lumo)
```

**Qualitative checks by system**:

| System | Expected HOMO Character | Expected LUMO Character |
|--------|------------------------|------------------------|
| H₂ | Bonding σ | Antibonding σ* |
| N₂ | Bonding σ (or degenerate π) | Antibonding π* |
| CO | Lone pair on C (or bonding σ) | Antibonding π* |
| H₂O | Lone pair on O (nonbonding) | Antibonding O-H σ* |
| Benzene | Degenerate π bonding | Degenerate π* antibonding |

---

## Density / Charge Parity

| Check | Method | Tolerance |
|-------|--------|-----------|
| Total charge integral | ∫ ρ(r) dr = N_electrons | 1e-4 |
| Bader charges (if available) | Sign consistent with electronegativity | Qualitative |
| Dipole moment | Direction: negative end toward more electronegative atom | Sign |
| Electron density at bond midpoint | ρ > 0, higher for multiple bonds | Order-of-magnitude |

---

## Qualitative Checks

| Property | Expected | Tolerance |
|----------|----------|-----------|
| Binding energy sign | Negative (stable molecule) | Must be negative |
| Bond order | H₂: 1, N₂: 3, CO: ~2.5 | Integer or half-integer |
| Planarity of sp² systems | All atoms coplanar | < 0.1 Å out-of-plane |
| Symmetry of wavefunctions | Matches point group | Visual / automatic |

---

## Rules

1. **Basis ordering is the silent killer of QM ports.** Verify orbital-by-orbital, not just total energy.
2. **SCF convergence is not guaranteed.** Provide damping, level shifting, and good initial guesses.
3. **Check orbital ordering, not just eigenvalues.** A wrong ordering with correct energies is still wrong chemistry.
4. **Start with minimal basis, then expand.** Do not test polarization functions until s and p are correct.
