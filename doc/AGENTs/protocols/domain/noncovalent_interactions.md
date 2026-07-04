---
description: Domain-specific test catalog for pairwise non-bonded interactions: dispersion, electrostatics, repulsion, and hydrogen bonding.
---

# Non-Covalent Interaction Validation Skill

## Scope

Pairwise and many-body interactions between non-bonded atoms: van der Waals, electrostatics, short-range repulsion, and hydrogen-bond corrections.

---

## Canonical Test Progression

| System | Interaction Types | Key Checks |
|--------|-------------------|------------|
| **Rare gas dimer** | Pure dispersion (Lennard-Jones) | Well depth, equilibrium distance, repulsive wall |
| **Ion pair** | Coulomb + repulsion | Binding energy, distance near ionic radii sum |
| **Water dimer** | LJ + Coulomb + H-bond | O···H distance ~1.8 Å, binding ~−5 kcal/mol |
| **Benzene dimer** | Dispersion + quadrupole | π–π stacking distance ~3.5 Å, binding ~−2 kcal/mol |
| **Polar molecule + ionic surface** | All terms + registry | Adsorption site, orientation, corrugation |

---

## Potential Form Parity

### Lennard-Jones (12-6)

| Check | Expected | Failure Mode |
|-------|----------|--------------|
| Minimum at r = 2^(1/6) σ | Yes | Wrong exponent, wrong mixing rule |
| E → 0 as r → ∞ | Yes | Missing cutoff tail correction |
| E → +∞ as r → 0 | Yes | Division by zero, NaN |
| Force = −dE/dr | Yes | Sign error in derivative |

### Morse Potential

| Check | Expected | Failure Mode |
|-------|----------|--------------|
| Minimum at r = r₀ | Yes | Wrong r₀ → shifted minimum |
| Well depth = Dₑ | Yes | Wrong dissociation energy |
| E → 0 as r → ∞ | Yes | Missing asymptote |
| E finite at r = 0 | Yes | Exponential prevents NaN |

### Damped Coulomb

| Check | Expected | Failure Mode |
|-------|----------|--------------|
| E → QᵢQⱼ/r as r → ∞ | Yes | Damping too strong → wrong long-range |
| E finite at r = 0 | Yes | Undamped → NaN |
| Force on like charges repulsive | Yes | Sign error → attraction |
| Force on opposite charges attractive | Yes | Sign error → repulsion |

---

## Exclusion Rules

Bonded neighbors must be excluded from non-bonded evaluation to avoid double-counting.

| Exclusion Type | Distance in Bond Graph | Typical Rule |
|----------------|------------------------|--------------|
| 1-2 (bonded) | 1 | Always excluded |
| 1-3 (angle) | 2 | Always excluded |
| 1-4 (dihedral) | 3 | Often scaled, sometimes excluded |

**Parity check**: Compare exclusion lists between reference and target atom-by-atom.

---

## Mixing Rules

For cross-species parameters between atom types A and B:

| Parameter | Arithmetic Mean | Geometric Mean |
|-----------|-----------------|---------------|
| σ (LJ radius) | (σₐ + σᵦ)/2 | √(σₐ·σᵦ) |
| ε (LJ well depth) | — | √(εₐ·εᵦ) |

The reference and target must use the **same mixing rule**. Switching between arithmetic and geometric mean silently changes all heterogeneous interactions.

---

## Qualitative Checks

| System | Expected Behavior | Tolerance |
|--------|-------------------|-----------|
| Like charges | Repulsive, E > 0 at close range | Positive energy |
| Opposite charges | Attractive, E < 0 | Negative energy |
| Rare gas dimer | Weakly attractive, E ~ −0.01 to −0.5 eV | Order-of-magnitude |
| Water dimer | H-bond: O···H < 2.0 Å, binding ~−5 kcal/mol | ±2 kcal/mol |
| π–π stacking | Planar, separation 3.2–3.8 Å | ±0.5 Å |

---

## Rules

1. **Sign errors in Coulomb are catastrophic.** Like charges must repel.
2. **Exclusion rules must match exactly.** A missing 1-3 exclusion changes angle equilibrium.
3. **Mixing rules are part of the method.** Do not "improve" them during porting.
4. **Test ion pair before neutral dimer.** Coulomb sign errors are easiest to catch with charged systems.
