---
description: Domain-specific test catalog for molecule-surface interactions: grid potentials, adsorption geometry, and long-range electrostatics.
---

# Molecule-Surface Interaction Validation Skill

## Scope

Interactions between adsorbate molecules and rigid substrates, represented via precomputed grids, compact analytic bases, or direct pairwise summation.

---

## Canonical Test Progression

| System | Substrate | Adsorbate | What It Tests |
|--------|-----------|-----------|---------------|
| **Single ion scan** | Ionic crystal | Cation / anion | Charge–lattice registry, Madelung potential |
| **Rare gas atom** | Ionic / metal | Ne, Ar | Physisorption well, no chemisorption |
| **Water on ionic surface** | NaCl-like | H₂O | Dipole orientation: O→cation, H→anion |
| **Planar aromatic** | Ionic / metal | Benzene | π–surface distance, registry, corrugation |
| **Flexible organic** | Variable | Formic acid, oxalate | H-bond to surface, conformational change |

---

## Potential Representation Parity

### Precomputed 3D Grid (Interpolated)

| Check | Expected | Failure Mode |
|-------|----------|--------------|
| Grid origin aligned with substrate | Yes | Shifted grid → wrong adsorption sites |
| Sampling parity (corners vs. centers) | Consistent with reference | Half-step shift → systematic offset |
| Interpolation smoothness | C² continuous | Non-smooth → energy drift in MD |
| PBC in-plane (x,y) | Yes, with cell replication | Missing PBC → edge artifacts |
| Clamping out-of-plane (z) | Yes | Periodic z → atoms fall through surface |

### Compact Analytic Basis

| Check | Expected | Failure Mode |
|-------|----------|--------------|
| Basis captures lateral periodicity | cos(kx) / sin(kx) terms | Missing harmonics → wrong corrugation |
| Basis captures vertical decay | exp(−az) or polynomial | Wrong decay → wrong well depth |
| Fit residual small in accessible region | RMSE < 1% of well depth | Large residual → wrong forces |

### Direct Pairwise Summation

| Check | Expected | Failure Mode |
|-------|----------|--------------|
| Converges with enough periodic images | Yes | Too few images → wrong Madelung energy |
| Coulomb converges for charged slab | Requires Ewald / FFT | Brute force → divergent |

---

## Channel Separation

Substrate interaction is typically split into independent channels:

| Channel | Physics | Short-Range Behavior | Long-Range Behavior |
|---------|---------|----------------------|---------------------|
| Pauli repulsion | Electron overlap | Steeply rising ~1/r¹² | Negligible |
| London dispersion | Induced dipole | Attractive ~1/r⁶ | Slow decay |
| Coulomb electrostatics | Charge interaction | Damped finite | 1/r tail |

**Parity check**: Sample each channel independently. Pauli must be positive everywhere. London must be negative everywhere. Coulomb must change sign with charge product.

---

## Adsorption Geometry Qualitative Checks

| Molecule | Surface Site | Expected Orientation | Distance |
|----------|--------------|----------------------|----------|
| Water | Cation site | Oxygen closest to cation, dipole pointing away | O–cation: 2.0–3.5 Å |
| Water | Anion site | Hydrogens closest to anion | H–anion: 2.0–3.5 Å |
| Benzene | Hollow / bridge | Planar, ring parallel to surface | 3.2–3.8 Å |
| Carboxylic acid | Cation–anion bridge | Both oxygens interacting | O–ion: 2.0–3.0 Å |

---

## Coulomb Convergence

For charged or polar substrates, the Coulomb potential does not converge by direct real-space summation. Two approaches:

1. **FFT-based Poisson solver**: Solve ∇²V = −4πρ in reciprocal space.
2. **Ewald summation**: Split into short-range (real-space cutoff) and long-range (reciprocal-space sum).

**Parity check**: Compare GridFFT / Ewald result against brute-force summation with very large cutoff. RMSE should be < 1e-5 energy units.

---

## Rules

1. **Grid origin conventions differ between platforms.** Always verify alignment with a known reference point (e.g., ion position).
2. **Coulomb on charged substrates requires special treatment.** Do not assume brute-force summation converges.
3. **Adsorption orientation is the most sensitive qualitative check.** Water binding O-down vs. H-down reveals force sign errors.
4. **Test each channel independently.** A bug in one channel may be hidden if others dominate.
