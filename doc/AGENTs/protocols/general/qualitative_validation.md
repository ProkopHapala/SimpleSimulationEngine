---
description: Validate that computed results obey physical laws, geometric constraints, and conservation principles regardless of numerical parity.
---

# Qualitative Validation Protocol

## Purpose

Not all correctness can be expressed as RMSE. A new method may legitimately differ from a reference while still being physically sound. This protocol checks:

- Conservation laws (energy, momentum)
- Symmetry invariance (translation, rotation, reflection)
- Geometric plausibility (bond lengths, angles, planarity)
- Physical intuition (adsorption orientation, charge distribution)

It runs after parity checks (if parity was required) and before performance optimization.

---

## Agentic Loop Integration

```
parity checking (if same-method port)
    |
    v
THIS PROTOCOL (Level 5)
    |
    v
performance optimization (Level 6)
    |
    v
human review (only if flagged)
```

---

## Checks

### 1. Numerical Sanity (Automatic)

Cheap, run on every iteration:

| Check | Assertion |
|-------|-----------|
| No NaN / Inf | `np.isfinite(arr).all()` |
| Energy in reasonable range | per system class (see domain skill) |
| Forces bounded | `max|F| < F_max` |
| Coordinates not exploded | `max displacement per step < threshold` |
| No atom overlap | `all pairwise distances > r_hardcore` |

---

### 2. Conservation Laws (Automatic)

| Law | Test | Tolerance |
|-----|------|-----------|
| Energy conservation (NVE) | `ΔE_total / E_total < 1e-4` over 100–1000 steps | 1e-4 |
| Momentum conservation (isolated) | `|sum(F_i)| < 1e-6` | 1e-6 |
| Angular momentum (isolated, no torque) | `|dL/dt| < 1e-6` | 1e-6 |
| Virial theorem (stationary state) | `2<T> + <Σ r·F> ≈ 0` | 1% |

**Note**: Energy conservation tolerance must be relaxed if the method itself is approximate (e.g., a forcefield with cutoff). Compare against the reference method's drift, not against zero.

---

### 3. Symmetry Invariance (Automatic)

| Operation | Test |
|-----------|------|
| Translation | `E(R + Δ) == E(R)` |
| Rotation | `E(R) == E(U·R)` for orthogonal U |
| Permutation of identical atoms | `E(swap identical) == E(original)` |
| Reflection (if applicable) | `E(mirror) == E(original)` |

For forces:
- `sum_i F_i == 0` for isolated systems
- `F_ij == -F_ji` for every interacting pair

---

### 4. Geometric Plausibility (Automatic + LLM)

After geometry relaxation or MD, check:

| Property | Check | Tolerance |
|----------|-------|-----------|
| Bond lengths | within known covalent range | ±10% of tabulated value |
| Bond angles | within VSEPR-predicted range | ±10° |
| Planarity (sp²) | max out-of-plane distance < 0.1 Å | 0.1 Å |
| Tetrahedral (sp³) | angle ~109.5° | ±5° |
| Ring closure | bonded path length matches expected | exact |
| No unphysical dissociation | max bond stretch < 2× equilibrium | 2× |

**LLM role**: When automatic checks are borderline (e.g., bond length 5% longer than expected but within tolerance), the LLM reads the geometry report and decides whether this is physically meaningful (e.g., strained ring) or a bug.

---

### 5. Physical Intuition (LLM + Human)

These are domain-specific and require reasoning. Examples by domain:

| Domain | Intuition Check |
|--------|----------------|
| Polar molecule + ionic surface | Oxygen (negative end of dipole) binds to cation; hydrogens (positive end) to anion |
| Hydrogen bond | H-bond donor–acceptor distance 1.5–2.0 Å, angle near 180° |
| Aromatic stacking | Planar rings parallel, interplanar distance 3.2–3.8 Å |
| Charge distribution | Partial charges sum to total charge; sign consistent with electronegativity |
| Orbital ordering | Bonding orbital lower in energy than antibonding |

The LLM reads the computed geometry/energies and issues a verdict:

```
Physical Intuition Report:
  O-cation distance: 2.31 Å  (expected 2.0–3.5 Å)  PASS
  H-anion distance:  2.45 Å  (expected 2.0–3.5 Å)  PASS
  Dipole orientation: O toward cation  PASS
  Verdict: QUALITATIVE BEHAVIOR CORRECT
```

If any check is borderline or FAIL, generate a `.xyz` snapshot and plot for human review.

---

## Artifact Naming Convention

```
debug/<date>_<task>/
  ├── conservation_NVE.png          # Energy vs time
  ├── conservation_momentum.txt   # Momentum drift
  ├── symmetry_rotation.txt       # Symmetry test results
  ├── geometry_relaxed.xyz        # Final structure
  ├── geometry_report.md          # Bond/angle/dihedral table
  └── qualitative_intuition.md    # LLM verdict
```

---

## Rules

1. **Conservation laws are non-negotiable.** If energy drifts or momentum is not conserved, the implementation is wrong, not "approximate."
2. **Geometric checks use known physical ranges, not reference values.** A new method may predict a slightly different bond length, but it must still be physically reasonable (e.g., C-C bond 1.2–1.6 Å, not 0.5 Å or 5 Å).
3. **LLM interprets borderline cases; human reviews only failures.** Do not escalate every geometric check to human.
4. **Symmetry checks catch index-mapping bugs** that topology checks miss (e.g., force on atom 5 equals force on atom 3 after permutation, revealing a swapped index).
