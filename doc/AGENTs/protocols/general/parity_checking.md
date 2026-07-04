---
description: Validate reimplemented methods against reference via exact numerical parity or qualitative landscape agreement.
---

# Parity Checking Protocol

## Purpose

Determine whether a reimplemented method (new language, platform, or optimized kernel) produces results equivalent to a trusted reference implementation.

Two modes exist:

- **Exact parity** — same algorithm, different code. Expect numerical identity within floating-point tolerance.
- **Qualitative parity** — different formulation or approximation. Expect physically equivalent behavior (same minima, barriers, trends) even if absolute energies differ.

---

## Agentic Loop Integration

This protocol is invoked after compilation and basic runtime sanity pass. It gates all higher-level qualitative and performance checks.

```
edit_code
    |
    v
compile (Level 0)
    |
    v
runtime sanity (Level 1 — no NaN, finite numbers)
    |
    v
THIS PROTOCOL (Level 2–4)
    |
    v
qualitative validation
    |
    v
performance optimization
```

If any parity check fails, the loop feeds the diff report back to the agent for diagnosis. Do not proceed to qualitative or performance checks until parity is resolved or intentionally waived.

---

## Test Hierarchy

### Level 2: Single-Point Quantitative Parity

Run the same input geometry through both reference and target implementations. Compare:

| Quantity | Exact Parity Tolerance (double) | Exact Parity Tolerance (single) |
|----------|-----------------------------------|---------------------------------|
| Total energy | `rtol = 1e-6` | `rtol = 1e-3` |
| Per-atom forces | `atol = 1e-5` energy/length | `atol = 1e-3` energy/length |
| Stress / virial | `atol = 1e-4` pressure units | `atol = 1e-2` pressure units |
| Derived properties (charges, densities) | `rtol = 1e-5` | `rtol = 1e-3` |

**Force comparison**: Use both magnitude RMSE and **direction cosine**.

```
cos(F_ref, F_test) = dot(F_ref, F_test) / (|F_ref| * |F_test|)
```

A direction cosine < 0.999 indicates a sign error or missing term even when magnitudes look reasonable.

**When exact parity is impossible** (different method): skip this level and proceed to Level 3.

---

### Level 3: Rigid Scan Parity

Compare energy/force profiles along controlled deformation coordinates. More discriminating than a single point because it reveals shifted minima, incorrect curvature, or missing barriers.

**Scan types**:

| Coordinate | Typical Range | Steps | What It Tests |
|------------|-------------|-------|---------------|
| Bond stretch | 0.5–2.0 × equilibrium | 20–50 | Bond term, dissociation limit |
| Bond angle | 60°–180° | 20–30 | Angular potential, equilibrium geometry |
| Dihedral / torsion | 0°–360° | 36–72 | Rotational barrier, conjugation |
| Molecule–surface distance | 2–8 Å | 30–50 | Adsorption well depth, orientation |
| Lateral position (x,y) over lattice | Unit cell + margin | 20×20 | Registry, corrugation |

**Metrics**:

| Metric | Exact Parity Pass | Qualitative Parity Pass |
|--------|-------------------|------------------------|
| Energy RMSE | < 1e-4 eV/atom | < 0.1 eV/atom (method-class dependent) |
| Minimum position shift | < 0.01 Å | < 0.1 Å |
| Barrier height error | < 1% | Same order of magnitude |
| Curvature at minimum (frequency) | < 1% | Within 10% |
| Curve shape correlation | Pearson r > 0.9999 | Pearson r > 0.99 |

**Artifact**: Generate overlay plot with both curves labeled. Reference should be distinguishable from target by line style (e.g., solid vs. dotted), **not** by color alone. Caption must include: scan type, RMSE, max deviation, minimum positions.

---

### Level 4: Systematic Shift Handling

When the target method is intentionally different (better basis, different functional, approximate Hamiltonian), absolute parity is meaningless. Instead:

1. **Compute correlation**: Pearson r between energy profiles. r > 0.99 required.
2. **Compare relative energies**: Energy differences between configurations must match, even if absolute offset differs.
3. **Check monotonicity**: If reference has barrier A > barrier B, target must preserve this ordering.
4. **Verify asymptotic behavior**: At large distances, both must converge to the same physical limit (e.g., dissociated atoms, no interaction).

The agent must report:

```
Exact parity:    SKIPPED (different method class)
Correlation:     0.997  PASS
Relative barrier: 0.12 eV (ref 0.11 eV)  PASS
Asymptote match:  < 1e-3 eV at r = 5.0 Å  PASS
```

---

## Tolerance Strategy Selection

The agent chooses the tolerance strategy based on the task:

| Strategy | Use When | Tolerance Source |
|----------|----------|----------------|
| `exact` | Direct port, same algorithm | Machine epsilon × problem scale |
| `relative` | Same method, different precision (double → single) | User-specified or 1e-3 default |
| `systematic_shift` | Better basis, improved approximation | Correlation + relative barriers |
| `method_class` | Entirely different method family | Domain skill file (e.g., forcefield vs. QM) |

---

## Artifact Naming Convention

All parity artifacts must be self-describing:

```
debug/<date>_<task>/
  ├── scan_bond_<label>.png          # Level 3 overlay plot
  ├── scan_angle_<label>.png
  ├── scan_surface_<label>.png
  ├── parity_singlepoint.txt          # Level 2 numeric comparison
  ├── parity_scans.txt                # Level 3 RMSE table
  └── REPORT.md                       # Human-readable summary
```

---

## Failure Modes & Diagnosis

| Symptom | Likely Cause | Agent Action |
|---------|--------------|--------------|
| Energy matches, forces wrong | Missing force derivative, sign error in gradient | Check force assembly pattern |
| Small system OK, large system fails | Buffer underallocation, uninitialized memory | Run topology + array-bound checks |
| Scan shape matches, minimum shifted | Equilibrium parameter mismatch (l0, r0) | Compare parameter tables |
| Chaotic scan (noise, spikes) | Race condition, atomic collision, bad indexing | Topology check + single-thread CPU reference |
| Systematic offset only | Different zero of energy (normal for different methods) | Switch to `systematic_shift` strategy |
| Forces finite but direction cosine < 0 | Total sign flip or coordinate system mismatch | Check coordinate handedness / index mapping |

---

## Rules

1. **Never skip parity because it is inconvenient.** If exact parity is impossible, explicitly switch to qualitative parity and document why.
2. **Always produce plots.** Numeric RMSE is necessary but not sufficient; human review of scan shapes catches issues numbers hide.
3. **Report relative barriers, not just absolute energies.** Chemistry is driven by energy differences.
4. **Cache reference scan outputs.** Do not regenerate reference data on every iteration.
