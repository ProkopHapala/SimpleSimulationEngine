---
description: Domain-specific skill for building and validating molecular topology: graph construction, hybridization assignment, and conversion to flat evaluation arrays.
---

# Topology Building Validation Skill

## Scope

The two-phase model of molecular simulation:
- **Phase 1**: Dynamic graph editing (atoms, bonds, rings, hybridization)
- **Phase 2**: Static flat arrays for force field evaluation (neighbor lists, parameter tables)

This skill validates the conversion boundary between phases.

---

## Canonical Graph Operations

### Bond Detection

| Check | Method | Failure Mode |
|-------|--------|--------------|
| All covalent bonds found | Distance < sum of covalent radii × factor | Missing bond → wrong hybridization |
| No spurious bonds | Distance > cutoff | Extra bond → wrong ring detection |
| Periodic image bonds handled | Minimum image convention | Bonds across cell boundaries missed |
| Bond order assigned | Single/double/triple from distance | Wrong order → wrong parameters |

### Hybridization (VSEPR)

| Domains | Hybridization | Geometry | Examples |
|---------|--------------|----------|----------|
| 2 | sp | Linear | CO₂, acetylene |
| 3 | sp² | Trigonal planar | Benzene, formaldehyde |
| 4 | sp³ | Tetrahedral | Methane, water (2 bonds + 2 lone pairs) |

**Check**: Count σ bonds + π bonds + lone pairs = total domains. Map domains → hybridization → geometry.

---

## Two-Phase Consistency

After converting from dynamic graph to flat arrays, verify:

| Check | Reference vs. Target |
|-------|-------------------|
| Same number of atoms | natoms_ref == natoms_new |
| Same number of bonds | nbonds_ref == nbonds_new |
| Same bond atom pairs | set(bonds_ref) == set(bonds_new) |
| Same neighbor count per atom | len(neigh_ref[i]) == len(neigh_new[i]) |
| Same angle definitions | set(angles_ref) == set(angles_new) |
| Same atom type assignments | type_ref[i] == type_new[i] |
| Same π-fragment IDs | pi_fragment_ref[i] == pi_fragment_new[i] |

---

## Node vs. Terminal Atom Separation

Many force fields distinguish "node" atoms (heavy, multiply bonded) from "terminal" atoms (H, lone pairs).

| Property | Node Atom | Terminal Atom |
|----------|-----------|---------------|
| Angular terms | Full (angles, dihedrals, π-alignment) | None (recoil only) |
| Neighbor count | ≥ 2 | 1 |
| Force computation | Independent thread | Inherits from host |

**Check**: After conversion, verify that every terminal atom has exactly one bonded neighbor and receives force only via recoil from that neighbor.

---

## Pi-System Detection

Conjugated and aromatic systems require special handling.

| Check | Algorithm | Failure Mode |
|-------|-----------|--------------|
| All sp² carbons in same ring found | DFS/BFS for connected π-network | Missed conjugation → no aromaticity |
| π-orbital directions computed | Normal to σ-plane | Wrong direction → wrong π–π interaction |
| Ring planarity enforced | All ring atoms coplanar within tolerance | Puckered ring → broken conjugation |

---

## Qualitative Checks

| Molecule | Expected Topology | Check |
|----------|-------------------|-------|
| Water | 2 O-H bonds, O has 2 lone pairs | 4 domains → sp3, bent |
| Benzene | 6 C-C bonds, each C has 1 π orbital | 6-membered ring, all sp2, planar |
| Pyridine | As benzene, but one N has lone pair instead of H | Same ring, N correctly typed |
| Carboxylic acid dimer | Two H-bonds between monomers | Dimer topology, not two separate molecules |

---

## Rules

1. **The graph is the source of truth.** If the graph is wrong, all subsequent physics is wrong.
2. **Hybridization follows from bond count, not vice versa.** Do not assign sp2 because the atom is carbon in a ring — count bonds and lone pairs.
3. **After any reordering, rebuild topology from scratch.** Do not incrementally patch — subtle index shifts accumulate.
4. **Test topology independently of force evaluation.** A correct graph with wrong parameters is easier to fix than a wrong graph with correct parameters.
