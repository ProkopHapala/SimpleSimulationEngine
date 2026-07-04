---
description: Verify that structural mappings (neighbor lists, angle definitions, index padding) are identical between reference and reimplemented code.
---

# Topology Verification Protocol

## Purpose

Topology errors are silent killers: the code compiles, numbers are finite, but forces point to the wrong atoms. This protocol catches mapping mismatches before any physics check runs.

It applies whenever:
- Porting from one language/platform to another
- Changing memory layout (AoS → SoA, adding/removing padding)
- Reordering atoms for performance (node-first vs. original order)
- Generating neighbor lists or bonded interaction lists from a dynamic graph

---

## Agentic Loop Integration

```
compile (Level 0)
    |
    v
runtime sanity (Level 1)
    |
    v
THIS PROTOCOL (Level 2)
    |
    v
parity checking (Level 3–4)
    |
    v
qualitative validation
```

If topology fails, abort immediately. Do not run scans or qualitative checks — wrong topology means wrong physics by definition.

---

## Checks

### 1. Neighbor List Isomorphism

For every atom, the set of neighbors must be identical between reference and target, ignoring order within the list.

```
for each atom i:
    ref_neighbors = set(neigh_ref[i] excluding padding_value)
    new_neighbors = set(neigh_new[i] excluding padding_value)
    assert ref_neighbors == new_neighbors
```

**Padding value**: Explicit sentinel (typically -1) must be used. `0` is never acceptable as padding because atom 0 is a valid index.

**Order tolerance**: Neighbors may be reordered between implementations, but the *set* must match. If order matters for your algorithm, enforce list equality, not just set equality.

---

### 2. Back-Neighbor Consistency

In a directed neighbor list (i → j), the back-neighbor is the index of i in j's neighbor list. It must exist and be correct for force recoil / Newton's third law.

```
for each atom i:
    for each neighbor j in neigh[i]:
        k = index_of(i) in neigh[j]
        assert k >= 0, "Missing back-neighbor"
        assert back_neigh[j][k] == i
```

**Common failure**: After reordering atoms (e.g., node-first for GPU), back-neighbor indices are not updated, causing recoil forces to be deposited on wrong atoms.

---

### 3. Angle & Dihedral Atom Index Mapping

For every angle (i,j,k) and dihedral (i,j,k,l):

```
for each angle in angles_ref:
    mapped = [atom_map[x] for x in angle]
    assert set(mapped) in angles_new
```

The `atom_map` accounts for reordering between reference and target. If no reordering occurred, `atom_map` is identity.

**Critical**: Verify that the *middle* atom of each angle (the vertex) is the same physical atom, not just any atom with the same element type.

---

### 4. Padding Pattern Match

Unused slots in padded arrays must have identical padding positions and identical sentinel values.

```
assert (arr_ref == padding_value) == (arr_new == padding_value)
assert no value outside [valid_range] exists in non-padding slots
```

**Common failure**: Partially initialized arrays where some padded slots contain garbage that is later read (e.g., in a `for` loop that does not check the sentinel).

---

### 5. Self-Neighbor Presence

For methods that store diagonal blocks or self-interactions in a sparse/block format, verify that each atom has a self-entry if and only if the reference does.

```
self_present_ref = any(neigh_j == i and cell_shift == 0 for ...)
self_present_new = any(neigh_j == i and cell_shift == 0 for ...)
assert self_present_ref == self_present_new
```

**Note**: The self-neighbor may not be at a fixed index; it must be detected by atom identity and cell shift, not by position in the list.

---

### 6. Exclusion List Consistency

Bonded neighbors excluded from non-bonded evaluation must match:

```
for each atom i:
    ref_excl = set(excl_ref[i])
    new_excl = set(excl_new[i])
    assert ref_excl == new_excl
```

**Failure modes**: Missing 1-3 exclusions (angle neighbors), incorrect handling of periodic images, or different bond-order thresholds for exclusion.

---

## Topology Diff Report

When a mismatch is found, generate a concise text report the LLM can diagnose:

```
TOPOLOGY MISMATCH — Atom 7
  Reference neighbors: {3, 5, 8}
  New neighbors:       {3, 5, 8, 8}     ← DUPLICATE

TOPOLOGY MISMATCH — Angle 12
  Reference atoms: (2, 5, 9)
  New atoms:       (2, 5, 10)         ← INDEX SHIFT

TOPOLOGY MISMATCH — Atom 3
  Reference padding: [7, 8, -1, -1]
  New padding:     [7, 8, 0, -1]       ← 0 USED AS PADDING
```

---

## Rules

1. **Never trust that topology "should be the same."** Prove it on every reimplementation.
2. **-1 is the only acceptable padding sentinel.** Zero is a valid atom index.
3. **Topology checks are automatic, not LLM-reviewed.** The diff report is text; the LLM diagnoses it, but the comparison itself is a deterministic script.
4. **After any reordering or memory-layout change, re-run topology.** Even "harmless" optimizations can break back-neighbors.
