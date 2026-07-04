---
name: numerical-parity
description: Comparing two calculations — reference vs conservation-law parity, bisect strategy, float32/64 tolerance
trigger:
  glob:
    - "**/tests/**/*"
    - "**/*test*.py"
    - "**/*parity*.py"
    - "**/*check*.py"
    - "**/*debug*.py"
    - "**/test_*.sh"
    - "**/run_*.sh"
    - "**/*benchmark*.py"
---

## Protocol

**Define correctness before coding:** Verify via parity checks against reference code, known analytical solutions, physical conservation laws, symmetry checks, or known physical limits.

**Two modes of parity:**
1. **Reference parity** — you have a trusted implementation (C++, Fortran). Match inputs, compare outputs component by component.
2. **Conservation-law parity** — no reference exists. Verify physical invariants instead:
   - `|ΔP| < tol` after one step with no external forces
   - `|ΔL| < tol` for rotation-capable systems
   - Energy bounded (MD) or monotonically decreasing (relaxation)
   - Symmetry preserved (e.g. mirror symmetry of a symmetric molecule)

1. **Input-first**
   - Compare buffers before physics: topology, parameters, indexing.
   - Check neighbor lists, bond connectivity, stiffness (k), equilibrium values (l0).
   - Watch for 0/1-based indexing, struct padding, sorting differences.

2. **Component isolation**
   - Disable ALL forces. Test one component at a time: bonds → angles → dihedrals → nonbonded.
   - Use switches to gate. Verify exclusion masks (1-2, 1-3, 1-4).

3. **Synchronized tracing**
   - If inputs match but outputs differ, inject identical printf format in both implementations.
   - Gate by atom ID and step to avoid flooding (e.g., `if(id==0 && step==0)`).

4. **Single-step**
   - Run exactly ONE iteration. Compare forces (F) first, then energies (E).

5. **Tolerance discipline**
   - GPU float32 vs CPU float64: expect 1e-6 for bonds/angles, up to 5e-5 for trig-heavy parts.
   - Don't chase post-parity noise.

## Bisect Strategy

When outputs differ, drill down systematically:
1. **Input & Topology Check:** Verify buffers *before* physics loops.
2. **Component Isolation:** Disable all forces via switches. Test one component at a time.
3. **Block-Level Parity:** Compare individual matrix blocks per neighbor slot.
4. **Dense Reassembly:** Rebuild dense matrices from blocks and compare to reference.
5. **Single-Step End-to-End Execution:** Run exactly ONE iteration. Compare forces first, then energies.

Stop and fix issues at each level before moving to the next.

## Structural vs. Numeric Divergence

Keep structural checks (indices, counts, topology maps) separate from floating-point parity:
- **Structural:** Expect exact integer matching and true zeros.
- **Numeric:** Expect Fp32 vs Fp64 noise. Do not chase noise past tolerance limits.

## Discrepancy Tracking

- **Range Sanity:** Place strategic checks throughout calculations to ensure values are within physical boundaries and are not `NaN`, infinity, or unexpected zeros.
- **Worst-Diff Tracking:** When debugging mismatched arrays, scan the entire buffer and print out the **worst-diff location (index) and value** to pinpoint exactly where the divergence starts.
