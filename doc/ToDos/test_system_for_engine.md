# TODO: Robust Test System for C++ Engine and Games

## Problem

Current C++ `test_*.cpp` files (47 of them) are interactive SDL/OpenGL apps — visual demos, not tests. No headless mode, no assertions, no regression detection. User must manually build, run, navigate, and judge each one. This is high-friction and unsustainable as the engine grows.

Python tests have a working pytest system (see skill:`test-runner`), but C++ has no equivalent.

## Goal

Build a layered test system with clearly separated levels of user interaction — from fully automatic to interactive GUI — minimizing user burden at each level.

## Test Interaction Layers

### Layer 0 — Fully Automatic (no human)

Headless assertion-based tests. Run on every commit. Exit code 0/1.

- **Conservation laws**: `|ΔE| < tol`, `|ΣF| < tol`, `|ΔL| < tol` after N steps
- **Numerical sanity**: no NaN/Inf, forces bounded, no atom overlap, coordinates in range
- **Topology/index checks**: neighbor list isomorphism, back-neighbor consistency, padding sentinels (see protocol:`topology_verification`)
- **Platform parity**: C++ vs OpenCL vs WebGPU, same inputs → same outputs within tolerance (see skill:`numerical-parity`)
- **Game logic assertions**: score, entity count, checkpoint reached, no clipping through walls
- **Deterministic replay**: record input sequence → replay with fixed seed → assert same outcome

Implementation:
- Add `--headless` mode to every `test_*.cpp` app
- Backend must be headlessly testable (no SDL/OpenGL dependency in test logic) — see `agentic_debugging_principles` §2.2
- `--headless` runs N steps, asserts invariants, dumps state, exits 0/1
- Add C++ invariant helpers to `testUtils.h`: `assertConservation()`, `assertFinite()`, `assertNoOverlap()`, `dumpState()`

### Layer 1 — LLM-Evaluated (no human, or human only on flag)

Automatic checks that require reasoning, not just numeric comparison.

- **Physical intuition**: "Oxygen (negative dipole end) binds to cation?" — LLM reads geometry report and issues verdict (see protocol:`qualitative_validation` §5)
- **Geometric plausibility**: bond lengths in covalent range, angles match VSEPR, ring closure — LLM interprets borderline cases
- **Scan shape assessment**: "does this energy curve have the expected minima/barriers?" — LLM reads numeric summary + sees plot
- **Anomaly detection**: "is this vortex pattern physically reasonable?" — LLM gets statistics + thumbnail

Implementation:
- Test dumps structured report (`.json` + `.md`) with key quantities
- LLM reads report, issues PASS/FAIL/BORDERLINE verdict
- Human reviews only BORDERLINE/FAIL

### Layer 2 — Generated Artifacts for Quick Visual Review (minimal human)

Still images and plots the user can browse in seconds, not minutes.

- **Conservation plots**: energy vs time, momentum drift
- **Trajectory plots**: position vs time, phase space
- **Geometry snapshots**: `.xyz`, `.obj`, `.vtk` rendered to `.png`
- **Field visualizations**: heatmap slices, contour plots
- **Comparison overlays**: reference vs current (see skill:`visual-debug` for plot conventions)

Implementation:
- `--headless` mode auto-generates artifacts to `debug/<test_name>/`
- **Artifact browser**: lightweight HTML gallery or Python viewer that shows all plots from a test run in one window
  - User scrolls through images, clicks OK/Flag per image
  - Generates review verdict file
- **3D viewer**: browser-based (WebGL) or Vispy viewer for `.xyz`/`.obj` snapshots
  - Auto-loads camera angle, atom colors, bonds
  - One command: `python view_snapshot.py debug/test_Mech2D/frame_0500.xyz`
- **State inspector**: export simulation state to standard formats, auto-open in appropriate viewer

### Layer 3 — Interactive GUI Review (human, one-click launch)

Full interactive simulation for cases where still images aren't enough.

- **One-click launch**: `python gui_launcher.py test_Mech2D --t=500`
  - Auto-builds if needed
  - Auto-loads scenario/inputs
  - Auto-positions camera at interesting moment
  - Pauses at step N (or at event trigger)
  - User just looks and presses SPACE to continue or ESC to flag
- **Prepared state**: GUI opens already at the right moment — user does NOT need to:
  - Run the sim from start
  - Load inputs manually
  - Navigate camera
  - Wait for the right frame
- **Review mode**: `--review` flag on test apps
  - Runs sim headlessly to step N
  - Then hands off to interactive SDL loop
  - User can inspect, manipulate, then flag OK/BAD

## Implementation Plan

### Phase 1: Headless C++ Test Harness
- [ ] Add `--headless` mode to 2-3 representative test apps (test_Mech2D, test_SpaceFlightODE, test_Fluid2D)
- [ ] Add invariant check helpers to `testUtils.h`:
  - `assertConservation(E0, E1, tol)` — energy drift
  - `assertFinite(arr, n)` — no NaN/Inf
  - `assertNoOverlap(pos, n, r_min)` — no atom/entity overlap
  - `assertMomentum(p, tol)` — net momentum ~0
  - `dumpState(path)` — export pos/vel/energy to `.csv`/`.xyz`
- [ ] Create `run_all_tests.sh` — runs all `test_* --headless` + pytest, generates `TEST_REPORT.md`

### Phase 2: Artifact Generation + Browser
- [ ] Headless mode auto-generates plots (conservation, trajectory, snapshot)
- [ ] Create `artifact_browser.py` — HTML gallery of all plots from a test run
- [ ] Create `view_snapshot.py` — 3D viewer for `.xyz`/`.obj` (Vispy or WebGL)

### Phase 3: One-Click GUI Review
- [ ] Create `gui_launcher.py` — builds, loads, positions, pauses
- [ ] Add `--review --t=N` flag to test apps
- [ ] Add `--replay input.log` for deterministic input replay

### Phase 4: Aggregated Report
- [ ] `run_all_tests.sh` generates `TEST_REPORT.md`:
  - Summary: X auto-passed, Y auto-failed, Z need visual review
  - Auto-failed: link to artifacts + error details
  - Needs review: link to plots + one-click GUI launch command
  - Auto-passed: compact list, no details

### Phase 5: Game-Specific Patterns
- [ ] Input recording/replay system (serialize SDL events to `.log`)
- [ ] State assertion helpers (score, entity count, checkpoint, collision)
- [ ] Boundary testing (edge of map, extreme input, simultaneous collisions)

## Cross-References

- skill:`test-runner` — Python pytest system, ref_data regression
- skill:`visual-debug` — plot conventions, output paths, shared utilities
- skill:`numerical-parity` — reference vs conservation-law parity, bisect strategy
- protocol:`parity_checking` — Level 2-4 test hierarchy
- protocol:`qualitative_validation` — Level 5 (conservation, symmetry, geometry, intuition)
- protocol:`topology_verification` — structural mapping checks
- `agentic_debugging_principles` — Layer 0-5 test construction, GUI backend/frontend split
- `testUtils.h` — existing C++ test utilities (checkDeriv, printArray, isnan for Vec/Quat)
