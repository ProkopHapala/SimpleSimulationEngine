---
name: test-runner
description: Running/writing tests — pytest classes, markers, ref_data regression, debug output paths
trigger:
  glob:
    - "**/tests/**/*"
    - "**/*test*.py"
    - "**/pytest.ini"
    - "**/conftest.py"
    - "**/doc/TEST_DESIGN.md"
---

## Test System Overview

The test system follows the design in `doc/TEST_DESIGN.md`. Scripts in `tests/` fall into **three classes**:

### Class 1: Pytest Tests (with assertions)

Scripts containing `def test_*` functions with `assert` statements. Run via `pytest`. These are the **primary test suite** — they verify numerical correctness, physical consistency, and regression parity.

- Collected and run by pytest automatically
- Use markers (`slow`, `gpu`) to control which tests run
- May produce debug plots as a side effect, but the pass/fail decision is always an `assert`
- Examples: `test_topology.py`, `test_forcefield.py`, `test_surface.py`, `test_folded_relax.py`, `SPM/test_afm_morse.py`, `SPM/test_afm_fdbm.py`

### Class 2: Standalone Diagnostic Scripts (no assertions)

Scripts run directly via `python tests/<script>.py`. They compute physics, generate plots, and print metrics — but make **no pass/fail assertions**. Their purpose is visual review and manual inspection by the user.

- **Not** collected by pytest (no `def test_*` functions)
- Run on-demand: `python tests/SPM/plot_fdbm_relax.py`
- Output goes to `debug/<script_name>/`
- Examples: `test_folded_surface_scan.py`, `test_tensor_parity.py`, `run_manipulation.py`, `SPM/plot_*.py`, `SPM/run_afm_morse_visual.py`

### Class 3: Utility Modules (no tests, no main)

Helper modules imported by Class 1 and Class 2. Never run directly.

- Examples: `conftest.py`, `helpers/parity.py`, `helpers/geometry.py`, `helpers/scan.py`, `helpers/folded_rigid.py`

### Structure

```
tests/
  conftest.py              # fixtures: data paths, molecule loader, --update-refs option
  test_topology.py          # Class 1: bond/angle/hybridization/type assignment
  test_forcefield.py        # Class 1: UFF/SPFF optimization, NVE conservation
  test_surface.py           # Class 1: Ewald vs brute, GridFF, folded function
  test_folded_relax.py      # Class 1: rigid body relaxation + manipulation on folded basis
  test_lingebra.py          # Class 1: linear algebra eigenvalue tests
  test_integration.py       # Class 1: relaxed scan (stubs — TODO)
  test_folded_surface_scan.py  # Class 2: folded basis fitting + plots
  test_tensor_parity.py     # Class 2: GPU vs CPU tensor kernel parity + plots
  run_manipulation.py       # Class 2: CLI run relaxed scan, export .xyz movie
  SPM/test_afm_morse.py     # Class 1: AFM Morse/LJ force field tests
  SPM/test_afm_fdbm.py      # Class 1: AFM Full Density-Based Model tests
  SPM/plot_density_projection.py  # Class 2: DFTB density visualization
  SPM/plot_fdbm_potentials.py     # Class 2: FDBM potential diagnostics
  SPM/plot_fdbm_relax.py          # Class 2: FDBM relaxation + AFM images
  SPM/run_afm_morse_visual.py     # Class 2: Morse AFM visualization
  ref_data/                 # git-tracked reference files for regression tests
    *.ref.json              # physical properties (z_rel, force, torque, distances)
    *.ref.xyz               # final geometry snapshots
  helpers/
    parity.py               # Class 3: rmse, correlation, overlay_plot, assert_parity
    geometry.py             # Class 3: bond_lengths, angles, planarity, distort
    scan.py                 # Class 3: 1D scan runner (z-scan, x-scan), compare_scans
    folded_rigid.py         # Class 3: rigid body setup, relaxation, manipulation, reference system
```

### Running Tests

```bash
pytest -m "not slow"               # default: fast tests only (< 1s each)
pytest -m "not slow and not gpu"   # fast, no GPU
pytest -m "gpu"                    # GPU tests only
pytest                             # everything (including slow)
pytest --durations=10              # show 10 slowest tests
```

Markers: `slow`, `gpu` (defined in `pytest.ini`).

- `slow` — tests taking > 1s (hard limit 5s); excluded from default runs
- `gpu` — tests requiring OpenCL GPU

### Reference Data System

**Reference files in `tests/ref_data/` are tracked in git** — they are permanent regression test data, not debug artifacts.

- Each test saves physical properties as `{ref_name}.ref.json` and final geometry as `{ref_name}.ref.xyz`
- JSON includes `test_func` and `test_module` fields for bidirectional linking
- Tests compare against references with physical tolerances (not exact float equality)
- Regenerate references after intentional physics changes:
```bash
pytest tests/test_folded_relax.py --update-refs
```
- Always commit updated `.ref.json` and `.ref.xyz` files to git
- See `tests/TEST_RESULTS.md` → "Reference Data System" section for details

### Debug Output Policy

**All diagnostic scripts and visual outputs go to:**

```
$REPO_ROOT/debug/<script_name>/
```

- `<script_name>` = generating script name without `.py`
- Subfolders allowed (e.g., `debug/plot_fdbm_relax/dftb_work/`)
- **Never** write to `/tmp/` or repository root
- Always report the exact paths of generated plots

### Execution Time Policy

**Default test runs must be fast.** The goal is to run as many tests as possible in seconds, not minutes.

- **Each individual test function** (not the whole script) must take **< 1 second** (hard limit: 5 seconds)
- Tests exceeding 1 second should be marked `@pytest.mark.slow`
- Tests exceeding 5 seconds **must** be marked `@pytest.mark.slow` — they are not run by default
- Default run: `pytest -m "not slow"` — fast tests only
- Full run: `pytest` — includes slow tests
- Use `pytest --durations=10` to identify the slowest tests
- When writing new tests, keep inputs small (few atoms, few steps, small grids) to stay under 1 second
- If a test inherently needs to be slow (e.g. convergence study, large molecule relaxation), mark it `@pytest.mark.slow` and document why in the docstring

### Class 2: Standalone Diagnostic Scripts

Standalone scripts (not pytest tests) produce PNG plots and numerical output for **human review**. They are run directly:

```bash
python tests/test_folded_surface_scan.py
python tests/test_tensor_parity.py
python tests/SPM/plot_fdbm_relax.py
python tests/SPM/run_afm_morse_visual.py
```

These scripts:
- Have `if __name__ == '__main__'` blocks, no `def test_*` functions
- Output plots/data to `debug/<script_name>/`
- May use existing plotting functions from `pyBall/SPM/AFM_utils.py` and `pyBall/SPM/AFM.py`
- Are **not** collected by pytest — run them on-demand when visual inspection is needed

### Helpers

- `tests/helpers/parity.py` — `rmse`, `correlation`, `max_err`, `overlay_plot`, `assert_parity`
- `tests/helpers/geometry.py` — `bond_lengths`, `bond_angle`, `planarity`, `distort`, `check_geometry`
- `tests/helpers/scan.py` — `z_scan`, `x_scan`, `compare_scans`
- `tests/helpers/folded_rigid.py` — `setup_rigid_folded`, `relax_folded`, `relaxed_scan`, `plot_relaxed_scan`, `plot_manipulation_trail`, `save_reference`, `compare_to_reference`, `replicate_substrate`

### Consolidation Principle

Before writing new test/plot functions, search existing utility modules. Reuse and generalize rather than reimplement.
