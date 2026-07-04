---
name: visual-debug
description: Diagnostic plots/visualizations — shared utilities, output paths, plot style conventions
trigger:
  glob:
    - "**/tests/**/*"
    - "**/*test*.py"
    - "**/*debug*.py"
    - "**/test_*.sh"
    - "**/run_*.sh"
    - "**/*benchmark*.py"
---

## Shared Utilities

Before writing ad-hoc debugging/plotting code, check these existing modules:

**Python Visualization:**
- `pyBall/plotUtils.py` - matplotlib utilities: 1D/2D function plots, geometry visualization, field slices, scan profiles, derivative plots
- `pyBall/VispyUtils.py` - GPU-accelerated 3D visualization: AtomScene class for interactive molecular viewing, bond visualization, force vectors

**Python Testing/Diagnostics:**
- `pyBall/DFTB/TestUtils.py` - RMS error computation, checkpoint management (save/load/compare), grid generation, eigenvec printing
- `pyBall/atomicUtils.py` - Atomic utilities: normalize, findAllBonds, graph preprocessing, adjacency lists

**C++ Testing/Diagnostics:**
- `cpp/common/testUtils.h` - Print arrays/vectors/matrices, compareVecs, derivative checking (checkDeriv, checkDeriv3d), timing (StopWatch), error macros (TEST_ERROR_PROC_N, SPEED_TEST_FUNC)

## Test Artifacts

- **Output location policy:** All diagnostic scripts and visual test outputs must be saved under `debug/<script_name>/` (e.g., `debug/plot_fdbm_relax/`). The `<script_name>` is the generating script's name without `.py`. Subfolders are allowed for organization. **Never** write to `/tmp/` or the repository root. This keeps artifacts persistent and easy to find.
- **Structured outputs:** Group all debugging, benchmarking, and testing outputs into organized, numbered directories (e.g., `tests/003_case_name/`). Do not clutter root directories. Explicitly report their location.
- **Foreground execution:** Run tests synchronously in the foreground with full output. Never hide output or use background commands (`&`, `| tail`, `| head`, or silent redirects). Full `stdout` must be visible.

## Visual Review

- **Python tests:** Generate diagnostic plots using `matplotlib` saved as `.png` files. Use shared helpers like `plotUtils.py` (e.g., `plot_scan_profile`, `plot_field_slice`, `plotGeometryWithForces`).
- **Optional plotting:** Make plotting optional via flags (e.g., `--noPlot`, `--saveFig`). Isolate `plt.show()` strictly to the CLI/main entry point.
- **Report paths:** Always report the exact paths/folders of generated plots.

## Diagnostics

- **Numerical range sanity:** Strategically place checks throughout calculations to ensure values are within reasonable limits and are not `NaN`, infinity, or unexpected zeros.
- **Checkpointing:** Use `pyBall/DFTB/TestUtils.py` checkpoint functions (`save_checkpoint`, `load_checkpoint`, `compare_checkpoint`) for parity testing and reproducible debugging.
- **RMS error:** Use `compute_rms_error` from `TestUtils.py` for array comparisons.

## Consolidation Principle

- **Reuse over reinvent:** Before writing new debug/plot/test functions, search existing utility modules. Generalize existing functions if they almost fit your needs.
- **Separate concerns:** Keep compute algorithms separate from plotting/diagnostics. Move ad-hoc plotting code from test scripts into shared utilities.
- **Zero-copy buffers:** For Python-C++ interop, use `np.ctypeslib.as_array` pattern (see `python_native_bindings` skill) instead of copying data.

## Plot Style Preferences

When creating diagnostic comparison plots (e.g., reference vs model curves):

- **Subtract constant offset (DC term):** Methods like Ewald summation cannot capture the constant (DC/zero-frequency) term, causing a constant shift between model and reference. Compute `dc = mean(model - ref)` and subtract it from the model so curves overlap. Center reference by subtracting its own mean for display. The residual difference after this shift reveals true numerical error.
- **Reference curve:** `ls=':'`, `lw=1.5` (dotted, thick — clearly visible as the reference)
- **Model curve(s):** `ls='-'`, `lw=0.5` (solid, thin — the thing being tested)
- **Residual difference on twin axis:** Plot `(model_shifted - ref) * 100` on a secondary y-axis (twinx) so small residual errors are visible. Label it `diff x100`.
- **RMSE/MaxErr box:** Show RMSE and max error in a text box (upper-left, monospace, semi-transparent background).
- **General rule:** The reference should always be more visually prominent (thicker, dotted) than the model being compared.
