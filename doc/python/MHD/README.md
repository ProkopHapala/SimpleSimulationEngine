# MHD — Magnetohydrodynamics & Flux-Conserving Coil Simulations

Axisymmetric (r,z) simulations of magnetic fields, flux conservation, and
plasma-coil interactions for fusion-relevant and pulsed-power scenarios.

## Physics

All simulations assume **axisymmetry** (φ-symmetry), reducing the problem to
the (r,z) meridional plane. The magnetic field is described via:

- **Vector potential** A_φ(r,z) for circular current loops and axial dipoles
- **Flux function** Ψ = 2πr·A_φ (contours of Ψ are field lines)
- **B-field** (B_r, B_z) from analytic elliptic-integral kernels
- **Mutual inductance** M(i,j) between coaxial loops via elliptic integrals

Key physical principles modeled:

- **Flux freezing**: In perfectly conducting loops, Φ_i = ∮B·dA is conserved.
  The inductance matrix K(t)·I(t) = Φ₀ is solved at each geometry step.
- **Diamagnetic plasma**: Modeled as a point dipole whose moment m is solved
  to enforce Ψ_total = 0 at the plasma boundary (separatrix condition).
- **Lorentz force**: F = I·(2πr)·B on plasma/current loops from external fields.
- **Energy gradient force**: F_z = -∂W/∂z for flux-conserving systems.

## Files

### Core library

- **`inductance_core.py`** — Shared library: mutual/self-inductance
  (`calc_mutual_inductance`, `calc_self_inductance`), inductance matrix
  (`build_inductance_matrix`), flux conservation solver
  (`solve_flux_conserving`), B-field kernels (`field_loop_rz`,
  `field_dipole_rz`), vector potential (`Aphi_loop_rz`, `Aphi_dipole_rz`),
  geometry generators (parabolic nozzle, spherical plasma, disk, tube).
- **`MHD_plots.py`** — Plotting utilities: coil geometry visualization
  with HSV/magnitude B-field background, streamlines, 1D B-profile overlays.

### Validation scripts (text output, no plots)

- **`check_B_consistency.py`** — Verifies numerical flux integration of
  `field_loop_rz` matches analytic mutual inductance. Checks on-axis B_z
  against textbook formula. **All checks PASS** (rel err < 1e-7).
- **`check_B_kernels.py`** — Sanity checks: outer-flux vs L_outer formula,
  on-axis B_z scan. On-axis check **PASS** (rel err ~1e-16).

### Demo scripts (generate plots)

- **`demo_dipole_gemini.py`** — Simplest separatrix demo: single seed coil +
  point dipole, solves for dipole moment to create diamagnetic bubble at
  specified radius. Uses Ψ = r·A_φ convention.
- **`demo_dipole_gemini2.py`** — Coupled nozzle+dipole: seed coil + parabolic
  cage rings + moving plasma dipole. Solves [I_cage, m_dipole] simultaneously
  for flux conservation + separatrix condition. Uses Ψ = 2πr·A_φ convention.
- **`demo_diamagnetic_dipole_bubble.py`** — Static field visualization of SC
  coil + axial dipole. Includes curl(A) vs direct B verification and radial
  profile tests.
- **`demo_coil_motion_flux.py`** — Flux-conserving coil motion: all coils
  move linearly from (R0,z0) to (R1,z1), currents solved to conserve initial
  flux. Default: parabolic nozzle with SC seed + cage + plasma armature.
- **`demo_plasma_sphere_flux.py`** — Spherical plasma shell approximated by
  many coaxial loops on sphere surface. Reuses `demo_coil_motion_flux`
  machinery. Includes B-field profile overlays.
- **`demo_plasma_dynamics_flux.py`** — Dynamic single-plasma-coil simulation
  with explicit Euler time integration. Includes gas pressure (adiabatic),
  magnetic self-pressure (hoop stress), and Lorentz force from external coils.
- **`demo_gauss_gun_flux.py`** — Coilgun accelerator: driver tube coils +
  projectile disk. Verlet integration with energy-gradient force. Sequential
  coil switching as projectile passes.
- **`vector_potential_coil.py`** — Vector potential A_φ of a single circular
  loop. Validates B = curl(A) via finite differences vs analytic B-field
  formulas.
- **`fast_eliptke_Abramowitz.py`** — Fast polynomial approximation (Abramowitz
  & Stegun 17.3.34/17.3.36) for elliptic integrals K(m), E(m). Compares
  resulting B-field vs scipy `ellipk`/`ellipe` reference.

## Usage

All demo scripts support `--noshow` to save PNG files instead of opening
a display window:

```bash
cd doc/python/MHD

# Run a demo interactively
python demo_coil_motion_flux.py

# Generate PNG files for review (no display needed)
python demo_coil_motion_flux.py --noshow

# Validation scripts (text output only)
python check_B_consistency.py
python check_B_kernels.py

# Custom parameters
python demo_dipole_gemini.py --r-seed 3.0 --r-plasma 1.0 --noshow
python demo_gauss_gun_flux.py --n-drive 6 --steps 2000 --noshow
```

## Dependencies

- Python 3.10+
- NumPy (>=1.20, compatible with both 1.x and 2.x)
- SciPy (for `ellipk`, `ellipe`)
- Matplotlib

## Code Reusability & Refactoring

All heavy physics and plotting code lives in two shared modules:

- **`inductance_core.py`** — All B-field/A_phi kernels, flux functions (Psi),
  inductance matrices, flux conservation solver, Lorentz force, numerical curl,
  and geometry generators (parabolic nozzle, sphere, disk, tube, driver,
  projectile).
- **`MHD_plots.py`** — `make_background_rgb` (HSV/magnitude field coloring),
  `plot_coil_geometry`, `plot_B_profiles`, and overlay helpers.

Demo scripts are thin CLI wrappers that import from these modules. No script
duplicates kernel or plotting code. The following duplications were removed:

1. `vector_potential_coil.py` had its own `vector_potential_loop_rz`,
   `field_loop_rz`, `compute_numerical_B_from_A` — all identical to
   `inductance_core`. Removed, now imports `Aphi_loop_rz`, `field_loop_rz`,
   `compute_numerical_B_from_A`.
2. `fast_eliptke_Abramowitz.py` had `field_loop_rz_ref` — exact copy of
   `inductance_core.field_loop_rz`. Removed, imports `field_loop_rz as
   field_loop_rz_ref`. Keeps its unique `fast_ellip_ke` +
   `field_loop_rz_fast`.
3. `demo_dipole_gemini.py` had inline `psi_loop`/`psi_dipole` — duplicates of
   `r * Aphi_loop_rz` / `r * Aphi_dipole_rz`. Removed, imports
   `psi_loop_rz`, `psi_dipole_rz` from core.
4. `demo_dipole_gemini2.py` had inline HSV background code (~15 lines) —
   duplicate of `MHD_plots.make_background_rgb`. Replaced with import.
5. `demo_diamagnetic_dipole_bubble.py` had its own `make_background_rgb` —
   duplicate of `MHD_plots.make_background_rgb`. Removed, imports from
   `MHD_plots`.
6. `demo_plasma_dynamics_flux.py` had `compute_force_on_plasma` — duplicate
   of the Lorentz force loop. Removed, imports `compute_lorentz_force_loop`
   from core.
7. `demo_gauss_gun_flux.py` had `build_driver_tube`, `build_projectile_disk`,
   `build_geometry_vectors` — geometry builders moved to `inductance_core`.
   Removed, imports from core.

## Bugs Fixed During Review

1. **`fast_eliptke_Abramowitz.py`**: Factor 2π error in `field_loop_rz_fast`
   and `field_loop_rz_ref` — used `4e-7·π·I` (= μ₀I) instead of
   `μ₀I/(2π)`. Fixed to use `MU0 * I / (2π·√denom)`.
2. **`check_B_consistency.py`**: Used `np.trapz` (deprecated in NumPy 2.0)
   / `np.trapezoid` (absent in NumPy <2.0). Replaced with compatibility
   shim: `_trapezoid = getattr(np, 'trapezoid', None) or np.trapz`.
3. **`check_B_kernels.py`**: Same `np.trapezoid` compatibility fix applied.
4. **`demo_coil_motion_flux.py`**: Double `plt.show()` (one inside
   `run_coil_motion_flux`, one at end of `__main__`) caused blocking when
   called from `demo_plasma_sphere_flux.py`. Removed trailing `plt.show()`,
   added `noshow` parameter.
5. **`vector_potentil_coil.py`** → renamed to **`vector_potential_coil.py`**
   (typo fix per audit recommendation).
6. **`demo_dipole_gemini.py`**, **`vector_potential_coil.py`**,
   **`fast_eliptke_Abramowitz.py`**: Added `argparse` + `if __name__ ==
   "__main__"` guard (previously had hardcoded parameters at module level).

## Known Limitations

- `demo_plasma_dynamics_flux.py` uses explicit Euler (1st order) — energy
  drift expected for long runs. Consider Verlet or RK4 for production use.
- `check_B_kernels.py` outer-flux test shows large discrepancy (rel err ~37)
  because the thin-wire self-inductance formula L_outer = μ₀a·ln(8a/r_core)
  is a logarithmic approximation, not exact for the numerical flux integral.
- `tmp.md` contains draft notes with embedded code snippets — should be
  formalized or removed.
- Psi convention differs between `demo_dipole_gemini.py` (Ψ=r·A_φ) and
  `demo_dipole_gemini2.py` (Ψ=2πr·A_φ). Both are valid for separatrix
  detection (Ψ=0 contour), but mixing them would cause scaling errors.
