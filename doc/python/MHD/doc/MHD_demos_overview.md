## USER

please look on all files in directory
`/home/prokophapala/git/SimpleSimulationEngine/doc/python/MHD`

I'm a bit confused about thatm (I lost track what is what), perhaps we can split demos, and modules which the demos are using. and documentation. Can you examine it and write in char what is each file good for, what it is doing, a bit of theory, basically some overview of this MHD dictionary

---

Here’s a concise map of everything under `doc/python/MHD`—what each file does, how it fits, and the core theory it touches.

### Core physics/module
- **inductance_core.py** — Kernel library for coaxial circular loops. Provides self/mutual inductance (`calc_self_inductance`, `calc_mutual_inductance`), matrix assembly (`build_inductance_matrix`), flux helpers (`compute_flux`, `init_flux`), flux-conserving solver with optional fixed currents, magnetic energy, Biot–Savart field kernel (`field_loop_rz` from `vector_potential_coil.py`). Used by all demos and plots.

### Plot helpers
- **MHD_plots.py** — Visualization utilities: coil geometry plot with field background/streamlines, axial and radial B profiles, overlays. Depends on `field_loop_rz` from `vector_potential_coil.py`.

### Documentation/notes
- **MHD_Plasma_Nozle_python.md** — Main design doc for the flux‑conserving magnetic nozzle prototype: physics assumptions, inductance matrix, flux conservation, solver structure, plotting, and how [demo_coil_motion_flux.py](cci:7://file:///home/prokophapala/git/SimpleSimulationEngine/doc/python/MHD/demo_coil_motion_flux.py:0:0-0:0) ties in.
- **MHD_GaussGun_demo.md** — Explains the Gauss‑gun accelerator scenario ([demo_gauss_gun_flux.py](cci:7://file:///home/prokophapala/git/SimpleSimulationEngine/doc/python/MHD/demo_gauss_gun_flux.py:0:0-0:0)): driver coil tube, projectile disk loops, active/inactive drivers, flux conservation, force computation.
- **MHD_vectro_potential.md** — Derivation and explanation of vector potential \(A_\phi\) for a circular loop, relation to mutual inductance, and curl‑test to recover \(B\).
- **MHD_check_problem.md** — (Not shown above, but in folder) Likely a checklist/problem statement doc for the MHD work; review for context when needed.

### Demos (flux-conserving circuit solvers)
- **demo_coil_motion_flux.py** — Generic flux‑conserving interpolation: moves each coil linearly from (R0,z0) to (R1,z1), solves \(K(t) I(t) = \Phi_0\) at each step, plots energy/currents and geometry. Can read a coil file (`TYPE R0 z0 R1 z1 I0`) or auto‑generate a parabolic nozzle.
- **demo_plasma_sphere_flux.py** — Builds many plasma surface loops (sphere meridional cross‑section) + seed + cage; reuses `run_coil_motion_flux` to follow expansion of a spherical plasma shell; optional B‑profile overlays.
- **demo_gauss_gun_flux.py** — Gauss‑gun/coil‑gun flux demo: driver coil tube + projectile disk rings. Computes initial flux, progressively deactivates driver coils, solves currents on active set, computes axial force by magnetic‑energy finite differences.
- **demo_plasma_dynamics_flux.py** — Dynamic single‑plasma‑coil model with inertia: seed (+optional cage) fixed; plasma loop moves under Lorentz force and gas‑pressure model while enforcing flux conservation each step.
- **demo_diamagnetic_dipole_bubble.py** — Field‑pattern demo combining a superconducting loop and a point dipole; visualizes Br/Bz, optional separatrix via \( \psi = r A_\phi \); includes consistency checks curl(A) vs B.
- **demo_plasma_dynamics_flux.py** (already listed) and **demo_plasma_sphere_flux.py** cover moving plasma; the former is single‑loop dynamics, the latter many‑loop quasi‑static expansion.

### Tests/checks
- **check_B_kernels.py** — Sanity checks for `field_loop_rz` function: numerical flux vs self‑inductance, on‑axis Bz vs textbook formula.
- **check_B_consistency.py** — Cross‑checks mutual inductance kernel vs integrated B‑field flux; also on‑axis field exactness.

### Vector potential validation
- **vector_potentil_coil.py** — Standalone script deriving/plotting \(A_\phi\) for a loop, computing B via curl numerics, and comparing to analytic `field_loop_rz`; mirrors the math in `MHD_vectro_potential.md`.

How to think about the directory
- **Reusable core:** `inductance_core.py`, `MHD_plots.py`.
- **Docs:** `MHD_Plasma_Nozle_python.md`, `MHD_GaussGun_demo.md`, `MHD_vectro_potential.md`, `MHD_check_problem.md`.
- **Demos:** `demo_*` scripts (nozzle motion, plasma sphere, Gauss gun, dynamic plasma coil, dipole bubble).
- **Checks/validation:** `check_B_*`, `vector_potentil_coil.py`.

If you want, I can propose a clearer split into `core/`, `demos/`, `docs/`, `checks/` folders and add README bullets per file.