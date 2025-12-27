# Gauss-Gun Coil Accelerator – Axisymmetric Flux/Force Demo

(This document describes the Gauss-gun style demo implemented in `demo_gauss_gun_flux.py`. It reuses the inductance/field machinery from `inductance_core.py` and follows the conventions of `MHD_Plasma_Nozle_python.md`.)

## 1. Scenario Overview

We model an axisymmetric **coil gun / Gauss gun** in cylindrical `(r,z)` coordinates:

- A tube of **driver coils** (superconducting or powered) along the `z` axis.
- A **metallic projectile disk**, represented as several coaxial loops at a common axial position.
- Driver coils are initially **pre-magnetized**; as the projectile moves, we **switch off** coils slightly ahead of it.
- Switching off a coil rapidly changes the external flux through the projectile, inducing currents in the disk loops.
- The projectile experiences an axial magnetic force and is accelerated along `+z`.

The demo keeps full axisymmetry and uses the same inductance matrix and Biot–Savart kernels as the plasma nozzle prototypes, but now includes explicit **dynamics of the projectile** (1D motion along `z`).


## 2. Geometry and Coils

### 2.1 Coordinate system

- Cylindrical coordinates `(r,z)` with symmetry axis at `r = 0`.
- All coils are circular loops centered on the axis.
- The projectile is constrained to move along `z`; its radial structure is fixed.

### 2.2 Driver coil tube

We model `N_drive` identical (or parameterized) driver coils arranged along `z`:

- Coil radii: `R_drive[k]` (often a constant `R_drive`).
- Axial positions: `Z_drive[k]`, typically evenly spaced between `z_start` and `z_end`.
- Currents: `I_drive[k]`.

Interpretation in the current implementation:

- At `t = 0` all driver coils are part of the **flux-conserving set**, together with the projectile rings.
- As the projectile moves, we gradually **disconnect** driver coils in front of it. A disconnected coil is:
  - Removed from the active flux-conserving set.
  - Held at `I = 0` afterwards (it no longer cares about its initial flux).
- Remaining active driver coils and the projectile loops continue to conserve their initial fluxes.

We track an `active_drive[k]` boolean mask indicating which driver coils still belong to the active flux-conserving set. Geometry is kept directly in arrays `(R_drive, Z_drive)` without using the textual `TYPE R0 Z0 R1 Z1 I0` format.

### 2.3 Projectile disk

The projectile is a thin conducting **disk** modeled as `N_disk` coaxial loops:

- Radii: `r_disk[i]` between `r_inner` and `r_outer`.
- Axial positions: `z_disk[i](t) = z_proj(t)` for all `i`.

Currents in projectile rings `I_disk[i](t)` are **unknowns** solved from flux conservation conditions (see below); they are not externally driven.

In code, projectile loops occupy indices `N_drive .. N_drive + N_disk - 1` in the global loop arrays.


## 3. Inductance Matrix and Flux Constraints

We reuse `inductance_core.py`:

- `build_inductance_matrix(rs, zs)`
- `calc_mutual_inductance`, `calc_self_inductance`
- `magnetic_energy(K, I)`

At any time step we construct arrays:

- `rs[0:N_drive]  = R_drive[:]`
- `zs[0:N_drive]  = Z_drive[:]`
- `rs[N_drive:]   = r_disk[:]`
- `zs[N_drive:]   = z_proj * np.ones(N_disk)`

Then:

```python
K = build_inductance_matrix(rs, zs)
```

### 3.1 Initial flux for all loops

At `t = 0`:

- Driver coil currents are set to a chosen `I_drive0` profile (typically all equal and non-zero).
- Projectile disk currents start from `0`.

We build the initial inductance matrix `K0` at `(z_proj0, Z_drive)` and currents

```python
I0 = np.concatenate([I_drive0, np.zeros(N_disk)])
Phi0 = K0 @ I0
```

Here `Phi0[i]` is the **initial flux** linked to loop `i` (drivers and projectile alike).

Physical interpretation:

- Initially, all loops are treated as perfectly conducting and flux-conserving.
- As the projectile moves and we disconnect driver coils, only the **active subset** of loops (remaining drivers + all projectile rings) continues to satisfy `Phi_active(t) = Phi0_active`.

### 3.2 Solving for currents in the active flux-conserving set

At a later time step with projectile position `z_proj(t)` we rebuild the inductance matrix for all loops at the current geometry, but solve only on the **active** set:

1. Form global geometry vectors for all loops

   ```python
   rs_full, zs_full = build_geometry_vectors(R_drive, Z_drive, r_disk, z_proj)
   K_full = build_inductance_matrix(rs_full, zs_full)
   ```

2. Define active indices:

   ```python
   idx_drive_active = np.where(active_drive)[0]
   idx_proj = np.arange(N_drive, N_total)
   idx_active = np.concatenate([idx_drive_active, idx_proj])
   ```

3. Solve the reduced system on the active subset:

   ```python
   K_act   = K_full[np.ix_(idx_active, idx_active)]
   Phi0_act = Phi0[idx_active]
   I_act   = np.linalg.solve(K_act, Phi0_act)
   ```

4. Scatter the solution back into a full current vector (inactive drivers set to zero):

   ```python
   I_full = np.zeros(N_total)
   I_full[idx_active] = I_act
   ```

This is equivalent in spirit to using `solve_flux_conserving` with a shrinking active set, but implemented explicitly for clarity and to make it easy to completely remove disconnected coils from the flux-conserving matrix.


## 4. Coil Switching Logic

To emulate a Gauss gun we implement a simple rule:

> At each time step, find the nearest **active** driver coil ahead of the projectile along `+z`. If its axial distance from the projectile is less than a threshold `a_switch`, mark that coil as disconnected. Disconnected coils are removed from the flux-conserving system and their currents remain zero afterwards.

Algorithm sketch:

1. Given `z_proj` and arrays `Z_drive[k]`, `active_drive[k]`:
   - Consider all `k` with `Z_drive[k] > z_proj` and `active_drive[k] == True`.
   - If none exist, no switching this step.
2. Find the closest such coil in `z`:

   ```text
   dz_k = Z_drive[k] - z_proj > 0
   k_next = argmin(dz_k)
   ```

3. If `0 < (Z_drive[k_next] - z_proj) < a_switch`, then mark it as disconnected:

   ```python
   active_drive[k_next] = False
   # coil k_next is then excluded from the active index set
   # and its current is forced to zero in I_full
   ```

This provides a simple, controllable switching pattern that follows the projectile as it advances while gradually shrinking the active flux-conserving matrix.


## 5. Force Model and Dynamics

We want the axial force `Fz` on the projectile and we integrate its motion along `z` using a simple **Verlet scheme**.

### 5.1 Magnetic force from energy gradient

Total magnetic energy of the system at a given configuration is

```python
W = magnetic_energy(K, I)
```

For a conservative magnetic system, the generalized force along a coordinate `q` is

```text
F_q = -∂W/∂q
```

For the projectile coordinate `q = z_proj` we approximate this derivative by a symmetric finite difference:

```python
Fz ≈ - (W_plus - W_minus) / (2 * dz_fd)
```

where

- `W_minus` is computed at `z_proj - dz_fd`.
- `W_plus` is computed at `z_proj + dz_fd`.

At each of these offset positions we recompute:

1. Projectile geometry (`zs` for projectile loops).
2. Inductance matrix `K`.
3. Projectile currents `I_disk` from the same flux constraint `Phi0_disk` and the **same** driver currents `I_drive` (no further switching inside the finite-difference evaluation).
4. Total energy `W = 0.5 * I^T K I`.

This yields a robust and relatively simple `Fz` consistent with the inductance formulation.

### 5.2 Equation of motion and Verlet integration

We treat the projectile as a rigid body of mass `m_proj` with position `z(t)` and velocity `v(t)`.

We use the standard position-Verlet update:

```text
z_{n+1} = 2 z_n - z_{n-1} + (F_n / m_proj) * dt^2
```

Implementation outline:

1. Initialize:
   - `z_0` (initial axial position of projectile).
   - `v_0` (initial velocity, often 0).
2. Compute the initial force `F_0` at `z_0`.
3. Set a “previous” position consistent with the initial velocity and force, e.g.

   ```text
   z_{-1} = z_0 - v_0 * dt + 0.5 * (F_0 / m_proj) * dt^2
   ```

4. For each time step `n = 0 .. N_steps-1`:
   - From `z_n` and the current active set, apply coil switching logic.
   - Solve for currents in all **active** loops via flux constraints (drivers + projectile).
   - Compute `F_n` via the energy finite-difference procedure.
   - Update `z_{n+1}` from the Verlet formula.
   - Optionally reconstruct velocity as

     ```text
     v_n ≈ (z_{n+1} - z_{n-1}) / (2 * dt)
     ```

The simulation stops when the projectile exits the driver tube (`z_proj > z_end + margin`) or a maximum number of steps is reached.


## 6. Implementation: `demo_gauss_gun_flux.py`

The script `demo_gauss_gun_flux.py` (in `doc/python/MHD`) implements the above model. It:

1. Parses CLI arguments to define:
   - Driver tube geometry (`n_drive`, `r_drive`, `z_start`, `z_end`).
   - Driver current (`i_drive`).
   - Projectile disk geometry (`n_disk`, `r_inner`, `r_outer`).
   - Projectile mass and initial state (`m_proj`, `z0`, `v0`).
   - Switching distance (`a_switch`).
   - Time stepping parameters (`dt`, `steps`).
2. Generates driver coil positions and projectile disk radii.
3. Builds `K0`, computes initial flux `Phi0`, and extracts `Phi0_disk`.
4. Steps forward in time with the Verlet integrator, updating:
   - The active set of driver coils based on switching logic.
   - Currents in all active loops from flux constraints.
   - Projectile position from the force.
5. Records time series of:
   - `t`, `z(t)`, `v(t)`, `Fz(t)`.
   - Magnetic energy `W(t)`.
   - Projectile kinetic energy `K(t) = 0.5 m v^2`.
   - Total energy `E_tot(t) = W(t) + K(t)`.
   - Currents in each driver coil and projectile ring.
6. Produces diagnostic plots for trajectory, forces, energies, and coil currents, plus a geometry/B-field snapshot using the shared `MHD_plots.plot_coil_geometry` utility.


## 7. CLI Sketch and Typical Run

A simple default run from `doc/python/MHD` is:

```bash
python demo_gauss_gun_flux.py
```

which (with current defaults) uses for example:

- `n_drive = 4`, `r_drive ≈ 1.0`, `z_start = 0.0`, `z_end = 2.0`.
- `n_disk = 1`, projectile radius ≈ `0.9 * r_drive`.
- `m_proj = 0.1`, `z0 = -0.05`, `v0 = 5.0`.
- `a_switch ≈ 0.08`, `dt ~ 5e-6`, `steps ~ 4000`.

The script prints basic diagnostics (number and positions of switch events, `z` range) and opens plots showing projectile motion, energies, and coil currents, plus a geometry/B-field snapshot.

You can of course override any of these via the CLI, e.g. to increase tube length or number of driver coils:

```bash
python demo_gauss_gun_flux.py \
  --n-drive 8 --r-drive 1.0 --z-start 0.0 --z-end 4.0 --i-drive 1e6 \
  --n-disk 3 --r-inner 0.5 --r-outer 0.9 --m-proj 0.1 \
  --z0 -0.05 --v0 5.0 --a-switch 0.2 --dt 5e-6 --steps 8000
```


## 8. Insights and Possible Extensions

Some observations from running the current demo:

- With all loops initially flux-conserving, the system behaves like a coupled electromagnetic–mechanical oscillator. **Total energy** `E_tot(t)` is approximately conserved (within numerical error), while magnetic and kinetic energy exchange as the projectile is accelerated and then interacts with the remaining coils.
- Switching driver coils off by **removing them from the active flux-conserving set** produces clear, discrete changes in the field structure and current distribution, visible both in the current time series and in the B-field geometry plot.
- The projectile’s final speed and the shape of the current pulses are sensitive to
  - the tube geometry (coil spacing, radius, number of coils),
  - the switching distance `a_switch`, and
  - the initial seeding current `I_drive0`.

Future refinements that fit into the same framework:

Future refinements that fit into the same framework:

- Smooth (time-ramped) coil switching instead of instant cutoff.
- Multiple projectile segments or a short “train” of disks.
- Inclusion of resistive damping in projectile loops.
- Using `plot_B_profiles` to overlay 1D axial/radial B-field profiles on top of the geometry.
