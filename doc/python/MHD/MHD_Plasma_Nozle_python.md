# MHD Plasma Nozzle – Python Prototype

High‑level documentation for the axisymmetric MHD / coil‑circuit prototype implemented in this directory. This document ties together two core components:

- `inductance_core.py` – core inductance and B‑field kernels
- `demo_coil_motion_flux.py` – a flux‑conserving coil‑motion demo

The goal is to document the physics assumptions, numerical choices, and practical issues encountered while building a minimal but robust Python prototype of a flux‑conserving magnetic nozzle, intended as a reference and sanity check for future C++ / WebGL implementations.


## 1. Physical Scenario and Design Intent

The target physical system is an axisymmetric magnetic nozzle for pulsed nuclear propulsion:

- **Superconducting (SC) seed coil** provides an initial bias magnetic flux. It sits outside and behind the hot region, shielded from radiation.
- **Normal‑metal cage** (e.g. tungsten / tantalum rings) forms a hot nozzle interior. During a pulse the cage is highly conducting on the microsecond timescale and can support large induced currents.
- **Diamagnetic plasma ball** is born inside the nozzle and rapidly expands, pushing magnetic flux outwards against the cage.

Engineering objectives:

- Protect the SC coil from current spikes and direct radiation.
- Drive most of the transient current in the cage and, secondarily, in the outer plasma shell.
- Let plasma expansion do mechanical work against magnetic pressure, increasing magnetic energy (the system behaves as an electric generator).

Simulation objectives for this Python prototype:

- Work in cylindrical coordinates `(r, z)` with coaxial circular loops.
- Treat coils and plasma “filaments” as 1D rings; no volumetric MHD grid.
- Use **flux conservation in closed loops** (Method 1 in the design doc) rather than control‑point fitting (Method 2), to get a clean reference
  solution.
- Provide tools to:
  - Build inductance matrices from geometry
  - Compute fluxes and magnetic energy
  - Evolve a set of loops between two geometries under strict
    per‑loop flux conservation
  - Visualize geometry, fields, currents and energies for inspection

This is deliberately a **quasi‑static** prototype: there is no time stepping of plasma dynamics here, only initial vs final geometries with flux‑based current redistribution. Full dynamics (forces, gas pressure, springs, remeshing, instabilities) are reserved for later.


## 2. Core Physics: Inductance Matrix and Flux Conservation

The central abstraction is the inductance matrix `K` relating loop currents `I` to loop fluxes `Phi`:

```python
Phi = K @ I
```

For `N` coaxial circular loops with radii `r[i]` and axial positions `z[i]`, elements of `K` are:

- **Diagonal** – self‑inductance of a single loop

  ```text
  L_i = 4e-7 * pi * R * (ln(8 R / r_wire) - 1.75)
  ```

  Thin circular loop, with a core radius `r_wire` used as a   regularization length.

- **Off‑diagonal** – mutual inductance between loops `i` and `j`, from   standard elliptic‑integral formulas for coaxial circles:

  ```text
  M_ij = M(r_i, z_i; r_j, z_j)
  ```

where `M` is implemented via `scipy.special.ellipk` / `ellipe`. See `calc_mutual_inductance` for the exact kernel.

### Flux conservation assumption

For a fast process (microseconds) in highly conducting loops and plasma, resistive effects are negligible. The relevant invariant is **magnetic flux** through each closed loop, not current. For loop `i` we enforce

```text
Phi_i(t1) = Phi_i(t0)
```

Even if the loop geometry changes between `t0` and `t1`, and even if other loops move or appear, the total flux linked to that loop is conserved as long as no resistive dissipation or reconnection is modeled.

In matrix form, if `K0` is the inductance matrix at `t0`, `I0` the initial currents, and `K1` the matrix at `t1`, we have

```python
Phi0 = K0 @ I0    # initial total flux per loop
K1 @ I1 = Phi0    # flux conservation constraint at t1
```

The unknowns are the currents `I1` at the final geometry. Solving this system (or a reduced version with some fixed‑current loops) gives a clean, physically interpretable answer: how the system redistributes currents to preserve per‑loop flux during geometry change.


## 3. Module `inductance_core.py`

This module holds the core physics building blocks used by all demos. It is designed to be small, testable, and reusable from C++/JS ports.

### Key functions

- **`calc_mutual_inductance(r1, z1, r2, z2)`**

  Mutual inductance for two coaxial circular loops using elliptic integrals.   Handles degenerate cases `r1 == 0` or `r2 == 0` by returning `0`.

- **`calc_self_inductance(R, r_wire=0.01)`**

  Self‑inductance of a thin loop of radius `R` with a specified wire radius. The constant `1.75` is the usual empirical correction for thin circular loops.

- **`build_inductance_matrix(rs, zs, self_r_wire=0.01)`**

  Assembles the full `N x N` inductance matrix `K` for a set of loops from vectors of radii and axial positions. This is the primary entry point for building `K0` and `K(t)` in the demos.

- **`compute_flux(K, I)`**

  Simple helper: `Phi = K @ I`. Used to compute initial fluxes and to verify flux conservation in diagnostics.

- **`init_flux(K0, I0, flux_override=None)`**

  Computes initial fluxes from `K0` and `I0` and optionally allows manual overrides per loop (e.g. for strictly diamagnetic loops where one wants to force `Phi0[i] = 0` rather than use the current‑generated flux). Currently the main demo uses default behavior without overrides.

- **`solve_flux_conserving(K1, Phi0, fixed_mask=None, I_fixed=None)`**

  General solver for flux‑conserving problems with optional fixed‑current loops. In the simplest case (all loops flux‑conserving), this reduces to

  ```python
  I1 = np.linalg.solve(K1, Phi0)
  ```

  For mixed systems (some loops fixed‑current, some flux‑conserving), this partitions `K1` and solves only for the unknown subset.

- **`magnetic_energy(K, I)`**

  Computes total magnetic energy

  ```text
  W = 0.5 * I^T K I
  ```

  including all mutual terms. This is the correct scalar energy of the coupled system. Note that when comparing “energy of coil A vs B”, the presence of mutual terms makes interpretation subtle; the demo focuses on total energy and currents rather than attempting to partition `W`.

- **`generate_parabolic_nozzle(...)`**

  Utility to generate a *textual* description of a simple parabolic cage + plasma + SC configuration, compatible with the `demo_coil_motion_flux` parser. Output is a multi‑line string with columns

  ```text
  TYPE  R0  Z0  R1  Z1  I0
  ```

  where `TYPE` is one of `SC`, `CAGE`, `PLASMA`. Geometry is defined by throat radius, exit radius, axial extent, and optional overrides for SC and plasma paths.

- **`field_loop_rz(a, z0, I, r_grid, z_grid, eps=1e-12)`**

  Biot–Savart magnetic field kernel for a circular loop at radius `a` and  axial position `z0`, evaluated on a 2D `(r,z)` grid. Returns the radial and axial components `(Br, Bz)` using the standard elliptic‑integral expressions. A small regularization `eps` is used near the axis to avoid division by zero; by symmetry `Br` is forced to zero on the axis.

This kernel is used solely for visualization in the Python demo but shares its mathematical structure with the mutual‑inductance formulas; agreement between visual field patterns and flux conservation is therefore a strong sanity check.


## 4. Script `demo_coil_motion_flux.py`

This script is the main **numerical experiment** in the Python prototype. It takes one or more coils defined by

```text
TYPE  R0  Z0  R1  Z1  I0
```

and simulates a linear motion from initial positions `(R0, Z0)` to final positions `(R1, Z1)` while conserving per‑loop flux.

### Parsing and geometry

- Coil definitions are read either from a text file or from the generated  parabolic nozzle string. The parser `parse_coils` accepts a 2D string array (as produced by `np.genfromtxt(..., dtype=str)`) and returns

  ```text
  types, R0, z0, R1, z1, I0
  ```

- Special types:
  - `SC` – seed superconducting coil (non‑special numerically; the label
    is used only for coloring and legend text).
  - `CAGE` – cage rings.
  - `PLASMA` – plasma loop.

Currently all loops are treated as flux‑conserving; implementing truly fixed‑current SC coils would be a small modification using
`solve_flux_conserving`.

### Flux‑conserving evolution

The main driver `run_coil_motion_flux` implements the following steps:

1. Parse coils → `types, R0, z0, R1, z1, I0`.
2. Build initial inductance matrix `K0` at `(R0, z0)`.
3. Compute initial flux vector

   ```python
   Phi0 = K0 @ I0
   ```

4. Sample an interpolation parameter `t` from `0 → 1` in `n_steps` steps.
5. For each `t`:
   - Compute interpolated geometry

     ```python
     rs = (1 - t) * R0 + t * R1
     zs = (1 - t) * z0 + t * z1
     ```

   - Build inductance matrix `K(t)` for `rs, zs`.
   - Solve `K(t) @ I(t) = Phi0` via `np.linalg.solve`.
   - Compute total magnetic energy `W(t)`.

6. Store `I(t)` and `W(t)` for diagnostic plots.

### Flux conservation check

At the end of the run, the script explicitly recomputes the inductance
matrix at the final geometry and prints a per‑coil flux comparison:

- `Phi0[i]` – initial flux of coil `i`
- `Phi_final[i]` – final flux computed from `K_final @ I_final`
- relative error

This allows direct verification that the numerical solve indeed preserved
per‑coil flux to machine precision. For well‑conditioned configurations the
relative error is typically at the level of floating‑point noise.

### Diagnostics and plotting

The script can produce two kinds of plots:

- **Time series plots** (vs interpolation parameter `t`):
  - Total magnetic energy `W(t)`
  - Per‑coil currents `I_i(t)` (in MA)

  These plots are essential for understanding whether the expanding plasma
  configuration behaves as an **energy generator** (`W` increases when the
  plasma does work pushing against the magnetic field) or as a sink.

- **Geometry + field plots** (initial vs final):

  The function `plot_geometry` produces a two‑panel figure:

  - Left: initial configuration (`t = 0`)
  - Right: final configuration (`t = 1`)

  Both panels show:

  - Coil positions in the `(z, r)` plane, mirrored to `±r` to emphasize
    axisymmetry.
  - Coil markers:
    - Color by type (`SC`, `CAGE`, `PLASMA`)
    - Small circular outlines plus central dots for visibility
    - Text labels with current in MA.
  - A background magnetic field magnitude map on a regular `(z, r)` grid,
    computed by summing contributions from all loops via `field_loop_rz`.
  - Streamlines (`streamplot`) of the `(Br, Bz)` field to visualize
    qualitative field line structure.

  Two background modes are supported:

  - `bg_mode="hsv"` – HSV colormap with hue encoding field direction
    (angle in `Br–Bz` plane) and value encoding magnitude.
  - `bg_mode="mag"` – magnitude‑only view (`magma` colormap) where
    bright regions correspond to large `|B|`.

  To make weaker structures visible next to a dominant seed coil, the
  background normalization deliberately uses a small reference scale

  ```text
  B_ref = 0.01 * max(|B|)
  ```

  so gradients in weaker regions are not completely washed out.

### CLI entry point

The bottom of `demo_coil_motion_flux.py` implements a small `argparse`
CLI that supports two usage modes:

1. **External coil definition file**

   ```bash
   python demo_coil_motion_flux.py coils.dat
   ```

   where `coils.dat` contains one coil per line with columns

   ```text
   TYPE  R0  Z0  R1  Z1  I0
   ```

2. **Generated parabolic nozzle** (no filename argument):

   ```bash
   python demo_coil_motion_flux.py
   ```

   In this case the script calls `generate_parabolic_nozzle` with parameters set by optional CLI flags (throat/exit radius, axial extent, SC current and position, plasma start/end positions).

Additional flags control the number of interpolation steps, whether to plot energy/current time series and/or the geometry, and which background B‑field mode to use.


## 5. Issues Encountered and How They Were Solved

During development several conceptual and technical issues arose. The key ones and their resolutions are summarized here.

### 5.1 Flux vs current conservation

Initial intuition phrased invariants in terms of “conserved currents”. Careful analysis clarified that for ideal conductors and superconductors in fast transients it is **flux** that is conserved, not current. The correct
constraint is

```text
Phi_i = const
```

not `I_i = const`.

**Resolution:**

- Implemented the inductance‑matrix formulation

  ```python
  Phi = K @ I
  ```

  and explicitly enforced `K1 @ I1 = Phi0` for moving loops.
- Verified that with this approach, in simple 2‑coil and 3‑coil toy  geometries, currents behave in a physically reasonable way (e.g. induced currents oppose flux change, cage shields SC, etc.).

### 5.2 Magnetic energy apparently decreasing

Early toy demos appeared to show total magnetic energy *decreasing* as the plasma loop expanded, which contradicted the expectation that plasma expansion should *do work* on the magnetic field and therefore increase `W`.

**Key insights:**

- Energy accounting must include mutual inductances and use the full quadratic form `0.5 * I^T K I`.
- The reference case for validating monotonicity is a controlled geometry where intuition is clear (e.g. plasma loop expanding inside a fixed SC loop, or a concentric cage configuration).

**Resolution:**

- Switched to the matrix‑based flux‑conserving formulation and recomputed energies consistently from `K(t)` and `I(t)`.
- Constructed simpler, more controlled geometries (2‑coil, 3‑coil with cage) and verified that energy changes align with physical expectations in those cases.

### 5.3 Geometry and plotting clarity

Initial geometry plots were visually confusing:

- Grid lines too strong relative to coil markers.
- Coil markers too small or obscured.
- Points near the plot boundary, making some coils appear cut off.
- Only one B‑field snapshot shown, not initial vs final.

**Resolution:**

- Lightened grid lines and increased marker visibility by drawing both circle outlines and central dots.
- Added symmetric r‑extent (`±r`) and dynamic margins so all coils are well inside the viewing window with equal aspect ratio.
- Plotted **both** initial and final B‑fields separately, using independent field evaluations for `t=0` and `t=1`.
- Added HSV vs magnitude background modes.

This made visual debugging much easier: one can now clearly see how field lines “pile up” between plasma and cage, how shielding develops, and whether flux appears conserved visually.

### 5.4 B‑field kernel correctness (elliptic integrals)

There was concern that mistakes in the elliptic‑integral formulas might break flux conservation and produce unphysical fields.

**Resolution:**

- Implemented `field_loop_rz` directly from standard reference formulas for a single circular loop, taking care with the parameter convention (`m = k^2`) used by `scipy.special.ellipk/ellipe`.
- Introduced a robust flux conservation check using the same inductance machinery that is built from elliptic integrals; agreement between `Phi0` and `Phi_final` provides an indirect but strong consistency test.
- Recommended simple analytic spot checks (not included in the demo):
  - Compare on‑axis `Bz` from `field_loop_rz` at `r=0` to the textbook

    ```text
    Bz_on_axis = mu0 * I * a^2 / (2 * (a^2 + z^2)^(3/2))
    ```

  - Check symmetry: `Br(r=0) = 0`, reflection symmetry in `z` for single
    loop, etc.

### 5.5 Input robustness (`np.genfromtxt` quirks)

While refactoring input handling, a couple of minor but practical issues showed up:

- Using `split_whitespace=True` with `np.genfromtxt` is not portable.
- When only a single coil line is present, `genfromtxt` may return a 1D array instead of 2D, which confused the parser.

**Resolution:**

- Simplified to plain `np.genfromtxt(..., dtype=str)` and treated the result as arbitrary shape; `parse_coils` now explicitly reshapes 1D arrays to `(1, -1)` and iterates safely.


## 6. Current Capabilities and Limitations

### What the Python prototype does now

- Builds inductance matrices for arbitrary sets of coaxial loops.
- Computes initial fluxes and total magnetic energy.
- Evolves coil geometries from `(R0, z0)` to `(R1, z1)` under strict per‑loop
  flux conservation, solving `K(t) @ I(t) = Phi0` at each step.
- Produces:
  - Time series of total magnetic energy and per‑coil currents.
  - Visualizations of initial and final coil geometry with magnetic field
    background and streamlines.
- Prints a per‑coil flux conservation table to verify numerical correctness.
- Generates simple parabolic nozzle configurations for quick experimentation
  without external input files.

### What it intentionally does *not* do yet

- No multi‑ring or full soft‑body plasma dynamics (no springs between many
  plasma segments, no remeshing, no shape change beyond a single loop).
- No resistive losses, radiation, or finite conductivity effects.
- No full MHD (pressure, density, temperature fields in the volume).
- No control‑point (Method 2) solver; everything is Method 1 (flux/circuit).

The main scripts described so far (especially `demo_coil_motion_flux.py`)
are therefore best viewed as **static flux‑compression circuit solvers and
visualizers**, not as full MHD engines.


## 7. Relation to the WebGL Design Document

The WebGL design document describes a more ambitious interactive demo with both **flux‑conservation** and **control‑point** solvers, UI, and 3D rendering. The Python code here corresponds to a subset of that design:

- It implements *Method 1: Flux Conservation* in a clean, inspectable way using NumPy and SciPy.
- It shares the same underlying B‑field kernel structure (elliptic integrals for circular loops).
- It focuses on numerical and physical sanity checks (flux conservation, energy trends, qualitative field shapes) that are harder to debug once the solver is buried inside GLSL shaders and interactive UI code.

Once confidence is gained in the Python prototype, the same inductance and flux‑conservation logic can be ported to JS/WebGL (or C++/OpenCL) following this module as a reference.


## 8. How to Use This Prototype in Practice

For day‑to‑day experimentation you generally do **not** import this module from Python. Instead you:

1. Change directory to `doc/python/MHD`.
2. Run `demo_coil_motion_flux.py` with either
   - a custom coils file, or
   - no file to use the parabolic nozzle generator.
3. Inspect
   - the printed flux conservation table
   - energy/current time series
   - geometry + field plots.

Example commands:

- Default parabolic nozzle with magnitude‑only B‑background:

  ```bash
  python demo_coil_motion_flux.py
  ```

- Same, but with HSV background to see direction of `B`:

  ```bash
  python demo_coil_motion_flux.py --bg-mode hsv
  ```

- Custom geometry from file and no geometry plot (only time series):

  ```bash
  python demo_coil_motion_flux.py coils.dat --no-geometry
  ```

- More interpolation steps for smoother curves:

  ```bash
  python demo_coil_motion_flux.py --steps 200
  ```


## 9. Dynamic Single‑Plasma‑Coil Demo

To start exploring dynamics while still keeping the system simple, there is an additional script:

- `demo_plasma_dynamics_flux.py`

This script reuses the inductance matrix machinery and flux conservation but replaces prescribed geometry with **explicit motion of a single plasma coil**:

- The seed SC coil (and optional cage coil) are fixed in space.
- The plasma is represented by one circular loop with:
  - mass and velocity in `(r, z)`
  - current determined at each time step from the flux‑conserving system
    `K(t) @ I(t) = Phi0` (all coils included).
- At each time step the code:
  1. Builds `K(rs, zs)` for the current geometry and solves for `I(t)`.
  2. Computes the Lorentz force on the plasma loop from the fields of the
     other coils (using the same Biot–Savart kernel).
  3. Adds two **expansion forces** on the plasma loop:
     - a gas pressure term with a simple cylindrical volume model
       `V ~ π r² L_eff` and adiabatic law `P ~ (V0/V)^γ`,
     - a magnetic self‑pressure (hoop stress) term proportional to `I_p²/r`.
  4. Integrates position and velocity forward in time with explicit Euler.

Diagnostics and plots include:

- Plasma trajectory `r(t)` and `z(t)`,
- kinetic, magnetic, and total energy vs time,
- coil currents vs time,
- a phase‑space style trajectory in the `(z, r)` plane.

This dynamic demo is still a **toy model**:

- Only a single plasma ring is mobile (no deformable surface),
- integration is explicit Euler, so time step must be chosen small enough
  for stability,
- no coupling to a full gas volume or to multiple plasma segments.

However, it forms a useful bridge between the purely static flux demos and the more ambitious soft‑body / MHD dynamics envisioned in the design documents. 
It lets you test how flux‑conserving currents, external fields, and simple pressure models interact to accelerate or decelerate an expanding plasma ring before moving to a full WebGL or C++ implementation.


## 10. Spherical Plasma Shell Demo

To better approximate a **3D plasma “ball”** and study shielding/leakage through a cage, an additional demo script was added:

- `demo_plasma_sphere_flux.py`

This script reuses the same inductance and flux‑conservation machinery but changes how the geometry is generated.

### 10.1 Geometry: parabolic cage + spherical plasma

The configuration is built in two steps, using helpers from `inductance_core.py`:

1. **Parabolic SC + cage** via `generate_parabolic_nozzle(...)`:

   - Generates an SC seed coil plus a set of `CAGE` rings following a parabolic wall defined by throat radius `r_throat`, exit radius `r_exit` and axial range `[z_start, z_end]`.
   - Output is a text block in the usual

     ```text
     TYPE  R0  Z0  R1  Z1  I0
     ```

     format. `demo_plasma_sphere_flux.py` discards the parabolic `PLASMA` line and keeps only the `SC` and `CAGE` entries, so the cage geometry is identical to the nozzle demo.

2. **Spherical plasma shell** via `generate_spherical_plasma_loops(...)`:

   - Implemented directly in `inductance_core.py`.
   - Produces multiple `PLASMA` loops approximating a spherical shell in the `(r, z)` meridional plane:

     ```text
     PLASMA  R0  Z0   R1  Z1   I0
     ```

   - The shell is centered at `z0` and expands **purely radially** from radius `r0` to `r1` (the center does not move).
   - A small number of loops (e.g. 5–8) is used so the geometry remains readable while still approximating a sphere.

The demo concatenates the parabolic `SC` + `CAGE` lines with the spherical `PLASMA` lines and feeds the result into `run_coil_motion_flux`, exactly like any other coils file.

### 10.2 CLI and typical use

`demo_plasma_sphere_flux.py` exposes CLI parameters for both the plasma ball and the parabolic wall, in a compact one‑line style:

- Plasma shell:
  - `--n-plasma` – number of plasma loops on the sphere
  - `--r0`, `--r1` – initial and final sphere radii
  - `--z0` – sphere center `z` (fixed for both initial and final state)

- Seed coil:
  - `--sc-r`, `--sc-z`, `--sc-current`

- Parabolic cage:
  - `--n-cage` – number of cage rings from throat to exit
  - `--r-throat`, `--r-exit` – throat and exit radii
  - `--z-start`, `--z-end` – axial extent of the parabola

- Solver / plotting:
  - `--steps` – number of interpolation steps
  - `--no-timeseries`, `--no-geometry` – disable plots
  - `--bg-mode` – HSV vs magnitude background for B‑field plots

The seed coil is typically centered at the parabola vertex (`z = 0`), and the plasma sphere expands from a tiny initial radius (e.g. `r0 = 0.01`) to a radius that just touches the throat (e.g. `r1 = r_throat`).

### 10.3 Insights from the spherical demo

- A **finite number of cage rings** and plasma loops means shielding is never perfectly “topological”: field lines can and do leak between rings, which is physically expected.
- When the SC, parabolic cage and spherical plasma are all built from the same inductance machinery and flux‑conserving solve, any apparent “field penetrating the ball” is due to geometry and discretization, not a kernel bug.
- The shared plotting utilities (`MHD_plots.py`) make it easier to visually compare different geometries (nozzle vs sphere) with consistent B‑field rendering.


## 11. B‑Field Kernel and Inductance Consistency Checks

Because the demos use elliptic‑integral formulas in **two different roles**:

- to build the inductance matrix (`calc_mutual_inductance`, `calc_self_inductance`), and
- to visualize fields and compute forces (`field_loop_rz`),

it is important to verify that these kernels are mutually consistent.

Two small standalone check scripts were added for this purpose:

- `check_B_kernels.py`
- `check_B_consistency.py`

### 11.1 On‑axis field check

`check_B_kernels.py` compares the on‑axis `Bz` from `field_loop_rz` against the textbook formula for a single loop of radius `a` and current `I`:

```text
Bz_on_axis(z) = mu0 * I * a^2 / (2 * (a^2 + (z - z0)^2)^(3/2))
```

This is a clean, non‑ambiguous test because:

- the geometry is simple (single loop, axis evaluation), and
- the analytic expression is exact.

Result: numerical and analytic `Bz` agree to machine precision (relative errors ~1e‑16), strongly confirming the correctness and normalization of `field_loop_rz`.

### 11.2 Mutual inductance vs integrated B: “apples to apples”

Attempting to compare a crude **self‑flux integral** directly to `calc_self_inductance` produced huge discrepancies, because

- `calc_self_inductance` assumes a finite wire radius and an internal field model, while
+- integrating `Bz` from a **filamentary** kernel over the loop area is not the same quantity and is formally divergent near the wire.

To avoid this mismatch, `check_B_consistency.py` focuses on **mutual inductance** between *distinct* loops, where the geometry is non‑singular and the physics is unambiguous:

1. Compute

   ```python
   M = calc_mutual_inductance(R_src, Z_src, R_tgt, Z_tgt)
   Phi_analytic = M * I
   ```

2. Numerically integrate the flux of `field_loop_rz` from the source loop through the target loop’s area:

   ```text
   Phi_numeric = ∫_0^{R_tgt} Bz(r, Z_tgt; R_src, Z_src) * 2π r dr
   ```

3. Compare `Phi_numeric` vs `Phi_analytic` for:
   - separated loops (e.g. cage vs plasma at different `z`),
   - concentric loops (e.g. plasma expanding inside a larger seed coil).

Result: relative errors are in the `1e-8 – 1e-9` range, i.e. essentially perfect agreement given floating‑point and quadrature. This shows that:

- `calc_mutual_inductance` and `field_loop_rz` are mathematically consistent,
- the B‑field used for visualization and forces is exactly the same field that underlies the inductance matrix used in the flux solver.

### 11.3 Takeaways for interpreting the demos

- The **B‑field kernel and inductance matrix are internally consistent**; discrepancies between intuition and plots are almost certainly due to geometry and discretization (finite number of cage rings / plasma loops), not kernel errors.
- Self‑inductance is inherently sensitive to **core radius and internal current distribution**. Naively integrating filamentary `Bz` over the loop area is not equivalent to `calc_self_inductance` and can give order‑of‑magnitude differences; mutual‑inductance tests avoid this ambiguity.
- For debugging and design work it is therefore better to:
  - trust the mutual‑inductance vs integrated‑B checks and the on‑axis tests for kernel correctness, and
  - use energy/flux diagnostics plus controlled geometries (simple 2‑coil/3‑coil setups, spherical plasma shell inside a parabolic cage) to reason about physical shielding and flux compression.
