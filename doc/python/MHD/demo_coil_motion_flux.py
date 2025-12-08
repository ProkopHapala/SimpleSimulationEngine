"""Flux-conserving coil motion demo.

This script reads a set of coaxial coils, linearly moves each coil from
(R0, z0) to (R1, z1), and at each interpolation step solves for currents
that conserve the initial total flux in every loop.

Usage
-----

1) Using a custom definition file:

    python demo_coil_motion_flux.py coils.dat

   where each non-empty line of `coils.dat` has the form:

    TYPE   R0   Z0   R1   Z1   I0

   Examples of TYPE: SC, CAGE, PLASMA (labels are used only for plotting).

2) Using an auto-generated parabolic nozzle (no file argument):

    python demo_coil_motion_flux.py

   In this case, a default parabolic nozzle configuration is generated
   via `generate_parabolic_nozzle()` from `inductance_core`.

Options
-------

    --no-timeseries   Disable energy/current vs interpolation plots
    --no-geometry     Disable geometry plot (initial vs final coils)
    --steps N         Number of interpolation steps (default: 50)
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb
from io import StringIO

from inductance_core import (
    build_inductance_matrix,
    compute_flux,
    magnetic_energy,
    generate_parabolic_nozzle,
    field_loop_rz,
)


def parse_coils(data):
    """Parse coil definitions from a token array.

    Expected columns per row:
        type  R0  z0  R1  z1  I0

    `data` is assumed to be an ndarray of shape (N, >=6) with string entries,
    e.g. produced by np.genfromtxt(..., dtype=str).
    """
    data = np.asarray(data)
    if data.ndim == 1:
        data = data[None, :]

    types = []
    R0_list = []
    z0_list = []
    R1_list = []
    z1_list = []
    I0_list = []

    for row in data:
        if row.size < 6:
            continue
        t, R0, z0, R1, z1, I0 = row[:6]
        types.append(str(t))
        R0_list.append(float(R0))
        z0_list.append(float(z0))
        R1_list.append(float(R1))
        z1_list.append(float(z1))
        I0_list.append(float(I0))
    return (
        np.array(types, dtype=object),
        np.array(R0_list, dtype=float),
        np.array(z0_list, dtype=float),
        np.array(R1_list, dtype=float),
        np.array(z1_list, dtype=float),
        np.array(I0_list, dtype=float),
    )


def plot_geometry(types, R0, z0, I0, R1, z1, I_final, title="Coil geometry (initial vs final)", bg_mode="hsv"):
    """Plot initial and final coil positions with current labels in (z,r) plane.

    types : array of strings (e.g. SC, CAGE, PLASMA)
    R0,z0 : initial radii/positions
    R1,z1 : final radii/positions
    I0    : initial currents
    I_final : final currents from solver (same order)
    """
    fig, axes = plt.subplots(1, 2, figsize=(10,10), sharey=True)
    states = [
        (R0, z0, I0,      "t=0"     ),
        (R1, z1, I_final, "t=final" ),
    ]

    # Determine global bounds for margins
    all_R = np.concatenate([R0, R1])
    all_Z = np.concatenate([z0, z1])
    r_span = max(1.0, np.max(np.abs(all_R)) - np.min(np.abs(all_R)))
    r_pad = 0.1 * r_span
    r_extent = np.max(np.abs(all_R)) + r_pad
    r_min_plot, r_max_plot = -r_extent, r_extent
    z_min, z_max = all_Z.min(), all_Z.max()
    dz = 0.1 * max(1.0, z_max - z_min)

    Nr, Nz = 120, 240
    r_grid = np.linspace(r_min_plot, r_max_plot, Nr)
    z_grid = np.linspace(z_min - dz, z_max + dz, Nz)
    Zg, Rg = np.meshgrid(z_grid, r_grid)

    def make_field_and_rgb(Rs, Zs, Is):
        Br = np.zeros_like(Rg)
        Bz = np.zeros_like(Rg)
        for a, zc, I in zip(Rs, Zs, Is):
            Br_i_abs, Bz_i = field_loop_rz(a, zc, I, np.abs(Rg), Zg)
            Br += np.sign(Rg) * Br_i_abs
            Bz += Bz_i
        Bmag = np.sqrt(Br * Br + Bz * Bz)
        B_ref = 0.01 * np.max(Bmag)  # damp dominant field to reveal weaker structures
        Bmag_norm = np.clip(Bmag / max(B_ref, 1e-30), 0.0, 1.0)
        phi = np.arctan2(Bz, Br)  # direction in (r,z)
        if bg_mode == "hsv":
            hue = (phi + np.pi) / (2 * np.pi)
            sat = np.ones_like(hue)
            val = Bmag_norm
            hsv = np.stack([hue, sat, val], axis=-1)
            rgb = hsv_to_rgb(hsv)
        else:
            cmap = plt.cm.magma
            rgb = cmap(Bmag_norm)[..., :3]
        return Br, Bz, rgb

    Br_init, Bz_init, rgb_initial = make_field_and_rgb(R0, z0, I0)
    Br_final, Bz_final, rgb_final = make_field_and_rgb(R1, z1, I_final)
    field_states = [(Br_init, Bz_init, rgb_initial), (Br_final, Bz_final, rgb_final)]

    for ax, (Rs, Zs, Is, subt), (Br_s, Bz_s, rgb) in zip(axes, states, field_states):
        ax.imshow( rgb, extent=(z_min - dz, z_max + dz, r_min_plot, r_max_plot), origin="lower", aspect="equal", alpha=0.8, )
        # Field lines (streamplot) in the (z,r) plane; note streamplot expects x=z_grid, y=r_grid
        ax.streamplot( z_grid, r_grid, Bz_s, Br_s, color="k", linewidth=0.6, density=1.2, arrowsize=0.7, )
        for typ, r, z, I in zip(types, Rs, Zs, Is):
            I_MA = I / 1e6
            if typ ==   "SC":
                color = "k"
            elif typ == "CAGE":
                color = "g"
            else:
                color = "r"
            # Circle outline
            radius_draw = 0.02 * max(1.0, r_extent)
            for sign in (+1, -1):
                r_pos = sign * r
                circle = plt.Circle((z, r_pos), radius_draw, fill=False, color=color, linewidth=1.0)
                ax.add_patch(circle)
                # Marker at coil center for visibility
                ax.plot(z, r_pos, marker="o", markersize=4, color=color)
            #ax.text(z, r, f"{typ}\n{I_MA:.6f} MA", ha="center", va="bottom", fontsize=8)
            ax.text(z, r, f"{I_MA:.6f} MA", ha="center", va="bottom", fontsize=8)

        ax.set_title(subt)
        ax.set_xlabel("z [m]")
        ax.set_xlim(z_min - dz, z_max + dz)
        ax.set_ylim(r_min_plot, r_max_plot)
        ax.set_aspect("equal", adjustable="box")
        ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.5)

    axes[0].set_ylabel("r [m]")
    fig.suptitle(title)
    plt.tight_layout()


def run_coil_motion_flux(data, n_steps=50, plot_timeseries=True, plot_geom=True, bg_mode="hsv"):
    """Simulate flux-conserving evolution of a set of coils moving in (R,z).

    Scenario:
    - At t = 0, each coil i is at (R0[i], z0[i]) and carries current I0[i].
    - We compute the initial inductance matrix K0 and flux vector Phi0 = K0 @ I0.
    - During a fast event (explosion / compression), each coil moves linearly
      from (R0,z0) to (R1,z1). We assume all coils are ideal (superconducting
      or perfectly conducting) so that Phi_i is conserved for each closed loop.
    - At each intermediate geometry, we solve K(t) @ I(t) = Phi0.

    This models typical flux-conserving machines (explosive generators,
    flux-compression, plasma expansion against seed coils, etc.) where
    there is no external power supply during the fast motion.
    """
    types, R0, z0, R1, z1, I0 = parse_coils(data)
    N = types.size
    if N == 0:
        raise ValueError("No coils loaded from dat file")
    # Initial geometry and flux
    rs0 = R0.copy()
    zs0 = z0.copy()
    K0 = build_inductance_matrix(rs0, zs0)
    Phi0 = compute_flux(K0, I0)
    ts = np.linspace(0.0, 1.0, n_steps)
    energies = []
    currents = []  # shape (n_steps, N)

    for t in ts:
        rs = (1.0 - t) * R0 + t * R1
        zs = (1.0 - t) * z0 + t * z1
        Kt = build_inductance_matrix(rs, zs)
        It = np.linalg.solve(Kt, Phi0)
        Wt = magnetic_energy(Kt, It)
        energies.append(Wt)
        currents.append(It)

    energies = np.array(energies)
    currents = np.array(currents)  # (n_steps, N)

    # Flux conservation check at final step
    K_final = build_inductance_matrix(R1, z1)
    Phi_final = compute_flux(K_final, currents[-1])
    print("\nFlux check (initial vs final):")
    print(" idx  type        Phi0 [Wb]           Phi_final [Wb]       rel_err")
    for i, typ in enumerate(types):
        phi0_i = Phi0[i]
        phif_i = Phi_final[i]
        rel = 0.0 if phi0_i == 0 else (phif_i - phi0_i) / phi0_i
        print(f" {i:2d}   {typ:>6}   {phi0_i: .6e}   {phif_i: .6e}   {rel: .3e}")

    if plot_timeseries:
        # Simple diagnostics: total energy and per-coil currents vs step index
        fig, ax = plt.subplots(2, 1, figsize=(7, 8), sharex=True)
        ax[0].plot(ts, energies)
        ax[0].set_ylabel("Magnetic energy [J]")
        ax[0].set_title("Total magnetic energy vs interpolation step")
        ax[0].grid(True)
        for i in range(N):  ax[1].plot(ts, currents[:, i] / 1e6, label=f"I_{i} ({types[i]}) [MA]")
        ax[1].set_xlabel("Interpolation parameter t (0 = initial, 1 = final)")
        ax[1].set_ylabel("Current [MA]")
        ax[1].legend()
        ax[1].grid(True)
        plt.tight_layout()
    if plot_geom:
        # Geometry figure: initial vs final positions and currents
        plot_geometry(types, R0, z0,   I0, R1, z1, currents[-1], title="Coil geometry (initial vs final)", bg_mode=bg_mode)
        plt.show()

    


'''
Examples
========

# 1) Default parabolic nozzle, 50 steps, with all plots
python demo_coil_motion_flux.py

# 2) Custom throat/exit radii and length
python demo_coil_motion_flux.py --r-throat 1.0 --r-exit 2.5 --z-start 0.0 --z-end 3.0

# 3) Shifted seed coil and stronger current
python demo_coil_motion_flux.py --sc-r 1.8 --sc-z -0.5 --sc-current 2.0e6

# 4) Plasma starting downstream and ending before exit
python demo_coil_motion_flux.py --plasma-r0 0.3 --plasma-r1 2.0 --plasma-z0 0.5 --plasma-z1 1.8

# 5) Use external file and disable geometry plot
python demo_coil_motion_flux.py coils.dat --no-geometry

'''

if __name__ == "__main__":
    """Entry point for CLI argument parsing and demo execution."""
    parser = argparse.ArgumentParser( description="Flux-conserving coil motion demo (see module docstring for details)."  )
    parser.add_argument( "coils_file", nargs="?", default=None, help="Path to coils.dat file (if omitted, generate a parabolic nozzle)", )
    parser.add_argument( "--steps", type=int, default=50, help="Number of interpolation steps (default: 50)", )
    # Parabolic nozzle geometry parameters (used when no coils_file is given)
    parser.add_argument( "--n-rings",    type=int,   default=5, help="Number of cage rings in generated nozzle (default: 12)", )
    parser.add_argument( "--r-throat",   type=float, default=0.2, help="Throat radius R_throat (default: 1.0)", )
    parser.add_argument( "--r-exit",     type=float, default=2.0, help="Exit radius R_exit (default: 3.0)", )
    parser.add_argument( "--z-start",    type=float, default=0.5, help="Axial start position z_start (default: 0.0)", )
    parser.add_argument( "--z-end",      type=float, default=1.8, help="Axial end position z_end (default: 2.0)", )

    parser.add_argument( "--sc-r",       type=float, default=0.3, help="Seed coil radius (default: 1.5 * r_throat)", )
    parser.add_argument( "--sc-z",       type=float, default=0.3, help="Seed coil axial position (default: z_start)", )
    parser.add_argument( "--sc-current", type=float, default=1.0e6, help="Seed coil current in A (default: 1e6)", )

    parser.add_argument( "--plasma-r0",  type=float, default=0.01, help="Plasma initial radius at throat (default: 0.2)", )
    parser.add_argument( "--plasma-r1",  type=float, default=1.0, help="Plasma final radius at exit (default: 0.9 * r_exit)", )
    parser.add_argument( "--plasma-z0",  type=float, default=1.0, help="Plasma initial z (default: z_start)", )
    parser.add_argument( "--plasma-z1",  type=float, default=1.0, help="Plasma final z (default: z_end)", )
    parser.add_argument( "--no-timeseries", action="store_true", help="Disable energy/current vs interpolation plots", )
    parser.add_argument( "--no-geometry",   action="store_true", help="Disable geometry plot (initial vs final coils)", )
    #parser.add_argument( "--bg-mode", choices=["hsv", "mag"], default="hsv", help="Background: hsv (direction hue, magnitude value) or mag (brightness only)" )
    parser.add_argument( "--bg-mode", choices=["hsv", "mag"], default="mag", help="Background: hsv (direction hue, magnitude value) or mag (brightness only)" )
    args = parser.parse_args()

    if args.coils_file is None:
        # Auto-generate a parabolic nozzle configuration using CLI geometry
        data_str = generate_parabolic_nozzle(
            n_rings=args.n_rings,
            r_throat=args.r_throat,
            r_exit=args.r_exit,
            z_start=args.z_start,
            z_end=args.z_end,
            sc_current=args.sc_current,
            plasma_r_start=args.plasma_r0,
            plasma_r_end=args.plasma_r1,
            sc_r=args.sc_r,
            sc_z=args.sc_z,
            plasma_z_start=args.plasma_z0,
            plasma_z_end=args.plasma_z1,
        )
        data = np.genfromtxt(StringIO(data_str), dtype=str)
    else:
        data = np.genfromtxt(args.coils_file, dtype=str)

    run_coil_motion_flux( data, n_steps=args.steps, plot_timeseries=not args.no_timeseries, plot_geom=not args.no_geometry, bg_mode=args.bg_mode )

    plt.show()
