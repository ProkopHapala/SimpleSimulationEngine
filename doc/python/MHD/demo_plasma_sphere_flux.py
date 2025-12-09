"""Flux-conserving spherical plasma shell demo.

This script constructs a configuration consisting of

- one seed superconducting (SC) coil, fixed in space,
- several cage coils between seed and plasma (stacked along z, similar in
  spirit to the cage in the dynamic demo), and
- several plasma coils sampling the surface of a sphere in the (r,z)
  meridional plane (axisymmetric plasma "ball").

The sphere center position and radius are linearly interpolated between an
initial state (t=0) and a final state (t=1). At each interpolation step we
solve for currents that conserve the initial total flux in every loop,
reusing the generic machinery from `demo_coil_motion_flux`.

The geometry is encoded in the same textual format as `demo_coil_motion_flux`:

    TYPE  R0  Z0  R1  Z1  I0

where we use

- TYPE = "SC"      for the seed coil,
- TYPE = "PLASMA"  for the surface coils on the plasma sphere.

Usage (from doc/python/MHD):

    python demo_plasma_sphere_flux.py

Key CLI options:

    --n-plasma     Number of coils on the plasma surface (default: 5)
    --n-cage       Number of cage coils (default: 6)
    --r0           Initial plasma sphere radius
    --r1           Final plasma sphere radius
    --z0           Initial plasma sphere center z
    --z1           Final plasma sphere center z
    --sc-r         Seed coil radius
    --sc-z         Seed coil axial position
    --sc-current   Seed coil current (A)
    --steps        Number of interpolation steps
    --bg-mode      Background mode for field plot ("hsv" or "mag")

The physics and plotting are identical in spirit to demo_coil_motion_flux,
except that we now have many plasma loops approximating a continuous
conductive shell.
"""

import argparse
import numpy as np
from io import StringIO
import matplotlib.pyplot as plt

from inductance_core import (
    generate_sc_seed_coil,
    generate_parabolic_rings,
    generate_spherical_plasma_loops,
)
from MHD_plots import plot_coil_geometry, plot_B_profiles
from demo_coil_motion_flux import run_coil_motion_flux

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=("Flux-conserving evolution of a spherical plasma shell represented by many coaxial loops on the sphere surface."))
    parser.add_argument("--n-plasma", type=int, default=16, help="Number of plasma coils on sphere surface (default: 8)")
    parser.add_argument("--n-cage",   type=int, default=10, help="Number of cage coils (default: 8)")
    parser.add_argument("--r0",       type=float, default=0.01, help="Initial plasma sphere radius")
    parser.add_argument("--r1",       type=float, default=0.9,  help="Final plasma sphere radius")
    parser.add_argument("--z0",       type=float, default=1.0,  help="Plasma sphere center z")
    parser.add_argument("--sc-r",     type=float, default=1.0,  help="Seed coil radius")
    parser.add_argument("--sc-z",     type=float, default=0.0,  help="Seed coil axial position")
    parser.add_argument("--sc-current", type=float, default=1.0e6, help="Seed coil current in A")
    parser.add_argument("--r-throat", type=float, default=0.05,  help="Parabola throat radius")
    parser.add_argument("--r-exit",   type=float, default=1.4,  help="Parabola exit radius")
    parser.add_argument("--z-start",  type=float, default=0.0,  help="Parabola start z")
    parser.add_argument("--z-end",    type=float, default=1.2,  help="Parabola end z")
    parser.add_argument("--steps",    type=int, default=50,     help="Number of interpolation steps")
    parser.add_argument("--no-timeseries", action="store_true", help="Disable energy/current vs interpolation plots")
    parser.add_argument("--no-geometry", action="store_true",   help="Disable geometry plot (initial vs final coils)")
    parser.add_argument("--no-Bprofiles", action="store_true",  help="Disable 1D B-field profiles (axial and radial)")
    parser.add_argument("--bg-mode", choices=["hsv", "mag"], default="mag", help="Background: hsv (direction hue, magnitude value) or mag (brightness only)")

    args = parser.parse_args()

    sc_line    = generate_sc_seed_coil(args.sc_r, args.sc_z, sc_current=args.sc_current, coil_type="SC")
    cage_block = generate_parabolic_rings(n_rings=args.n_cage, r_throat=args.r_throat, r_exit=args.r_exit, z_start=args.z_start, z_end=args.z_end, coil_type="CAGE")
    sphere_str = generate_spherical_plasma_loops(n_plasma=args.n_plasma, r0=args.r0, r1=args.r1, z0=args.z0)
    data_str   = "\n".join([sc_line, cage_block, sphere_str])

    # Reuse the existing demo machinery by feeding it a string buffer.
    from numpy import genfromtxt

    data = genfromtxt(StringIO(data_str), dtype=str)

    # Run flux-conserving solve without geometry; we reconstruct geometry here to overlay profiles.
    res = run_coil_motion_flux(data, n_steps=args.steps, plot_timeseries=not args.no_timeseries, plot_geom=False, bg_mode=args.bg_mode)

    types = res["types"]; R0 = res["R0"]; z0 = res["z0"]; I0 = res["I0"]; R1 = res["R1"]; z1 = res["z1"]; I_final = res["I_final"]

    fig = None; axes = None
    if not args.no_geometry:
        fig, axes = plot_coil_geometry(types, R0, z0, I0, R1, z1, I_final, title="Coil geometry (initial vs final)", bg_mode=args.bg_mode)

    # Overlay 1D B-field profiles and create separate profile figure if requested.
    if not args.no_Bprofiles:
        if axes is not None:
            # Overlay on initial panel (no separate figure) and final panel (with separate figure)
            plot_B_profiles(R0, z0, I0, float(args.z0), axes=axes[0], show_separate=False)
            plot_B_profiles(R1, z1, I_final, float(args.z0), axes=axes[1], show_separate=True)
        else:
            # No geometry: just overlay into a fresh axes using final state
            plot_B_profiles(R1, z1, I_final, float(args.z0), axes=None, show_separate=True)

    plt.tight_layout(); plt.show()