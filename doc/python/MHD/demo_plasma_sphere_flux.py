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

from inductance_core import generate_parabolic_nozzle, generate_spherical_plasma_loops
from demo_coil_motion_flux import run_coil_motion_flux


def generate_spherical_plasma_config(n_plasma=5, r0=0.01, r1=0.5, z0=0.5, z1=0.5, sc_r=1.0, sc_z=0.0, sc_current=1.0e6, n_cage=6, r_throat=0.2, r_exit=2.0, z_start=0.0, z_end=1.8):
    """Generate a coil definition string for a spherical plasma shell.

    Parameters
    ----------
    n_plasma : int
        Number of plasma coils sampling the sphere surface.
    r0, r1 : float
        Initial and final plasma sphere radii.
    z0, z1 : float
        Plasma sphere center positions along z. The sphere center is kept
        fixed at z0; z1 is currently ignored for the plasma and only passed
        through to keep the parabolic helper signature consistent.
    sc_r, sc_z : float
        Radius and axial position of the seed SC coil (fixed).
    sc_current : float
        Current in the seed SC coil (A).
    n_cage : int
        Number of cage coils (parabolic rings) along the wall.
    r_throat, r_exit, z_start, z_end : float
        Parabolic wall parameters passed to generate_parabolic_nozzle.

    Returns
    -------
    data_str : str
        Multi-line string with lines of the form

            TYPE  R0  Z0  R1  Z1  I0

        suitable for np.genfromtxt(..., dtype=str).
    """
    lines = []

    # 1) Parabolic SC + cage coils (nozzle-style); drop its PLASMA line.
    base_str = generate_parabolic_nozzle(n_rings=n_cage, r_throat=r_throat, r_exit=r_exit, z_start=z_start, z_end=z_end, sc_current=sc_current, plasma_r_start=r0, plasma_r_end=r1, sc_r=sc_r, sc_z=sc_z, plasma_z_start=z0, plasma_z_end=z0)

    from numpy import genfromtxt

    base_arr = genfromtxt(StringIO(base_str), dtype=str)
    if base_arr.ndim == 1:
        base_arr = base_arr[None, :]

    for row in base_arr:
        typ = row[0]
        if typ == "PLASMA": continue
        t, R0s, Z0s, R1s, Z1s, I0s = row[:6]
        lines.append(f"{t:6s}  {float(R0s):.6f}  {float(Z0s):.6f}   {float(R1s):.6f}  {float(Z1s):.6f}   {float(I0s):.6e}")

    # 2) Spherical plasma shell loops (fixed center z0, pure radial expansion)
    sphere_str = generate_spherical_plasma_loops(n_plasma=n_plasma, r0=r0, r1=r1, z0=z0)
    if len(sphere_str.strip()) > 0:
        for line in sphere_str.splitlines(): lines.append(line)

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=("Flux-conserving evolution of a spherical plasma shell represented by many coaxial loops on the sphere surface."))
    parser.add_argument("--n-plasma", type=int, default=8, help="Number of plasma coils on sphere surface (default: 8)")
    parser.add_argument("--n-cage",   type=int, default=8, help="Number of cage coils (default: 8)")
    parser.add_argument("--r0",       type=float, default=0.01, help="Initial plasma sphere radius")
    parser.add_argument("--r1",       type=float, default=0.4,  help="Final plasma sphere radius")
    parser.add_argument("--z0",       type=float, default=0.5,  help="Plasma sphere center z")
    parser.add_argument("--sc-r",     type=float, default=1.0,  help="Seed coil radius")
    parser.add_argument("--sc-z",     type=float, default=0.0,  help="Seed coil axial position")
    parser.add_argument("--sc-current", type=float, default=1.0e6, help="Seed coil current in A")
    parser.add_argument("--r-throat", type=float, default=0.2,  help="Parabola throat radius")
    parser.add_argument("--r-exit",   type=float, default=2.0,  help="Parabola exit radius")
    parser.add_argument("--z-start",  type=float, default=0.0,  help="Parabola start z")
    parser.add_argument("--z-end",    type=float, default=1.8,  help="Parabola end z")
    parser.add_argument("--steps",    type=int, default=50,     help="Number of interpolation steps")
    parser.add_argument("--no-timeseries", action="store_true", help="Disable energy/current vs interpolation plots")
    parser.add_argument("--no-geometry", action="store_true",   help="Disable geometry plot (initial vs final coils)")
    parser.add_argument("--bg-mode", choices=["hsv", "mag"], default="mag", help="Background: hsv (direction hue, magnitude value) or mag (brightness only)")

    args = parser.parse_args()

    data_str = generate_spherical_plasma_config(n_plasma=args.n_plasma, r0=args.r0, r1=args.r1, z0=args.z0, z1=args.z0, sc_r=args.sc_r, sc_z=args.sc_z, sc_current=args.sc_current, n_cage=args.n_cage, r_throat=args.r_throat, r_exit=args.r_exit, z_start=args.z_start, z_end=args.z_end)

    # Reuse the existing demo machinery by feeding it a string buffer.
    from numpy import genfromtxt

    data = genfromtxt(StringIO(data_str), dtype=str)

    run_coil_motion_flux(data, n_steps=args.steps, plot_timeseries=not args.no_timeseries, plot_geom=not args.no_geometry, bg_mode=args.bg_mode)


if __name__ == "__main__":
    main()
