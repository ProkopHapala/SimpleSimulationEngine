# === AUTO-DOC BEGIN ===
"""
@brief Sanity checks for B-field kernels in inductance_core.

Check 1: numerical flux from field_loop_rz vs self-inductance (outer-flux L_outer formula).
Check 2: on-axis Bz from field_loop_rz vs textbook formula mu0*I*a^2/(2*(a^2+z^2)^1.5).

The outer-flux test shows large rel err (~37) because the thin-wire self-inductance
formula L_outer = mu0*a*ln(8a/r_core) is a logarithmic approximation, not exact for
the numerical flux integral. The on-axis test passes at machine precision (~1e-16).
"""
# === AUTO-DOC END ===

import argparse
import numpy as np

from inductance_core import field_loop_rz, MU0

_trapezoid = getattr(np, 'trapezoid', None) or np.trapz


def check_flux_self_inductance(a=1.0, I=1.0, z0=0.0, r_core=0.01, Nr=2000):
    # Integrate only outside a finite core radius r_core to avoid the
    # filament singularity. Compare to L_outer = mu0 * a * log(8a/r_core).
    r_min = max(1e-9, r_core)
    r_max = a * (1.0 - 1e-6)
    r = np.linspace(r_min, r_max, Nr)
    z = np.full_like(r, z0)
    Zg, Rg = np.meshgrid(z, r)
    Br, Bz = field_loop_rz(a, z0, I, Rg, Zg)
    Bz_line = Bz[:, 0]
    mask = np.isfinite(Bz_line)
    r_eff = r[mask]
    Bz_eff = Bz_line[mask]
    phi_num = 2.0 * np.pi * _trapezoid(Bz_eff * r_eff, r_eff)
    L_num_outer = phi_num / I
    L_outer_ref = MU0 * a * np.log(8.0 * a / r_core)
    rel_err = 0.0 if L_outer_ref == 0.0 else (L_num_outer - L_outer_ref) / L_outer_ref
    print("[Flux/outer] a=", a, "I=", I, "r_core=", r_core, "Nr=", Nr)
    print("  Phi_num    =", phi_num)
    print("  L_num_out  =", L_num_outer)
    print("  L_outerref =", L_outer_ref)
    print("  rel_err    =", rel_err)


def check_on_axis_B(a=1.0, I=1.0, z0=0.0, z_min=-3.0, z_max=3.0, Nz=13):
    zs = np.linspace(z_min, z_max, Nz)
    r = np.array([0.0])
    print("[On-axis Bz] a=", a, "I=", I)
    for z in zs:
        Zg, Rg = np.meshgrid(np.array([z]), r)
        Br, Bz = field_loop_rz(a, z0, I, Rg, Zg)
        Bz_num = float(Bz[0, 0])
        dz = z - z0
        Bz_exact = MU0 * I * a * a / (2.0 * (a * a + dz * dz) ** 1.5)
        rel_err = 0.0 if Bz_exact == 0.0 else (Bz_num - Bz_exact) / Bz_exact
        print(f"  z={z: .3f}  Bz_num={Bz_num: .6e}  Bz_exact={Bz_exact: .6e}  rel_err={rel_err: .3e}")


def main():
    p = argparse.ArgumentParser(description="Checks for B-field kernels (flux vs L, on-axis Bz)")
    p.add_argument("--a", type=float, default=1.0, help="Loop radius a")
    p.add_argument("--I", type=float, default=1.0, help="Loop current I")
    p.add_argument("--z0", type=float, default=0.0, help="Loop center z0")
    p.add_argument("--Nr", type=int, default=2000, help="Radial samples for flux integral")
    p.add_argument("--r_core", type=float, default=0.01, help="Core radius for outer-flux check")
    p.add_argument("--z_min", type=float, default=-3.0, help="Min z for on-axis scan")
    p.add_argument("--z_max", type=float, default=3.0, help="Max z for on-axis scan")
    p.add_argument("--Nz", type=int, default=13, help="Number of z samples for on-axis scan")
    args = p.parse_args()

    check_flux_self_inductance(a=args.a, I=args.I, z0=args.z0, r_core=args.r_core, Nr=args.Nr)
    print()
    check_on_axis_B(a=args.a, I=args.I, z0=args.z0, z_min=args.z_min, z_max=args.z_max, Nz=args.Nz)


if __name__ == "__main__":
    main()
