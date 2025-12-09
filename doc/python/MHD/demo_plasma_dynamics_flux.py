"""Dynamic plasma-coil demo with flux-conserving currents.

This script is a dynamic counterpart to ``demo_coil_motion_flux.py``.
Instead of prescribing a linear motion of all coils between (R0,Z0) and
(R1,Z1), it treats a single plasma coil as a moving mass point with
inertia in (r,z). Seed and (optionally) cage coils remain fixed.

At each time step:

- Geometry of all loops (seed, optional cage, plasma) is known.
- We build the inductance matrix K(r,z) for the current geometry.
- We enforce per-loop flux conservation using the initial flux vector
  Phi0 = K0 @ I0, by solving

      K(t) @ I(t) = Phi0

- Using the resulting currents, we compute the Lorentz force on the
  plasma coil from the magnetic field of the other coils.
- We integrate plasma position and velocity forward in time.

Outputs:

- Plasma trajectory r(t), z(t)
- Kinetic, magnetic, and total energy vs time
- Coil currents vs time
- Optional check of flux conservation over time

This is still a toy model:

- Single moving loop for the plasma surface
- No gas pressure, springs, or volume evolution
- No resistivity or losses

But it combines the flux-conserving circuit picture with a minimal
explicit dynamics of the plasma coil.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt

from inductance_core import (
    build_inductance_matrix,
    compute_flux,
    magnetic_energy,
    field_loop_rz,
)


def setup_system(with_cage: bool = True):
    """Set up a simple Seed(+Cage)+Plasma configuration.

    Returns
    -------
    types : ndarray of str
        Coil labels ("SC", "CAGE", "PLASMA").
    rs, zs : ndarray
        Radii and axial positions.
    I0 : ndarray
        Initial currents.
    movable_idx : int
        Index of the movable plasma coil in the arrays.
    mass : float
        Mass associated with the plasma coil.
    v0 : ndarray, shape (2,)
        Initial velocity components (vr, vz) of the plasma coil.
    """
    types = []
    rs = []
    zs = []
    I0 = []

    # Seed SC coil: fixed position, fixed initial current
    types.append("SC")
    rs.append(1.0)
    zs.append(0.0)
    I0.append(1.0e6)

    # Optional cage coil: fixed, initially zero current
    if with_cage:
        types.append("CAGE")
        rs.append(0.8)
        zs.append(0.5)
        I0.append(0.0)

    # Plasma coil: small initial radius, zero current, movable
    movable_idx = len(types)
    types.append("PLASMA")
    rs.append(0.2)
    zs.append(0.5)
    I0.append(0.0)

    types = np.array(types, dtype=object)
    rs = np.array(rs, dtype=float)
    zs = np.array(zs, dtype=float)
    I0 = np.array(I0, dtype=float)

    mass = 0.1  # arbitrary mass for visualization [kg]
    v0 = np.array([1000.0, 0.0])  # initial (vr, vz) [m/s]
    return types, rs, zs, I0, movable_idx, mass, v0


def compute_force_on_plasma(idx_plasma, types, rs, zs, I):
    """Compute Lorentz force on plasma coil from other coils.

    We use the standard thin-loop expression:

        F_r = I_plasma * (2*pi*r_plasma) * Bz_total
        F_z = I_plasma * (2*pi*r_plasma) * (-Br_total)

    where (Br_total, Bz_total) is the magnetic field at the plasma loop
    center due to all *other* coils.
    """
    r_p = rs[idx_plasma]
    z_p = zs[idx_plasma]
    I_p = I[idx_plasma]

    if I_p == 0.0 or r_p <= 0.0:
        return 0.0, 0.0

    Br_tot = 0.0
    Bz_tot = 0.0
    for j in range(len(types)):
        if j == idx_plasma:
            continue
        a = rs[j]
        zc = zs[j]
        Ij = I[j]
        if Ij == 0.0 or a <= 0.0:
            continue
        # field_loop_rz works on arrays; call it on scalars via 0-d arrays
        Br_j, Bz_j = field_loop_rz(a, zc, Ij, np.array([[r_p]]), np.array([[z_p]]))
        Br_tot += float(Br_j[0, 0])
        Bz_tot += float(Bz_j[0, 0])

    L_len = 2.0 * np.pi * r_p
    Fr = I_p * L_len * Bz_tot
    Fz = I_p * L_len * (-Br_tot)
    return Fr, Fz


def run_plasma_dynamics(
    with_cage: bool = True,
    dt: float = 5e-6,
    n_steps: int = 400,
    r_min: float = 0.01,
    P0: float = 1e6,
    gamma_gas: float = 5.0 / 3.0,
    L_eff: float = 1.0,
    mag_self_coeff: float = 1e-7,
):
    """Run dynamic simulation of a single plasma coil.

    Parameters
    ----------
    with_cage : bool
        If True, include a fixed cage coil between SC and plasma.
    dt : float
        Time step [s].
    n_steps : int
        Number of integration steps.
    r_min : float
        Minimum allowed radius for the plasma coil (simple axis guard).
    """
    types, rs, zs, I0, idx_p, mass, v0 = setup_system(with_cage=with_cage)

    # Initial inductance matrix and flux (Phi0)
    # Plasma may carry initial magnetization from the seed field; we keep
    # whatever external flux is present at t=0 and conserve it, while
    # additional gas and magnetic self-pressure terms act to drive
    # expansion and prevent collapse.
    K0 = build_inductance_matrix(rs, zs)
    Phi0 = compute_flux(K0, I0)

    # State variables for plasma coil
    vr, vz = float(v0[0]), float(v0[1])

    # Gas pressure model: simple cylindrical volume V ~ pi * r^2 * L_eff
    r0 = float(rs[idx_p])
    V0 = np.pi * r0 * r0 * L_eff

    # Histories
    times = []
    r_hist = []
    z_hist = []
    E_kin_hist = []
    E_mag_hist = []
    E_tot_hist = []
    I_hist = []  # shape (steps, N)

    t = 0.0
    I = I0.copy()

    for step in range(n_steps):
        # 1) Solve flux-conserving currents at current geometry
        K = build_inductance_matrix(rs, zs)
        I = np.linalg.solve(K, Phi0)

        # 2) Compute energies
        E_mag = magnetic_energy(K, I)
        E_kin = 0.5 * mass * (vr * vr + vz * vz)

        # 3) Log state
        times.append(t)
        r_hist.append(rs[idx_p])
        z_hist.append(zs[idx_p])
        E_mag_hist.append(E_mag)
        E_kin_hist.append(E_kin)
        E_tot_hist.append(E_mag + E_kin)
        I_hist.append(I.copy())

        # 4) Compute forces on plasma coil
        # 4a) Lorentz force from external fields
        Fr_ext, Fz_ext = compute_force_on_plasma(idx_p, types, rs, zs, I)

        # 4b) Gas pressure force (outward), adiabatic P ~ (V0/V)^gamma
        r_p = float(rs[idx_p])
        V = np.pi * max(r_p, r_min) ** 2 * L_eff
        P = P0 * (V0 / V) ** gamma_gas
        Fr_gas = P * (2.0 * np.pi * r_p * L_eff)

        # 4c) Magnetic self-pressure (hoop stress) ~ I_p^2 / r
        I_p = float(I[idx_p])
        Fr_self = mag_self_coeff * I_p * I_p / max(r_p, r_min)

        Fr = Fr_ext + Fr_gas + Fr_self
        Fz = Fz_ext

        # 5) Integrate one time step (simple explicit Euler)
        ar = Fr / mass
        az = Fz / mass
        vr += ar * dt
        vz += az * dt
        rs[idx_p] += vr * dt
        zs[idx_p] += vz * dt

        # Prevent radius from collapsing through axis
        if rs[idx_p] < r_min:
            rs[idx_p] = r_min
            vr = -0.5 * vr  # crude bounce

        t += dt

    times      = np.array(times)
    r_hist     = np.array(r_hist)
    z_hist     = np.array(z_hist)
    E_mag_hist = np.array(E_mag_hist)
    E_kin_hist = np.array(E_kin_hist)
    E_tot_hist = np.array(E_tot_hist)
    I_hist     = np.array(I_hist)

    return {
        "types": types,
        "times": times,
        "r": r_hist,
        "z": z_hist,
        "E_mag": E_mag_hist,
        "E_kin": E_kin_hist,
        "E_tot": E_tot_hist,
        "I_hist": I_hist,
    }


def plot_results(result):
    """Plot trajectory, energies, and currents from the simulation result."""
    types = result["types"]
    t_ms = result["times"] * 1e3
    r = result["r"]
    z = result["z"]
    E_mag = result["E_mag"]
    E_kin = result["E_kin"]
    E_tot = result["E_tot"]
    I_hist = result["I_hist"]  # (steps, N)

    N = I_hist.shape[1]

    fig, axes = plt.subplots(4, 1, figsize=(8, 12), sharex=True)

    # 1) Plasma trajectory
    ax = axes[0]
    ax.plot(t_ms, r, label="r (radius)")
    ax.plot(t_ms, z, label="z (axial)")
    ax.set_ylabel("Position [m]")
    ax.set_title("Plasma coil trajectory")
    ax.legend()
    ax.grid(True)

    # 2) Energies
    ax = axes[1]
    ax.plot(t_ms, E_mag, label="Magnetic")
    ax.plot(t_ms, E_kin, label="Kinetic")
    ax.plot(t_ms, E_tot, label="Total", linestyle="--", linewidth=1.5)
    ax.set_ylabel("Energy [J]")
    ax.set_title("Energies vs time")
    ax.legend()
    ax.grid(True)

    # 3) Currents
    ax = axes[2]
    for i in range(N):
        ax.plot(t_ms, I_hist[:, i] / 1e6, label=f"I_{i} ({types[i]}) [MA]")
    ax.set_ylabel("Current [MA]")
    ax.set_title("Coil currents vs time")
    ax.legend()
    ax.grid(True)

    # 4) r–z phase plot of plasma path
    ax = axes[3]
    ax.plot(z, r, "-o", markersize=2)
    ax.set_xlabel("z [m]")
    ax.set_ylabel("r [m]")
    ax.set_title("Plasma trajectory in (z,r) plane")
    ax.grid(True)
    ax.set_aspect("equal", adjustable="box")

    plt.tight_layout()


def main():
    parser = argparse.ArgumentParser(description="Dynamic plasma-coil demo with flux-conserving currents.")
    parser.add_argument( "--no-cage", action="store_true", help="Disable cage coil (simulate only seed + plasma)" )
    parser.add_argument( "--dt", type=float, default=5e-6, help="Time step [s] (default: 5e-6)")
    parser.add_argument( "--steps", type=int, default=400, help="Number of integration steps (default: 400)")
    parser.add_argument( "--P0", type=float, default=1e6, help="Initial effective plasma pressure [Pa] (default: 1e6)")
    parser.add_argument( "--gamma-gas", type=float, default=5.0 / 3.0, help="Polytropic index gamma for gas pressure (default: 5/3)", )
    parser.add_argument( "--L-eff", type=float, default=1.0, help="Effective axial length of plasma ring for pressure model (default: 1.0)", )
    parser.add_argument( "--mag-self-coeff", type=float, default=1e-7, help="Coefficient for magnetic self-pressure term ~ coeff * I^2 / r (default: 1e-7)", )
    args = parser.parse_args()

    res = run_plasma_dynamics( with_cage=not args.no_cage, dt=args.dt, n_steps=args.steps, P0=args.P0, gamma_gas=args.gamma_gas, L_eff=args.L_eff, mag_self_coeff=args.mag_self_coeff, )
    plot_results(res)
    plt.show()


if __name__ == "__main__":
    main()
