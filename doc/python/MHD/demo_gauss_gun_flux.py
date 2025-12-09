import argparse
import numpy as np
import matplotlib.pyplot as plt

from inductance_core import build_inductance_matrix, magnetic_energy
from MHD_plots import plot_coil_geometry


def build_driver_tube(n_drive=4, r_drive=1.0, z_start=0.0, z_end=None, i_drive=1.0e6):
    """Generate driver tube coil positions and initial currents.

    By default, coils are spaced roughly by one barrel radius along z.

    Returns
    -------
    R_drive, Z_drive, I_drive : 1D ndarrays of length n_drive
    """
    if z_end is None:
        # Space coils by approximately one barrel radius
        length = r_drive * max(1, (n_drive - 1))
        z_end = z_start + length
    Z_drive = np.linspace(z_start, z_end, n_drive)
    R_drive = np.full_like(Z_drive, float(r_drive))
    I_drive = np.full_like(Z_drive, float(i_drive))
    return R_drive, Z_drive, I_drive


def build_projectile_disk(n_disk=1, r_inner=0.9, r_outer=0.9):
    """Generate projectile disk ring radii.

    By default, outer projectile radius is ~0.9 * barrel radius (assuming
    r_drive=1.0). Radii are ordered from outermost to innermost so proj_1 is
    the largest ring.

    Returns
    -------
    r_disk : 1D ndarray of length n_disk
    """
    # Start from outer radius and go inward so the first index is the largest ring
    if n_disk <= 1:
        return np.array([r_outer], dtype=float)
    return np.linspace(r_outer, r_inner, n_disk)


def build_geometry_vectors(R_drive, Z_drive, r_disk, z_proj):
    """Assemble global radius/axial arrays for driver + projectile loops."""
    n_drive = R_drive.size
    n_disk = r_disk.size
    rs = np.empty(n_drive + n_disk, dtype=float)
    zs = np.empty_like(rs)
    rs[:n_drive] = R_drive
    zs[:n_drive] = Z_drive
    rs[n_drive:] = r_disk
    zs[n_drive:] = z_proj
    return rs, zs


def compute_initial_flux(R_drive, Z_drive, I_drive0, r_disk, z_proj0):
    """Compute initial flux vector Phi0 for all loops (drivers + projectile).

    At t=0, driver coils carry I_drive0, projectile loops carry zero current.
    """
    rs0, zs0 = build_geometry_vectors(R_drive, Z_drive, r_disk, z_proj0)
    K0 = build_inductance_matrix(rs0, zs0)
    I0 = np.zeros_like(rs0)
    I0[: R_drive.size] = I_drive0
    Phi0 = K0 @ I0
    return Phi0, K0, I0


def solve_active_currents(R_drive, Z_drive, r_disk, z_proj, Phi0, active_drive):
    """Solve for currents in the active flux-conserving set.

    Drivers with active_drive[k] == True and all projectile loops are included
    in the system K_active @ I_active = Phi0_active. Switched-off drivers are
    removed from this system and have I=0.
    """
    n_drive = R_drive.size
    n_disk = r_disk.size
    n_total = n_drive + n_disk

    # Global geometry vectors (for convenience)
    rs_full, zs_full = build_geometry_vectors(R_drive, Z_drive, r_disk, z_proj)

    # Active indices: all still-connected drivers + all projectile loops
    idx_drive_active = np.where(active_drive)[0]
    idx_proj = np.arange(n_drive, n_total)
    idx_active = np.concatenate([idx_drive_active, idx_proj])

    # Build reduced system on active subset
    rs_act = rs_full[idx_active]
    zs_act = zs_full[idx_active]
    K_act = build_inductance_matrix(rs_act, zs_act)
    Phi0_act = Phi0[idx_active]
    I_act = np.linalg.solve(K_act, Phi0_act)

    # Scatter back to full vector and build full K for energy/diagnostics
    I_full = np.zeros(n_total, dtype=float)
    I_full[idx_active] = I_act
    K_full = build_inductance_matrix(rs_full, zs_full)
    return I_full, K_full


def compute_force_z_energy_gradient(R_drive, Z_drive, r_disk, z_proj, Phi0, active_drive, dz_fd=1e-4):
    """Compute axial force on projectile via energy finite difference.

    Fz ≈ -[W(z+dz) - W(z-dz)] / (2 dz)
    with global flux invariants Phi0 and a shrinking active set of flux-conserving coils.
    """
    # z - dz
    I_minus, K_minus = solve_active_currents(R_drive, Z_drive, r_disk, z_proj - dz_fd, Phi0, active_drive)
    W_minus = magnetic_energy(K_minus, I_minus)

    # z + dz
    I_plus, K_plus = solve_active_currents(R_drive, Z_drive, r_disk, z_proj + dz_fd, Phi0, active_drive)
    W_plus = magnetic_energy(K_plus, I_plus)

    Fz = -(W_plus - W_minus) / (2.0 * dz_fd)
    return Fz, (W_minus, W_plus)


def apply_switching_logic(z_proj, Z_drive, active_drive, a_switch):
    """Mark nearest driver coil ahead of projectile as switched-off if within a_switch.

    Updates active_drive in-place; deactivated coils are removed from the
    flux-conserving system and their currents remain zero afterwards.

    Returns index of switched-off coil or None if no switch occurred.
    """
    n_drive = Z_drive.size
    dz = Z_drive - z_proj
    mask_ahead = (dz > 0.0) & active_drive
    if not np.any(mask_ahead):
        return None
    idx = np.where(mask_ahead)[0]
    k_next = idx[np.argmin(dz[idx])]
    if 0.0 < dz[k_next] < a_switch:
        active_drive[k_next] = False
        return int(k_next)
    return None


def run_gauss_gun_sim(args):
    # Geometry and initial currents
    R_drive, Z_drive, I_drive0 = build_driver_tube(n_drive=args.n_drive, r_drive=args.r_drive, z_start=args.z_start, z_end=args.z_end, i_drive=args.i_drive)
    r_disk = build_projectile_disk(n_disk=args.n_disk, r_inner=args.r_inner, r_outer=args.r_outer)

    # Global flux invariants for all loops at t=0
    Phi0, K0, I0 = compute_initial_flux(R_drive, Z_drive, I_drive0, r_disk, args.z0)

    # State variables
    n_drive = R_drive.size
    n_disk = r_disk.size
    n_total = n_drive + n_disk
    active_drive = np.ones(n_drive, dtype=bool)

    z_prev = args.z0 - args.v0 * args.dt  # simple bootstrap without initial force term
    z_curr = args.z0

    n_steps = args.steps
    ts = np.zeros(n_steps + 1, dtype=float)
    zs = np.zeros(n_steps + 1, dtype=float)
    vs = np.zeros(n_steps + 1, dtype=float)
    Fs = np.zeros(n_steps + 1, dtype=float)
    Ws = np.zeros(n_steps + 1, dtype=float)
    KEs = np.zeros(n_steps + 1, dtype=float)   # kinetic energy of projectile
    Etot = np.zeros(n_steps + 1, dtype=float)  # total = magnetic + kinetic

    disk_currents = np.zeros((n_steps + 1, n_disk), dtype=float)
    drive_currents = np.zeros((n_steps + 1, n_drive), dtype=float)
    switch_events = []  # list of (t, z_proj, coil_index)

    # Initial step: solve for currents in the full active set at z_curr
    I_init, K_init = solve_active_currents(R_drive, Z_drive, r_disk, z_curr, Phi0, active_drive)
    W0 = magnetic_energy(K_init, I_init)
    F0, _ = compute_force_z_energy_gradient(R_drive, Z_drive, r_disk, z_curr, Phi0, active_drive, dz_fd=args.dz_fd)

    ts[0] = 0.0
    zs[0] = z_curr
    vs[0] = args.v0
    Fs[0] = F0
    Ws[0] = W0
    KEs[0] = 0.5 * args.m_proj * args.v0 * args.v0
    Etot[0] = Ws[0] + KEs[0]
    disk_currents[0, :] = I_init[n_drive:]
    drive_currents[0, :] = I_init[:n_drive]

    # Correct bootstrap for z_prev using initial acceleration
    z_prev = z_curr - args.v0 * args.dt + 0.5 * (F0 / args.m_proj) * args.dt * args.dt

    for i in range(1, n_steps + 1):
        t = i * args.dt

        # Apply coil switching logic based on current projectile position
        k_switched = apply_switching_logic(z_curr, Z_drive, active_drive, args.a_switch)
        if k_switched is not None:
            switch_events.append((t, float(z_curr), int(k_switched)))

        # Solve for currents in active set and total energy
        I_all, K = solve_active_currents(R_drive, Z_drive, r_disk, z_curr, Phi0, active_drive)
        W = magnetic_energy(K, I_all)

        # Force from energy gradient
        Fz, _ = compute_force_z_energy_gradient(R_drive, Z_drive, r_disk, z_curr, Phi0, active_drive, dz_fd=args.dz_fd)

        # Verlet position update
        z_next = 2.0 * z_curr - z_prev + (Fz / args.m_proj) * args.dt * args.dt

        # Approximate velocity for diagnostics
        v_curr = (z_next - z_prev) / (2.0 * args.dt)
        KE_curr = 0.5 * args.m_proj * v_curr * v_curr

        ts[i] = t
        zs[i] = z_curr
        vs[i] = v_curr
        Fs[i] = Fz
        Ws[i] = W
        KEs[i] = KE_curr
        Etot[i] = W + KE_curr
        disk_currents[i, :] = I_all[n_drive:]
        drive_currents[i, :] = I_all[:n_drive]

        # Advance steps
        z_prev, z_curr = z_curr, z_next

        # Stop if projectile has left the driver region by a margin
        if z_curr > (Z_drive.max() + args.exit_margin):
            ts = ts[: i + 1]
            zs = zs[: i + 1]
            vs = vs[: i + 1]
            Fs = Fs[: i + 1]
            Ws = Ws[: i + 1]
            KEs = KEs[: i + 1]
            Etot = Etot[: i + 1]
            disk_currents = disk_currents[: i + 1, :]
            drive_currents = drive_currents[: i + 1, :]
            break

    return {
        "t": ts,
        "z": zs,
        "v": vs,
        "Fz": Fs,
        "W": Ws,
        "K": KEs,
        "E_tot": Etot,
        "disk_currents": disk_currents,
        "drive_currents": drive_currents,
        "R_drive": R_drive,
        "Z_drive": Z_drive,
        "r_disk": r_disk,
        "switch_events": switch_events,
    }


def plot_gauss_gun_results(res, args):
    t = res["t"]
    z = res["z"]
    v = res["v"]
    Fz = res["Fz"]
    W = res["W"]
    K = res.get("K")
    E_tot = res.get("E_tot")
    disk_currents = res["disk_currents"]
    drive_currents = res["drive_currents"]
    R_drive = res["R_drive"]
    Z_drive = res["Z_drive"]
    r_disk = res["r_disk"]

    fig1, ax = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
    ax[0].plot(t, z)
    ax[0].set_ylabel("z_proj [m]")
    ax[0].grid(True, alpha=0.3)

    ax[1].plot(t, v)
    ax[1].set_ylabel("v_proj [m/s]")
    ax[1].grid(True, alpha=0.3)

    ax[2].plot(t, Fz)
    ax[2].set_ylabel("Fz [N]")
    ax[2].set_xlabel("t [s]")
    ax[2].grid(True, alpha=0.3)

    fig1.tight_layout()

    # Energy (magnetic, kinetic, total) + combined coil currents (drivers + projectile)
    fig2, ax2 = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    ax2[0].plot(t, W, label="W_mag")
    if K is not None:
        ax2[0].plot(t, K, label="K_proj")
    if E_tot is not None:
        ax2[0].plot(t, E_tot, label="E_tot")
    ax2[0].set_ylabel("Energy [J]")
    ax2[0].grid(True, alpha=0.3)
    ax2[0].legend(loc="best")

    ax2[1].set_xlabel("t [s]")
    ax2[1].set_ylabel("I [MA]")
    ax2[1].grid(True, alpha=0.3)

    # Driver coils: ring_i
    n_drive = drive_currents.shape[1]
    for k in range(n_drive):
        ax2[1].plot(t, drive_currents[:, k] * 1e-6, label=f"ring_{k+1}")

    # Projectile coils: proj_j
    n_disk = r_disk.size
    for j in range(n_disk):
        ax2[1].plot(t, disk_currents[:, j] * 1e-6, linestyle="--", label=f"proj_{j+1}")

    if n_drive + n_disk <= 16:
        ax2[1].legend(loc="best")

    fig2.tight_layout()

    # Geometry + B-field snapshot using shared MHD_plots utility
    n_drive = drive_currents.shape[1]
    n_disk = r_disk.size
    n_total = n_drive + n_disk

    types = np.empty(n_total, dtype=object)
    types[:n_drive] = "CAGE"   # driver tube rings
    types[n_drive:] = "PLASMA"  # projectile disk rings

    R0 = np.concatenate([R_drive, r_disk])
    z0_arr = np.concatenate([Z_drive, np.full(n_disk, z[0], dtype=float)])
    I0 = np.concatenate([drive_currents[0, :], disk_currents[0, :]])

    R1 = R0.copy()  # radii do not change in this demo
    z1_arr = np.concatenate([Z_drive, np.full(n_disk, z[-1], dtype=float)])
    I_final = np.concatenate([drive_currents[-1, :], disk_currents[-1, :]])

    plot_coil_geometry(types, R0, z0_arr, I0, R1, z1_arr, I_final, bg_mode="mag", title="Gauss-gun: initial vs final geometry")

    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=("Gauss-gun style coil accelerator: projectile disk accelerated by pre-magnetized driver tube with flux-conserving projectile loops."))
    parser.add_argument("--n-drive", type=int,   default=4, help="Number of driver coils along z")
    parser.add_argument("--r-drive", type=float, default=0.5, help="Radius of driver coils")
    parser.add_argument("--z-start", type=float, default=0.0, help="Start z of driver tube")
    parser.add_argument("--z-end",   type=float, default=2.0, help="End z of driver tube")
    parser.add_argument("--i-drive", type=float, default=1.0e6, help="Current in active driver coils [A]")

    parser.add_argument("--n-disk", type=int, default=1, help="Number of rings in projectile disk")
    parser.add_argument("--r-inner", type=float, default=0.1, help="Inner radius of projectile disk")
    parser.add_argument("--r-outer", type=float, default=0.45, help="Outer radius of projectile disk")
    parser.add_argument("--m-proj", type=float, default=0.1, help="Projectile mass [kg]")

    parser.add_argument("--z0", type=float, default=-0.05, help="Initial projectile z position")
    parser.add_argument("--v0", type=float, default=5.0, help="Initial projectile axial velocity")

    parser.add_argument("--a-switch", type=float, default=0.08, help="Distance ahead of projectile to switch off next driver coil")
    parser.add_argument("--dt", type=float, default=5e-6, help="Time step [s]")
    parser.add_argument("--steps", type=int, default=4000, help="Maximum number of time steps")
    parser.add_argument("--dz-fd", type=float, default=1e-4, help="Finite-difference step for force evaluation")
    parser.add_argument("--exit-margin", type=float, default=0.5, help="Margin beyond last driver coil to stop simulation")

    args = parser.parse_args()

    res = run_gauss_gun_sim(args)

    # Debug summary: switching and trajectory range
    switches = res.get("switch_events", [])
    print(f"Number of switch events: {len(switches)}")
    for (t_sw, z_sw, k_sw) in switches:
        print(f"  switch: t={t_sw:.6e} s, z={z_sw:.6e} m, coil_index={k_sw}")

    z = res["z"]
    print(f"z range: z_min={z.min():.6e} m, z_max={z.max():.6e} m")

    plot_gauss_gun_results(res, args)
