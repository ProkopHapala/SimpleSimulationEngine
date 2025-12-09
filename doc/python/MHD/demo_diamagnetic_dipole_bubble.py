import argparse
import numpy as np
import matplotlib.pyplot as plt

from inductance_core import field_loop_rz, field_dipole_rz, Aphi_loop_rz, Aphi_dipole_rz, MU0


def compute_fields_coil_plus_dipole(r_grid, z_grid, a_coil, z_coil, I_coil, m_dipole, z_dipole):
    """Compute Br,Bz for a single circular coil plus an axial point dipole.

    Parameters
    ----------
    r_grid, z_grid : ndarray
        Meshgrid in (r,z) where the field is evaluated.
    a_coil, z_coil : float
        Radius and axial position of the superconducting coil.
    I_coil : float
        Coil current (A).
    m_dipole, z_dipole : float
        Dipole moment (A m^2) and axial position of the point dipole.
    """
    Br_c, Bz_c = field_loop_rz(a_coil, z_coil, I_coil, r_grid, z_grid)
    Br_d, Bz_d = field_dipole_rz(m_dipole, z_dipole, r_grid, z_grid)
    Br = Br_c + Br_d
    Bz = Bz_c + Bz_d
    return Br, Bz, Br_c, Bz_c, Br_d, Bz_d


def make_background_rgb(Br, Bz, bg_mode="hsv", B_ref=None):
    """Construct RGB background from Br,Bz as in MHD_plots.plot_coil_geometry."""
    Bmag = np.sqrt(Br * Br + Bz * Bz)
    if B_ref is None or B_ref <= 0.0:
        B_ref = 0.01 * np.max(Bmag)
        if B_ref <= 0.0:
            B_ref = 1.0
    Bmag_norm = np.clip(Bmag / B_ref, 0.0, 1.0)
    phi = np.arctan2(Bz, Br)

    if bg_mode.lower() == "hsv":
        from matplotlib.colors import hsv_to_rgb

        hue = (phi + np.pi) / (2 * np.pi)
        sat = np.ones_like(hue)
        val = Bmag_norm
        hsv = np.stack([hue, sat, val], axis=-1)
        rgb = hsv_to_rgb(hsv)
    else:
        cmap = plt.cm.magma
        rgb = cmap(Bmag_norm)[..., :3]
    return rgb


def check_Aphi_vs_B(a_coil, z_coil, I_coil, m_dipole, z_dipole):
    """Numerically verify that curl(Aphi) reproduces Br,Bz at one sample point."""
    r0 = 0.7 * a_coil if a_coil > 0.0 else 0.5
    z0 = z_coil + 0.5 * a_coil
    dr = 1e-3 * max(a_coil, 1.0)
    dz = 1e-3 * max(a_coil, 1.0)

    def Aphi_tot(r, z):
        A_c = Aphi_loop_rz(a_coil, z_coil, I_coil, r, z)
        A_d = Aphi_dipole_rz(m_dipole, z_dipole, r, z)
        return A_c + A_d

    A0 = Aphi_tot(r0, z0)
    Azp = Aphi_tot(r0, z0 + dz)
    Azm = Aphi_tot(r0, z0 - dz)
    Arp = Aphi_tot(r0 + dr, z0)
    Arm = Aphi_tot(r0 - dr, z0)

    Br_num = -(Azp - Azm) / (2.0 * dz)
    Bz_num = ( ((r0 + dr) * Arp - (r0 - dr) * Arm) / (2.0 * dr) ) / r0

    Br_c, Bz_c = field_loop_rz(a_coil, z_coil, I_coil, np.array([[r0]]), np.array([[z0]]))
    Br_d, Bz_d = field_dipole_rz(m_dipole, z_dipole, np.array([[r0]]), np.array([[z0]]))
    Br_ref = (Br_c + Br_d)[0, 0]
    Bz_ref = (Bz_c + Bz_d)[0, 0]

    err_Br = abs(Br_num - Br_ref) / (abs(Br_ref) + 1e-12)
    err_Bz = abs(Bz_num - Bz_ref) / (abs(Bz_ref) + 1e-12)
    print(f"[check_Aphi_vs_B] sample (r,z)=({r0:.3g},{z0:.3g}), |rel err Br|={err_Br:.3e}, |rel err Bz|={err_Bz:.3e}")


def run_Aphi_profile_tests(a_coil, z_coil, I_coil, m_dipole, z_dipole, r_max, n_r):
    """Systematic radial scans: compare Br,Bz from direct formulas vs curl(Aphi)."""
    # Radial line (avoid r=0 exactly)
    r = np.linspace(1e-4 * max(a_coil, 1.0), r_max, n_r)
    dr = r[1] - r[0]
    dz = 1e-3 * max(a_coil, 1.0)

    # Choose an offset in front of the coil for the test plane
    z0 = z_coil + 0.5 * a_coil

    # --- Coil-only ---
    A_c = Aphi_loop_rz(a_coil, z_coil, I_coil, r, z0)
    A_c_zp = Aphi_loop_rz(a_coil, z_coil, I_coil, r, z0 + dz)
    A_c_zm = Aphi_loop_rz(a_coil, z_coil, I_coil, r, z0 - dz)
    A_c_rp = Aphi_loop_rz(a_coil, z_coil, I_coil, r + dr, z0)
    A_c_rm = Aphi_loop_rz(a_coil, z_coil, I_coil, r - dr, z0)

    Br_c_num = -(A_c_zp - A_c_zm) / (2.0 * dz)
    Bz_c_num = (((r + dr) * A_c_rp - (r - dr) * A_c_rm) / (2.0 * dr)) / r

    Rg_c = r[None, :]
    Zg_c = np.full_like(Rg_c, z0)
    Br_c_ref, Bz_c_ref = field_loop_rz(a_coil, z_coil, I_coil, Rg_c, Zg_c)
    Br_c_ref = Br_c_ref[0, :]
    Bz_c_ref = Bz_c_ref[0, :]

    # --- Dipole-only ---
    A_d = Aphi_dipole_rz(m_dipole, z_dipole, r, z0)
    A_d_zp = Aphi_dipole_rz(m_dipole, z_dipole, r, z0 + dz)
    A_d_zm = Aphi_dipole_rz(m_dipole, z_dipole, r, z0 - dz)
    A_d_rp = Aphi_dipole_rz(m_dipole, z_dipole, r + dr, z0)
    A_d_rm = Aphi_dipole_rz(m_dipole, z_dipole, r - dr, z0)

    Br_d_num = -(A_d_zp - A_d_zm) / (2.0 * dz)
    Bz_d_num = (((r + dr) * A_d_rp - (r - dr) * A_d_rm) / (2.0 * dr)) / r

    Rg_d = r[None, :]
    Zg_d = np.full_like(Rg_d, z0)
    Br_d_ref, Bz_d_ref = field_dipole_rz(m_dipole, z_dipole, Rg_d, Zg_d)
    Br_d_ref = Br_d_ref[0, :]
    Bz_d_ref = Bz_d_ref[0, :]

    fig, axes = plt.subplots(2, 2, figsize=(10, 6), sharex=True)

    axes[0, 0].plot(r, Br_c_ref, "k-", label="Br coil (direct)")
    axes[0, 0].plot(r, Br_c_num, "r--", label="Br coil (curl A)")
    axes[0, 0].set_ylabel("Br [T]")
    axes[0, 0].legend(); axes[0, 0].grid(True, linestyle=":", alpha=0.5)

    axes[1, 0].plot(r, Bz_c_ref, "k-", label="Bz coil (direct)")
    axes[1, 0].plot(r, Bz_c_num, "r--", label="Bz coil (curl A)")
    axes[1, 0].set_xlabel("r [m]")
    axes[1, 0].set_ylabel("Bz [T]")
    axes[1, 0].legend(); axes[1, 0].grid(True, linestyle=":", alpha=0.5)

    axes[0, 1].plot(r, Br_d_ref, "k-", label="Br dipole (direct)")
    axes[0, 1].plot(r, Br_d_num, "r--", label="Br dipole (curl A)")
    axes[0, 1].set_ylabel("Br [T]")
    axes[0, 1].legend(); axes[0, 1].grid(True, linestyle=":", alpha=0.5)

    axes[1, 1].plot(r, Bz_d_ref, "k-", label="Bz dipole (direct)")
    axes[1, 1].plot(r, Bz_d_num, "r--", label="Bz dipole (curl A)")
    axes[1, 1].set_xlabel("r [m]")
    axes[1, 1].set_ylabel("Bz [T]")
    axes[1, 1].legend(); axes[1, 1].grid(True, linestyle=":", alpha=0.5)

    fig.suptitle(f"Aphi consistency test at z0={z0:.3f}")
    plt.tight_layout()
    plt.show()


def run_demo(args):
    # Construct grid in (r,z)
    r = np.linspace(-args.r_max, args.r_max, args.n_r)
    z = np.linspace(args.z_min, args.z_max, args.n_z)
    Zg, Rg = np.meshgrid(z, r)

    # Compute fields
    Br, Bz, Br_c, Bz_c, Br_d, Bz_d = compute_fields_coil_plus_dipole(
        np.abs(Rg), Zg,
        a_coil=args.r_coil,
        z_coil=args.z_coil,
        I_coil=args.i_coil,
        m_dipole=args.m_dipole,
        z_dipole=args.z_dipole,
    )

    # Background coloring
    rgb = make_background_rgb(Br, Bz, bg_mode=args.bg_mode, B_ref=args.B_ref)

    fig, ax = plt.subplots(figsize=(8, 6))
    extent = (args.z_min, args.z_max, -args.r_max, args.r_max)
    ax.imshow(rgb, extent=extent, origin="lower", aspect="equal", alpha=0.8)

    # Streamlines of total field
    ax.streamplot(z, r, Bz, Br * np.sign(Rg), color="gray", linewidth=0.6, density=1.8, arrowsize=0.6)

    # Optional: Aphi-based separatrix via psi = r Aphi contour
    if getattr(args, "highlight_separatrix", False):
        Aphi_c = Aphi_loop_rz(args.r_coil, args.z_coil, args.i_coil, np.abs(Rg), Zg)
        Aphi_d = Aphi_dipole_rz(args.m_dipole, args.z_dipole, np.abs(Rg), Zg)
        Aphi = Aphi_c + Aphi_d
        psi = np.abs(Rg) * Aphi

        ir0 = int(np.argmin(np.abs(r)))
        Bz_axis = Bz[ir0, :]
        mask = z > args.z_coil
        z_fwd = z[mask]
        Bz_fwd = Bz_axis[mask]
        z_sep = None
        if z_fwd.size > 1:
            for k in range(z_fwd.size - 1):
                b1, b2 = Bz_fwd[k], Bz_fwd[k + 1]
                if b1 == 0.0 or b1 * b2 < 0.0:
                    z1, z2 = z_fwd[k], z_fwd[k + 1]
                    t = -b1 / (b2 - b1) if b2 != b1 else 0.5
                    z_sep = z1 + t * (z2 - z1)
                    break
        if z_sep is not None:
            r_seed = 0.1 * args.r_coil
            ir_seed = int(np.argmin(np.abs(r - r_seed)))
            iz_seed = int(np.argmin(np.abs(z - z_sep)))
            psi_sep = psi[ir_seed, iz_seed]
            ax.contour(Zg, Rg, psi, levels=[psi_sep], colors="cyan", linewidths=1.5)

    # Mark coil (SC loop) at (r=±a_coil, z=z_coil)
    for sign in (+1, -1):
        ax.plot(args.z_coil, sign * args.r_coil, marker="o", color="k", markersize=5)
    I_MA = args.i_coil / 1e6
    ax.text(args.z_coil, args.r_coil, f"SC: {I_MA:.3f} MA", ha="center", va="bottom")

    # Mark dipole location on axis
    ax.plot(args.z_dipole, 0.0, marker="^", color="tab:red", markersize=6)
    ax.text(args.z_dipole, 0.0, "dipole", ha="center", va="bottom", color="tab:red")

    ax.set_xlabel("z [m]")
    ax.set_ylabel("r [m]")
    ax.set_xlim(args.z_min, args.z_max)
    ax.set_ylim(-args.r_max, args.r_max)
    ax.set_aspect("equal", "box")
    ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.5)
    ax.set_title("SC coil + diamagnetic dipole bubble (total B field)")

    # Optional: contour where Bz changes sign (proxy for plasma boundary)
    if args.show_Bz_zero:
        levels = [0.0]
        ax.contour(Zg, Rg, Bz, levels=levels, colors="white", linewidths=1.0)
        ax.plot([], [], color="white", linewidth=1.0, label="Bz = 0")
        ax.legend(loc="best")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=(
        "Static diamagnetic bubble demo: field of one superconducting coil plus an axial magnetic dipole, "
        "visualized in (r,z) using the same style as MHD_plots."
    ))
    parser.add_argument("--r-coil",   type=float, default=1.0, help="Radius of superconducting coil")
    parser.add_argument("--z-coil",   type=float, default=0.0, help="Axial position of superconducting coil")
    parser.add_argument("--i-coil",   type=float, default=1.0e6, help="Current in superconducting coil [A]")
    parser.add_argument("--m-dipole", type=float, default=-0.5e6, help="Magnetic dipole moment [A m^2] (along +z)")
    parser.add_argument("--z-dipole", type=float, default=1.0, help="Axial position of dipole")
    parser.add_argument("--z-min",    type=float, default=-1.0, help="Minimum z for plot window")
    parser.add_argument("--z-max",    type=float, default=3.0, help="Maximum z for plot window")
    parser.add_argument("--r-max",    type=float, default=2.0, help="Maximum |r| for plot window")
    parser.add_argument("--n-r",      type=int,   default=160, help="Number of radial samples")
    parser.add_argument("--n-z",      type=int,   default=320, help="Number of axial samples")

    parser.add_argument("--bg-mode", choices=["hsv", "mag"], default="mag", help="Background coloring mode")
    parser.add_argument("--B-ref", type=float, default=1.0, help="Reference |B| for background normalization [T]; if <=0, use auto (0.01*max)")
    parser.add_argument("--show-Bz-zero",         type=int, default=1, help="Overlay contour where Bz=0 as proxy for plasma boundary")
    parser.add_argument("--highlight-separatrix", type=int, default=1, help="Highlight separatrix via psi=r*Aphi contour from on-axis Bz=0 in front of coil")
    parser.add_argument("--verify-A",             type=int, default=1, help="Numerically verify that curl(Aphi) matches Br,Bz at one sample point")
    parser.add_argument("--test-Aphi-profiles",   type=int, default=1, help="Run radial profile tests comparing direct B vs curl(Aphi) for coil and dipole")

    args = parser.parse_args()
    if args.verify_A:
        check_Aphi_vs_B(args.r_coil, args.z_coil, args.i_coil, args.m_dipole, args.z_dipole)
    if args.test_Aphi_profiles:
        run_Aphi_profile_tests(args.r_coil, args.z_coil, args.i_coil, args.m_dipole, args.z_dipole, args.r_max, args.n_r)
    run_demo(args)
