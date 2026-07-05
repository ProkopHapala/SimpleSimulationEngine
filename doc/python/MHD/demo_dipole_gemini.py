# === AUTO-DOC BEGIN ===
"""
@brief Simplest separatrix demo: single seed coil + point dipole = diamagnetic bubble.

Solves for the dipole moment m that creates Psi_total = 0 at a desired plasma radius,
then plots Psi contours (field lines). Uses the Psi = r*A_phi convention (not 2*pi*r),
so the separatrix is the Psi=0 contour. This is the minimal building block for
understanding diamagnetic plasma confinement before adding cage coils (see demo_dipole_gemini2).
"""
# === AUTO-DOC END ===

import argparse
import numpy as np
import matplotlib.pyplot as plt

from inductance_core import psi_loop_rz, psi_dipole_rz, MU0

# --- 2. SOLVER ---

def solve_dipole_strength(R_boundary, Z_boundary, I_seed, R_seed, Z_seed):
    """
    Finds dipole moment 'm' such that Psi_total(R_bound, Z_bound) = 0.
    Psi_seed + Psi_dipole = 0
    """
    psi_s = psi_loop_rz(R_seed, Z_seed, I_seed, np.array([R_boundary]), np.array([Z_boundary]))[0]
    psi_d_unit = psi_dipole_rz(1.0, 0.0, np.array([R_boundary]), np.array([Z_boundary]))[0]
    m_solution = -psi_s / psi_d_unit
    return m_solution

# --- 3. EXECUTION ---

def run_demo(args):
    R_SEED = args.r_seed
    Z_SEED = args.z_seed
    I_SEED = args.i_seed

    R_PLASMA = args.r_plasma
    Z_PLASMA = args.z_plasma

    m_sol = solve_dipole_strength(R_PLASMA, Z_PLASMA, I_SEED, R_SEED, Z_SEED)
    print(f"Required Dipole Moment: {m_sol:.2e} A*m^2")

    r_vals = np.linspace(0, 3.0, 200)
    z_vals = np.linspace(-2.0, 2.0, 200)
    R, Z = np.meshgrid(r_vals, z_vals)

    Psi_S = psi_loop_rz(R_SEED, Z_SEED, I_SEED, R, Z)
    Psi_D = psi_dipole_rz(m_sol, 0.0, R, Z)
    Psi_Tot = Psi_S + Psi_D

    plt.figure(figsize=(8, 10))

    plt.plot(R_SEED, Z_SEED, 'ro', markersize=10, label='Seed Coil')
    plt.plot(0, 0, 'b*', markersize=15, label='Plasma Dipole')

    levels = np.linspace(np.min(Psi_Tot), np.max(Psi_Tot), 50)
    plt.contour(R, Z, Psi_Tot, levels=levels, colors='gray', linewidths=0.5, alpha=0.5)

    plt.contour(R, Z, Psi_Tot, levels=[0], colors='blue', linewidths=3, linestyles='solid')

    plt.title(f"Rankine Body Method (Plasma Dipole)\nSeparatrix at R={R_PLASMA}")
    plt.xlabel("Radius R [m]")
    plt.ylabel("Axial Z [m]")
    plt.axis('equal')
    plt.grid(True)
    plt.legend()

    plt.text(0.1, 0.2, "Plasma Region\n(Diamagnetic)", color='blue', fontweight='bold')
    plt.text(1.5, 0.5, "External Field", color='gray')

    plt.tight_layout()
    if args.noshow:
        for n in plt.get_fignums():
            plt.figure(n); plt.savefig(f'demo_dipole_gemini_{n}.png', dpi=150, bbox_inches='tight')
        plt.close('all')
    else:
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Diamagnetic dipole bubble demo (Rankine body method)")
    parser.add_argument("--r-seed", type=float, default=2.0, help="Seed coil radius")
    parser.add_argument("--z-seed", type=float, default=0.0, help="Seed coil axial position")
    parser.add_argument("--i-seed", type=float, default=1.0e6, help="Seed coil current [A]")
    parser.add_argument("--r-plasma", type=float, default=0.8, help="Desired plasma bubble radius")
    parser.add_argument("--z-plasma", type=float, default=0.0, help="Plasma bubble axial position")
    parser.add_argument("--noshow", action="store_true", help="Save figures to PNG instead of displaying")
    args = parser.parse_args()
    run_demo(args)