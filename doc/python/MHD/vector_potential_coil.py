# === AUTO-DOC BEGIN ===
"""
@brief Vector potential A_phi of a single circular loop: validates B = curl(A) numerically.

Computes A_phi via inductance_core.Aphi_loop_rz, then derives Br/Bz via finite-difference
curl (compute_numerical_B_from_A from inductance_core) and compares against the analytic
field_loop_rz formula. Plots A_phi(r) and overlaid analytic vs numerical B components at
a fixed z. This is the fundamental validation that A_phi and B kernels are consistent.
"""
# === AUTO-DOC END ===

import argparse
import numpy as np
import matplotlib.pyplot as plt

from inductance_core import Aphi_loop_rz, field_loop_rz, compute_numerical_B_from_A

# ==========================================
# RUN TEST
# ==========================================

def run_demo(args):
    R_COIL = args.r_coil
    I_COIL = args.i_coil
    Z_COIL = args.z_coil
    Z_EVAL = args.z_eval

    r_eval = np.linspace(0.0, 3.0, 100)

    r_mesh, z_mesh = np.meshgrid(r_eval, [Z_EVAL])
    Br_ana, Bz_ana = field_loop_rz(R_COIL, Z_COIL, I_COIL, r_mesh, z_mesh)
    Br_ana = Br_ana[0]
    Bz_ana = Bz_ana[0]

    Br_num, Bz_num = compute_numerical_B_from_A(R_COIL, Z_COIL, I_COIL, r_eval, Z_EVAL)

    A_val = Aphi_loop_rz(R_COIL, Z_COIL, I_COIL, r_mesh, z_mesh)[0]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10), sharex=True)

    ax1.plot(r_eval, A_val, 'k-', lw=2, label=r'$A_\phi$')
    ax1.set_ylabel(r'Vector Potential $A_\phi$ [T$\cdot$m]')
    ax1.set_title(f'Vector Potential at Z = {Z_EVAL} (Coil R={R_COIL})')
    ax1.grid(True)
    ax1.legend()

    ax2.plot(r_eval, Bz_ana, 'b-', lw=3, alpha=0.5, label=r'$B_z$ (Analytic)')
    ax2.plot(r_eval, Bz_num, 'b--', label=r'$B_z$ (Num Curl)')

    ax2.plot(r_eval, Br_ana, 'r-', lw=3, alpha=0.5, label=r'$B_r$ (Analytic)')
    ax2.plot(r_eval, Br_num, 'r--', label=r'$B_r$ (Num Curl)')

    ax2.set_xlabel('Radius r [m]')
    ax2.set_ylabel('Magnetic Field [T]')
    ax2.set_title('Validation: B = curl(A) vs Analytic B')
    ax2.legend()
    ax2.grid(True)

    err_z = np.mean(np.abs(Bz_num - Bz_ana))
    err_r = np.mean(np.abs(Br_num - Br_ana))
    print(f"Mean Absolute Error Bz: {err_z:.2e} T")
    print(f"Mean Absolute Error Br: {err_r:.2e} T")

    plt.tight_layout()
    if args.noshow:
        for n in plt.get_fignums():
            plt.figure(n); plt.savefig(f'vector_potential_coil_{n}.png', dpi=150, bbox_inches='tight')
        plt.close('all')
    else:
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Vector potential of a circular coil: validation B = curl(A)")
    parser.add_argument("--r-coil", type=float, default=1.0, help="Coil radius")
    parser.add_argument("--i-coil", type=float, default=1.0e6, help="Coil current [A]")
    parser.add_argument("--z-coil", type=float, default=0.0, help="Coil axial position")
    parser.add_argument("--z-eval", type=float, default=0.5, help="Z position for evaluation line")
    parser.add_argument("--noshow", action="store_true", help="Save figures to PNG instead of displaying")
    args = parser.parse_args()
    run_demo(args)