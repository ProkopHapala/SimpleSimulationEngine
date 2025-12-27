"""
Hybrid Simulation: Superconducting Nozzle with a Diamagnetic Dipole Bubble.

This script simulates a plasma "bubble" (modeled as a Rankine body/Dipole separatrix)
moving through a flux-conserving metallic nozzle (Cage).

Physics Model:
--------------
1. Seed Coil (SC): Fixed current, creates background field.
2. Cage Coils: Passive rings. They trap the initial flux from the SC.
   As the plasma enters, they induce currents to maintain Phi = const.
3. Plasma: Modeled as a point dipole 'm' at position Z_p.
   The magnitude 'm' is NOT fixed. It is an unknown solved for at every step
   to ensure that the Magnetic Flux Function Psi_total = 0 at the desired
   plasma radius R_p (the 'Bubble Surface').

System of Equations (Coupled):
------------------------------
We solve for vector x = [I_cage_1, ..., I_cage_N, m_dipole].

[ K_cage      C_dip_cage ] [ I_cage ]   [ Phi_0 - Phi_sc_cage ]
[                        ] [        ] = [                     ]
[ C_cage_surf    C_m     ] [   m    ]   [   -Psi_sc_surf      ]

Where:
- K_cage: Mutual inductance matrix of cage rings.
- C_dip_cage: Flux through cage ring due to unit dipole.
- C_cage_surf: Flux function Psi at plasma surface due to unit cage current.
- C_m: Flux function Psi at plasma surface due to unit dipole.

Usage:
------
    python demo_nozzle_dipole_bubble.py

"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, hsv_to_rgb

from inductance_core import build_inductance_matrix, calc_mutual_inductance, Aphi_loop_rz, Aphi_dipole_rz, generate_parabolic_rings, MU0

# --- HELPER: Flux Function Psi ---
# Psi = 2 * pi * r * A_phi
# (Note: In axisymmetry, contours of Psi are field lines. Difference in Psi is Flux.)

def calc_Psi_loop(r_eval, z_eval, r_src, z_src, I_src):
    """Calculate Flux Function Psi at (r,z) due to a loop source."""
    # Aphi returns Vector Potential. Psi = 2*pi*r*Aphi
    A = Aphi_loop_rz(r_src, z_src, I_src, np.array([r_eval]), np.array([z_eval]))
    return 2 * np.pi * r_eval * A[0]

def calc_Psi_dipole(r_eval, z_eval, m_moment, z_dip):
    """Calculate Flux Function Psi at (r,z) due to a dipole source."""
    A = Aphi_dipole_rz(m_moment, z_dip, np.array([r_eval]), np.array([z_eval]))
    return 2 * np.pi * r_eval * A[0]

# --- SIMULATION CLASS ---

class DipoleNozzleSimulation:
    def __init__(self, args):
        self.args = args
        
        # 1. Generate Geometry
        # --------------------
        # Seed Coil
        self.sc = {'r': args.sc_r, 'z': args.sc_z, 'I': args.sc_current}
        
        # Cage Coils (Parabolic Nozzle) – follow the same wall profile as other nozzle demos
        cage_str = generate_parabolic_rings(n_rings=args.n_cage, r_throat=args.r_throat, r_exit=args.r_exit, z_start=args.z_start, z_end=args.z_end, coil_type="CAGE")
        self.cage_coils = []
        for line in cage_str.split('\n'):
            parts = line.split()
            if not parts: continue
            # Format: CAGE  R0  Z0  R1  Z1  0.0
            self.cage_coils.append({
                'r': float(parts[1]),
                'z': float(parts[2]),
                'I': 0.0,            # To be solved
                'phi_target': 0.0    # Conserved flux
            })
        
        self.n_cage = len(self.cage_coils)
        
        # 2. Precompute Cage Inductance Matrix (Fixed Geometry)
        # -----------------------------------------------------
        cage_rs = [c['r'] for c in self.cage_coils]
        cage_zs = [c['z'] for c in self.cage_coils]
        self.K_cage = build_inductance_matrix(cage_rs, cage_zs)
        
        # 3. Initialize Flux Conservation (t=0)
        # -------------------------------------
        # Assume at t=0, Cage currents are 0 (or response to SC).
        # We assume "Flux Freezing": The total flux through the cage rings
        # is whatever the SC puts there initially.
        for i in range(self.n_cage):
            # Flux from SC -> Cage Ring i
            M_sc_i = calc_mutual_inductance(self.sc['r'], self.sc['z'], cage_rs[i], cage_zs[i])
            phi_init = M_sc_i * self.sc['I']
            self.cage_coils[i]['phi_target'] = phi_init

        # History storage
        self.history = {
            'z_plasma': [], 'r_plasma': [],
            'm_dipole': [], 'I_cage_sum': [],
            'energies': [] 
        }

    def solve_step(self, z_p, r_p):
        """
        Solve for [I_cage, m_dipole] given the plasma position/radius.
        Constraints:
        1. Flux(Cage_i) = Phi_target_i
        2. Psi_total(at R_p, Z_p) = 0   (Separatrix condition)
        """
        N = self.n_cage
        
        # --- Build System Matrix [N+1 x N+1] ---
        # Block structure:
        # [ K_cage (NxN)      |  Col_dipole (Nx1) ]
        # [ Row_cage_surf (1xN)|  Val_dip_surf (1x1)]
        
        Matrix = np.zeros((N + 1, N + 1))
        RHS = np.zeros(N + 1)
        
        # Block 1: Cage-Cage couplings (Precomputed)
        Matrix[:N, :N] = self.K_cage
        
        # Block 2: Dipole -> Cage couplings (Flux through ring due to unit dipole)
        # Flux = Psi_dipole(r_ring, z_ring) * 2pi (already in calc_Psi)
        for i in range(N):
            c = self.cage_coils[i]
            # Coupling factor: Flux per unit moment
            C_dip_cage = calc_Psi_dipole(c['r'], c['z'], 1.0, z_p)
            Matrix[i, N] = C_dip_cage
            
            # RHS for Cage: Target Flux - Flux from SC
            # Phi_target - M_sc_cage * I_sc
            M_sc_c = calc_mutual_inductance(self.sc['r'], self.sc['z'], c['r'], c['z'])
            RHS[i] = c['phi_target'] - M_sc_c * self.sc['I']

        # Block 3: Cage -> Plasma Surface (Psi at bubble equator due to unit cage I)
        for j in range(N):
            c = self.cage_coils[j]
            # Psi per unit current
            C_cage_surf = calc_Psi_loop(r_p, z_p, c['r'], c['z'], 1.0)
            Matrix[N, j] = C_cage_surf
            
        # Block 4: Dipole -> Plasma Surface (Psi at bubble equator due to unit moment)
        # Note: Psi is singular at (0,0), but we evaluate at (r_p, 0) relative to dipole
        C_m_surf = calc_Psi_dipole(r_p, z_p, 1.0, z_p)
        Matrix[N, N] = C_m_surf
        
        # RHS for Plasma: -Psi from SC
        # We want Psi_total = 0 => Psi_cage + Psi_dip = -Psi_sc
        Psi_sc_surf = calc_Psi_loop(r_p, z_p, self.sc['r'], self.sc['z'], self.sc['I'])
        RHS[N] = -Psi_sc_surf
        
        # --- SOLVE ---
        Solution = np.linalg.solve(Matrix, RHS)
        
        # Distribute results
        I_cage_res = Solution[:N]
        m_dip_res = Solution[N]
        
        for i in range(N):
            self.cage_coils[i]['I'] = I_cage_res[i]
            
        return I_cage_res, m_dip_res

    def run(self):
        # Generate trajectory
        z_steps = np.linspace(self.args.z_p_start, self.args.z_p_end, self.args.steps)
        r_steps = np.linspace(self.args.r_p_start, self.args.r_p_end, self.args.steps)
        
        print(f"Running Simulation: {self.args.steps} steps...")
        
        for idx, (z, r) in enumerate(zip(z_steps, r_steps)):
            I_cage, m_dip = self.solve_step(z, r)
            
            # Log
            self.history['z_plasma'].append(z)
            self.history['r_plasma'].append(r)
            self.history['m_dipole'].append(m_dip)
            self.history['I_cage_sum'].append(np.sum(I_cage))
            
            # Visualize last step or requested frames
            if idx == self.args.steps - 1 and not self.args.no_plot:
                self.plot_state(z, r, m_dip, I_cage)

    def plot_state(self, z_p, r_p, m_dip, I_cage):
        """Visualize the magnetic streamlines (Psi contours) and geometry."""
        fig, ax = plt.subplots(figsize=(10, 7))
        
        # 1. Set dynamic extents from geometry/trajectory
        max_r_geom = max([c['r'] for c in self.cage_coils] + [self.sc['r'], r_p, self.args.r_p_end])
        r_max = 1.1 * max_r_geom
        z_min_geom = min(self.args.z_start, self.args.z_p_start, self.sc['z'])
        z_max_geom = max(self.args.z_end, self.args.z_p_end, self.sc['z'])
        dz_pad = 0.1 * (z_max_geom - z_min_geom + 1e-6)
        z_min = z_min_geom - dz_pad
        z_max = z_max_geom + dz_pad
        
        # 2. Grid (symmetric r for nicer fieldlines)
        r_grid = np.linspace(-r_max, r_max, 220)
        z_grid = np.linspace(z_min, z_max, 260)
        Zg, Rg = np.meshgrid(z_grid, r_grid)  # shape (Nr, Nz)
        Rabs = np.abs(Rg)
        
        # 3. Compute Flux Function Psi on grid
        Psi_tot = np.zeros_like(Rg)
        Psi_tot += calc_Psi_loop(Rabs, Zg, self.sc['r'], self.sc['z'], self.sc['I'])
        for c in self.cage_coils:
            Psi_tot += calc_Psi_loop(Rabs, Zg, c['r'], c['z'], c['I'])
        Psi_tot += calc_Psi_dipole(Rabs, Zg, m_dip, z_p)
        
        # 4. Derive B-field from Psi (B_r = -(1/r)dPsi/dz, B_z = (1/r)dPsi/dr)
        dPsi_dr, dPsi_dz = np.gradient(Psi_tot, r_grid, z_grid, edge_order=2)
        eps = 1e-9
        denom = np.where(np.abs(Rg) < eps, eps, Rg)
        Br = -dPsi_dz / denom
        Bz = dPsi_dr / denom
        Bmag = np.sqrt(Br * Br + Bz * Bz)
        B_ref = np.percentile(Bmag, 95)
        if B_ref <= 0.0:
            B_ref = 1.0
        B_norm = np.clip(Bmag / B_ref, 0.0, 1.0)
        hue = (np.arctan2(Bz, Br) + np.pi) / (2 * np.pi)
        sat = np.ones_like(hue)
        val = B_norm
        hsv = np.stack([hue, sat, val], axis=-1)
        rgb = hsv_to_rgb(hsv)
        
        # 5. Plot background and streamlines like the existing nozzle demo
        ax.imshow(rgb, extent=(z_min, z_max, -r_max, r_max), origin="lower", aspect="equal", alpha=0.8)
        ax.streamplot(z_grid, r_grid, Bz, Br, color="gray", linewidth=0.5, density=1.8, arrowsize=0.6)
        
        # Separatrix (Psi=0)
        ax.contour(Zg, Rg, Psi_tot, levels=[0], colors='blue', linewidths=2.0)
        
        # 6. Geometry overlays (mirror to both sides)
        cage_z = [c['z'] for c in self.cage_coils]
        cage_r = [c['r'] for c in self.cage_coils]
        norm = Normalize(vmin=min(I_cage), vmax=max(I_cage))
        cmap = plt.cm.plasma
        for sign in (+1, -1):
            sc = ax.scatter(cage_z, [sign * r for r in cage_r], c=I_cage, cmap=cmap, norm=norm, s=35, label='Cage Rings' if sign == 1 else None, zorder=5)
        plt.colorbar(sc, ax=ax, label='Induced Current [A]')
        for sign in (+1, -1):
            ax.plot(self.sc['z'], sign * self.sc['r'], 'ro', markersize=7, label='Seed Coil' if sign == 1 else None)
        ax.plot(z_p, 0, 'b*', markersize=12, label='Dipole Center')
        ax.plot([z_p], [r_p], 'bx', label='Target Radius')
        ax.plot([z_p], [-r_p], 'bx')
        
        # Trace the parabolic wall
        ax.plot(cage_z, cage_r, 'k--', alpha=0.4, lw=1.0, label='Nozzle wall (parabolic)')
        ax.plot(cage_z, [-r for r in cage_r], 'k--', alpha=0.4, lw=1.0)
        
        ax.set_title(f"Plasma Dipole Bubble in Flux-Conserving Nozzle\nZ_p={z_p:.2f}, R_p={r_p:.2f}, m={m_dip:.2e}")
        ax.set_xlabel("Z [m]")
        ax.set_ylabel("R [m]")
        ax.set_xlim(z_min, z_max)
        ax.set_ylim(-r_max, r_max)
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.25)
        ax.legend(loc='upper right', fontsize=8)
        
        plt.tight_layout()
        plt.show()

# --- MAIN ---

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulation of Diamagnetic Dipole Bubble in Flux-Conserving Nozzle")
    
    # Geometry
    # Match parabolic profile defaults to demo_plasma_sphere_flux for consistency
    parser.add_argument("--n-cage", type=int, default=10, help="Number of cage rings")
    parser.add_argument("--r-throat", type=float, default=0.05, help="Nozzle throat radius")
    parser.add_argument("--r-exit", type=float, default=1.4, help="Nozzle exit radius")
    parser.add_argument("--z-start", type=float, default=0.0, help="Nozzle start Z")
    parser.add_argument("--z-end", type=float, default=1.2, help="Nozzle end Z")
    
    # SC Coil
    parser.add_argument("--sc-r", type=float, default=1.5, help="Seed coil radius")
    parser.add_argument("--sc-z", type=float, default=0.0, help="Seed coil Z")
    parser.add_argument("--sc-current", type=float, default=1.0e6, help="Seed coil current")
    
    # Plasma Trajectory
    parser.add_argument("--z-p-start", type=float, default=0.0, help="Plasma start Z")
    parser.add_argument("--z-p-end",   type=float, default=0.2,  help="Plasma end Z")
    parser.add_argument("--r-p-start", type=float, default=0.3,  help="Plasma start Radius")
    parser.add_argument("--r-p-end",   type=float, default=1.5,  help="Plasma end Radius (expansion)")
    
    # Sim settings
    parser.add_argument("--steps", type=int, default=50,  help="Simulation steps")
    parser.add_argument("--no-plot", action="store_true", help="Skip plotting")
    
    args = parser.parse_args()
    
    sim = DipoleNozzleSimulation(args)
    sim.run()