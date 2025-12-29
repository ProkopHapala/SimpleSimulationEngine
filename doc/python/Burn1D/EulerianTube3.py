import numpy as np
import matplotlib.pyplot as plt
import argparse
import time

# --- Configuration & Physics Constants ---
GAMMA = 1.4      # Adiabatic index (Air)
R_GAS = 287.0    # Specific gas constant (Air)
CFL = 0.5        # Courant-Friedrichs-Lewy number (Stability < 1.0)

class CompressibleSolver:
    def __init__(self, nx, length, geom_segments, u0, u_amb, p_init, rho_init, bc_mode="ambient"):
        self.nx = nx
        self.dx = length / nx
        self.t = 0.0
        self.u_amb = u_amb
        self.bc_mode = bc_mode
        
        # Grid generation (Cell centers)
        self.x = np.linspace(0.5*self.dx, length - 0.5*self.dx, nx)
        
        # Geometry generation
        self.A, self.A_faces, self.dA_dx_term, self.vol = self._generate_geometry(length, geom_segments)

        # State Vector U: [rho*A, rho*u*A, E*A]
        # We solve for Conserved Variables per unit length
        self.U = np.zeros((3, nx))
        
        # Initialization
        e_init = p_init / ((GAMMA - 1) * rho_init) # Internal energy
        E_init = rho_init * (e_init + 0.5 * u0**2) # Total energy density
        
        self.U[0, :] = rho_init * self.A
        self.U[1, :] = rho_init * u0 * self.A
        self.U[2, :] = E_init * self.A

        # Ambient conditions (for boundary conditions)
        self.p_amb = p_init
        self.rho_amb = rho_init
        self.T_amb = p_init / (R_GAS * rho_init)

    @staticmethod
    def _fill_slopes(xs, ys, ms):
        ms = np.array(ms, dtype=float)
        for i in range(len(ms)):
            if np.isnan(ms[i]):
                if i == 0: ms[i] = (ys[1] - ys[0]) / (xs[1] - xs[0])
                elif i == len(ms) - 1: ms[i] = (ys[-1] - ys[-2]) / (xs[-1] - xs[-2])
                else: ms[i] = (ys[i + 1] - ys[i - 1]) / (xs[i + 1] - xs[i - 1])
        return ms

    @staticmethod
    def _hermite_eval(xq, xs, ys, ms):
        xs = np.asarray(xs); ys = np.asarray(ys); ms = np.asarray(ms)
        idx = np.searchsorted(xs, xq, side='right') - 1
        idx = np.clip(idx, 0, len(xs) - 2)
        x0 = xs[idx]; x1 = xs[idx + 1]
        y0 = ys[idx]; y1 = ys[idx + 1]
        m0 = ms[idx]; m1 = ms[idx + 1]
        dx = x1 - x0
        t = (xq - x0) / dx
        t = np.clip(t, 0.0, 1.0)
        h00 = 2 * t**3 - 3 * t**2 + 1
        h10 = t**3 - 2 * t**2 + t
        h01 = -2 * t**3 + 3 * t**2
        h11 = t**3 - t**2
        y = h00 * y0 + h10 * dx * m0 + h01 * y1 + h11 * dx * m1
        return y

    def _generate_geometry(self, length, segments):
        xs = [s[0] for s in segments]
        ds = [s[1] for s in segments]
        ms = [s[2] if len(s) > 2 else np.nan for s in segments]
        ms = self._fill_slopes(xs, ds, ms)

        d_center = self._hermite_eval(self.x, xs, ds, ms)
        x_edges = np.linspace(0.0, length, self.nx + 1)
        d_edge = self._hermite_eval(x_edges, xs, ds, ms)

        A_center = np.pi * 0.25 * d_center**2
        A_faces = np.pi * 0.25 * d_edge**2
        
        # Precompute source term geometry factor
        # dA/dx term for the source must match the flux difference discretization exactly
        dA_dx_term = (A_faces[1:] - A_faces[:-1]) / self.dx
        
        return A_center, A_faces, dA_dx_term, A_center * self.dx

    def get_primitive(self):
        """Decode U back to physical variables: rho, u, p, T."""
        rho = self.U[0] / self.A
        u   = self.U[1] / self.U[0]
        E_vol = self.U[2] / self.A
        p = (GAMMA - 1) * (E_vol - 0.5 * rho * u**2)
        T = p / (R_GAS * rho)
        return rho, u, p, T

    def compute_flux_rusanov_well_balanced(self, U_L, U_R, A_L, A_R, A_face):
        """
        Well-Balanced Rusanov Flux.
        1. Fluxes are computed using primitives + Interface Area.
        2. Dissipation is computed using Conserved Variables SCALED to Interface Area.
        """
        # 1. Decode Primitives (using Source Cell Areas)
        rho_L = U_L[0] / A_L
        u_L   = U_L[1] / U_L[0]
        p_L   = (GAMMA - 1) * ((U_L[2]/A_L) - 0.5*rho_L*u_L**2)
        E_vol_L = U_L[2] / A_L

        rho_R = U_R[0] / A_R
        u_R   = U_R[1] / U_R[0]
        p_R   = (GAMMA - 1) * ((U_R[2]/A_R) - 0.5*rho_R*u_R**2)
        E_vol_R = U_R[2] / A_R

        # 2. Compute Physical Fluxes at Interface Area
        # If the flow is static, F = [0, p*A_face, 0]
        F_L = np.array([
            rho_L * u_L,
            rho_L * u_L**2 + p_L,
            u_L * (E_vol_L + p_L)
        ]) * A_face
        
        F_R = np.array([
            rho_R * u_R,
            rho_R * u_R**2 + p_R,
            u_R * (E_vol_R + p_R)
        ]) * A_face

        # 3. Compute Wave Speeds
        c_L = np.sqrt(GAMMA * abs(p_L) / (abs(rho_L) + 1e-12))
        c_R = np.sqrt(GAMMA * abs(p_R) / (abs(rho_R) + 1e-12))
        S_max = max(abs(u_L) + c_L, abs(u_R) + c_R)

        # 4. Well-Balanced Dissipation Term
        # We rescale the conservative variables U to what they would be 
        # if the cell area was exactly A_face.
        # This ensures that if rho, u, E are constant, U*_L == U*_R, 
        # and dissipation becomes zero.
        U_star_L = U_L * (A_face / A_L)
        U_star_R = U_R * (A_face / A_R)

        # Rusanov Flux Formula
        F_face = 0.5 * (F_L + F_R) - 0.5 * S_max * (U_star_R - U_star_L)
        
        return F_face, S_max

    def step(self, dt_fixed=None, prop_params=None):
        rho, u, p, T = self.get_primitive()
        c = np.sqrt(GAMMA * abs(p) / (abs(rho) + 1e-12))
        
        max_speed = np.max(np.abs(u) + c)
        if max_speed < 1e-5: max_speed = 1e-5
        dt = CFL * self.dx / max_speed
        if dt_fixed: dt = dt_fixed

        self.t += dt

        # --- 1. Reconstruction (Ghost Cells) ---
        # We need to handle areas for ghost cells correctly to maintain equilibrium.
        # Simple 'edge' padding for A is sufficient for zero-slope boundaries.
        A_padded = np.pad(self.A, (1, 1), 'edge')
        
        # We pad U manually to ensure primitives are consistent
        U_padded = np.zeros((3, self.nx + 2))
        U_padded[:, 1:-1] = self.U
        
        e_amb = self.p_amb / ((GAMMA-1)*self.rho_amb)
        if self.bc_mode == "ambient":
            # Left Ghost: always ambient free-stream entering from left
            U_padded[0, 0] = self.rho_amb * A_padded[0]
            U_padded[1, 0] = self.rho_amb * self.u_amb * A_padded[0]
            U_padded[2, 0] = self.rho_amb * (e_amb + 0.5*self.u_amb**2) * A_padded[0]

            # Right Ghost: ambient free-stream entering from right (reverse sign)
            u_in = -self.u_amb
            U_padded[0, -1] = self.rho_amb * A_padded[-1]
            U_padded[1, -1] = self.rho_amb * u_in * A_padded[-1]
            U_padded[2, -1] = self.rho_amb * (e_amb + 0.5*u_in**2) * A_padded[-1]
        else:
            # Original: inflow uses ambient stagnation, outflow extrapolates / fixes pressure
            if u[0] > 0: # Inlet
                U_padded[0, 0] = self.rho_amb * A_padded[0]
                U_padded[1, 0] = self.rho_amb * self.u_amb * A_padded[0]
                U_padded[2, 0] = self.rho_amb * (e_amb + 0.5*self.u_amb**2) * A_padded[0]
            else: # Exhaust
                U_padded[:, 0] = self.U[:, 0] * (A_padded[0] / self.A[0])

            if u[-1] < 0: # Inlet from right
                u_in = -self.u_amb
                U_padded[0, -1] = self.rho_amb * A_padded[-1]
                U_padded[1, -1] = self.rho_amb * u_in * A_padded[-1]
                U_padded[2, -1] = self.rho_amb * (e_amb + 0.5*u_in**2) * A_padded[-1]
            else: # Outlet: fix ambient pressure, extrapolate rho,u
                rho_g = rho[-1]
                u_g = u[-1]
                p_g = self.p_amb
                e_g = p_g / ((GAMMA-1)*rho_g)
                U_padded[0, -1] = rho_g * A_padded[-1]
                U_padded[1, -1] = rho_g * u_g * A_padded[-1]
                U_padded[2, -1] = rho_g * (e_g + 0.5*u_g**2) * A_padded[-1]


        # --- 2. Flux Calculation ---
        Flux_interface = np.zeros((3, self.nx + 1))
        
        for i in range(self.nx + 1):
            # i=0 is left face of first cell. 
            # In padded arrays, Left State is index i, Right State is index i+1
            U_L = U_padded[:, i]
            U_R = U_padded[:, i+1]
            A_L = A_padded[i]
            A_R = A_padded[i+1]
            A_face = self.A_faces[i]
            
            Flux_interface[:, i], _ = self.compute_flux_rusanov_well_balanced(
                U_L, U_R, A_L, A_R, A_face
            )

        # Net Flux
        dF = Flux_interface[:, :-1] - Flux_interface[:, 1:]

        # --- 3. Source Terms ---
        S = np.zeros_like(self.U)
        
        # 3a. Geometric Source (Pressure acting on walls)
        # S_mom = p * dA/dx
        # To maintain equilibrium, this term must exactly match the flux difference p*(A_R - A_L)
        # self.dA_dx_term is precalculated as (A_face_R - A_face_L)/dx
        # This matches the flux difference formulation.
        S[1, :] = p * self.dA_dx_term

        # 3b. Propeller Source
        if prop_params and prop_params['force'] != 0:
            enable = True
            if prop_params['pulse_mode']:
                cycle_time = self.t % (prop_params['period_on'] + prop_params['period_off'])
                if cycle_time > prop_params['period_on']:
                    enable = False
            
            if enable:
                mask = (self.x > prop_params['x_start']) & (self.x < prop_params['x_end'])
                if np.any(mask):
                    # Force per unit length
                    f_density = prop_params['force'] / (prop_params['x_end'] - prop_params['x_start'])
                    S[1, mask] += f_density
                    S[2, mask] += f_density * u[mask]

        # --- 4. Update ---
        self.U += (dt / self.dx) * dF + dt * S
        
        return dt

def main():
    parser = argparse.ArgumentParser(description="1D Euler Solver (Well-Balanced)")
    parser.add_argument("--cells",       type=int,   default=100,    help="Number of cells")
    parser.add_argument("--length",      type=float, default=2.0,    help="Tube length (m)")
    parser.add_argument("--dt",          type=float, default=None,   help="Fixed time step")
    parser.add_argument("--prop_force",  type=float, default=10.0,    help="Propeller force (N)")
    parser.add_argument("--pulse",       type=int,   default=0,      help="Enable pulse mode")
    parser.add_argument("--period_on",   type=float, default=0.005,  help="Pulse ON time")
    parser.add_argument("--period_off",  type=float, default=0.05,   help="Pulse OFF time")
    parser.add_argument("--u0",          type=float, default=1.0,    help="Initial velocity inside tube at t=0 (m/s)")
    parser.add_argument("--u_amb",       type=float, default=10.0,   help="Ambient/free-stream velocity at boundaries (m/s)")
    parser.add_argument("--verbosity",   type=int,   default=1,      help="Verbosity level for summary output (0 = off)")
    parser.add_argument("--bc_mode",     type=str,   choices=["ambient", "original"], default="ambient",
                        help="Boundary condition mode: 'ambient' = fixed free-stream both ends; 'original' = inlet uses ambient, outlet extrap/ambient pressure")
    parser.add_argument("--autoscale",   type=int,   default=1,      help="Enable dynamic autoscaling of plots (default: fixed ranges)")
    args = parser.parse_args()

    # Geometry: Wide -> Throat -> Wide
    geom = [
        (0.0,               0.1, 0.0),
        (0.3 * args.length, 0.1, 0.0),
        (0.5 * args.length, 0.05, 0.0), # Narrow throat
        (0.7 * args.length, 0.1, 0.0),
        (      args.length, 0.1, 0.0)
    ]

    solver = CompressibleSolver(
        nx=args.cells, 
        length=args.length, 
        geom_segments=geom, 
        u0=args.u0,
        u_amb=args.u_amb,
        p_init=101325.0,
        rho_init=1.225
        , bc_mode=args.bc_mode
    )
    rho0, u0, p0, T0 = solver.get_primitive()
    
    prop_config = {
        'force': args.prop_force,
        'x_start': 0.1 * args.length,
        'x_end': 0.2 * args.length,
        'pulse_mode': args.pulse,
        'period_on': args.period_on,
        'period_off': args.period_off
    }

    plt.ion()
    fig = plt.figure(figsize=(12, 8))
    gs = fig.add_gridspec(3, 2)
    
    ax_geom = fig.add_subplot(gs[0, 0])
    ax_p = fig.add_subplot(gs[1, 0], sharex=ax_geom)
    ax_u = fig.add_subplot(gs[2, 0], sharex=ax_geom)
    ax_rho = fig.add_subplot(gs[1, 1], sharex=ax_geom)
    ax_hist = fig.add_subplot(gs[2, 1])

    # Static Geometry plot
    ax_geom.plot(solver.x, np.sqrt(solver.A/np.pi)*2, 'k-', lw=2)
    ax_geom.fill_between(solver.x, np.sqrt(solver.A/np.pi)*2, color='gray', alpha=0.3)
    ax_geom.set_ylabel("Diameter (m)")
    ax_geom.set_title("Tube Geometry")
    
    # Propeller region marker
    if args.prop_force != 0:
        ax_geom.axvspan(prop_config['x_start'], prop_config['x_end'], color='red', alpha=0.3, label='Propeller')
        ax_geom.legend()

    # Lines for updating
    line_p, = ax_p.plot([], [], 'b-')
    line_u, = ax_u.plot([], [], 'g-')
    line_rho, = ax_rho.plot([], [], 'r-')
    
    # Conservation history
    history_t = []
    history_mass = []
    history_energy = []
    line_mass, = ax_hist.plot([], [], 'k-', label='Total Mass')
    ax_hist2 = ax_hist.twinx()
    line_eng, = ax_hist2.plot([], [], 'r--', label='Total Energy')
    
    # Plot settings
    ax_p.set_ylabel("Pressure (Pa)")
    ax_u.set_ylabel("Velocity (m/s)")
    ax_rho.set_ylabel("Density (kg/m3)")
    ax_u.set_xlabel("Position (m)")
    ax_hist.set_xlabel("Time (s)")
    ax_hist.set_ylabel("Mass (kg)")
    ax_hist2.set_ylabel("Energy (J)")

    # Reference lines
    ax_p.axhline(solver.p_amb, color='gray', linestyle=':', linewidth=1, label='Ambient p')
    ax_rho.axhline(solver.rho_amb, color='gray', linestyle=':', linewidth=1, label='Ambient rho')
    ax_u.axhline(0.0, color='gray', linestyle=':', linewidth=1, label='u=0')
    
    # Legend for history
    lines = [line_mass, line_eng]
    labels = [l.get_label() for l in lines]
    ax_hist.legend(lines, labels, loc='center right')

    def fixed_limits():
        def padded_range(arr, frac=0.05, min_span=1e-3):
            a_min = float(np.min(arr))
            a_max = float(np.max(arr))
            span = a_max - a_min
            pad = max(frac * max(abs(a_min), abs(a_max), 1.0), min_span)
            if span < min_span:
                a_min -= pad
                a_max += pad
            else:
                a_min -= pad
                a_max += pad
            return a_min, a_max
        p_lim = padded_range(p0, frac=0.01, min_span=100.0)
        rho_lim = padded_range(rho0, frac=0.01, min_span=0.1)
        u_span = max(np.max(np.abs(u0)), 0.5)
        u_lim = (-1.2 * u_span, 1.2 * u_span)
        ax_p.set_ylim(*p_lim)
        ax_rho.set_ylim(*rho_lim)
        ax_u.set_ylim(*u_lim)

    if not args.autoscale:
        fixed_limits()

    def print_profiles(tag):
        rho_f, u_f, p_f, T_f = solver.get_primitive()
        d_f = 2.0 * np.sqrt(solver.A / np.pi)
        data = np.column_stack([solver.x, d_f, solver.A, rho_f, u_f, p_f, T_f])
        header = "x[m]   D[m]   A[m2]   rho[kg/m3]   u[m/s]   p[Pa]   T[K]"
        print(f"\nProfiles ({tag}):")
        print(header)
        np.set_printoptions(linewidth=200)
        print(data)
        np.set_printoptions()

    try:
        frame = 0
        while True:
            dt_sim = solver.step(dt_fixed=args.dt, prop_params=prop_config)
            
            # Monitoring Drift
            if args.u0 == 0 and args.prop_force == 0 and frame % 50 == 0:
                v_max = np.max(np.abs(solver.get_primitive()[1]))
                if v_max > 1e-2:
                    print(f"WARNING: Drift detected! Max Vel = {v_max:.5f} m/s")
            
            if args.verbosity > 0:
                rho_v, u_v, _, _ = solver.get_primitive()
                total_mom = np.sum(solver.U[1, :] * solver.dx)
                total_E = np.sum(solver.U[2, :] * solver.dx)
                total_ke = np.sum(0.5 * rho_v * u_v**2 * solver.A * solver.dx)
                print(f"t={solver.t:.6f}s  mom={total_mom:.6e} kg·m/s  E={total_E:.6e} J  KE={total_ke:.6e} J")
            
            # Visualization (every 5 steps)
            if frame % 5 == 0:
                rho, u, p, T = solver.get_primitive()
                
                line_p.set_data(solver.x, p)
                line_u.set_data(solver.x, u)
                line_rho.set_data(solver.x, rho)
                
                # Auto-scale axes (optional)
                if args.autoscale:
                    ax_p.set_ylim(np.min(p)*0.95, np.max(p)*1.05)
                    ax_u.set_ylim(np.min(u)-10, np.max(u)+10)
                    ax_rho.set_ylim(np.min(rho)*0.95, np.max(rho)*1.05)
                
                # Update history
                history_t.append(solver.t)
                history_mass.append(np.sum(solver.U[0, :] * solver.dx))
                history_energy.append(np.sum(solver.U[2, :] * solver.dx))
                
                # Trim history to keep plot fast
                if len(history_t) > 200:
                    history_t.pop(0)
                    history_mass.pop(0)
                    history_energy.pop(0)

                line_mass.set_data(history_t, history_mass)
                line_eng.set_data(history_t, history_energy)
                
                ax_hist.set_xlim(min(history_t), max(history_t) + dt_sim)
                ax_hist.set_ylim(min(history_mass)*0.999, max(history_mass)*1.001)
                ax_hist2.set_ylim(min(history_energy)*0.95, max(history_energy)*1.05)

                ax_geom.set_title(f"t={solver.t:.4f}s | Propeller: {'ON' if (solver.t % (prop_config['period_on']+prop_config['period_off']) < prop_config['period_on']) and args.pulse else 'OFF'}")

                plt.pause(0.001)
                if not plt.fignum_exists(fig.number):
                    break
            
            frame += 1

    except KeyboardInterrupt:
        print("\nSimulation stopped by user.")
    finally:
        print_profiles("final")

if __name__ == "__main__":
    main()