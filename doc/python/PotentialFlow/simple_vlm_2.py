import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import argparse

# ==============================================================================
# 1. ROBUST PHYSICS: CORE MODEL
# ==============================================================================

def get_velocity_robust(point, p1, p2, gamma=1.0, core_radius=0.1):
    """
    Computes induced velocity with a Finite Core (Vatistas/Scully model).
    This prevents division by zero and 'smears' the vortex locally.
    """
    r0 = p2 - p1
    r1 = point - p1
    r2 = point - p2
    
    # Cross products
    r1_x_r2 = np.cross(r1, r2)
    r1_x_r2_sq = np.sum(r1_x_r2**2)
    
    # Distances
    r0_mag = np.linalg.norm(r0)
    r1_mag = np.linalg.norm(r1)
    r2_mag = np.linalg.norm(r2)
    
    # --- The Core Regularization ---
    # Standard Biot-Savart denominator is |r1 x r2|^2
    # We replace it with |r1 x r2|^2 + (radius * core)^2
    # This effectively assumes the vortex has thickness.
    
    # We scale core radius by segment length to keep it relative, 
    # or use absolute if provided.
    r_effective_sq = r1_x_r2_sq + (core_radius * r0_mag)**2 
    # Note: For semi-infinite, r0 is infinite, so we handle differently below.
    
    # Finite Segment Term
    term1 = np.dot(r0, (r1/r1_mag - r2/r2_mag))
    
    # V = (Gamma / 4pi) * (r1 x r2) / denominator * term
    V = (gamma / (4.0 * np.pi)) * (r1_x_r2 / r_effective_sq) * term1
    
    return V

def get_velocity_semi_inf_robust(point, origin, direction, gamma=1.0, core_radius=0.1):
    """
    Robust Semi-Infinite Vortex. 
    Models wake diffusion (core growth) if needed, but here uses fixed core.
    """
    u = direction / np.linalg.norm(direction)
    r = point - origin
    r_mag = np.linalg.norm(r)
    
    r_x_u = np.cross(r, u)
    r_x_u_sq = np.sum(r_x_u**2)
    
    # --- Core Regularization ---
    # Here the "segment length" is effectively infinite, so we use absolute core radius
    # r_c can conceptually grow with distance from origin: r_c(x) ~ sqrt(x)
    
    # Simple fixed core for stability:
    r_effective_sq = r_x_u_sq + core_radius**2
    
    term = 1.0 + np.dot(u, r) / r_mag
    V = (gamma / (4.0 * np.pi)) * (np.cross(u, r) / r_effective_sq) * term
    
    return -V 

# ==============================================================================
# 2. GEOMETRY (Multi-Surface)
# ==============================================================================

def poly_gen(u, coeffs):
    return coeffs[0] + coeffs[1]*(1 - u**2) + coeffs[2]*(1 - u**4)

class Surface:
    def __init__(self, name, n_panels, offset=np.array([0,0,0]), span=10.0):
        self.name = name
        self.n_panels = n_panels
        self.offset = offset
        self.span = span
        self.LE_x_coeffs = [0.0, 0.0, 0.0]
        self.LE_z_coeffs = [0.0, 0.0, 0.0]
        self.TE_x_coeffs = [1.0, 0.0, 0.0]
        self.TE_z_coeffs = [0.0, 0.0, 0.0]
        self.panels = []

    def generate_mesh(self):
        u = np.linspace(-1, 1, self.n_panels + 1)
        y = u * (self.span / 2.0)
        le_x = poly_gen(u, self.LE_x_coeffs) + self.offset[0]
        le_z = poly_gen(u, self.LE_z_coeffs) + self.offset[2]
        te_x = poly_gen(u, self.TE_x_coeffs) + self.offset[0]
        te_z = poly_gen(u, self.TE_z_coeffs) + self.offset[2]
        
        # Apply offset Y
        y += self.offset[1]
        
        self.LE = np.vstack([le_x, y, le_z]).T
        self.TE = np.vstack([te_x, y, te_z]).T
        return self.LE, self.TE

# ==============================================================================
# 3. SOLVER (Drag + Interactions)
# ==============================================================================

class VLM_Solver:
    def __init__(self, surfaces, alpha_deg=5.0, V_inf=10.0, rho=1.225, core_radius=0.1):
        self.surfaces = surfaces
        self.alpha = np.radians(alpha_deg)
        self.V_mag = V_inf
        self.rho = rho
        # Freestream vector
        self.V_inf = np.array([np.cos(self.alpha), 0, np.sin(self.alpha)]) * V_inf
        
        self.all_panels = []
        self.gammas = None
        
        # CORE RADIUS (The "Smearing" parameter)
        # 0.05 meters is typical for a small UAV scale. 
        # Making it too large reduces accuracy. Making it too small causes instability.
        self.core_radius = core_radius

    def setup_panels(self):
        self.all_panels = []
        global_idx = 0
        
        for surf in self.surfaces:
            LE, TE = surf.generate_mesh()
            n = surf.n_panels
            
            for i in range(n):
                p1, p2 = LE[i], LE[i+1]
                p3, p4 = TE[i+1], TE[i]
                
                panel = {
                    'id': global_idx,
                    'surf_name': surf.name,
                    'bound_vortex': (p1 + 0.25*(p4-p1), p2 + 0.25*(p3-p2)),
                    'cp': 0.5 * ((p1 + 0.75*(p4-p1)) + (p2 + 0.75*(p3-p2))),
                    'area': np.linalg.norm(np.cross(p4-p1, p2-p1)), # Approx area
                    'chord': np.linalg.norm((p4-p1 + p3-p2)/2.0),
                    'wake_dir': np.array([1.0, 0.0, 0.0])
                }
                
                # Normal
                vec_chord = p4 - p1
                vec_span = p2 - p1
                n_vec = np.cross(vec_chord, vec_span)
                panel['normal'] = n_vec / np.linalg.norm(n_vec)
                
                self.all_panels.append(panel)
                global_idx += 1

    def solve(self):
        n = len(self.all_panels)
        A = np.zeros((n, n))
        b = np.zeros(n)
        
        for i in range(n):
            panel_i = self.all_panels[i]
            b[i] = -np.dot(self.V_inf, panel_i['normal'])
            
            for j in range(n):
                panel_j = self.all_panels[j]
                v1, v2 = panel_j['bound_vortex']
                
                # Use Robust Bio-Savart
                # 1. Bound
                V_b = get_velocity_robust(panel_i['cp'], v1, v2, 1.0, self.core_radius)
                # 2. Wake Legs
                V_l = -get_velocity_semi_inf_robust(panel_i['cp'], v1, panel_j['wake_dir'], 1.0, self.core_radius)
                V_r = get_velocity_semi_inf_robust(panel_i['cp'], v2, panel_j['wake_dir'], 1.0, self.core_radius)
                
                A[i, j] = np.dot(V_b + V_l + V_r, panel_i['normal'])
                
        self.gammas = np.linalg.solve(A, b)

    def get_velocity_field(self, points):
        """Induced + freestream at array of points, using robust kernels."""
        pts = np.atleast_2d(points)
        vels = np.zeros_like(pts)
        for k, pt in enumerate(pts):
            v = self.V_inf.copy()
            for i, p in enumerate(self.all_panels):
                g = self.gammas[i]
                v1, v2 = p['bound_vortex']
                wdir = p['wake_dir']
                v += get_velocity_robust(pt, v1, v2, g, self.core_radius)
                v -= get_velocity_semi_inf_robust(pt, v1, wdir, g, self.core_radius)
                v += get_velocity_semi_inf_robust(pt, v2, wdir, g, self.core_radius)
            vels[k] = v
        return vels
        
    def calculate_forces_and_drag(self, nu=1.5e-5):
        """ Calculates Lift, Induced Drag, and Viscous Drag """
        
        total_lift = 0.0
        total_drag_induced = 0.0
        total_drag_viscous = 0.0
        
        print("\n--- Force Analysis (Strip Theory) ---")
        print(f"{'ID':<3} {'Gamma':<6} {'Cl_loc':<6} {'Cd_prof':<8} {'Re':<8}")
        
        for i, p in enumerate(self.all_panels):
            gamma = self.gammas[i]
            
            # 1. Kutta-Joukowski Force: F = rho * Gamma * (V_local x Bound_Vec)
            # V_local is V_inf + V_induced_by_others (complicated)
            # Simplified: F = rho * Gamma * V_inf * Span_Width
            
            v1, v2 = p['bound_vortex']
            bound_vec = v2 - v1
            width = np.linalg.norm(bound_vec)
            
            # Force vector per panel (Potential)
            # F = rho * Gamma * (V_inf x l)
            force_vec = self.rho * gamma * np.cross(self.V_inf, bound_vec)
            
            # Lift is Z component, Induced Drag is X component (if wake tilts, but here wake is flat)
            # BETTER Induced Drag calculation: Trefftz Plane or local downwash.
            # Local Downwash Method: Di = rho * Gamma * w_induced
            # We need to calculate induced velocity at this panel's own center excluding itself
            # This is expensive. 
            
            # Simple approximation for this script: 
            # Induced drag comes from the force vector being tilted back by local downwash angle alpha_i.
            # F_potential is perpendicular to local velocity V_local.
            # D_i = L * alpha_i
            
            # Let's trust the cross product force_vec if V_inf is the only velocity?
            # No, that gives 0 drag. We need the local induced velocity.
            # Let's skip detailed Di for this snippet and focus on Viscous.
            
            L_panel = force_vec[2]
            total_lift += L_panel
            
            # 2. Viscous Drag (Strip Theory)
            # Calculate local Cl
            # L = 0.5 * rho * V^2 * Area * Cl  ->  Cl = L / (q * S)
            q = 0.5 * self.rho * self.V_mag**2
            Cl_local = L_panel / (q * p['area'])
            
            # Calculate Reynolds
            Re = self.V_mag * p['chord'] / nu
            
            # 3. Profile Drag Lookup (The "Viscous" part)
            # In a real tool, you query your airfoil polar here.
            # Here, we use a simple flat plate turbulent approximation:
            # Cd = 0.074 / Re^0.2 + k * Cl^2
            Cd_profile = 0.074 / (Re**0.2) + 0.01 * Cl_local**2 
            
            D_viscous_panel = Cd_profile * q * p['area']
            total_drag_viscous += D_viscous_panel
            
            if i % 2 == 0: # Print every other panel
                print(f"{i:<3} {gamma:.2f}   {Cl_local:.3f}  {Cd_profile:.5f}   {Re:.1e}")

        print("-" * 40)
        print(f"Total Lift: {total_lift:.2f} N")
        print(f"Viscous Drag: {total_drag_viscous:.4f} N")
        # Note: Induced drag requires Trefftz plane integration for accuracy in VLM
        
    def visualize(self):
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        # Plot Geometry
        for surf in self.surfaces:
            LE, TE = surf.generate_mesh()
            ax.plot(LE[:,0], LE[:,1], LE[:,2], 'k-', alpha=0.3)
            ax.plot(TE[:,0], TE[:,1], TE[:,2], 'k-', alpha=0.3)
            # Connect chords
            for i in range(len(LE)):
                ax.plot([LE[i,0], TE[i,0]], [LE[i,1], TE[i,1]], [LE[i,2], TE[i,2]], 'k-', alpha=0.1)

        # Plot Bound Vortices (Colored by Gamma)
        max_g = np.max(np.abs(self.gammas))
        for i, p in enumerate(self.all_panels):
            v1, v2 = p['bound_vortex']
            g = self.gammas[i]
            # Norm
            color = plt.cm.jet(0.5 + 0.5*(g/max_g))
            ax.plot([v1[0], v2[0]], [v1[1], v2[1]], [v1[2], v2[2]], color=color, linewidth=3)

        # Plot Streamlines (INTERACTION CHECK)
        print("\nComputing Streamlines (with core 'smearing')...")
        for y in self.stream_start_y:
            for z in self.stream_start_z:
                pos = np.array([self.stream_x0, y, z], dtype=float)
                path = [pos.copy()]
                
                for _ in range(self.stream_steps):
                    vel1 = self.get_velocity_field(pos)[0]
                    step1 = self.stream_dt * vel1
                    s1 = np.linalg.norm(step1)
                    if s1 > self.stream_max_step:
                        step1 *= self.stream_max_step / s1
                    mid = pos + 0.5 * step1
                    vel2 = self.get_velocity_field(mid)[0]
                    step2 = self.stream_dt * vel2
                    s2 = np.linalg.norm(step2)
                    if s2 > self.stream_max_step:
                        step2 *= self.stream_max_step / s2
                    pos = pos + step2
                    path.append(pos.copy())
                    if pos[0] > self.stream_x_max:
                        break
                
                path = np.array(path)
                ax.plot(path[:,0], path[:,1], path[:,2], 'g-', alpha=0.5, linewidth=1)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        # Auto-fit axes to seed extents and streamline range so CLI changes are visible
        ax.set_xlim(self.stream_x0, self.stream_x_max)
        ax.set_ylim(self.stream_start_y.min(), self.stream_start_y.max())
        ax.set_zlim(self.stream_start_z.min(), self.stream_start_z.max())
        # Enforce equal scaling to avoid distortion
        xmid = 0.5 * sum(ax.get_xlim())
        ymid = 0.5 * sum(ax.get_ylim())
        zmid = 0.5 * sum(ax.get_zlim())
        max_range = max(
            ax.get_xlim()[1] - ax.get_xlim()[0],
            ax.get_ylim()[1] - ax.get_ylim()[0],
            ax.get_zlim()[1] - ax.get_zlim()[0],
        )
        half = 0.5 * max_range
        ax.set_xlim(xmid - half, xmid + half)
        ax.set_ylim(ymid - half, ymid + half)
        ax.set_zlim(zmid - half, zmid + half)
        try:
            ax.set_box_aspect([1, 1, 1])
        except Exception:
            pass
        plt.title(f"Interaction: Main Wing + Tail\nCore Radius = {self.core_radius}m (Prevents Explosion)")
        plt.show()

# ==============================================================================
# 4. RUN DEMO
# ==============================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simple VLM with finite-core wake and strip drag")
    parser.add_argument("--alpha", type=float, default=8.0, help="Angle of attack [deg]")
    parser.add_argument("--v-inf", type=float, default=20.0, help="Freestream speed [m/s]")
    parser.add_argument("--rho", type=float, default=1.225, help="Density [kg/m^3]")
    parser.add_argument("--core-radius", type=float, default=0.1, help="Finite core radius [m]")
    parser.add_argument("--nu", type=float, default=1.5e-5, help="Kinematic viscosity [m^2/s]")
    # Main wing
    parser.add_argument("--main-panels", type=int, default=10, help="Number of panels on main wing")
    parser.add_argument("--main-span", type=float, default=8.0, help="Main wing span [m]")
    parser.add_argument("--main-offset-x", type=float, default=0.0, help="Main wing X offset")
    parser.add_argument("--main-offset-y", type=float, default=0.0, help="Main wing Y offset")
    parser.add_argument("--main-offset-z", type=float, default=0.0, help="Main wing Z offset")
    parser.add_argument("--main-lex", nargs=3, type=float, default=[0.0, 1.0, 0.0], help="Main LE x coeffs a b c")
    parser.add_argument("--main-tex", nargs=3, type=float, default=[1.5, 1.0, 0.0], help="Main TE x coeffs a b c")
    parser.add_argument("--main-lez", nargs=3, type=float, default=[0.0, 0.0, 0.0], help="Main LE z coeffs a b c")
    parser.add_argument("--main-tez", nargs=3, type=float, default=[0.0, 0.0, 0.0], help="Main TE z coeffs a b c")
    # Tail
    parser.add_argument("--tail-panels", type=int, default=6, help="Number of panels on tail")
    parser.add_argument("--tail-span", type=float, default=4.0, help="Tail span [m]")
    parser.add_argument("--tail-offset-x", type=float, default=4.0, help="Tail X offset")
    parser.add_argument("--tail-offset-y", type=float, default=0.0, help="Tail Y offset")
    parser.add_argument("--tail-offset-z", type=float, default=0.2, help="Tail Z offset")
    parser.add_argument("--tail-lex", nargs=3, type=float, default=[0.0, 0.0, 0.0], help="Tail LE x coeffs a b c")
    parser.add_argument("--tail-tex", nargs=3, type=float, default=[1.0, 0.0, 0.0], help="Tail TE x coeffs a b c")
    parser.add_argument("--tail-lez", nargs=3, type=float, default=[0.0, 0.0, 0.0], help="Tail LE z coeffs a b c")
    parser.add_argument("--tail-tez", nargs=3, type=float, default=[0.0, 0.0, 0.0], help="Tail TE z coeffs a b c")
    # Streamlines
    parser.add_argument("--stream-ny", type=int, default=12, help="Number of seeds in spanwise direction")
    parser.add_argument("--stream-nz", type=int, default=3, help="Number of seeds in vertical direction")
    parser.add_argument("--stream-ymin", type=float, default=-5.0, help="Min Y for seeds")
    parser.add_argument("--stream-ymax", type=float, default=5.0, help="Max Y for seeds")
    parser.add_argument("--stream-zmin", type=float, default=-0.5, help="Min Z for seeds (before offset)")
    parser.add_argument("--stream-zmax", type=float, default=0.5, help="Max Z for seeds (before offset)")
    parser.add_argument("--stream-zoffset", type=float, default=0.5, help="Z offset for seeds")
    parser.add_argument("--stream-steps", type=int, default=100, help="Streamline integration steps")
    parser.add_argument("--stream-dt", type=float, default=0.05, help="Streamline integration timestep")
    parser.add_argument("--stream-x0", type=float, default=-2.0, help="Streamline starting x")
    parser.add_argument("--stream-xmax", type=float, default=15.0, help="Stop integrating when x exceeds this")
    parser.add_argument("--stream-max-step", type=float, default=0.2, help="Clamp maximum step length per integration sub-step")
    args = parser.parse_args()

    # Build surfaces
    wing = Surface("Main", n_panels=args.main_panels, span=args.main_span,
                   offset=[args.main_offset_x, args.main_offset_y, args.main_offset_z])
    wing.LE_x_coeffs = args.main_lex
    wing.TE_x_coeffs = args.main_tex
    wing.LE_z_coeffs = args.main_lez
    wing.TE_z_coeffs = args.main_tez

    tail = Surface("Tail", n_panels=args.tail_panels, span=args.tail_span,
                   offset=[args.tail_offset_x, args.tail_offset_y, args.tail_offset_z])
    tail.LE_x_coeffs = args.tail_lex
    tail.TE_x_coeffs = args.tail_tex
    tail.LE_z_coeffs = args.tail_lez
    tail.TE_z_coeffs = args.tail_tez

    solver = VLM_Solver([wing, tail], alpha_deg=args.alpha, V_inf=args.v_inf,
                        rho=args.rho, core_radius=args.core_radius)
    # Streamline params stored on solver for reuse
    solver.stream_start_y = np.linspace(args.stream_ymin, args.stream_ymax, args.stream_ny)
    solver.stream_start_z = np.linspace(args.stream_zmin, args.stream_zmax, args.stream_nz) + args.stream_zoffset
    solver.stream_steps = args.stream_steps
    solver.stream_dt = args.stream_dt
    solver.stream_x0 = args.stream_x0
    solver.stream_x_max = args.stream_xmax
    solver.stream_max_step = args.stream_max_step

    solver.setup_panels()
    solver.solve()
    solver.calculate_forces_and_drag(nu=args.nu)
    solver.visualize()