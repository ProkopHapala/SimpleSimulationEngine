import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate  # For trapezoid
import argparse

from SlenderBodyTheory2 import AxialSolver

# ==============================================================================
# 1. ADVANCED BOUNDARY LAYER SOLVER (with better separation & drag)
# ==============================================================================

class AdvancedBoundaryLayer:
    def __init__(self, nu=1.0e-6, rho=1025.0):
        self.nu = nu
        self.rho = rho
        self.viscous_force = 0.0
        self.pressure_force = 0.0
        self.total_force_axial = 0.0
        self.total_force_crossflow = 0.0
        self.Cd_profile_final = 0.0
        self.Cd_crossflow_final = 0.0

    def solve(self, s, Ue, V_inf, alpha_deg, L_sub, R_sub, Ref_Area):
        n = len(s)
        theta = np.zeros(n)
        H = np.zeros(n)
        Cf = np.zeros(n)
        
        dUds = np.gradient(Ue, s)
        
        is_turbulent = False
        self.separated = False
        self.sep_index = -1
        
        # Laminar Init
        if dUds[0] > 1e-6:
            theta[0] = np.sqrt(0.075 * self.nu / dUds[0])
        else:
            theta[0] = 0.0
        H[0] = 2.6

        for i in range(n - 1):
            if self.separated:
                # Wake Model: Momentum Thickness stays constant
                theta[i+1] = theta[i]
                H[i+1] = H[i]
                Cf[i+1] = 0.0
                continue

            ds = s[i+1] - s[i]
            U = Ue[i]
            dU = dUds[i]
            
            # --- LAMINAR STEP ---
            if not is_turbulent:
                term1 = theta[i]**2 * (U**6)
                term2 = 0.45 * self.nu * (U**5) * ds
                theta_sq = (term1 + term2) / (Ue[i+1]**6 + 1e-12)
                theta[i+1] = np.sqrt(theta_sq)
                
                lam = (theta[i+1]**2 / self.nu) * dUds[i+1]
                lam = max(-0.09, min(0.25, lam))
                
                if lam >= 0:
                    H[i+1] = 2.61 - 3.75*lam + 5.24*lam**2
                else:
                    H[i+1] = 2.088 + 0.0731/(lam + 0.14)
                
                l = 0.22 + 1.57*lam - 1.8*lam**2 if lam >= 0 else 0.22 + 1.402*lam + (0.018*lam)/(lam+0.107)
                Cf[i+1] = (2 * self.nu * l) / (Ue[i+1] * theta[i+1] + 1e-12)

                Re_th = U * theta[i+1] / self.nu
                if Re_th > 400 or lam < -0.08:
                    is_turbulent = True
                    H[i+1] = 1.4
            
            # --- TURBULENT STEP (Head's Method) ---
            else:
                th = theta[i]
                Sh = H[i]
                Re_th = U * th / self.nu
                if Re_th < 100: Re_th = 100
                Cf_turb = 0.246 * (10**(-0.678*Sh)) * (Re_th**-0.268)
                Cf[i+1] = Cf_turb
                
                if Sh < 1.1: Sh = 1.1
                H1 = 1.535 * (Sh - 0.7)**(-2.715) + 3.3
                Ce = 0.0306 * (H1 - 3.0)**(-0.6169)
                
                dTh_ds = (Cf_turb / 2.0) - (Sh + 2.0) * (th / (U+1e-9)) * dU
                theta[i+1] = th + dTh_ds * ds
                
                # Simplified H evolution (assume H ~ 1.4 for turbulent)
                # A full solver would integrate H1 based on Ce and other terms
                H[i+1] = 1.4 

                # --- SEPARATION CHECK ---
                if H[i+1] > 2.5:
                    self.separated = True
                    self.sep_index = i+1

        self.s = s
        self.theta = theta
        self.H = H
        self.Cf = Cf
        
        # --- Drag Calculation ---
        # Friction Drag Force (integrate Cf over wetted area)
        # dA = perimeter * ds. Perimeter ~ 2*pi*radius(s)
        perimeter = 2 * np.pi * get_radius_func(L_sub, R_sub)(s) # Assuming circular cross-section
        wetted_area = scipy.integrate.trapezoid(perimeter, s)
        
        Cf_vals = Cf
        if self.separated:
            # Integrate Cf only up to separation
            Cf_vals = Cf[:self.sep_index]
            s_active = s[:self.sep_index]
            perimeter_active = perimeter[:self.sep_index]
        else:
            s_active = s
            perimeter_active = perimeter
            
        friction_force = scipy.integrate.trapezoid(Cf_vals * perimeter_active, s_active) * (0.5 * self.rho * V_inf**2)
        self.viscous_force = friction_force
        
        # Pressure Drag (from wake momentum deficit)
        # Based on Squire-Young at separation or TE
        theta_wake = theta[self.sep_index] if self.separated else theta[-1]
        H_wake = H[self.sep_index] if self.separated else H[-1]
        
        # Velocity ratio in wake (approximate)
        U_wake_ratio = 0.95 # Can be refined
        
        # Wake momentum thickness scaled by reference Length
        # This gives a Cd_profile value first
        # Cd_profile_temp = 2 * theta_wake * (U_wake_ratio)**((H_wake + 5)/2) / L_sub # Length based
        
        # Force = Integral( Pressure Difference * Area element )
        # DeltaP ~ rho * U_inf^2 * theta_wake * (H_wake+2)/U_inf^2 
        # Force ~ rho * U_inf^2 * theta_wake * (H_wake+2) * (Projected Area) ??
        # Let's use a simpler Base Drag model for separation:
        # Pressure_drag = (Pressure_behind - Pressure_inf) * Area_base
        # Simplified: Pressure_behind ~ stagnation pressure
        # DeltaP = 0.5 * rho * U_inf^2 (assuming potential flow recovery at front)
        # Base Area ~ pi * R_sep^2
        
        # Simplified Wake drag: Use momentum thickness
        # Force ~ Wake Area Momentum deficit
        # Area of wake ~ pi * R_wake^2. Momentum deficit ~ U_inf * Theta_wake
        # Force ~ 0.5 * rho * U_inf * Theta_wake * (2*pi*R_wake) * U_inf ?
        # This is complex. A simpler approximation for slender bodies is:
        # Pressure Drag ~ Friction Drag * (H-1)
        # OR Base Drag ~ Cd_base * 0.5 * rho * U_inf^2 * (pi*R_sep^2)
        # Let's use a relation to Cd_profile
        
        # Re-estimating Cd_profile using area
        # Drag_profile = Cd_profile_temp * Ref_Area (if using temp Cd)
        # Let's use a simpler force approach:
        # Force_pressure = Integral ( -dP * Area )
        # dP ~ -rho * Ue * dUe -> Integral dP ~ -rho * Integral(Ue dUe) = -0.5*rho*Ue^2
        # The pressure recovery is NOT to U_inf at the tail, it's to U_wake.
        # Force_pressure = 0.5 * rho * U_inf^2 * ( (U_wake/U_inf)^2 - 1 ) * Area_reference
        
        # A pragmatic approach: Base Drag from turbulent separation
        # Cd_base ~ 0.1 to 0.5 for bluff bodies, much lower for streamlined
        # Let's use the theta/H relationship directly for force
        # Effective 'wake displacement thickness' ~ theta * (H+1)/2
        wake_disp_thickness = theta_wake * (H_wake + 1) / 2.0
        
        # Force = Pressure Diff * Area ~ (0.5 * rho * U_inf^2) * (wake_disp_thickness / L_sub) * (Projection Area)
        # Approx Project Area: L * 2*R_sub
        # Force_pressure = 0.5 * rho * U_inf**2 * (wake_disp_thickness / L_sub) * (L_sub * 2*R_sub)
        # Force_pressure = rho * U_inf**2 * wake_disp_thickness * R_sub
        
        # Let's use the Squire-Young force relation: F ~ rho * U^2 * Theta * (H+1)/2
        # This is a dimensionally correct way to estimate the momentum deficit force
        wake_momentum_deficit_area = theta_wake * (H_wake + 1.0) / 2.0
        pressure_force = (0.5 * self.rho * V_inf**2) * wake_momentum_deficit_area * (2 * np.pi * R_sub) # Scaling force by perimeter
        self.pressure_force = pressure_force
        
        # Total Axial Force
        self.total_force_axial = self.viscous_force + self.pressure_force
        
        # Final Cd based on Reference Area
        # For slender bodies, Reference Area often is L*D or just L^2
        # Let's use Lateral Area = L_sub * 2*R_sub
        self.Cd_profile_final = self.total_force_axial / (0.5 * self.rho * V_inf**2 * (L_sub * 2*R_sub))
        
        return theta, H, Cf

# ==============================================================================
# 2. AXIAL SOLVER (With Pressure Field Visualization)
# ==============================================================================

# --- Assume AxialSolver class is defined here ---
# (Including get_point_source_vel, get_point_dipole_vel, solve, geometry_ellipsoid, etc.)

class AxialSolver:
    def __init__(self, length=2.0, radius=0.2, n_basis=20):
        self.L = length
        self.R = radius
        self.n_basis = n_basis
        
        t = np.linspace(0, np.pi, n_basis)
        self.basis_x = (length / 2.0) * (1 - np.cos(t)) - (length/2.0)
        self.basis_locs = np.zeros((n_basis, 3))
        self.basis_locs[:, 0] = self.basis_x + length/2.0
        
        self.sigma_weights = np.zeros(n_basis)
        self.mu_weights = np.zeros(n_basis)
        self.surface_pts = None
        self.normals = None

    def get_point_source_vel(self, field_points, source_loc):
        r_vec = field_points - source_loc
        r_sq = np.sum(r_vec**2, axis=1)
        r_mag = np.sqrt(r_sq)
        r_mag3 = r_mag**3 + 1e-9
        return (1.0 / (4 * np.pi)) * r_vec / r_mag3[:, np.newaxis]

    def get_point_dipole_vel(self, field_points, dipole_loc, axis_dir=np.array([0,0,1])):
        r_vec = field_points - dipole_loc
        r_sq = np.sum(r_vec**2, axis=1)
        r_mag = np.sqrt(r_sq)
        r_5 = r_mag**5 + 1e-9
        r_3 = r_mag**3 + 1e-9
        m_dot_r = np.dot(r_vec, axis_dir)
        term1 = 3 * m_dot_r[:, np.newaxis] * r_vec / r_5[:, np.newaxis]
        term2 = axis_dir / r_3[:, np.newaxis]
        return (1.0 / (4 * np.pi)) * (term1 - term2)

    def geometry_ellipsoid(self, n_points=100):
        theta = np.linspace(0, 2*np.pi, n_points)
        x = (self.L/2.0) * np.cos(theta) + self.L/2.0
        z = self.R * np.sin(theta)
        y = np.zeros_like(x)
        points = np.vstack([x, y, z]).T
        
        nx = 2*(x - self.L/2.0) / (self.L/2.0)**2
        nz = 2*z / self.R**2
        ny = np.zeros_like(nx)
        normals = np.vstack([nx, ny, nz]).T
        norm = np.linalg.norm(normals, axis=1, keepdims=True)
        normals = normals / (norm + 1e-9)
        return points, normals

    def solve(self, V_inf_mag, alpha_deg):
        alpha = np.radians(alpha_deg)
        V_freestream = V_inf_mag * np.array([np.cos(alpha), 0, np.sin(alpha)])
        self.V_freestream = V_freestream
        
        pts, normals = self.geometry_ellipsoid(n_points=200)
        self.surface_pts = pts
        self.normals = normals
        
        n_pts = len(pts)
        n_dof = 2 * self.n_basis
        A = np.zeros((n_pts, n_dof))
        
        for j in range(self.n_basis):
            loc = self.basis_locs[j]
            V_src = self.get_point_source_vel(pts, loc)
            A[:, j] = np.sum(V_src * normals, axis=1)
            V_dip = self.get_point_dipole_vel(pts, loc, axis_dir=np.array([0,0,1]))
            A[:, self.n_basis + j] = np.sum(V_dip * normals, axis=1)

        b = -np.sum(V_freestream * normals, axis=1)
        x, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)
        
        self.sigma_weights = x[0:self.n_basis]
        self.mu_weights    = x[self.n_basis:]

    def get_total_velocity(self, grid_pts):
        V_tot = np.zeros_like(grid_pts) + self.V_freestream
        for j in range(self.n_basis):
            loc = self.basis_locs[j]
            V_tot += self.sigma_weights[j] * self.get_point_source_vel(grid_pts, loc)
            V_tot += self.mu_weights[j] * self.get_point_dipole_vel(grid_pts, loc, np.array([0,0,1]))
        return V_tot

# Helper function for get_radius (used by BL solver)
def get_radius_func(L_sub, R_sub):
    def radius_at_x(x_val):
        a = L_sub / 2.0
        # Avoid sqrt of negative numbers from x_val outside [0, L_sub]
        term = 1 - ((x_val - a) / a)**2
        term = np.maximum(term, 0) # Clamp to zero
        return R_sub * np.sqrt(term)
    return radius_at_x

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(description="SlenderBodyTheory5-a potential flow + BL solver")
    parser.add_argument("--alpha",   type=float, default=15.0, help="Angle of attack in degrees")
    parser.add_argument("--n-basis", type=int,   default=8,    help="Number of axial basis points (sources+dipoles)")
    parser.add_argument("--n-arc",   type=int,   default=100,  help="Number of surface samples along meridian for BL")
    parser.add_argument("--nx",      type=int,   default=80,   help="Grid resolution in x for plots")
    parser.add_argument("--nz",      type=int,   default=80,   help="Grid resolution in z for plots")
    parser.add_argument("--cp-vmin", type=float, default=None, help="Explicit Cp vmin for plotting (default: min of finite Cp)")
    parser.add_argument("--cp-vmax", type=float, default=0.2,  help="Explicit Cp vmax for plotting (default: max of finite Cp)")
    parser.add_argument("--length",  type=float, default=5.0,  help="Body length [m]")
    parser.add_argument("--radius",  type=float, default=0.5,  help="Max radius [m]")
    parser.add_argument("--v-inf",   type=float, default=10.0, help="Freestream speed [m/s]")
    parser.add_argument("--nu",      type=float, default=1e-6, help="Kinematic viscosity")
    parser.add_argument("--rho",     type=float, default=1025.0, help="Density")
    args = parser.parse_args()

    L_sub = args.length
    R_sub = args.radius
    V_inf = args.v_inf
    alpha = args.alpha
    nu = args.nu
    rho = args.rho

    Ref_Area_Lateral = L_sub * (2 * R_sub)

    solver = AxialSolver(length=L_sub, radius=R_sub, n_basis=args.n_basis)
    solver.solve(V_inf_mag=V_inf, alpha_deg=alpha)

    n_pts_arc = args.n_arc
    x_arc = np.linspace(0.001, L_sub - 0.001, n_pts_arc)
    radius_at_x = get_radius_func(L_sub, R_sub)
    z_arc = radius_at_x(x_arc)
    query_pts_arc = np.vstack([x_arc, np.zeros(n_pts_arc), z_arc]).T

    V_vec_arc = solver.get_total_velocity(query_pts_arc)
    Ue = np.linalg.norm(V_vec_arc, axis=1)

    deltas = np.diff(query_pts_arc, axis=0)
    s = np.concatenate(([0], np.cumsum(np.linalg.norm(deltas, axis=1))))

    bl = AdvancedBoundaryLayer(nu=nu, rho=rho)
    theta, H, Cf = bl.solve(s, Ue, V_inf, alpha, L_sub, R_sub, Ref_Area_Lateral)

    F_friction = bl.viscous_force
    F_pressure = bl.pressure_force
    F_drag_axial = bl.total_force_axial

    alpha_rad = np.radians(alpha)
    V_cross_flow_avg = V_inf * np.sin(alpha_rad)
    x_geom = np.linspace(0, L_sub, 50)
    radii_geom = radius_at_x(x_geom)
    diameter_geom = 2 * radii_geom
    projected_area_side = scipy.integrate.trapezoid(diameter_geom, x_geom)

    Cd_cyl_turb = 0.7
    F_cross_force = Cd_cyl_turb * (0.5 * rho * V_cross_flow_avg**2) * projected_area_side
    F_drag_total_opposing = F_drag_axial + F_cross_force * np.sin(alpha_rad)

    # Cd based on lateral area
    Cd_axial_final = F_drag_axial / (0.5 * rho * V_inf**2 * Ref_Area_Lateral)
    Cd_crossflow_final = F_cross_force / (0.5 * rho * V_inf**2 * Ref_Area_Lateral)
    Cd_total_final = Cd_axial_final + Cd_crossflow_final

    print(f"--- RESULTS (Alpha={alpha} deg) ---")
    print(f"Separated:      {bl.separated} at x = {x_arc[bl.sep_index] if bl.separated else 'None'}")
    print(f"Viscous Friction Force: {F_friction:.2f} N")
    print(f"Pressure (Wake) Force:  {F_pressure:.2f} N")
    print(f"Total Axial Drag Force: {F_drag_axial:.2f} N")
    print(f"Cross-Flow Force (Lateral): {F_cross_force:.2f} N")
    print(f"Total Drag Force Opposing Flow: {F_drag_total_opposing:.2f} N")
    print(f"--- Coefficients (Ref Area = {Ref_Area_Lateral:.2f} m^2) ---")
    print(f"Cd_Axial (Profile+Pressure): {Cd_axial_final:.5f}")
    print(f"Cd_Crossflow (Lateral Component): {Cd_crossflow_final:.5f}")
    print(f"Total Cd (Simplified Sum): {Cd_total_final:.5f}")

    fig, ax = plt.subplots(3, 1, figsize=(10, 12), sharex=True)

    margin = 0.5
    nx, nz = args.nx, args.nz
    x_grid = np.linspace(-margin, L_sub + margin, nx)
    z_grid = np.linspace(-L_sub/2, L_sub/2, nz)
    X, Z = np.meshgrid(x_grid, z_grid)
    flat_pts_grid = np.vstack([X.ravel(), np.zeros_like(X.ravel()), Z.ravel()]).T

    V_field = solver.get_total_velocity(flat_pts_grid)
    U = V_field[:, 0].reshape(X.shape)
    W = V_field[:, 2].reshape(X.shape)
    Vel_Mag = np.sqrt(U**2 + W**2)

    Cp = 1.0 - (Vel_Mag / V_inf)**2
    Cp_masked = Cp[np.isfinite(Cp)]
    if Cp_masked.size > 0:
        if args.cp_vmax is not None and args.cp_vmin is None:
            cmax = float(args.cp_vmax)
            cmin = -cmax
        else:
            cmin = args.cp_vmin if args.cp_vmin is not None else float(np.nanmin(Cp_masked))
            cmax = args.cp_vmax if args.cp_vmax is not None else float(np.nanmax(Cp_masked))
        if cmax <= cmin:
            cmax = cmin + 1e-6
    else:
        cmin, cmax = -1.0, 1.0
    Cp_plot = np.clip(Cp, cmin, cmax)
    im = ax[0].imshow(Cp_plot, extent=(x_grid.min(), x_grid.max(), z_grid.min(), z_grid.max()), origin="lower", cmap='coolwarm', vmin=cmin, vmax=cmax, aspect='equal')
    plt.colorbar(im, ax=ax[0], label='Cp')

    stream_density = 1.5
    st = ax[0].streamplot(X, Z, U, W, color='black', cmap='viridis', density=stream_density, linewidth=0.7)

    surf_x = solver.surface_pts[:,0]
    surf_z = solver.surface_pts[:,2]
    surf_x_closed = np.append(surf_x, surf_x[0])
    surf_z_closed = np.append(surf_z, surf_z[0])
    ax[0].plot(surf_x_closed, surf_z_closed, 'k-', linewidth=2.5, label='Submarine Hull')

    basis_x = solver.basis_locs[:,0]
    basis_z = solver.basis_locs[:,2]
    ax[0].plot(basis_x, basis_z, 'x', color='orange', markersize=6, label='Internal Basis')
    
    # Separation region highlight
    if bl.separated:
        sep_x_val = x_arc[bl.sep_index]
        sep_z_val = get_radius_func(L_sub, R_sub)(x_arc[bl.sep_index])
        # Draw a rectangle or line indicating separation
        ax[0].fill_between([sep_x_val, L_sub], -sep_z_val*1.1, sep_z_val*1.1, color='red', alpha=0.1, label='Separated Wake Region')
        
    ax[0].set_title(f"Flow Field and Pressure (Cp) - AoA={alpha}°")
    ax[0].set_aspect('equal')
    ax[0].set_xlabel("X [m]")
    ax[0].set_ylabel("Z [m]")
    ax[0].legend(loc='lower left', fontsize=8)
    
    # --- Plot 2: BL Parameters ---
    ax[1].plot(s, Ue, 'b-', label='Ue Potential')
    ax[1].set_ylabel('Ue [m/s]')
    ax[1].grid(True, alpha=0.3)
    
    ax2 = ax[1].twinx()
    color_theta = 'r-'
    ax2.plot(s, theta*1000, color_theta, label='Theta (mm)')
    ax2.set_ylabel('Theta [mm]', color=color_theta.strip('-'))
    ax2.tick_params(axis='y', color=color_theta.strip('-'), labelcolor=color_theta.strip('-'))
    
    if bl.separated:
        ax[1].axvline(x=s[bl.sep_index], color='k', linestyle='--', label='Separation')
        
    ax[1].legend(loc='upper left', fontsize=8)
    ax2.legend(loc='upper right', fontsize=8)
    ax[1].set_title('Boundary Layer Growth (Suction Side)')
    
    # --- Plot 3: Shape Factor and Friction ---
    ax[2].plot(s, H, 'g-', label='Shape Factor H')
    ax[2].axhline(y=2.5, color='gray', linestyle=':', label='Separation Threshold')
    ax[2].set_ylabel('H')
    ax[2].set_xlabel('Arc Length s [m]')
    ax[2].legend(loc='upper left', fontsize=8)
    
    ax3 = ax[2].twinx()
    color_cf = 'purple'
    ax3.plot(s, Cf, color=color_cf, label='Skin Friction Cf')
    ax3.set_ylabel('Cf', color=color_cf)
    ax3.tick_params(axis='y', color=color_cf, labelcolor=color_cf)
    ax3.legend(fontsize=8)
    ax3.grid(True, alpha=0.3)
    ax[2].set_title('Shape Factor and Skin Friction')

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()