import numpy as np
import matplotlib.pyplot as plt
from scipy.special import legendre
import argparse
from scipy.interpolate import UnivariateSpline

# ==============================================================================
# 1. GEOMETRY & POTENTIAL FLOW (Your InternalBasisSolver)
# ==============================================================================

def get_naca_coords(naca='2412', n_points=200):
    """ Generates NACA 4-digit coordinates (Surface Sample Points) """
    m = int(naca[0]) / 100.0
    p = int(naca[1]) / 10.0
    t = int(naca[2:]) / 100.0

    # Cluster points at LE and TE
    beta = np.linspace(0, np.pi, n_points)
    x = 0.5 * (1.0 - np.cos(beta))
    
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1015 * x**4)
    
    yc = np.zeros_like(x)
    dyc_dx = np.zeros_like(x)
    for i, xi in enumerate(x):
        if xi < p:
            yc[i] = m / p**2 * (2*p*xi - xi**2)
            dyc_dx[i] = 2*m / p**2 * (p - xi)
        else:
            yc[i] = m / (1-p)**2 * ((1-2*p) + 2*p*xi - xi**2)
            dyc_dx[i] = 2*m / (1-p)**2 * (p - xi)
            
    theta = np.arctan(dyc_dx)
    xu = x - yt * np.sin(theta)
    yu = yc + yt * np.cos(theta)
    xl = x + yt * np.sin(theta)
    yl = yc - yt * np.cos(theta)
    
    X = np.concatenate((xl[::-1], xu[1:]))
    Y = np.concatenate((yl[::-1], yu[1:]))
    
    dx = np.gradient(X)
    dy = np.gradient(Y)
    norms = np.sqrt(dx**2 + dy**2)
    nx = -dy / norms
    ny = dx / norms
    
    return X, Y, nx, ny

GL_X, GL_W = np.polynomial.legendre.leggauss(64) 

class InternalBasisSolver:
    def __init__(self, surf_x, surf_y, surf_nx, surf_ny, n_basis=12):
        self.sx = surf_x
        self.sy = surf_y
        self.nx = surf_nx
        self.ny = surf_ny
        self.N_samples = len(surf_x)
        self.N_basis = n_basis
        
        self.integ_x = 0.5 * (GL_X + 1.0)
        self.integ_w = 0.5 * GL_W
        self.integ_N = len(GL_X)
        
        self.basis_vals_src = np.zeros((n_basis, self.integ_N))
        self.basis_vals_vtx = np.zeros((n_basis, self.integ_N))
        map_x = 2*self.integ_x - 1.0
        
        for k in range(n_basis):
            poly = legendre(k)
            self.basis_vals_src[k, :] = poly(map_x)
            self.basis_vals_vtx[k, :] = poly(map_x)

    def _get_induced_velocity_matrix(self, target_x, target_y):
        n_pts = len(target_x)
        Tx = target_x[:, np.newaxis]
        Ty = target_y[:, np.newaxis]
        Ix = self.integ_x[np.newaxis, :]
        Iw = self.integ_w[np.newaxis, :] 
        Iy = 0.0 
        Dx = Tx - Ix
        Dy = Ty - Iy
        R2 = Dx**2 + Dy**2
        InvR2 = 1.0 / (R2 + 1e-12)
        fact = (1.0 / (2*np.pi)) * Iw * InvR2
        K_src_u = fact * Dx
        K_src_v = fact * Dy
        K_vtx_u = fact * (-Dy)
        K_vtx_v = fact * Dx
        M_src_u = K_src_u @ self.basis_vals_src.T
        M_src_v = K_src_v @ self.basis_vals_src.T
        M_vtx_u = K_vtx_u @ self.basis_vals_vtx.T
        M_vtx_v = K_vtx_v @ self.basis_vals_vtx.T
        return M_src_u, M_src_v, M_vtx_u, M_vtx_v

    def solve(self, V_inf, alpha_deg):
        alpha = np.radians(alpha_deg)
        V_free = V_inf * np.array([np.cos(alpha), np.sin(alpha)])
        S_u, S_v, V_u, V_v = self._get_induced_velocity_matrix(self.sx, self.sy)
        Nx = self.nx[:, np.newaxis]
        Ny = self.ny[:, np.newaxis]
        A_src = S_u * Nx + S_v * Ny
        A_vtx = V_u * Nx + V_v * Ny
        A = np.hstack([A_src, A_vtx])
        b = - (V_free[0] * self.nx + V_free[1] * self.ny)
        
        row_closure = np.zeros(2 * self.N_basis)
        basis_integrals = np.sum(self.basis_vals_src * self.integ_w, axis=1)
        row_closure[0:self.N_basis] = basis_integrals * 10.0 
        
        row_kutta = np.zeros(2 * self.N_basis)
        row_kutta[self.N_basis:] = 1.0 * 10.0 
        
        A_aug = np.vstack([A, row_closure, row_kutta])
        b_aug = np.concatenate([b, [0, 0]])
        
        res = np.linalg.lstsq(A_aug, b_aug, rcond=None)
        self.coeffs = res[0]
        self.src_coeffs = self.coeffs[0:self.N_basis]
        self.vtx_coeffs = self.coeffs[self.N_basis:]
        self.V_inf = V_free

    def get_velocity_field(self, X, Y):
        flat_shape = X.ravel().shape
        x_flat = X.ravel()
        y_flat = Y.ravel()
        S_u, S_v, V_u, V_v = self._get_induced_velocity_matrix(x_flat, y_flat)
        u_ind = S_u @ self.src_coeffs + V_u @ self.vtx_coeffs
        v_ind = S_v @ self.src_coeffs + V_v @ self.vtx_coeffs
        U = u_ind + self.V_inf[0]
        V = v_ind + self.V_inf[1]
        return U.reshape(X.shape), V.reshape(X.shape)

# ==============================================================================
# 2. VISCOUS SOLVER (Drela / XFOIL style IBL)
# ==============================================================================

class ViscousSolver:
    def __init__(self, nu):
        self.nu = nu
        
    def solve_side(self, s, Ue):
        """
        Solves one side of the airfoil (stagnation -> TE).
        s: arc length (must start at 0)
        Ue: Edge velocity (must be positive)
        """
        n = len(s)
        theta = np.zeros(n) # Momentum thickness
        H = np.zeros(n)     # Shape factor (delta*/theta)
        Cf = np.zeros(n)    # Skin friction coeff
        
        # State
        is_turbulent = False
        is_separated = False
        
        # --- LAMINAR INITIALIZATION (Thwaites Stagnation) ---
        # dU/ds at stagnation
        if s[1] > 0:
            dUds0 = Ue[1] / s[1]
        else:
            dUds0 = 0
            
        if dUds0 > 1e-6:
            theta[0] = np.sqrt(0.075 * self.nu / dUds0)
        else:
            theta[0] = 0.0
        H[0] = 2.6 # Laminar stagnation shape factor
        
        for i in range(n - 1):
            if is_separated:
                # Wake approximation for separated flow
                theta[i+1] = theta[i] * 1.05 # Grow wake
                H[i+1] = 3.0 # Separated H
                Cf[i+1] = 0.0
                continue

            ds = s[i+1] - s[i]
            U = Ue[i]
            U_next = Ue[i+1]
            U_avg = 0.5 * (U + U_next)
            
            # --- LAMINAR REGION (Thwaites Method) ---
            if not is_turbulent:
                # 1. Integration: Theta^2
                # theta^2(x) = 0.45*nu/U^6 * int(U^5 dx)
                # Differential form: d(theta^2)/ds = 0.45*nu/U - 6/U * theta^2 * dU/ds
                # Robust step:
                term1 = theta[i]**2 * (U**6)
                term2 = 0.45 * self.nu * (U_avg**5) * ds
                theta[i+1] = np.sqrt((term1 + term2) / (U_next**6 + 1e-20))
                
                # 2. Thwaites Parameter lambda
                dUds = (U_next - U) / ds
                lam = (theta[i+1]**2 / self.nu) * dUds
                
                # 3. Shape Factor & Cf Correlations (Thwaites)
                lam = max(-0.1, min(0.25, lam)) # Clamp limits
                
                if lam >= 0:
                    H[i+1] = 2.61 - 3.75*lam + 5.24*lam**2
                    l_shear = 0.22 + 1.57*lam - 1.8*lam**2
                else:
                    H[i+1] = 2.088 + 0.0731/(lam + 0.14)
                    l_shear = 0.22 + 1.402*lam + (0.018*lam)/(lam+0.107)
                    
                Cf[i+1] = 2 * self.nu * l_shear / (U_next * theta[i+1] + 1e-12)
                
                # 4. TRANSITION CHECKS
                Re_theta = U_next * theta[i+1] / self.nu
                Re_x = U_next * s[i+1] / self.nu
                
                # A) Laminar Separation Bubble Check (Drela)
                # If H > 3.5 (approx), laminar layer separates.
                # If Re_theta is low, it forms a bubble.
                if H[i+1] > 3.5:
                    # LAMINAR SEPARATION DETECTED
                    # Model: Instant transition with momentum loss penalty
                    # In XFOIL this is complex (e^N). Here we simplify.
                    is_turbulent = True
                    # Jump conditions:
                    # Theta roughly conserved across bubble (plus drag penalty)
                    # H drops to turbulent value
                    H[i+1] = 1.4 
                    # Add empirical bubble drag penalty to theta
                    # delta_theta_bubble approx width of bubble
                    theta[i+1] += 0.05 * theta[i+1] 
                
                # B) Natural Transition (Michel's Criterion)
                # Re_theta_tr = 1.174 * (1 + 22400/Re_x) * Re_x^0.46
                elif Re_x > 10000: # avoid singularity at 0
                    Re_tr_michel = 1.174 * (1 + 22400.0/Re_x) * (Re_x**0.46)
                    if Re_theta > Re_tr_michel:
                        is_turbulent = True
                        H[i+1] = 1.4 # Jump to turbulent H

            # --- TURBULENT REGION (Head's Entrainment Method) ---
            else:
                th = theta[i]
                Sh = H[i]
                
                # 1. Skin Friction (Ludwieg-Tillmann)
                Re_th = U * th / self.nu
                if Re_th < 50: Re_th = 50
                Cf_turb = 0.246 * (10**(-0.678*Sh)) * (Re_th**-0.268)
                Cf[i+1] = Cf_turb
                
                # 2. Entrainment Shape Factor H1 (Head's correlation)
                if Sh < 1.1: Sh = 1.1
                H1 = 1.535 * (Sh - 0.7)**(-2.715) + 3.3
                
                # 3. Entrainment Coefficient Ce
                Ce = 0.0306 * (H1 - 3.0)**(-0.6169)
                
                # 4. Momentum Integral Eq: dTheta/ds
                # dTheta/ds = Cf/2 - (H+2)*(Theta/U)*dU/ds
                dUds = (U_next - U) / ds
                dTheta_ds = (Cf_turb / 2.0) - (Sh + 2.0) * (th / (U+1e-9)) * dUds
                
                theta[i+1] = th + dTheta_ds * ds
                
                # 5. Entrainment Eq: d(U*Theta*H1)/ds = U*Ce
                # Differential form for H1:
                # dH1/ds = (1/Theta) * [ Ce - H1 * ( dTheta/ds + (Theta/U)*dU/ds ) ]
                dH1_ds = (1.0/(th+1e-9)) * ( Ce - H1 * ( dTheta_ds + (th/(U+1e-9))*dUds ) )
                
                H1_next = H1 + dH1_ds * ds
                
                # Invert H1(H) to get new H
                term = (H1_next - 3.3) / 1.535
                if term > 0:
                    H[i+1] = 0.7 + term**(-1.0/2.715)
                else:
                    H[i+1] = 2.5 # Safety clamp
                    
                # 6. Turbulent Separation Check
                if H[i+1] > 2.8: # Typical separation shape factor
                    is_separated = True
        
        return theta, H, Cf

# ==============================================================================
# 3. MAIN DRIVER
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(description="Airfoil Solver with Drela-style Viscous Drag")
    parser.add_argument("--naca", default="2412", type=str)
    parser.add_argument("--alpha", default=4.0, type=float, help="Angle of attack (deg)")
    parser.add_argument("--re", default=1e6, type=float, help="Reynolds Number")
    parser.add_argument("--field", default="cp", choices=["vel", "cp"], help="Field to plot in flow map: velocity magnitude or Cp")
    parser.add_argument("--n-points", default=100, type=int, help="Surface sample points for geometry")
    parser.add_argument("--n-basis", default=8, type=int, help="Number of internal polynomial basis terms")
    parser.add_argument("--grid-nx", default=120, type=int, help="Flow-field grid resolution in x")
    parser.add_argument("--grid-ny", default=100, type=int, help="Flow-field grid resolution in y")
    parser.add_argument("--x-margin", default=0.2, type=float, help="Margin before/after chord for flow plot")
    parser.add_argument("--y-extent", default=0.5, type=float, help="Half-height of flow plot window")
    parser.add_argument("--v-inf", default=50.0, type=float, help="Freestream speed [m/s]")
    parser.add_argument("--chord", default=1.0, type=float, help="Reference chord length [m]")
    parser.add_argument("--surface-scalar", default="cf", choices=["none", "cp", "cf", "h", "theta", "vel"],
                        help="Scalar to plot along surface normals on airfoil plot")
    parser.add_argument("--scalar-scale", default=0.05, type=float, help="Scale factor for normal offsets")
    args = parser.parse_args()

    # Settings
    V_inf = args.v_inf
    alpha = args.alpha
    chord = args.chord
    nu = (V_inf * chord) / args.re # Calc viscosity from Desired Re
    
    print(f"--- Running Analysis: NACA {args.naca}, Alpha={alpha}, Re={args.re:.1e}, V_inf={V_inf} m/s ---")

    # 1. Geometry & Potential Flow
    sx, sy, nx, ny = get_naca_coords(args.naca, n_points=args.n_points)
    solver = InternalBasisSolver(sx, sy, nx, ny, n_basis=args.n_basis)
    solver.solve(V_inf, alpha)
    
    # Get Surface Velocity
    # We evaluate slightly off surface? No, Internal Basis is stable on surface.
    U_surf, V_surf = solver.get_velocity_field(sx, sy)
    V_mag_surf = np.sqrt(U_surf**2 + V_surf**2)
    
    # 2. Process Surface Data for BL Solver
    # Need to split into Upper and Lower surfaces starting from Stagnation Point
    idx_stag = np.argmin(V_mag_surf)
    
    # Identify split indices (0 to Stag = Lower, Stag to End = Upper)
    # Note: Points are usually TE_Lower -> Nose -> TE_Upper
    # Lower Surface: Index idx_stag down to 0
    # Upper Surface: Index idx_stag up to end
    
    def extract_surface(indices):
        x = sx[indices]
        y = sy[indices]
        u = V_mag_surf[indices]
        # Calculate Arc Length
        dist = np.sqrt(np.diff(x)**2 + np.diff(y)**2)
        s = np.concatenate(([0], np.cumsum(dist)))
        return s, u

    idx_lower = np.arange(idx_stag, -1, -1)
    idx_upper = np.arange(idx_stag, len(sx))
    
    s_l, ue_l = extract_surface(idx_lower)
    s_u, ue_u = extract_surface(idx_upper)
    
    # 3. Solve Boundary Layers
    viscous = ViscousSolver(nu)
    
    theta_u, H_u, Cf_u = viscous.solve_side(s_u, ue_u)
    theta_l, H_l, Cf_l = viscous.solve_side(s_l, ue_l)
    
    # 4. Drag Calculation (Squire-Young)
    # Cd = 2 * theta_te * (U_te / V_inf) ^ ((H_te + 5) / 2)
    # Apply to both sides
    
    def get_cd_side(theta, H, Ue):
        th_te = theta[-1]
        H_te = H[-1]
        U_te = Ue[-1]
        
        # Safety for separated flow
        if H_te > 4.0: H_te = 4.0 
        
        return 2 * th_te * (U_te / V_inf)**((H_te + 5)/2)

    Cd_u = get_cd_side(theta_u, H_u, ue_u)
    Cd_l = get_cd_side(theta_l, H_l, ue_l)
    Cd_total = (Cd_u + Cd_l) / chord
    
    # Lift (Kutta-Joukowski on Potential)
    # This assumes potential lift is maintained (no stall)
    # Circulation = sum(gamma * ds)? 
    # InternalBasisSolver doesn't output total Gamma directly easily.
    # Integrate Cp: Cl = Integral (Cp_lower - Cp_upper) dx
    Cp = 1.0 - (V_mag_surf / V_inf)**2
    # Simple trapezoidal integration of Cp around the foil projected on X
    Cl_integ = -np.trapz(Cp * nx, sx) + np.trapz(Cp * ny, sy) # Generalized force
    # Actually just -trapz(Cp * nx) is usually approx for Cl if alpha small
    # Let's do simplified:
    Cl_total = 0.0
    for i in range(len(sx)-1):
        dx = sx[i+1] - sx[i]
        dy = sy[i+1] - sy[i]
        cp_avg = 0.5*(Cp[i]+Cp[i+1])
        # Force vector = -Cp * n * ds
        # Lift is component perpendicular to V_inf
        # n = (nx, ny)
        # ds = sqrt(dx^2+dy^2)
        # Fx = -Cp * nx * ds
        # Fy = -Cp * ny * ds
        # Lift = Fy*cos(a) - Fx*sin(a)
        ds = np.sqrt(dx**2 + dy**2)
        nx_i = 0.5*(nx[i]+nx[i+1])
        ny_i = 0.5*(ny[i]+ny[i+1])
        
        Fx = -cp_avg * nx_i * ds
        Fy = -cp_avg * ny_i * ds
        
        Cl_total += Fy * np.cos(np.radians(alpha)) - Fx * np.sin(np.radians(alpha))
        
    Cl_total /= chord

    print(f"Total Cl: {Cl_total:.4f}")
    print(f"Total Cd: {Cd_total:.5f} ({Cd_total*10000:.1f} counts)")
    print(f"  Upper Cd: {Cd_u/chord:.5f}")
    print(f"  Lower Cd: {Cd_l/chord:.5f}")

    # Flow field for density-style visualization (like airfoil3/SlenderBodyTheory5-a)
    xg = np.linspace(-args.x_margin, 1.0 + args.x_margin, args.grid_nx)
    yg = np.linspace(-args.y_extent, args.y_extent, args.grid_ny)
    XG, YG = np.meshgrid(xg, yg)
    UG, VG = solver.get_velocity_field(XG, YG)
    Vel_Mag = np.sqrt(UG**2 + VG**2)
    if args.field == "cp":
        Cp_field = 1.0 - (Vel_Mag / V_inf)**2
    # Mask inside airfoil approximation
    def get_thick(x):
        t = int(args.naca[2:]) / 100.0
        return 5*t*(0.2969*np.sqrt(np.abs(x))-0.1260*x-0.3516*x**2+0.2843*x**3-0.1015*x**4)
    inside = (np.abs(YG) < get_thick(XG)) & (XG > 0) & (XG < 1)
    Vel_Mag[inside] = np.nan
    if args.field == "cp":
        Cp_field[inside] = np.nan

    fig_flow = plt.figure(figsize=(12, 5))
    axf = fig_flow.add_subplot(1, 1, 1)
    if args.field == "cp":
        im = axf.imshow(Cp_field, extent=(xg.min(), xg.max(), yg.min(), yg.max()),
                        origin="lower", cmap='coolwarm', aspect='equal')
        axf.set_title(f"Cp Field (AoA={alpha}°)")
        cbar_label = "Cp"
    else:
        im = axf.imshow(Vel_Mag, extent=(xg.min(), xg.max(), yg.min(), yg.max()),
                        origin="lower", cmap='plasma', aspect='equal')
        axf.set_title(f"Velocity Magnitude Field (AoA={alpha}°)")
        cbar_label = "|V|"
    axf.streamplot(XG, YG, UG, VG, color='k', density=1.4, linewidth=0.7)
    axf.plot(sx, sy, 'k-', linewidth=1.5, label='Airfoil')
    # Basis function quadrature nodes (placement along chord)
    node_stride = max(1, len(solver.integ_x) // max(6, solver.N_basis * 2))
    axf.plot(solver.integ_x[::node_stride], np.zeros_like(solver.integ_x[::node_stride]), 'rx', label='Internal Basis Nodes')
    # Optional scalar plotted along normals
    if args.surface_scalar != "none":
        scalar_plot = np.full_like(sx, np.nan, dtype=float)
        if args.surface_scalar == "cp":
            scalar_plot = Cp
        elif args.surface_scalar == "cf":
            scalar_plot[idx_upper] = Cf_u
            scalar_plot[idx_lower] = Cf_l
        elif args.surface_scalar == "h":
            scalar_plot[idx_upper] = H_u
            scalar_plot[idx_lower] = H_l
        elif args.surface_scalar == "theta":
            scalar_plot[idx_upper] = theta_u
            scalar_plot[idx_lower] = theta_l
        elif args.surface_scalar == "vel":
            scalar_plot = V_mag_surf
        finite_mask = np.isfinite(scalar_plot)
        if np.any(finite_mask):
            maxabs = np.nanmax(np.abs(scalar_plot))
            if maxabs <= 0:
                maxabs = 1.0
            scale = args.scalar_scale / maxabs
            x2 = sx + nx * scalar_plot * scale
            y2 = sy + ny * scalar_plot * scale
            for i in range(len(sx)):
                if finite_mask[i]:
                    axf.plot([sx[i], x2[i]], [sy[i], y2[i]], color='magenta', linewidth=0.6, alpha=0.7)
    axf.set_xlabel("x/c")
    axf.set_ylabel("y/c")
    axf.legend(loc='upper right')
    plt.colorbar(im, ax=axf, label=cbar_label)
    plt.tight_layout()

    # ================= VISUALIZATION =================
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(3, 2)
    
    # 1. Cp Distribution
    ax1 = fig.add_subplot(gs[0, :])
    ax1.plot(sx, Cp, 'k-', linewidth=1, alpha=0.5)
    ax1.plot(sx[idx_upper], Cp[idx_upper], 'b-', linewidth=2, label='Upper')
    ax1.plot(sx[idx_lower], Cp[idx_lower], 'r-', linewidth=2, label='Lower')
    ax1.invert_yaxis()
    ax1.set_ylabel("Cp")
    ax1.set_title(f"Pressure Distribution (Cl={Cl_total:.2f})")
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. Skin Friction
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.plot(s_u, Cf_u, 'b-', label='Upper')
    ax2.plot(s_l, Cf_l, 'r-', label='Lower')
    ax2.set_ylabel("Cf")
    ax2.set_xlabel("Arc Length s")
    ax2.set_title("Skin Friction (Separation when Cf -> 0)")
    ax2.set_ylim(0, 0.01)
    ax2.grid(True)
    
    # 3. Shape Factor (H) - Shows Laminar/Turbulent/Separation
    ax3 = fig.add_subplot(gs[1, 1])
    ax3.plot(s_u, H_u, 'b-', label='Upper')
    ax3.plot(s_l, H_l, 'r-', label='Lower')
    ax3.axhline(2.6, color='gray', linestyle=':', label='Laminar (2.6)')
    ax3.axhline(1.4, color='k', linestyle='--', label='Turbulent (1.4)')
    ax3.set_ylabel("Shape Factor H")
    ax3.set_title("Boundary Layer State (H)")
    ax3.legend()
    ax3.grid(True)
    
    # 4. Momentum Thickness
    ax4 = fig.add_subplot(gs[2, :])
    ax4.plot(s_u, theta_u*1000, 'b-', label='Upper')
    ax4.plot(s_l, theta_l*1000, 'r-', label='Lower')
    ax4.set_ylabel("Theta [mm]")
    ax4.set_xlabel("Arc Length s")
    ax4.set_title(f"Momentum Thickness (Drag Proportional to Theta_TE) - Cd={Cd_total:.5f}")
    ax4.grid(True)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()