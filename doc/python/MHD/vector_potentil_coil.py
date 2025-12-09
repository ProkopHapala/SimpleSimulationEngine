import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ellipk, ellipe

MU0 = 4 * np.pi * 1e-7

# ==========================================
# 1. KERNELS (A and B)
# ==========================================

def vector_potential_loop_rz(a, z0, I, r_grid, z_grid, eps=1e-12):
    """
    Calculates Azimuthal Vector Potential A_phi.
    Relation: Phi = A_phi * 2 * pi * r
    """
    r = np.asarray(r_grid, dtype=float)
    z = np.asarray(z_grid, dtype=float)
    z_rel = z - z0
    
    # Avoid singularity at r=0
    # A_phi is 0 on axis, but formula divides by sqrt(r)
    r_safe = np.where(r < eps, eps, r)
    
    # Elliptic parameter m = k^2
    # m = 4 a r / ((a+r)^2 + z^2)
    denom = (a + r_safe)**2 + z_rel**2
    m = (4 * a * r_safe) / denom
    k = np.sqrt(m)
    
    K = ellipk(m)
    E = ellipe(m)
    
    # Formula: A_phi = (mu0 I / pi) * sqrt(a/r) * (1/k) * [ (1-m/2)K - E ]
    # Note: 1/k = sqrt(denom) / (2 sqrt(a r))
    # sqrt(a/r) * (1/k) = sqrt(a/r) * sqrt(denom) / (2 sqrt(a) sqrt(r))
    #                   = sqrt(denom) / (2 r)
    
    prefactor = (MU0 * I) / (2 * np.pi * r_safe) * np.sqrt(denom)
    term_bracket = (1 - m/2)*K - E
    
    A_phi = prefactor * term_bracket
    
    # Enforce 0 on axis
    A_phi = np.where(r < eps, 0.0, A_phi)
    
    return A_phi

def field_loop_rz(a, z0, I, r_grid, z_grid, eps=1e-12):
    """Analytic B-field (Br, Bz) for comparison."""
    r = np.asarray(r_grid, dtype=float)
    z = np.asarray(z_grid, dtype=float)
    z_rel = z - z0
    r_safe = np.where(r < eps, eps, r)

    denom = (a + r_safe)**2 + z_rel**2
    m = (4 * a * r_safe) / denom
    K = ellipk(m)
    E = ellipe(m)
    
    factor = MU0 * I / (2 * np.pi * np.sqrt(denom))
    alpha2 = a**2 + r_safe**2 + z_rel**2
    beta2  = (a - r_safe)**2 + z_rel**2
    
    # B_r
    term1_r = -K
    term2_r = (alpha2 / beta2) * E
    Br = factor * (z_rel / r_safe) * (term1_r + term2_r)
    
    # B_z
    term1_z = K
    term2_z = ((a**2 - r_safe**2 - z_rel**2) / beta2) * E
    Bz = factor * (term1_z + term2_z)
    
    # Fix on axis
    Br = np.where(r < eps, 0.0, Br)
    
    return Br, Bz

# ==========================================
# 2. NUMERICAL CURL & TESTING
# ==========================================

def compute_numerical_B_from_A(a, z0, I, r_vals, z_fixed, dr=1e-5, dz=1e-5):
    """
    Computes B = Curl(A) using finite differences centered at z_fixed.
    """
    # Create a small 2D mesh around the line of interest
    # We need neighbors in R and Z to compute derivatives
    
    r_grid, z_grid = np.meshgrid(r_vals, [z_fixed - dz, z_fixed, z_fixed + dz])
    
    # Compute A on this strip
    A_phi = vector_potential_loop_rz(a, z0, I, r_grid, z_grid)
    
    # 1. Compute Br = - dA_phi / dz
    # Central difference in Z direction (axis 0)
    dA_dz = (A_phi[2, :] - A_phi[0, :]) / (2 * dz)
    Br_num = -dA_dz
    
    # 2. Compute Bz = (1/r) * d(r*A)/dr
    # Central difference in R direction is tricky because R varies.
    # Let's use np.gradient along axis 1 (R-axis)
    r_A = r_grid * A_phi
    
    # Slice out the middle row (z_fixed)
    r_A_mid = r_A[1, :]
    r_vals_mid = r_grid[1, :]
    
    # Gradient of (r*A) w.r.t r
    d_rA_dr = np.gradient(r_A_mid, r_vals_mid)
    
    # Bz = (1/r) * d_rA_dr
    # Handle r=0 for division
    Bz_num = np.divide(d_rA_dr, r_vals_mid, out=np.zeros_like(d_rA_dr), where=r_vals_mid!=0)
    
    # Special limit for Bz on axis: Bz(0) = 2 * dA/dr (L'Hopital)
    # But usually the finite diff handles it okay if grid is fine.
    # Let's overwrite index 0 if r[0]==0 with the analytical limit check later.
    
    return Br_num, Bz_num

# ==========================================
# 3. RUN TEST
# ==========================================

# Parameters
R_COIL = 1.0
I_COIL = 1.0e6
Z_COIL = 0.0
Z_EVAL = 0.5  # Scan line offset from coil plane

# Grid
r_eval = np.linspace(0.0, 3.0, 100)

# A. Analytic B
r_mesh, z_mesh = np.meshgrid(r_eval, [Z_EVAL])
Br_ana, Bz_ana = field_loop_rz(R_COIL, Z_COIL, I_COIL, r_mesh, z_mesh)
Br_ana = Br_ana[0]
Bz_ana = Bz_ana[0]

# B. Numerical B from A
Br_num, Bz_num = compute_numerical_B_from_A(R_COIL, Z_COIL, I_COIL, r_eval, Z_EVAL)

# C. Vector Potential for Plotting
A_val = vector_potential_loop_rz(R_COIL, Z_COIL, I_COIL, r_mesh, z_mesh)[0]

# --- PLOTTING ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10), sharex=True)

# Plot 1: Vector Potential Profile
ax1.plot(r_eval, A_val, 'k-', lw=2, label=r'$A_\phi$')
ax1.set_ylabel(r'Vector Potential $A_\phi$ [T$\cdot$m]')
ax1.set_title(f'Vector Potential at Z = {Z_EVAL} (Coil R={R_COIL})')
ax1.grid(True)
ax1.legend()

# Plot 2: Comparison of B fields
ax2.plot(r_eval, Bz_ana, 'b-', lw=3, alpha=0.5, label=r'$B_z$ (Analytic)')
ax2.plot(r_eval, Bz_num, 'b--', label=r'$B_z$ (Num Curl)')

ax2.plot(r_eval, Br_ana, 'r-', lw=3, alpha=0.5, label=r'$B_r$ (Analytic)')
ax2.plot(r_eval, Br_num, 'r--', label=r'$B_r$ (Num Curl)')

ax2.set_xlabel('Radius r [m]')
ax2.set_ylabel('Magnetic Field [T]')
ax2.set_title('Validation: B = curl(A) vs Analytic B')
ax2.legend()
ax2.grid(True)

# Calculate error
err_z = np.mean(np.abs(Bz_num - Bz_ana))
err_r = np.mean(np.abs(Br_num - Br_ana))
print(f"Mean Absolute Error Bz: {err_z:.2e} T")
print(f"Mean Absolute Error Br: {err_r:.2e} T")

plt.tight_layout()
plt.show()