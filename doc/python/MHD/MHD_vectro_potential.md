Yes, calculating the **Vector Potential $\mathbf{A}$** is often better for plotting (flux lines are simply contours of $r A_\phi$) and for coupling with other physics solvers.

### 1. Derivation of Vector Potential $A_\phi$

For a circular current loop of radius $a$ in the $xy$-plane carrying current $I$:
*   Due to cylindrical symmetry, $\mathbf{A}$ only has an azimuthal component $A_\phi(r, z)$.
*   There is a direct relationship between **Mutual Inductance ($M$)** and **Vector Potential**:
    The flux $\Phi$ through a circle of radius $r$ at height $z$ is:
    $$ \Phi = \oint \mathbf{A} \cdot d\mathbf{l} = A_\phi(r, z) \cdot 2\pi r $$
    We also know $\Phi = M(a, r, z) \cdot I$.
    Therefore:
    $$ A_\phi(r, z) = \frac{I \cdot M(a, r, z)}{2\pi r} $$

Using the elliptic integral formula for $M$ derived previously:
$$ m = k^2 = \frac{4 a r}{(a + r)^2 + z^2} $$
$$ M = \mu_0 \sqrt{a r} \left[ \left( \frac{2}{k} - k \right) K(m) - \frac{2}{k} E(m) \right] $$

Substituting this into the equation for $A_\phi$:
$$ A_\phi = \frac{\mu_0 I}{\pi} \sqrt{\frac{a}{r}} \frac{1}{k} \left[ (1 - \frac{m}{2}) K(m) - E(m) \right] $$

### 2. Derivation of B from A (The Curl Test)

The magnetic field is $\mathbf{B} = \nabla \times \mathbf{A}$. In cylindrical coordinates for $\mathbf{A} = (0, A_\phi, 0)$:

1.  **Radial Field ($B_r$):**
    $$ B_r = -\frac{\partial A_\phi}{\partial z} $$

2.  **Axial Field ($B_z$):**
    $$ B_z = \frac{1}{r} \frac{\partial (r A_\phi)}{\partial r} = \frac{A_\phi}{r} + \frac{\partial A_\phi}{\partial r} $$

---

### 3. The Python Script

This script implements `vector_potential_loop_rz`, computes it on a grid, calculates the curl numerically, and plots it against the analytic `field_loop_rz` results.

```python
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
```

### Explanation of the Code

1.  **`vector_potential_loop_rz`**: This is the new kernel. It handles the $1/\sqrt{r}$ term carefully by setting $A=0$ where $r < \epsilon$.
2.  **`compute_numerical_B_from_A`**:
    *   It creates a tiny 3-row grid around the target $Z$.
    *   **$B_r$** is calculated via central difference in $Z$: `(A_up - A_down) / (2*dz)`.
    *   **$B_z$** is calculated via `np.gradient` in $R$. We strictly follow the cylindrical curl formula: $\frac{1}{r} \frac{\partial (r A)}{\partial r}$.

The same vector-potential machinery is now used in the static **diamagnetic dipole bubble** demo (`demo_diamagnetic_dipole_bubble.py`). The physical picture is:

- **Coil**: a superconducting ring of radius $a$ at $(r=a, z=z_\text{coil})$ carrying current $I_\text{coil}$.
- **Dipole**: an on-axis magnetic dipole at $(r=0, z=z_\text{dip})$ with moment $m$ anti-aligned so that it creates a cavity (bubble) of opposite polarity field.
- **Plasma bubble surface**: approximated by a **magnetic separatrix** between field lines that close through the dipole region and those that close through the coil.

This is analogous to **potential flow theory**:

- In 2D potential flow, **streamlines** are level sets of a scalar stream function $\psi_\text{flow}(x,z)$.
- Here, in an axisymmetric magnetostatic configuration with $\mathbf{A} = A_\phi(r,z)\,\hat{\boldsymbol{\phi}}$, the **poloidal flux function**
  $$ \psi(r,z) = r\,A_\phi(r,z) $$
  plays the same role: magnetic field lines in the $(r,z)$ plane lie on contours of constant $\psi$.

#### 4.1 Two separatrix diagnostics

In the bubble demo we use two different, complementary diagnostics for the bubble boundary:

1. **Bz=0 contour + B-field streamlines (proxy)**
   - Compute the total field $B_r,B_z$ from `field_loop_rz` and `field_dipole_rz`.
   - Plot grey **streamlines** of $(B_r,B_z)$ in the $(r,z)$ plane.
   - Overlay a white contour where $B_z = 0$ as a simple indicator of where the axial component flips sign.
   - Optionally, highlight special streamlines starting near the on-axis stagnation point in front of the coil.

   This method is cheap and useful for a quick topological picture, but the $B_z=0$ curve is only a **proxy** for the true separatrix.

2. **Flux-function `psi = r Aphi` contour (more rigorous)**
   - Use the validated kernels `Aphi_loop_rz` and `Aphi_dipole_rz` (derived above and tested numerically) to compute
     $$ A_\phi^{\text{tot}} = A_\phi^{\text{coil}} + A_\phi^{\text{dipole}}, \qquad \psi = r\,A_\phi^{\text{tot}}. $$
   - Locate the **on-axis stagnation point** (or Bz=0 crossing) in front of the coil.
   - Take $\psi_\text{sep}$ as the value of $\psi$ on a small off-axis seed point near that stagnation location.
   - Plot the contour $\psi(r,z) = \psi_\text{sep}$ in a distinct color (cyan) as an approximation of the bubble surface / separatrix.

Because the underlying `Aphi_loop_rz` and `Aphi_dipole_rz` have been validated against `field_loop_rz` and `field_dipole_rz` via the curl test and radial profiles, these $\psi$-contours are consistent with the magnetic field and provide a more faithful picture of the magnetic topology.

#### 4.2 Coil and dipole vector potentials in the demo code

The actual implementations live in `inductance_core.py`:

- **Point dipole on axis** (`Aphi_dipole_rz`):
  - Based directly on the standard dipole potential in Coulomb gauge,
    $$ \mathbf{A}(\mathbf{R}) = \frac{\mu_0}{4\pi}\,\frac{\mathbf{m}\times\mathbf{R}}{|\mathbf{R}|^3}, $$
    which for an on-axis dipole with $\mathbf{m}\parallel\hat{z}$ yields a purely azimuthal component
    $$ A_\phi = \frac{\mu_0 m r}{4\pi R^3}, \quad R^2 = r^2 + (z-z_m)^2. $$
  - `curl(Aphi_dipole_rz)` matches `field_dipole_rz` along 1D scans to numerical precision.

- **Circular loop** (`Aphi_loop_rz`):
  - Derived from the mutual inductance expression and expressed via complete elliptic integrals using the formula documented in this file.
  - Implemented as
    $$ A_\phi(r,z) = \frac{\mu_0 I}{2\pi r}\,\sqrt{(a+r)^2 + (z-z_0)^2}\,\Big[(1-\tfrac{m}{2})K(m) - E(m)\Big], $$
    with
    $$ m = \frac{4 a r}{(a+r)^2 + (z-z_0)^2}. $$
  - In the demo, this kernel is validated by:
    - a local curl test at a single point (`check_Aphi_vs_B`), and
    - radial profile comparisons of `B_r,B_z` from `curl(Aphi)` versus the analytic `field_loop_rz` (coil-only) and `field_dipole_rz` (dipole-only).

Once both Aphi kernels pass these tests, the **separatrix drawn as a $\psi$-contour** in the dipole bubble demo becomes a reliable diagnostic of the bubble shape and its evolution when parameters (coil current, dipole moment, positions) are varied.