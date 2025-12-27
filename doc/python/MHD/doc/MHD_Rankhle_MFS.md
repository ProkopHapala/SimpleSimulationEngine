
No I'm trhinking that it is better to represent the plasma ball as a dipole such that magnetic filed is zero on the surface of the ball. That is I do not create explicit coils, I just vary simple dipole inside plasma (which is fake not true magnetic field inside plasma), just as a cheap way how to have same separatrix on the surface of plasma. I saw this trick when simulating e.g. ships or airfoils using potential flow theory, that I embed linear conbination of dipoles or sources/drain into the flow and vary thyr amplitude to get zero velocity perpendicular to the surface of the airfoild. Now with plasma we can have the B-field of the seed coil (instead of constant velocity field), and then plasma can either
1) expand by increasing the value of its dipole (we do not know where exactly is the surface, but the surface will grow as well somehow)
2) put some constran points at the surface of the sphere on which we want to have surface of the plasma and solve the system of equations to match it.

I do not know it it is easy to incorporate the constant flux condition for the metalic cage coils into method (1), with method (2) it is kinda natural choice to use the RHS method where surface of the plasma and surface of the metalic cage are the two surfaces where we want to match the B-field  

---

## Gemini

This is a standard and very elegant technique in mathematical physics, often called the **Method of Fundamental Solutions (MFS)** or the **Rankine Body** approach in fluid dynamics.

In Magnetics, this is based on the **Magnetic Flux Function** (Stream Function) $\Psi$.

### The Concept: The "Magnetic Bubble"

In 2D Axisymmetric coordinates ($r, z$), the magnetic field lines are exactly the contour lines of the scalar function:
$$ \Psi(r, z) = r \cdot A_\phi(r, z) $$

*   **Property:** The difference in $\Psi$ between two contours tells you exactly how much Magnetic Flux is between them.
*   **The Separatrix:** If you sum a uniform background field (lines going straight up) and a Dipole (lines looping around), there exists a specific contour $\Psi = C$ that separates the "closed" dipole lines from the "open" background lines.
*   **The Plasma Surface:** This separatrix **IS** the surface of your plasma. Inside it, the field topology is different (or zero, if treated as a perfect exclusion zone).

### Why this is better for you
1.  **Speed:** You calculate 1 dipole instead of 50 rings.
2.  **Smoothness:** You get a naturally smooth, physically realistic "teardrop" or spherical shape without defining mesh points.
3.  **No Leaks:** The separatrix is mathematically closed. Flux cannot cross a $\Psi$ contour.

---

### The Math: Integrating Dipoles with the Cage

You need to couple the **Dipole Moment ($m$)** with the **Cage Currents ($I_{cage}$)**.

#### 1. The Dipole Formula (The "Fake" Source)
For a magnetic dipole moment $m$ located at $(0, z_{dip})$, the Vector Potential $A_\phi$ at point $(r, z)$ is:
$$ A_{\phi, dipole} = \frac{\mu_0 m}{4\pi} \frac{r}{(r^2 + (z-z_{dip})^2)^{3/2}} $$

And the Flux Function $\Psi = r A_\phi$:
$$ \Psi_{dipole}(r, z) = \frac{\mu_0 m}{4\pi} \frac{r^2}{(r^2 + (z-z_{dip})^2)^{3/2}} $$

#### 2. The Cage Interaction (Mutual Inductance)
To solve the system, you need the flux through a Cage Ring (Radius $R_c$, Position $Z_c$) caused by the Dipole. This is simply $\Psi$ evaluated at the wire!
$$ \Phi_{ring \leftarrow dipole} = 2\pi \cdot \Psi_{dipole}(R_c, Z_c) $$
$$ M_{ring, dipole} = \frac{\mu_0 R_c^2}{2 (R_c^2 + (Z_c-z_{dip})^2)^{3/2}} $$

*(Note: This assumes the ring is centered on the axis. If you check units, $m$ has units $A \cdot m^2$. Flux is $T \cdot m^2$. Formula is correct.)*

### 3. The New Algorithm (Method 2: Fixed Surface)

You proposed **Method 2**: *Put control points on the surface and solve.*

Let's say you want the plasma to have a specific radius $R_{p}$ at the equator ($z=z_{dip}$).

**The System of Equations:**
1.  **Cage Physics:** Flux conservation.
    $$ L_{cage} I_{cage} + M_{cd} \cdot m = \Phi_{initial} - M_{cs} I_{seed} $$
    *(Where $m$ is the dipole magnitude, treated as a current-like variable)*

2.  **Plasma Physics:** The plasma surface is defined by the separatrix where the **Total Vector Potential is Zero** (assuming $A_\phi=0$ on axis).
    Actually, valid diamagnetism usually implies $\mathbf{B}_{total} = 0$ *inside*.
    For the "Stream Function" method, we set the condition that the **Net Flux** at the plasma boundary radius $R_p$ is zero (cancellation of external and internal sources).
    $$ \Psi_{total}(R_p, z_{dip}) = 0 $$
    $$ \Psi_{seed} + \Psi_{cage}(I_{cage}) + \Psi_{dipole}(m) = 0 $$

**Result:**
You have a linear system for unknowns $[I_{cage}, m]$.
Once solved, you plot the contour $\Psi_{total}(r, z) = 0$. That contour is your plasma shape.

---

### Python Proof-of-Concept

This script:
1.  Creates a background field (Seed Coil).
2.  Places a **Dipole** at the center.
3.  Calculates the dipole strength $m$ required to create a bubble of radius $R_{plasma}$.
4.  Plots the **Streamlines ($\Psi$)** showing the separatrix (the plasma skin).

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ellipk, ellipe

MU0 = 4 * np.pi * 1e-7

# --- 1. KERNELS ---

def psi_loop(r, z, a, z0, I):
    """
    Magnetic Stream Function Psi = r * A_phi for a current loop.
    Flux between axis and point (r,z) is 2*pi*Psi.
    """
    # Standard Vector Potential Kernel
    r = np.atleast_1d(r)
    z = np.atleast_1d(z)
    eps = 1e-12
    r_safe = np.where(r < eps, eps, r)
    z_rel = z - z0
    
    denom = (a + r_safe)**2 + z_rel**2
    m_param = (4 * a * r_safe) / denom
    k = np.sqrt(m_param)
    K = ellipk(m_param)
    E = ellipe(m_param)
    
    # A_phi formula
    prefactor = (MU0 * I) / (2 * np.pi * r_safe) * np.sqrt(denom)
    A_phi = prefactor * ((1 - m_param/2)*K - E)
    
    Psi = r_safe * A_phi
    return Psi

def psi_dipole(r, z, m_moment, z_dip):
    """
    Stream Function for a Magnetic Dipole aligned with Z-axis.
    Psi = (mu0 * m / 4pi) * r^2 / (r^2 + dz^2)^(3/2)
    """
    z_rel = z - z_dip
    denom = (r**2 + z_rel**2)**1.5
    # Avoid singularity
    denom = np.where(denom < 1e-12, 1e-12, denom)
    
    Psi = (MU0 * m_moment / (4 * np.pi)) * (r**2 / denom)
    return Psi

# --- 2. SOLVER ---

def solve_dipole_strength(R_boundary, Z_boundary, I_seed, R_seed, Z_seed):
    """
    Finds dipole moment 'm' such that Psi_total(R_bound, Z_bound) = 0.
    Psi_seed + Psi_dipole = 0
    """
    # 1. Calc Psi from Seed at boundary
    psi_s = psi_loop(R_boundary, Z_boundary, R_seed, Z_seed, I_seed)[0]
    
    # 2. Calc Unit Psi from Dipole (m=1)
    psi_d_unit = psi_dipole(R_boundary, Z_boundary, 1.0, 0.0) # Dipole at 0,0
    
    # 3. Solve: psi_s + m * psi_d_unit = 0
    m_solution = -psi_s / psi_d_unit
    return m_solution

# --- 3. EXECUTION ---

# Geometry
R_SEED = 2.0
Z_SEED = 0.0
I_SEED = 1.0e6  # 1 MA

# Desired Plasma Bubble
R_PLASMA = 0.8  # We want the plasma surface to be at r=0.8
Z_PLASMA = 0.0  # At the equator

# Solve
m_sol = solve_dipole_strength(R_PLASMA, Z_PLASMA, I_SEED, R_SEED, Z_SEED)
print(f"Required Dipole Moment: {m_sol:.2e} A*m^2")

# Grid for Plotting
r_vals = np.linspace(0, 3.0, 200)
z_vals = np.linspace(-2.0, 2.0, 200)
R, Z = np.meshgrid(r_vals, z_vals)

# Compute Total Psi
Psi_S = psi_loop(R, Z, R_SEED, Z_SEED, I_SEED).reshape(R.shape)
Psi_D = psi_dipole(R, Z, m_sol, 0.0).reshape(R.shape)
Psi_Tot = Psi_S + Psi_D

# --- PLOTTING ---
plt.figure(figsize=(8, 10))

# Plot Seed Coil
plt.plot(R_SEED, Z_SEED, 'ro', markersize=10, label='Seed Coil')
plt.plot(0, 0, 'b*', markersize=15, label='Plasma Dipole')

# Plot Contours (Streamlines)
# We specifically want the contour Psi=0 (The Separatrix)
levels = np.linspace(np.min(Psi_Tot), np.max(Psi_Tot), 50)
plt.contour(R, Z, Psi_Tot, levels=levels, colors='gray', linewidths=0.5, alpha=0.5)

# Highlight the Plasma Boundary (Psi=0)
plt.contour(R, Z, Psi_Tot, levels=[0], colors='blue', linewidths=3, linestyles='solid')

plt.title(f"Rankine Body Method (Plasma Dipole)\nSeparatrix at R={R_PLASMA}")
plt.xlabel("Radius R [m]")
plt.ylabel("Axial Z [m]")
plt.axis('equal')
plt.grid(True)
plt.legend()

# Annotation
plt.text(0.1, 0.2, "Plasma Region\n(Diamagnetic)", color='blue', fontweight='bold')
plt.text(1.5, 0.5, "External Field", color='gray')

plt.show()
```

### Interpretation
1.  **The Blue Line:** This is the contour $\Psi=0$.
    *   Everything **inside** this line is dominated by the dipole. In a real diamagnet, the total B-field here would be zero (perfect cancellation).
    *   Everything **outside** is the external field flowing around the object.
2.  **The Shape:** Notice it forms a perfect "Bubble" or spherical shape naturally. You didn't have to define nodes or segments.
3.  **Expansion:** To simulate expansion, you simply increase the "Desired R" in the solver. The dipole moment $m$ will increase, pushing the Blue Line outward.
4.  **Cage Integration:** If you add a Cage, you just add `Psi_Cage` to the sum. The Cage current reacts to keep its own flux constant, and the Dipole reacts to keep the Plasma boundary at $R$.

### Challenges with this approach
*   **Shape Control:** A single dipole creates a sphere/apple shape. If your nozzle is long and thin, a single dipole won't fit well.
*   **Solution:** Use a **Linear Array of Dipoles** along the Z-axis (e.g., 5-10 dipoles).
    *   You can tune the weights of these dipoles to fit a specific long cigar shape (like fitting a spline).
    *   This is still much faster and more robust than surface current rings.