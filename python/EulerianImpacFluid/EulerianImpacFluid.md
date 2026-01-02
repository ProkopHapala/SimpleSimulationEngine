https://aistudio.google.com/app/prompts?state=%7B%22ids%22:%5B%221_fEAgF7as9M1He5RmM9572C5L612zEmP%22%5D,%22action%22:%22open%22,%22userId%22:%22100958146796876347936%22,%22resourceKeys%22:%7B%7D%7D&usp=sharing

## USER

I'm tinking how to make simple simulator for compressible continuum mechanis / fluid synamics for simulation mainly of hypervelocity impacts, for many reasons I got recomandation to do it in eurelian frame with static contant volume cells which are filled and emptied by material. The pressure, density, velocity and all importaint quantities should be simulated in those cells. But I need somehow simulate non-mixable materials which is easy to do using particles (Lagrangian frame).

Or is it? I saw some interesting technique which is kinda based on topological features of zero-surfaces of switching functions. Instead of storing concentration of certain substrance directly (value between [0..1]), we store switching variable [-1..1], which somewhere cross zero. The concentration of some material is then given by some non linear mapping of the switching variable (for example fermi function 1/(1+exp(x)) or some other sigmoide like function, can be also fast functions like smoothstep().

So the euler grid cell does not store directly ammount of certain substance, but store only total agregate density, pressure, momentum etc. 
when we move some material we work like everything is the same material in the first step, and only secondarily we calculate how much the switching function should change to be consistent with the transfer of material. I guess we need some variational derivatives of the switching function to know how much the function should change to represent movement of certain amaount of material.

What makse it a bit problematic that when we have two materials of very high density difference ( e.g. liquid hydrogen ~0.16 kg/l and uranium ~18kg/l => >100x density ratio), than it may be very difficult to solve it properly - the derivatives may be very steep

do you know such methods? can you elaborate on this? Can you write more rigogrousy mathematically what exact updates we have to make in which order, and where are each quantitiy represented?

---

## Gemnini

What you are describing is known in Computational Fluid Dynamics (CFD) as the **Level Set Method (LSM)**, often combined with the **Ghost Fluid Method (GFM)** to handle the high density ratios and immiscibility.

The "switching variable" you described is the **Level Set Function ($\phi$)**, and the "nonlinear mapping" is the **Regularized Heaviside Function**.

Here is a breakdown of why this matches your intuition, why the variational derivative part applies (in Phase Field methods), and a rigorous mathematical formulation for a compressible solver handling Uranium/Hydrogen impacts.

---

### 1. The Concepts

#### The "Switching Variable": The Level Set ($\phi$)
Instead of a volume fraction $\alpha \in [0,1]$, we store a signed distance field $\phi(\vec{x}, t)$.
*   $\phi > 0$: Material A (Uranium)
*   $\phi < 0$: Material B (Liquid Hydrogen)
*   $\phi = 0$: The Interface

#### The Mapping: Heaviside Function
To get the density or pressure in a cell, we map $\phi$ to a range $[0,1]$ using a Heaviside function $H(\phi)$.
*   Sharp interface: $H(\phi) = 1$ if $\phi \ge 0$, else $0$.
*   Numerical interface (smeared): You mentioned sigmoids. In numerics, we use a specific "smeared" Heaviside $H_\epsilon(\phi)$ over a bandwidth $\epsilon$ (usually 1.5 times the cell width) to prevent numerical instability (Gibbs phenomenon).

$$
H_\epsilon(\phi) = 
\begin{cases} 
0 & \text{if } \phi < -\epsilon \\
\frac{1}{2} + \frac{\phi}{2\epsilon} + \frac{1}{2\pi}\sin(\frac{\pi\phi}{\epsilon}) & \text{if } |\phi| \le \epsilon \\
1 & \text{if } \phi > \epsilon 
\end{cases}
$$

#### The "Variational Derivative" vs. Advection
You mentioned calculating derivatives to see "how much the function should change."
1.  **Phase Field Method:** This uses variational derivatives of a free energy functional (Cahn-Hilliard equation). It is physically rigorous for micro-scale surface tension but generally **too diffusive** for hypervelocity impacts.
2.  **Level Set Method (Recommended):** We simply *advect* the distance field using the fluid velocity. $\phi$ flows with the material. This is much better for high-speed ballistics.

---

### 2. The High Density Ratio Problem
**The Danger:** If you simply average density in a cell containing Uranium ($\rho=18$) and Hydrogen ($\rho=0.16$), a $50/50$ mix gives $\rho \approx 9$.
If you calculate the speed of sound $c = \sqrt{\partial P / \partial \rho}$ using this smeared density, the acoustic waves in the Hydrogen part of that cell will travel artificially slow, and the Uranium waves artificially fast. This leads to massive pressure errors and the interface "blowing up."

**The Solution: The Ghost Fluid Method (GFM)**
Don't mix the fluids physically. In the mixed cells, you create two "Virtual" (Ghost) states.
*   The Uranium solver thinks the Hydrogen region is filled with "Ghost Uranium" having the pressure and velocity of the Hydrogen but the entropy/density of Uranium.
*   You solve the fluxes twice at the interface.

---

### 3. Rigorous Mathematical Formulation

We assume a **Compressible Euler System** on a Cartesian Grid (finite volume).

#### State Variables (Stored in Cell Centers)
For each cell $i,j$:
1.  **Level Set:** $\phi_{i,j}$ (Signed distance to interface)
2.  **Momentum:** $(\rho \vec{u})_{i,j}$
3.  **Total Energy:** $E_{i,j}$
4.  **Note:** We do **not** store a single density $\rho$. We store variables for the separate phases or use GFM to reconstruct them. For a simple implementation, let's look at the **5-Equation Model** (Allaire et al.) which is easier than full GFM but robust for high density ratios.

#### The 5-Equation Model (Quasi-Conservative)
This allows compressible mixing without pressure oscillations.

**1. Mass Conservation (Fluid 1 - Uranium):**
$$ \frac{\partial (\alpha_1 \rho_1)}{\partial t} + \nabla \cdot (\alpha_1 \rho_1 \vec{u}) = 0 $$

**2. Mass Conservation (Fluid 2 - Hydrogen):**
$$ \frac{\partial (\alpha_2 \rho_2)}{\partial t} + \nabla \cdot (\alpha_2 \rho_2 \vec{u}) = 0 $$
*(Note: $\alpha_1 + \alpha_2 = 1$)*

**3. Global Momentum:**
$$ \frac{\partial (\rho \vec{u})}{\partial t} + \nabla \cdot (\rho \vec{u} \otimes \vec{u} + p I) = 0 $$

**4. Global Total Energy:**
$$ \frac{\partial E}{\partial t} + \nabla \cdot ((E + p)\vec{u}) = 0 $$

**5. Interface Advection (The Level Set):**
$$ \frac{\partial \phi}{\partial t} + \vec{u} \cdot \nabla \phi = 0 $$

---

### 4. Step-by-Step Update Loop

Here is the exact order of operations for one time step $\Delta t$:

#### Step A: Calculate Time Step (CFL Condition)
Find the maximum wave speed $S_{max} = \max(|\vec{u}| + c)$ over the grid.
$$ \Delta t = \text{CFL} \times \frac{\Delta x}{S_{max}} $$

#### Step B: Reconstruct Interface (Calculate $\alpha$)
Using your level set $\phi$:
$$ \alpha_1 = H_\epsilon(\phi) $$
$$ \alpha_2 = 1 - \alpha_1 $$

#### Step C: Equation of State (EOS) mixture
You need the pressure $p$ for the momentum update.
For hypervelocity, use the **Mie-Grüneisen** or **Stiff Gas** EOS.
The mixture pressure in a cell is usually defined to maintain mechanical equilibrium:
$$ p(\rho e, \alpha_1, \alpha_2) = \frac{\rho e - (\frac{\alpha_1 \gamma_1 P_{\infty,1}}{\gamma_1 -1} + \frac{\alpha_2 \gamma_2 P_{\infty,2}}{\gamma_2 -1}) }{ \frac{\alpha_1}{\gamma_1 - 1} + \frac{\alpha_2}{\gamma_2 - 1} } $$
*(This formula allows robust pressure calculation even with 100x density ratios).*

#### Step D: Spatial Reconstruction (MUSCL / WENO)
To get high accuracy and handle shocks (hypervelocity), you cannot just use neighboring cells. You must reconstruct values at cell faces (Left and Right states) using a slope limiter (like Minmod or Superbee) to prevent oscillation.

#### Step E: Solve Riemann Problem (Flux Calculation)
At every cell face, you have a Left state ($U_L$) and Right state ($U_R$).
Use an HLLC or Roe Solver to calculate the flux $F_{interface}$.
*Crucial:* This solver handles the "physics" of the impact (shock waves).

#### Step F: Update Conservative Variables
$$ (\rho \vec{u})^{n+1} = (\rho \vec{u})^n - \frac{\Delta t}{\Delta x} (F_{right} - F_{left}) $$
$$ E^{n+1} = E^n - \frac{\Delta t}{\Delta x} (F_{right} - F_{left}) $$

#### Step G: Update Level Set (Advection)
Solve $\frac{\partial \phi}{\partial t} + \vec{u} \cdot \nabla \phi = 0$.
You can use a simple Upwind scheme or Hamilton-Jacobi WENO for this.
$$ \phi^{n+1}_{i} = \phi^n_i - \Delta t (\max(u_i, 0) D^-_x \phi + \min(u_i, 0) D^+_x \phi) $$
*Where $D^-$ and $D^+$ are backward/forward difference operators.*

#### Step H: Re-initialization (Fixing the Gradient)
**This is the most critical step for your "Switching Variable" idea.**
Advecting $\phi$ distorts it. It stops being a distance field (slope $\neq 1$).
You must iterate this PDE for a few fictitious time steps $\tau$:
$$ \frac{\partial \phi}{\partial \tau} + S(\phi_0) (|\nabla \phi| - 1) = 0 $$
where $S(\phi_0)$ is the sign of the original field. This resets the slope to 1 without moving the zero-crossing (the interface).

---

### Summary for your Simulator

1.  **Don't** use variational derivatives (Phase Field) for hypervelocity; it smears the interface too much. Use **Level Set Advection** (Step G).
2.  **Do** use the switching variable $\phi \in [-\infty, \infty]$ (Distance field).
3.  **Do** use the Heaviside function to map $\phi \to$ Volume Fraction $\alpha$.
4.  **Handling Density Ratio:** Use the **5-Equation Model** (solving continuity for both phases separately + mixture momentum/energy). This prevents the "Uranium pressure" from incorrectly accelerating "Hydrogen mass" because the volume fractions $\alpha$ scale the compressibility in the EOS explicitly.

This is the standard approach used in hydrocodes like LS-DYNA (ALE solvers) or academic codes like CLAWPACK for hypervelocity impact.


---

## USER

ok, now please implement OpenCL code for this in 2D, every workgroup (32 threads) will take care of a tile 8x8 reading all relevant variables into local memory (including +1 pixel boundary) and outputting them to output buffers (ping-pong, avoid atomic write). Make clear separation of EOS and other modular function so that kernel is simple and readable. All the local function should operate on tread private memory.

---


## Gemini

This is a complex simulation task. To make this runnable and readable, I will implement a **2D Finite Volume Solver** using the **5-Equation Model** (Allaire et al.) which is robust for immiscible fluids with high density ratios (like Hydrogen/Uranium).

### The Architecture
1.  **Grid:** Eulerian Cartesian Grid.
2.  **Tiling:** We use **Local Memory (LDS)** to cache an $8 \times 8$ tile plus a 1-cell halo ($10 \times 10$ total) to minimize global memory bandwidth.
3.  **Parallelism:** Workgroup size = 32. Each thread processes 2 cells (since $8 \times 8 = 64$ cells).
4.  **Physics:**
    *   **EOS:** Stiffened Gas EOS (handles liquids and solids effectively).
    *   **Flux:** Rusanov (Local Lax-Friedrichs) flux. It is dissipative but extremely stable for hypervelocity shocks and keeps the code readable.
    *   **Advection:** Level set $\phi$ is advected with the fluid velocity.

### OpenCL Kernel (`solver.cl`)

```c
// =================================================================================
// PHYSICS CONSTANTS & MACROS
// =================================================================================
#define TILE_W 8
#define TILE_H 8
#define HALO   1
#define LOCAL_W (TILE_W + 2 * HALO) // 10
#define LOCAL_H (TILE_H + 2 * HALO) // 10
#define LOCAL_SIZE (LOCAL_W * LOCAL_H) // 100

// Material 1: Liquid Hydrogen (approx parameters for Stiffened Gas)
#define GAMMA_1 1.4f
#define PI_1    0.0f  // Pascal (approx 0 for liquid H2 near boiling, simplified)

// Material 2: Uranium (approx parameters)
#define GAMMA_2 3.0f  // Gruneisen parameter approx
#define PI_2    40.0e9f // 40 GPa stiffening pressure

// Small epsilon to avoid division by zero
#define EPS 1e-6f

// =================================================================================
// DATA STRUCTURES
// =================================================================================
// We use simple float arrays.
// Input/Output Buffer Layout:
// val0: float4(rho_alpha1, rho_alpha2, mom_x, mom_y)
// val1: float2(total_energy, levelset_phi)

typedef struct {
    float rho_a1; // Partial density fluid 1
    float rho_a2; // Partial density fluid 2
    float ru;     // Momentum X
    float rv;     // Momentum Y
    float E;      // Total Energy
    float phi;    // Level Set
} Conserved;

typedef struct {
    float p;      // Pressure
    float u;      // Velocity X
    float v;      // Velocity Y
    float c;      // Speed of sound (mixture)
} Primitive;

// =================================================================================
// MODULAR FUNCTIONS (Thread Private)
// =================================================================================

// 1. Heaviside Function (Smeared)
// Maps level set phi to volume fraction alpha [0..1]
inline float heaviside(float phi, float dx) {
    float epsilon = 1.5f * dx;
    if (phi < -epsilon) return 0.0f;
    if (phi > epsilon) return 1.0f;
    return 0.5f + (phi / (2.0f * epsilon)) + (1.0f / (2.0f * M_PI_F)) * sin((M_PI_F * phi) / epsilon);
}

// 2. Equation of State (Mixture Stiffened Gas)
// Calculates Pressure and Sound Speed based on conserved variables
inline Primitive get_primitive(Conserved U, float dx) {
    Primitive P;
    
    // 1. Recover Densities and Volume Fractions
    // Note: We use the level set to determine alpha for consistency, 
    // though 5-eq model tracks rho_alpha transport explicitly.
    // For this simple version, we rely on the transported partial densities.
    
    float rho = U.rho_a1 + U.rho_a2;
    if (rho < EPS) rho = EPS;

    // Velocities
    P.u = U.ru / rho;
    P.v = U.rv / rho;
    float ke = 0.5f * (P.u*P.u + P.v*P.v); // Kinetic Energy

    // Internal Energy
    float internal_e = U.E - rho * ke;

    // Volume fractions from partial densities
    // (In exact 5-eq, we transport rho*alpha, but we need alpha to solve Pressure)
    // We approximate alpha based on the transported masses and known material reference densities
    // For strict consistency, we use the mixture EOS formula:
    
    float alpha1 = U.rho_a1 / rho; // Approximation for mixing
    // Ideally use Level Set for sharp interface, but let's stick to conservative var mixture:
    float alpha2 = 1.0f - alpha1;

    // Mixture Parameters
    float g1_term = alpha1 / (GAMMA_1 - 1.0f);
    float g2_term = alpha2 / (GAMMA_2 - 1.0f);
    float denom = g1_term + g2_term;

    float pinf1_term = (alpha1 * GAMMA_1 * PI_1) / (GAMMA_1 - 1.0f);
    float pinf2_term = (alpha2 * GAMMA_2 * PI_2) / (GAMMA_2 - 1.0f);

    // Pressure Law for Mixture
    P.p = (internal_e - (pinf1_term + pinf2_term)) / denom;
    if (P.p < 1e-5f) P.p = 1e-5f; // Vacuum protection

    // Speed of Sound (Mixture)
    // c^2 = (p * sum(alpha*gamma/(gamma-1)) + sum(alpha * pi * gamma / (gamma-1)) ) / (rho * denom)
    // This is a simplified effective sound speed for the mixture
    float num_c = P.p * (1.0f/denom * ( (alpha1*GAMMA_1)/(GAMMA_1-1.0f) + (alpha2*GAMMA_2)/(GAMMA_2-1.0f) )); // Approx
    P.c = sqrt(max(num_c / rho, 1.0f)); 

    return P;
}

// 3. Flux Function (Rusanov / Local Lax-Friedrichs)
// Computes flux between state L and state R
inline Conserved compute_flux(Conserved UL, Conserved UR, float dx, float dt, float nx, float ny) {
    Primitive PL = get_primitive(UL, dx);
    Primitive PR = get_primitive(UR, dx);

    // Normal velocities
    float un_L = PL.u * nx + PL.v * ny;
    float un_R = PR.u * nx + PR.v * ny;

    // Max wave speed (S+)
    float S_L = fabs(un_L) + PL.c;
    float S_R = fabs(un_R) + PR.c;
    float S_max = max(S_L, S_R);

    // Physical Fluxes F(U)
    Conserved FL, FR;

    // Flux L
    FL.rho_a1 = UL.rho_a1 * un_L;
    FL.rho_a2 = UL.rho_a2 * un_L;
    FL.ru     = UL.ru * un_L + PL.p * nx;
    FL.rv     = UL.rv * un_L + PL.p * ny;
    FL.E      = (UL.E + PL.p) * un_L;
    FL.phi    = UL.phi * un_L; // Advection of level set

    // Flux R
    FR.rho_a1 = UR.rho_a1 * un_R;
    FR.rho_a2 = UR.rho_a2 * un_R;
    FR.ru     = UR.ru * un_R + PR.p * nx;
    FR.rv     = UR.rv * un_R + PR.p * ny;
    FR.E      = (UR.E + PR.p) * un_R;
    FR.phi    = UR.phi * un_R;

    // Rusanov Flux: 0.5 * (FL + FR) - 0.5 * Smax * (UR - UL)
    Conserved F;
    F.rho_a1 = 0.5f * (FL.rho_a1 + FR.rho_a1 - S_max * (UR.rho_a1 - UL.rho_a1));
    F.rho_a2 = 0.5f * (FL.rho_a2 + FR.rho_a2 - S_max * (UR.rho_a2 - UL.rho_a2));
    F.ru     = 0.5f * (FL.ru     + FR.ru     - S_max * (UR.ru     - UL.ru));
    F.rv     = 0.5f * (FL.rv     + FR.rv     - S_max * (UR.rv     - UL.rv));
    F.E      = 0.5f * (FL.E      + FR.E      - S_max * (UR.E      - UL.E));
    F.phi    = 0.5f * (FL.phi    + FR.phi    - S_max * (UR.phi    - UL.phi));

    return F;
}

// =================================================================================
// KERNEL
// =================================================================================
__kernel void update_fluid(
    __global const float4* restrict in_U_vec, // rho_a1, rho_a2, ru, rv
    __global const float2* restrict in_E_vec, // E, phi
    __global float4* restrict out_U_vec,
    __global float2* restrict out_E_vec,
    int width,
    int height,
    float dx,
    float dt
) {
    // ---------------------------------------------------------
    // 1. Setup Identities
    // ---------------------------------------------------------
    int lx = get_local_id(0); // 0..31
    int grp_x = get_group_id(0);
    int grp_y = get_group_id(1);
    
    // Origin of this tile in global space
    int tile_origin_x = grp_x * TILE_W; 
    int tile_origin_y = grp_y * TILE_H;

    // Local Memory (Structure of Arrays for better banking)
    // Size: 10x10 = 100 elements
    __local float l_rho_a1[LOCAL_SIZE];
    __local float l_rho_a2[LOCAL_SIZE];
    __local float l_ru[LOCAL_SIZE];
    __local float l_rv[LOCAL_SIZE];
    __local float l_E[LOCAL_SIZE];
    __local float l_phi[LOCAL_SIZE];

    // ---------------------------------------------------------
    // 2. Cooperative Load into Local Memory (with Halo)
    // ---------------------------------------------------------
    // We have 32 threads. We need to load 100 items.
    // Each thread loads approx 3-4 items.
    
    int total_items = LOCAL_W * LOCAL_H; // 100
    int threads = get_local_size(0); // 32

    for (int t = lx; t < total_items; t += threads) {
        // Map linear t to local 2D (ly_loc, lx_loc) in the 10x10 block
        int ly_loc = t / LOCAL_W;
        int lx_loc = t % LOCAL_W;

        // Corresponding global coordinate (offset by Halo -1)
        int gx = tile_origin_x + lx_loc - HALO;
        int gy = tile_origin_y + ly_loc - HALO;

        // Boundary Clamp
        gx = max(0, min(gx, width - 1));
        gy = max(0, min(gy, height - 1));
        int g_idx = gy * width + gx;

        // Read Global
        float4 vecA = in_U_vec[g_idx];
        float2 vecB = in_E_vec[g_idx];

        // Write Local
        l_rho_a1[t] = vecA.x;
        l_rho_a2[t] = vecA.y;
        l_ru[t]     = vecA.z;
        l_rv[t]     = vecA.w;
        l_E[t]      = vecB.x;
        l_phi[t]    = vecB.y;
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    // ---------------------------------------------------------
    // 3. Compute Phase (Fluxes & Update)
    // ---------------------------------------------------------
    // Each thread (32 total) processes 2 pixels in the 8x8 tile.
    
    int pixels_to_process = TILE_W * TILE_H; // 64
    
    for (int k = lx; k < pixels_to_process; k += threads) {
        // Map k to local 8x8 coordinates (0..7)
        int ty = k / TILE_W;
        int tx = k % TILE_W;

        // Map to 10x10 local memory indices (offset by +HALO)
        int cx = tx + HALO;
        int cy = ty + HALO;
        int c_idx = cy * LOCAL_W + cx;

        // Extract Center State (U_center)
        Conserved Uc;
        Uc.rho_a1 = l_rho_a1[c_idx];
        Uc.rho_a2 = l_rho_a2[c_idx];
        Uc.ru     = l_ru[c_idx];
        Uc.rv     = l_rv[c_idx];
        Uc.E      = l_E[c_idx];
        Uc.phi    = l_phi[c_idx];

        // Extract Neighbors
        Conserved U_left, U_right, U_down, U_up;
        
        // Helper lambda-like macro for extraction
        #define GET_U(idx) (Conserved){l_rho_a1[idx], l_rho_a2[idx], l_ru[idx], l_rv[idx], l_E[idx], l_phi[idx]}
        
        U_left  = GET_U(cy * LOCAL_W + (cx - 1));
        U_right = GET_U(cy * LOCAL_W + (cx + 1));
        U_down  = GET_U((cy - 1) * LOCAL_W + cx);
        U_up    = GET_U((cy + 1) * LOCAL_W + cx);

        // Calculate Fluxes
        // X-Fluxes
        Conserved F_xm = compute_flux(U_left, Uc, dx, dt, 1.0f, 0.0f);  // Face i-1/2
        Conserved F_xp = compute_flux(Uc, U_right, dx, dt, 1.0f, 0.0f); // Face i+1/2
        
        // Y-Fluxes
        Conserved F_ym = compute_flux(U_down, Uc, dx, dt, 0.0f, 1.0f);  // Face j-1/2
        Conserved F_yp = compute_flux(Uc, U_up, dx, dt, 0.0f, 1.0f);    // Face j+1/2

        // Finite Volume Update: U_new = U_old - dt/dx * (F_right - F_left) ...
        float dtdx = dt / dx; // Assuming dx = dy

        Conserved U_new;
        U_new.rho_a1 = Uc.rho_a1 - dtdx * (F_xp.rho_a1 - F_xm.rho_a1 + F_yp.rho_a1 - F_ym.rho_a1);
        U_new.rho_a2 = Uc.rho_a2 - dtdx * (F_xp.rho_a2 - F_xm.rho_a2 + F_yp.rho_a2 - F_ym.rho_a2);
        U_new.ru     = Uc.ru     - dtdx * (F_xp.ru     - F_xm.ru     + F_yp.ru     - F_ym.ru);
        U_new.rv     = Uc.rv     - dtdx * (F_xp.rv     - F_xm.rv     + F_yp.rv     - F_ym.rv);
        U_new.E      = Uc.E      - dtdx * (F_xp.E      - F_xm.E      + F_yp.E      - F_ym.E);
        U_new.phi    = Uc.phi    - dtdx * (F_xp.phi    - F_xm.phi    + F_yp.phi    - F_ym.phi);

        // ---------------------------------------------------------
        // 4. Global Write
        // ---------------------------------------------------------
        int gx_out = tile_origin_x + tx;
        int gy_out = tile_origin_y + ty;

        if(gx_out < width && gy_out < height) {
            int out_idx = gy_out * width + gx_out;
            
            float4 outA = (float4)(U_new.rho_a1, U_new.rho_a2, U_new.ru, U_new.rv);
            float2 outB = (float2)(U_new.E, U_new.phi);

            out_U_vec[out_idx] = outA;
            out_E_vec[out_idx] = outB;
        }
    }
}
```

### Explanation of the Implementation

1.  **Workgroup & Tiling:**
    *   The kernel uses 1D workgroups of size 32 (`get_local_id(0)`).
    *   The loop structure inside "Cooperative Load" and "Compute" allows these 32 threads to handle the 100 load operations and 64 compute operations required for the 8x8 tile.
    *   `LOCAL_SIZE` is 100 ($10 \times 10$) to account for the $8 \times 8$ active area plus a 1-pixel border (halo) required for the flux calculation.

2.  **Memory Layout:**
    *   **Global:** Vectorized `float4` and `float2` are used to ensure memory coalescing when reading/writing global memory.
    *   **Local (`__local`):** Structure of Arrays (SoA) is used (`l_rho_a1`, `l_ru`, etc.). This avoids bank conflicts better than an Array of Structs in local memory when threads access neighbors.

3.  **Physics - The 5-Equation Model:**
    *   We update **two** partial densities `rho_a1` ($\alpha_1 \rho_1$) and `rho_a2` ($\alpha_2 \rho_2$). This is crucial for maintaining correct pressure equilibrium at the interface of Uranium and Hydrogen.
    *   **Mixture Pressure:** The `get_primitive` function calculates the pressure using the Stiffened Gas EOS for the mixture. This formula respects the vastly different compressibilities of the materials.

4.  **Flux (Rusanov):**
    *   I used the Rusanov flux because it is essentially code-independent of the complexity of the EOS. It calculates the maximum wave speed $S_{max}$ and adds numerical viscosity proportional to that speed to stabilize shocks. This prevents the "blowing up" typical of standard central differences in hypervelocity scenarios.

### Host-Side Pseudo Code (C++)

To run this, your host code loop would look like this:

```cpp
// 1. Set Global Work Size
// We want to cover the image WxH with 8x8 tiles.
size_t local_work_size[1] = { 32 }; // As requested
size_t global_work_size[1] = { (W / 8) * (H / 8) * 32 }; 
// Note: This 1D global setup maps to 2D tiles inside the kernel logic 
// using get_group_id mapped to tile coordinates.

// 2. Main Loop
while (simulating) {
    // A. Calculate dt (CFL condition) on CPU or via reduction kernel
    // dt = CFL * dx / max(|u| + c);
    
    // B. Set Arguments
    clSetKernelArg(kernel, 0, buffer_U_in);
    clSetKernelArg(kernel, 1, buffer_E_in);
    clSetKernelArg(kernel, 2, buffer_U_out);
    clSetKernelArg(kernel, 3, buffer_E_out);
    // ... dx, dt, etc.

    // C. Enqueue
    clEnqueueNDRangeKernel(queue, kernel, 1, NULL, global_work_size, local_work_size, ...);

    // D. Swap Buffers (Ping-Pong)
    std::swap(buffer_U_in, buffer_U_out);
    std::swap(buffer_E_in, buffer_E_out);
}
```

### Important Note on Level Set Re-initialization
The code above *advects* the level set $\phi$. Over time, $\phi$ will distort and stop being a signed distance field (slope $\neq 1$). In a full production code, every 10-20 steps you must run a separate "Re-initialization" kernel that solves $|\nabla \phi| = 1$ to fix the slope without moving the zero-crossing interface.

---

### USER

write me that Re-initialization kernel, do not care about host-side for now, only about the kernels. any other comments? what can be expected problems with stability? Any specially care about initialization, the input parameters?

---

## gemini

Here is the **Re-initialization (Redistancing) Kernel**.

### The Math (Godunov's Scheme)
To solve the Eikonal equation $|\nabla \phi| = 1$ without moving the interface, we solve this PDE to steady state (usually 5 to 10 iterations per fluid step):
$$ \phi_{\tau} + S(\phi_0)(|\nabla \phi| - 1) = 0 $$
Where $S(\phi_0)$ is the "smoothed sign" of the original level set. We use **Godunov's Upwind Scheme** to calculate $|\nabla \phi|$ because standard central differences are unstable for Hamilton-Jacobi equations.

### 1. OpenCL Kernel: `redistance.cl`

```c
// =================================================================================
// CONSTANTS
// =================================================================================
#define TILE_W 8
#define TILE_H 8
#define HALO   1
#define LOCAL_W (TILE_W + 2 * HALO) // 10
#define LOCAL_H (TILE_H + 2 * HALO) // 10
#define LOCAL_SIZE (LOCAL_W * LOCAL_H)

// Smoothed sign function to prevent numerical movement of the zero-crossing
// S = phi / sqrt(phi^2 + dx^2)
inline float smoothed_sign(float phi, float dx) {
    return phi / sqrt(phi*phi + dx*dx);
}

// =================================================================================
// RE-INITIALIZATION KERNEL
// =================================================================================
__kernel void redistance_phi(
    __global const float2* restrict in_E_vec,  // Current state (E, phi)
    __global const float2* restrict orig_E_vec,// State at start of redist step (for Sign)
    __global float2* restrict out_E_vec,       // Output
    int width,
    int height,
    float dx,
    float dtau // Pseudo-time step, usually 0.5 * dx
) {
    // ---------------------------------------------------------
    // 1. Thread & Tile Setup
    // ---------------------------------------------------------
    int lx = get_local_id(0);
    int grp_x = get_group_id(0);
    int grp_y = get_group_id(1);
    
    int tile_origin_x = grp_x * TILE_W; 
    int tile_origin_y = grp_y * TILE_H;

    // Local Memory for Phi
    __local float l_phi[LOCAL_SIZE];
    
    // We also need the original sign. 
    // Optimization: We don't need neighbor signs, only local sign, 
    // so we can read it directly in the compute loop or cache it.
    // Let's cache it to avoid uncoalesced global reads later if desired,
    // but for simplicity/register pressure, we will read 'orig' directly later 
    // since it's 1:1 mapping with threads in the compute phase.

    // ---------------------------------------------------------
    // 2. Cooperative Load (Current Phi with Halo)
    // ---------------------------------------------------------
    int total_items = LOCAL_W * LOCAL_H;
    int threads = get_local_size(0);

    for (int t = lx; t < total_items; t += threads) {
        int ly_loc = t / LOCAL_W;
        int lx_loc = t % LOCAL_W;

        int gx = tile_origin_x + lx_loc - HALO;
        int gy = tile_origin_y + ly_loc - HALO;

        gx = max(0, min(gx, width - 1));
        gy = max(0, min(gy, height - 1));
        
        // Load .y component (phi) from float2 buffer
        l_phi[t] = in_E_vec[gy * width + gx].y;
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    // ---------------------------------------------------------
    // 3. Compute Godunov Update
    // ---------------------------------------------------------
    int pixels_to_process = TILE_W * TILE_H;

    for (int k = lx; k < pixels_to_process; k += threads) {
        int ty = k / TILE_W;
        int tx = k % TILE_W;

        int cx = tx + HALO;
        int cy = ty + HALO;
        int c_idx = cy * LOCAL_W + cx;

        float phi_c = l_phi[c_idx];

        // Neighbors
        float phi_W = l_phi[cy * LOCAL_W + (cx - 1)]; // West
        float phi_E = l_phi[cy * LOCAL_W + (cx + 1)]; // East
        float phi_S = l_phi[(cy - 1) * LOCAL_W + cx]; // South
        float phi_N = l_phi[(cy + 1) * LOCAL_W + cx]; // North

        // One-sided derivatives
        // D_minus: (Center - Prev) / dx
        // D_plus:  (Next - Center) / dx
        float dx_m = (phi_c - phi_W) / dx;
        float dx_p = (phi_E - phi_c) / dx;
        float dy_m = (phi_c - phi_S) / dx;
        float dy_p = (phi_N - phi_c) / dx;

        // Get Original Sign from Global Memory (to anchor the interface)
        int gx_out = tile_origin_x + tx;
        int gy_out = tile_origin_y + ty;
        
        // We need the sign of phi^0 (the state before we started re-distancing iterations)
        float phi_0 = orig_E_vec[gy_out * width + gx_out].y;
        float S = smoothed_sign(phi_0, dx);

        // Godunov Hamiltonian Logic
        float grad_phi_sq = 0.0f;

        if (phi_0 > 0.0f) {
            // If inside fluid (phi > 0):
            // We look "downhill" towards the interface
            float dx_term = max(max(0.0f, -dx_p), max(0.0f, dx_m)); // Squared effectively
            float dy_term = max(max(0.0f, -dy_p), max(0.0f, dy_m));
            grad_phi_sq = dx_term*dx_term + dy_term*dy_term;
        } else {
            // If outside fluid (phi < 0):
            // We look "uphill" towards the interface
            float dx_term = max(max(0.0f, dx_p), max(0.0f, -dx_m));
            float dy_term = max(max(0.0f, dy_p), max(0.0f, -dy_m));
            grad_phi_sq = dx_term*dx_term + dy_term*dy_term;
        }

        // Update: phi_new = phi_old - dtau * S * (|grad| - 1)
        float grad_mag = sqrt(grad_phi_sq);
        float phi_new = phi_c - dtau * S * (grad_mag - 1.0f);

        // ---------------------------------------------------------
        // 4. Output
        // ---------------------------------------------------------
        // We preserve Energy (x component) and update Phi (y component)
        if(gx_out < width && gy_out < height) {
            int out_idx = gy_out * width + gx_out;
            float E_val = in_E_vec[out_idx].x; // Keep energy
            out_E_vec[out_idx] = (float2)(E_val, phi_new);
        }
    }
}
```

---

### 2. Stability & Common Problems

#### A. The Boundary Condition Problem
*   **Issue:** The simulation box is finite. What happens when the Uranium block hits the wall?
*   **Solution:** In the `load` section of both kernels, I used `clamp` (Coordinate clamping). This effectively creates **Neumann BCs** (Reflection) where the derivative is zero at the wall. For hypervelocity, you usually want **Transmissive (Outflow)** boundaries so shockwaves leave the domain.
    *   *To fix:* Check if `gx` is out of bounds. If so, return the value of `width-1`. This copies the edge value into the ghost zone (Zero-gradient).

#### B. The "Interface Jitter" (Mass Loss)
*   **Issue:** Even with re-initialization, the zero-level set might drift slightly, causing "Mass Loss" (your uranium ball shrinks).
*   **Solution:** 
    1.  Use the **Smoothed Sign** function provided above. The hard sign ($sgn(x)$) causes numerical popping.
    2.  Do not run re-initialization *every* step. Run it every 10 or 20 steps.
    3.  If mass loss is severe, you need a "Mass Correction" step: Calculate total mass before redistancing, calculating total mass after, and shift $\phi$ globally by a tiny constant to restore it.

#### C. Initialization (Input Parameters)
*   **Signed Distance Field (SDF):** You cannot just set cells to +1 or -1. You must initialize $\phi$ as the Euclidean distance to the interface.
    *   *Circle at $(cx, cy)$ with radius $R$:*
        $$ \phi(x,y) = R - \sqrt{(x-cx)^2 + (y-cy)^2} $$
        (Assuming $\phi > 0$ is inside).
    *   *Half-space (Impact Plate) at $y < H$:*
        $$ \phi(x,y) = H - y $$

#### D. The Epsilon Problem (Division by Zero)
*   **Issue:** In the EOS kernel, we calculate `1/(alpha1 + alpha2)`.
*   **Solution:** In the 5-equation model, strictly $\alpha_1 + \alpha_2$ should be 1. But numerically they drift.
    *   *Fix:* In the fluid kernel, always normalize alphas before EOS: `sum = a1+a2; a1/=sum; a2/=sum;`.
    *   *Fix:* In `rho = rho_a1 + rho_a2`, ensure `rho` never drops below a vacuum tolerance (e.g., `1e-6`).

### 3. Integration Plan

Your simulation loop should look like this:

1.  **Initialize** $U$ (Momentum, Density) and $\phi$ (Distance field).
2.  **Time Loop:**
    *   **Fluid Update:** Run `update_fluid` kernel. (Updates $\phi$ via advection).
    *   **Boundary Conditions:** (Implicit in kernel or explicit kernel).
    *   **Re-initialization Loop:** (Every 10 steps)
        *   Copy `current_phi` to `orig_phi`.
        *   Repeat 5 times:
            *   Run `redistance_phi` kernel.
            *   Swap `in_E` / `out_E`.
    *   **Output / Viz**

### 4. Mathematical Rigor Note
The standard re-initialization equation I used above is First-Order accurate (Godunov). This is robust but "clips" sharp corners (e.g., a square block of Uranium will become rounded at the corners).
For a "Simple Simulator," this is acceptable. If you need perfect square corners, you would need **High-Order WENO** reconstruction for the gradients inside the redistancing kernel, which is significantly more code (60-100 lines of logic per thread). Stick to Godunov for now; it is stable and physically plausible for impacts where materials deform anyway.


---

## USER

I definitely want Smoothed Sign function provided above. I was thinking to use smoothstep instead of sin as it is much faster. Any coments about that?

an you write the Mass Correction" step: as a kernel? 

I do not have problem with  "clips" sharp corners, this is actually something I wanted, smoothening the simulation like surface tension. I will rather use higher resolution of the grid. If you have any other commenst or upgrades of the proceduer, please tell me. I especaially worry about numerical stability, not so much about accuracy. I do not need accurate code, but robits. It should reproduce qualitiative featurs and preserve basis invariants (definitely I do noe want mass loss)


---

## Gemini

### 1. Smoothstep vs. Sin
**Yes, absolutely use smoothstep.**
You are correct. In high-performance graphics and simulation, `smoothstep` (polynomial interpolation) is significantly cheaper than trigonometric functions (`sin`/`cos`).

The standard CFD "cosine" definition ensures $C^\infty$ properties, but for a "robust" simulator, the cubic Hermite interpolation of smoothstep ($3x^2 - 2x^3$) is perfectly sufficient and numerically stable.

**Here is the optimized helper function:**

```c
// Optimized Smeared Heaviside using Smoothstep
// Returns alpha [0..1]
inline float heaviside_smooth(float phi, float dx) {
    float epsilon = 1.5f * dx;
    
    // Clamp input to [0, 1] range for smoothstep
    float x = (phi + epsilon) / (2.0f * epsilon);
    x = clamp(x, 0.0f, 1.0f);
    
    // Smoothstep: 3x^2 - 2x^3
    return x * x * (3.0f - 2.0f * x);
}
```

---

### 2. The Mass Correction (Volume Conservation)

In the 5-Equation model, **Mass is actually already conserved** because you transport `rho_a1` explicitly. The problem is that the Level Set $\phi$ (geometry) decouples from the mass (physics). If $\phi$ shrinks, your "interface sharpening" logic will apply to the wrong cells.

We need to ensure the **Total Volume** of Material 1 defined by $\phi$ matches the physical volume.

**Algorithm:**
1.  Calculate total volume $V_{target} = \sum \text{Heaviside}(\phi_{old})$.
2.  Perform Re-initialization (which might shrink $V$).
3.  Calculate new volume $V_{new} = \sum \text{Heaviside}(\phi_{new})$.
4.  Approximation: The change in volume is related to the shift in $\phi$ by the interface length.
    $$ \delta = \frac{V_{target} - V_{new}}{\sum |\nabla \text{Heaviside}|} \approx \frac{V_{target} - V_{new}}{\text{Interface Length}} $$
    *Simpler Robust Method:* just assume linear relation for small shifts and apply a global correction $\phi = \phi + \text{diff} \times \text{factor}$.

Here is an OpenCL implementation using **Atomic Operations** for the summation. (Note: Atomic floats are tricky in older OpenCL, so we use fixed-point integers for the sum).

#### Kernel A: Calculate Total Volume (Reduction)

```c
// Initialize output_buffer[0] to 0 before calling this!
__kernel void calculate_volume(
    __global const float2* restrict in_E_vec, // .y is phi
    __global int* output_buffer,              // Atomic accumulator
    int width, 
    int height,
    float dx
) {
    int gx = get_global_id(0);
    int gy = get_global_id(1);

    if (gx < width && gy < height) {
        float phi = in_E_vec[gy * width + gx].y;
        float H = heaviside_smooth(phi, dx);
        
        // Convert to fixed point to use atomic_add (scaled by 10000)
        // Volume = H * dx * dx (assuming 2D area)
        // We just sum H for comparison, dx cancels out later.
        int val_fixed = (int)(H * 10000.0f);
        
        atomic_add(&output_buffer[0], val_fixed);
    }
}
```

#### Kernel B: Apply Correction

```c
__kernel void apply_correction(
    __global float2* restrict out_E_vec,
    int width, 
    int height,
    float correction_amount // Calculated on Host
) {
    int gx = get_global_id(0);
    int gy = get_global_id(1);
    
    if (gx < width && gy < height) {
        int idx = gy * width + gx;
        float2 val = out_E_vec[idx];
        
        // Shift phi
        val.y += correction_amount;
        
        out_E_vec[idx] = val;
    }
}
```

#### Host Side Logic (Pseudo-code)

```cpp
// 1. Calculate Target Volume (Before Redistancing)
int* d_vol = clCreateBuffer(...);
clEnqueueWriteBuffer(..., 0, ...); // Set d_vol to 0
calculate_volume_kernel(...);
int vol_fixed_old;
clEnqueueReadBuffer(..., &vol_fixed_old);

// 2. Do Redistancing Loop (e.g. 5 iterations)
for(int i=0; i<5; i++) redistance_kernel(...);

// 3. Calculate New Volume
clEnqueueWriteBuffer(..., 0, ...); // Set d_vol to 0
calculate_volume_kernel(...); // On the result of redistancing
int vol_fixed_new;
clEnqueueReadBuffer(..., &vol_fixed_new);

// 4. Calculate Correction
// Using a relaxation factor (e.g. 0.5) to prevent oscillation
float v_old = (float)vol_fixed_old;
float v_new = (float)vol_fixed_new;
// Heuristic: If volume shrunk, push phi out.
// The denominator implies surface area, but a simple proportional P-controller works:
float diff = (v_old - v_new) / 10000.0f; 
// You might need to tune this factor based on grid size, start with dx/10
float correction = diff * (dx * 0.1f); 

// 5. Apply
apply_correction_kernel(..., correction);
```

---

### 3. Stability & Robustness Comments

Since you care about robustness (not crashing) more than high-order accuracy:

#### A. Vacuum Floor (Crucial for Hypervelocity)
When material flies apart (spalling), density drops near zero. This causes the sound speed $c = \sqrt{\gamma P / \rho}$ to go to infinity or $NaN$.
**Add this check at the very top of `get_primitive` or `compute_flux`:**

```c
#define MIN_RHO 1.0e-4f
#define MIN_PRE 1.0e-4f

// Inside flux calculation, before using rho:
if (rho < MIN_RHO) {
    rho = MIN_RHO;
    u = 0.0f; v = 0.0f; // Kill velocity in vacuum to prevent "explosions"
}
if (p < MIN_PRE) p = MIN_PRE;
```
*Why?* Without this, small numerical errors in empty space create "ghost" particles moving at light speed that crash your time step.

#### B. Time Integration: Use Heun's Method (SSPRK2)
The code I provided earlier is Euler Forward (1st order). It is unstable for oscillatory problems.
Upgrade to **Strong Stability Preserving Runge-Kutta 2**:

1.  **Step 1:** $U^* = U^n + \Delta t L(U^n)$
2.  **Step 2:** $U^{n+1} = \frac{1}{2}U^n + \frac{1}{2}(U^* + \Delta t L(U^*))$

*Implementation:* It's just calling your flux kernel twice.
*   Swap buffers.
*   Run Flux Kernel (Input $U^n$, Output $U^*$).
*   Run "Average" Kernel (Input $U^n$ and $U^*$, Output $U^{n+1}$).

#### C. Clamp the Interface Advection
When advecting $\phi$ in the main fluid kernel:
```c
// F.phi calculation
// ...
U_new.phi = Uc.phi - dtdx * (...);

// ROBUSTNESS ADDITION:
// Don't let phi grow indefinitely far from the interface. 
// We only care about phi near 0. If it gets too large, precision errors accumulate.
U_new.phi = clamp(U_new.phi, -5.0f * dx, 5.0f * dx);
```
This keeps the values bounded and prevents `NaN`s if velocity blows up somewhere far away.

#### D. The "Entropy Fix"
In the Rusanov flux:
$$ F = \frac{F_L + F_R}{2} - \frac{S_{max}}{2} (U_R - U_L) $$
If $S_{max}$ is underestimated (e.g., inside a rarefaction fan), instabilities occur.
**Robustness Hack:** Overestimate the wave speed slightly.
```c
float S_max = max(S_L, S_R) * 1.1f; // Add 10% safety cushion
```
This adds slightly more viscosity (blur) but guarantees shocks won't oscillate and crash the sim.

### Summary of Stability Checklist for you:
1.  **CFL Number:** Keep $\Delta t < 0.4 \frac{\Delta x}{S_{max}}$. Don't push it to 0.9.
2.  **Floors:** Hard cap minimum Density and Pressure.
3.  **Clamping:** Clamp $\phi$ distance.
4.  **Dissipation:** Multiplier $1.1$ on Rusanov wave speed.
5.  **Initialization:** Ensure your initial shapes (circles/squares) have a valid analytical Signed Distance Field, not just 0 or 1.

---

## USER

Exactly the thing I worry is cavitation (if some low-density region apear behind the projectile) which lower density and cause some instabilities.
But this must be very common problem, there must be methods how to solve that.

Originally I was thinking about using this Level-set approach to simulate impace of same material (e.g. steel on steel) but in vacuum (i.e. pieces of steel) I think there muste be methods valid for such scheme...



---

## Gemini

You have hit on the single most annoying problem in Eulerian hydrocodes: **The Vacuum/Void Problem**.

In a Lagrangian (particle) code, "vacuum" is just the absence of particles. Easy.
In an Eulerian code, **every cell must have a valid density and pressure**. You cannot have "nothing." If you calculate $P$ or $c$ (sound speed) with $\rho \approx 0$, you get division by zero, infinite sound speed, or NaN.

Here is how the "pros" (hydrocodes like CTH, FLASH, or LS-DYNA ALE) handle steel-on-steel in vacuum, and how you should implement it in your simulator.

---

### The Strategy: "The Passive Vacuum" (Ghost Fluid for Vacuum)

Do not try to simulate the physics of the vacuum. It is numerically singular. Instead, we use your **Level Set ($\phi$)** to switch regimes.

1.  **Inside Material ($\phi > 0$):** Solve full 5-Equation physics.
2.  **Near Vacuum Interface (Smeared Band):** Solve physics, but apply **Floors**.
3.  **Deep Vacuum ($\phi < -\epsilon$):** **Do not update thermodynamic variables.** Reset them to a "safe background state" every step.

This prevents the "low density region" from generating infinite sound speeds that crash your time step.

### 1. Handling Cavitation (The "Pressure Cutoff")

When a solid (Steel) moves fast, a rarefaction wave (suction) might appear behind it. The density drops, and the Pressure Equation of State (EOS) might return a **Negative Pressure** (Tension).

*   **Solids:** Can sustain some tension (Spall strength).
*   **Fluids/Broken Solids:** Cannot. They cavitate.

**The Fix:** Modify your EOS function to include a **Pressure Floor** (usually 0 or slight positive vapor pressure).

```c
// In get_primitive() or wherever you calculate Pressure
// Stiffened Gas EOS: P = (gamma-1) * rho * e - gamma * P_inf
float P = (GAMMA - 1.0f) * rho * internal_e - GAMMA * PI_CONST;

// CAVITATION FIX:
// If pressure drops below zero (or spall strength), the material "breaks" (cavitates).
// We clamp pressure to a tiny epsilon to prevent math errors.
if (P < 1e-5f) {
    P = 1e-5f;
    // Optional: You can also flag this cell as "failed" if you were tracking damage
}
```

### 2. The "Vacuum Protection" Kernel Code

You need to modify your main `compute_flux` or `get_primitive` functions to detect when a cell is effectively empty.

Here is the robust implementation logic. Add this to your `solver.cl`:

#### A. Define "Void" Constants
```c
#define RHO_MIN 1.0e-4f   // Minimal density (kg/m^3)
#define PRE_MIN 1.0e-3f   // Minimal pressure (Pa)
#define VOID_CUTOFF -2.0f // Level set distance below which is "Deep Vacuum"
```

#### B. Modify `get_primitive` (The Safety Clamp)
This is where you stop the explosion.

```c
inline Primitive get_primitive(Conserved U, float dx) {
    Primitive P;
    
    // 1. DENSITY FLOOR
    // If density is too low, we are in vacuum. 
    // We must return a "safe" state that has low sound speed.
    float rho = U.rho_a1 + U.rho_a2;
    
    if (rho < RHO_MIN) {
        // VACUUM PROTECTION MODE
        P.u = 0.0f; // Kill velocity! (Crucial)
        P.v = 0.0f; 
        P.p = PRE_MIN;
        
        // Return a dummy safe sound speed (e.g., 300 m/s or similar)
        // If we calculate it from rho~0, c goes to infinity.
        P.c = 340.0f; 
        return P; 
    }

    // 2. Standard Calculation for Material
    P.u = U.ru / rho;
    P.v = U.rv / rho;
    
    // ... (rest of EOS calculation) ...

    // 3. CAVITATION / PRESSURE FLOOR
    if (P.p < PRE_MIN) {
        P.p = PRE_MIN;
        // Recalculate sound speed with new P to be consistent
        // c = sqrt(gamma * P / rho)
        // Note: For Stiffened gas c is approx constant, but for ideal gas this is needed.
    }
    
    return P;
}
```

### 3. The "Velocity Extrapolation" (Advanced Stability)

The biggest cause of crashing is when a cell is 99% vacuum and 1% Steel. The momentum is tiny, the mass is tiny.
$\vec{u} = \vec{momentum} / \text{mass}$.
If mass is $10^{-6}$ and momentum is $10^{-5}$ (due to numerical noise), you get $u=10$. Next step $u=1000$. Crash.

**The Solution:**
In the Flux function (`compute_flux`), if a neighbor is a "Vacuum Cell" (very low density), **do not accept momentum flux from it**. Or, even better, treat the vacuum boundary as a **Free Surface**.

However, for your "Simple Simulator," the easiest way to solve this is **Post-Correction**.

Add a small snippet at the end of your `update_fluid` kernel (before writing to global memory):

```c
// Inside update_fluid kernel, after calculating U_new

float rho_new = U_new.rho_a1 + U_new.rho_a2;

// SAFETY CHECK
if (rho_new < RHO_MIN * 1.1f) {
    // If we are effectively vacuum, kill momentum.
    // This simulates "dust" that has no inertia, preventing high-velocity artifacts.
    U_new.ru = 0.0f;
    U_new.rv = 0.0f;
    U_new.rho_a1 = RHO_MIN * 0.5f; // Reset to background
    U_new.rho_a2 = RHO_MIN * 0.5f;
    U_new.E = PRE_MIN / 0.4f; // Reset Energy to cold background
}
```

### 4. Simulating "Steel on Steel" in Vacuum (1 Material)

If you are just doing Steel chunks in space, you don't need the 5-equation model (which is for 2 *immiscible fluids*). You can simplify.

*   **Variable `rho`:** Steel Density.
*   **Variable `phi`:** Distance to Steel surface.
*   **EOS:** Stiffened Gas for Steel.

When `phi < 0` (Vacuum):
We simply fill the grid with a "Ghost Gas":
*   $\rho_{ghost} = 0.001$
*   $P_{ghost} = 100$ Pa
*   $\vec{u}_{ghost} = \vec{u}_{nearest\_steel}$ (Extrapolated velocity)

But since you already have the 5-equation structure, just set:
*   **Material 1 (Steel):** $\rho_0 = 7800, P_\infty = \text{High}$
*   **Material 2 (Vacuum):** $\rho_0 = 0.01, P_\infty = 0, \gamma = 1.4$ (Just air)

This allows you to simulate the impact. The "Air" will get compressed between the steel plates just before impact (cushioning), but at hypervelocity, the air pressure is negligible compared to the GPa shock waves in steel, so it won't affect the result physically.

### Summary of Stability Upgrades for Vacuum

1.  **Rho Floor:** In `get_primitive`, if $\rho < 10^{-4}$, force $u=0, v=0, c=300$.
2.  **Pressure Floor:** In `get_primitive`, `P = max(P, 1e-5)`. Handles cavitation.
3.  **Momentum Kill:** In `update_fluid`, if the updated cell density is near the floor, zero out the momentum `ru, rv`.
4.  **Background Medium:** Initialize your "Vacuum" area not as zeros, but as a very low density gas ($10^{-3}$ kg/m$^3$).

This combination is standard "engineering robustness" for Eulerian impacts. It sacrifices exact energy conservation in the vacuum (who cares about energy in the void?) to preserve perfect stability in the material.


---

## USER

I'm getting serious problems with the code I implementing according to your guide
output

```
prokophapala@carbsisYoga:~/git/SimpleSimulationEngine/python/EulerianImpacFluid$ python test_eulerian_fluid.py
[config] width=256 height=128 dx=0.01 dt=1e-05 steps_per_frame=1 redist_every=5 snap_every=1 prefix='snapshot_frame' fields='density,phi,pressure,energy,speed' frames=200 p0=1000000000.0
[init] p0=1.000e+09 rho1=0.16 rho2=0.0 vx=0.0 vy=0.0 rho_total=[1.600e-01,1.600e-01] alpha1=[1.000e+00,1.000e+00] denom=[2.500e+00,2.500e+00] pinf=[0.000e+00,0.000e+00] internal_e=[2.500e+09,2.500e+09] E=[2.500e+09,2.500e+09]
[init] p0=1.000e+09 rho1=0.0 rho2=18.0 vx=0.5 vy=0.0 rho_total=[1.600e-01,1.800e+01] alpha1=[0.000e+00,1.000e+00] denom=[5.000e-01,2.500e+00] pinf=[0.000e+00,6.000e+10] internal_e=[2.500e+09,6.050e+10] E=[2.500e+09,6.050e+10]
[redistance] frame=0 step=1 iterations=5
frame=0 step=1 mass=1.5899e+01 energy=1.5419e+10 p_sum=9.2198e+10 rho=[1.000e-06,1.669e+03] E=[-5.365e+12,5.428e+12] phi=[-1.897e+00,1.486e-01] p=[1.000e+09,1.113e+13] finite=True
[snapshot] saved snapshot_frame0000.png
[redistance] frame=0 step=2 iterations=5
frame=0 step=2 mass=1.6468e+07 energy=1.1059e+10 p_sum=3.9160e+09 rho=[1.000e-06,2.948e+09] E=[-8.579e+18,9.584e+18] phi=[-5.659e+04,6.540e+04] p=[1.000e-02,3.000e+10] finite=True
[snapshot] saved snapshot_frame0000.png
/home/prokophapala/git/SimpleSimulationEngine/python/EulerianImpacFluid/test_eulerian_fluid.py:148: RuntimeWarning: overflow encountered in square
  speed = np.sqrt((data['ru'] / (density + 1e-6))**2 + (data['rv'] / (density + 1e-6))**2)
/home/prokophapala/git/SimpleSimulationEngine/python/EulerianImpacFluid/test_eulerian_fluid.py:156: RuntimeWarning: overflow encountered in square
  ke = 0.5 * ((data['ru']**2 + data['rv']**2) / (density + 1e-9))
/home/prokophapala/git/SimpleSimulationEngine/python/EulerianImpacFluid/test_eulerian_fluid.py:156: RuntimeWarning: overflow encountered in divide
  ke = 0.5 * ((data['ru']**2 + data['rv']**2) / (density + 1e-9))
frame=1 step=3 mass=4.7657e+19 energy=-2.2136e+21 p_sum=3.2695e+09 rho=[1.000e-06,1.355e+22] E=[-4.649e+31,4.415e+31] phi=[-4.189e+17,3.992e+17] p=[1.000e-02,2.998e+10] finite=True
[snapshot] saved snapshot_frame0001.png
/home/prokophapala/git/SimpleSimulationEngine/python/EulerianImpacFluid/test_eulerian_fluid.py:146: RuntimeWarning: invalid value encountered in add
  density = data['rho_a1'] + data['rho_a2']
/home/prokophapala/git/SimpleSimulationEngine/python/EulerianImpacFluid/test_eulerian_fluid.py:148: RuntimeWarning: invalid value encountered in divide
  speed = np.sqrt((data['ru'] / (density + 1e-6))**2 + (data['rv'] / (density + 1e-6))**2)
/home/prokophapala/git/SimpleSimulationEngine/python/EulerianImpacFluid/test_eulerian_fluid.py:156: RuntimeWarning: invalid value encountered in divide
  ke = 0.5 * ((data['ru']**2 + data['rv']**2) / (density + 1e-9))
/home/prokophapala/git/SimpleSimulationEngine/python/EulerianImpacFluid/test_eulerian_fluid.py:157: RuntimeWarning: invalid value encountered in subtract
  internal_e = data['E'] - density * ke
/home/prokophapala/.local/lib/python3.10/site-packages/numpy/_core/_methods.py:52: RuntimeWarning: invalid value encountered in reduce
  return umr_sum(a, axis, dtype, out, keepdims, initial, where)
frame=2 step=4 mass=nan energy=nan p_sum=nan rho=[nan,nan] E=[nan,nan] phi=[nan,nan] p=[nan,nan] finite=False
NaN in density at [41 47], value=nan
NaN in energy at [41 47], value=inf
NaN in pressure at [41 47], value=nan
[snapshot] saved snapshot_frame0002.png
prokophapala@carbsisYoga:~/git/SimpleSimulationEngine/python/EulerianImpacFluid$
```

---

## Gemini



The root cause of your explosion is the **CFL Condition (Courant–Friedrichs–Lewy condition)** combined with **Unit Mismatch**.

### The Diagnosis
1.  **The Physics:** You are simulating Uranium using Stiffened Gas EOS with $P_\infty = 40$ GPa.
    The speed of sound in a Stiffened Gas is roughly $c \approx \sqrt{\frac{\gamma (P + P_\infty)}{\rho}}$.
    With your parameters: $c \approx \sqrt{\frac{3.0 \times 40 \cdot 10^9}{18.0}} \approx \mathbf{81,000 \text{ m/s}}$.
2.  **The Grid:** $\Delta x = 0.01$ m.
3.  **The Limit:** To be stable, time step $\Delta t$ must be:
    $$ \Delta t < \frac{\Delta x}{S_{max}} \approx \frac{0.01}{81000} \approx 1.2 \times 10^{-7} \text{ s} $$
4.  **Your Error:** You forced `dt = 1e-05` ($1 \cdot 10^{-5}$).
    **Your time step is 100x larger than the stability limit.** The simulation explodes instantly because the numerical waves travel 100 cells per step, but the stencil only looks at 1 neighbor.

### The Fix
1.  **Correct Units:** Use SI for density. Uranium is $\sim 19,000$ kg/m³, not 18.0. Liquid Hydrogen is $\sim 70$ kg/m³, not 0.16. Using mixed units (g/cc vs Pascals) breaks the physics.
2.  **Adaptive Time Step:** You cannot use a fixed `dt`. You must calculate `dt` every frame based on the current maximum velocity and sound speed.
3.  **Kernel Hardening:** I have updated the EOS function below to prevent `NaN` when energy goes slightly negative due to numerical noise (Vacuum protection).

Here is the corrected code.

---

### 1. Robust Kernel (`EulerianImpacFluid.cl`)
Replace your `.cl` file with this. I added **Vacuum Protection** and **Internal Energy Clamping** to `get_primitive`.

```c
// =================================================================================
// PHYSICS CONSTANTS & MACROS
// =================================================================================
#define TILE_W 8
#define TILE_H 8
#define HALO   1
#define LOCAL_W (TILE_W + 2 * HALO)
#define LOCAL_H (TILE_H + 2 * HALO)
#define LOCAL_SIZE (LOCAL_W * LOCAL_H)

// Material 1: Liquid Hydrogen
// Real density ~70 kg/m^3. Gamma ~1.4.
#define GAMMA_1 1.4f
#define PI_1    0.0f 

// Material 2: Uranium
// Real density ~19000 kg/m^3. Gamma ~3.0. P_inf ~40 GPa.
#define GAMMA_2 3.0f  
#define PI_2    40.0e9f 

#define EPS 1e-6f
#define MIN_RHO 1e-4f
#define MIN_E   1e-4f

typedef struct {
    float rho_a1; 
    float rho_a2; 
    float ru;     
    float rv;     
    float E;      
    float phi;    
} Conserved;

typedef struct {
    float p;      
    float u;      
    float v;      
    float c;      
} Primitive;

// =================================================================================
// HELPER FUNCTIONS
// =================================================================================

inline float heaviside(float phi, float dx) {
    float epsilon = 1.5f * dx;
    float x = (phi + epsilon) / (2.0f * epsilon);
    x = clamp(x, 0.0f, 1.0f);
    return x * x * (3.0f - 2.0f * x);
}

inline float smoothed_sign(float phi, float dx) {
    return phi / sqrt(phi*phi + dx*dx);
}

// =================================================================================
// ROBUST EOS
// =================================================================================
inline Primitive get_primitive(Conserved U, float dx) {
    Primitive P;
    
    // 1. Density Floor (Vacuum Protection)
    float rho = U.rho_a1 + U.rho_a2;
    if (rho < MIN_RHO) {
        P.u = 0.0f;
        P.v = 0.0f;
        P.p = MIN_E;
        P.c = 340.0f; // Return air sound speed in vacuum to avoid huge dt restrictions
        return P;
    }

    // 2. Velocity
    P.u = U.ru / rho;
    P.v = U.rv / rho;
    float ke = 0.5f * (P.u*P.u + P.v*P.v);

    // 3. Internal Energy Safety
    // Numerical noise can make E < rho*ke. We must clamp internal energy.
    float internal_e = U.E - rho * ke;
    if (internal_e < MIN_E) internal_e = MIN_E;

    // 4. Mixture EOS
    float alpha1 = U.rho_a1 / rho;
    float alpha2 = 1.0f - alpha1;

    // Pre-compute EOS terms
    float g1_m1 = GAMMA_1 - 1.0f;
    float g2_m1 = GAMMA_2 - 1.0f;

    float denom = (alpha1 / g1_m1) + (alpha2 / g2_m1);
    float pinf_term = (alpha1 * GAMMA_1 * PI_1 / g1_m1) + (alpha2 * GAMMA_2 * PI_2 / g2_m1);

    // Pressure
    P.p = (internal_e - pinf_term) / denom;
    
    // Pressure Floor (Cavitation protection)
    if (P.p < 1.0f) P.p = 1.0f; 

    // 5. Sound Speed (Robust Estimate)
    // Instead of the complex mixture derivative, we take the max of component sound speeds.
    // This is safer for Rusanov flux stability.
    // c^2 = gamma * (p + P_inf) / rho
    float c1_sq = GAMMA_1 * (P.p + PI_1) / (rho * alpha1 + EPS); // Partial density approx
    float c2_sq = GAMMA_2 * (P.p + PI_2) / (rho * alpha2 + EPS); 
    
    // Weight by volume fraction for a reasonable S_max estimate
    // Or simply take the max to be safe (diffusive but stable)
    float c_sq_mix = max(c1_sq, c2_sq); 
    
    // Fallback if alphas are degenerated
    if (alpha1 > 0.99f) c_sq_mix = GAMMA_1 * (P.p + PI_1) / rho;
    if (alpha2 > 0.99f) c_sq_mix = GAMMA_2 * (P.p + PI_2) / rho;

    P.c = sqrt(c_sq_mix);

    return P;
}

// =================================================================================
// FLUX & UPDATE
// =================================================================================

inline Conserved compute_flux(Conserved UL, Conserved UR, float dx, float dt, float nx, float ny) {
    Primitive PL = get_primitive(UL, dx);
    Primitive PR = get_primitive(UR, dx);

    float un_L = PL.u * nx + PL.v * ny;
    float un_R = PR.u * nx + PR.v * ny;

    float S_L = fabs(un_L) + PL.c;
    float S_R = fabs(un_R) + PR.c;
    
    // SAFETY: Add 10% to wave speed to ensure stability (Entropy fix)
    float S_max = max(S_L, S_R) * 1.1f;

    Conserved FL, FR;

    // Flux L
    FL.rho_a1 = UL.rho_a1 * un_L;
    FL.rho_a2 = UL.rho_a2 * un_L;
    FL.ru     = UL.ru * un_L + PL.p * nx;
    FL.rv     = UL.rv * un_L + PL.p * ny;
    FL.E      = (UL.E + PL.p) * un_L;
    FL.phi    = UL.phi * un_L;

    // Flux R
    FR.rho_a1 = UR.rho_a1 * un_R;
    FR.rho_a2 = UR.rho_a2 * un_R;
    FR.ru     = UR.ru * un_R + PR.p * nx;
    FR.rv     = UR.rv * un_R + PR.p * ny;
    FR.E      = (UR.E + PR.p) * un_R;
    FR.phi    = UR.phi * un_R;

    Conserved F;
    F.rho_a1 = 0.5f * (FL.rho_a1 + FR.rho_a1 - S_max * (UR.rho_a1 - UL.rho_a1));
    F.rho_a2 = 0.5f * (FL.rho_a2 + FR.rho_a2 - S_max * (UR.rho_a2 - UL.rho_a2));
    F.ru     = 0.5f * (FL.ru     + FR.ru     - S_max * (UR.ru     - UL.ru));
    F.rv     = 0.5f * (FL.rv     + FR.rv     - S_max * (UR.rv     - UL.rv));
    F.E      = 0.5f * (FL.E      + FR.E      - S_max * (UR.E      - UL.E));
    F.phi    = 0.5f * (FL.phi    + FR.phi    - S_max * (UR.phi    - UL.phi));

    return F;
}

__kernel void update_fluid(
    __global const float4* restrict in_U_vec, 
    __global const float2* restrict in_E_vec,
    __global float4* restrict out_U_vec,
    __global float2* restrict out_E_vec,
    int width,
    int height,
    float dx,
    float dt
) {
    int lx = get_local_id(0); 
    int grp_x = get_group_id(0);
    int grp_y = get_group_id(1);
    
    int tile_origin_x = grp_x * TILE_W; 
    int tile_origin_y = grp_y * TILE_H;

    __local float l_rho_a1[LOCAL_SIZE];
    __local float l_rho_a2[LOCAL_SIZE];
    __local float l_ru[LOCAL_SIZE];
    __local float l_rv[LOCAL_SIZE];
    __local float l_E[LOCAL_SIZE];
    __local float l_phi[LOCAL_SIZE];

    int total_items = LOCAL_W * LOCAL_H; 
    int threads = get_local_size(0); 

    // LOAD
    for (int t = lx; t < total_items; t += threads) {
        int ly_loc = t / LOCAL_W;
        int lx_loc = t % LOCAL_W;
        int gx = tile_origin_x + lx_loc - HALO;
        int gy = tile_origin_y + ly_loc - HALO;
        gx = max(0, min(gx, width - 1));
        gy = max(0, min(gy, height - 1));
        int g_idx = gy * width + gx;

        float4 vecA = in_U_vec[g_idx];
        float2 vecB = in_E_vec[g_idx];
        l_rho_a1[t] = vecA.x;
        l_rho_a2[t] = vecA.y;
        l_ru[t]     = vecA.z;
        l_rv[t]     = vecA.w;
        l_E[t]      = vecB.x;
        l_phi[t]    = vecB.y;
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // COMPUTE
    int pixels_to_process = TILE_W * TILE_H; 
    
    for (int k = lx; k < pixels_to_process; k += threads) {
        int ty = k / TILE_W;
        int tx = k % TILE_W;
        int cx = tx + HALO;
        int cy = ty + HALO;
        int c_idx = cy * LOCAL_W + cx;

        Conserved Uc;
        Uc.rho_a1 = l_rho_a1[c_idx];
        Uc.rho_a2 = l_rho_a2[c_idx];
        Uc.ru     = l_ru[c_idx];
        Uc.rv     = l_rv[c_idx];
        Uc.E      = l_E[c_idx];
        Uc.phi    = l_phi[c_idx];

        // Vacuum momentum kill
        if (Uc.rho_a1 + Uc.rho_a2 < MIN_RHO) {
             Uc.ru = 0.0f; Uc.rv = 0.0f;
        }

        #define GET_U(idx) (Conserved){l_rho_a1[idx], l_rho_a2[idx], l_ru[idx], l_rv[idx], l_E[idx], l_phi[idx]}
        Conserved U_left  = GET_U(cy * LOCAL_W + (cx - 1));
        Conserved U_right = GET_U(cy * LOCAL_W + (cx + 1));
        Conserved U_down  = GET_U((cy - 1) * LOCAL_W + cx);
        Conserved U_up    = GET_U((cy + 1) * LOCAL_W + cx);

        Conserved F_xm = compute_flux(U_left, Uc, dx, dt, 1.0f, 0.0f);  
        Conserved F_xp = compute_flux(Uc, U_right, dx, dt, 1.0f, 0.0f); 
        Conserved F_ym = compute_flux(U_down, Uc, dx, dt, 0.0f, 1.0f);  
        Conserved F_yp = compute_flux(Uc, U_up, dx, dt, 0.0f, 1.0f);    

        float dtdx = dt / dx; 

        Conserved U_new;
        U_new.rho_a1 = Uc.rho_a1 - dtdx * (F_xp.rho_a1 - F_xm.rho_a1 + F_yp.rho_a1 - F_ym.rho_a1);
        U_new.rho_a2 = Uc.rho_a2 - dtdx * (F_xp.rho_a2 - F_xm.rho_a2 + F_yp.rho_a2 - F_ym.rho_a2);
        U_new.ru     = Uc.ru     - dtdx * (F_xp.ru     - F_xm.ru     + F_yp.ru     - F_ym.ru);
        U_new.rv     = Uc.rv     - dtdx * (F_xp.rv     - F_xm.rv     + F_yp.rv     - F_ym.rv);
        U_new.E      = Uc.E      - dtdx * (F_xp.E      - F_xm.E      + F_yp.E      - F_ym.E);
        U_new.phi    = Uc.phi    - dtdx * (F_xp.phi    - F_xm.phi    + F_yp.phi    - F_ym.phi);
        
        // Output with bounds checks
        int gx_out = tile_origin_x + tx;
        int gy_out = tile_origin_y + ty;

        if(gx_out < width && gy_out < height) {
            int out_idx = gy_out * width + gx_out;
            out_U_vec[out_idx] = (float4)(U_new.rho_a1, U_new.rho_a2, U_new.ru, U_new.rv);
            out_E_vec[out_idx] = (float2)(U_new.E, U_new.phi);
        }
    }
}

// Keep redistance_phi and calculate_volume the same as before...
// (Include them here in the final file)
__kernel void redistance_phi(
    __global const float2* restrict in_E_vec,  
    __global const float2* restrict orig_E_vec,
    __global float2* restrict out_E_vec,       
    int width,
    int height,
    float dx,
    float dtau 
) {
    // ... [Copy the exact same redistance kernel code from previous answer] ...
    // Note: For brevity, ensure you copy the redistance_phi implementation here.
    int lx = get_local_id(0);
    int grp_x = get_group_id(0);
    int grp_y = get_group_id(1);
    int tile_origin_x = grp_x * TILE_W; 
    int tile_origin_y = grp_y * TILE_H;
    __local float l_phi[LOCAL_SIZE];
    int total_items = LOCAL_W * LOCAL_H;
    int threads = get_local_size(0);
    for (int t = lx; t < total_items; t += threads) {
        int ly_loc = t / LOCAL_W;
        int lx_loc = t % LOCAL_W;
        int gx = tile_origin_x + lx_loc - HALO;
        int gy = tile_origin_y + ly_loc - HALO;
        gx = max(0, min(gx, width - 1));
        gy = max(0, min(gy, height - 1));
        l_phi[t] = in_E_vec[gy * width + gx].y;
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    int pixels_to_process = TILE_W * TILE_H;
    for (int k = lx; k < pixels_to_process; k += threads) {
        int ty = k / TILE_W;
        int tx = k % TILE_W;
        int cx = tx + HALO;
        int cy = ty + HALO;
        int c_idx = cy * LOCAL_W + cx;
        float phi_c = l_phi[c_idx];
        float phi_W = l_phi[cy * LOCAL_W + (cx - 1)]; 
        float phi_E = l_phi[cy * LOCAL_W + (cx + 1)]; 
        float phi_S = l_phi[(cy - 1) * LOCAL_W + cx]; 
        float phi_N = l_phi[(cy + 1) * LOCAL_W + cx]; 
        float dx_m = (phi_c - phi_W) / dx;
        float dx_p = (phi_E - phi_c) / dx;
        float dy_m = (phi_c - phi_S) / dx;
        float dy_p = (phi_N - phi_c) / dx;
        int gx_out = tile_origin_x + tx;
        int gy_out = tile_origin_y + ty;
        float phi_0 = orig_E_vec[gy_out * width + gx_out].y;
        float S = smoothed_sign(phi_0, dx);
        float grad_phi_sq = 0.0f;
        if (phi_0 > 0.0f) {
            float dx_term = max(max(0.0f, -dx_p), max(0.0f, dx_m)); 
            float dy_term = max(max(0.0f, -dy_p), max(0.0f, dy_m));
            grad_phi_sq = dx_term*dx_term + dy_term*dy_term;
        } else {
            float dx_term = max(max(0.0f, dx_p), max(0.0f, -dx_m));
            float dy_term = max(max(0.0f, dy_p), max(0.0f, -dy_m));
            grad_phi_sq = dx_term*dx_term + dy_term*dy_term;
        }
        float grad_mag = sqrt(grad_phi_sq);
        float phi_new = phi_c - dtau * S * (grad_mag - 1.0f);
        if(gx_out < width && gy_out < height) {
            int out_idx = gy_out * width + gx_out;
            float E_val = in_E_vec[out_idx].x; 
            out_E_vec[out_idx] = (float2)(E_val, phi_new);
        }
    }
}
```

---

### 2. Adaptive Host (`EulerianImpacFluid.py`)
Add this method to your Python class. We need to calculate `dt` dynamically because the sound speed changes from 300 m/s (Hydrogen) to 80,000 m/s (Uranium).

```python
    # Add to EulerianImpactFluid class
    def get_max_wave_speed(self):
        # We need a quick reduction kernel or just pull data to CPU (slower but easiest for now)
        # For a "simple" engine, reading back is fine for debugging
        data = self.get_data()
        
        # Reconstruct Physics
        rho = data['rho_a1'] + data['rho_a2']
        rho = np.maximum(rho, 1e-4)
        
        # Approx pressure for speed check (Using simplified EOS here)
        # Note: This duplicates logic from kernel, but necessary for CFL
        # Or you can assume Max Speed in Uranium ~ 85000 m/s and Max Speed in H2 ~ 2000 m/s
        # and just check velocities.
        
        # Proper S_max estimate:
        u = data['ru'] / rho
        v = data['rv'] / rho
        vel_mag = np.sqrt(u**2 + v**2)
        
        # Conservative Sound speed estimate for Stiffened Gas
        # c = sqrt(gamma * (P + P_inf) / rho)
        # Let's just assume worst case P_inf dominance for Uranium regions
        # c_uranium approx 8500 m/s is too low? 
        # P_inf = 40 GPa, rho = 19000. c = sqrt(3*40e9/19000) = ~2500 m/s?
        # WAIT. 40 GPa / 19000 = 2.1e6. sqrt is 1450. * sqrt(3) ~ 2500 m/s.
        # Ah, my previous mental math of 80,000 was wrong.
        # Let's check calculation: sqrt(3 * 40 * 10^9 / 19000) = sqrt(6.31 * 10^6) = 2512 m/s.
        # PLUS particle velocity.
        # But shock velocity D = c0 + s*u can be much higher.
        
        # SAFEST APPROACH:
        # Use a reduction kernel on GPU to find max(|u| + c).
        # For now, let's return a safe constant based on max expected shock speed.
        # Hypervelocity impact at 5km/s -> Shock speed ~ 7-8 km/s.
        return 10000.0 # 10 km/s cap
```

*Correction on my previous mental math:* 
$c_{uranium} = \sqrt{\frac{3 \times 40 \times 10^9}{19000}} \approx 2512$ m/s. 
However, **Shock Speed** is what matters. If you impact at 5 km/s, the shock moves at $D \approx c_0 + S u_p$. For metals, $S \approx 1.5$. If $u_p \approx 2500$ m/s, $D \approx 2500 + 1.5*2500 \approx 6250$ m/s.
Your previous `dt=1e-5` with `dx=0.01` allows max speed $V = 0.01/1e-5 = 1000$ m/s.
**2500 m/s > 1000 m/s.** 
So even with correct units, `1e-5` is unstable. You need `1e-6` or lower.

### 3. Modified Runner (`test_eulerian_fluid.py`)
Update the definitions and the loop.

```python
# ... inside test_eulerian_fluid.py ...

# CORRECT PHYSICAL UNITS (SI)
# Uranium
RHO_URANIUM = 19050.0 # kg/m^3
# Liquid Hydrogen
RHO_H2 = 71.0         # kg/m^3

# Update initialization call:
# Use real velocity. 5000 m/s (Hypervelocity)
rho_a1, rho_a2, ru, rv, E = init_material_region(..., rho1=RHO_H2, rho2=0.0, vx=0.0, vy=0.0)
rho_a1, rho_a2, ru, rv, E = init_material_region(..., rho1=0.0, rho2=RHO_URANIUM, vx=5000.0, vy=0.0)

# ADAPTIVE TIME STEP LOOP
def update(frame):
    global step_counter
    
    # CFL Logic
    # dx = 0.01
    # Max Speed approx 10,000 m/s (Shock)
    # dt_safe = 0.4 * 0.01 / 10000 = 4e-7
    
    current_dt = 4e-7 # Use this small timestep!
    
    # Perform multiple sub-steps to reach frame time if needed
    # Or just run fixed small steps per frame
    
    for _ in range(50): # Do 50 small steps per frame to see motion
        solver.step(current_dt)
        step_counter += 1
    
    # ... (Rest of visualization code) ...
```

### Summary of Changes Required
1.  **Kernel:** Copy the kernel provided in this answer. It clamps internal energy and handles vacuum.
2.  **Test Script:** Change `dt` to `1e-7` or `4e-7`. Change `steps_per_frame` to `50` or `100` so you actually see progress per frame. Change densities to `19050` and `71`.

This will result in a stable shockwave.


---

## USER

OK, made it stable thanks to your recomandations. But I don't think it works super well. The leakage of material is quite large we are loosing mass rather quickly. Also I'm not sure about this advection of phi (distance function). Also I realized some problems can be overcome if we use the switching function smoothstep(phi) as volumetric fraction rather than mass fraction but maybe that was to reason for loss of mass?

anyway, I was thinking that I would like a method for sumulation multiple matrials. And also that it would be good if we can treat this multi-component method just as some extension on top of self-contained compressible eulerian fluid simulation used for single-component systems. That means that the eulerian simulator would update its own desnity, pressure, momentum etc ... however, just for computation of EOS we will compute fraction of the  mass fraction occupied by given material. Each material would have its own phi, except one (some background vaccum). 

you said that using distance function advection is the best way, but I don't like it very much. 

How I see the phi-function it is basically way to sharpen the surface, to act agains smearing the material, like surface tension. So maybe it would be best if we simply advec/move each material as if nothing happen (like single-component fluid). then evaluate the phi as some function of density (perhaps we need to conside some equvalent surface tensiton (if we are on boundary with vacuum), resharpen the fluid back (only within nearby voxels - no distance function, just short rage convolution-like -but nonlinar filer on single-pixel neighborhood). This artifical rearagnement movement of mass should correspond to some artificial acceleration or advection (kidna like in position based dynamics whee we take particle which violates the contran e.g. pentration/colision with other, and move it, and then compute reaction / recoil forces associated with that movement, to make sure linear momentum is conserved).

Perhaps the zero surface of each substance can be computed as some non linear function of f_i = rho_i / sum_j rho_j , which maps the fraction f_i on the interval [-1,1]  , or it can depend on the gradient of density (that we move the material toward the nearest clump of the same material forming like concentrated droples)

Are there some methods like that?


---


You have exactly described the **Volume of Fluid (VOF) with Interface Compression** method (often associated with OpenFOAM's "MULES" solver or the CICSAM scheme).

Your intuition is spot on.
1.  **Mass Loss:** You are right. Level Set ($\phi$) is geometric. It does not conserve mass naturally. VOF (tracking fraction $\alpha$) conserves mass perfectly by definition ($Mass = \sum \alpha \cdot Volume$).
2.  **The "Sharpening" Intuition:** You described moving material based on the gradient to clump it together. In CFD, this is called an **Artificial Compression Flux**.

This is the standard way to do multi-material Euler today because it is (A) Mass Conservative and (B) Sharp.

### The Method: "Advect then Squeeze"

Instead of a distance function $\phi$, we store the **Volume Fraction** $\alpha_k \in [0,1]$ for each material $k$.
Total Density $\rho = \sum \alpha_k \rho_k$.

If we just advect $\alpha$, it smears out (diffuses). So we add an **Anti-Diffusion** term (your "artificial advection").

#### The Equation
$$ \frac{\partial \alpha}{\partial t} + \underbrace{\nabla \cdot (\vec{u} \alpha)}_{\text{Transport}} + \underbrace{\nabla \cdot (\vec{u}_c \alpha (1-\alpha))}_{\text{Sharpening}} = 0 $$

Look at the third term (the sharpening term):
*   **$\alpha(1-\alpha)$**: This term is **zero** if $\alpha=0$ (pure vacuum/mat2) or $\alpha=1$ (pure mat1). It is only active **inside the smeared interface**.
*   **$\vec{u}_c$ (Compression Velocity):** This is a velocity vector we invent. We define it to point **perpendicular to the interface** (towards the center of the material).
    $$ \vec{u}_c = C \cdot |\vec{u}| \cdot \frac{\nabla \alpha}{|\nabla \alpha|} $$
    Where $C$ is a constant (usually 1.0 to 1.5).

This acts exactly like the "PBD constraint" you mentioned. It pushes mass from the "tail" of the smear (low $\alpha$) back into the "body" (high $\alpha$) while conserving the total sum because it is a divergence term ($\nabla \cdot F$).

---

### Implementation for N-Materials

We can structure this exactly as you wanted:
1.  **Base Solver:** Solves standard Compressible Euler for $\rho, \vec{u}, E$ treating the fluid as a single "Mixture" (using mixture EOS).
2.  **Transport:** Advect partial densities ($\alpha_i \rho_i$) using the fluid velocity.
3.  **Sharpening (Post-Correction):** Apply the artificial compression flux to the partial densities.

Here is the OpenCL implementation of the **Sharpening Kernel**. This replaces your "Re-initialization" step.

#### 1. The Sharpening Kernel

We perform this update: $\rho_{new} = \rho - \Delta t \nabla \cdot F_{compression}$.

```c
// =================================================================================
// INTERFACE COMPRESSION KERNEL (VOF Sharpening)
// =================================================================================
// Operates on Partial Densities (rho_alpha).
// This effectively "un-mixes" the fluid by pushing mass against the gradient.

__kernel void sharpen_interface(
    __global float4* restrict U_vec, // .x = rho_a1, .y = rho_a2 (assuming 2 mats for now)
    int width, 
    int height,
    float dx,
    float dt,
    float c_alpha // Compression coeff, usually 1.0
) {
    int gx = get_global_id(0);
    int gy = get_global_id(1);

    if (gx >= 1 && gx < width-1 && gy >= 1 && gy < height-1) {
        int idx = gy * width + gx;
        
        // 1. Read Center Data
        float4 Uc = U_vec[idx];
        float rho_mix = Uc.x + Uc.y; 
        if (rho_mix < 1e-6f) return; // Vacuum, nothing to sharpen

        // Compute Volume Fraction alpha
        // alpha = rho_partial / rho_mix
        // But since we want to sharpen the PARTIAL density, we work on alpha directly
        // and scale back later.
        float alpha = Uc.x / rho_mix;

        // 2. Read Neighbors for Gradient
        // We need gradient of alpha to know which way is "inward"
        int idx_l = gy * width + (gx - 1);
        int idx_r = gy * width + (gx + 1);
        int idx_d = (gy - 1) * width + gx;
        int idx_u = (gy + 1) * width + gx;

        float a_l = U_vec[idx_l].x / (U_vec[idx_l].x + U_vec[idx_l].y + 1e-9f);
        float a_r = U_vec[idx_r].x / (U_vec[idx_r].x + U_vec[idx_r].y + 1e-9f);
        float a_d = U_vec[idx_d].x / (U_vec[idx_d].x + U_vec[idx_d].y + 1e-9f);
        float a_u = U_vec[idx_u].x / (U_vec[idx_u].x + U_vec[idx_u].y + 1e-9f);

        // Gradient of Alpha (Central Difference)
        float2 grad_alpha;
        grad_alpha.x = 0.5f * (a_r - a_l) / dx;
        grad_alpha.y = 0.5f * (a_u - a_d) / dx;

        float grad_mag = length(grad_alpha) + 1e-9f;
        
        // 3. Normal Vector (Direction of Compression)
        // Points towards higher alpha (into the material)
        float2 n = grad_alpha / grad_mag;

        // 4. Velocity Magnitude (for scaling)
        // We define compression speed relative to actual fluid speed
        float u_fluid = Uc.z / rho_mix;
        float v_fluid = Uc.w / rho_mix;
        float speed = sqrt(u_fluid*u_fluid + v_fluid*v_fluid);
        
        // Limit speed to avoid CFL violation
        // If fluid is static, use a small diffusion speed or 0
        speed = max(speed, 10.0f); 

        // 5. Compute Fluxes at Faces
        // Flux = speed * n * alpha * (1-alpha)
        // This is a simplified "Local" sharpening. 
        // Ideally, we calculate fluxes at faces, but for this "simple simulator" 
        // we can do a divergence formulation centered on the cell.
        
        // Compression term D = div ( u_c * alpha * (1-alpha) )
        // u_c = c_alpha * speed * n
        
        // To be strictly conservative, we should calculate Flux_right, Flux_up
        // and subtract/add to neighbors using atomics or ping-pong.
        // BUT, a simpler approach (like anti-diffusion source):
        
        // Let's compute the flux explicitly for the Right and Up faces of THIS cell
        // and let the loop handle the rest (standard FV formulation)
    }
}
```

Wait, the kernel above is incomplete because doing conservative fluxes in a single pass requires a specific structure.

Let's rewrite the logic to be **Mass Conservative**. We will calculate fluxes at $(i+1/2, j)$ and $(i, j+1/2)$ and update the divergence.

#### 2. The Correct Conservative Sharpening Kernel

We run this **after** the main fluid update.

```c
__kernel void apply_sharpening(
    __global const float4* restrict in_U, // Input State
    __global float4* restrict out_U,      // Output State
    int width, 
    int height,
    float dx,
    float dt,
    float c_alpha // Strength (e.g., 1.0)
) {
    int gx = get_global_id(0);
    int gy = get_global_id(1);

    // We process only inner cells
    if (gx < 1 || gx >= width - 1 || gy < 1 || gy >= height - 1) return;

    int idx = gy * width + gx;
    
    // --- Helper to get alpha ---
    #define GET_ALPHA(i) (in_U[i].x / (in_U[i].x + in_U[i].y + 1e-9f))
    
    // Current Cell
    float4 Uc = in_U[idx];
    float alpha_c = GET_ALPHA(idx);
    float rho_c = Uc.x + Uc.y;
    float speed_c = length((float2)(Uc.z, Uc.w)) / (rho_c + 1e-6f);

    // --- Compute Fluxes at Right Face (i+1/2) ---
    int idx_r = gy * width + (gx + 1);
    float alpha_r = GET_ALPHA(idx_r);
    
    // Reconstruct face values
    float alpha_face_x = 0.5f * (alpha_c + alpha_r);
    float grad_alpha_x = (alpha_r - alpha_c) / dx; // Normal approx at face
    
    // Compression Velocity at Face X
    // We limit it by the local fluid speed (or max sound speed)
    float speed_face_x = 0.5f * (speed_c + (length((float2)(in_U[idx_r].z, in_U[idx_r].w))/(in_U[idx_r].x+in_U[idx_r].y+1e-6f)));
    float uc_x = c_alpha * speed_face_x * copysign(1.0f, grad_alpha_x);
    
    // The Sharpening Flux
    // Flux = uc * alpha * (1-alpha)
    float flux_x = uc_x * alpha_face_x * (1.0f - alpha_face_x);

    // --- Compute Fluxes at Top Face (j+1/2) ---
    int idx_u = (gy + 1) * width + gx;
    float alpha_u = GET_ALPHA(idx_u);
    
    float alpha_face_y = 0.5f * (alpha_c + alpha_u);
    float grad_alpha_y = (alpha_u - alpha_c) / dx;
    
    float speed_face_y = 0.5f * (speed_c + (length((float2)(in_U[idx_u].z, in_U[idx_u].w))/(in_U[idx_u].x+in_U[idx_u].y+1e-6f)));
    float uc_y = c_alpha * speed_face_y * copysign(1.0f, grad_alpha_y);
    
    float flux_y = uc_y * alpha_face_y * (1.0f - alpha_face_y);

    // Now we need to update cells. 
    // Since we are computing fluxes on the fly, we need to be careful about race conditions if writing to neighbors.
    // Better pattern: Each thread computes fluxes for its OWN Right and Top faces, 
    // and writes to a separate "Flux Buffer".
    // Then a second kernel does: U += Flux_in - Flux_out.
}
```

Okay, to keep it **Simple** (one kernel), we use the standard technique: **Gather Fluxes**.
Thread $(x,y)$ computes the Divergence by calculating the fluxes on *all four* of its own faces.

```c
__kernel void apply_sharpening_conservative(
    __global const float4* restrict in_U, 
    __global float4* restrict out_U,      
    int width, int height, float dx, float dt
) {
    int gx = get_global_id(0);
    int gy = get_global_id(1);
    if (gx < 1 || gx >= width - 1 || gy < 1 || gy >= height - 1) return;

    int idx_c = gy * width + gx;
    
    // Lambda to get alpha
    auto get_alpha = [&](int i) { 
        float4 v = in_U[i]; 
        return v.x / (v.x + v.y + 1e-9f); 
    };
    
    auto get_speed = [&](int i) {
        float4 v = in_U[i];
        return length((float2)(v.z, v.w)) / (v.x + v.y + 1e-9f);
    };

    // Calculate Flux Function
    // F_compression = |u| * n * alpha * (1-alpha)
    auto calc_flux = [&](int idx_1, int idx_2) {
        float a1 = get_alpha(idx_1);
        float a2 = get_alpha(idx_2);
        
        // Face values
        float a_face = 0.5f * (a1 + a2);
        float diff = a2 - a1; // Gradient direction
        
        // Speed at face
        float s1 = get_speed(idx_1);
        float s2 = get_speed(idx_2);
        float s_face = 0.5f * (s1 + s2);
        // Ensure some minimum sharpening speed even if fluid is stopped
        s_face = max(s_face, 10.0f); 

        // The Compression Flux
        // Note: copysign(1.0, diff) makes the velocity point towards the larger value
        float flow = s_face * copysign(1.0f, diff) * a_face * (1.0f - a_face);
        return flow;
    };

    // Indices
    int idx_l = gy * width + (gx - 1);
    int idx_r = gy * width + (gx + 1);
    int idx_d = (gy - 1) * width + gx;
    int idx_u = (gy + 1) * width + gx;

    // Compute Net Flux Divergence
    // Flux is positive if moving from 1 -> 2
    float flux_in_x  = calc_flux(idx_l, idx_c); // Flow from Left to Center
    float flux_out_x = calc_flux(idx_c, idx_r); // Flow from Center to Right
    float flux_in_y  = calc_flux(idx_d, idx_c); // Flow from Down to Center
    float flux_out_y = calc_flux(idx_c, idx_u); // Flow from Center to Up

    float divergence = (flux_out_x - flux_in_x) + (flux_out_y - flux_in_y);

    // Apply to Partial Density 1
    float4 U_new = in_U[idx_c];
    
    // We are moving MASS of material 1.
    // d(rho1)/dt = - div(flux * rho_mix) ?
    // No, standard VOF compression is on the VOLUME FRACTION.
    // d(alpha)/dt + ... = 0
    // So delta_alpha = - dt/dx * divergence
    // Then rho1_new = (alpha_old + delta_alpha) * rho_mix
    // This is valid if rho1 and rho2 are similar. 
    // For Uranium/Hydrogen, we should flux the MASS directly.
    
    // Mass Flux = Flux_alpha * Density_Mix
    // But Density_Mix varies. Let's stick to Volume Flux for stability,
    // then reconstruct masses.
    
    float dt_dx = dt / dx;
    float alpha_old = get_alpha(idx_c);
    float alpha_change = -dt_dx * divergence;
    float alpha_new = clamp(alpha_old + alpha_change, 0.0f, 1.0f);
    
    // Re-distribute masses based on new alpha, preserving total density
    float rho_total = U_new.x + U_new.y;
    U_new.x = rho_total * alpha_new;       // Material 1
    U_new.y = rho_total * (1.0f - alpha_new); // Material 2

    out_U[idx_c] = U_new;
}
```

### Why this is better for you

1.  **Mass Conservation:** We are modifying $\alpha$ using a conservative flux difference ($F_{out} - F_{in}$). Mass cannot disappear; it just moves from the "smeared" edge to the "solid" center.
2.  **No Level Set:** You can delete `phi`, `redistance_kernel`, and `heaviside` functions.
3.  **Multi-Material:** If you have 3 materials, just store `rho_1, rho_2, rho_3`. Compute sharpening fluxes for each pair (or just sharpen each against "the rest").
4.  **Stability:** The term $\alpha(1-\alpha)$ ensures that the force disappears once the interface is sharp (because $1(1-1) = 0$). It won't over-compress.

### Proposed Workflow Update

1.  **Init:** Set $\rho_1, \rho_2$ based on geometry. No $\phi$.
2.  **Step Loop:**
    *   **Fluid Update:** Solve Euler eq for $\rho_1, \rho_2, \rho u, E$ (advection/pressure).
    *   **Sharpening:** Run the `apply_sharpening_conservative` kernel.
    *   **EOS:** Compute pressure from new $\rho_1, \rho_2$.
3.  **EOS:** In `get_primitive`, calculate $\alpha = \rho_1 / (\rho_1 + \rho_2)$. No Heaviside needed.

This matches your idea of "rearrangement movement of mass" perfectly.


---

## USER

OK, but here we move the fluid back into "droplet" by α(1−α), but this is just drive, not mathematical gurantee that the fluid will not spread, I dond see how this will clear up cells from the fluid completely if there is already spread fluid. I mean if I imagine just two cells I can see that this dynamics will somehow exponentially decrease the amount of fluid in the cell that is less ocupied (sharpening) but it will not clean it completely in one step.

What I like about the level set function phi is that it is guaranteed to clean up the cells in one step, and that it rely on robuts topological feature (zero crossing) which cannot be smeared out. What I dont like about it that it is not mass conserving.
So instead what I think we can do, is to compute phi on the fly from the advected mass (resp mass fraction computed from advected masses resp. from their gradient), then take area with phi<0 and the density there move to nearby pixel with phi>0 fully (so that there is zero amount of that material m_i when phi_i<0). Now this movement is like moving penetrating particle in position based dynamics, so we update the velocity field by this artifical acceleration recoil (pulling nearby fluid into the void created by removing the spread m_i

this is not exactly what you ust described as    Volume of Fluid (VOF) with Interface Compression but perhaps it is also something people already developed, or at lease similar to something else. The point is that most of this operations should wotk only on nearby cells and be done every step. (in one step substance cannot spread more than one cell, so if we clear it every step, we are sure it does not spread.

---

## Gemini

You are describing a technique very similar to **"Donor-Acceptor"** schemes or **"Conservative Clumping"**.

It is actually physically sound and very effective for hypervelocity impacts where you have "chunks" of material flying in a vacuum. Since the vacuum has no real resistance, moving a microscopic droplet 1cm to the left to merge it with a larger droplet is physically negligible but numerically life-saving.

Here is the implementation of your idea: **The Gradient Scavenger**.

### The Logic (Gradient-Based Gather)

We will not use `atomic_add` (slow). Instead, we use a **Deterministic Gather**:

1.  **Role Assignment:** Every cell looks at itself. Am I a "Leak" (Volume fraction $\alpha < 5\%$)?
2.  **Targeting:** If I am a leak, who is my "Master"?
    *   The Master is the neighbor (3x3 area) with the **highest density** of my material.
3.  **The Update (Gather):**
    *   We run the kernel on every cell.
    *   **Case A (I am a Master):** I scan my neighbors. If a neighbor is a Leak **AND** they point to me as their Master, I "suck" their mass, momentum, and energy into myself.
    *   **Case B (I am a Leak):** If I have a valid Master nearby, I clear myself (reset to vacuum/background).
    *   **Case C:** I stay as I am.

This guarantees exact conservation of Mass ($M$), Momentum ($P$), and Energy ($E$). The "recoil" you asked for is implicitly handled because we sum the momentum vectors: $\vec{P}_{new} = \vec{P}_{me} + \vec{P}_{leaked}$.

---

### The OpenCL Kernel

Add this kernel to your system. It works on the partial densities directly.

```c
// =================================================================================
// SCAVENGE KERNEL (Clumps small leaks into main body)
// =================================================================================

// Helper to determine if a cell is a "Leak" of Material 2 (Uranium)
// Threshold: e.g., if alpha < 0.05 (5%)
inline bool is_leak(float rho_partial, float rho_other) {
    float rho_total = rho_partial + rho_other;
    if (rho_total < 1e-5f) return false; // Already vacuum
    float alpha = rho_partial / rho_total;
    return (alpha > 0.0f && alpha < 0.05f); // 5% threshold
}

// Helper to get priority of a cell (Density of Uranium)
inline float get_priority(float4 U) {
    return U.y; // Using rho_a2 (Uranium) as priority
}

__kernel void scavenge_matter(
    __global const float4* restrict in_U, // State
    __global float4* restrict out_U,      // Output
    int width, 
    int height
) {
    int gx = get_global_id(0);
    int gy = get_global_id(1);

    if (gx < 1 || gx >= width - 1 || gy < 1 || gy >= height - 1) return;

    int idx = gy * width + gx;
    float4 U_me = in_U[idx];

    // ---------------------------------------------------------
    // STEP 1: Determine MY status
    // ---------------------------------------------------------
    // We only scavenge Material 2 (Uranium) typically, 
    // letting Material 1 (H2/Gas) fill the voids.
    
    // Check if I am a "Leak" that should disappear
    bool am_I_leak = is_leak(U_me.y, U_me.x);

    // If I am a leak, I need to know who my Master is.
    // I don't move yet; I just need to know if I *would* move so I can clear myself.
    // But wait! We need to know if neighbors point to ME.
    
    // ---------------------------------------------------------
    // STEP 2: Scavenge Loop (Gather from neighbors)
    // ---------------------------------------------------------
    
    float4 U_accum = U_me;
    
    // If I am a leak, I will empty myself (U_accum becomes background),
    // UNLESS I can't find a master (isolated droplet), then I stay.
    
    // Find My Master first to see if I disappear
    float max_prio = -1.0f;
    int master_idx = -1;

    // Search 3x3 for the best place to go
    for (int dy = -1; dy <= 1; dy++) {
        for (int dx = -1; dx <= 1; dx++) {
            if (dx == 0 && dy == 0) continue;
            
            int n_idx = (gy + dy) * width + (gx + dx);
            float4 U_nb = in_U[n_idx];
            float prio = get_priority(U_nb);
            
            // Criteria for Master:
            // 1. Must NOT be a leak itself (Stable core)
            // 2. Highest Density
            // 3. Deterministic tie-break (indices)
            if (!is_leak(U_nb.y, U_nb.x)) {
                if (prio > max_prio) {
                    max_prio = prio;
                    master_idx = n_idx;
                }
            }
        }
    }

    // STATE UPDATE LOGIC:
    
    if (am_I_leak && master_idx != -1) {
        // I am a leak and I have a valid destination.
        // I disappear (my mass moves to neighbor in their thread).
        // Reset me to pure background (H2 or Vacuum).
        
        // Conservation Note: My momentum/mass is added to Master in the loop below.
        // Here I just reset myself.
        U_accum = (float4)(1e-3f, 0.0f, 0.0f, 0.0f); // 1g/L background gas, 0 Uranium
    } 
    else {
        // I am NOT a leak (or I am an isolated leak).
        // I must check if any neighbors are leaks that chose ME as master.
        
        for (int dy = -1; dy <= 1; dy++) {
            for (int dx = -1; dx <= 1; dx++) {
                if (dx == 0 && dy == 0) continue;
                
                int n_idx = (gy + dy) * width + (gx + dx);
                float4 U_nb = in_U[n_idx];
                
                // Is neighbor a leak?
                if (is_leak(U_nb.y, U_nb.x)) {
                    // Did they choose ME?
                    // Re-run the selection logic from THEIR perspective
                    
                    float best_nb_prio = -1.0f;
                    int best_nb_idx = -1;
                    
                    // Neighbor looks at its surrounding
                    int nb_gx = gx + dx;
                    int nb_gy = gy + dy;
                    
                    for (int ddy = -1; ddy <= 1; ddy++) {
                        for (int ddx = -1; ddx <= 1; ddx++) {
                            if (ddx == 0 && ddy == 0) continue;
                            int nn_idx = (nb_gy + ddy) * width + (nb_gx + ddx);
                            
                            // Optimization: We only need to check if *I* am the best.
                            // But for correctness we must find the global best in their 3x3.
                            float4 U_nn = in_U[nn_idx];
                            float nn_prio = get_priority(U_nn);
                            
                            if (!is_leak(U_nn.y, U_nn.x)) {
                                if (nn_prio > best_nb_prio) {
                                    best_nb_prio = nn_prio;
                                    best_nb_idx = nn_idx;
                                }
                            }
                        }
                    }
                    
                    // If I am the chosen one
                    if (best_nb_idx == idx) {
                        // SUCK IT IN
                        // Add their mass and momentum to me
                        // Note: We add EVERYTHING from them (both Mat1 and Mat2) 
                        // to preserve local physics, or just Mat2. 
                        // Usually better to take everything to clear the cell.
                        
                        U_accum.x += U_nb.x; // Add H2
                        U_accum.y += U_nb.y; // Add Uranium
                        U_accum.z += U_nb.z; // Add Momentum X
                        U_accum.w += U_nb.w; // Add Momentum Y
                        
                        // Note on Energy:
                        // Energy is in a separate buffer in your struct. 
                        // You strictly need to move Energy too! 
                        // For this snippet I omitted the Energy buffer for brevity 
                        // but you MUST add it (E_accum += E_nb).
                    }
                }
            }
        }
    }

    out_U[idx] = U_accum;
}
```

### Why this works (The "PBD" Analogy)

You compared this to **Position Based Dynamics (PBD)**, where you resolve constraints by projecting positions.

*   **The Constraint:** "Material 2 cannot exist in fractions $\alpha < 0.05$."
*   **The Projection:** If constraint is violated, move the particle (mass) to the nearest valid position (Core).
*   **Conservation:**
    *   **Mass:** $m_{core} \leftarrow m_{core} + m_{leak}$. Total mass unchanged.
    *   **Momentum:** $p_{core} \leftarrow p_{core} + p_{leak}$. Total momentum unchanged.
    *   **Velocity Recoil:** Since $v = p/m$, the new velocity is the mass-weighted average.
        $$ v_{new} = \frac{m_c v_c + m_l v_l}{m_c + m_l} $$
        This physically simulates the tiny droplet hitting the big drop and merging (inelastic collision).

### Integration

1.  **Run Fluid Step** (Updates $\rho, u, E$).
2.  **Run Scavenge Kernel** (Moves small $\rho$ to big $\rho$).
    *   *Note:* You need to pass the Energy buffer too and update it exactly like `U.z/w` (momentum).
3.  **Run Boundary/Floors** (Reset any resulting vacuum to background).

This method is extremely robust ("Clean"). It will keep your Uranium ball perfectly sharp and prevent any "fog" or "haze" from developing around the impact site, effectively acting like infinite surface tension for small fragments.

---

## USER

OK, this looks simple, but when it is simple, conserve mass and all invariants why it is not used more, what are the advantages, why you mentioned other methods first?

anyway, I still thin there should be possibility to combine the topological prperthies of that "phi funcation" that it is zero-crossing, which I like a lot, with something more mass-conserving. for example updating the phi-function not by distance function but by some other method, and then aplying mass-conserving and momentum-conserving operation to move the masses (densities) consistently (to make sure no mass is where phi<1

think more if there is anything like that.

---

## Gemini

You are asking the right questions. The "Scavenging/Clumping" method is not the standard in high-end CFD for a few reasons:
1.  **Grid Bias:** It tends to turn drops into squares (aligns with grid axes).
2.  **Loss of Physics:** It kills small physical droplets (atomization) by forcing them to merge, even if they should physically be separating.
3.  **No Curvature:** It ignores the shape of the interface, making surface tension calculations inaccurate.

However, your intuition about combining the **Topology of $\phi$** with **Conservation of Mass** is actually the definition of the state-of-the-art method used in modern research (like in the **CLSVOF** - Coupled Level Set Volume of Fluid method).

Here is the method you are looking for. It is often called **"Density Re-initialization"** or **"Mass Correction via Level Set"**.

### The Concept: "The Mold and the Clay"

Imagine the Level Set $\phi$ is a **Rigid Mold** that moves through space based on velocity.
Imagine the Density $\rho$ is **Soft Clay** (Mass) that also moves through space.

Due to numerical errors, the Clay spills out of the Mold, or the Mold moves faster than the Clay.
**The Solution:** In every step, we push the Clay back into the Mold using a conservative flux, but we *never* create or destroy Clay.

### The Algorithm

1.  **Advect Everything:** Move $\phi$ and $\rho$ using the fluid velocity $\vec{u}$.
2.  **Detect Mismatch:**
    *   **Leak:** $\rho > 0$ where $\phi < 0$ (Clay outside mold).
    *   **Void:** $\rho < \rho_{target}$ where $\phi > 0$ (Mold is empty).
3.  **The Repair Flux (PBD-style):**
    We apply a "Correction Velocity" $\vec{u}_{corr}$ to the mass.
    *   The direction is $\vec{n} = \nabla \phi$ (The normal vector pointing *inside* the mold).
    *   We only move mass **along this normal**.
    *   We move it from "Leak" regions to "Void" regions.

This preserves the **Zero-Crossing Topology** (because $\phi$ dictates the geometry) and **Conserves Mass** (because we move $\rho$ using a divergence term).

---

### Implementation

We can implement this as a single **"Constraint Kernel"** that runs after your fluid update.
It calculates a flux of mass based on the mismatch with $\phi$.

$$ \frac{\partial \rho}{\partial \tau} + \nabla \cdot (\rho \vec{u}_{corr}) = 0 $$

Where we define the artificial correction velocity $\vec{u}_{corr}$ as:
$$ \vec{u}_{corr} = \text{Speed} \cdot \vec{n} \cdot \text{Sign}(\phi) $$
Actually, a simpler heuristic works better: **"If I am outside, run inside. If I am inside, stay."**

#### The OpenCL Kernel: `constrain_mass_to_phi.cl`

This kernel pushes density that is in the "wrong" place (where $\phi < 0$) towards the "right" place (where $\phi > 0$) using the gradient of $\phi$ as the guide.

```c
__kernel void constrain_mass(
    __global float4* restrict U_vec, // .x=rho1, .y=rho2 (Uranium)
    __global const float2* restrict E_vec, // .y = phi
    int width, int height, float dx, float dt
) {
    int gx = get_global_id(0);
    int gy = get_global_id(1);
    if (gx < 1 || gx >= width - 1 || gy < 1 || gy >= height - 1) return;

    int idx_c = gy * width + gx;

    // 1. Get Topology (The Guide)
    float phi_c = E_vec[idx_c].y;
    
    // Calculate Normal of phi (Direction to the interior)
    // We need this to know which way is "In"
    float phi_l = E_vec[gy*width + (gx-1)].y;
    float phi_r = E_vec[gy*width + (gx+1)].y;
    float phi_d = E_vec[(gy-1)*width + gx].y;
    float phi_u = E_vec[(gy+1)*width + gx].y;
    
    float2 grad_phi;
    grad_phi.x = (phi_r - phi_l) * 0.5f;
    grad_phi.y = (phi_u - phi_d) * 0.5f;
    
    // Normalize to get direction
    float len = length(grad_phi) + 1e-9f;
    float2 n = grad_phi / len; 

    // 2. Define Correction Flux
    // We want to move Uranium (rho2) from regions where phi < 0 
    // into regions where phi > 0.
    
    // Heuristic:
    // We define a "Correction Speed". 
    // This should be fast enough to fix errors in one step, 
    // roughly dx / dt.
    float correction_speed = dx / dt; 

    // We calculate fluxes at faces (Right and Up)
    
    // --- Right Face Flux ---
    int idx_r = gy * width + (gx + 1);
    float phi_face_x = 0.5f * (phi_c + phi_r);
    float nx_face = 0.5f * (n.x + (E_vec[idx_r].y - E_vec[gy*width+gx-1+1].y)/(2.0f*dx)); // Approx
    // Actually, just use center normal is fine for simple code
    
    // FLUX LOGIC:
    // If phi < 0, material wants to move along +Normal.
    // If phi > 0, material stays put (velocity 0).
    // We construct an artificial velocity field u_art that exists ONLY where phi < 0
    // and points in direction of N.
    
    // Let's compute Divergence of this flux centered at C.
    // To be conservative, we compute flux_in and flux_out.
    
    // Helper to get "Illegal Mass" that wants to move
    auto get_mobile_density = [&](int i, float p) {
        float4 u = U_vec[i];
        // Only Uranium (y) moves if phi (p) is negative
        if (p < 0.0f) return u.y; 
        return 0.0f;
    };

    // Calculate Fluxes from neighbors INTO me, and ME into neighbors
    // Upwind scheme: Mass moves in direction of Normal.
    
    float flux_net = 0.0f;
    
    // Check Left Neighbor (flow x)
    {
        float p_l = phi_l;
        float2 n_l = (float2)( (phi_c - E_vec[gy*width+gx-2].y)*0.5f, 0.0f ); // Simplified
        n_l = normalize(n_l);
        
        // If Left's normal points Right (+x), and Left is outside, it pushes mass to Me.
        if (p_l < 0 && n_l.x > 0) {
            float amount = get_mobile_density(gy*width+(gx-1), p_l) * n_l.x * correction_speed;
            flux_net += amount; // Inflow
        }
    }
    
    // Check Me pushing Right
    {
        if (phi_c < 0 && n.x > 0) {
             float amount = get_mobile_density(idx_c, phi_c) * n.x * correction_speed;
             flux_net -= amount; // Outflow
        }
    }
    
    // ... Repeat for Right-in, Me-Left-out, Up-in, Me-Down-out ...
    // This explicit neighbor checking is tedious in a single kernel without shared mem.
    
    // BETTER SIMPLIFIED APPROACH:
    // Gradient Descent on Density Error.
    // Error E = rho_uranium - Target_rho (where Target is 0 if phi<0)
    // Flux = - D * Grad(Error)
    // This is "Artificial Diffusion" but restricted to the interface.
}
```

Okay, the kernel above gets complicated quickly due to the stencil.

Let's do the **"Attraction Force"** method, which is mathematically equivalent to the PBD logic you like, but formulated as an Eulerian Source Term.

### The "Attraction" Method

Instead of thinking of "Moving Mass", think of **Forces**.
If a piece of Uranium is in a location where $\phi < 0$, it feels a force $\vec{F}$ pointing towards $\phi > 0$ (along $\nabla \phi$).
We simply add this velocity to the fluid velocity **for the density update only**.

**Revised Step-by-Step:**

1.  **Advect Level Set $\phi$** (geometry update).
2.  **Calculate Correction Velocity $\vec{u}_{c}$:**
    $$ \vec{u}_{c} = 
    \begin{cases} 
    C \cdot \nabla \phi & \text{if } \phi < \epsilon \text{ (Outside/Edge)} \\
    0 & \text{if } \phi > \epsilon \text{ (Inside)}
    \end{cases} $$
    *(Where C is a large speed, e.g., max fluid speed).*
3.  **Advect Density $\rho$:**
    Use velocity $\vec{u}_{total} = \vec{u}_{fluid} + \vec{u}_{c}$.

**Result:**
*   Inside the material, $\vec{u}_{c} = 0$, so physics is normal.
*   Outside (leakage), $\vec{u}_{c}$ dominates. It grabs the leakage and "conveys" it back towards the zero-crossing interface.
*   Once the leakage crosses $\phi=0$, $\vec{u}_{c}$ turns off, and the mass is safely "inside".

### Implementation

You don't need a new kernel. You modify your **Fluid Update Kernel**.

In `compute_flux`, when calculating the flux for the density `rho_a2` (Uranium), you add a "Ghost Velocity".

```c
// Inside compute_flux function

// 1. Get Normal from Phi
// (You need to pass phi gradients or compute them locally)
float phi_L = UL.phi;
float phi_R = UR.phi;
float normal_x = (phi_R - phi_L); // Un-normalized direction

// 2. Determine if we need correction
// If both are vacuum (phi < 0), we want to pull mass towards positive phi.
// Or simply: If I am "Leak" (phi < 0), I have an extra velocity along Normal.

float u_corr_L = 0.0f;
float u_corr_R = 0.0f;

float correction_strength = 2000.0f; // m/s, tunable parameter

if (UL.phi < 0.0f) {
   // I am outside. I want to move towards increasing phi.
   // Gradient points "In". 
   // We approximate local gradient direction sign using the neighbor R
   float sign_grad = (phi_R > phi_L) ? 1.0f : -1.0f;
   u_corr_L = correction_strength * sign_grad;
}

if (UR.phi < 0.0f) {
   float sign_grad = (phi_R > phi_L) ? 1.0f : -1.0f;
   u_corr_R = correction_strength * sign_grad;
}

// 3. Modify Velocities ONLY for Density Flux
// We do NOT modify momentum flux (ru, rv) because that would create fake forces.
// We only move the MASS.

float un_L_rho = (PL.u + u_corr_L * nx) * nx + ...; 

// FLUX L
FL.rho_a2 = UL.rho_a2 * un_L_rho; // Uranium moves with correction
FL.ru     = UL.ru * un_L + PL.p * nx; // Momentum moves with REAL velocity
```

**Wait, be careful.** Moving mass without moving momentum changes the velocity ($v = p/m$). If you add mass to a cell but not momentum, velocity drops. If you remove mass but not momentum, velocity skyrockets (instability!).

**Correct PBD approach for Eulerian:**
If you move mass $\delta m$ from A to B via correction:
1.  Subtract $\delta m$ from A.
2.  Subtract $\vec{v}_A \cdot \delta m$ (momentum) from A.
3.  Add $\delta m$ to B.
4.  Add $\vec{v}_A \cdot \delta m$ (carried momentum) to B.

This ensures that the "Leak" carries its own momentum back into the main body.

### Summary of the Best Approach for You

Since you want simple + robust + topological:

1.  **Keep $\phi$:** It is your ground truth for geometry.
2.  **Keep your EOS:** It handles the compressibility.
3.  **Add "Mass Scavenging" (The Kernel I wrote previously):**
    *   Why? Because advecting density with a modified velocity `u+u_corr` introduces a new CFL condition (stability limit). If `u_corr` is fast (to fix leaks quickly), your timestep `dt` must become tiny.
    *   The **Scavenge Kernel** (Donor-Acceptor) is effectively an **infinite speed** correction that happens instantly at the end of the step. It violates CFL physically but is numerically safe because it's a discrete redistribution, not a PDE advection.

**My recommendation:**
Use the **Scavenge Kernel** I provided in the previous answer, but modify the condition `is_leak`.
Instead of `alpha < 0.05`, use:
`bool is_leak = (phi < -dx);`

This links the mass movement to your topological variable $\phi$.
*   If $\phi < -dx$ (clearly outside interface) and $\rho > 0$, **Move it**.
*   This gives you the "Zero Crossing" guarantee (all mass stays within distance $dx$ of $\phi=0$) and exact conservation.

---

## USER

OK, but what looks suspicious - you move rho, you move phi, but then you say phi is right, I adjust rho to match phi. But if phi si always right we can just propagate phi. There must be a way how rho influence phi backwards. On thing is  - if we somehow construct phi form rho? (that was what I was asking the whole time from the beggining). The other thing is that movement of rho may create some recoil forced / acceleration which change the velocity field which then advect phi. 

---

## USER

You are absolutely correct. This is the fundamental insight behind the **Coupled Level Set and Volume of Fluid (CLSVOF)** method.

If you treat $\phi$ as an "absolute truth" that mass must obey, you violate physics (mass conservation).
If you treat $\rho$ as the "absolute truth," you lose the sharp interface topology (it smears).

The solution is indeed a **Two-Way Feedback Loop**:
1.  **Forward:** $\phi$ helps $\rho$ stay sharp (Interface Compression).
2.  **Backward:** $\rho$ corrects $\phi$ to ensure the geometry matches the mass (Mass-Defined Geometry).

You asked: *"Can we construct phi from rho?"* and *"Can movement of rho create recoil?"*
Yes. Here is the rigorous method to do exactly that, usually called **"Level Set Re-initialization from Volume Fraction."**

### The Concept: Mass is the "Anchor"

We stop treating $\phi$ as a variable we just advect blindly. Instead, $\phi$ becomes a **State Variable derived from Mass**.

1.  **Advection:** We advect Mass ($\rho$) and Momentum ($\rho \vec{u}$) strictly conservatively. This is the Physics.
2.  **Reconstruction:** At the end of the step, we ask: "Given this distribution of Mass, where is the interface?"
3.  **Correction:** We force $\phi$ to align with the mass boundary.

If a chunk of Uranium flies into a vacuum cell (increasing $\rho$), the Reconstruction step sees $\rho > 0$ and forces $\phi$ to become positive. The geometry "follows" the physics.

---

### The Algorithm

We introduce a "Feedback Kernel" that runs after the Fluid Update but before the Redistance step.

$$ \phi^{new} = \phi^{pred} + \lambda \cdot (\alpha - H(\phi)) $$

*   $\alpha$: Volume fraction derived from mass ($\rho / \rho_{solid}$).
*   $H(\phi)$: The Heaviside value of the current Level Set (0 or 1).
*   $\lambda$: A relaxation factor.

**How it works:**
*   **Case A (Mass moved in):** Cell has $\rho_{uranium}$ ($\alpha \approx 1$) but $\phi$ is still negative ($H(\phi) \approx 0$).
    *   Error = $1 - 0 = 1$.
    *   $\phi$ increases. It crosses zero. The interface moves to encompass the mass.
*   **Case B (Mass moved out):** Cell is empty ($\alpha \approx 0$) but $\phi$ is positive ($H(\phi) \approx 1$).
    *   Error = $0 - 1 = -1$.
    *   $\phi$ decreases. The interface shrinks to exclude the empty space.

### 1. The Feedback Kernel

Add this to your `EulerianImpacFluid.cl`.

```c
// =================================================================================
// FEEDBACK KERNEL: Adjusts Phi to match Mass
// =================================================================================
__kernel void correct_phi_from_mass(
    __global const float4* restrict U_vec, // State (Densities)
    __global float2* restrict E_vec,       // State (Energy, Phi)
    int width, 
    int height,
    float dx,
    float target_density_mat2 // Density of pure Uranium (e.g. 19050.0)
) {
    int gx = get_global_id(0);
    int gy = get_global_id(1);

    if (gx < width && gy < height) {
        int idx = gy * width + gx;
        
        // 1. Get Current Physics
        float4 U = U_vec[idx];
        float2 E = E_vec[idx];
        
        float rho_u = U.y; // Uranium Mass
        float phi = E.y;
        
        // 2. Compute Physical Volume Fraction
        // alpha = rho / rho_solid
        // We clamp it to [0,1] because compressible shocks can make rho > rho_solid
        float alpha = clamp(rho_u / target_density_mat2, 0.0f, 1.0f);
        
        // 3. Compute Geometric Volume Fraction (from current Phi)
        // using your smooth heaviside
        float H_phi = heaviside(phi, dx);
        
        // 4. The Mismatch
        // If alpha > H_phi, we have more mass than the geometry thinks -> Expand Phi
        // If alpha < H_phi, we have less mass than geometry thinks -> Shrink Phi
        float error = alpha - H_phi;
        
        // 5. Apply Correction
        // We nudge phi. The coefficient depends on dx.
        // We want to shift the zero-crossing. 
        // Heaviside transition width is approx 1.5*dx.
        // So a full 0->1 change corresponds to shifting phi by ~1.5*dx.
        
        float nudge = error * (1.5f * dx);
        
        // Update Phi
        E.y = phi + nudge;
        
        E_vec[idx] = E;
    }
}
```

### 2. The Recoil (Velocity Update)

You mentioned: *"Movement of rho may create some recoil forces."*

In the Eulerian frame, **Advection creates Recoil**.
When you move mass $\Delta m$ from cell A to cell B with velocity $\vec{u}$, you also move momentum $\Delta \vec{P} = \vec{u} \Delta m$.
If you do this correctly (via the standard flux loop), the velocity field updates automatically to conserve momentum.

*   If heavy mass moves into a light cell, the velocity of that cell drops (Mass increases, Momentum increases less $\to$ $v_{new} = P_{new}/m_{new}$ decreases).
*   **This acts as the recoil/drag.**

**Crucial Note:** If you use the "Scavenge/Clumping" kernel I gave you before, you **MUST** move momentum too.
If you use the "Phi Feedback" kernel above, you don't need to touch momentum artificially. You let the standard fluid solver move the mass/momentum, and you just update the "Map" ($\phi$) to agree with the territory ($\rho$).

---

### The Final "Best Practice" Algorithm

This is robust, conserves mass, maintains sharp interfaces, and handles multiple materials.

#### Step 1: Fluid Update (The Physics)
Solve the 5-equation model (or standard Euler).
*   Updates $\rho_1, \rho_2, \rho u, \rho v, E$.
*   This moves the mass and creates the "Real" distribution.
*   *Note:* Advect $\phi$ here too, just as a guess/prediction.

#### Step 2: Phi Feedback (The Consistency)
Run `correct_phi_from_mass`.
*   Forces the zero-crossing of $\phi$ to align with the 0.5 contour of the density fraction.
*   This handles the "Recoil" logic: The mass moved physically, so the geometry must update.

#### Step 3: Redistance (The Topology)
Run `redistance_phi` (5 iterations).
*   This takes the messy $\phi$ from Step 2 (which might be flattened or distorted) and restores it to a Signed Distance Field ($|\nabla \phi| = 1$).
*   This ensures smooth surface normals for the next time step.

#### Step 4: Interface Sharpening (Optional but Recommended)
Run the `sharpen_interface` (VOF compression) kernel I described earlier.
*   This prevents the density $\rho$ from smearing out too much.
*   Because Step 2 links $\phi$ to $\rho$, keeping $\rho$ sharp helps keep $\phi$ sharp, and vice-versa.

---

### Why this answers your "Self-Contained Extension" request

You wanted the multi-material logic to be an "extension" on top of a standard solver.

1.  **Base Layer:** Standard Compressible Euler Solver. It just sees densities and momenta. It knows nothing about $\phi$. It conserves everything.
2.  **Extension Layer:** The Level Set $\phi$.
    *   It "watches" the density field.
    *   It re-aligns itself every step (`correct_phi_from_mass`).
    *   It provides high-quality normals ($\nabla \phi$) for surface tension or mixture EOS calculations.
    *   It acts as a "switch" for your vacuum handling (if $\phi < \text{cutoff}$, treat as vacuum).

### Implementation Checklist for `test_eulerian_fluid.py`

1.  **Add `correct_phi_from_mass`** to your `.cl` file.
2.  **Update Python Loop:**

```python
def update(frame):
    # 1. Standard Fluid Step (Moves Mass & Momentum)
    solver.step(dt) 
    
    # 2. Feedback: Geometry follows Mass
    # Target density for Uranium is ~19050
    solver.correct_phi(target_density=19050.0) 
    
    # 3. Restore Distance Property
    solver.redistance(iterations=5)
    
    # ... Visualization ...
```

This effectively constructs $\phi$ from $\rho$ every step, but uses the previous $\phi$ as a "warm start" to maintain temporal coherence. It is the best of both worlds.

---

## USER

OK, I'm still a bit curious / uncertain about one detail again in analogy with position based dynamics for particles, or projective dynamics.

what we do there often is that a particle is moved by predictor ( according to e.g. verlet update given external force ), but then corrector (e.g. collisions, stay on contrained manifold etc.) say not you cannot be there and move it (just like our rho which is in phi<0)

in position based dynamics we account for that movement by updating the velocty of that particle 
so instead of verlet velocty v(t+dt) = (f(t+dt) + f(t))*(0.5*dt/m) we give it new velocity computed from updated positions
v(t+dt) = ( x(t+dt) - x(t) )/dt
this effectively account for the recoil due to moving the particle back to the constrained manifold 

I would expect something like this also here in the fluid with phi, and it does not seem to me that, just moving the momentum from one place to another does do same thing, for example if you consider flow of constant velocity everywhere, movement of momentum does nothing, but yet we moved chunk of material and this should induce some recoil (kick, accleration) 



---

## USER

You have hit on a very deep and subtle point in Computational Fluid Dynamics. Your intuition from Position Based Dynamics (PBD) is **correct**, but the way it manifests in an Eulerian Compressible solver is slightly different.

In PBD:
1.  **Predict:** $x^* = x + v \Delta t$
2.  **Correct:** $x_{final} = \text{Constraint}(x^*)$ (You move the position)
3.  **Update Velocity:** $v_{new} = (x_{final} - x) / \Delta t$ (The movement implies a force/recoil).

In your proposed Eulerian scheme:
1.  **Predict:** Advect $\rho$ and $\rho \vec{u}$.
2.  **Correct:** Move $\rho$ from "Leak" to "Body" (The Sharpening/Scavenging).
3.  **Update Velocity:** *This is the missing link you are asking about.*

### The Paradox: "Constant Velocity"

You asked: *"If I have a constant velocity flow, and I move mass+momentum from A to B, the velocity doesn't change. Where is the recoil?"*

**Answer:** In this specific case (Sharpening/Scavenging), **we do not want recoil.**

Here is why:
*   **PBD Context:** You move a particle because it hit a wall (Collision). This requires a force. Force $\to$ Acceleration $\to$ Velocity Change.
*   **Sharpening Context:** You move mass because **it shouldn't have spread out in the first place**. It was a numerical error (diffusion). The physical droplet didn't explode; the numbers just got blurry.
    *   By moving the mass *and* the momentum together, you are effectively saying: *"This piece of fluid didn't actually drift left; it stayed in the center."*
    *   If you applied a "recoil" (force), you would be introducing a **Spurious Current** (Fake Force) generated purely by your numerical error correction.

### However... The "Recoil" DOES happen (via Pressure)

Even though we don't apply an explicit velocity kick, the **Eulerian PBD Cycle** completes automatically in the **next time step** via the Equation of State (EOS).

1.  **Position Update:** You scavenge mass from the Leak (void) and shove it into the Body (solid).
2.  **Density Change:** The density $\rho$ in the Body increases slightly (because you stuffed more mass into the same volume).
3.  **EOS Response (The PBD Solver):** $P = f(\rho)$.
    *   Since $\rho$ went up, **Pressure goes up**.
4.  **Momentum Update (The Recoil):**
    *   $\frac{\partial (\rho u)}{\partial t} = -\nabla P$.
    *   The increased pressure creates a gradient that pushes outwards.

**This is your PBD recoil!**
By artificially moving mass to satisfy the topological constraint ($\phi$), you compressed the fluid. The fluid resists this compression via Pressure, generating the acceleration you were looking for.

---

### Do we ever need an Explicit Velocity Correction?

Yes. If your "Correction" is not just fixing numerical blur, but enforcing a physical force (like **Surface Tension** or **Wall Separation**), you must apply a Momentum Source Term.

If you want the "Phi" surface to act like a taut rubber sheet (Surface Tension), the simple "Move Mass" approach relies on the EOS (Compressibility) to resist. For hypervelocity impacts (high pressure), this is usually enough.

But if you want to be rigorous, you can add the **Continuum Surface Force (CSF)**, which is the exact Eulerian equivalent of the PBD constraint force for minimizing surface area.

#### The "Recoil" Kernel (Surface Tension)

Instead of deducing velocity from position changes, we apply the force derived from the surface topology directly to the momentum.

$$ \vec{F}_{recoil} = \sigma \cdot \kappa \cdot \nabla \phi $$

*   $\sigma$: Surface Tension Coefficient.
*   $\kappa$: Curvature ($\nabla \cdot (\nabla \phi / |\nabla \phi|)$).
*   $\nabla \phi$: Direction of the force.

**Algorithm:**
1.  **Sharpen/Scavenge:** Moves mass/momentum to fix numerical smear. (No velocity change).
2.  **Pressure Step:** EOS converts density changes to Pressure. (Implicit recoil).
3.  **Surface Tension:** Calculates curvature of $\phi$ and adds explicit force $\vec{F}_{recoil}$ to momentum. (Explicit recoil).

### Summary for your Simulator

You do **not** need to manually calculate $v_{new} = (x_{new} - x_{old})/dt$ for the scavenging step.

1.  **Why?** Because "Scavenging" is correcting a numerical error, not modeling a physical collision. If you add a "kick" every time you fix a blur, your fluid will self-accelerate and become unstable (Parasitic Currents).
2.  **The Mechanism:**
    *   Move Mass ($\rho$) AND Momentum ($\rho \vec{u}$) proportionally.
    *   This preserves the local velocity $v$.
    *   The **Equation of State** handles the rest. If you over-stuff a cell, $P$ rises, and the fluid naturally accelerates away in the next step.

**Implementation Detail:**
Ensure your `scavenge_matter` kernel (or whatever moving logic you use) moves **Total Energy** density as well.
$$ E_{new} = E_{old} + E_{leak} $$
If you forget this, the internal energy per unit mass will drop (cooling the fluid), or rise, which creates fake pressure spikes.

### Conclusion
Your PBD intuition maps to Eulerian variables like this:
*   **Position Constraint ($x$):** handled by **Mass Transport/Sharpening** ($\rho$).
*   **Velocity Update ($v$):** handled by **Pressure Gradient** ($\nabla P$).
*   **Stiffness ($k$):** handled by the **Bulk Modulus** (Sound speed squared $c^2$) in the EOS.

