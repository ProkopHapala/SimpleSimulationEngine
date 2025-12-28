import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.collections import PolyCollection
from matplotlib import cm
from numba import njit

parser = argparse.ArgumentParser(description="Spherical Implosion v4 (Numba) with adjustable parameters")
parser.add_argument("--zones",        type=int,   default=30,     help="Number of radial zones")
parser.add_argument("--t-end",        type=float, default=20.0e-6, help="End time (s)")
parser.add_argument("--cfl",          type=float, default=0.3,    help="CFL safety factor")
parser.add_argument("--max-e-change", type=float, default=0.1,    help="Max fractional energy change per step")
parser.add_argument("--visc-lin",     type=float, default=0.5,    help="Linear artificial viscosity coeff")
parser.add_argument("--visc-quad",    type=float, default=2.0  ,  help="Quadratic artificial viscosity coeff")
parser.add_argument("--driver-vel",   type=float, default=5000.0, help="Initial inward velocity for tamper nodes (m/s)")
parser.add_argument("--substeps",     type=int,   default=20,     help="Physics substeps per animation frame")
parser.add_argument("--dt-initial",   type=float, default=1e-10,  help="Initial timestep (s)")
args = parser.parse_args()

# Extract CLI to module-level scalars (avoid passing argparse into hot loops)
DRIVER_VEL = args.driver_vel
SUBSTEPS = args.substeps
DT_INITIAL = args.dt_initial

# ==========================================
# PART 1: CONFIGURATION & CONSTANTS
# ==========================================

# Physics
KB = 1.38e-23
NA = 6.022e23
MEV_TO_J = 1.602e-13
A_RAD = 7.56e-16
C_LIGHT = 3.0e8

# Simulation Setup
# REDUCED RESOLUTION for clear layer visualization
N_ZONES = args.zones
T_END = args.t_end

# Stability Controls
CFL_SAFETY = args.cfl        # Lower is more stable
MAX_E_CHANGE = args.max_e_change      # Max 10% energy change per step (prevents fusion crashes)
VISC_LIN = args.visc_lin          # Von Neumann Viscosity linear
VISC_QUAD = args.visc_quad         # Von Neumann Viscosity quadratic

# Material IDs
MAT_GAS = 0
MAT_PUSHER = 1
MAT_TAMPER = 2

# ==========================================
# PART 2: ACCELERATED PHYSICS KERNELS
# ==========================================

@njit(fastmath=True)
def get_eos(rho, e_int, gamma, cv):
    # Temperature from internal energy
    T = e_int / cv
    if T < 300.0: T = 300.0
    
    # Pressure = Ideal/Stiff Gas + Radiation
    p_mat = (gamma - 1.0) * rho * e_int
    p_rad = (1.0/3.0) * A_RAD * (T**4)
    
    return p_mat + p_rad, T

@njit(fastmath=True)
def calc_fusion_rate(T):
    # Bosch-Hale approximation for DT fusion
    keV = 1.16e7
    T_keV = T / keV
    if T_keV < 1.0: return 0.0
    
    # Reactivity <sigma v>
    sig_v = 3.68e-18 * (T_keV**(-2/3)) * np.exp(-19.94 * (T_keV**(-1/3)))
    return sig_v

@njit
def physics_step_staggered(
    r, v, rho, mass, e_int, 
    mat_gamma, mat_cv, mat_type, mat_A,
    neutron_flux, dt_prev
):
    """
    STAGGERED GRID SCHEME (Von Neumann-Richtmyer)
    v is defined at half-steps (n+1/2)
    r, rho, e are defined at full-steps (n)
    """
    N = len(rho)
    
    # 1. Calculate Geometry & Density (Step n)
    # ----------------------------------------
    vol = np.zeros(N)
    r_clip = np.copy(r)
    # Prevent crossing (negative volume)
    for k in range(1, N+1):
        if r_clip[k] <= r_clip[k-1] + 1e-7:
            r_clip[k] = r_clip[k-1] + 1e-7

    for i in range(N):
        vol[i] = (4.0/3.0) * np.pi * (r_clip[i+1]**3 - r_clip[i]**3)
        rho[i] = mass[i] / vol[i]

    # 2. Calculate Pressure & Viscosity (Step n)
    # ------------------------------------------
    P_total = np.zeros(N)
    q = np.zeros(N)
    Te = np.zeros(N)
    
    for i in range(N):
        p_eos, t_val = get_eos(rho[i], e_int[i], mat_gamma[i], mat_cv[i])
        Te[i] = t_val
        
        # Artificial Viscosity (Shock Handling)
        # Using previous velocity to estimate compression
        # div_v approx
        u_diff = v[i+1] - v[i]
        if u_diff < 0:
            c_s = np.sqrt(mat_gamma[i] * p_eos / rho[i])
            q[i] = rho[i] * (VISC_LIN * c_s * abs(u_diff) + VISC_QUAD * u_diff**2)
        
        P_total[i] = p_eos + q[i]

    # 3. Time Step Calculation (Stability Check)
    # ------------------------------------------
    # Hydro CFL
    dr_min = np.min(r_clip[1:] - r_clip[:-1])
    cs_max = np.max(np.sqrt(mat_gamma * P_total / rho))
    dt_hydro = CFL_SAFETY * dr_min / (cs_max + 1e-5)
    
    # Clamp dt to prevent massive jumps
    dt = min(dt_hydro, 1.1 * dt_prev) 
    dt = max(dt, 1e-13)
    if dt > 1e-9: dt = 1e-9

    # 4. Accelerate (Update v to n+1/2)
    # ---------------------------------
    # a = Force / NodeMass
    # Force = -Area * Grad(P)
    new_v = np.copy(v)
    for i in range(1, N): # Inner nodes
        area = 4.0 * np.pi * r_clip[i]**2
        m_node = 0.5 * (mass[i-1] + mass[i])
        
        force = -area * (P_total[i] - P_total[i-1]) # P_total includes q
        accel = force / m_node
        new_v[i] = v[i] + accel * dt

    new_v[0] = 0.0 # Center symmetry
    # Outer boundary is free (no pressure outside)
    
    # 5. Move (Update r to n+1)
    # -------------------------
    new_r = r_clip + new_v * dt
    
    # 6. Energy Update (Work + Sources)
    # ---------------------------------
    new_e_int = np.copy(e_int)
    new_neutron_flux = np.copy(neutron_flux)
    fusion_power = np.zeros(N)  # W/m3
    
    # A. PdV Work (Mechanical Heating)
    # Using v(n+1/2) * dt * P(n) is the symplectic choice
    for i in range(N):
        # dV/dt approx
        dv_vol = 4.0 * np.pi * r_clip[i]**2 * new_v[i] # Inner surface velocity
        dv_vol_next = 4.0 * np.pi * r_clip[i+1]**2 * new_v[i+1] # Outer surface
        dV = (dv_vol_next - dv_vol) * dt
        
        work = -P_total[i] * dV
        new_e_int[i] += work / mass[i]

    # B. Fusion & Neutronics (Explicit Source)
    # Limit max energy change to avoid instability
    for i in range(N):
        if mat_type[i] == MAT_GAS:
            n_ion = (rho[i] * 1000.0 / mat_A[i]) * NA
            rate = calc_fusion_rate(Te[i])
            
            # Reactions/m3
            # Limit reaction rate if it causes >10% energy jump
            raw_yield = (0.25 * n_ion**2) * rate * dt
            
            alpha_E = raw_yield * (3.5 * MEV_TO_J) / rho[i]
            
            # STABILITY CLAMP: Don't let energy explode in one step
            if alpha_E > MAX_E_CHANGE * new_e_int[i]:
                alpha_E = MAX_E_CHANGE * new_e_int[i]
                # Scale back effective yield
                raw_yield = alpha_E * rho[i] / (3.5 * MEV_TO_J)
            
            new_e_int[i] += alpha_E
            
            # Neutrons (14.1 MeV) go to flux
            # Simplified coupling: Source term
            new_neutron_flux[i] += raw_yield * C_LIGHT # Roughly flux contribution
            # Fusion power density (total 17.6 MeV per reaction)
            fusion_power[i] = raw_yield * (17.6 * MEV_TO_J) / dt
            
    # C. Neutron Diffusion (Simplified Explicit)
    # Very rough approx to spread the heat
    flux_diff = np.zeros(N)
    for i in range(1, N-1):
        flux_diff[i] = (new_neutron_flux[i+1] - 2*new_neutron_flux[i] + new_neutron_flux[i-1])
    
    new_neutron_flux += flux_diff * 0.1 # Artificial diffusion constant
    new_neutron_flux *= 0.95 # Leakage/Absorption decay
    
    return new_r, new_v, rho, new_e_int, Te, P_total, new_neutron_flux, fusion_power, dt

# ==========================================
# PART 3: SIMULATION & VISUALIZATION
# ==========================================

class Simulation:
    def __init__(self):
        self.reset()
        
    def reset(self):
        # Grid Init
        self.r = np.linspace(0, 0.08, N_ZONES + 1)
        self.v = np.zeros(N_ZONES + 1)
        
        # Physics State
        self.rho = np.zeros(N_ZONES)
        self.e_int = np.zeros(N_ZONES)
        self.mass = np.zeros(N_ZONES)
        self.neutron_flux = np.zeros(N_ZONES)
        self.Te = np.zeros(N_ZONES)
        self.P_total = np.zeros(N_ZONES)
        self.fusion_power = np.zeros(N_ZONES)
        self.step_count = 0
        
        # Material Properties
        self.mat_gamma = np.zeros(N_ZONES)
        self.mat_cv = np.zeros(N_ZONES)
        self.mat_A = np.zeros(N_ZONES)
        self.mat_type = np.zeros(N_ZONES, dtype=np.int32)
        
        # ---------------------------
        # DEFINE LAYERS (The "Segments")
        # ---------------------------
        # Using 30 zones total:
        # 0-5: Gas (6 zones)
        # 6-15: Pusher (10 zones)
        # 16-29: Tamper (14 zones)
        
        BARN = 1e-28
        
        for i in range(N_ZONES):
            vol = (4.0/3.0) * np.pi * (self.r[i+1]**3 - self.r[i]**3)
            
            if i < 6: # GAS
                self.mat_type[i] = MAT_GAS
                rho0 = 225.0
                self.mat_gamma[i] = 1.66
                self.mat_cv[i] = 3e4
                self.mat_A[i] = 2.5
            elif i < 16: # PUSHER (Pu)
                self.mat_type[i] = MAT_PUSHER
                rho0 = 19800.0
                self.mat_gamma[i] = 3.0
                self.mat_cv[i] = 130.0
                self.mat_A[i] = 239.0
            else: # TAMPER (U)
                self.mat_type[i] = MAT_TAMPER
                rho0 = 18900.0
                self.mat_gamma[i] = 3.0
                self.mat_cv[i] = 116.0
                self.mat_A[i] = 238.0
            
            self.rho[i] = rho0
            self.mass[i] = rho0 * vol
            self.e_int[i] = self.mat_cv[i] * 300.0
            
        # DRIVER: Kick the tamper
        # Assign inward velocity to Tamper nodes
        for k in range(16, N_ZONES+1):
            self.v[k] = -abs(DRIVER_VEL) 

        self.time = 0.0
        self.dt = DT_INITIAL

    def step(self):
        # Substep for animation smoothness
        for _ in range(SUBSTEPS):
            if self.time > T_END: break
            
            (self.r, self.v, self.rho, self.e_int, 
             self.Te, self.P_total, self.neutron_flux, self.fusion_power, self.dt) = physics_step_staggered(
                self.r, self.v, self.rho, self.mass, self.e_int,
                self.mat_gamma, self.mat_cv, self.mat_type, self.mat_A,
                self.neutron_flux, self.dt
            )
            self.time += self.dt
            self.step_count += 1
            if self.step_count % 100 == 0:
                print(f"[diag] t={self.time*1e6:7.3f} us | rho0={self.rho[0]:.3e} | T0={self.Te[0]:.3e} K | flux_max={np.max(self.neutron_flux):.3e} | fusionP_max={np.max(self.fusion_power):.3e} W/m3")

# ==========================================
# PART 4: VISUALIZATION (THE WEDGE)
# ==========================================

sim = Simulation()
fig = plt.figure(figsize=(12, 10))
gs = fig.add_gridspec(3, 1, height_ratios=[3, 1, 1])

# PLOT 1: THE WEDGE (Physical Layers)
ax_sim = fig.add_subplot(gs[0])
ax_sim.set_xlim(0, 0.08)
ax_sim.set_ylim(-0.04, 0.04)
ax_sim.set_aspect('equal')
ax_sim.set_facecolor('black')
ax_sim.set_title("Nuclear Implosion (Material Layers)")
ax_sim.set_xlabel("Radius (m)")

# Create Polygons for each zone
# We will color them by Material Type + Temperature Tint
polys = []
for i in range(N_ZONES):
    polys.append(np.zeros((4, 2))) # Placeholder
    
coll = PolyCollection(polys, edgecolors='white', linewidths=0.5)
ax_sim.add_collection(coll)

# Material Legend
ax_sim.text(0.01, 0.035, "CYAN: Gas", color='cyan', fontsize=10, fontweight='bold')
ax_sim.text(0.03, 0.035, "RED: Pusher", color='red', fontsize=10, fontweight='bold')
ax_sim.text(0.05, 0.035, "GRAY: Tamper", color='gray', fontsize=10, fontweight='bold')

# PLOT 2: Pressure Piston
ax_p = fig.add_subplot(gs[1])
ax_p.set_xlim(0, 0.08)
ax_p.set_yscale('log')
ax_p.set_ylim(1e9, 1e16)
ax_p.set_ylabel("Pressure (Pa)")
ax_p.grid(True, which='both', alpha=0.3)
line_p, = ax_p.plot([], [], 'y-', lw=2)

# PLOT 3: Temperature/Neutrons/Fusion Power
ax_t = fig.add_subplot(gs[2])
ax_t.set_xlim(0, 0.08)
ax_t.set_yscale('log')
ax_t.set_ylim(300, 1e9)
ax_t.set_ylabel("Temperature (K)")
ax_t.set_xlabel("Radius (m)")
ax_t.grid(True, which='both', alpha=0.3)
line_t, = ax_t.plot([], [], 'r-', lw=2, label='Temp')
# Twin axis for flux
ax_n = ax_t.twinx()
ax_n.set_ylim(1e10, 1e35)
ax_n.set_yscale('log')
ax_n.set_ylabel("Neutrons / Fusion Power", color='cyan')
line_n, = ax_n.plot([], [], 'c--', lw=2, label='Neutrons')
line_fp, = ax_n.plot([], [], 'y:', lw=2, label='Fusion Power')

def get_layer_color(mat_type, temp):
    # Base Colors
    if mat_type == MAT_GAS:
        base = np.array([0.0, 1.0, 1.0]) # Cyan
    elif mat_type == MAT_PUSHER:
        base = np.array([1.0, 0.0, 0.0]) # Red
    else:
        base = np.array([0.5, 0.5, 0.5]) # Gray
        
    # Thermal Glow Logic
    # As T goes from 300 to 1e7, add White/Yellow
    glow = np.clip(np.log10(temp)/8.0, 0, 1) # 0 to 1 scale roughly
    
    # Mix base with White (1,1,1) based on glow
    color = base * (1 - glow*0.5) + np.array([1.0, 1.0, 0.8]) * (glow*0.8)
    return np.clip(color, 0, 1)

def update(frame):
    sim.step()
    
    # 1. Update Polygons
    verts = []
    colors = []
    
    for i in range(N_ZONES):
        r_inner = sim.r[i]
        r_outer = sim.r[i+1]
        
        # Wedge coordinates (Trapezoid)
        # y = +/- 0.5 * r to make a 45-degree wedge appearance
        xi = [r_inner, r_outer, r_outer, r_inner]
        yi = [r_inner*0.5, r_outer*0.5, -r_outer*0.5, -r_inner*0.5]
        verts.append(list(zip(xi, yi)))
        
        # Color Logic
        c = get_layer_color(sim.mat_type[i], sim.Te[i])
        colors.append(c)
        
    coll.set_verts(verts)
    coll.set_facecolors(colors)
    
    # 2. Update Lines
    centers = 0.5 * (sim.r[:-1] + sim.r[1:])
    line_p.set_data(centers, sim.P_total)
    line_t.set_data(centers, sim.Te)
    line_n.set_data(centers, sim.neutron_flux + 1e10) # offset 0
    line_fp.set_data(centers, sim.fusion_power + 1e-6) # fusion power density (offset to stay >0)
    
    ax_sim.set_title(f"Time: {sim.time*1e6:.3f} µs | Center Density: {sim.rho[0]/1000:.1f} g/cc")
    
    return coll, line_p, line_t, line_n, line_fp

ani = animation.FuncAnimation(fig, update, frames=600, interval=20, blit=False)
plt.show()