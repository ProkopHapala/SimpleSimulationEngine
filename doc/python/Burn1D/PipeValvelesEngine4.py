import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.collections import LineCollection

parser = argparse.ArgumentParser(description="1D Lagrangian pulsejet with reservoirs")
parser.add_argument("--cfl",             type=float, default=0.4,  help="CFL safety factor")
parser.add_argument("--dt-min",          type=float, default=1e-8, help="Minimum timestep clamp")
parser.add_argument("--dt-max",          type=float, default=1e-4, help="Maximum timestep clamp")
parser.add_argument("--steps-per-frame", type=int,   default=50,   help="Physics substeps per frame")
parser.add_argument("--dx-pipe",         type=float, default=0.02, help="Target dx inside pipe")
parser.add_argument("--dx-res",          type=float, default=0.05, help="Target dx inside reservoirs")
parser.add_argument("--burn-time",       type=float, default=0.01, help="Time to trigger burn (s)")
parser.add_argument("--burn-multiplier", type=float, default=4.0,  help="Temperature multiplier on burn")
parser.add_argument("--jet-damping",     type=float, default=5.0,  help="Damping when exhausting into reservoirs")
parser.add_argument("--wall-friction",   type=float, default=0.5,  help="Wall friction coefficient")
args = parser.parse_args()

# --- CONFIGURATION ---
# Physical Constants
Gamma = 1.4
P_atm = 101325.0
Rho_atm = 1.225
C_sound_atm = np.sqrt(Gamma * P_atm / Rho_atm)

# Simulation Stability
CFL_SAFETY = args.cfl
VISC_COEFF = 2.5       # Shock thickening
WALL_FRICTION = args.wall_friction    # Drag inside the pipe
JET_DAMPING = args.jet_damping      # Energy loss when exiting to atmosphere (Crucial for valveless!)

# Combustion
T_burn_start = args.burn_time
Burn_Multiplier = args.burn_multiplier

# --- 1. LOCKWOOD-HILLER GEOMETRY ---
# Defined by control points (x, radius)
# x < 0 and x > 1.6 are the "Atmosphere Reservoirs"
GEOM_POINTS = [
    (-2.0, 1.00),  # Left Reservoir start (huge)
    (-0.1, 0.08),  # Funnel into intake
    ( 0.0, 0.025),  # Intake Start
    ( 0.3, 0.025), # Intake End / Chamber Start
    ( 0.40, 0.1), # Chamber Max Width
    ( 0.50, 0.1), # Chamber Max Width
    ( 0.6, 0.035), # Chamber End / Exhaust Start
    ( 1.6, 0.05),  # Exhaust End
    ( 3.5, 1.00)   # Right Reservoir (huge)
]

def get_radius(x):
    x_arr = np.asarray(x)
    pts = np.array(GEOM_POINTS)
    xp, yp = pts[:,0], pts[:,1]
    # linear interp for robustness; endpoints clipped
    r = np.interp(x_arr, xp, yp)
    if np.isscalar(x):
        return float(r)
    return r

def get_area(x):
    r = get_radius(x)
    return np.pi * r**2

# --- INITIALIZATION ---
# We use a fixed high number of elements now, because the reservoirs handle the BCs
target_dx_pipe = args.dx_pipe
target_dx_res = args.dx_res

# Generate non-uniform grid (dense in pipe, sparse in reservoir)
xs = []
x = GEOM_POINTS[0][0]
end = GEOM_POINTS[-1][0]

while x < end:
    xs.append(x)
    # Variable resolution based on where we are
    if 0.0 <= x <= 1.6:
        step = target_dx_pipe
    else:
        step = target_dx_res
    x += step
xs = np.array(xs)

N = len(xs)
pos = xs.copy()
vel = np.zeros(N)
mass = np.zeros(N-1)
entropy = np.ones(N-1)

# Initialize Air Mass
for i in range(N-1):
    x_mid = 0.5 * (pos[i] + pos[i+1])
    area = get_area(x_mid)
    vol = area * (pos[i+1] - pos[i])
    mass[i] = vol * Rho_atm

sim_time = 0.0
burned = False

# --- PLOTTING ---
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12))

# Background Geometry
gx = np.linspace(GEOM_POINTS[0][0], GEOM_POINTS[-1][0], 1000)
gr = np.array([get_radius(x) for x in gx])
ax1.plot(gx, gr, 'k', lw=2)
ax1.plot(gx, -gr, 'k', lw=2)
# Highlight Engine Body
ax1.fill_between(gx, gr, 0.6, where=(gx>0)&(gx<1.6), color='gray', alpha=0.2)
ax1.fill_between(gx, -gr, -0.6, where=(gx>0)&(gx<1.6), color='gray', alpha=0.2)

lc = LineCollection([], linewidths=1)
ax1.add_collection(lc)
ax1.set_xlim(-0.5, 2.0) # Zoom in on engine
ax1.set_ylim(-0.3, 0.3)
ax1.set_title("Pulsejet Simulation (Zoomed)")

line_p, = ax2.plot([], [], 'r-')
ax2.set_xlim(GEOM_POINTS[0][0], GEOM_POINTS[-1][0])
ax2.set_ylim(0, P_atm * 3.5)
ax2.axhline(P_atm, color='gray', ls=':')
ax2.set_ylabel("Pressure (Pa)")

line_v, = ax3.plot([], [], 'g-')
ax3.set_xlim(GEOM_POINTS[0][0], GEOM_POINTS[-1][0])
ax3.set_ylim(-400, 400)
ax3.axhline(0, color='gray', ls=':')
ax3.set_ylabel("Velocity (m/s)")

# --- PHYSICS LOOP ---
def physics_step():
    global pos, vel, mass, entropy, sim_time, burned
    
    # 1. Geometry Update
    dx = pos[1:] - pos[:-1]
    
    # Robustness: Fix collapsed cells
    min_dx = 1e-5
    if np.any(dx < min_dx):
        bad_idx = np.where(dx < min_dx)[0]
        for i in bad_idx:
            mid = 0.5 * (pos[i] + pos[i+1])
            pos[i] = mid - 0.51 * min_dx
            pos[i+1] = mid + 0.51 * min_dx
        dx = pos[1:] - pos[:-1] # Recalc

    # Area & Volume
    A_nodes = np.array([get_area(x) for x in pos])
    A_L, A_R = A_nodes[:-1], A_nodes[1:]
    vol = (dx / 3.0) * (A_L + A_R + np.sqrt(A_L * A_R))
    
    # 2. Thermodynamics
    rho = mass / vol
    P_eos = P_atm * entropy * (rho / Rho_atm)**Gamma
    
    # 3. Artificial Viscosity (Shock Handling)
    du = vel[1:] - vel[:-1]
    q = np.zeros_like(P_eos)
    compressing = du < 0
    q[compressing] = (VISC_COEFF**2) * rho[compressing] * (du[compressing]**2)
    
    P_total = P_eos + q
    
    # 4. Combustion Trigger
    if not burned and sim_time > T_burn_start:
        burned = True
        centers = 0.5 * (pos[:-1] + pos[1:])
        # Burn only in chamber (0.3 to 0.6)
        mask = (centers > 0.3) & (centers < 0.6)
        entropy[mask] *= Burn_Multiplier
        print(f"COMBUSTION at {sim_time*1000:.1f}ms")

    # 5. Time Step (CFL)
    c_sound = np.sqrt(Gamma * P_total / rho)
    max_vel = np.max(c_sound + np.abs(vel[:-1])) + 1e-5
    current_dt = CFL_SAFETY * np.min(dx) / max_vel
    current_dt = np.clip(current_dt, args.dt_min, args.dt_max)
    sim_time += current_dt
    
    # 6. Forces (Pressure Gradient)
    # Boundary Conditions: Fixed Pressure at ends of huge reservoirs
    P_bc = np.concatenate(([P_atm], P_total, [P_atm]))
    Force = (P_bc[:-1] - P_bc[1:]) * A_nodes
    
    # Node Inertia
    node_mass = np.zeros(len(pos))
    node_mass[1:-1] = 0.5 * (mass[:-1] + mass[1:])
    node_mass[0] = mass[0] * 10 # Heavy boundary anchor
    node_mass[-1] = mass[-1] * 10
    
    # 7. Integration
    accel = Force / node_mass
    
    # --- PHYSICAL REALISM ADDITIONS ---
    
    # A. Wall Friction (Inside Engine Only)
    # f = -C * v * |v| / Diameter
    in_engine = (pos > 0.0) & (pos < 1.6)
    diams = 2 * np.sqrt(A_nodes / np.pi)
    friction = -WALL_FRICTION * vel * np.abs(vel) / (diams + 1e-6)
    accel[in_engine] += friction[in_engine]
    
    # B. Jet Damping (The Valveless Mechanism)
    # If gas is moving rapidly into a wider area (Exhausting), kill KE
    # If gas is moving rapidly into a narrower area (Intaking), keep KE
    # We detect "Exhausting" by checking if we are in the reservoir zones and moving away from engine
    
    in_left_res = (pos < 0.0) & (vel < 0) # Blowing out intake
    in_right_res = (pos > 1.6) & (vel > 0) # Blowing out exhaust
    
    # Apply heavy damping to "Jet" regions to simulate turbulence loss
    accel[in_left_res] -= JET_DAMPING * vel[in_left_res]
    accel[in_right_res] -= JET_DAMPING * vel[in_right_res]

    # Update
    vel += accel * current_dt
    pos += vel * current_dt

def update(frame):
    # Run multiple physics sub-steps per frame
    for _ in range(args.steps_per_frame):
        physics_step()
        
    # Visualization
    centers = 0.5 * (pos[:-1] + pos[1:])
    radii = np.sqrt(np.array([get_area(x) for x in pos]) / np.pi)
    
    # 1. Elements (Color by Temp)
    segs = []
    colors = []
    norm_ent = np.clip((entropy - 1.0)/(Burn_Multiplier - 1.0), 0, 1)
    
    # Only draw lines if they are within view
    view_mask = (pos > -0.6) & (pos < 2.1)
    indices = np.where(view_mask)[0]
    
    for i in indices:
        if i >= len(entropy): continue
        x = pos[i]
        r = radii[i]
        segs.append([[x, -r], [x, r]])
        # Red = Hot, Blue = Cold
        c = norm_ent[i]
        colors.append((c, 0, 1-c, 0.6))
        
    lc.set_segments(segs)
    lc.set_color(colors)
    
    # 2. Pressure
    dx = np.diff(pos)
    dx = np.maximum(dx, 1e-6)
    vol = dx * get_area(centers)
    P = P_atm * entropy * (mass/vol / Rho_atm)**Gamma
    line_p.set_data(centers, P)
    
    # 3. Velocity
    v_mid = 0.5 * (vel[:-1] + vel[1:])
    line_v.set_data(centers, v_mid)
    
    return lc, line_p, line_v

ani = animation.FuncAnimation(fig, update, frames=200, interval=20, blit=True)
plt.show()