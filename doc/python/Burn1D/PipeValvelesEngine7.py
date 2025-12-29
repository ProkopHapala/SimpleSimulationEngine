import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.collections import LineCollection

# --- CLI ARGUMENTS ---
def parse_args():
    p = argparse.ArgumentParser(description="1D pulsejet toy model with combustion and diffusion")
    # Time stepping / numerics
    p.add_argument("--cfl",             type=float, default=0.4,    help="CFL number")
    p.add_argument("--dt-max",          type=float, default=5e-5,   help="Max timestep [s]")
    p.add_argument("--steps-per-frame", type=int,   default=20,     help="Physics steps per animation frame")
    # Grid (constant mass/volume elements)
    p.add_argument("--target-mass",     type=float, default=20e-5,   help="Target mass per element [kg]")
    p.add_argument("--length",          type=float, default=1.0,    help="Pipe length (overridden by geometry end unless set here)")
    # Physics / mixing
    p.add_argument("--friction",        type=float, default=1.5,    help="Wall friction coefficient")
    p.add_argument("--viscosity",       type=float, default=2.0,    help="Artificial viscosity magnitude")
    p.add_argument("--thermal-cond",    type=float, default=0.1,    help="Thermal mixing coefficient")
    p.add_argument("--species-diff",    type=float, default=0.1,    help="Species mixing coefficient")
    # Relax
    p.add_argument("--relax-time",      type=float, default=0.003,  help="Pre-injection relaxation time [s]")
    # Injection / combustion
    p.add_argument("--inject-start",    type=float, default=0.005,  help="Injection start time [s]")
    p.add_argument("--inject-duration", type=float, default=0.010, help="Injection duration [s]")
    p.add_argument("--fuel-mass",       type=float, default=3e-4,   help="Total fuel mass to inject [kg]")
    p.add_argument("--ignite-time",     type=float, default=None,   help="Ignition time [s] (default: end of injection)")
    return p.parse_args()

args = parse_args()
# Derived / compatibility names
args.target_mass = args.target_mass
args.total_fuel_mass = args.fuel_mass
if args.ignite_time is None:
    args.ignite_time = args.inject_start + args.inject_duration

# --- CHEMICAL CONSTANTS ---
# Species Indices: 0:N2, 1:O2, 2:Fuel, 3:CO2, 4:H2O
MOLAR_MASSES = np.array([28.0, 32.0, 58.1, 44.0, 18.0]) * 1e-3 # kg/mol
R_UNIV = 8.314
R_SPECIES = R_UNIV / MOLAR_MASSES # J/(kg K)

# Air Composition (Mass Fraction)
Y_AIR = np.array([0.77, 0.23, 0.0, 0.0, 0.0]) # N2, O2

# Thermo
GAMMA = 1.4

# Stoichiometry (Butane C4H10)
# 2 C4H10 + 13 O2 -> 8 CO2 + 10 H2O
# Mass Balance per 1kg of Fuel:
# M_fuel = 58.1 * 2 = 116.2
# M_ox   = 32.0 * 13 = 416.0
# Ratio Ox/Fuel = 3.58
O2_PER_FUEL = 3.58
# Products per 1kg Fuel:
# CO2 = (8 * 44) / 116.2 = 3.03 kg
# H2O = (10 * 18) / 116.2 = 1.55 kg
Y_PROD_CO2 = 3.03
Y_PROD_H2O = 1.55

HEAT_OF_COMBUSTION = 45e6 # J/kg (Butane LHV approx)

P_ATM = 101325.0
T_ATM = 293.0
# Initial atmospheric density depends on air composition R
R_AIR = np.sum(Y_AIR * R_SPECIES)
RHO_ATM = P_ATM / (R_AIR * T_ATM)

# --- 1. ENGINE GEOMETRY ---
GEOM_NODES = [
    (0.15, 0.030), # Intake
    (0.20, 0.030), 
    (0.30, 0.080), # Chamber Start
    (0.55, 0.080), # Chamber End
    (0.65, 0.045), # Nozzle
    (1.50, 0.080)  # Exhaust Flare
]

# Shift geometry so simulation domain starts at the intake position
geom_x0 = GEOM_NODES[0][0]
geom_xs = [p[0] - geom_x0 for p in GEOM_NODES]
geom_rs = [p[1] for p in GEOM_NODES]
# Override length to span only the defined geometry after shift
args.length = geom_xs[-1]

def get_radius(x):
    return np.interp(x, geom_xs, geom_rs)

def get_area(x):
    return np.pi * get_radius(x)**2

# --- 2. INITIALIZATION ---
def generate_initial_grid():
    # Constant mass per element target
    target_mass = args.target_mass
    xs = [0.0]
    curr_x = 0.0
    while curr_x < args.length:
        area = get_area(curr_x)
        dx = target_mass / (area * RHO_ATM)
        curr_x += dx
        if curr_x >= args.length:
            xs.append(args.length)
            break
        xs.append(curr_x)
    return np.array(xs), target_mass

nodes_x, TARGET_ELEMENT_MASS = generate_initial_grid()
nodes_v = np.zeros(len(nodes_x))

# Species Mass Matrix: [N_elements, 5 species]
N_elem = len(nodes_x) - 1
elem_species = np.zeros((N_elem, 5))

# Initialize with Air
for i in range(N_elem):
    # Calculate volume just to double check, but we assign mass based on target
    # to keep density uniform initially
    elem_species[i, :] = Y_AIR * TARGET_ELEMENT_MASS

elem_entropy = np.ones(N_elem)

sim_time = 0.0
ignited = False
hist_time = []
hist_momentum = []

# --- PLOTTING SETUP ---
fig = plt.figure(figsize=(14, 12))
gs = fig.add_gridspec(5, 2)

# 1. Geometry
ax_geo = fig.add_subplot(gs[0, :])
ax_geo.set_title(f"1D Pulsejet: Composition & Combustion (t=0.00s)")
ax_geo.set_xlim(-0.2, args.length + 0.2)
ax_geo.set_ylim(-0.15, 0.15)
ax_geo.set_aspect('equal')
gx = np.linspace(0, args.length, 300)
gr = np.array([get_radius(x) for x in gx])
ax_geo.plot(gx, gr, 'k', lw=2)
ax_geo.plot(gx, -gr, 'k', lw=2)
lc_elements = LineCollection([], linewidths=0.5) 
lc_centers = LineCollection([], linewidths=3.0)  # temperature-colored center lines
ax_geo.add_collection(lc_elements)
ax_geo.add_collection(lc_centers)

# 2. Composition (New!)
ax_spec = fig.add_subplot(gs[1, :])
ax_spec.set_title("Chemical Composition (Mass Fraction)")
ax_spec.set_xlim(-0.2, args.length + 0.2)
ax_spec.set_ylim(0, 1.0)
ax_spec.grid(True, alpha=0.3)
# Stackplot requires recreating collections usually, so we use lines for speed
l_fuel, = ax_spec.plot([], [], 'orange', label='Fuel', lw=2)
l_n2,   = ax_spec.plot([], [], 'gray', label='N2', lw=1.5)
l_o2,   = ax_spec.plot([], [], 'b', label='O2', lw=1.5)
l_prod, = ax_spec.plot([], [], 'g', label='Products (CO2+H2O)', lw=1.5)
ax_spec.legend(loc='upper right')

# 3. Thermodynamics
ax_thermo = fig.add_subplot(gs[2, :])
ax_thermo.set_xlim(-0.2, args.length + 0.2)
ax_thermo.set_title("Thermodynamics")
line_p, = ax_thermo.plot([], [], 'r-', lw=1, label='Pressure (atm)')
line_t, = ax_thermo.plot([], [], 'm-', lw=1, label='Temp (x300K)')
line_rho, = ax_thermo.plot([], [], 'b-', lw=1, label='Density (kg/m3)')
ax_thermo.legend(loc='upper right')
ax_thermo.grid(True, alpha=0.3)

# 4. Velocity
ax_vel = fig.add_subplot(gs[3, :])
ax_vel.set_xlim(-0.2, args.length + 0.2)
ax_vel.set_ylim(-400, 400)
ax_vel.set_ylabel("Velocity (m/s)")
line_v, = ax_vel.plot([], [], 'k-', lw=1)
ax_vel.grid(True, alpha=0.3)

# 5. Momentum
ax_mom = fig.add_subplot(gs[4, :])
ax_mom.set_title("Total Momentum")
ax_mom.set_xlabel("Time (s)")
line_mom, = ax_mom.plot([], [], 'k-', lw=1.5)
ax_mom.grid(True)

plt.tight_layout()

# --- PHYSICS LOGIC ---

def physics_step():
    global nodes_x, nodes_v, elem_species, elem_entropy, sim_time, ignited
    
    # --- 1. DYNAMIC TOPOLOGY (INFLOW/OUTFLOW) ---
    # Left Boundary (Intake)
    gap_left = nodes_x[0]
    area_intake = get_area(0)
    mass_capacity = gap_left * area_intake * RHO_ATM
    
    if mass_capacity > TARGET_ELEMENT_MASS:
        # Spawn Ambient Air
        dx = TARGET_ELEMENT_MASS / (area_intake * RHO_ATM)
        nodes_x = np.insert(nodes_x, 0, nodes_x[0] - dx)
        nodes_v = np.insert(nodes_v, 0, 0.0)
        
        new_species = Y_AIR * TARGET_ELEMENT_MASS
        elem_species = np.vstack([new_species, elem_species])
        elem_entropy = np.insert(elem_entropy, 0, 1.0)
        
    elif nodes_x[0] < -0.1: # Outflow
        nodes_x, nodes_v = nodes_x[1:], nodes_v[1:]
        elem_species = elem_species[1:]
        elem_entropy = elem_entropy[1:]

    # Right Boundary (Exhaust)
    gap_right = args.length - nodes_x[-1]
    area_ex = get_area(args.length)
    mass_capacity = gap_right * area_ex * RHO_ATM
    
    if mass_capacity > TARGET_ELEMENT_MASS:
        dx = TARGET_ELEMENT_MASS / (area_ex * RHO_ATM)
        nodes_x = np.append(nodes_x, nodes_x[-1] + dx)
        nodes_v = np.append(nodes_v, 0.0)
        
        new_species = Y_AIR * TARGET_ELEMENT_MASS
        elem_species = np.vstack([elem_species, new_species])
        elem_entropy = np.append(elem_entropy, 1.0)
        
    elif nodes_x[-1] > args.length + 0.1:
        nodes_x, nodes_v = nodes_x[:-1], nodes_v[:-1]
        elem_species = elem_species[:-1]
        elem_entropy = elem_entropy[:-1]

    # --- 2. CALCULATE PROPERTIES ---
    dx = np.diff(nodes_x)
    # Anti-collapse
    min_dx = 1e-6
    if np.any(dx < min_dx):
        for i in range(len(dx)):
            if nodes_x[i+1] <= nodes_x[i] + min_dx:
                nodes_x[i+1] = nodes_x[i] + min_dx + 1e-9
        dx = np.diff(nodes_x)

    centers = 0.5 * (nodes_x[:-1] + nodes_x[1:])
    areas = get_area(centers)
    # Conic volumes
    A_L, A_R = get_area(nodes_x[:-1]), get_area(nodes_x[1:])
    vols = (dx/3.0)*(A_L + A_R + np.sqrt(A_L*A_R))
    
    # Mass & Composition
    elem_masses = np.sum(elem_species, axis=1)
    Y_fracs = elem_species / elem_masses[:, None]
    
    rho = elem_masses / vols
    
    # Variable R based on mixture
    # R_mix = sum(Y_i * R_i)
    R_mix = np.sum(Y_fracs * R_SPECIES, axis=1)
    
    # Equation of State
    # P = P_atm * S * (rho/rho_atm)^gamma
    # Note: This implies Gamma is constant. 
    # For a toy model, this is fine. Real combustion engines var gamma 1.4->1.3.
    pressures = P_ATM * elem_entropy * (rho / RHO_ATM)**GAMMA
    temps = pressures / (rho * R_mix)

    # --- 3. INJECTION & COMBUSTION ---
    # Determine adaptive DT first for rate calculations
    c_sound = np.sqrt(GAMMA * pressures / rho)
    max_v = np.max(c_sound + np.abs(nodes_v[1:]))
    dt = args.cfl * np.min(dx) / (max_v + 1e-5)
    dt = np.clip(dt, 1e-9, args.dt_max)
    sim_time += dt

    # A. INJECTION (single-point in chamber center) — only after relax_time
    if sim_time >= args.relax_time and args.inject_start < sim_time < (args.inject_start + args.inject_duration):
        chamber_center = 0.4
        chamber_mask = (centers > 0.3) & (centers < 0.5)
        if np.any(chamber_mask):
            idx = np.argmin(np.abs(centers - chamber_center))
            # Rate of mass injection per second
            mdot = args.total_fuel_mass / args.inject_duration
            d_mass = mdot * dt
            elem_species[idx, 2] += d_mass
            # (Note: adding mass changes rho and P in next step naturally)

    # B. ONE-SHOT IGNITION (at ignite_time) — only after relax_time
    if (not ignited) and sim_time >= max(args.ignite_time, args.relax_time):
        ignited = True
        # Ignite chamber slice; burn stoich-limited mass instantly
        chamber_mask = (centers > 0.3) & (centers < 0.5)
        idxs = np.where(chamber_mask)[0]
        for i in idxs:
            m_fuel = elem_species[i, 2]
            m_ox = elem_species[i, 1]
            burn_amount = min(m_fuel, m_ox / O2_PER_FUEL)
            if burn_amount <= 1e-12:
                continue
            elem_species[i, 2] -= burn_amount
            elem_species[i, 1] -= burn_amount * O2_PER_FUEL
            elem_species[i, 3] += burn_amount * Y_PROD_CO2
            elem_species[i, 4] += burn_amount * Y_PROD_H2O
            dQ = burn_amount * HEAT_OF_COMBUSTION
            dS = dQ / temps[i]
            elem_entropy[i] += dS * 0.8

    # --- 4. DIFFUSION (Mass & Heat) ---
    # Mixing between neighbors
    
    # A. Heat Diffusion
    # Nodes connect elements i-1 and i. 
    dT = np.zeros(len(nodes_x))
    # Pad temps with boundaries
    T_bound = np.concatenate(([T_ATM], temps, [T_ATM]))
    dx_nodes = np.concatenate(([dx[0]], (dx[:-1]+dx[1:])*0.5, [dx[-1]]))
    dT = T_bound[1:] - T_bound[:-1]
    
    heat_flux = -args.thermal_cond * get_area(nodes_x) * dT / dx_nodes
    # dS = (Flux_in - Flux_out)*dt / T
    net_heat = (heat_flux[:-1] - heat_flux[1:]) * dt
    elem_entropy += net_heat / (elem_masses * temps)
    
    # B. Species Diffusion (Fick's Law)
    # Flux J = -D * A * dC/dx
    # Concentration C = Mass_Fraction (simplification)
    # We iterate over 5 species
    for sp in range(5):
        Y = elem_species[:, sp] / elem_masses
        Y_bound = np.concatenate(([Y_AIR[sp]], Y, [Y_AIR[sp]])) # Boundary is air
        dY = Y_bound[1:] - Y_bound[:-1]
        
        mass_flux = -args.species_diff * get_area(nodes_x) * dY / dx_nodes * 10.0 # Multiplier for visible effect
        
        # Apply flux (Mass transfer)
        net_mass_flow = (mass_flux[:-1] - mass_flux[1:]) * dt
        elem_species[:, sp] += net_mass_flow

    # --- 5. MECHANICAL INTEGRATION ---
    
    # Artificial Viscosity
    du = nodes_v[1:] - nodes_v[:-1]
    q = np.zeros_like(pressures)
    compressing = du < 0
    q[compressing] = (args.viscosity**2) * rho[compressing] * (du[compressing]**2)
    
    P_total = pressures + q
    
    # Forces
    P_bc = np.concatenate(([P_ATM], P_total, [P_ATM]))
    delta_P = P_bc[:-1] - P_bc[1:]
    force = delta_P * get_area(nodes_x)
    
    # Node Inertia
    node_m = np.zeros(len(nodes_x))
    node_m[1:-1] = 0.5 * (elem_masses[:-1] + elem_masses[1:])
    node_m[0] = elem_masses[0]
    node_m[-1] = elem_masses[-1]
    
    friction = -args.friction * nodes_v * np.abs(nodes_v) * node_m
    force += friction
    
    accel = force / node_m
    nodes_v += accel * dt
    nodes_x += nodes_v * dt
    
    # Diagnostics
    total_mom = np.sum(node_m * nodes_v)
    hist_time.append(sim_time)
    hist_momentum.append(total_mom)

# --- ANIMATION ---
def update_plot(frame):
    for _ in range(args.steps_per_frame):
        physics_step()
        
    centers = 0.5 * (nodes_x[:-1] + nodes_x[1:])
    
    # Recalc thermo for plotting
    m_t = np.sum(elem_species, axis=1)
    v_t = (np.diff(nodes_x)/3.0)*(get_area(nodes_x[:-1]) + get_area(nodes_x[1:]) + np.sqrt(get_area(nodes_x[:-1])*get_area(nodes_x[1:])))
    rho = m_t/v_t
    Y = elem_species / m_t[:, None]
    R_mix = np.sum(Y * R_SPECIES, axis=1)
    p = P_ATM * elem_entropy * (rho/RHO_ATM)**GAMMA
    t = p / (rho * R_mix)
    
    # 1. Elements (Color by Temp)
    segs = []
    colors = []
    center_segs = []
    center_colors = []
    norm_t = np.clip((t - 300) / 400, 0, 1)
    radii = get_radius(nodes_x)
    for i in range(len(nodes_x)):
        x = nodes_x[i]
        r = radii[i]
        segs.append([[x, -r], [x, r]])
        # Color nodes
        if i < len(t): val = norm_t[i]
        elif i > 0: val = norm_t[i-1]
        else: val = 0
        colors.append((val, 0, 1-val, 1)) # Blue -> Red
        # if i < len(centers):
        #     cx = centers[i]
        #     cr = get_radius(cx)
        #     center_segs.append([[cx, -cr], [cx, cr]])
        #     center_colors.append((val, 0, 1-val, 1))
    lc_elements.set_segments(segs)
    lc_elements.set_color(colors)
    # lc_centers.set_segments(center_segs)
    # lc_centers.set_color(center_colors)
    ax_geo.set_title(f"Time: {sim_time*1000:.1f} ms | Elements: {len(nodes_x)}")
    
    # 2. Composition (Mass Fraction)
    # 0:N2, 1:O2, 2:Fuel, 3:CO2, 4:H2O
    # l_n2.set_data(centers, Y[:, 0])
    l_fuel.set_data(centers, Y[:, 2])
    l_o2.set_data(centers, Y[:, 1])
    l_prod.set_data(centers, Y[:, 3] + Y[:, 4])
    
    # 3. Thermo
    line_p.set_data(centers, p / P_ATM)
    line_t.set_data(centers, t / 300.0)
    line_rho.set_data(centers, rho)
    # Autoscale thermo axes to keep traces visible
    y_min = np.min([np.min(p / P_ATM), np.min(t / 300.0), np.min(rho)])
    y_max = np.max([np.max(p / P_ATM), np.max(t / 300.0), np.max(rho)])
    pad = 0.1 * max(1e-3, y_max - y_min)
    ax_thermo.set_ylim(y_min - pad, y_max + pad)
    
    # 4. Velocity
    v_elem = 0.5 * (nodes_v[:-1] + nodes_v[1:])
    line_v.set_data(centers, v_elem)
    
    # 5. Momentum
    if len(hist_time) > 0:
        line_mom.set_data(hist_time, hist_momentum)
        ax_mom.set_xlim(0, max(hist_time[-1], 0.05))
        vals = np.array(hist_momentum)
        if len(vals) > 5:
            ax_mom.set_ylim(np.min(vals), np.max(vals))

    return lc_elements, lc_centers, line_p

ani = animation.FuncAnimation(fig, update_plot, interval=10, blit=False)
plt.show()