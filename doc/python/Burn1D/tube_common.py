import numpy as np

# Shared physical constants
GAMMA = 1.4
P_ATM = 101325.0
RHO_ATM = 1.225
T_ATM = 293.0
R_GAS = 287.0

# --- Geometry helpers ---

def get_area_from_radius(get_radius_fn, x):
    return np.pi * get_radius_fn(x)**2

def get_radius_nozzle(x, length):
    x_arr = np.clip(np.asarray(x), 0.0, length)
    base_r = 0.05
    t = np.clip((x_arr - 1.4) / 0.6, 0.0, 1.0)
    return base_r - 0.03 * (1 - np.cos(t * np.pi / 2.0))

def get_radius_cone(x, length):
    x_arr = np.clip(np.asarray(x), 0.0, length)
    r_start, r_end = 0.02, 0.08
    t = x_arr / length
    return r_start + (r_end - r_start) * t

def make_get_radius(geo, length):
    if geo == "nozzle":
        return lambda x: get_radius_nozzle(x, length)
    elif geo == "cone":
        return lambda x: get_radius_cone(x, length)
    else:
        raise ValueError(f"Unknown geometry: {geo}")

def make_interp_geometry(geom_nodes):
    """Linear interpolation geometry from list of (x, r) nodes."""
    x0 = geom_nodes[0][0]
    xs = [p[0] - x0 for p in geom_nodes]
    rs = [p[1] for p in geom_nodes]
    length = xs[-1]
    def get_radius(x):
        return np.interp(np.clip(np.asarray(x), 0, length), xs, rs)
    def get_area(x):
        return np.pi * get_radius(x)**2
    return get_radius, get_area, length, x0

# --- Grid generation ---

def generate_equal_mass_grid(get_area_fn, length, n_cells=None, target_mass=None):
    """Generate equal-mass Lagrangian grid. Specify either n_cells or target_mass."""
    if target_mass is None:
        avg_area = 0.5 * (get_area_fn(0) + get_area_fn(length))
        target_mass = avg_area * (length / n_cells) * RHO_ATM
    xs = [0.0]
    curr_x = 0.0
    while curr_x < length:
        area_loc = get_area_fn(curr_x)
        dx_needed = target_mass / (RHO_ATM * area_loc)
        curr_x += dx_needed
        if curr_x >= length:
            xs.append(length)
            break
        xs.append(curr_x)
    return np.array(xs), target_mass

# --- Volume and mass ---

def conic_frustum_vol(dx, A_L, A_R):
    return (dx / 3.0) * (A_L + A_R + np.sqrt(A_L * A_R))

def node_mass_from_cell_mass(cell_mass):
    node_m = np.zeros(len(cell_mass) + 1)
    node_m[1:-1] = 0.5 * (cell_mass[:-1] + cell_mass[1:])
    node_m[0] = 0.5 * cell_mass[0]
    node_m[-1] = 0.5 * cell_mass[-1]
    return node_m

# --- CFL ---

def cfl_dt(dx, c_sound, v, cfl, dt_max=1e9):
    max_v = np.max(c_sound + np.abs(v))
    dt = cfl * np.min(dx) / (max_v + 1e-5)
    return np.clip(dt, 1e-9, dt_max)
