import numpy as np

GAMMA = 1.4
P_ATM = 101325.0
RHO_ATM = 1.225

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

def get_area(get_radius_fn, x):
    return np.pi * get_radius_fn(x)**2

def generate_grid(get_radius_fn, length, n_cells):
    xs = [0.0]
    curr_x = 0.0
    avg_area = 0.5 * (get_area(get_radius_fn, 0) + get_area(get_radius_fn, length))
    target_dx_approx = length / n_cells
    target_mass = avg_area * target_dx_approx * RHO_ATM
    while curr_x < length:
        area_loc = get_area(get_radius_fn, curr_x)
        dx_needed = target_mass / (RHO_ATM * area_loc)
        curr_x += dx_needed
        if curr_x >= length:
            xs.append(length)
            break
        xs.append(curr_x)
    nodes_x = np.array(xs)
    dxs = np.diff(nodes_x)
    A1 = get_area(get_radius_fn, nodes_x[:-1])
    A2 = get_area(get_radius_fn, nodes_x[1:])
    vols = (dxs / 3.0) * (A1 + A2 + np.sqrt(A1 * A2))
    cell_mass = vols * RHO_ATM
    node_mass = np.zeros(len(nodes_x))
    node_mass[1:-1] = 0.5 * (cell_mass[:-1] + cell_mass[1:])
    node_mass[0] = 0.5 * cell_mass[0]
    node_mass[-1] = 0.5 * cell_mass[-1]
    return nodes_x, cell_mass, node_mass, target_mass

class WaveTubeSolver:
    def __init__(self, geo="nozzle", length=2.0, n_cells=1000, kick_vel=50.0, dt=2e-6):
        self.geo = geo
        self.length = length
        self.n_cells = n_cells
        self.kick_vel = kick_vel
        self.dt = dt
        self.get_radius = make_get_radius(geo, length)
        self.nodes_x, self.cell_mass, self.node_mass, self.target_mass = generate_grid(
            self.get_radius, length, n_cells)
        center_x = length / 2.0
        self.nodes_v = kick_vel * np.exp(-(self.nodes_x - center_x)**2 / (2 * 0.05**2))
        self.nodes_f = np.zeros_like(self.nodes_v)
        self.sim_t = 0.0
        self.hist_t = []
        self.hist_mom = []
        self.hist_E = []
        self.nodes_f, _ = self._get_forces(self.nodes_x)

    def _get_forces(self, x):
        dx = x[1:] - x[:-1]
        A_nodes = get_area(self.get_radius, x)
        A_L = A_nodes[:-1]
        A_R = A_nodes[1:]
        vol = (dx / 3.0) * (A_L + A_R + np.sqrt(A_L * A_R))
        rho = self.cell_mass / vol
        p = P_ATM * (rho / RHO_ATM)**GAMMA
        f = np.zeros(len(x))
        delta_p = p[:-1] - p[1:]
        f[1:-1] = delta_p * A_nodes[1:-1]
        return f, p

    def step(self):
        dt = self.dt
        self.nodes_v += 0.5 * (self.nodes_f / self.node_mass) * dt
        self.nodes_v[0] = 0.0
        self.nodes_v[-1] = 0.0
        self.nodes_x += self.nodes_v * dt
        self.nodes_f, p = self._get_forces(self.nodes_x)
        self.nodes_v += 0.5 * (self.nodes_f / self.node_mass) * dt
        self.nodes_v[0] = 0.0
        self.nodes_v[-1] = 0.0
        self.sim_t += dt
        mom = np.sum(self.node_mass * self.nodes_v)
        ke = 0.5 * np.sum(self.node_mass * self.nodes_v**2)
        dx = self.nodes_x[1:] - self.nodes_x[:-1]
        A_nodes = get_area(self.get_radius, self.nodes_x)
        vol = (dx / 3.0) * (A_nodes[:-1] + A_nodes[1:] + np.sqrt(A_nodes[:-1] * A_nodes[1:]))
        rho = self.cell_mass / vol
        p_exact = P_ATM * (rho / RHO_ATM)**GAMMA
        ie = np.sum(p_exact * vol / (GAMMA - 1))
        self.hist_t.append(self.sim_t)
        self.hist_mom.append(mom)
        self.hist_E.append(ke + ie)
        return p

    def get_centers(self):
        return 0.5 * (self.nodes_x[1:] + self.nodes_x[:-1])

    def get_vel_elements(self):
        return 0.5 * (self.nodes_v[:-1] + self.nodes_v[1:])
