import numpy as np

GAMMA = 1.4
P_ATM = 101325.0
T_ATM = 293.0

MOLAR_MASSES = np.array([28.0, 32.0, 58.1, 44.0, 18.0]) * 1e-3
R_UNIV = 8.314
R_SPECIES = R_UNIV / MOLAR_MASSES
Y_AIR = np.array([0.77, 0.23, 0.0, 0.0, 0.0])
R_AIR = np.sum(Y_AIR * R_SPECIES)
RHO_ATM = P_ATM / (R_AIR * T_ATM)

O2_PER_FUEL = 3.58
Y_PROD_CO2 = 3.03
Y_PROD_H2O = 1.55
HEAT_OF_COMBUSTION = 45e6

GEOM_NODES_ENGINE = [
    (-0.30, 0.030), (-0.20, 0.030), (-0.10, 0.080),
    (0.00, 0.080), (0.15, 0.045), (1.00, 0.080)
]
GEOM_NODES_INERTIA = [(0.0, 0.03), (2.0, 0.03)]

def make_geometry(mode="ENGINE"):
    if mode == "INERTIA":
        nodes = GEOM_NODES_INERTIA
    else:
        nodes = GEOM_NODES_ENGINE
    geom_x0 = nodes[0][0]
    geom_xs = [p[0] - geom_x0 for p in nodes]
    geom_rs = [p[1] for p in nodes]
    length = geom_xs[-1]
    def get_radius(x):
        x_clamped = np.clip(np.asarray(x), 0, length if mode == "ENGINE" else 999)
        return np.interp(x_clamped, geom_xs, geom_rs)
    def get_area(x):
        return np.pi * get_radius(x)**2
    return get_radius, get_area, length, geom_x0

class PulsejetSolver:
    def __init__(self, mode="ENGINE", target_mass=20e-5, cfl=0.4, dt_max=5e-5,
                 friction=0.1, viscosity=0.1, thermal_cond=0.1, species_diff=0.1,
                 relax_time=0.003, inject_period=0.020, inject_start=0.005,
                 inject_duration=0.005, inject_point=0.0, fuel_mass=3e-4,
                 ignite_delay=0.001, burn_time=0.02, burn_mult=4.0):
        self.mode = mode
        self.target_mass = target_mass
        self.cfl = cfl
        self.dt_max = dt_max
        self.friction = friction
        self.viscosity = viscosity
        self.thermal_cond = thermal_cond
        self.species_diff = species_diff
        self.relax_time = relax_time
        self.inject_period = inject_period
        self.inject_start = inject_start
        self.inject_duration = inject_duration
        self.inject_point = inject_point
        self.fuel_mass = fuel_mass
        self.ignite_delay = ignite_delay
        self.burn_time = burn_time
        self.burn_mult = burn_mult

        self.get_radius, self.get_area, self.length, self.geom_x0 = make_geometry(mode)
        self._init_grid()
        self.sim_time = 0.0
        self.last_ignite_cycle = -1
        self.burned = False
        self.hist_time = []
        self.hist_momentum = []

    def _init_grid(self):
        xs = [0.0]
        curr_x = 0.0
        while curr_x < self.length:
            area = self.get_area(curr_x)
            dx = self.target_mass / (area * RHO_ATM)
            curr_x += dx
            if curr_x >= self.length:
                xs.append(self.length)
                break
            xs.append(curr_x)
        self.nodes_x = np.array(xs)
        self.nodes_v = np.zeros(len(self.nodes_x))
        if self.mode == "INERTIA":
            self.nodes_v[:] = 200.0
        self.N_elem = len(self.nodes_x) - 1
        if self.mode == "ENGINE":
            self.elem_species = np.zeros((self.N_elem, 5))
            for i in range(self.N_elem):
                self.elem_species[i, :] = Y_AIR * self.target_mass
        else:
            self.elem_species = None
        self.elem_entropy = np.ones(self.N_elem)

    def step(self):
        if self.mode == "ENGINE" and self.elem_species is not None:
            self._step_engine()
        else:
            self._step_inertia()

    def _step_engine(self):
        self._dynamic_topology_engine()
        dx = np.diff(self.nodes_x)
        min_dx = 1e-6
        if np.any(dx < min_dx):
            for i in range(len(dx)):
                if self.nodes_x[i + 1] <= self.nodes_x[i] + min_dx:
                    self.nodes_x[i + 1] = self.nodes_x[i] + min_dx + 1e-9
            dx = np.diff(self.nodes_x)
        centers = 0.5 * (self.nodes_x[:-1] + self.nodes_x[1:])
        A_L = self.get_area(self.nodes_x[:-1])
        A_R = self.get_area(self.nodes_x[1:])
        vols = (dx / 3.0) * (A_L + A_R + np.sqrt(A_L * A_R))
        elem_masses = np.sum(self.elem_species, axis=1)
        Y_fracs = self.elem_species / elem_masses[:, None]
        rho = elem_masses / vols
        R_mix = np.sum(Y_fracs * R_SPECIES, axis=1)
        pressures = P_ATM * self.elem_entropy * (rho / RHO_ATM) ** GAMMA
        temps = pressures / (rho * R_mix)
        c_sound = np.sqrt(GAMMA * pressures / rho)
        max_v = np.max(c_sound + np.abs(self.nodes_v[1:]))
        dt = self.cfl * np.min(dx) / (max_v + 1e-5)
        dt = np.clip(dt, 1e-9, self.dt_max)
        self.sim_time += dt
        if self.sim_time >= self.relax_time and self.sim_time >= self.inject_start:
            phase = (self.sim_time - self.inject_start) % self.inject_period
            cycle_idx = int((self.sim_time - self.inject_start) // self.inject_period)
            if 0 <= phase < self.inject_duration:
                inj_x = self.inject_point - self.geom_x0
                idx = np.argmin(np.abs(centers - inj_x))
                self.elem_species[idx, 2] += (self.fuel_mass / self.inject_duration) * dt
            ign_start = self.inject_duration + self.ignite_delay
            if ign_start <= phase < ign_start + dt * 1.1 and cycle_idx != self.last_ignite_cycle:
                self.last_ignite_cycle = cycle_idx
                inj_x = self.inject_point - self.geom_x0
                chamber_mask = np.abs(centers - inj_x) < 0.05
                for i in np.where(chamber_mask)[0]:
                    burn = min(self.elem_species[i, 2], self.elem_species[i, 1] / O2_PER_FUEL)
                    if burn > 1e-12:
                        self.elem_species[i, 2] -= burn
                        self.elem_species[i, 1] -= burn * O2_PER_FUEL
                        self.elem_species[i, 3] += burn * Y_PROD_CO2
                        self.elem_species[i, 4] += burn * Y_PROD_H2O
                        self.elem_entropy[i] += (burn * HEAT_OF_COMBUSTION / temps[i]) * 0.8
        dT_bound = np.concatenate(([T_ATM], temps, [T_ATM]))
        dx_nodes = np.concatenate(([dx[0]], (dx[:-1] + dx[1:]) * 0.5, [dx[-1]]))
        dT = dT_bound[1:] - dT_bound[:-1]
        heat_flux = -self.thermal_cond * self.get_area(self.nodes_x) * dT / dx_nodes
        net_heat = (heat_flux[:-1] - heat_flux[1:]) * dt
        self.elem_entropy += net_heat / (elem_masses * temps)
        for sp in range(5):
            Y = self.elem_species[:, sp] / elem_masses
            Y_bound = np.concatenate(([Y_AIR[sp]], Y, [Y_AIR[sp]]))
            dY = Y_bound[1:] - Y_bound[:-1]
            mass_flux = -self.species_diff * self.get_area(self.nodes_x) * dY / dx_nodes * 10.0
            net_mass_flow = (mass_flux[:-1] - mass_flux[1:]) * dt
            self.elem_species[:, sp] += net_mass_flow
        du = self.nodes_v[1:] - self.nodes_v[:-1]
        q = np.zeros_like(pressures)
        compressing = du < 0
        q[compressing] = (self.viscosity ** 2) * rho[compressing] * (du[compressing] ** 2)
        P_total = pressures + q
        P_bc = np.concatenate(([P_ATM], P_total, [P_ATM]))
        delta_P = P_bc[:-1] - P_bc[1:]
        force = delta_P * self.get_area(self.nodes_x)
        node_m = np.zeros(len(self.nodes_x))
        node_m[1:-1] = 0.5 * (elem_masses[:-1] + elem_masses[1:])
        node_m[0] = elem_masses[0]
        node_m[-1] = elem_masses[-1]
        force += -self.friction * self.nodes_v * np.abs(self.nodes_v) * node_m
        accel = force / node_m
        self.nodes_v += accel * dt
        self.nodes_x += self.nodes_v * dt
        total_mom = np.sum(node_m * self.nodes_v)
        self.hist_time.append(self.sim_time)
        self.hist_momentum.append(total_mom)

    def _step_inertia(self):
        A_in = self.get_area(self.nodes_x[0])
        spawn_threshold = self.target_mass / (A_in * RHO_ATM)
        if self.nodes_x[0] > spawn_threshold:
            dx_new = spawn_threshold
            x_new = self.nodes_x[0] - dx_new
            v_new = self.nodes_v[0]
            self.nodes_x = np.insert(self.nodes_x, 0, x_new)
            self.nodes_v = np.insert(self.nodes_v, 0, v_new)
            self.elem_entropy = np.insert(self.elem_entropy, 0, 1.0)
        elif self.nodes_x[0] < -0.2:
            self.nodes_x = self.nodes_x[1:]
            self.nodes_v = self.nodes_v[1:]
            self.elem_entropy = self.elem_entropy[1:]
        A_out = self.get_area(self.nodes_x[-1])
        spawn_threshold_r = self.target_mass / (A_out * RHO_ATM)
        if self.nodes_x[-1] < self.length - spawn_threshold_r:
            x_new = self.nodes_x[-1] + spawn_threshold_r
            v_new = self.nodes_v[-1]
            self.nodes_x = np.append(self.nodes_x, x_new)
            self.nodes_v = np.append(self.nodes_v, v_new)
            self.elem_entropy = np.append(self.elem_entropy, 1.0)
        elif self.nodes_x[-1] > self.length + 0.2:
            self.nodes_x = self.nodes_x[:-1]
            self.nodes_v = self.nodes_v[:-1]
            self.elem_entropy = self.elem_entropy[:-1]
        dx = np.diff(self.nodes_x)
        min_dx = 1e-5
        if np.any(dx < min_dx):
            dx = np.maximum(dx, min_dx)
            self.nodes_x[1:] = self.nodes_x[0] + np.cumsum(dx)
        centers = 0.5 * (self.nodes_x[:-1] + self.nodes_x[1:])
        A_L = self.get_area(self.nodes_x[:-1])
        A_R = self.get_area(self.nodes_x[1:])
        vols = (dx / 3.0) * (A_L + A_R + np.sqrt(A_L * A_R))
        n_cells = len(dx)
        mass = np.full(n_cells, self.target_mass)
        rho = mass / vols
        P_thermo = P_ATM * self.elem_entropy * (rho / RHO_ATM) ** GAMMA
        du = self.nodes_v[1:] - self.nodes_v[:-1]
        q = np.zeros_like(P_thermo)
        if self.viscosity > 0:
            compressing = du < 0
            q[compressing] = (self.viscosity ** 2) * rho[compressing] * (du[compressing] ** 2)
        P_tot = P_thermo + q
        if self.mode == "INERTIA":
            P_ghost_L = P_tot[0]
            P_ghost_R = P_tot[-1]
        else:
            P_ghost_L = P_ATM
            P_ghost_R = P_ATM
        P_bcs = np.concatenate(([P_ghost_L], P_tot, [P_ghost_R]))
        delta_P = P_bcs[:-1] - P_bcs[1:]
        A_nodes = self.get_area(self.nodes_x)
        Force = delta_P * A_nodes
        node_m = np.zeros(len(self.nodes_x))
        node_m[1:-1] = 0.5 * (mass[:-1] + mass[1:])
        node_m[0] = 0.5 * mass[0]
        node_m[-1] = 0.5 * mass[-1]
        if self.friction > 0:
            Force -= self.friction * self.nodes_v * np.abs(self.nodes_v) * node_m
        accel = Force / node_m
        c_s = np.sqrt(GAMMA * P_tot / rho)
        max_v = np.max(c_s + np.abs(self.nodes_v[1:]))
        dt = self.cfl * np.min(dx) / (max_v + 1e-5)
        dt = np.clip(dt, 1e-9, self.dt_max)
        if self.mode == "ENGINE":
            self.sim_time += dt
            if not self.burned and self.sim_time > self.burn_time:
                self.burned = True
                mask = (centers > 0.35) & (centers < 0.55)
                self.elem_entropy[mask] *= self.burn_mult
        self.nodes_v += accel * dt
        self.nodes_x += self.nodes_v * dt
        total_mom = np.sum(node_m * self.nodes_v)
        self.hist_time.append(self.sim_time)
        self.hist_momentum.append(total_mom)

    def _dynamic_topology_engine(self):
        gap_left = self.nodes_x[0]
        area_intake = self.get_area(0)
        mass_capacity = gap_left * area_intake * RHO_ATM
        if mass_capacity > self.target_mass:
            dx = self.target_mass / (area_intake * RHO_ATM)
            self.nodes_x = np.insert(self.nodes_x, 0, self.nodes_x[0] - dx)
            self.nodes_v = np.insert(self.nodes_v, 0, 0.0)
            new_species = Y_AIR * self.target_mass
            self.elem_species = np.vstack([new_species, self.elem_species])
            self.elem_entropy = np.insert(self.elem_entropy, 0, 1.0)
        elif self.nodes_x[0] < -0.1:
            self.nodes_x = self.nodes_x[1:]
            self.nodes_v = self.nodes_v[1:]
            self.elem_species = self.elem_species[1:]
            self.elem_entropy = self.elem_entropy[1:]
        gap_right = self.length - self.nodes_x[-1]
        area_ex = self.get_area(self.length)
        mass_capacity = gap_right * area_ex * RHO_ATM
        if mass_capacity > self.target_mass:
            dx = self.target_mass / (area_ex * RHO_ATM)
            self.nodes_x = np.append(self.nodes_x, self.nodes_x[-1] + dx)
            self.nodes_v = np.append(self.nodes_v, 0.0)
            new_species = Y_AIR * self.target_mass
            self.elem_species = np.vstack([self.elem_species, new_species])
            self.elem_entropy = np.append(self.elem_entropy, 1.0)
        elif self.nodes_x[-1] > self.length + 0.1:
            self.nodes_x = self.nodes_x[:-1]
            self.nodes_v = self.nodes_v[:-1]
            self.elem_species = self.elem_species[:-1]
            self.elem_entropy = self.elem_entropy[:-1]

    def get_state(self):
        dx = np.diff(self.nodes_x)
        centers = 0.5 * (self.nodes_x[:-1] + self.nodes_x[1:])
        A_L = self.get_area(self.nodes_x[:-1])
        A_R = self.get_area(self.nodes_x[1:])
        vols = (dx / 3.0) * (A_L + A_R + np.sqrt(A_L * A_R))
        if self.elem_species is not None:
            m_t = np.sum(self.elem_species, axis=1)
        else:
            m_t = np.full(len(dx), self.target_mass)
        rho = m_t / vols
        if self.elem_species is not None:
            Y = self.elem_species / m_t[:, None]
            R_mix = np.sum(Y * R_SPECIES, axis=1)
        else:
            Y = None
            R_mix = np.full(len(dx), R_AIR)
        p = P_ATM * self.elem_entropy * (rho / RHO_ATM) ** GAMMA
        t = p / (rho * R_mix)
        return centers, rho, p, t, Y, vols, m_t
