import numpy as np
from tube_common import (GAMMA, P_ATM, RHO_ATM, T_ATM,
    make_get_radius, make_interp_geometry, get_area_from_radius,
    generate_equal_mass_grid, conic_frustum_vol, node_mass_from_cell_mass, cfl_dt)

# --- Species constants (pulsejet) ---
MOLAR_MASSES = np.array([28.0, 32.0, 58.1, 44.0, 18.0]) * 1e-3
R_UNIV = 8.314
R_SPECIES = R_UNIV / MOLAR_MASSES
Y_AIR = np.array([0.77, 0.23, 0.0, 0.0, 0.0])
R_AIR = np.sum(Y_AIR * R_SPECIES)
O2_PER_FUEL = 3.58
Y_PROD_CO2 = 3.03
Y_PROD_H2O = 1.55
HEAT_OF_COMBUSTION = 45e6

GEOM_NODES_ENGINE = [
    (-1.0, 0.50),    # Left reservoir (far)
    (-0.40, 0.20),   # Left funnel
    (-0.30, 0.030),  # Intake start (engine)
    (-0.20, 0.030),  # Intake end
    (-0.10, 0.080),  # Chamber start
    (0.00, 0.080),   # Chamber end
    (0.15, 0.045),   # Throat
    (1.00, 0.080),   # Exhaust end (engine)
    (1.20, 0.20),    # Right funnel
    (2.0, 0.50)      # Right reservoir (far)
]
ENGINE_START_X = -0.30
ENGINE_END_X = 1.00
GEOM_NODES_INERTIA = [(0.0, 0.03), (2.0, 0.03)]

# ===========================================================================
# Base class: shared Lagrangian machinery
# ===========================================================================

class LagrangianTube1D:
    """Base class for 1D Lagrangian tube solvers (moving mesh, constant mass per element)."""

    def compute_pressure(self, rho, entropy):
        p = P_ATM * entropy * (rho / RHO_ATM) ** GAMMA
        return np.maximum(p, 100.0)  # floor to prevent NaN in c_sound

    def compute_forces(self, nodes_x, P_total, P_ghost_L, P_ghost_R):
        P_bcs = np.concatenate(([P_ghost_L], P_total, [P_ghost_R]))
        delta_P = P_bcs[:-1] - P_bcs[1:]
        A_nodes = self.get_area(nodes_x)
        return delta_P * A_nodes

    def compute_wall_forces(self, nodes_x, P_total, P_ghost_L, P_ghost_R):
        """Wall force from pressure on changing cross-section: 0.5 * P * dA per node.
        This is the discrete equivalent of the p*dA/dx source term in Eulerian framework.
        Without this, momentum is not conserved through variable-area sections."""
        A_nodes = self.get_area(nodes_x)
        dA = np.diff(A_nodes)
        f_wall = np.zeros(len(nodes_x))
        f_wall[:-1] += 0.5 * P_total * dA
        f_wall[1:]  += 0.5 * P_total * dA
        return f_wall

    def artificial_viscosity(self, du, rho, viscosity):
        q = np.zeros_like(du)
        if viscosity > 0:
            compressing = du < 0
            q[compressing] = (viscosity ** 2) * rho[compressing] * (du[compressing] ** 2)
        return q

    def get_centers(self):
        return 0.5 * (self.nodes_x[1:] + self.nodes_x[:-1])

    def get_vel_elements(self):
        return 0.5 * (self.nodes_v[:-1] + self.nodes_v[1:])

    def get_volumes(self):
        dx = np.diff(self.nodes_x)
        A_L = self.get_area(self.nodes_x[:-1])
        A_R = self.get_area(self.nodes_x[1:])
        return conic_frustum_vol(dx, A_L, A_R), dx

# ===========================================================================
# WaveTubeSolver — isentropic acoustic wave, velocity Verlet, fixed walls
# ===========================================================================

class WaveTubeSolver(LagrangianTube1D):
    def __init__(self, geo="nozzle", length=2.0, n_cells=1000, kick_vel=50.0, dt=2e-6):
        self.geo = geo
        self.length = length
        self.n_cells = n_cells
        self.kick_vel = kick_vel
        self.dt = dt
        self.get_radius = make_get_radius(geo, length)
        self.get_area = lambda x: get_area_from_radius(self.get_radius, x)
        self.nodes_x, self.target_mass = generate_equal_mass_grid(self.get_area, length, n_cells=n_cells)
        dxs = np.diff(self.nodes_x)
        A1 = self.get_area(self.nodes_x[:-1])
        A2 = self.get_area(self.nodes_x[1:])
        vols = conic_frustum_vol(dxs, A1, A2)
        self.cell_mass = vols * RHO_ATM
        self.node_mass = node_mass_from_cell_mass(self.cell_mass)
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
        A_nodes = self.get_area(x)
        A_L = A_nodes[:-1]
        A_R = A_nodes[1:]
        vol = conic_frustum_vol(dx, A_L, A_R)
        rho = self.cell_mass / vol
        p = P_ATM * (rho / RHO_ATM)**GAMMA  # isentropic (entropy=1)
        f = np.zeros(len(x))
        delta_p = p[:-1] - p[1:]
        f[1:-1] = delta_p * A_nodes[1:-1]
        return f, p

    def step(self):
        dt = self.dt
        self.nodes_v += 0.5 * (self.nodes_f / self.node_mass) * dt
        self.nodes_v[0] = 0.0; self.nodes_v[-1] = 0.0
        self.nodes_x += self.nodes_v * dt
        self.nodes_f, p = self._get_forces(self.nodes_x)
        self.nodes_v += 0.5 * (self.nodes_f / self.node_mass) * dt
        self.nodes_v[0] = 0.0; self.nodes_v[-1] = 0.0
        self.sim_t += dt
        mom = np.sum(self.node_mass * self.nodes_v)
        ke = 0.5 * np.sum(self.node_mass * self.nodes_v**2)
        vols, dx = self.get_volumes()
        rho = self.cell_mass / vols
        p_exact = P_ATM * (rho / RHO_ATM)**GAMMA
        ie = np.sum(p_exact * vols / (GAMMA - 1))
        self.hist_t.append(self.sim_t)
        self.hist_mom.append(mom)
        self.hist_E.append(ke + ie)
        return p

# ===========================================================================
# PulsejetSolver — entropy EOS, species, combustion, dynamic topology
# ===========================================================================

class PulsejetSolver(LagrangianTube1D):
    def __init__(self, mode="ENGINE", target_mass=200e-5, cfl=0.4, dt_max=5e-5,
                 friction=0.1, viscosity=0.1, thermal_cond=0.1, species_diff=0.1,
                 relax_time=0.003, inject_period=0.020, inject_start=0.005,
                 inject_duration=0.005, inject_point=0.0, fuel_mass=3e-5,
                 ignite_delay=0.001, burn_time=0.02, burn_mult=4.0,
                 jet_damping=5.0):
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
        self.jet_damping = jet_damping
        self.auto_inject = True
        self.manual_inject_x = None
        geom_nodes = GEOM_NODES_INERTIA if mode == "INERTIA" else GEOM_NODES_ENGINE
        self.get_radius, self.get_area, self.length, self.geom_x0 = make_interp_geometry(geom_nodes)
        if mode == "ENGINE":
            self.engine_start_local = ENGINE_START_X - self.geom_x0
            self.engine_end_local = ENGINE_END_X - self.geom_x0
        else:
            self.engine_start_local = 0.0
            self.engine_end_local = self.length
        self._init_grid()
        self.sim_time = 0.0
        self.last_ignite_cycle = -1
        self.burned = False
        self.hist_time = []
        self.hist_momentum = []
        self.hist_thrust = []
        self.mom_out_left = 0.0
        self.mom_out_right = 0.0
        self.prev_accel = None

    def _init_grid(self):
        self.nodes_x, _ = generate_equal_mass_grid(self.get_area, self.length, target_mass=self.target_mass)
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
        if not np.all(np.isfinite(self.nodes_x)) or not np.all(np.isfinite(self.nodes_v)):
            raise RuntimeError(f"NaN/Inf in solver state at t={self.sim_time:.6f}s — simulation diverged")

    def _step_engine(self):
        self._dynamic_topology_engine()
        dx = np.diff(self.nodes_x)
        min_dx = 1e-6
        if np.any(dx < min_dx):
            for i in range(len(dx)):
                if self.nodes_x[i + 1] <= self.nodes_x[i] + min_dx:
                    self.nodes_x[i + 1] = self.nodes_x[i] + min_dx + 1e-9
            dx = np.diff(self.nodes_x)
        centers = self.get_centers()
        A_L = self.get_area(self.nodes_x[:-1])
        A_R = self.get_area(self.nodes_x[1:])
        vols = conic_frustum_vol(dx, A_L, A_R)
        vols = np.maximum(vols, 1e-12)
        elem_masses = np.sum(self.elem_species, axis=1)
        Y_fracs = self.elem_species / elem_masses[:, None]
        rho = elem_masses / vols
        R_mix = np.sum(Y_fracs * R_SPECIES, axis=1)
        pressures = self.compute_pressure(rho, self.elem_entropy)
        temps = pressures / (rho * R_mix)
        c_sound = np.sqrt(GAMMA * pressures / rho)
        dt = cfl_dt(dx, c_sound, self.nodes_v[1:], self.cfl, self.dt_max)
        self.sim_time += dt
        node_m = np.zeros(len(self.nodes_x))
        node_m[1:-1] = 0.5 * (elem_masses[:-1] + elem_masses[1:])
        node_m[0] = elem_masses[0]
        node_m[-1] = elem_masses[-1]
        # --- A4: Verlet first half-kick (disabled for engine mode — unstable with dynamic topology) ---
        if self.prev_accel is not None and len(self.prev_accel) == len(self.nodes_v) and self.mode != "ENGINE":
            self.nodes_v += 0.5 * self.prev_accel * dt
        # --- Drift ---
        self.nodes_x += self.nodes_v * dt
        dx = np.diff(self.nodes_x)
        if np.any(dx < min_dx):
            for i in range(len(dx)):
                if self.nodes_x[i + 1] <= self.nodes_x[i] + min_dx:
                    self.nodes_x[i + 1] = self.nodes_x[i] + min_dx + 1e-9
            dx = np.diff(self.nodes_x)
        # --- Recompute thermodynamics at new positions ---
        centers = self.get_centers()
        A_L = self.get_area(self.nodes_x[:-1])
        A_R = self.get_area(self.nodes_x[1:])
        vols = conic_frustum_vol(dx, A_L, A_R)
        vols = np.maximum(vols, 1e-12)
        elem_masses = np.sum(self.elem_species, axis=1)
        Y_fracs = self.elem_species / elem_masses[:, None]
        rho = elem_masses / vols
        R_mix = np.sum(Y_fracs * R_SPECIES, axis=1)
        pressures = self.compute_pressure(rho, self.elem_entropy)
        temps = pressures / (rho * R_mix)
        # --- Combustion ---
        self._entropy_before_combustion = self.elem_entropy.copy()
        if self.auto_inject and self.sim_time >= self.relax_time and self.sim_time >= self.inject_start:
            phase = (self.sim_time - self.inject_start) % self.inject_period
            cycle_idx = int((self.sim_time - self.inject_start) // self.inject_period)
            if 0 <= phase < self.inject_duration:
                inj_x = self.inject_point - self.geom_x0
                chamber_mask = np.abs(centers - inj_x) < 0.15
                chamber_idx = np.where(chamber_mask)[0]
                if len(chamber_idx) > 0:
                    fuel_per_elem = (self.fuel_mass / self.inject_duration) * dt / len(chamber_idx)
                    for idx in chamber_idx:
                        self.elem_species[idx, 2] += fuel_per_elem
                else:
                    idx = np.argmin(np.abs(centers - inj_x))
                    self.elem_species[idx, 2] += (self.fuel_mass / self.inject_duration) * dt
            ign_start = self.inject_duration + self.ignite_delay
            burn_window = 0.001  # burn over 1ms window
            if ign_start <= phase < ign_start + burn_window:
                self._do_ignite(centers, temps, burn_fraction=0.3)
        # --- Manual injection (one-shot, spread over chamber) ---
        if self.manual_inject_x is not None:
            inj_x = self.manual_inject_x - self.geom_x0
            chamber_mask = np.abs(centers - inj_x) < 0.15
            chamber_idx = np.where(chamber_mask)[0]
            if len(chamber_idx) > 0:
                fuel_per_elem = self.fuel_mass / len(chamber_idx)
                for idx in chamber_idx:
                    self.elem_species[idx, 2] += fuel_per_elem
            else:
                idx = np.argmin(np.abs(centers - inj_x))
                self.elem_species[idx, 2] += self.fuel_mass
            self.manual_inject_x = None
        # --- Reset Verlet prev_accel if entropy changed significantly (combustion discontinuity) ---
        if hasattr(self, '_entropy_before_combustion'):
            dS_max = np.max(np.abs(self.elem_entropy - self._entropy_before_combustion))
            if dS_max > 0.5:
                self.prev_accel = None
            del self._entropy_before_combustion
        # --- Heat conduction ---
        dT_bound = np.concatenate(([T_ATM], temps, [T_ATM]))
        dx_nodes = np.concatenate(([dx[0]], (dx[:-1] + dx[1:]) * 0.5, [dx[-1]]))
        dx_nodes = np.maximum(dx_nodes, 1e-12)
        dT = dT_bound[1:] - dT_bound[:-1]
        heat_flux = -self.thermal_cond * self.get_area(self.nodes_x) * dT / dx_nodes
        net_heat = (heat_flux[:-1] - heat_flux[1:]) * dt
        self.elem_entropy += net_heat / (elem_masses * temps)
        # --- Species diffusion ---
        for sp in range(5):
            Y = self.elem_species[:, sp] / elem_masses
            Y_bound = np.concatenate(([Y_AIR[sp]], Y, [Y_AIR[sp]]))
            dY = Y_bound[1:] - Y_bound[:-1]
            mass_flux = -self.species_diff * self.get_area(self.nodes_x) * dY / dx_nodes * 10.0
            net_mass_flow = (mass_flux[:-1] - mass_flux[1:]) * dt
            self.elem_species[:, sp] += net_mass_flow
        # --- Recompute thermodynamics after source terms ---
        elem_masses = np.sum(self.elem_species, axis=1)
        Y_fracs = self.elem_species / elem_masses[:, None]
        rho = elem_masses / vols
        R_mix = np.sum(Y_fracs * R_SPECIES, axis=1)
        pressures = self.compute_pressure(rho, self.elem_entropy)
        temps = pressures / (rho * R_mix)
        # --- Artificial viscosity ---
        du = self.nodes_v[1:] - self.nodes_v[:-1]
        q = self.artificial_viscosity(du, rho, self.viscosity)
        P_total = pressures + q
        # --- L1: Forces (pressure gradient + wall force) ---
        force = self.compute_forces(self.nodes_x, P_total, P_ATM, P_ATM)
        force += self.compute_wall_forces(self.nodes_x, P_total, P_ATM, P_ATM)
        node_m = np.zeros(len(self.nodes_x))
        node_m[1:-1] = 0.5 * (elem_masses[:-1] + elem_masses[1:])
        node_m[0] = elem_masses[0]
        node_m[-1] = elem_masses[-1]
        force += -self.friction * self.nodes_v * np.abs(self.nodes_v) * node_m
        # --- L3: Jet damping (valveless mechanism) ---
        if self.jet_damping > 0:
            in_left_res = (self.nodes_x < self.engine_start_local) & (self.nodes_v < 0)
            in_right_res = (self.nodes_x > self.engine_end_local) & (self.nodes_v > 0)
            force[in_left_res] -= self.jet_damping * self.nodes_v[in_left_res] * node_m[in_left_res]
            force[in_right_res] -= self.jet_damping * self.nodes_v[in_right_res] * node_m[in_right_res]
        accel = force / node_m
        # --- A4: Verlet second half-kick (full kick for engine mode = semi-implicit Euler) ---
        kick = 1.0 if self.mode == "ENGINE" else 0.5
        self.nodes_v += kick * accel * dt
        self.prev_accel = accel
        # --- L4: Diagnostics ---
        total_mom = np.sum(node_m * self.nodes_v)
        self.hist_time.append(self.sim_time)
        self.hist_momentum.append(total_mom)
        self.hist_thrust.append(self._compute_thrust(P_total, centers, rho))

    def _step_inertia(self):
        A_in = self.get_area(self.nodes_x[0])
        spawn_threshold = self.target_mass / (A_in * RHO_ATM)
        if self.nodes_x[0] > spawn_threshold:
            x_new = self.nodes_x[0] - spawn_threshold
            self.nodes_x = np.insert(self.nodes_x, 0, x_new)
            self.nodes_v = np.insert(self.nodes_v, 0, self.nodes_v[0])
            self.elem_entropy = np.insert(self.elem_entropy, 0, 1.0)
        elif self.nodes_x[0] < -0.2:
            self.nodes_x = self.nodes_x[1:]
            self.nodes_v = self.nodes_v[1:]
            self.elem_entropy = self.elem_entropy[1:]
        A_out = self.get_area(self.nodes_x[-1])
        spawn_threshold_r = self.target_mass / (A_out * RHO_ATM)
        if self.nodes_x[-1] < self.length - spawn_threshold_r:
            x_new = self.nodes_x[-1] + spawn_threshold_r
            self.nodes_x = np.append(self.nodes_x, x_new)
            self.nodes_v = np.append(self.nodes_v, self.nodes_v[-1])
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
        centers = self.get_centers()
        vols, _ = self.get_volumes()
        n_cells = len(dx)
        mass = np.full(n_cells, self.target_mass)
        rho = mass / vols
        P_thermo = self.compute_pressure(rho, self.elem_entropy)
        du = self.nodes_v[1:] - self.nodes_v[:-1]
        q = self.artificial_viscosity(du, rho, self.viscosity)
        P_tot = P_thermo + q
        if self.mode == "INERTIA":
            P_ghost_L = P_tot[0]; P_ghost_R = P_tot[-1]
        else:
            P_ghost_L = P_ATM; P_ghost_R = P_ATM
        force = self.compute_forces(self.nodes_x, P_tot, P_ghost_L, P_ghost_R)
        force += self.compute_wall_forces(self.nodes_x, P_tot, P_ghost_L, P_ghost_R)
        node_m = node_mass_from_cell_mass(mass)
        if self.friction > 0:
            force -= self.friction * self.nodes_v * np.abs(self.nodes_v) * node_m
        accel = force / node_m
        c_s = np.sqrt(GAMMA * P_tot / rho)
        dt = cfl_dt(dx, c_s, self.nodes_v[1:], self.cfl, self.dt_max)
        if self.mode == "ENGINE":
            self.sim_time += dt
            if not self.burned and self.sim_time > self.burn_time:
                self.burned = True
                mask = (centers > 0.35) & (centers < 0.55)
                self.elem_entropy[mask] *= self.burn_mult
        # --- A4: Verlet ---
        if self.prev_accel is not None and len(self.prev_accel) == len(self.nodes_v):
            self.nodes_v += 0.5 * self.prev_accel * dt
        self.nodes_x += self.nodes_v * dt
        # Recompute forces at new positions (simplified for inertia mode)
        dx = np.diff(self.nodes_x)
        if np.any(dx < min_dx):
            dx = np.maximum(dx, min_dx)
            self.nodes_x[1:] = self.nodes_x[0] + np.cumsum(dx)
        vols, _ = self.get_volumes()
        rho = mass / vols
        P_thermo = self.compute_pressure(rho, self.elem_entropy)
        du = self.nodes_v[1:] - self.nodes_v[:-1]
        q = self.artificial_viscosity(du, rho, self.viscosity)
        P_tot = P_thermo + q
        if self.mode == "INERTIA":
            P_ghost_L = P_tot[0]; P_ghost_R = P_tot[-1]
        else:
            P_ghost_L = P_ATM; P_ghost_R = P_ATM
        force = self.compute_forces(self.nodes_x, P_tot, P_ghost_L, P_ghost_R)
        force += self.compute_wall_forces(self.nodes_x, P_tot, P_ghost_L, P_ghost_R)
        node_m = node_mass_from_cell_mass(mass)
        if self.friction > 0:
            force -= self.friction * self.nodes_v * np.abs(self.nodes_v) * node_m
        accel = force / node_m
        self.nodes_v += 0.5 * accel * dt
        self.prev_accel = accel
        total_mom = np.sum(node_m * self.nodes_v)
        self.hist_time.append(self.sim_time)
        self.hist_momentum.append(total_mom)
        self.hist_thrust.append(self._compute_thrust(P_tot, centers, rho))

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
            # L4: Track momentum of culled node (thrust measurement)
            node_m_cull = 0.5 * np.sum(self.elem_species[0])
            self.mom_out_left += node_m_cull * self.nodes_v[0]
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
            # L4: Track momentum of culled node (thrust measurement)
            node_m_cull = 0.5 * np.sum(self.elem_species[-1])
            self.mom_out_right += node_m_cull * self.nodes_v[-1]
            self.nodes_x = self.nodes_x[:-1]
            self.nodes_v = self.nodes_v[:-1]
            self.elem_species = self.elem_species[:-1]
            self.elem_entropy = self.elem_entropy[:-1]

    def _compute_thrust(self, P_total, centers, rho):
        """Instantaneous thrust at engine boundaries: momentum flux + pressure thrust."""
        exit_elems = np.where(centers >= self.engine_end_local)[0]
        intake_elems = np.where(centers <= self.engine_start_local)[0]
        if len(exit_elems) == 0 or len(intake_elems) == 0:
            return 0.0
        i_exit = exit_elems[0]
        i_intake = intake_elems[-1]
        p_exit = P_total[i_exit]
        p_intake = P_total[i_intake]
        A_exit = self.get_area(self.engine_end_local)
        A_intake = self.get_area(self.engine_start_local)
        v_exit = 0.5 * (self.nodes_v[i_exit] + self.nodes_v[i_exit + 1])
        v_intake = 0.5 * (self.nodes_v[i_intake] + self.nodes_v[i_intake + 1])
        rho_exit = rho[i_exit]
        rho_intake = rho[i_intake]
        mdot_exit = rho_exit * v_exit * A_exit
        mdot_intake = rho_intake * v_intake * A_intake
        thrust = (mdot_exit * v_exit - mdot_intake * v_intake
                  + (p_exit - P_ATM) * A_exit - (p_intake - P_ATM) * A_intake)
        return thrust

    def _do_ignite(self, centers, temps, burn_fraction=0.3):
        inj_x = self.inject_point - self.geom_x0
        chamber_mask = np.abs(centers - inj_x) < 0.15
        T_floor = 200.0
        max_dS = 5.0  # cap entropy increase per element per step
        for i in np.where(chamber_mask)[0]:
            burn = min(self.elem_species[i, 2], self.elem_species[i, 1] / O2_PER_FUEL) * burn_fraction
            if burn > 1e-12:
                self.elem_species[i, 2] -= burn
                self.elem_species[i, 1] -= burn * O2_PER_FUEL
                self.elem_species[i, 3] += burn * Y_PROD_CO2
                self.elem_species[i, 4] += burn * Y_PROD_H2O
                T_safe = max(temps[i], T_floor)
                dS = (burn * HEAT_OF_COMBUSTION / T_safe) * 0.8
                self.elem_entropy[i] += min(dS, max_dS)

    def manual_inject(self, x=None):
        """Queue a one-shot fuel injection at position x (absolute coords). Defaults to inject_point."""
        self.manual_inject_x = x if x is not None else self.inject_point

    def manual_ignite(self):
        """Immediately burn fuel in chamber region (iterative to consume most fuel)."""
        centers = self.get_centers()
        dx = np.diff(self.nodes_x)
        A_L = self.get_area(self.nodes_x[:-1])
        A_R = self.get_area(self.nodes_x[1:])
        vols = conic_frustum_vol(dx, A_L, A_R)
        vols = np.maximum(vols, 1e-12)
        elem_masses = np.sum(self.elem_species, axis=1)
        Y_fracs = self.elem_species / elem_masses[:, None]
        rho = elem_masses / vols
        R_mix = np.sum(Y_fracs * R_SPECIES, axis=1)
        pressures = self.compute_pressure(rho, self.elem_entropy)
        temps = pressures / (rho * R_mix)
        for _ in range(10):
            self._do_ignite(centers, temps, burn_fraction=0.1)
            elem_masses = np.sum(self.elem_species, axis=1)
            Y_fracs = self.elem_species / elem_masses[:, None]
            rho = elem_masses / vols
            R_mix = np.sum(Y_fracs * R_SPECIES, axis=1)
            pressures = self.compute_pressure(rho, self.elem_entropy)
            temps = pressures / (rho * R_mix)

    def get_avg_thrust(self):
        """Average thrust from cumulative momentum of culled nodes."""
        if self.sim_time > 0:
            return (self.mom_out_right - self.mom_out_left) / self.sim_time
        return 0.0

    def get_state(self):
        dx = np.diff(self.nodes_x)
        centers = 0.5 * (self.nodes_x[:-1] + self.nodes_x[1:])
        A_L = self.get_area(self.nodes_x[:-1])
        A_R = self.get_area(self.nodes_x[1:])
        vols = conic_frustum_vol(dx, A_L, A_R)
        vols = np.maximum(vols, 1e-12)
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
        p = self.compute_pressure(rho, self.elem_entropy)
        t = p / (rho * R_mix)
        return centers, rho, p, t, Y, vols, m_t
