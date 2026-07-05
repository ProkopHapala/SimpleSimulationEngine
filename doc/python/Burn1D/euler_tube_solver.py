import numpy as np
from tube_common import GAMMA, R_GAS, P_ATM, RHO_ATM, T_ATM

CFL = 0.5

# Combustion constants (shared with Lagrangian)
HEAT_OF_COMBUSTION = 45e6
O2_PER_FUEL = 3.58
Y_PROD_CO2 = 3.03
Y_PROD_H2O = 1.55
T_IGNITE = 600.0


def minmod(a, b):
    """Minmod slope limiter for MUSCL reconstruction."""
    return np.where(a * b > 0, np.sign(a) * np.minimum(np.abs(a), np.abs(b)), 0.0)

class CompressibleSolver:
    def __init__(self, nx, length, geom_segments, u0, u_amb, p_init, rho_init, bc_mode="ambient",
                 combustion=False, jet_damping=0.0, engine_start=0.0, engine_end=None,
                 inject_point=None, inject_period=0.020, inject_duration=0.005,
                 fuel_mass=3e-4, ignite_delay=0.001, muscl=True, flux_type="hllc"):
        self.nx = nx
        self.dx = length / nx
        self.t = 0.0
        self.u_amb = u_amb
        self.bc_mode = bc_mode
        self.combustion = combustion
        self.jet_damping = jet_damping
        self.engine_start = engine_start
        self.engine_end = engine_end if engine_end is not None else length
        self.inject_point = inject_point if inject_point is not None else 0.5 * length
        self.inject_period = inject_period
        self.inject_duration = inject_duration
        self.fuel_mass = fuel_mass
        self.ignite_delay = ignite_delay
        self.muscl = muscl
        self.flux_type = flux_type
        self.x = np.linspace(0.5 * self.dx, length - 0.5 * self.dx, nx)
        self.A, self.A_faces, self.dA_dx_term, self.vol = self._generate_geometry(length, geom_segments)
        n_vars = 4 if combustion else 3
        self.U = np.zeros((n_vars, nx))
        e_init = p_init / ((GAMMA - 1) * rho_init)
        E_init = rho_init * (e_init + 0.5 * u0**2)
        self.U[0, :] = rho_init * self.A
        self.U[1, :] = rho_init * u0 * self.A
        self.U[2, :] = E_init * self.A
        if combustion:
            self.U[3, :] = 0.0
        self.p_amb = p_init
        self.rho_amb = rho_init
        self.T_amb = p_init / (R_GAS * rho_init)
        self.hist_t = []
        self.hist_thrust = []
        self.last_inject_cycle = -1

    @staticmethod
    def _fill_slopes(xs, ys, ms):
        ms = np.array(ms, dtype=float)
        for i in range(len(ms)):
            if np.isnan(ms[i]):
                if i == 0: ms[i] = (ys[1] - ys[0]) / (xs[1] - xs[0])
                elif i == len(ms) - 1: ms[i] = (ys[-1] - ys[-2]) / (xs[-1] - xs[-2])
                else: ms[i] = (ys[i + 1] - ys[i - 1]) / (xs[i + 1] - xs[i - 1])
        return ms

    @staticmethod
    def _hermite_eval(xq, xs, ys, ms):
        xs = np.asarray(xs); ys = np.asarray(ys); ms = np.asarray(ms)
        idx = np.searchsorted(xs, xq, side='right') - 1
        idx = np.clip(idx, 0, len(xs) - 2)
        x0 = xs[idx]; x1 = xs[idx + 1]
        y0 = ys[idx]; y1 = ys[idx + 1]
        m0 = ms[idx]; m1 = ms[idx + 1]
        dx = x1 - x0
        t = (xq - x0) / dx
        t = np.clip(t, 0.0, 1.0)
        h00 = 2 * t**3 - 3 * t**2 + 1
        h10 = t**3 - 2 * t**2 + t
        h01 = -2 * t**3 + 3 * t**2
        h11 = t**3 - t**2
        return h00 * y0 + h10 * dx * m0 + h01 * y1 + h11 * dx * m1

    def _generate_geometry(self, length, segments):
        xs = [s[0] for s in segments]
        ds = [s[1] for s in segments]
        ms = [s[2] if len(s) > 2 else np.nan for s in segments]
        ms = self._fill_slopes(xs, ds, ms)
        d_center = self._hermite_eval(self.x, xs, ds, ms)
        x_edges = np.linspace(0.0, length, self.nx + 1)
        d_edge = self._hermite_eval(x_edges, xs, ds, ms)
        A_center = np.pi * 0.25 * d_center**2
        A_faces = np.pi * 0.25 * d_edge**2
        dA_dx_term = (A_faces[1:] - A_faces[:-1]) / self.dx
        return A_center, A_faces, dA_dx_term, A_center * self.dx

    def get_primitive(self):
        rho = self.U[0] / self.A
        u = self.U[1] / self.U[0]
        E_vol = self.U[2] / self.A
        p = (GAMMA - 1) * (E_vol - 0.5 * rho * u**2)
        T = p / (R_GAS * rho)
        Y_fuel = self.U[3] / self.U[0] if self.combustion else None
        return rho, u, p, T, Y_fuel

    def compute_flux_rusanov_well_balanced(self, U_L, U_R, A_L, A_R, A_face):
        rho_L = U_L[0] / A_L
        u_L = U_L[1] / U_L[0]
        p_L = (GAMMA - 1) * ((U_L[2] / A_L) - 0.5 * rho_L * u_L**2)
        E_vol_L = U_L[2] / A_L
        rho_R = U_R[0] / A_R
        u_R = U_R[1] / U_R[0]
        p_R = (GAMMA - 1) * ((U_R[2] / A_R) - 0.5 * rho_R * u_R**2)
        E_vol_R = U_R[2] / A_R
        n = len(U_L)
        F_L = np.zeros(n)
        F_R = np.zeros(n)
        F_L[0] = rho_L * u_L; F_L[1] = rho_L * u_L**2 + p_L; F_L[2] = u_L * (E_vol_L + p_L)
        F_R[0] = rho_R * u_R; F_R[1] = rho_R * u_R**2 + p_R; F_R[2] = u_R * (E_vol_R + p_R)
        if n > 3:
            F_L[3] = rho_L * u_L * (U_L[3] / U_L[0])
            F_R[3] = rho_R * u_R * (U_R[3] / U_R[0])
        F_L *= A_face; F_R *= A_face
        c_L = np.sqrt(GAMMA * abs(p_L) / (abs(rho_L) + 1e-12))
        c_R = np.sqrt(GAMMA * abs(p_R) / (abs(rho_R) + 1e-12))
        S_max = max(abs(u_L) + c_L, abs(u_R) + c_R)
        U_star_L = U_L * (A_face / A_L)
        U_star_R = U_R * (A_face / A_R)
        F_face = 0.5 * (F_L + F_R) - 0.5 * S_max * (U_star_R - U_star_L)
        return F_face, S_max

    def compute_flux_hllc(self, U_L, U_R, A_L, A_R, A_face):
        """HLLC flux — restores contact wave for sharper species interfaces."""
        rho_L = U_L[0] / A_L
        u_L = U_L[1] / U_L[0]
        p_L = (GAMMA - 1) * ((U_L[2] / A_L) - 0.5 * rho_L * u_L**2)
        E_vol_L = U_L[2] / A_L
        rho_R = U_R[0] / A_R
        u_R = U_R[1] / U_R[0]
        p_R = (GAMMA - 1) * ((U_R[2] / A_R) - 0.5 * rho_R * u_R**2)
        E_vol_R = U_R[2] / A_R
        c_L = np.sqrt(GAMMA * abs(p_L) / (abs(rho_L) + 1e-12))
        c_R = np.sqrt(GAMMA * abs(p_R) / (abs(rho_R) + 1e-12))
        S_L = min(u_L - c_L, u_R - c_R)
        S_R = max(u_L + c_L, u_R + c_R)
        S_M = (p_R - p_L + rho_L * u_L * (S_L - u_L) - rho_R * u_R * (S_R - u_R)) / (
               rho_L * (S_L - u_L) - rho_R * (S_R - u_R) + 1e-12)
        p_star = p_L + rho_L * (S_L - u_L) * (S_M - u_L)
        n = len(U_L)
        F_L = np.zeros(n)
        F_R = np.zeros(n)
        F_L[0] = rho_L * u_L; F_L[1] = rho_L * u_L**2 + p_L; F_L[2] = u_L * (E_vol_L + p_L)
        F_R[0] = rho_R * u_R; F_R[1] = rho_R * u_R**2 + p_R; F_R[2] = u_R * (E_vol_R + p_R)
        if n > 3:
            F_L[3] = rho_L * u_L * (U_L[3] / U_L[0])
            F_R[3] = rho_R * u_R * (U_R[3] / U_R[0])
        F_L *= A_face; F_R *= A_face
        S_max = max(abs(S_L), abs(S_R))
        if S_L >= 0:
            return F_L, S_max
        if S_R <= 0:
            return F_R, S_max
        # Compute star states
        rho_star_L = rho_L * (S_L - u_L) / (S_L - S_M)
        rho_star_R = rho_R * (S_R - u_R) / (S_R - S_M)
        U_star_L = np.zeros(n)
        U_star_R = np.zeros(n)
        U_star_L[0] = rho_star_L * A_face
        U_star_L[1] = rho_star_L * S_M * A_face
        e_star_L = p_star / (GAMMA - 1) + 0.5 * rho_star_L * S_M**2
        U_star_L[2] = e_star_L * A_face
        U_star_R[0] = rho_star_R * A_face
        U_star_R[1] = rho_star_R * S_M * A_face
        e_star_R = p_star / (GAMMA - 1) + 0.5 * rho_star_R * S_M**2
        U_star_R[2] = e_star_R * A_face
        if n > 3:
            Y_fuel_L = U_L[3] / U_L[0]
            Y_fuel_R = U_R[3] / U_R[0]
            U_star_L[3] = rho_star_L * Y_fuel_L * A_face
            U_star_R[3] = rho_star_R * Y_fuel_R * A_face
        if S_M >= 0:
            F_star = F_L + S_L * (U_star_L - U_L * (A_face / A_L))
            return F_star, S_max
        else:
            F_star = F_R + S_R * (U_star_R - U_R * (A_face / A_R))
            return F_star, S_max

    def _apply_bc(self, U_padded, A_padded, rho, u, e_amb):
        """Fill ghost cells based on bc_mode."""
        nv = U_padded.shape[0]
        if self.bc_mode == "nonreflecting":
            c_amb = np.sqrt(GAMMA * self.p_amb / self.rho_amb)
            # Left boundary: J- outgoing (from interior), J+ incoming (from ambient)
            c0 = np.sqrt(GAMMA * abs(self.p_amb) / (abs(rho[0]) + 1e-12))
            c0 = max(c0, 1e-6)
            J_minus = u[0] - 2*c0/(GAMMA-1)          # outgoing, from interior
            J_plus_amb = self.u_amb + 2*c_amb/(GAMMA-1)  # incoming, from ambient
            u_gL = 0.5*(J_plus_amb + J_minus)
            c_gL = max(0.25*(GAMMA-1)*(J_plus_amb - J_minus), 1e-6)
            p_gL = self.p_amb * (c_gL/c_amb)**(2*GAMMA/(GAMMA-1))
            rho_gL = GAMMA * p_gL / c_gL**2
            e_gL = p_gL / ((GAMMA-1)*rho_gL)
            U_padded[0, 0] = rho_gL * A_padded[0]
            U_padded[1, 0] = rho_gL * u_gL * A_padded[0]
            U_padded[2, 0] = rho_gL * (e_gL + 0.5*u_gL**2) * A_padded[0]
            if nv > 3: U_padded[3, 0] = 0.0
            # Right boundary: J+ outgoing (from interior), J- incoming (from ambient)
            cN = np.sqrt(GAMMA * abs(self.p_amb) / (abs(rho[-1]) + 1e-12))
            cN = max(cN, 1e-6)
            J_plus = u[-1] + 2*cN/(GAMMA-1)           # outgoing, from interior
            J_minus_amb = -self.u_amb - 2*c_amb/(GAMMA-1)  # incoming, from ambient
            u_gR = 0.5*(J_plus + J_minus_amb)
            c_gR = max(0.25*(GAMMA-1)*(J_plus - J_minus_amb), 1e-6)
            p_gR = self.p_amb * (c_gR/c_amb)**(2*GAMMA/(GAMMA-1))
            rho_gR = GAMMA * p_gR / c_gR**2
            e_gR = p_gR / ((GAMMA-1)*rho_gR)
            U_padded[0, -1] = rho_gR * A_padded[-1]
            U_padded[1, -1] = rho_gR * u_gR * A_padded[-1]
            U_padded[2, -1] = rho_gR * (e_gR + 0.5*u_gR**2) * A_padded[-1]
            if nv > 3: U_padded[3, -1] = 0.0
        elif self.bc_mode == "ambient":
            U_padded[0, 0] = self.rho_amb * A_padded[0]
            U_padded[1, 0] = self.rho_amb * self.u_amb * A_padded[0]
            U_padded[2, 0] = self.rho_amb * (e_amb + 0.5 * self.u_amb**2) * A_padded[0]
            if nv > 3: U_padded[3, 0] = 0.0
            u_in = -self.u_amb
            U_padded[0, -1] = self.rho_amb * A_padded[-1]
            U_padded[1, -1] = self.rho_amb * u_in * A_padded[-1]
            U_padded[2, -1] = self.rho_amb * (e_amb + 0.5 * u_in**2) * A_padded[-1]
            if nv > 3: U_padded[3, -1] = 0.0
        else:
            if u[0] > 0:
                U_padded[0, 0] = self.rho_amb * A_padded[0]
                U_padded[1, 0] = self.rho_amb * self.u_amb * A_padded[0]
                U_padded[2, 0] = self.rho_amb * (e_amb + 0.5 * self.u_amb**2) * A_padded[0]
                if nv > 3: U_padded[3, 0] = 0.0
            else:
                U_padded[:, 0] = self.U[:, 0] * (A_padded[0] / self.A[0])
            if u[-1] < 0:
                u_in = -self.u_amb
                U_padded[0, -1] = self.rho_amb * A_padded[-1]
                U_padded[1, -1] = self.rho_amb * u_in * A_padded[-1]
                U_padded[2, -1] = self.rho_amb * (e_amb + 0.5 * u_in**2) * A_padded[-1]
                if nv > 3: U_padded[3, -1] = 0.0
            else:
                rho_g = rho[-1]; u_g = u[-1]; p_g = self.p_amb
                e_g = p_g / ((GAMMA - 1) * rho_g)
                U_padded[0, -1] = rho_g * A_padded[-1]
                U_padded[1, -1] = rho_g * u_g * A_padded[-1]
                U_padded[2, -1] = rho_g * (e_g + 0.5 * u_g**2) * A_padded[-1]
                if nv > 3: U_padded[3, -1] = self.U[3, -1] * (A_padded[-1] / self.A[-1])

    def _U_to_W(self, U, A):
        """Convert conserved variables to primitive (rho, u, p, [Y_fuel])."""
        nv = U.shape[0]
        rho = U[0] / A
        u = U[1] / (U[0] + 1e-12)
        E_vol = U[2] / A
        p = (GAMMA - 1) * (E_vol - 0.5 * rho * u**2)
        if nv > 3:
            Y = U[3] / (U[0] + 1e-12)
            return np.array([rho, u, p, Y])
        return np.array([rho, u, p])

    def _W_to_U(self, W, A):
        """Convert primitive variables to conserved."""
        nv = W.shape[0]
        rho = W[0]; u = W[1]; p = W[2]
        E_vol = p / (GAMMA - 1) + 0.5 * rho * u**2
        U = np.zeros(nv)
        U[0] = rho * A
        U[1] = rho * u * A
        U[2] = E_vol * A
        if nv > 3:
            U[3] = W[3] * rho * A
        return U

    def _compute_fluxes(self, U_padded, A_padded):
        """Compute fluxes at all interfaces, with optional MUSCL reconstruction."""
        nv = U_padded.shape[0]
        Flux_interface = np.zeros((nv, self.nx + 1))
        flux_fn = self.compute_flux_hllc if self.flux_type == "hllc" else self.compute_flux_rusanov_well_balanced
        if self.muscl:
            # MUSCL on primitive variables (smooth across area changes)
            W_padded = np.zeros((nv if nv > 3 else 3, self.nx + 2))
            for j in range(self.nx + 2):
                W_padded[:, j] = self._U_to_W(U_padded[:, j], A_padded[j])
            slopes = np.zeros_like(W_padded)
            for v in range(W_padded.shape[0]):
                dL = W_padded[v, 1:-1] - W_padded[v, :-2]
                dR = W_padded[v, 2:]  - W_padded[v, 1:-1]
                slopes[v, 1:-1] = minmod(dL, dR)
            for i in range(self.nx + 1):
                W_L = W_padded[:, i] + 0.5 * slopes[:, i]
                W_R = W_padded[:, i + 1] - 0.5 * slopes[:, i + 1]
                A_L = A_padded[i]
                A_R = A_padded[i + 1]
                A_face = self.A_faces[i]
                U_L = self._W_to_U(W_L, A_L)
                U_R = self._W_to_U(W_R, A_R)
                Flux_interface[:, i], _ = flux_fn(U_L, U_R, A_L, A_R, A_face)
        else:
            for i in range(self.nx + 1):
                U_L = U_padded[:, i]
                U_R = U_padded[:, i + 1]
                A_L = A_padded[i]
                A_R = A_padded[i + 1]
                A_face = self.A_faces[i]
                Flux_interface[:, i], _ = flux_fn(U_L, U_R, A_L, A_R, A_face)
        return Flux_interface

    def _combustion_source(self, dt):
        """E2: Combustion source — fuel injection + heat release (cycle-based ignition)."""
        rho, u, p, T, Y_fuel = self.get_primitive()
        S = np.zeros_like(self.U)
        t_old = self.t - dt
        # Fuel injection
        phase = self.t % self.inject_period
        phase_old = t_old % self.inject_period
        cycle_idx = int(self.t // self.inject_period)
        if phase < self.inject_duration:
            chamber_mask = np.abs(self.x - self.inject_point) < 0.1
            chamber_idx = np.where(chamber_mask)[0]
            inject_rate = self.fuel_mass / self.inject_duration
            if len(chamber_idx) > 0:
                for idx in chamber_idx:
                    S[3, idx] += inject_rate / (self.dx * len(chamber_idx))
                    S[0, idx] += inject_rate / (self.dx * len(chamber_idx))
            else:
                idx = np.argmin(np.abs(self.x - self.inject_point))
                S[3, idx] += inject_rate / self.dx
                S[0, idx] += inject_rate / self.dx
        # Ignition: detect phase crossing ign_start boundary
        ign_start = self.inject_duration + self.ignite_delay
        crossed = (phase_old < ign_start <= phase) or (cycle_idx != int(t_old // self.inject_period) and ign_start <= phase + self.inject_period)
        ignite_now = crossed and (cycle_idx != self.last_inject_cycle)
        if ignite_now:
            self.last_inject_cycle = cycle_idx
        if Y_fuel is not None:
            if ignite_now:
                chamber_mask = np.abs(self.x - self.inject_point) < 0.1
                burn_mask = chamber_mask & (Y_fuel > 1e-10)
                if np.any(burn_mask):
                    burn_rate = Y_fuel[burn_mask] * rho[burn_mask] / dt
                    S[3, burn_mask] -= burn_rate * self.A[burn_mask]
                    S[2, burn_mask] += burn_rate * HEAT_OF_COMBUSTION * self.A[burn_mask]
            else:
                burn_mask = (T > T_IGNITE) & (Y_fuel > 1e-10)
                if np.any(burn_mask):
                    burn_rate = np.minimum(Y_fuel[burn_mask] * rho[burn_mask] / dt, rho[burn_mask] * 0.1)
                    S[3, burn_mask] -= burn_rate * self.A[burn_mask]
                    S[2, burn_mask] += burn_rate * HEAT_OF_COMBUSTION * self.A[burn_mask]
        return S

    def _jet_damping_source(self, u):
        """E3: Jet damping — momentum/energy sink in reservoir zones for outflow."""
        if self.jet_damping <= 0:
            return np.zeros_like(self.U)
        S = np.zeros_like(self.U)
        left_res = (self.x < self.engine_start) & (u < 0)
        right_res = (self.x > self.engine_end) & (u > 0)
        S[1, left_res]  -= self.jet_damping * self.U[1, left_res]
        S[1, right_res] -= self.jet_damping * self.U[1, right_res]
        S[2, left_res]  -= self.jet_damping * 0.5 * self.U[1, left_res]**2 / (self.U[0, left_res] + 1e-12)
        S[2, right_res] -= self.jet_damping * 0.5 * self.U[1, right_res]**2 / (self.U[0, right_res] + 1e-12)
        return S

    def _compute_thrust(self, rho, u, p):
        """E4: Thrust at engine exit: momentum flux + pressure thrust."""
        i_exit = np.argmin(np.abs(self.x - self.engine_end))
        i_intake = np.argmin(np.abs(self.x - self.engine_start))
        A_exit = self.A[i_exit]
        A_intake = self.A[i_intake]
        thrust = (rho[i_exit] * u[i_exit]**2 * A_exit
                  - rho[i_intake] * u[i_intake]**2 * A_intake
                  + (p[i_exit] - self.p_amb) * A_exit
                  - (p[i_intake] - self.p_amb) * A_intake)
        return thrust

    def step(self, dt_fixed=None, prop_params=None):
        rho, u, p, T, Y_fuel = self.get_primitive()
        c = np.sqrt(GAMMA * abs(p) / (abs(rho) + 1e-12))
        max_speed = np.max(np.abs(u) + c)
        if max_speed < 1e-5: max_speed = 1e-5
        dt = CFL * self.dx / max_speed
        if dt_fixed: dt = dt_fixed
        self.t += dt
        nv = self.U.shape[0]
        A_padded = np.pad(self.A, (1, 1), 'edge')
        U_padded = np.zeros((nv, self.nx + 2))
        U_padded[:, 1:-1] = self.U
        e_amb = self.p_amb / ((GAMMA - 1) * self.rho_amb)
        self._apply_bc(U_padded, A_padded, rho, u, e_amb)
        Flux_interface = self._compute_fluxes(U_padded, A_padded)
        dF = Flux_interface[:, :-1] - Flux_interface[:, 1:]
        S = np.zeros_like(self.U)
        S[1, :] = p * self.dA_dx_term
        if prop_params and prop_params['force'] != 0:
            enable = True
            if prop_params['pulse_mode']:
                cycle_time = self.t % (prop_params['period_on'] + prop_params['period_off'])
                if cycle_time > prop_params['period_on']:
                    enable = False
            if enable:
                mask = (self.x > prop_params['x_start']) & (self.x < prop_params['x_end'])
                if np.any(mask):
                    f_density = prop_params['force'] / (prop_params['x_end'] - prop_params['x_start'])
                    S[1, mask] += f_density
                    S[2, mask] += f_density * u[mask]
        if self.combustion:
            S += self._combustion_source(dt)
        S += self._jet_damping_source(u)
        self.U += (dt / self.dx) * dF + dt * S
        # E4: Track thrust
        rho_n, u_n, p_n, _, _ = self.get_primitive()
        thrust = self._compute_thrust(rho_n, u_n, p_n)
        self.hist_t.append(self.t)
        self.hist_thrust.append(thrust)
        return dt

    def print_profiles(self, tag="final"):
        rho_f, u_f, p_f, T_f, Y_fuel = self.get_primitive()
        d_f = 2.0 * np.sqrt(self.A / np.pi)
        if Y_fuel is not None:
            data = np.column_stack([self.x, d_f, self.A, rho_f, u_f, p_f, T_f, Y_fuel])
            header = "x[m]   D[m]   A[m2]   rho[kg/m3]   u[m/s]   p[Pa]   T[K]   Y_fuel"
        else:
            data = np.column_stack([self.x, d_f, self.A, rho_f, u_f, p_f, T_f])
            header = "x[m]   D[m]   A[m2]   rho[kg/m3]   u[m/s]   p[Pa]   T[K]"
        print(f"\nProfiles ({tag}):")
        print(header)
        np.set_printoptions(linewidth=200)
        print(data)
        np.set_printoptions()
