import numpy as np
from tube_common import GAMMA, R_GAS

CFL = 0.5

class CompressibleSolver:
    def __init__(self, nx, length, geom_segments, u0, u_amb, p_init, rho_init, bc_mode="ambient"):
        self.nx = nx
        self.dx = length / nx
        self.t = 0.0
        self.u_amb = u_amb
        self.bc_mode = bc_mode
        self.x = np.linspace(0.5 * self.dx, length - 0.5 * self.dx, nx)
        self.A, self.A_faces, self.dA_dx_term, self.vol = self._generate_geometry(length, geom_segments)
        self.U = np.zeros((3, nx))
        e_init = p_init / ((GAMMA - 1) * rho_init)
        E_init = rho_init * (e_init + 0.5 * u0**2)
        self.U[0, :] = rho_init * self.A
        self.U[1, :] = rho_init * u0 * self.A
        self.U[2, :] = E_init * self.A
        self.p_amb = p_init
        self.rho_amb = rho_init
        self.T_amb = p_init / (R_GAS * rho_init)

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
        return rho, u, p, T

    def compute_flux_rusanov_well_balanced(self, U_L, U_R, A_L, A_R, A_face):
        rho_L = U_L[0] / A_L
        u_L = U_L[1] / U_L[0]
        p_L = (GAMMA - 1) * ((U_L[2] / A_L) - 0.5 * rho_L * u_L**2)
        E_vol_L = U_L[2] / A_L
        rho_R = U_R[0] / A_R
        u_R = U_R[1] / U_R[0]
        p_R = (GAMMA - 1) * ((U_R[2] / A_R) - 0.5 * rho_R * u_R**2)
        E_vol_R = U_R[2] / A_R
        F_L = np.array([rho_L * u_L, rho_L * u_L**2 + p_L, u_L * (E_vol_L + p_L)]) * A_face
        F_R = np.array([rho_R * u_R, rho_R * u_R**2 + p_R, u_R * (E_vol_R + p_R)]) * A_face
        c_L = np.sqrt(GAMMA * abs(p_L) / (abs(rho_L) + 1e-12))
        c_R = np.sqrt(GAMMA * abs(p_R) / (abs(rho_R) + 1e-12))
        S_max = max(abs(u_L) + c_L, abs(u_R) + c_R)
        U_star_L = U_L * (A_face / A_L)
        U_star_R = U_R * (A_face / A_R)
        F_face = 0.5 * (F_L + F_R) - 0.5 * S_max * (U_star_R - U_star_L)
        return F_face, S_max

    def step(self, dt_fixed=None, prop_params=None):
        rho, u, p, T = self.get_primitive()
        c = np.sqrt(GAMMA * abs(p) / (abs(rho) + 1e-12))
        max_speed = np.max(np.abs(u) + c)
        if max_speed < 1e-5: max_speed = 1e-5
        dt = CFL * self.dx / max_speed
        if dt_fixed: dt = dt_fixed
        self.t += dt
        A_padded = np.pad(self.A, (1, 1), 'edge')
        U_padded = np.zeros((3, self.nx + 2))
        U_padded[:, 1:-1] = self.U
        e_amb = self.p_amb / ((GAMMA - 1) * self.rho_amb)
        if self.bc_mode == "ambient":
            U_padded[0, 0] = self.rho_amb * A_padded[0]
            U_padded[1, 0] = self.rho_amb * self.u_amb * A_padded[0]
            U_padded[2, 0] = self.rho_amb * (e_amb + 0.5 * self.u_amb**2) * A_padded[0]
            u_in = -self.u_amb
            U_padded[0, -1] = self.rho_amb * A_padded[-1]
            U_padded[1, -1] = self.rho_amb * u_in * A_padded[-1]
            U_padded[2, -1] = self.rho_amb * (e_amb + 0.5 * u_in**2) * A_padded[-1]
        else:
            if u[0] > 0:
                U_padded[0, 0] = self.rho_amb * A_padded[0]
                U_padded[1, 0] = self.rho_amb * self.u_amb * A_padded[0]
                U_padded[2, 0] = self.rho_amb * (e_amb + 0.5 * self.u_amb**2) * A_padded[0]
            else:
                U_padded[:, 0] = self.U[:, 0] * (A_padded[0] / self.A[0])
            if u[-1] < 0:
                u_in = -self.u_amb
                U_padded[0, -1] = self.rho_amb * A_padded[-1]
                U_padded[1, -1] = self.rho_amb * u_in * A_padded[-1]
                U_padded[2, -1] = self.rho_amb * (e_amb + 0.5 * u_in**2) * A_padded[-1]
            else:
                rho_g = rho[-1]; u_g = u[-1]; p_g = self.p_amb
                e_g = p_g / ((GAMMA - 1) * rho_g)
                U_padded[0, -1] = rho_g * A_padded[-1]
                U_padded[1, -1] = rho_g * u_g * A_padded[-1]
                U_padded[2, -1] = rho_g * (e_g + 0.5 * u_g**2) * A_padded[-1]
        Flux_interface = np.zeros((3, self.nx + 1))
        for i in range(self.nx + 1):
            U_L = U_padded[:, i]
            U_R = U_padded[:, i + 1]
            A_L = A_padded[i]
            A_R = A_padded[i + 1]
            A_face = self.A_faces[i]
            Flux_interface[:, i], _ = self.compute_flux_rusanov_well_balanced(U_L, U_R, A_L, A_R, A_face)
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
        self.U += (dt / self.dx) * dF + dt * S
        return dt

    def print_profiles(self, tag="final"):
        rho_f, u_f, p_f, T_f = self.get_primitive()
        d_f = 2.0 * np.sqrt(self.A / np.pi)
        data = np.column_stack([self.x, d_f, self.A, rho_f, u_f, p_f, T_f])
        header = "x[m]   D[m]   A[m2]   rho[kg/m3]   u[m/s]   p[Pa]   T[K]"
        print(f"\nProfiles ({tag}):")
        print(header)
        np.set_printoptions(linewidth=200)
        print(data)
        np.set_printoptions()
