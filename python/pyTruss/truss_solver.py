"""
Truss dynamics solver module.

Provides implicit Euler time-stepping with pluggable solver backends:
- Vertex Block Descent (VBD)
- Projective Dynamics (dense/iterative)
- Jacobi/Gauss-Seidel
"""

import numpy as np
from typing import List, Tuple, Callable, Optional, Dict, Any


VERBOSITY = 0


def set_verbosity(level: int) -> None:
    global VERBOSITY
    VERBOSITY = int(level)


def _apply_fixed(x: np.ndarray, fixed_idx: np.ndarray, fixed_vals: np.ndarray) -> None:
    if fixed_idx.size == 0:
        return
    x[fixed_idx] = fixed_vals

def _print_state(tag: str, positions: np.ndarray, velocities: np.ndarray) -> None:
    print(f"{tag} positions:\n{positions}")
    print(f"{tag} velocities:\n{velocities}")

class TrussSolver:
    """
    Implicit Euler dynamics driver with pluggable solver callbacks.

    Stores positions, velocities, and simulation parameters to avoid repeatedly
    passing large arrays between functions.
    """

    def __init__(self, truss, *, dt: float, gravity: np.ndarray,
                 solver: Callable,
                 solver_config: Optional[Dict[str, Any]] = None,
                 fixed_points: Optional[List[int]] = None,
                 track_indices: Optional[List[int]] = None,
                 verbose: int = 0) -> None:
        print("TrussSolver.__init__()")
        self.truss = truss
        self.dt = float(dt)
        self.gravity = np.asarray(gravity, dtype=np.float64)
        self.solver_callback = solver
        self.solver_config = dict(solver_config) if solver_config is not None else {}
        self.fixed = np.array(sorted(truss.fixed if fixed_points is None else fixed_points), dtype=int)
        self.verbose = int(verbose)

        self.x = truss.points.astype(np.float64, copy=True)
        self.v = np.zeros_like(self.x)

        if track_indices is not None and len(track_indices) > 0:
            self.track_idx = np.asarray(track_indices, dtype=int)
            self.trajectory: Optional[list[np.ndarray]] = [self.x[self.track_idx].copy()]
        else:
            self.track_idx = None
            self.trajectory = None

        self.n_points = self.x.shape[0]
        self.ndof = self.n_points * 3
        self.bonds = np.asarray(truss.bonds, dtype=np.int32)
        if self.bonds.size == 0:
            self.bonds = self.bonds.reshape(0, 2)
        self.ks = truss.ks.astype(np.float64, copy=False)
        if self.ks.ndim == 0:
            self.ks = np.full(self.bonds.shape[0], float(self.ks))
        self.masses = truss.masses.astype(np.float64, copy=False)
        self.rest_positions = self.x.copy()
        if self.bonds.size > 0:
            rest_i = self.rest_positions[self.bonds[:, 0]]
            rest_j = self.rest_positions[self.bonds[:, 1]]
            self.rest_vectors = rest_j - rest_i
            self.rest_lengths = np.linalg.norm(self.rest_vectors, axis=1)
        else:
            self.rest_vectors = np.zeros((0, 3), dtype=np.float64)
            self.rest_lengths = np.zeros(0, dtype=np.float64)

        inv_dt2 = 1.0 / (self.dt * self.dt)
        self.mass_diag = self.masses * inv_dt2
        mass_eps = 1e-12
        self.mass_diag = np.where(self.mass_diag > 0.0, self.mass_diag, mass_eps)

        if self.fixed.size > 0:
            self.fixed_positions = self.x[self.fixed].copy()
        else:
            self.fixed_positions = None

        # Lazy caches for fly solvers
        self._max_neighs: Optional[int] = None
        self._neigh_indices: Optional[np.ndarray] = None
        self._neigh_k: Optional[np.ndarray] = None
        self._neigh_l0: Optional[np.ndarray] = None
        self._linear_matrix: Optional[np.ndarray] = None
        self._linear_diag: Optional[np.ndarray] = None
        self._cholesky_factor: Optional[np.ndarray] = None

        self.solver_state: Dict[str, Any] = {'owner': self}
        self.ps_pred = np.zeros_like(self.x)
        self.ps_cor = np.zeros_like(self.x)
        self.forces = np.zeros_like(self.x)

    def _build_linear_diff_system(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        self._ensure_fly_neighbors()
        neighs = self._neigh_indices
        kngs = self._neigh_k
        l0ngs = self._neigh_l0
        max_neighs = self._max_neighs or 0

        rhs = np.zeros_like(self.ps_pred)
        diag = self.mass_diag.copy()
        for vid in range(self.n_points):
            pi = self.ps_pred[vid]
            bi = np.zeros(3, dtype=np.float64)
            Aii = self.mass_diag[vid]
            for jj in range(max_neighs):
                j = neighs[vid, jj]
                if j < 0:
                    break
                k = kngs[vid, jj]
                rest = l0ngs[vid, jj]
                pj = self.ps_pred[j]
                dij = pi - pj
                length = np.linalg.norm(dij)
                safe = max(length, 1e-9)
                bi += dij * k * (rest / safe - 1.0)
                Aii += k
            rhs[vid] = bi
            diag[vid] = Aii

        fixed_mask = np.zeros(self.n_points, dtype=bool)
        if self.fixed.size > 0:
            fixed_mask[self.fixed] = True
            rhs[fixed_mask] = 0.0
            diag[fixed_mask] = 1.0

        return rhs, diag, fixed_mask

    def reset_state(self, *, positions: Optional[np.ndarray] = None,
                    velocities: Optional[np.ndarray] = None) -> None:
        print("TrussSolver.reset_state()")
        if positions is not None:
            self.x = positions.astype(np.float64, copy=True)
        if velocities is not None:
            self.v = velocities.astype(np.float64, copy=True)
        if self.trajectory is not None and self.track_idx is not None:
            self.trajectory = [self.x[self.track_idx].copy()]
        if self.fixed.size > 0:
            self.fixed_positions = self.x[self.fixed].copy()
        else:
            self.fixed_positions = None
        self._max_neighs = None
        self._neigh_indices = None
        self._neigh_k = None
        self._neigh_l0 = None
        self._invalidate_linear_cache()
        self.solver_state.clear()
        self.solver_state['owner'] = self

    def _invalidate_linear_cache(self) -> None:
        print("TrussSolver._invalidate_linear_cache()")
        self._linear_matrix = None
        self._linear_diag = None
        self._cholesky_factor = None

    def _get_linear_matrix(self) -> np.ndarray:
        print("TrussSolver._get_linear_matrix()")
        if self._linear_matrix is not None:
            return self._linear_matrix

        n = self.ndof
        A = np.zeros((n, n), dtype=np.float64)
        eye3 = np.eye(3, dtype=np.float64)

        for vid in range(self.n_points):
            sl = slice(3 * vid, 3 * vid + 3)
            A[sl, sl] += self.mass_diag[vid] * eye3

        for idx, (i, j) in enumerate(self.bonds):
            k = float(self.ks[idx]) if self.ks.size > 0 else 0.0
            if k == 0.0:
                continue
            block = k * eye3
            si = slice(3 * i, 3 * i + 3)
            sj = slice(3 * j, 3 * j + 3)
            A[si, si] += block
            A[sj, sj] += block
            A[si, sj] -= block
            A[sj, si] -= block

        if self.fixed.size > 0:
            for vid in self.fixed:
                sl = slice(3 * vid, 3 * vid + 3)
                A[sl, :] = 0.0
                A[:, sl] = 0.0
                A[sl, sl] = eye3

        reg = 1e-9
        A[np.diag_indices_from(A)] += reg

        self._linear_matrix = A
        self._linear_diag = np.diag(A).copy()
        return self._linear_matrix

    def _get_linear_diag(self) -> np.ndarray:
        print("TrussSolver._get_linear_diag()")
        if self._linear_diag is None:
            self._get_linear_matrix()
        return self._linear_diag

    def _get_cholesky_factor(self) -> np.ndarray:
        print("TrussSolver._get_cholesky_factor()")
        if self._cholesky_factor is None:
            matrix = self._get_linear_matrix()
            self._cholesky_factor = np.linalg.cholesky(matrix)
        return self._cholesky_factor

    def _assemble_linear_rhs(self, y: np.ndarray) -> np.ndarray:
        print("TrussSolver._assemble_linear_rhs()")
        rhs = np.zeros(self.ndof, dtype=np.float64)
        
        for vid in range(self.n_points):
            sl = slice(3 * vid, 3 * vid + 3)
            rhs[sl] += self.mass_diag[vid] * y[vid]

        for idx, (i, j) in enumerate(self.bonds):
            k = float(self.ks[idx]) if self.ks.size > 0 else 0.0
            if k == 0.0:
                continue
            # Match C++ rhs_ProjectiveDynamics_i: d = pnew[i] - pnew[j]; bi += k * l0 * d / |d|
            d = y[i] - y[j]
            length = np.linalg.norm(d)
            l0 = np.linalg.norm(self.rest_vectors[idx])
            if length > 1e-12:
                direction = d / length
                spring_contrib = k * l0 * direction
            else:
                spring_contrib = np.zeros(3)
            si = slice(3 * i, 3 * i + 3)
            sj = slice(3 * j, 3 * j + 3)
            rhs[si] += spring_contrib
            rhs[sj] -= spring_contrib
        
        if self.fixed.size > 0 and self.fixed_positions is not None:
            for pos_idx, vid in enumerate(self.fixed):
                sl = slice(3 * vid, 3 * vid + 3)
                rhs[sl] = self.fixed_positions[pos_idx]

        return rhs

    def _enforce_fixed_flat(self, flat: np.ndarray) -> None:
        print("TrussSolver._enforce_fixed_flat()")
        if self.fixed.size == 0 or self.fixed_positions is None:
            return
        for pos_idx, vid in enumerate(self.fixed):
            sl = slice(3 * vid, 3 * vid + 3)
            flat[sl] = self.fixed_positions[pos_idx]

    def _zero_fixed_flat(self, flat: np.ndarray) -> None:
        print("TrussSolver._zero_fixed_flat()")
        if self.fixed.size == 0:
            return
        for vid in self.fixed:
            sl = slice(3 * vid, 3 * vid + 3)
            flat[sl] = 0.0

    def _ensure_fly_neighbors(self) -> None:
        if self._neigh_indices is not None:
            return
        if self.bonds.size == 0:
            shape = (self.n_points, 0)
            self._max_neighs = 0
            self._neigh_indices = np.zeros(shape, dtype=np.int32)
            self._neigh_k = np.zeros(shape, dtype=np.float64)
            self._neigh_l0 = np.zeros(shape, dtype=np.float64)
            return

        degrees = np.bincount(self.bonds.flatten(), minlength=self.n_points)
        max_neighs = int(degrees.max()) if degrees.size > 0 else 0
        self._max_neighs = max_neighs
        neigh_shape = (self.n_points, max_neighs)
        neigh_indices = np.full(neigh_shape, -1, dtype=np.int32)
        neigh_k = np.zeros(neigh_shape, dtype=np.float64)
        neigh_l0 = np.zeros(neigh_shape, dtype=np.float64)

        for ib, (i, j) in enumerate(self.bonds):
            k_val = float(self.ks[ib])
            l0_val = float(self.rest_lengths[ib]) if self.rest_lengths.size > 0 else 0.0
            for slot in range(max_neighs):
                if neigh_indices[i, slot] < 0:
                    neigh_indices[i, slot] = j
                    neigh_k[i, slot] = k_val
                    neigh_l0[i, slot] = l0_val
                    break
            for slot in range(max_neighs):
                if neigh_indices[j, slot] < 0:
                    neigh_indices[j, slot] = i
                    neigh_k[j, slot] = k_val
                    neigh_l0[j, slot] = l0_val
                    break

        self._neigh_indices = neigh_indices
        self._neigh_k = neigh_k
        self._neigh_l0 = neigh_l0

    def _prepare_step(self) -> None:
        dt = self.dt
        dt_sq = dt * dt
        np.copyto(self.ps_pred, self.x + dt * self.v + dt_sq * self.gravity)
        if self.fixed.size > 0:
            self.ps_pred[self.fixed] = self.x[self.fixed]
        self.forces.fill(0.0)

    def _correct_step(self, damping: float) -> None:
        inv_dt = 1.0 / self.dt
        cdamp = max(1.0 - damping * self.dt, 0.0)
        if self.fixed.size > 0:
            fixed_mask = np.zeros(self.n_points, dtype=bool)
            fixed_mask[self.fixed] = True
            free_mask = ~fixed_mask
            delta = self.ps_cor - self.x
            self.v[free_mask] = delta[free_mask] * inv_dt * cdamp
            self.x[free_mask] = self.ps_cor[free_mask]
            self.v[self.fixed] = 0.0
        else:
            delta = self.ps_cor - self.x
            self.v[:] = delta * inv_dt * cdamp
            self.x[:] = self.ps_cor


    def run(self, nsteps: int) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]:
        print("TrussSolver.run()")
        solver = self.solver_callback
        config = self.solver_config
        trajectory = self.trajectory
        track_idx = self.track_idx
        damping = float(config.get('damping', 0.0))

        for step in range(nsteps):
            self._prepare_step()
            solver(self, config)
            self._correct_step(damping)
            if trajectory is not None and track_idx is not None:
                trajectory.append(self.x[track_idx].copy())
            if self.verbose and (step == 0 or (step + 1) % max(1, nsteps // 10) == 0):
                print(f"  Step {step + 1}/{nsteps}, t={self.dt * (step + 1):.3f} s")
            if VERBOSITY >= 1:
                _print_state(f"[CPU][step {step + 1}/{nsteps}]", self.x, self.v)

        traj_out = None
        if trajectory is not None:
            traj_out = np.stack(trajectory, axis=0)

        return self.x, self.v, traj_out


def solve_vbd(solver: TrussSolver, config: dict) -> None:
    """
    Vertex Block Descent solver callback for one implicit Euler step.

    Config keys:
        niter: Number of VBD iterations (default 50)
        det_eps: Determinant threshold (default 1e-6)
        verbose: Print iteration diagnostics (default 0)
    """
    print("truss_solver.solve_vbd()")
    niter = config.get('niter', 50)
    det_eps = config.get('det_eps', 1e-6)
    verbose = config.get('verbose', 0)

    solver._ensure_fly_neighbors()

    neighs = solver._neigh_indices
    kngs = solver._neigh_k
    l0ngs = solver._neigh_l0
    max_neighs = solver._max_neighs or 0

    n_points = solver.n_points
    mass_diag = solver.mass_diag
    fixed_mask = np.zeros(n_points, dtype=bool)
    if solver.fixed.size > 0:
        fixed_mask[solver.fixed] = True

    x_curr = solver.x.copy()
    y = solver.ps_pred.copy()

    I3 = np.eye(3, dtype=np.float64)

    for itr in range(niter):
        if solver.fixed.size > 0:
            x_curr[fixed_mask] = solver.x[fixed_mask]

        max_dx = 0.0

        for vid in range(n_points):
            if fixed_mask[vid]:
                continue

            pi = x_curr[vid]
            yi = y[vid]
            mass_term = mass_diag[vid]

            grad = mass_term * (pi - yi)
            H = mass_term * I3.copy()

            for jj in range(max_neighs):
                j = neighs[vid, jj]
                if j < 0:
                    break
                k = kngs[vid, jj]
                rest = l0ngs[vid, jj]
                pj = x_curr[j]
                d = pi - pj
                length = np.linalg.norm(d)
                if length > 1e-9:
                    dir_vec = d / length
                    stretch = length - rest
                else:
                    dir_vec = np.zeros(3, dtype=np.float64)
                    length = 1e-9
                    stretch = -rest
                grad += k * stretch * dir_vec
                coeff_iso = k * (1.0 - rest / length)
                coeff_dir = k * (rest / length)
                H += coeff_iso * I3 + coeff_dir * np.outer(dir_vec, dir_vec)

            det = np.linalg.det(H)
            if abs(det) < det_eps:
                continue

            dx = -np.linalg.solve(H, grad)
            x_curr[vid] += dx
            step = np.linalg.norm(dx)
            if step > max_dx:
                max_dx = step

        if verbose:
            print(f"    VBD iter {itr}: max |dx| = {max_dx:.3e}")
        if VERBOSITY >= 2:
            vel_iter = (x_curr - solver.x) / solver.dt
            if solver.fixed.size > 0:
                vel_iter[fixed_mask] = 0.0
            _print_state(f"[CPU][VBD iter {itr + 1}/{niter}]", x_curr, vel_iter)

    if solver.fixed.size > 0:
        x_curr[fixed_mask] = solver.x[fixed_mask]

    solver.ps_cor[:] = x_curr


def _require_owner(state: dict) -> TrussSolver:
    print("truss_solver._require_owner()")
    owner = state.get('owner')
    if owner is None or not isinstance(owner, TrussSolver):
        raise ValueError("solver state missing TrussSolver owner; ensure TrussSolver.run() provided the state")
    return owner


def solve_direct(solver: TrussSolver, config: dict) -> None:
    print("truss_solver.solve_direct()")
    rhs = solver._assemble_linear_rhs(solver.ps_pred)
    matrix = solver._get_linear_matrix()
    flat = np.linalg.solve(matrix, rhs)
    solver._enforce_fixed_flat(flat)
    solver.ps_cor[:] = flat.reshape(-1, 3)
    solver.solver_state['last_flat'] = flat.copy()


def solve_cholesky(solver: TrussSolver, config: dict) -> None:
    print("truss_solver.solve_cholesky()")
    rhs = solver._assemble_linear_rhs(solver.ps_pred)
    L = solver._get_cholesky_factor()
    tmp = np.linalg.solve(L, rhs)
    flat = np.linalg.solve(L.T, tmp)
    solver._enforce_fixed_flat(flat)
    solver.ps_cor[:] = flat.reshape(-1, 3)
    solver.solver_state['last_flat'] = flat.copy()


def solve_iterative_momentum(solver: TrussSolver, config: dict) -> None:
    print("truss_solver.solve_iterative_momentum()")
    matrix = solver._get_linear_matrix()
    diag = solver._get_linear_diag()
    rhs = solver._assemble_linear_rhs(solver.ps_pred)

    state = solver.solver_state
    flat_prev = state.get('last_flat')
    if flat_prev is None or flat_prev.shape[0] != solver.ndof:
        flat = solver.ps_pred.reshape(-1).astype(np.float64, copy=True)
    else:
        flat = flat_prev.copy()
    solver._enforce_fixed_flat(flat)

    momentum = state.get('momentum')
    if momentum is None or momentum.shape[0] != solver.ndof:
        momentum = np.zeros(solver.ndof, dtype=np.float64)

    niter = int(config.get('niter', 10))
    beta = float(config.get('momentum_beta', config.get('momentum', 0.85)))
    beta = max(0.0, min(beta, 0.999))
    relax = float(config.get('relaxation', 1.0))
    relax = max(1e-6, relax)

    for itr in range(max(1, niter)):
        residual = rhs - matrix @ flat
        solver._zero_fixed_flat(residual)
        delta = (relax * residual) / diag
        momentum = beta * momentum + delta
        solver._zero_fixed_flat(momentum)
        flat += momentum
        solver._enforce_fixed_flat(flat)
        if config.get('tol'):
            tol = float(config['tol'])
            if np.linalg.norm(residual, ord=np.inf) < tol:
                break

    state['momentum'] = momentum
    state['last_flat'] = flat.copy()
    solver.ps_cor[:] = flat.reshape(-1, 3)

def solve_jacobi_diff(solver: TrussSolver, config: dict) -> None:
    """Jacobi-Diff: precomputed linear RHS, iterate on displacements."""
    rhs, diag, fixed_mask = solver._build_linear_diff_system()
    niter = int(config.get('niter', 20))
    solver._ensure_fly_neighbors()
    neighs = solver._neigh_indices
    kngs = solver._neigh_k
    max_neighs = solver._max_neighs or 0

    dp_in = np.zeros_like(solver.ps_pred)
    dp_out = np.zeros_like(dp_in)

    for it in range(niter):
        max_step = 0.0
        for vid in range(solver.n_points):
            if fixed_mask[vid]:
                dp_out[vid] = 0.0
                continue
            sum_j = np.zeros(3, dtype=np.float64)
            for jj in range(max_neighs):
                j = neighs[vid, jj]
                if j < 0:
                    break
                sum_j += kngs[vid, jj] * dp_in[j]
            new_dp = (rhs[vid] + sum_j) / max(diag[vid], 1e-12)
            step = np.linalg.norm(new_dp - dp_in[vid])
            if step > max_step:
                max_step = step
            dp_out[vid] = new_dp
        dp_in, dp_out = dp_out, dp_in
        if VERBOSITY >= 2:
            print(f"[CPU][JacobiDiff iter {it + 1}/{niter}] max |Δdp| = {max_step:.3e}")

    solver.ps_cor[:] = solver.ps_pred + dp_in
    if solver.fixed.size > 0:
        solver.ps_cor[solver.fixed] = solver.x[solver.fixed]


def solve_jacobi_fly(solver: TrussSolver, config: dict) -> None:
    """Jacobi-Fly: Copy of TrussDynamics_d::updateJacobi_fly() and OpenCL jacobi_fly kernel."""
    niter = int(config.get('niter', 10))
    solver._ensure_fly_neighbors()
    max_neighs = solver._max_neighs or 0
    neighs = solver._neigh_indices
    kngs = solver._neigh_k
    l0ngs = solver._neigh_l0

    ps_in = solver.ps_pred.copy()
    ps_out = ps_in.copy()
    
    # Match C++ updateJacobi_fly: for i in nPoint: sum_j = Ii*pi; for jj in nNeighMax: sum_j += k*l0*dij/|dij| + k*pj; ps_out[i] = sum_j/Aii
    for it in range(niter):
        max_step = 0.0
        for i in range(solver.n_points):
            pi = ps_in[i]
            Ii = solver.mass_diag[i]
            sum_j = pi * Ii  # bi = M_i/dt^2 * p_i
            Aii = Ii
            for jj in range(max_neighs):
                jG = neighs[i, jj]
                if jG < 0:
                    break
                k = kngs[i, jj]
                l0_ij = l0ngs[i, jj]
                pj = ps_in[jG]
                dij = pi - pj
                length = np.linalg.norm(dij)
                safe = max(length, 1e-8)
                sum_j += dij * (k * l0_ij / safe)  # sum_j.add_mul(dij, k*par.x/l)
                sum_j += pj * k                    # sum_j.add_mul(pj, k)
                Aii += k
            new_pi = sum_j / max(Aii, 1e-12)
            step = np.linalg.norm(new_pi - pi)
            if step > max_step:
                max_step = step
            ps_out[i] = new_pi
        if solver.fixed.size > 0:
            ps_out[solver.fixed] = solver.x[solver.fixed]
        ps_in, ps_out = ps_out, ps_in
        if VERBOSITY >= 2:
            print(f"[CPU][JacobiFly iter {it + 1}/{niter}] max |Δp| = {max_step:.3e}")
    
    solver.ps_cor[:] = ps_in


def solve_gs_diff(solver: TrussSolver, config: dict) -> None:
    """GS-Diff: precomputed linear RHS, in-place Gauss-Seidel on displacements."""
    rhs, diag, fixed_mask = solver._build_linear_diff_system()
    niter = int(config.get('niter', 20))
    solver._ensure_fly_neighbors()
    neighs = solver._neigh_indices
    kngs = solver._neigh_k
    max_neighs = solver._max_neighs or 0

    dp = np.zeros_like(solver.ps_pred)

    for it in range(niter):
        max_step = 0.0
        for vid in range(solver.n_points):
            if fixed_mask[vid]:
                dp[vid] = 0.0
                continue
            sum_j = np.zeros(3, dtype=np.float64)
            for jj in range(max_neighs):
                j = neighs[vid, jj]
                if j < 0:
                    break
                sum_j += kngs[vid, jj] * dp[j]
            new_dp = (rhs[vid] + sum_j) / max(diag[vid], 1e-12)
            step = np.linalg.norm(new_dp - dp[vid])
            if step > max_step:
                max_step = step
            dp[vid] = new_dp
        if VERBOSITY >= 2:
            print(f"[CPU][GSDiff iter {it + 1}/{niter}] max |Δdp| = {max_step:.3e}")

    solver.ps_cor[:] = solver.ps_pred + dp
    if solver.fixed.size > 0:
        solver.ps_cor[solver.fixed] = solver.x[solver.fixed]


def solve_gs_fly(solver: TrussSolver, config: dict) -> None:
    niter = int(config.get('niter', 10))
    solver._ensure_fly_neighbors()
    max_neighs = solver._max_neighs or 0
    neighs = solver._neigh_indices
    kngs = solver._neigh_k
    l0ngs = solver._neigh_l0

    ps = solver.ps_pred.copy()
    
    # Match C++ updateGaussSeidel_fly: in-place update
    for it in range(niter):
        max_step = 0.0
        for i in range(solver.n_points):
            pi = ps[i]
            Ii = solver.mass_diag[i]
            sum_j = pi * Ii
            Aii = Ii
            for jj in range(max_neighs):
                jG = neighs[i, jj]
                if jG < 0:
                    break
                k = kngs[i, jj]
                l0_ij = l0ngs[i, jj]
                pj = ps[jG]
                dij = pi - pj
                length = np.linalg.norm(dij)
                safe = max(length, 1e-8)
                sum_j += dij * (k * l0_ij / safe)
                sum_j += pj * k
                Aii += k
            new_pi = sum_j / max(Aii, 1e-12)
            step = np.linalg.norm(new_pi - pi)
            if step > max_step:
                max_step = step
            ps[i] = new_pi
        if VERBOSITY >= 2:
            print(f"[CPU][GSFly iter {it + 1}/{niter}] max |Δp| = {max_step:.3e}")
    
    solver.ps_cor[:] = ps
    if solver.fixed.size > 0:
        solver.ps_cor[solver.fixed] = solver.x[solver.fixed]


# Solver registry
SOLVERS = {
    'vbd': solve_vbd,
    'vbd_serial': solve_vbd,
    'direct': solve_direct,
    'cholesky': solve_cholesky,
    'momentum': solve_iterative_momentum,
    'jacobi_diff': solve_jacobi_diff,
    'jacobi_fly': solve_jacobi_fly,
    'gs_diff': solve_gs_diff,
    'gs_fly': solve_gs_fly,
}


def get_solver(name: str) -> Callable:
    """Retrieve solver callback by name."""
    print("truss_solver.get_solver()")
    if name not in SOLVERS:
        raise ValueError(f"Unknown solver '{name}'; available: {list(SOLVERS.keys())}")
    return SOLVERS[name]
