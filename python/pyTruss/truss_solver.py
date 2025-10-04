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
        else:
            self.rest_vectors = np.zeros((0, 3), dtype=np.float64)

        inv_dt2 = 1.0 / (self.dt * self.dt)
        self.mass_diag = self.masses * inv_dt2
        mass_eps = 1e-12
        self.mass_diag = np.where(self.mass_diag > 0.0, self.mass_diag, mass_eps)

        if self.fixed.size > 0:
            self.fixed_positions = self.x[self.fixed].copy()
        else:
            self.fixed_positions = None

        self._linear_matrix: Optional[np.ndarray] = None
        self._linear_diag: Optional[np.ndarray] = None
        self._cholesky_factor: Optional[np.ndarray] = None

        self.solver_state: Dict[str, Any] = {'owner': self}

    def reset_state(self, *, positions: Optional[np.ndarray] = None,
                    velocities: Optional[np.ndarray] = None) -> None:
        """Override stored positions/velocities for a fresh run."""
        print("TrussSolver.reset_state()")
        if positions is not None:
            self.x = positions.astype(np.float64, copy=True)
        if velocities is not None:
            self.v = velocities.astype(np.float64, copy=True)
        if self.trajectory is not None and self.track_idx is not None:
            self.trajectory = [self.x[self.track_idx].copy()]
        if self.fixed.size > 0:
            self.fixed_positions = self.x[self.fixed].copy()
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
            rest_vec = self.rest_vectors[idx]
            si = slice(3 * i, 3 * i + 3)
            sj = slice(3 * j, 3 * j + 3)
            rhs[si] += k * rest_vec
            rhs[sj] -= k * rest_vec

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

    def run(self, nsteps: int) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]:
        print("TrussSolver.run()")
        dt = self.dt
        gravity = self.gravity
        solver = self.solver_callback
        config = self.solver_config
        truss = self.truss

        x = self.x
        v = self.v
        fixed = self.fixed
        solver_state = self.solver_state
        trajectory = self.trajectory
        track_idx = self.track_idx

        for step in range(nsteps):
            x_new = solver(truss, state=solver_state, dt=dt, gravity=gravity,
                           x=x, v=v, fixed=fixed, config=config)
            v = (x_new - x) / dt
            x = x_new
            if len(fixed) > 0:
                v[fixed] = 0.0
            if trajectory is not None and track_idx is not None:
                trajectory.append(x[track_idx].copy())
            if self.verbose and (step == 0 or (step + 1) % max(1, nsteps // 10) == 0):
                print(f"  Step {step + 1}/{nsteps}, t={dt * (step + 1):.3f} s")
            if VERBOSITY >= 1:
                _print_state(f"[CPU][step {step + 1}/{nsteps}]", x, v)

        self.x = x
        self.v = v

        traj_out = None
        if trajectory is not None:
            traj_out = np.stack(trajectory, axis=0)

        return x, v, traj_out


def solve_vbd(truss, *, state: dict, dt: float, gravity: np.ndarray,
              x: np.ndarray, v: np.ndarray, fixed: np.ndarray,
              config: dict) -> np.ndarray:
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

    bonds = np.array(truss.bonds, dtype=np.int32)
    masses = truss.masses.astype(np.float64, copy=False)
    ks = truss.ks.astype(np.float64, copy=False)
    l0s = truss.get_rest_lengths().astype(np.float64, copy=False)

    n_points = x.shape[0]
    dt_sq = dt * dt
    inv_h2 = 1.0 / dt_sq
    I3 = np.eye(3, dtype=np.float64)

    fixed_mask = np.zeros(n_points, dtype=bool)
    if fixed.size > 0:
        fixed_mask[fixed] = True

    x_curr = x.copy()
    y = x + dt * v + dt_sq * gravity

    if fixed.size > 0:
        x_curr[fixed_mask] = x[fixed_mask]
        y[fixed_mask] = x[fixed_mask]

    eye_mass = np.einsum('i,ab->iab', masses * inv_h2, I3)

    for itr in range(niter):
        if fixed.size > 0:
            x_curr[fixed_mask] = x[fixed_mask]

        xi = x_curr[bonds[:, 0]]
        xj = x_curr[bonds[:, 1]]
        d = xi - xj
        L = np.linalg.norm(d, axis=1)
        safe_L = np.maximum(L, 1e-9)
        dir_vec = np.zeros_like(d)
        mask_nonzero = L > 1e-9
        dir_vec[mask_nonzero] = d[mask_nonzero] / safe_L[mask_nonzero, None]

        stretch = L - l0s
        grad_edge = ks[:, None] * stretch[:, None] * dir_vec

        grad = (masses * inv_h2)[:, None] * (x_curr - y)
        np.add.at(grad, bonds[:, 0], grad_edge)
        np.add.at(grad, bonds[:, 1], -grad_edge)

        coeff_iso = ks * (1.0 - l0s / safe_L)
        coeff_dir = ks * (l0s / safe_L)
        coeff_iso[~mask_nonzero] = ks[~mask_nonzero]
        coeff_dir[~mask_nonzero] = ks[~mask_nonzero]
        outer = dir_vec[:, :, None] * dir_vec[:, None, :]
        H_edges = coeff_iso[:, None, None] * I3 + coeff_dir[:, None, None] * outer

        H = eye_mass.copy()
        np.add.at(H, bonds[:, 0], H_edges)
        np.add.at(H, bonds[:, 1], H_edges)

        grad[fixed_mask] = 0.0
        if fixed.size > 0:
            H[fixed_mask] = eye_mass[fixed_mask]

        max_dx = 0.0
        for vid in range(n_points):
            if fixed_mask[vid]:
                continue
            Hi = H[vid]
            gi = grad[vid]
            det = np.linalg.det(Hi)
            if abs(det) < det_eps:
                continue
            dx = -np.linalg.solve(Hi, gi)
            x_curr[vid] += dx
            step = np.linalg.norm(dx)
            if step > max_dx:
                max_dx = step

        if verbose:
            print(f"    VBD iter {itr}: max |dx| = {max_dx:.3e}")
        if VERBOSITY >= 2:
            vel_iter = (x_curr - x) / dt
            if fixed.size > 0:
                vel_iter[fixed_mask] = 0.0
            _print_state(f"[CPU][VBD iter {itr + 1}/{niter}]", x_curr, vel_iter)

    if fixed.size > 0:
        x_curr[fixed_mask] = x[fixed_mask]

    return x_curr


def _require_owner(state: dict) -> TrussSolver:
    print("truss_solver._require_owner()")
    owner = state.get('owner')
    if owner is None or not isinstance(owner, TrussSolver):
        raise ValueError("solver state missing TrussSolver owner; ensure TrussSolver.run() provided the state")
    return owner


def solve_direct(truss, *, state: dict, dt: float, gravity: np.ndarray,
                 x: np.ndarray, v: np.ndarray, fixed: np.ndarray,
                 config: dict) -> np.ndarray:
    print("truss_solver.solve_direct()")
    solver = _require_owner(state)
    y = x + dt * v + (dt * dt) * gravity
    rhs = solver._assemble_linear_rhs(y)
    matrix = solver._get_linear_matrix()
    flat = np.linalg.solve(matrix, rhs)
    solver._enforce_fixed_flat(flat)
    result = flat.reshape(-1, 3)
    state['last_flat'] = flat.copy()
    return result


def solve_cholesky(truss, *, state: dict, dt: float, gravity: np.ndarray,
                   x: np.ndarray, v: np.ndarray, fixed: np.ndarray,
                   config: dict) -> np.ndarray:
    print("truss_solver.solve_cholesky()")
    solver = _require_owner(state)
    y = x + dt * v + (dt * dt) * gravity
    rhs = solver._assemble_linear_rhs(y)
    L = solver._get_cholesky_factor()
    tmp = np.linalg.solve(L, rhs)
    flat = np.linalg.solve(L.T, tmp)
    solver._enforce_fixed_flat(flat)
    result = flat.reshape(-1, 3)
    state['last_flat'] = flat.copy()
    return result


def solve_iterative_momentum(truss, *, state: dict, dt: float, gravity: np.ndarray,
                             x: np.ndarray, v: np.ndarray, fixed: np.ndarray,
                             config: dict) -> np.ndarray:
    print("truss_solver.solve_iterative_momentum()")
    solver = _require_owner(state)
    matrix = solver._get_linear_matrix()
    diag = solver._get_linear_diag()
    y = x + dt * v + (dt * dt) * gravity
    rhs = solver._assemble_linear_rhs(y)

    flat_prev = state.get('last_flat')
    if flat_prev is None or flat_prev.shape[0] != solver.ndof:
        flat = x.reshape(-1).astype(np.float64, copy=True)
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
    result = flat.reshape(-1, 3)
    return result


# Solver registry
SOLVERS = {
    'vbd': solve_vbd,
    'direct': solve_direct,
    'cholesky': solve_cholesky,
    'momentum': solve_iterative_momentum,
}


def get_solver(name: str) -> Callable:
    """Retrieve solver callback by name."""
    print("truss_solver.get_solver()")
    if name not in SOLVERS:
        raise ValueError(f"Unknown solver '{name}'; available: {list(SOLVERS.keys())}")
    return SOLVERS[name]
