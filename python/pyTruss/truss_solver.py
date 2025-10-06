"""
Truss dynamics solver module.

Provides implicit Euler time-stepping with pluggable solver backends:
- Vertex Block Descent (VBD)
- Projective Dynamics (dense/iterative)
- Jacobi/Gauss-Seidel
"""

import csv
import math

import numpy as np
from typing import List, Tuple, Callable, Optional, Dict, Any, Sequence


VERBOSITY = 0


def set_verbosity(level: int) -> None:
    global VERBOSITY
    VERBOSITY = int(level)


def _iter_logging_enabled(config: dict) -> bool:
    return bool(config.get('capture_iters'))


def _ensure_log_base(solver: 'TrussSolver') -> None:
    if 'iter_logs' not in solver.solver_state:
        solver.solver_state['iter_logs'] = {}
    if 'log_base' not in solver.solver_state:
        solver.solver_state['log_base'] = solver.ps_pred.copy()


def _log_iteration_state(solver: 'TrussSolver', config: dict, data: np.ndarray,
                         max_step: Optional[float], *, is_displacement: bool) -> None:
    if not _iter_logging_enabled(config):
        return
    logs = solver.solver_state.setdefault('iter_logs', {})
    label = config.get('log_label', 'solver')
    entry = logs.setdefault(label, {'positions': [], 'max_step': []})
    base = solver.solver_state.get('log_base')
    if base is None:
        base = solver.ps_pred
    positions = base + data if is_displacement else data
    entry['positions'].append(np.array(positions, copy=True))
    if max_step is not None:
        entry['max_step'].append(float(max_step))


def estimate_iterative_spectral_radius(solver: 'TrussSolver', method: str = 'jacobi_diff') -> float:
    """Compute spectral radius of the selected linear iteration matrix (debug helper).

    Parameters
    ----------
    solver : TrussSolver
        Solver instance providing assembled linear system via `_get_linear_matrix()`.
    method : str
        Iterative scheme to analyze. Supported values: 'jacobi', 'jacobi_diff',
        'gs', 'gs_diff'.

    Returns
    -------
    float
        Spectral radius (max |λ|) of the iteration matrix.
    """
    matrix = solver._get_linear_matrix()
    n = matrix.shape[0]
    diag = np.diag(matrix)
    if diag.size != n:
        raise ValueError("diagonal extraction failed for spectral radius estimate")

    if method in ('jacobi', 'jacobi_diff'):
        inv_diag = 1.0 / diag
        iteration = np.eye(n) - inv_diag[:, None] * matrix
    elif method in ('gs', 'gs_diff'):
        D = np.diag(diag)
        L = np.tril(matrix, k=-1)
        U = np.triu(matrix, k=1)
        DL = D + L
        iteration = -np.linalg.solve(DL, U)
    else:
        raise ValueError(f"unsupported method '{method}' for spectral radius estimate")

    eigvals = np.linalg.eigvals(iteration)
    if eigvals.size == 0:
        raise ValueError("eigenvalue computation returned empty result")
    return float(np.max(np.abs(eigvals)))


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

    capture_iters = _iter_logging_enabled(config)
    if capture_iters:
        _ensure_log_base(solver)

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

    if capture_iters:
        _log_iteration_state(solver, config, x_curr, 0.0, is_displacement=False)

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

        if capture_iters:
            _log_iteration_state(solver, config, x_curr, max_dx, is_displacement=False)

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


def _momentum_params(config: dict) -> dict:
    params = {
        'b_start': float(config.get('b_start', 0.55)),
        'b_end': float(config.get('b_end', 0.75)),
        'b_last': float(config.get('b_last', 0.0)),
        'istart': int(config.get('istart', 3)),
        'iend': config.get('iend'),
        'niter': int(config.get('niter', 10)),
    }
    if params['iend'] is None:
        params['iend'] = max(params['istart'] + 1, params['niter'] - 1)
    else:
        params['iend'] = int(params['iend'])
    params['istart'] = max(0, params['istart'])
    params['iend'] = max(params['istart'], params['iend'])
    params['niter'] = max(1, params['niter'])
    return params


def _momentum_mix(params: dict, itr: int) -> float:
    niter = params['niter']
    if itr == 0 or itr >= niter - 1:
        return params['b_last']
    if itr < params['istart']:
        return 0.0
    if itr >= params['iend']:
        return params['b_end']
    span = params['iend'] - params['istart']
    if span <= 0:
        return params['b_end']
    t = (itr - params['istart']) / span
    return params['b_start'] + t * (params['b_end'] - params['b_start'])


def solve_iterative_momentum(solver: TrussSolver, config: dict) -> None:
    """Momentum-accelerated iterative solver. Matches C++ updateIterativeMomentum architecture."""
    print("truss_solver.solve_iterative_momentum()")
    
    # Sub-solver selection
    sub_method = config.get('sub_method', 'jacobi_diff')
    if sub_method not in ['jacobi_diff', 'gs_diff', 'jacobi_fly', 'gs_fly']:
        raise ValueError(f"Invalid sub_method '{sub_method}' for momentum solver")
    
    params = _momentum_params(config)
    niter = params['niter']
    
    state = solver.solver_state
    inertia = state.get('inertia')
    if inertia is None or inertia.shape != solver.ps_pred.shape:
        inertia = np.zeros_like(solver.ps_pred)
    
    solver._ensure_fly_neighbors()

    capture_iters = _iter_logging_enabled(config)
    if capture_iters:
        _ensure_log_base(solver)

    # Pre-compute RHS ONCE for *_lin methods (matches C++ updatePD_dRHS call)
    rhs = diag = fixed_mask = None
    if sub_method in ['jacobi_diff', 'gs_diff']:
        rhs, diag, fixed_mask = solver._build_linear_diff_system()
    
    # For *_diff methods: iterate on displacements dp (psa/psb are dp, not positions)
    # For *_fly methods: iterate on positions directly
    if sub_method in ['jacobi_diff', 'gs_diff']:
        psa = np.zeros_like(solver.ps_pred)  # dp_in
        psb = np.zeros_like(psa)             # dp_out
    else:
        psa = solver.ps_pred.copy()
        psb = np.zeros_like(psa)
    
    is_diff = sub_method in ['jacobi_diff', 'gs_diff']
    if capture_iters:
        _log_iteration_state(solver, config, psa, 0.0, is_displacement=is_diff)

    for i in range(niter):
        # Single iteration (matches C++ updateJacobi_lin, updateGaussSeidel_lin, etc.)
        if sub_method == 'jacobi_diff':
            _update_jacobi_lin(solver, psa, psb, rhs, diag, fixed_mask)
        elif sub_method == 'gs_diff':
            psb[:] = psa  # GS needs copy for in-place update
            _update_gs_lin(solver, psb, rhs, diag, fixed_mask)
        elif sub_method == 'jacobi_fly':
            _update_jacobi_fly(solver, psa, psb)
        elif sub_method == 'gs_fly':
            psb[:] = psa
            _update_gs_fly(solver, psb)
        
        bmix = _momentum_mix(params, i)
        
        # Momentum mixing: p_{k+1} = p'_k + bmix * d_k
        p = psb + bmix * inertia
        np.subtract(p, psa, out=inertia)  # d_{k+1} = p_{k+1} - p_k
        
        max_step = float(np.abs(inertia).max())
        if VERBOSITY >= 2:
            print(f"[CPU][Momentum iter {i+1}/{niter}] bmix={bmix:.2f} max |Δp| = {max_step:.3e}")

        psa[:] = p                # p_k = p_{k+1}
        if capture_iters:
            _log_iteration_state(solver, config, psa, max_step, is_displacement=is_diff)
    
    # Final result
    if sub_method in ['jacobi_diff', 'gs_diff']:
        solver.ps_cor[:] = solver.ps_pred + psa  # position = predictor + displacement
    else:
        solver.ps_cor[:] = psa  # psa already holds positions
    
    if solver.fixed.size > 0:
        solver.ps_cor[solver.fixed] = solver.x[solver.fixed]
    state['inertia'] = inertia


def _ensure_state_array(state: dict, key: str, ref: np.ndarray) -> np.ndarray:
    arr = state.get(key)
    if arr is None or arr.shape != ref.shape:
        arr = ref.copy()
        state[key] = arr
    return arr


def solve_iterative_chebyshev(solver: TrussSolver, config: dict) -> None:
    """Chebyshev-accelerated iterative solver."""
    print("truss_solver.solve_iterative_chebyshev()")

    sub_method = config.get('sub_method', 'jacobi_diff')
    if sub_method not in ['jacobi_diff', 'gs_diff', 'jacobi_fly', 'gs_fly']:
        raise ValueError(f"Invalid sub_method '{sub_method}' for chebyshev solver")

    rho = float(config.get('rho', 0.95))
    delayed = max(0, int(config.get('delayed_start', 5)))
    gamma = float(config.get('gamma', 1.0))
    niter = max(1, int(config.get('niter', 10)))

    state = solver.solver_state

    solver._ensure_fly_neighbors()

    capture_iters = _iter_logging_enabled(config)
    if capture_iters:
        _ensure_log_base(solver)

    rhs = diag = fixed_mask = None
    if sub_method in ['jacobi_diff', 'gs_diff']:
        rhs, diag, fixed_mask = solver._build_linear_diff_system()

    if sub_method in ['jacobi_diff', 'gs_diff']:
        psa = np.zeros_like(solver.ps_pred)
        psb = np.zeros_like(psa)
    else:
        psa = solver.ps_pred.copy()
        psb = np.zeros_like(psa)

    prev = _ensure_state_array(state, 'cheb_prev', psa)
    prev2 = _ensure_state_array(state, 'cheb_prev2', psa)

    omega = 1.0
    is_diff = sub_method in ['jacobi_diff', 'gs_diff']

    if capture_iters:
        _log_iteration_state(solver, config, psa, 0.0, is_displacement=is_diff)

    for itr in range(niter):
        if sub_method == 'jacobi_diff':
            _update_jacobi_lin(solver, psa, psb, rhs, diag, fixed_mask)
        elif sub_method == 'gs_diff':
            psb[:] = psa
            _update_gs_lin(solver, psb, rhs, diag, fixed_mask)
        elif sub_method == 'jacobi_fly':
            _update_jacobi_fly(solver, psa, psb)
        else:  # gs_fly
            psb[:] = psa
            _update_gs_fly(solver, psb)

        if gamma != 1.0:
            # under-relaxation: psa holds previous iterate
            psb[:] = gamma * (psb - psa) + psa

        if itr < delayed:
            omega = 1.0
            p_next = psb.copy()
        else:
            if itr == delayed:
                omega = 2.0 / max(2.0 - rho * rho, 1e-9)
            else:
                denom = 4.0 - rho * rho * omega
                omega = 4.0 / max(denom, 1e-9)
            p_next = prev + omega * (psb - prev2)

        delta = float(np.linalg.norm(p_next - psa, ord=np.inf))
        if VERBOSITY >= 2:
            print(f"[CPU][Chebyshev iter {itr+1}/{niter}] omega={omega:.3f} max |Δp| = {delta:.3e}")

        prev2[:] = prev
        prev[:] = p_next
        psa[:] = p_next
        if capture_iters:
            _log_iteration_state(solver, config, psa, delta, is_displacement=is_diff)

    if sub_method in ['jacobi_diff', 'gs_diff']:
        solver.ps_cor[:] = solver.ps_pred + psa
    else:
        solver.ps_cor[:] = psa

    if solver.fixed.size > 0:
        solver.ps_cor[solver.fixed] = solver.x[solver.fixed]

    state['cheb_prev'] = prev
    state['cheb_prev2'] = prev2

def _update_jacobi_lin(solver: TrussSolver, dp_in: np.ndarray, dp_out: np.ndarray, 
                       rhs: np.ndarray, diag: np.ndarray, fixed_mask: np.ndarray) -> float:
    """Single Jacobi iteration on linearized system. Matches C++ updateJacobi_lin."""
    neighs = solver._neigh_indices
    kngs = solver._neigh_k
    max_neighs = solver._max_neighs or 0
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
    return max_step


def solve_jacobi_diff(solver: TrussSolver, config: dict) -> None:
    """Jacobi-Diff: precomputed linear RHS, iterate on displacements."""
    rhs, diag, fixed_mask = solver._build_linear_diff_system()
    niter = int(config.get('niter', 20))
    solver._ensure_fly_neighbors()

    dp_in = np.zeros_like(solver.ps_pred)
    dp_out = np.zeros_like(dp_in)

    capture_iters = _iter_logging_enabled(config)
    if capture_iters:
        _ensure_log_base(solver)
        _log_iteration_state(solver, config, dp_in, 0.0, is_displacement=True)

    for it in range(niter):
        max_step = _update_jacobi_lin(solver, dp_in, dp_out, rhs, diag, fixed_mask)
        dp_in, dp_out = dp_out, dp_in
        if VERBOSITY >= 2:
            print(f"[CPU][JacobiDiff iter {it + 1}/{niter}] max |Δdp| = {max_step:.3e}")
        if capture_iters:
            _log_iteration_state(solver, config, dp_in, max_step, is_displacement=True)

    solver.ps_cor[:] = solver.ps_pred + dp_in
    if solver.fixed.size > 0:
        solver.ps_cor[solver.fixed] = solver.x[solver.fixed]


def _update_jacobi_fly(solver: TrussSolver, ps_in: np.ndarray, ps_out: np.ndarray) -> float:
    """Single Jacobi-fly iteration. Matches C++ updateJacobi_fly."""
    neighs = solver._neigh_indices
    kngs = solver._neigh_k
    l0ngs = solver._neigh_l0
    max_neighs = solver._max_neighs or 0
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
            sum_j += dij * (k * l0_ij / safe)
            sum_j += pj * k
            Aii += k
        new_pi = sum_j / max(Aii, 1e-12)
        step = np.linalg.norm(new_pi - pi)
        if step > max_step:
            max_step = step
        ps_out[i] = new_pi
    
    if solver.fixed.size > 0:
        ps_out[solver.fixed] = solver.x[solver.fixed]
    return max_step


def solve_jacobi_fly(solver: TrussSolver, config: dict) -> None:
    """Jacobi-Fly: Copy of TrussDynamics_d::updateJacobi_fly() and OpenCL jacobi_fly kernel."""
    niter = int(config.get('niter', 10))
    solver._ensure_fly_neighbors()

    ps_in = solver.ps_pred.copy()
    ps_out = ps_in.copy()
    
    for it in range(niter):
        max_step = _update_jacobi_fly(solver, ps_in, ps_out)
        ps_in, ps_out = ps_out, ps_in
        if VERBOSITY >= 2:
            print(f"[CPU][JacobiFly iter {it + 1}/{niter}] max |Δp| = {max_step:.3e}")
    
    solver.ps_cor[:] = ps_in


def _update_gs_lin(solver: TrussSolver, dp: np.ndarray, 
                   rhs: np.ndarray, diag: np.ndarray, fixed_mask: np.ndarray) -> float:
    """Single Gauss-Seidel iteration on linearized system. Matches C++ updateGaussSeidel_lin."""
    neighs = solver._neigh_indices
    kngs = solver._neigh_k
    max_neighs = solver._max_neighs or 0
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
    return max_step


def solve_gs_diff(solver: TrussSolver, config: dict) -> None:
    """GS-Diff: precomputed linear RHS, in-place Gauss-Seidel on displacements."""
    rhs, diag, fixed_mask = solver._build_linear_diff_system()
    niter = int(config.get('niter', 20))
    solver._ensure_fly_neighbors()

    dp = np.zeros_like(solver.ps_pred)

    capture_iters = _iter_logging_enabled(config)
    if capture_iters:
        _ensure_log_base(solver)
        _log_iteration_state(solver, config, dp, 0.0, is_displacement=True)

    for it in range(niter):
        max_step = _update_gs_lin(solver, dp, rhs, diag, fixed_mask)
        if VERBOSITY >= 2:
            print(f"[CPU][GSDiff iter {it + 1}/{niter}] max |Δdp| = {max_step:.3e}")
        if capture_iters:
            _log_iteration_state(solver, config, dp, max_step, is_displacement=True)

    solver.ps_cor[:] = solver.ps_pred + dp
    if solver.fixed.size > 0:
        solver.ps_cor[solver.fixed] = solver.x[solver.fixed]


def _update_gs_fly(solver: TrussSolver, ps: np.ndarray) -> float:
    """Single Gauss-Seidel-fly iteration. Matches C++ updateGaussSeidel_fly."""
    neighs = solver._neigh_indices
    kngs = solver._neigh_k
    l0ngs = solver._neigh_l0
    max_neighs = solver._max_neighs or 0
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
    return max_step


def solve_gs_fly(solver: TrussSolver, config: dict) -> None:
    """GS-Fly: Gauss-Seidel with on-the-fly force computation."""
    niter = int(config.get('niter', 10))
    solver._ensure_fly_neighbors()

    ps = solver.ps_pred.copy()
    capture_iters = _iter_logging_enabled(config)
    if capture_iters:
        _ensure_log_base(solver)
        _log_iteration_state(solver, config, ps, 0.0, is_displacement=False)
    
    for it in range(niter):
        max_step = _update_gs_fly(solver, ps)
        if VERBOSITY >= 2:
            print(f"[CPU][GSFly iter {it + 1}/{niter}] max |Δp| = {max_step:.3e}")
        if capture_iters:
            _log_iteration_state(solver, config, ps, max_step, is_displacement=False)
    
    solver.ps_cor[:] = ps
    if solver.fixed.size > 0:
        solver.ps_cor[solver.fixed] = solver.x[solver.fixed]


# Solver registry
SOLVERS = {
    'vbd': solve_vbd,
    'vbd_serial': solve_vbd,
    'direct': solve_direct,
    'cholesky': solve_cholesky,
    'jacobi_diff': solve_jacobi_diff,
    'jacobi_fly': solve_jacobi_fly,
    'gs_diff': solve_gs_diff,
    'gs_fly': solve_gs_fly,
    'momentum': solve_iterative_momentum,
    'chebyshev': solve_iterative_chebyshev,
}

def get_solver(name: str) -> Callable:
    """Retrieve solver callback by name."""
    print("truss_solver.get_solver()")
    if name not in SOLVERS:
        raise ValueError(f"Unknown solver '{name}'; available: {list(SOLVERS.keys())}")
    return SOLVERS[name]


def select_spectral_radius_method(solver_name: str, sub_method: str) -> Optional[str]:
    """Return compatible spectral radius estimation scheme for the solver."""
    name = solver_name.lower()
    if name in {'jacobi_diff', 'jacobi'}:
        return 'jacobi_diff'
    if name in {'gs_diff', 'gauss_seidel'}:
        return 'gs_diff'
    if name in {'momentum', 'chebyshev'}:
        sub = str(sub_method).lower()
        if sub in {'jacobi_diff', 'gs_diff'}:
            return sub
        return None
    return None


def _execute_solver_method(
    solver: 'TrussSolver',
    *,
    solver_name: str,
    label: str,
    base_positions: np.ndarray,
    displacement: np.ndarray,
    config: Dict[str, Any],
    estimate_rho: bool,
    spectral_method: Optional[str],
) -> Dict[str, Any]:
    """Run a single solver method on the provided solver instance."""

    solver.solver_callback = get_solver(solver_name)
    solver.solver_config = dict(config)
    solver.verbose = int(solver.solver_config.get('verbose', solver.verbose))

    base_copy = base_positions.astype(np.float64, copy=True)
    initial = base_copy + displacement
    zero_vel = np.zeros_like(base_copy)

    solver.reset_state(positions=base_copy, velocities=zero_vel)
    solver.truss.points[:] = base_copy

    solver.ps_pred[:] = initial
    solver.ps_cor[:] = initial
    if solver.fixed.size > 0:
        solver.ps_pred[solver.fixed] = base_copy[solver.fixed]
        solver.ps_cor[solver.fixed] = base_copy[solver.fixed]

    state = solver.solver_state
    state['iter_logs'] = {}
    state['log_base'] = solver.ps_pred.copy()
    state['diagnostics'] = {}

    diagnostics = state['diagnostics']
    if estimate_rho:
        method = spectral_method or select_spectral_radius_method(solver_name, solver.solver_config.get('sub_method', ''))
        if method is not None:
            try:
                rho = estimate_iterative_spectral_radius(solver, method=method)
                diagnostics['spectral_method'] = method
                diagnostics['spectral_radius'] = rho
            except Exception as err:
                diagnostics['spectral_method'] = method
                diagnostics['spectral_radius_error'] = str(err)
        else:
            diagnostics['spectral_radius_error'] = 'unsupported solver for estimation'

    solver.solver_callback(solver, solver.solver_config)

    logs = state.get('iter_logs', {}).get(label)
    if logs and logs.get('positions'):
        positions = np.stack(logs['positions'], axis=0)
        max_entries = logs.get('max_step', [])
        max_step = np.array(max_entries, dtype=np.float64)
        if max_step.size < positions.shape[0]:
            pad = positions.shape[0] - max_step.size
            max_step = np.pad(max_step, (0, pad))
        deltas = np.linalg.norm(positions[1:] - positions[:-1], axis=2)
        residual = deltas.max(axis=1) if deltas.size else np.array([], dtype=np.float64)
    else:
        final_pos = solver.ps_cor.copy()
        positions = np.stack([initial, final_pos], axis=0)
        diff = final_pos - initial
        max_delta = float(np.max(np.abs(diff))) if diff.size else 0.0
        max_step = np.array([0.0, max_delta], dtype=np.float64)
        residual = np.array([max_delta], dtype=np.float64) if diff.size else np.array([], dtype=np.float64)

    return {
        'label': label,
        'solver': solver_name,
        'positions': positions,
        'max_step': max_step,
        'residual': residual,
        'base': base_copy,
        'initial': initial.copy(),
        'final': positions[-1],
        'diagnostics': diagnostics.copy(),
    }


def run_solver_suite(
    truss_factory: Callable[[], Any],
    dt: float,
    displacement: np.ndarray,
    solver_requests: List[Any],
    *,
    base_positions: Optional[np.ndarray] = None,
    fixed_points: Optional[Sequence[int]] = None,
    gravity: Optional[np.ndarray] = None,
    estimate_rho: bool = False,
    existing_solver: Optional['TrussSolver'] = None,
) -> Dict[str, Dict[str, Any]]:
    """Run identical initial state across multiple solvers and collect traces."""

    if existing_solver is not None:
        template_truss = existing_solver.truss
    else:
        template_truss = truss_factory()
    template_points = template_truss.points.astype(np.float64, copy=True)
    base = template_points if base_positions is None else np.asarray(base_positions, dtype=np.float64)
    if base.shape != template_points.shape:
        raise ValueError("base_positions shape mismatch with template truss")

    disp_array = np.asarray(displacement, dtype=np.float64)
    if disp_array.shape != base.shape:
        raise ValueError("displacement shape mismatch")

    if fixed_points is None:
        fixed_seq = sorted(template_truss.fixed)
    else:
        fixed_seq = list(int(idx) for idx in fixed_points)

    gravity_vec = np.zeros(3, dtype=np.float64) if gravity is None else np.asarray(gravity, dtype=np.float64)
    if gravity_vec.shape != (3,):
        raise ValueError(f"gravity vector must be length 3, got {gravity_vec.shape}")

    prepared: List[Dict[str, Any]] = []
    for idx, spec in enumerate(solver_requests):
        if isinstance(spec, dict):
            solver_name = spec['solver']
            config = dict(spec.get('config', {}))
            label = spec.get('label', solver_name)
            request_estimate = bool(spec.get('estimate_rho', estimate_rho))
            spectral_method = spec.get('spectral_method')
        elif isinstance(spec, tuple):
            if len(spec) == 3:
                solver_name, config, label = spec
            elif len(spec) == 2:
                solver_name, config = spec
                label = f"solver_{idx}"
            else:
                raise ValueError("solver request tuple must have 2 or 3 entries")
            config = dict(config)
            request_estimate = estimate_rho
            spectral_method = None
        else:
            raise TypeError("solver_requests entries must be dict or tuple")

        if not label:
            label = f"solver_{idx}"

        config.setdefault('capture_iters', True)
        config.setdefault('log_label', label)

        prepared.append({
            'solver': solver_name,
            'config': config,
            'label': label,
            'estimate_rho': request_estimate,
            'spectral_method': spectral_method,
        })

    if not prepared:
        return {}

    if existing_solver is not None:
        shared_solver = existing_solver
    else:
        shared_solver = TrussSolver(
            truss_factory(),
            dt=dt,
            gravity=gravity_vec,
            solver=get_solver(prepared[0]['solver']),
            solver_config=dict(prepared[0]['config']),
            fixed_points=fixed_seq,
            verbose=int(prepared[0]['config'].get('verbose', 0)),
        )

    results: Dict[str, Dict[str, Any]] = {}
    for entry in prepared:
        result = _execute_solver_method(
            shared_solver,
            solver_name=entry['solver'],
            label=entry['label'],
            base_positions=base,
            displacement=disp_array,
            config=entry['config'],
            estimate_rho=entry['estimate_rho'],
            spectral_method=entry['spectral_method'],
        )
        results[entry['label']] = result

    return results


def iter_convergence_rows(results: Dict[str, Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Convert solver results to per-iteration convergence rows."""
    rows: List[Dict[str, Any]] = []
    for data in results.values():
        solver_name = data['solver']
        label = data.get('label', solver_name)
        max_step = np.asarray(data.get('max_step', []), dtype=np.float64)
        residual = np.asarray(data.get('residual', []), dtype=np.float64)
        for itr, step_val in enumerate(max_step):
            if itr == 0:
                resid_val = float('nan')
            else:
                idx = itr - 1
                resid_val = residual[idx] if idx < residual.size else float('nan')
            rows.append({
                'solver': solver_name,
                'label': label,
                'iteration': itr,
                'max_step': float(step_val),
                'residual': float(resid_val),
            })
    return rows


def write_convergence_csv(results: Dict[str, Dict[str, Any]], csv_path: str) -> None:
    """Persist convergence rows into a CSV file."""
    fieldnames = ['solver', 'label', 'iteration', 'max_step', 'residual']
    rows = iter_convergence_rows(results)
    with open(csv_path, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            output = dict(row)
            if math.isnan(output['residual']):
                output['residual'] = ''
            writer.writerow(output)
