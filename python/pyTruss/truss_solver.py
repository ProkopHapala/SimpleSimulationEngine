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


def _build_csr(n_points: int, edges_i: np.ndarray, edges_j: np.ndarray,
               weights: np.ndarray, rest_lengths: Optional[np.ndarray] = None
               ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Optional[np.ndarray]]:
    deg = np.zeros(n_points, dtype=np.int32)
    if edges_i.size > 0:
        np.add.at(deg, edges_i, 1)
        np.add.at(deg, edges_j, 1)

    row_ptr = np.empty(n_points + 1, dtype=np.int32)
    row_ptr[0] = 0
    np.cumsum(deg, out=row_ptr[1:])
    total = int(row_ptr[-1])

    col_idx = np.empty(total, dtype=np.int32)
    data = np.empty(total, dtype=np.float64)
    rest_data = np.empty(total, dtype=np.float64) if rest_lengths is not None else None

    if edges_i.size > 0:
        cursor = row_ptr[:-1].copy()
        for e in range(edges_i.shape[0]):
            i = int(edges_i[e])
            j = int(edges_j[e])
            w = float(weights[e])
            pos_i = cursor[i]
            col_idx[pos_i] = j
            data[pos_i] = w
            if rest_data is not None:
                rest_data[pos_i] = float(rest_lengths[e])
            cursor[i] += 1
            pos_j = cursor[j]
            col_idx[pos_j] = i
            data[pos_j] = w
            if rest_data is not None:
                rest_data[pos_j] = float(rest_lengths[e])
            cursor[j] += 1

    return row_ptr, col_idx, data, rest_data


def _prepare_fixed_data(solver: "TrussSolver", x_ref: np.ndarray, fixed: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    fixed_mask = np.zeros(solver.n_points, dtype=bool)
    if fixed.size == 0:
        return fixed_mask, np.zeros((0, 3), dtype=np.float64), np.zeros((0, 3), dtype=np.float64)
    fixed_mask[fixed] = True
    if solver.fixed_positions is not None and solver.fixed_positions.shape[0] == fixed.size:
        fixed_pos = solver.fixed_positions.copy()
    else:
        fixed_pos = x_ref[fixed].copy()
    fixed_dp = fixed_pos - x_ref[fixed]
    return fixed_mask, fixed_pos, fixed_dp


def _apply_fixed(x: np.ndarray, fixed_idx: np.ndarray, fixed_vals: np.ndarray) -> None:
    if fixed_idx.size == 0:
        return
    x[fixed_idx] = fixed_vals


def _compute_diff_rhs(solver: "TrussSolver", y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    bonds = solver.bonds
    ks_edges = solver.ks
    rest_vec = solver.rest_vectors
    diag = solver.mass_diag.copy()
    b = solver.mass_diag[:, None] * y
    if bonds.size > 0:
        i_idx = bonds[:, 0]
        j_idx = bonds[:, 1]
        np.add.at(diag, i_idx, ks_edges)
        np.add.at(diag, j_idx, ks_edges)
        np.add.at(b, i_idx, ks_edges[:, None] * rest_vec)
        np.add.at(b, j_idx, -ks_edges[:, None] * rest_vec)
    return b, diag


def _get_initial_guess(state: dict, key: str, reference: np.ndarray) -> np.ndarray:
    guess = state.get(key)
    if guess is None or guess.shape != reference.shape:
        return reference.copy()
    return guess.copy()


def _jacobi_diff(solver: "TrussSolver", dp_init: np.ndarray, b: np.ndarray,
                 diag: np.ndarray, fixed_idx: np.ndarray, fixed_dp: np.ndarray,
                 niter: int, tol: Optional[float]) -> np.ndarray:
    inv_diag = 1.0 / np.where(diag > 1e-12, diag, 1e-12)
    bonds = solver.bonds
    ks_edges = solver.ks

    dp_prev = dp_init.copy()
    dp_next = dp_prev.copy()
    sum_term = np.zeros_like(dp_prev)

    for _ in range(max(1, niter)):
        sum_term.fill(0.0)
        if bonds.size > 0:
            i_idx = bonds[:, 0]
            j_idx = bonds[:, 1]
            np.add.at(sum_term, i_idx, ks_edges[:, None] * dp_prev[j_idx])
            np.add.at(sum_term, j_idx, ks_edges[:, None] * dp_prev[i_idx])
        dp_next[:] = (b + sum_term) * inv_diag[:, None]
        _apply_fixed(dp_next, fixed_idx, fixed_dp)
        if tol is not None and np.linalg.norm(dp_next - dp_prev, ord=np.inf) < tol:
            dp_prev = dp_next.copy()
            break
        dp_prev, dp_next = dp_next, dp_prev

    _apply_fixed(dp_prev, fixed_idx, fixed_dp)
    return dp_prev


def _gs_diff(solver: "TrussSolver", dp_init: np.ndarray, b: np.ndarray,
             diag: np.ndarray, fixed_mask: np.ndarray, fixed_idx: np.ndarray,
             fixed_dp: np.ndarray, niter: int, tol: Optional[float]) -> np.ndarray:
    bonds = solver.bonds
    ks_edges = solver.ks
    row_ptr, col_idx, weights, _ = _build_csr(solver.n_points, bonds[:, 0] if bonds.size > 0 else np.zeros(0, dtype=np.int32),
                                              bonds[:, 1] if bonds.size > 0 else np.zeros(0, dtype=np.int32), ks_edges)
    inv_diag = 1.0 / np.where(diag > 1e-12, diag, 1e-12)

    dp_curr = dp_init.copy()

    for _ in range(max(1, niter)):
        max_delta = 0.0
        for vid in range(solver.n_points):
            if fixed_mask[vid]:
                continue
            start = row_ptr[vid]
            end = row_ptr[vid + 1]
            if start == end:
                sum_term = np.zeros(3, dtype=np.float64)
            else:
                neigh = col_idx[start:end]
                w = weights[start:end]
                sum_term = (w[:, None] * dp_curr[neigh]).sum(axis=0)
            new_val = (b[vid] + sum_term) * inv_diag[vid]
            delta = np.linalg.norm(new_val - dp_curr[vid], ord=np.inf)
            if delta > max_delta:
                max_delta = delta
            dp_curr[vid] = new_val
        _apply_fixed(dp_curr, fixed_idx, fixed_dp)
        if tol is not None and max_delta < tol:
            break

    return dp_curr


def _jacobi_fly(solver: "TrussSolver", dp_init: np.ndarray, y: np.ndarray, diag: np.ndarray,
                fixed_idx: np.ndarray, fixed_dp: np.ndarray,
                niter: int, tol: Optional[float]) -> np.ndarray:
    inv_diag = 1.0 / np.where(diag > 1e-12, diag, 1e-12)
    mass_term = solver.mass_diag[:, None] * y
    bonds = solver.bonds
    ks_edges = solver.ks
    rest_len_edge = np.linalg.norm(solver.rest_vectors, axis=1) if solver.rest_vectors.size > 0 else np.zeros(0, dtype=np.float64)

    dp_prev = dp_init.copy()
    dp_next = dp_prev.copy()

    for _ in range(max(1, niter)):
        bi = mass_term.copy()
        if bonds.size > 0:
            i_idx = bonds[:, 0]
            j_idx = bonds[:, 1]
            diff = dp_prev[i_idx] - dp_prev[j_idx]
            lengths = np.linalg.norm(diff, axis=1)
            safe = np.maximum(lengths, 1e-12)
            coeff = ks_edges * (rest_len_edge / safe - 1.0)
            contrib_i = coeff[:, None] * diff + ks_edges[:, None] * dp_prev[j_idx]
            contrib_j = -coeff[:, None] * diff + ks_edges[:, None] * dp_prev[i_idx]
            np.add.at(bi, i_idx, contrib_i)
            np.add.at(bi, j_idx, contrib_j)
        dp_next[:] = bi * inv_diag[:, None]
        _apply_fixed(dp_next, fixed_idx, fixed_dp)
        if tol is not None and np.linalg.norm(dp_next - dp_prev, ord=np.inf) < tol:
            dp_prev = dp_next.copy()
            break
        dp_prev, dp_next = dp_next, dp_prev

    _apply_fixed(dp_prev, fixed_idx, fixed_dp)
    return dp_prev


def _gs_fly(solver: "TrussSolver", dp_init: np.ndarray, y: np.ndarray,
            fixed_mask: np.ndarray, fixed_idx: np.ndarray, fixed_dp: np.ndarray,
            niter: int, tol: Optional[float]) -> np.ndarray:
    bonds = solver.bonds
    ks_edges = solver.ks
    rest_len_edge = np.linalg.norm(solver.rest_vectors, axis=1) if solver.rest_vectors.size > 0 else np.zeros(0, dtype=np.float64)
    row_ptr, col_idx, weights, rest_l0 = _build_csr(solver.n_points,
                                                    bonds[:, 0] if bonds.size > 0 else np.zeros(0, dtype=np.int32),
                                                    bonds[:, 1] if bonds.size > 0 else np.zeros(0, dtype=np.int32),
                                                    ks_edges, rest_lengths=rest_len_edge)
    mass_diag = solver.mass_diag

    dp_curr = dp_init.copy()

    for _ in range(max(1, niter)):
        max_delta = 0.0
        for vid in range(solver.n_points):
            if fixed_mask[vid]:
                continue
            bi = mass_diag[vid] * y[vid]
            Aii = mass_diag[vid]
            start = row_ptr[vid]
            end = row_ptr[vid + 1]
            if start != end:
                neigh = col_idx[start:end]
                w = weights[start:end]
                l0 = rest_l0[start:end]
                diff = dp_curr[vid] - dp_curr[neigh]
                lengths = np.linalg.norm(diff, axis=1)
                safe = np.maximum(lengths, 1e-12)
                coeff = w * (l0 / safe - 1.0)
                bi += (coeff[:, None] * diff).sum(axis=0)
                bi += (w[:, None] * dp_curr[neigh]).sum(axis=0)
                Aii += w.sum()
            new_val = bi / max(Aii, 1e-12)
            delta = np.linalg.norm(new_val - dp_curr[vid], ord=np.inf)
            if delta > max_delta:
                max_delta = delta
            dp_curr[vid] = new_val
        _apply_fixed(dp_curr, fixed_idx, fixed_dp)
        if tol is not None and max_delta < tol:
            break

    return dp_curr


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

        # Predictor/corrector buffers for per-step solve
        self.ps_pred = np.zeros_like(self.x)
        self.ps_cor = np.zeros_like(self.x)
        self.forces = np.zeros_like(self.x)

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

    bonds = solver.bonds
    masses = solver.masses
    ks = solver.ks
    rest_vec = solver.rest_vectors
    l0s = np.linalg.norm(rest_vec, axis=1) if rest_vec.size > 0 else np.zeros(0, dtype=np.float64)

    n_points = solver.n_points
    dt_sq = solver.dt * solver.dt
    inv_h2 = 1.0 / dt_sq
    I3 = np.eye(3, dtype=np.float64)

    fixed_mask = np.zeros(n_points, dtype=bool)
    if solver.fixed.size > 0:
        fixed_mask[solver.fixed] = True

    x_curr = solver.x.copy()
    y = solver.ps_pred.copy()

    eye_mass = np.einsum('i,ab->iab', masses * inv_h2, I3)

    for itr in range(niter):
        if solver.fixed.size > 0:
            x_curr[fixed_mask] = solver.x[fixed_mask]

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
        if solver.fixed.size > 0:
            H[fixed_mask] = eye_mass[fixed_mask]

        if VERBOSITY >= 2 and n_points > 1:
            ivert = 1
            gi = grad[ivert]
            Hi = H[ivert]
            print(f"[CPU][VBD iter {itr}] grad[{ivert}] = {gi}")
            print(f"[CPU][VBD iter {itr}] H[{ivert}] =\n{Hi}")

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


def solve_jacobi_diff_cpu(solver: TrussSolver, config: dict) -> None:
    print("truss_solver.solve_jacobi_diff_cpu()")
    state = solver.solver_state
    fixed_mask, fixed_pos, fixed_dp = _prepare_fixed_data(solver, solver.x, solver.fixed)
    b, diag = _compute_diff_rhs(solver, solver.ps_pred)
    dp0 = _get_initial_guess(state, 'last_dp', np.zeros_like(solver.x))
    _apply_fixed(dp0, solver.fixed, fixed_dp)
    niter = int(config.get('niter', 20))
    tol = config.get('tol')
    tol = float(tol) if tol is not None else None
    dp = _jacobi_diff(solver, dp0, b, diag, solver.fixed, fixed_dp, niter, tol)
    _apply_fixed(dp, solver.fixed, fixed_dp)
    solver.ps_cor[:] = solver.ps_pred + dp
    if solver.fixed.size > 0:
        solver.ps_cor[solver.fixed] = fixed_pos
    state['last_dp'] = dp.copy()


def solve_gs_diff_cpu(solver: TrussSolver, config: dict) -> None:
    print("truss_solver.solve_gs_diff_cpu()")
    state = solver.solver_state
    fixed_mask, fixed_pos, fixed_dp = _prepare_fixed_data(solver, solver.x, solver.fixed)
    b, diag = _compute_diff_rhs(solver, solver.ps_pred)
    dp0 = _get_initial_guess(state, 'last_dp', np.zeros_like(solver.x))
    _apply_fixed(dp0, solver.fixed, fixed_dp)
    niter = int(config.get('niter', 20))
    tol = config.get('tol')
    tol = float(tol) if tol is not None else None
    dp = _gs_diff(solver, dp0, b, diag, fixed_mask, solver.fixed, fixed_dp, niter, tol)
    _apply_fixed(dp, solver.fixed, fixed_dp)
    solver.ps_cor[:] = solver.ps_pred + dp
    if solver.fixed.size > 0:
        solver.ps_cor[solver.fixed] = fixed_pos
    state['last_dp'] = dp.copy()


def solve_jacobi_fly_cpu(solver: TrussSolver, config: dict) -> None:
    print("truss_solver.solve_jacobi_fly_cpu()")
    state = solver.solver_state
    fixed_mask, fixed_pos, fixed_dp = _prepare_fixed_data(solver, solver.x, solver.fixed)
    _, diag = _compute_diff_rhs(solver, solver.ps_pred)
    dp0 = _get_initial_guess(state, 'last_dp', np.zeros_like(solver.x))
    _apply_fixed(dp0, solver.fixed, fixed_dp)
    niter = int(config.get('niter', 10))
    tol = config.get('tol')
    tol = float(tol) if tol is not None else None
    dp = _jacobi_fly(solver, dp0, solver.ps_pred, diag, solver.fixed, fixed_dp, niter, tol)
    _apply_fixed(dp, solver.fixed, fixed_dp)
    solver.ps_cor[:] = solver.ps_pred + dp
    if solver.fixed.size > 0:
        solver.ps_cor[solver.fixed] = fixed_pos
    state['last_dp'] = dp.copy()


def solve_gs_fly_cpu(solver: TrussSolver, config: dict) -> None:
    print("truss_solver.solve_gs_fly_cpu()")
    state = solver.solver_state
    fixed_mask, fixed_pos, fixed_dp = _prepare_fixed_data(solver, solver.x, solver.fixed)
    dp0 = _get_initial_guess(state, 'last_dp', np.zeros_like(solver.x))
    _apply_fixed(dp0, solver.fixed, fixed_dp)
    niter = int(config.get('niter', 10))
    tol = config.get('tol')
    tol = float(tol) if tol is not None else None
    dp = _gs_fly(solver, dp0, solver.ps_pred, fixed_mask, solver.fixed, fixed_dp, niter, tol)
    _apply_fixed(dp, solver.fixed, fixed_dp)
    solver.ps_cor[:] = solver.ps_pred + dp
    if solver.fixed.size > 0:
        solver.ps_cor[solver.fixed] = fixed_pos
    state['last_dp'] = dp.copy()


# Solver registry
SOLVERS = {
    'vbd': solve_vbd,
    'vbd_serial': solve_vbd,
    'direct': solve_direct,
    'cholesky': solve_cholesky,
    'momentum': solve_iterative_momentum,
    'jacobi_diff_cpu': solve_jacobi_diff_cpu,
    'jacobi_diff': solve_jacobi_diff_cpu,
    'gs_diff_cpu': solve_gs_diff_cpu,
    'gs_diff': solve_gs_diff_cpu,
    'jacobi_fly_cpu': solve_jacobi_fly_cpu,
    'jacobi_fly': solve_jacobi_fly_cpu,
    'gs_fly_cpu': solve_gs_fly_cpu,
    'gs_fly': solve_gs_fly_cpu,
}


def get_solver(name: str) -> Callable:
    """Retrieve solver callback by name."""
    print("truss_solver.get_solver()")
    if name not in SOLVERS:
        raise ValueError(f"Unknown solver '{name}'; available: {list(SOLVERS.keys())}")
    return SOLVERS[name]
