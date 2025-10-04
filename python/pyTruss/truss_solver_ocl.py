import os
from typing import Any, Callable, Dict, List, Optional

import numpy as np
import pyopencl as cl

from sparse import build_neighbor_list, neigh_stiffness, neighs_to_dense_arrays
from truss import Truss


VERBOSITY = 0


def set_verbosity(level: int) -> None:
    global VERBOSITY
    VERBOSITY = int(level)


def _print_state(tag: str, positions: np.ndarray, velocities: np.ndarray) -> None:
    print(f"{tag} positions:\n{positions}")
    print(f"{tag} velocities:\n{velocities}")


def build_rest_length_dense(neigh_lists: List[List[int]], bonds: np.ndarray, l0s: np.ndarray, n_max: int) -> np.ndarray:
    """Dense rest-length table matching the padded neighbor layout."""
    print("truss_solver_ocl.build_rest_length_dense()")
    assert len(bonds) == len(l0s), "bonds and l0s must have identical length"
    pair_rest: Dict[tuple[int, int], float] = {}
    for (ia, ib), rest in zip(bonds, l0s):
        ia = int(ia)
        ib = int(ib)
        val = float(rest)
        pair_rest[(ia, ib)] = val
        pair_rest[(ib, ia)] = val
    rest_dense = np.zeros((len(neigh_lists), n_max), dtype=np.float32)
    for vid, neigh in enumerate(neigh_lists):
        for jj, nj in enumerate(neigh):
            rest_dense[vid, jj] = pair_rest[(vid, int(nj))]
    return rest_dense


class OCLRuntime:
    def __init__(self) -> None:
        print("truss_solver_ocl.OCLRuntime.__init__()")
        platforms = cl.get_platforms()
        if not platforms:
            raise RuntimeError("No OpenCL platforms are available")
        self.ctx = cl.Context(
            dev_type=cl.device_type.ALL,
            properties=[(cl.context_properties.PLATFORM, platforms[0])],
        )
        self.queue = cl.CommandQueue(self.ctx)
        kernel_path = os.path.join(os.path.dirname(__file__), "truss.cl")
        with open(kernel_path, "r", encoding="utf-8") as f:
            kernel_src = f.read()
        self.prg = cl.Program(self.ctx, kernel_src).build()
        self.mf = cl.mem_flags


class TrussSolverOCL:
    """Projective Dynamics time-stepper driven by OpenCL kernels."""

    def __init__(
        self,
        truss: Truss,
        *,
        dt: float,
        gravity: np.ndarray,
        solver: Callable[[Truss, Any], np.ndarray],
        solver_config: Optional[Dict[str, Any]] = None,
        fixed_points: Optional[List[int]] = None,
        track_indices: Optional[List[int]] = None,
        verbose: int = 0,
    ) -> None:
        print("truss_solver_ocl.TrussSolverOCL.__init__()")
        self.truss = truss
        self.dt = float(dt)
        self.gravity = np.asarray(gravity, dtype=np.float64)
        self.runtime = OCLRuntime()
        self.solver_callback = solver
        self.solver_config = dict(solver_config) if solver_config is not None else {}
        self.verbose = int(verbose)

        self.x = truss.points.astype(np.float64, copy=True)
        self.v = np.zeros_like(self.x)

        if fixed_points is None:
            fixed_seq = sorted(truss.fixed)
        else:
            fixed_seq = sorted(fixed_points)
        self.fixed = np.asarray(fixed_seq, dtype=int) if len(fixed_seq) > 0 else np.zeros(0, dtype=int)

        if track_indices is not None and len(track_indices) > 0:
            self.track_idx = np.asarray(track_indices, dtype=int)
            self.trajectory: Optional[List[np.ndarray]] = [self.x[self.track_idx].copy()]
        else:
            self.track_idx = None
            self.trajectory = None

        self.bonds = np.asarray(truss.bonds, dtype=np.int32)
        if self.bonds.size == 0:
            self.bonds = self.bonds.reshape(0, 2)
        self.ks = truss.ks.astype(np.float64, copy=False)
        if self.ks.ndim == 0:
            self.ks = np.full(self.bonds.shape[0], float(self.ks))
        self.masses = truss.masses.astype(np.float64, copy=False)
        self.rest_lengths = truss.get_rest_lengths().astype(np.float64, copy=False)

        self.n_points = self.x.shape[0]
        if self.n_points == 0:
            raise ValueError("Truss must contain at least one vertex for GPU solver")

        inv_dt2 = 1.0 / (self.dt * self.dt)
        self.mass_diag = np.where(self.masses > 0.0, self.masses * inv_dt2, 1e-12)

        self.fixed_positions = self.x[self.fixed].copy() if self.fixed.size > 0 else None

        self.solver_state: Dict[str, Any] = {"owner": self}

        self._prepare_static_pd_data()

    def _prepare_static_pd_data(self) -> None:
        print("truss_solver_ocl.TrussSolverOCL._prepare_static_pd_data()")
        if self.bonds.size == 0:
            raise ValueError("GPU PD solver requires at least one bond (no springs found)")
        neighbs = build_neighbor_list(self.bonds, self.n_points)
        neighs, kngs, n_max = neigh_stiffness(neighbs, self.bonds, self.ks)
        padded_neighs, padded_kngs, _ = neighs_to_dense_arrays(neighs, kngs, n_max)
        self.neigh_lists = neighs
        self.n_max = n_max
        self.neighs_dense = padded_neighs.astype(np.int32, copy=True)
        self.kngs_dense = padded_kngs.astype(np.float32, copy=True)
        self.rest_dense = build_rest_length_dense(neighs, self.bonds, self.rest_lengths, n_max).astype(np.float32, copy=False)

    def _create_vbd_serial_workspace(self) -> Dict[str, Any]:
        print("truss_solver_ocl.TrussSolverOCL._create_vbd_serial_workspace()")
        runtime = self.runtime
        mf = runtime.mf
        x_host = np.zeros((self.n_points, 4), dtype=np.float32)
        y_host = np.zeros_like(x_host)
        workspace = {
            "x_host": x_host,
            "y_host": y_host,
            "x_buf": cl.Buffer(runtime.ctx, mf.READ_WRITE, x_host.nbytes),
            "y_buf": cl.Buffer(runtime.ctx, mf.READ_WRITE, y_host.nbytes),
            "neighs_buf": cl.Buffer(runtime.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.neighs_dense),
            "kngs_buf": cl.Buffer(runtime.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.kngs_dense),
            "l0_buf": cl.Buffer(runtime.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.rest_dense),
        }
        return workspace

    def _get_vbd_serial_workspace(self) -> Dict[str, Any]:
        workspace = self.solver_state.get("vbd_serial_workspace")
        if workspace is None:
            workspace = self._create_vbd_serial_workspace()
            self.solver_state["vbd_serial_workspace"] = workspace
        return workspace

    def reset_state(
        self,
        *,
        positions: Optional[np.ndarray] = None,
        velocities: Optional[np.ndarray] = None,
    ) -> None:
        print("truss_solver_ocl.TrussSolverOCL.reset_state()")
        if positions is not None:
            self.x = positions.astype(np.float64, copy=True)
        if velocities is not None:
            self.v = velocities.astype(np.float64, copy=True)
        else:
            self.v = np.zeros_like(self.x)
        if self.fixed.size > 0:
            self.fixed_positions = self.x[self.fixed].copy()
        if self.track_idx is not None and self.trajectory is not None:
            self.trajectory = [self.x[self.track_idx].copy()]
        self.solver_state.clear()
        self.solver_state["owner"] = self

    def run(self, nsteps: int) -> tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]:
        print("truss_solver_ocl.TrussSolverOCL.run()")
        if nsteps <= 0:
            raise ValueError("nsteps must be positive")

        dt = self.dt
        gravity = self.gravity
        solver = self.solver_callback
        config = self.solver_config

        x = self.x
        v = self.v
        fixed = self.fixed
        trajectory = self.trajectory
        track_idx = self.track_idx

        for step in range(nsteps):
            x_new = solver(
                self.truss,
                state=self.solver_state,
                dt=dt,
                gravity=gravity,
                x=x,
                v=v,
                fixed=fixed,
                config=config,
            )
            if x_new.shape != x.shape:
                raise ValueError("Solver callback returned array with unexpected shape")
            v = (x_new - x) / dt
            if fixed.size > 0:
                x_new[fixed] = x[fixed]
                v[fixed] = 0.0
            x = x_new
            if trajectory is not None and track_idx is not None:
                trajectory.append(x[track_idx].copy())
            if self.verbose and (step == 0 or (step + 1) % max(1, nsteps // 10) == 0):
                print(f"  [GPU] Step {step + 1}/{nsteps}, t={dt * (step + 1):.3f} s")
            if VERBOSITY >= 1:
                _print_state(f"[GPU][step {step + 1}/{nsteps}]", x, v)

        self.x = x
        self.v = v
        self.truss.points = x.copy()
        if self.fixed.size > 0 and self.fixed_positions is not None:
            self.fixed_positions = x[self.fixed].copy()

        traj_out = None
        if trajectory is not None:
            traj_out = np.stack(trajectory, axis=0)

        return x, v, traj_out

    # Workspace accessors -------------------------------------------------

    def get_vbd_serial_workspace(self) -> Dict[str, Any]:
        return self._get_vbd_serial_workspace()


def _require_owner(state: Dict[str, Any]) -> TrussSolverOCL:
    owner = state.get("owner")
    if owner is None or not isinstance(owner, TrussSolverOCL):
        raise ValueError("solver state missing TrussSolverOCL owner; ensure solver was constructed correctly")
    return owner


def solve_vbd_serial_gpu(
    truss: Truss,
    *,
    state: Dict[str, Any],
    dt: float,
    gravity: np.ndarray,
    x: np.ndarray,
    v: np.ndarray,
    fixed: np.ndarray,
    config: Dict[str, Any],
) -> np.ndarray:
    print("truss_solver_ocl.solve_vbd_serial_gpu()")
    solver = _require_owner(state)
    workspace = solver.get_vbd_serial_workspace()
    runtime = solver.runtime
    queue = runtime.queue

    niter = int(config.get("niter", 1))
    det_eps = np.float32(config.get("det_eps", 1e-6))
    if not bool(config.get("serial", True)):
        raise NotImplementedError("Parallel VBD kernel is not yet ported in the refactored solver")

    n_points = solver.n_points
    n_max = solver.n_max
    inv_h2 = np.float32(1.0 / (dt * dt))

    x32 = workspace["x_host"]
    y32 = workspace["y_host"]
    x32.fill(0.0)
    y32.fill(0.0)

    # Load current state and predicted positions into pinned host buffers.
    pos_np = x.astype(np.float32, copy=False)
    y_pred = x + dt * v + (dt * dt) * gravity
    y_pred = y_pred.astype(np.float32, copy=False)
    if fixed.size > 0:
        y_pred[fixed] = x[fixed].astype(np.float32, copy=False)

    mass32 = solver.masses.astype(np.float32, copy=False)
    x32[:, :3] = pos_np
    x32[:, 3] = mass32
    y32[:, :3] = y_pred

    cl.enqueue_copy(queue, workspace["x_buf"], x32)
    cl.enqueue_copy(queue, workspace["y_buf"], y32)

    verbose_iter = int(config.get("verbose", 0))
    n_iter = max(1, niter)

    buf_a = workspace["x_buf"]

    need_logging = bool(verbose_iter) or VERBOSITY >= 2
    if need_logging:
        lerp_prev = x.astype(np.float64, copy=True)
        lerp_curr = np.zeros_like(lerp_prev)
    else:
        lerp_prev = lerp_curr = None

    for itr_idx in range(n_iter):
        event = runtime.prg.vbd_vertex_serial(
            queue,
            (1,),
            None,
            buf_a,
            workspace["y_buf"],
            workspace["neighs_buf"],
            workspace["kngs_buf"],
            workspace["l0_buf"],
            np.int32(n_points),
            np.int32(n_max),
            inv_h2,
            det_eps,
        )
        event.wait()

        if need_logging:
            cl.enqueue_copy(queue, x32, buf_a)
            queue.finish()
            lerp_curr[:, :] = x32[:, :3].astype(np.float64)
            dx_iter = lerp_curr - lerp_prev
            max_dx = np.linalg.norm(dx_iter, axis=1).max()
            if verbose_iter:
                print(f"    VBD iter {itr_idx}: max |dx| = {max_dx:.3e}")
            if VERBOSITY >= 2:
                vel_iter = (lerp_curr - x) / dt
                if fixed.size > 0:
                    vel_iter[fixed] = 0.0
                _print_state(f"[GPU][VBD iter {itr_idx + 1}/{n_iter}]", lerp_curr, vel_iter)
            lerp_prev[:, :] = lerp_curr

    cl.enqueue_copy(queue, x32, buf_a)
    queue.finish()

    x_out = x32[:, :3].astype(np.float64, copy=True)
    if fixed.size > 0:
        x_out[fixed] = x[fixed]
    return x_out


SOLVERS: Dict[str, Callable[..., np.ndarray]] = {
    "vbd_serial": solve_vbd_serial_gpu,
    "vbd": solve_vbd_serial_gpu,
}


def get_solver(name: str) -> Callable[..., np.ndarray]:
    print("truss_solver_ocl.get_solver()")
    if name not in SOLVERS:
        raise ValueError(f"Unknown OCL solver '{name}'; available: {list(SOLVERS.keys())}")
    return SOLVERS[name]
