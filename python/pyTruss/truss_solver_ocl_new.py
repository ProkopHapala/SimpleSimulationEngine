import os, sys
import numpy as np
import pyopencl as cl
from typing import Any, Callable, Dict, List, Optional, Tuple

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'pyMolecular'))
from OCL.OpenCLBase import OpenCLBase
from OCL import clUtils as clu  # ensures package initialisation for other modules using OpenCLBase
from sparse import build_neighbor_list, neigh_stiffness, neighs_to_dense_arrays
from truss import Truss

VERBOSITY = 0


def set_verbosity(level: int) -> None:
    global VERBOSITY
    VERBOSITY = int(level)


def roundup_global_size(n, nloc): return ((n + nloc - 1) // nloc) * nloc

class TrussSolverOCL(OpenCLBase):
    """GPU-accelerated truss dynamics solver. Follows MolecularDynamics.py design: clear separation of buffer allocation and kernel execution, pre-generated kernel args, Python-side MD loop with GPU linear solvers."""
    
    def __init__(self, truss: Truss, dt: float = 0.01, gravity: Optional[np.ndarray] = None, fixed_points: Optional[List[int]] = None, nloc: int = 32, device_index: int = 0):
        super().__init__(nloc=nloc, device_index=device_index)
        self.truss    = truss
        self.dt       = float(dt)
        self.gravity  = np.array(gravity if gravity is not None else [0.0, 0.0, -9.8], dtype=np.float64)
        self.n_points = truss.points.shape[0]
        self.bonds    = np.asarray(truss.bonds, dtype=np.int32).reshape(-1, 2) if len(truss.bonds) > 0 else np.zeros((0, 2), dtype=np.int32)
        self.ks = truss.ks.astype(np.float32, copy=False)
        if self.ks.ndim == 0: self.ks = np.full(self.bonds.shape[0], float(self.ks), dtype=np.float32)
        self.masses   = truss.masses.astype(np.float64, copy=False)
        self.rest_lengths = truss.get_rest_lengths().astype(np.float32, copy=False)
        fixed_seq   = sorted(fixed_points if fixed_points is not None else truss.fixed)
        self.fixed  = np.asarray(fixed_seq, dtype=np.int32) if len(fixed_seq) > 0 else np.zeros(0, dtype=np.int32)
        self.x      = truss.points.astype(np.float64, copy=True)
        self.v      = np.zeros_like(self.x)
        self.x0     = self.x.copy()
        kernel_path = os.path.join(os.path.dirname(__file__), "truss.cl")
        if not self.load_program(kernel_path=kernel_path, bPrint=False): raise RuntimeError(f"Failed to load: {kernel_path}")
        self._prepare_topology()
        self._allocate_buffers()
        self._upload_static_data()
        self._setup_kernels()
        if VERBOSITY >= 1:
            print(f"TrussSolverOCL init: {self.n_points} points, {len(self.bonds)} bonds, n_max={self.n_max}")

    def _prepare_topology(self):
        if len(self.bonds) == 0: raise ValueError("Truss must have at least one bond")
        neighbs             = build_neighbor_list(self.bonds, self.n_points)
        neighs, kngs, n_max = neigh_stiffness(neighbs, self.bonds, self.ks)
        padded_neighs, padded_kngs, _ = neighs_to_dense_arrays(neighs, kngs, n_max)
        self.neigh_lists, self.n_max  = neighs, n_max
        self.neighs_dense             = padded_neighs.astype(np.int32, copy=True)
        self.kngs_dense               = padded_kngs.astype(np.float32, copy=True)
        pair_rest = {}
        for (ia, ib), rest in zip(self.bonds, self.rest_lengths): 
            pair_rest[(int(ia), int(ib))] = pair_rest[(int(ib), int(ia))] = float(rest)
        self.rest_dense = np.zeros((self.n_points, n_max), dtype=np.float32)
        for vid, neigh in enumerate(neighs):
            for jj, nj in enumerate(neigh): 
                self.rest_dense[vid, jj] = pair_rest[(vid, int(nj))]
        self.Aii = np.array([np.sum(k) for k in kngs], dtype=np.float32)

    def _allocate_buffers(self):
        mf = cl.mem_flags
        fs, ds, i_s = np.float32().itemsize, np.float64().itemsize, np.int32().itemsize
        self.create_buffer('pos',       self.n_points * 4 * fs, mf.READ_WRITE)
        self.buffer_dict['pos_f'] = self.buffer_dict['pos']
        self.create_buffer('pos_d',     self.n_points * 4 * ds, mf.READ_WRITE)
        self.create_buffer('y_pred',    self.n_points * 4 * fs, mf.READ_WRITE)
        self.create_buffer('neighs',    self.n_points * self.n_max * i_s,  mf.READ_ONLY)
        self.create_buffer('kngs',      self.n_points * self.n_max * fs, mf.READ_ONLY)
        self.create_buffer('l0s',       self.n_points * self.n_max * fs, mf.READ_ONLY)
        self.buffer_dict['l0ngs'] = self.buffer_dict['l0s']
        self.buffer_dict['rest_lengths'] = self.buffer_dict['l0s']
        self.create_buffer('Aii_buf',   self.n_points * fs, mf.READ_WRITE)
        self.create_buffer('bvec',      self.n_points * 4 * fs, mf.READ_WRITE)
        self.create_buffer('x_sol',     self.n_points * 4 * fs, mf.READ_WRITE)
        self.create_buffer('x_sol2',    self.n_points * 4 * fs, mf.READ_WRITE)
        self.create_buffer('residual',  self.n_points * 4 * fs, mf.READ_WRITE)
        self.create_buffer('dps_store', self.n_points * 4 * fs, mf.READ_WRITE)
        self.create_buffer('dp',        self.n_points * 4 * fs, mf.READ_WRITE)
        self.create_buffer('kfix',      self.n_points * fs, mf.READ_WRITE)
        self.pos_host = np.zeros((self.n_points, 4), dtype=np.float32)
        self.vel_host = np.zeros((self.n_points, 4), dtype=np.float32)
        self.buffer_dict['x'] = self.buffer_dict['pos']
        self.buffer_dict['y'] = self.buffer_dict['y_pred']

    def _upload_static_data(self):
        self.toGPU('neighs', self.neighs_dense.flatten())
        self.toGPU('kngs', self.kngs_dense.flatten())
        self.toGPU('l0s', self.rest_dense.flatten())
        self.toGPU('Aii_buf', self.Aii)
        kfix = np.zeros(self.n_points, dtype=np.float32)
        if len(self.fixed) > 0: kfix[self.fixed] = 1e10
        self.toGPU('kfix', kfix)
        self.queue.finish()

    def _setup_kernels(self):
        self.kernel_params = {
            'nmax': np.int32(self.n_max),
            'nverts': np.int32(self.n_points),
            'n': np.int32(self.n_points),
            'npoint': np.int32(self.n_points),
            'dt': np.float32(self.dt),
            'inv_h2': np.float32(1.0/(self.dt*self.dt)),
            'det_eps': np.float32(1e-6),
            'nIters': np.int32(1),
            'bmix': np.float32(0.1),
            'nmax_neigh': np.int32(self.n_max),
            'n_points': np.int32(self.n_points),
        }
        self.sz_loc    = (self.nloc,)
        self.sz_global = (roundup_global_size(self.n_points, self.nloc),)
        self.sz_single = (1,)
        required = ['pos', 'pos_f', 'pos_d', 'y_pred', 'neighs', 'kngs', 'l0s', 'Aii_buf', 'bvec', 'x_sol', 'x_sol2', 'residual', 'dps_store', 'dp', 'kfix']
        for name in required:
            if name not in self.buffer_dict:
                self.buffer_dict[name] = None

    def step(self, solver_callback: Callable, solver_config: Optional[Dict[str, Any]] = None):
        config  = solver_config if solver_config is not None else {}
        x_prev = self.x.copy()
        x_new  = solver_callback(self, config)
        v_new  = (x_new - x_prev) / self.dt
        if len(self.fixed) > 0:
            x_new[self.fixed] = x_prev[self.fixed]
            v_new[self.fixed] = 0.0
        self.x = x_new
        self.v = v_new
    
    def run(self, nsteps: int, solver_callback: Callable, solver_config: Optional[Dict[str, Any]] = None, track_indices: Optional[List[int]] = None, verbose: bool = False):
        trajectory = [] if track_indices is not None else None
        track_idx = np.asarray(track_indices, dtype=int) if track_indices is not None else None
        if trajectory is not None and track_idx is not None: trajectory.append(self.x[track_idx].copy())
        for step in range(nsteps):
            self.step(solver_callback, solver_config)
            if trajectory is not None and track_idx is not None: 
                trajectory.append(self.x[track_idx].copy())
            if verbose and (step == 0 or (step + 1) % max(1, nsteps // 10) == 0): 
                print(f"  Step {step + 1}/{nsteps}, t={self.dt * (step + 1):.3f}")
            if VERBOSITY >= 1:
                print(f"[OCL][step {step + 1}/{nsteps}] min|x|={np.linalg.norm(self.x, axis=1).min():.3e}")
        self.truss.points = self.x.copy()
        return self.x, self.v, (np.stack(trajectory, axis=0) if trajectory else None)

    def run_precompute_dRHS(self, use_double: bool = True):
        if use_double:
            pos_d = np.zeros((self.n_points, 4), dtype=np.float64); 
            pos_d[:, :3] = self.x; 
            pos_d[:, 3]  = self.masses
            self.toGPU('pos_d', pos_d)
            args = self.generate_kernel_args('precompute_dRHS_d', bPrint=False)
            self.prg.precompute_dRHS_d(self.queue, self.sz_global, self.sz_loc, *args)
        else:
            args = self.generate_kernel_args('precompute_dRHS', bPrint=False)
            self.prg.precompute_dRHS(self.queue, self.sz_global, self.sz_loc, *args)
        self.queue.finish()

    def run_vbd_serial(self, niter: int = 1):
        for _ in range(niter): 
            args = self.generate_kernel_args('vbd_vertex_serial', bPrint=False)
            self.prg.vbd_vertex_serial(self.queue, self.sz_single, None, *args); 
            self.queue.finish()

    def run_jacobi_fly(self, niter: int = 10):
        self.kernel_params['nIters'] = np.int32(niter); 
        args = self.generate_kernel_args('jacobi_fly', bPrint=False)
        self.prg.jacobi_fly(self.queue, self.sz_global, self.sz_loc, *args); self.queue.finish()

    def run_jacobi_diff(self, niter: int = 10):
        self.kernel_params['nIters'] = np.int32(niter); 
        args = self.generate_kernel_args('jacobi_iteration_diff_serial', bPrint=False)
        self.prg.jacobi_iteration_diff_serial(self.queue, self.sz_global, self.sz_loc, *args); self.queue.finish()

    def download_positions(self):
        self.fromGPU('pos_f', self.pos_host); return self.pos_host[:, :3].astype(np.float64)
    
    def upload_positions(self, x: np.ndarray):
        self.x               = x.astype(np.float64, copy=True)
        self.pos_host[:, :3] = x.astype(np.float32); 
        self.pos_host[:, 3]  =  self.masses.astype(np.float32)
        self.toGPU('pos_f', self.pos_host)
        pos_d = np.zeros((self.n_points, 4), dtype=np.float64); pos_d[:, :3] = x; pos_d[:, 3] = self.masses
        if 'pos_d' in self.buffer_dict: self.toGPU('pos_d', pos_d)
        self.queue.finish()
    
    def reset_state(self, positions: Optional[np.ndarray] = None, velocities: Optional[np.ndarray] = None):
        if positions  is not None: self.upload_positions(positions)
        if velocities is not None: self.v = velocities.astype(np.float64, copy=True)
        else: self.v = np.zeros_like(self.x)

# Solver callbacks
def solve_vbd_serial(solver: TrussSolverOCL, config: Dict[str, Any]) -> np.ndarray:
    niter, dt = int(config.get('niter', 10)), solver.dt
    y_pred = solver.x + dt * solver.v + (dt * dt) * solver.gravity
    if len(solver.fixed) > 0: y_pred[solver.fixed] = solver.x[solver.fixed]
    pos32, y32                = np.zeros((solver.n_points, 4), dtype=np.float32), np.zeros((solver.n_points, 4), dtype=np.float32)
    pos32[:, :3]              = solver.x.astype(np.float32)
    pos32[:,  3]              = solver.masses.astype(np.float32)
    y32  [:, :3]              = y_pred.astype(np.float32)
    solver.toGPU('pos_f', pos32); solver.toGPU('y_pred', y32)
    solver.kernel_params['inv_h2'] = np.float32(1.0/(dt*dt)); solver.kernel_params['det_eps'] = np.float32(config.get('det_eps', 1e-6))
    solver.run_vbd_serial(niter=niter)
    x_out = solver.download_positions()
    if len(solver.fixed) > 0: x_out[solver.fixed] = solver.x[solver.fixed]
    return x_out

def solve_jacobi_fly(solver: TrussSolverOCL, config: Dict[str, Any]) -> np.ndarray:
    niter        = int(config.get('niter', 10))
    dt           = solver.dt
    y_pred       = solver.x + dt * solver.v + (dt * dt) * solver.gravity
    pos32        = np.zeros((solver.n_points, 4), dtype=np.float32); 
    pos32[:, :3] = y_pred.astype(np.float32)
    pos32[:,  3] = solver.masses.astype(np.float32)
    solver.toGPU('pos_f', pos32)
    solver.kernel_params['dt'] = np.float32(dt); 
    solver.run_jacobi_fly(niter=niter)
    x_out  = solver.download_positions()
    if len(solver.fixed) > 0: x_out[solver.fixed] = solver.x[solver.fixed]
    return x_out

def solve_jacobi_diff(solver: TrussSolverOCL, config: Dict[str, Any]) -> np.ndarray:
    niter        = int (config.get('niter', 10))
    dt           = solver.dt
    use_double   = bool(config.get('use_double', True))
    y_pred       = solver.x + dt * solver.v + (dt * dt) * solver.gravity
    solver.upload_positions(y_pred)
    solver.kernel_params['dt'] = np.float32(dt); 
    solver.run_precompute_dRHS(use_double=use_double)
    dp_init      = np.zeros((solver.n_points, 4), dtype=np.float32); 
    solver.toGPU('dp', dp_init)
    solver.run_jacobi_diff(niter=niter)
    dp_out       = np.zeros((solver.n_points, 4), dtype=np.float32); 
    solver.fromGPU('dp', dp_out)
    x_out        = y_pred + dp_out[:, :3].astype(np.float64)
    if len(solver.fixed) > 0: x_out[solver.fixed] = solver.x[solver.fixed]
    return x_out

SOLVERS = {
    'vbd_serial': solve_vbd_serial,
    'vbd': solve_vbd_serial,
    'jacobi_fly': solve_jacobi_fly,
    'jacobi_diff': solve_jacobi_diff,
}
def get_solver(name: str) -> Callable: 
    if name not in SOLVERS: raise ValueError(f"Unknown solver '{name}'; available: {list(SOLVERS.keys())}")
    return SOLVERS[name]
