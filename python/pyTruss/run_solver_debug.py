#!/usr/bin/env python3
"""CPU-only solver debugging utility for truss solvers.

Builds a static truss configuration, applies a controlled perturbation, and
runs selected solver(s) with per-iteration logging enabled. Produces trajectory
and convergence plots to aid debugging of solver behavior.
"""

import argparse
from typing import Dict, List, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt

from truss import Truss
import truss_solver as ts


def _build_truss(args: argparse.Namespace) -> Truss:
    truss = Truss()
    ny_effective = 0 if bool(args.chain) else args.ny
    truss.build_grid_2d(
        nx=args.nx,
        ny=ny_effective,
        m=args.mass,
        m_end=args.anchor_mass,
        l=args.edge,
        k=args.stiffness,
        k_diag=args.diag_stiffness,
    )

    anchor_mode = args.anchor_mode.lower()
    if anchor_mode not in {"none", "left", "right", "both"}:
        raise ValueError(f"unsupported anchor mode: {anchor_mode}")

    if anchor_mode == "none":
        fixed = set()
    elif anchor_mode == "left":
        fixed = {0}
    elif anchor_mode == "right":
        fixed = {args.nx}
    else:
        fixed = {0, args.nx}

    truss.masses[:] = args.mass
    truss.fixed = fixed
    if truss.fixed:
        idx = np.array(sorted(truss.fixed), dtype=int)
        truss.masses[idx] = args.anchor_mass

    return truss


def _parse_vertices(spec: str, n_points: int) -> List[int]:
    if not spec or spec.strip().lower() == "all":
        return list(range(n_points))
    items = [tok.strip() for tok in spec.split(',') if tok.strip()]
    indices = [int(tok) for tok in items]
    for idx in indices:
        if idx < 0 or idx >= n_points:
            raise ValueError(f"plot vertex index {idx} out of range 0..{n_points-1}")
    return indices


def _make_displacement(base: np.ndarray, fixed: List[int], args: argparse.Namespace) -> np.ndarray:
    if args.perturb_path:
        disp = np.load(args.perturb_path)
        if disp.shape != base.shape:
            raise ValueError(f"perturb file shape {disp.shape} != base {base.shape}")
    else:
        rng = np.random.default_rng(args.seed)
        disp = rng.normal(scale=args.perturb_sigma, size=base.shape)
    if fixed:
        disp[np.array(fixed, dtype=int)] = 0.0
    return disp


def _solver_config(args: argparse.Namespace, label: str) -> Dict[str, object]:
    cfg: Dict[str, object] = {
        'niter': args.niter,
        'sub_method': args.sub_method,
        'capture_iters': True,
        'log_label': label,
        'rho': args.rho,
        'delayed_start': args.delayed_start,
        'gamma': args.gamma,
        'b_start': args.b_start,
        'b_end': args.b_end,
        'b_last': args.b_last,
        'istart': args.istart,
        'iend': args.iend if args.iend >= 0 else None,
    }
    return cfg


def _spectral_radius_method(solver_name: str, sub_method: str) -> Optional[str]:
    name = solver_name.lower()
    if name in {'jacobi_diff', 'jacobi'}:
        return 'jacobi_diff'
    if name in {'gs_diff', 'gauss_seidel'}:
        return 'gs_diff'
    if name in {'momentum', 'chebyshev'}:
        sub = sub_method.lower()
        if sub in {'jacobi_diff', 'gs_diff'}:
            return sub
        return None
    return None


def _run_solver_once(args: argparse.Namespace, solver_name: str, label: str,
                     base_positions: np.ndarray, displacement: np.ndarray,
                     fixed_points: List[int]) -> Dict[str, object]:
    truss = _build_truss(args)
    solver_callback = ts.get_solver(solver_name)
    cfg = _solver_config(args, label)

    solver = ts.TrussSolver(
        truss,
        dt=args.dt,
        gravity=np.zeros(3, dtype=np.float64),
        solver=solver_callback,
        solver_config=cfg,
        fixed_points=fixed_points,
        verbose=0,
    )

    base = base_positions.copy()
    initial = base + displacement
    initial_copy = initial.copy()

    solver.x[:] = base
    solver.v[:] = 0.0
    solver.ps_pred[:] = initial
    solver.ps_cor[:] = initial
    if solver.fixed.size > 0:
        solver.ps_pred[solver.fixed] = base[solver.fixed]
        solver.ps_cor[solver.fixed] = base[solver.fixed]

    solver.solver_state['iter_logs'] = {}
    solver.solver_state['log_base'] = solver.ps_pred.copy()

    if getattr(args, 'estimate_rho', False):
        method = _spectral_radius_method(solver_name, args.sub_method)
        if method is not None:
            try:
                rho = ts.estimate_iterative_spectral_radius(solver, method=method)
                print(f"[spectral-radius] {solver_name} ({method}) = {rho:.6f}")
            except Exception as err:
                print(f"[spectral-radius] {solver_name} unavailable: {err}")
        else:
            print(f"[spectral-radius] {solver_name} not supported for estimation")

    solver.solver_callback(solver, cfg)

    logs = solver.solver_state.get('iter_logs', {}).get(label)
    if logs and logs.get('positions'):
        positions = np.stack(logs['positions'], axis=0)
        max_entries = logs.get('max_step', [])
        if max_entries:
            max_step = np.array(max_entries, dtype=np.float64)
        else:
            max_step = np.zeros(positions.shape[0], dtype=np.float64)
        if max_step.size < positions.shape[0]:
            pad = positions.shape[0] - max_step.size
            max_step = np.pad(max_step, (0, pad))
        deltas = np.linalg.norm(positions[1:] - positions[:-1], axis=2)
        residual = deltas.max(axis=1) if deltas.size else np.array([], dtype=np.float64)
    else:
        final_pos = solver.ps_cor.copy()
        positions = np.stack([initial_copy, final_pos], axis=0)
        diff = final_pos - initial_copy
        max_delta = float(np.max(np.abs(diff))) if diff.size else 0.0
        max_step = np.array([0.0, max_delta], dtype=np.float64)
        residual = np.array([max_delta], dtype=np.float64) if diff.size else np.array([], dtype=np.float64)

    return {
        'label': label,
        'solver': solver_name,
        'positions': positions,
        'max_step': max_step,
        'residual': residual,
        'base': base,
        'initial': initial_copy,
        'final': positions[-1],
    }


def _plot_results(base: np.ndarray, displacement: np.ndarray, vertices: List[int], zoom: float,  truss: Truss, results: Dict[str, Dict[str, object]], args: argparse.Namespace) -> None:
    if args.no_plot:
        return

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Skeleton overlay on left subplot
    bonds = np.asarray(getattr(truss, 'bonds', np.zeros((0, 2), dtype=int)), dtype=int)
    if bonds.ndim != 2 or bonds.shape[1] != 2:
        bonds = np.zeros((0, 2), dtype=int)

    def _plot_skeleton(points: np.ndarray, *, color: str, linestyle: str = '-', linewidth: float = 1.0, marker: str = 'o', label: Optional[str] = None, size: float = 15) -> None:
        if bonds.size:
            for i, j in bonds:
                axes[0].plot([points[i, 0], points[j, 0]], [points[i, 1], points[j, 1]],  color=color, linestyle=linestyle, linewidth=linewidth, alpha=0.9)
        axes[0].scatter(points[:, 0], points[:, 1], s=size, color=color, marker=marker, label=label, alpha=0.9)

    initial = base
    perturbed = base + displacement * zoom
    _plot_skeleton(initial, color='0.3', label='base')
    _plot_skeleton(perturbed, color='tab:orange', linestyle='--', label='perturbed')

    for idx, (key, data) in enumerate(results.items()):
        positions = data['positions']  # (iters, n, 3)
        disp = (positions - base[np.newaxis]) * zoom
        color = colors[idx % len(colors)]

        final = data['final']
        final_scaled = base + (final - base) * zoom
        _plot_skeleton(final_scaled, color=color, linewidth=1.5, marker='x', label=f"{data['solver']} final", size=25)

        for j, vid in enumerate(vertices):
            lbl = data['solver'] if j == 0 else None
            traj = base[vid] + disp[:, vid, :]
            axes[0].plot(traj[:, 0], traj[:, 1], marker='o', ms=3, color=color, label=lbl, alpha=0.7)

        steps = np.arange(len(data['max_step']))
        axes[1].semilogy(steps, data['max_step'], marker='o', color=color, label=data['solver'])
        if data['residual'].size:
            axes[1].semilogy(steps[1:], data['residual'], linestyle='--', color=color, label=f"{data['solver']} Δp")

    axes[0].set_title("Trajectories & skeleton overlay")
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("y")
    axes[0].axis('equal')
    axes[0].legend(loc='best')

    axes[1].set_title("Per-iteration convergence")
    axes[1].set_xlabel("Iteration")
    axes[1].set_ylabel("max step / residual")
    axes[1].legend()
    axes[1].grid(True, which='both', linestyle=':')

    fig.tight_layout()
    if args.savefig:
        fig.savefig(args.savefig, dpi=150)
    plt.show()


def _print_summary(results: Dict[str, Dict[str, object]]) -> None:
    for key, data in results.items():
        solver = data['solver']
        steps = len(data['max_step']) - 1
        final_step = data['max_step'][-1] if data['max_step'].size else 0.0
        print(f"[{solver}] iterations={steps} final max |Δp|={final_step:.3e}")

'''
python run_solver_debug.py --nx 2 --ny 0 --niter 10 --solver momentum --sub-method jacobi_diff --ref-solver gs_diff  --verb 2
python run_solver_debug.py --nx 2 --ny 0 --niter 10 --solver jacobi_diff --ref-solver gs_diff  --verb 2

python run_solver_debug.py --nx 5 --ny 5 --anchor-mode both --niter 20 --solver momentum --sub-method gs_diff --ref-solver vbd  --verb 2

python run_solver_debug.py --nx 5 --ny 5 --anchor-mode both --niter 20 --solver gs_fly --ref-solver gs_diff  --verb 2

python run_solver_debug.py --nx 5 --ny 5 --anchor-mode both --niter 20 --solver gs_diff --ref-solver vbd  --verb 2

'''

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CPU-only solver debugging harness.")
    parser.add_argument("--nx",            type=int,   default=2,       help="Grid cells in X direction (nx+1 vertices).")
    parser.add_argument("--ny",            type=int,   default=0,      help="Grid cells in Y direction (ny+1 vertices).")
    parser.add_argument("--edge",          type=float, default=1.0,    help="Rest edge length (m).")
    parser.add_argument("--stiffness",     type=float, default=1000.0, help="Spring stiffness for axis edges.")
    parser.add_argument("--diag-stiffness",type=float, default=100.0,    help="Spring stiffness for diagonal edges.")
    parser.add_argument("--mass",          type=float, default=1.0,    help="Mass for non-anchored vertices.")
    parser.add_argument("--anchor-mass",   type=float, default=1e6,    help="Mass for anchored vertices.")
    parser.add_argument("--anchor-mode",   type=str,   default="left", choices=["none", "left", "right", "both"], help="Anchoring pattern.")
    parser.add_argument("--chain",         type=int,   default=0,      help="Force ny=0 to simulate a chain.")
    parser.add_argument("--dt",            type=float, default=0.05,   help="Pseudo time-step used for solver mass scaling.")
    parser.add_argument("--solver",        type=str,   default="jacobi_diff", help="Test solver name (CPU).")
    parser.add_argument("--ref-solver",    type=str,   default="none", help="Optional reference solver name (CPU).")
    parser.add_argument("--sub-method",    type=str,   default="jacobi_diff", choices=["jacobi_diff", "gs_diff", "jacobi_fly", "gs_fly"], help="Sub-solver for iterative accelerators.")
    parser.add_argument("--niter",         type=int,   default=10,    help="Maximum iterations per solver call.")
    parser.add_argument("--rho",           type=float, default=0.3,   help="Chebyshev spectral radius estimate.")
    parser.add_argument("--delayed-start", type=int,   default=5,     dest="delayed_start", help="Chebyshev delay (iterations).")
    parser.add_argument("--gamma",         type=float, default=1.0,   help="Chebyshev under-relaxation factor.")
    
    parser.add_argument("--b-start",       type=float, default=0.2,   help="Momentum mixing start value.")
    parser.add_argument("--b-end",         type=float, default=0.2,   help="Momentum mixing end value.")

    #parser.add_argument("--b-start",       type=float, default=0.15,   help="Momentum mixing start value.")
    #parser.add_argument("--b-end",         type=float, default=0.25,   help="Momentum mixing end value.")

    parser.add_argument("--b-last",        type=float, default=0.0,   help="Momentum mixing value for final iteration.")
    parser.add_argument("--istart",        type=int,   default=5,     help="Momentum mixing start iteration (inclusive).")
    parser.add_argument("--iend",          type=int,   default=-1,    help="Momentum mixing end iteration (inclusive, -1 for auto).")
    parser.add_argument("--seed",          type=int,   default=42,    help="Random seed for perturbation generation.")
    parser.add_argument("--perturb-sigma", type=float, default=0.001,  help="Std-dev of random perturbation (in meters).")
    parser.add_argument("--perturb-path",  type=str,   default="",    help="Optional .npy file with explicit displacement array.")
    parser.add_argument("--zoom",          type=float, default=50.0,  help="Visualization zoom multiplier for displacement trajectories.")
    parser.add_argument("--plot-vertices", type=str,   default="all", help="Comma-separated vertex indices to plot trajectories for (or 'all').")
    parser.add_argument("--savefig",       type=str,   default="",    help="Optional output path for the plot figure.")
    parser.add_argument("--verb",          type=int,   default=0,     help="Set truss solver verbosity (0 quiet, 1 steps, 2 per-iter).")
    parser.add_argument("--no-plot", action="store_true",             help="Disable plotting.")
    parser.add_argument("--estimate-rho", action="store_true",        help="Estimate and print spectral radius for supported solvers.")

    args = parser.parse_args()

    ts.set_verbosity(args.verb)

    base_truss = _build_truss(args)
    base_positions = base_truss.points.copy()
    fixed_points = sorted(base_truss.fixed)

    displacement = _make_displacement(base_positions, fixed_points, args)
    vertices = _parse_vertices(args.plot_vertices, base_positions.shape[0])

    results: Dict[str, Dict[str, object]] = {}

    if args.ref_solver.lower() != "none":
        ref_label = "reference"
        ref_result = _run_solver_once(args, args.ref_solver, ref_label, base_positions, displacement, fixed_points)
        results['reference'] = ref_result

    test_label = "solver"
    test_result = _run_solver_once(args, args.solver, test_label, base_positions, displacement, fixed_points)
    results['test'] = test_result

    _print_summary(results)
    _plot_results(base_positions, displacement, vertices, args.zoom, base_truss, results, args)
