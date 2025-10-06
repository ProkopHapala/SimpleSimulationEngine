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


def _parse_solver_suite(spec: str) -> List[Tuple[str, Optional[str], Optional[str]]]:
    entries: List[Tuple[str, Optional[str], Optional[str]]] = []
    if not spec:
        return entries
    tokens = [tok.strip() for tok in spec.split(',') if tok.strip()]
    for token in tokens:
        label: Optional[str] = None
        sub_method: Optional[str] = None
        solver_spec = token
        if ':' in solver_spec:
            solver_part, label_part = solver_spec.split(':', 1)
            solver_spec = solver_part.strip()
            label = label_part.strip() or None
        if '[' in solver_spec:
            if not solver_spec.endswith(']'):
                raise ValueError(f"Malformed solver suite entry '{token}' (missing ']')")
            solver_name, sub_part = solver_spec.split('[', 1)
            solver_name = solver_name.strip()
            sub_method = sub_part[:-1].strip()
        else:
            solver_name = solver_spec.strip()
        if not solver_name:
            raise ValueError(f"Empty solver name in suite entry '{token}'")
        entries.append((solver_name, sub_method, label))
    return entries


def _run_solver_once(solver: ts.TrussSolver, args: argparse.Namespace, solver_name: str, label: str,
                     base_positions: np.ndarray, displacement: np.ndarray,
                     fixed_points: List[int]) -> Dict[str, object]:
    request = {
        'solver': solver_name,
        'config': _solver_config(args, label),
        'label': label,
        'estimate_rho': getattr(args, 'estimate_rho', False),
    }
    results = ts.run_solver_suite(
        lambda: solver.truss,
        dt=args.dt,
        displacement=displacement,
        solver_requests=[request],
        base_positions=base_positions,
        fixed_points=fixed_points,
        estimate_rho=getattr(args, 'estimate_rho', False),
        existing_solver=solver,
    )
    return results[label]


def _maybe_export_linear_matrix(args: argparse.Namespace, solver: ts.TrussSolver,
                                base_positions: np.ndarray) -> None:
    if not (args.print_linear_matrix or args.save_linear_matrix_npy or args.save_linear_matrix_txt):
        return

    zero_vel = np.zeros_like(base_positions)
    solver.reset_state(positions=base_positions, velocities=zero_vel)
    solver.ps_pred[:] = base_positions
    solver.ps_cor[:] = base_positions

    matrix = solver._get_linear_matrix().copy()
    shape_info = f"Linearized system matrix shape: {matrix.shape}"

    if args.print_linear_matrix:
        print(shape_info)
        print(matrix)
    else:
        print(shape_info)

    if args.save_linear_matrix_npy:
        np.save(args.save_linear_matrix_npy, matrix)
        print(f"Linear matrix saved to {args.save_linear_matrix_npy} (NumPy binary)")

    if args.save_linear_matrix_txt:
        np.savetxt(args.save_linear_matrix_txt, matrix)
        print(f"Linear matrix saved to {args.save_linear_matrix_txt} (text)")


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

prokop@GTX3090:~/git/SimpleSimulationEngine/python/pyTruss$ python run_solver_debug.py --solver-suite "vbd:VBD,gs_diff:GS-diff,jacobi_diff:Jacobi-diff,jacobi_fly:Jacobi-fly,gs_fly:GS-fly" --csv-path solver_suite.csv --nx 4 --ny 4 --niter 20 

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
    parser.add_argument("--solver-suite",  type=str,   default="",    help="Comma-separated solvers to compare; entries may use 'solver[sub_method]:label'. Overrides --solver/--ref-solver.")
    parser.add_argument("--csv-path",      type=str,   default="",    help="Optional output CSV path for convergence traces.")
    parser.add_argument("--print-linear-matrix", action="store_true", help="Print the assembled linearized solver matrix (Jacobi/GS diff).")
    parser.add_argument("--save-linear-matrix-npy", type=str, default="", help="Path to save the linearized solver matrix via np.save.")
    parser.add_argument("--save-linear-matrix-txt", type=str, default="", help="Path to save the linearized solver matrix via np.savetxt.")

    args = parser.parse_args()

    ts.set_verbosity(args.verb)

    plot_truss = _build_truss(args)

    solver_truss = _build_truss(args)
    base_positions = solver_truss.points.copy()
    fixed_points = sorted(solver_truss.fixed)

    suite_entries = _parse_solver_suite(args.solver_suite)
    if not suite_entries:
        suite_entries = []
        if args.ref_solver.lower() != "none":
            suite_entries.append((args.ref_solver, None, "reference"))
        suite_entries.append((args.solver, None, "solver"))

    matrix_solver = suite_entries[0][0] if suite_entries else args.solver
    shared_solver = ts.TrussSolver(
        solver_truss,
        dt=args.dt,
        gravity=np.zeros(3, dtype=np.float64),
        solver=ts.get_solver(matrix_solver),
        solver_config=_solver_config(args, matrix_solver),
        fixed_points=fixed_points,
        verbose=args.verb,
    )

    _maybe_export_linear_matrix(args, shared_solver, base_positions)

    displacement = _make_displacement(base_positions, fixed_points, args)
    vertices = _parse_vertices(args.plot_vertices, base_positions.shape[0])

    solver_requests: List[Dict[str, object]] = []
    used_labels: Dict[str, int] = {}
    for idx, (solver_name, sub_method_override, label) in enumerate(suite_entries):
        label_final = label or solver_name or f"solver_{idx}"
        count = used_labels.get(label_final, 0)
        if count > 0:
            label_final = f"{label_final}_{count}"
        used_labels[label_final] = count + 1

        cfg = _solver_config(args, label_final)
        if sub_method_override:
            cfg['sub_method'] = sub_method_override

        solver_requests.append({
            'solver': solver_name,
            'config': cfg,
            'label': label_final,
            'estimate_rho': args.estimate_rho,
        })

    def truss_factory() -> Truss:
        return shared_solver.truss

    results = ts.run_solver_suite(
        truss_factory,
        dt=args.dt,
        displacement=displacement,
        solver_requests=solver_requests,
        base_positions=base_positions,
        fixed_points=fixed_points,
        estimate_rho=args.estimate_rho,
        existing_solver=shared_solver,
    )

    if args.csv_path:
        ts.write_convergence_csv(results, args.csv_path)

    _print_summary(results)
    _plot_results(base_positions, displacement, vertices, args.zoom, plot_truss, results, args)

