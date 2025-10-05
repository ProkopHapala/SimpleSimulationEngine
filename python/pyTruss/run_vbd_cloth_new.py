import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from truss import Truss

from truss_solver import TrussSolver, get_solver
import truss_solver as truss_solver_module

# from   truss_solver_bak3 import TrussSolver, get_solver
# import truss_solver_bak3 as     truss_solver_module

import truss_solver_ocl_new as truss_solver_ocl_mod
from truss_solver_ocl_new import TrussSolverOCL as TrussSolverOCLGPU, get_solver as get_solver_ocl
from plot_utils import plot_truss


CPU_SOLVER_CHOICES = set(truss_solver_module.SOLVERS.keys())
GPU_SOLVER_CHOICES = set(truss_solver_ocl_mod.SOLVERS.keys())
ALL_SOLVER_CHOICES = tuple(sorted(CPU_SOLVER_CHOICES | GPU_SOLVER_CHOICES))


def plot_graph_coloring(truss: Truss, vertex_colors: np.ndarray, partitions, ax=None, cmap_name: str = "tab20"):
    if vertex_colors.size == 0:
        return ax
    if ax is None:
        ax = plt.gca()
    n_colors = int(vertex_colors.max()) + 1
    palette = cm.get_cmap(cmap_name, n_colors)(np.linspace(0.0, 1.0, n_colors))
    node_colors = palette[vertex_colors.astype(int)]
    ax = plot_truss(truss.points, truss.bonds, ax=ax, edge_color='0.3', edge_alpha=0.5, point_size=60, node_colors=node_colors)
    for color_id, verts in enumerate(partitions):
        p = truss.points[verts[0], :2]
        ax.text(p[0], p[1], str(color_id), color='k', fontsize=10, ha='center', va='center')
    ax.set_title("Graph coloring (color id annotated)")
    return ax

'''
python run_vbd_cloth.py --nx 2 --ny 0 --nsteps 100 --niter 10 --stiffness 10000.0 --cpu 1 --anchor-mode left
python run_vbd_cloth.py --nx 4 --ny 0 --nsteps 100 --niter 10 --stiffness 10000.0 --cpu 1 --anchor-mode lef

python run_vbd_cloth.py --nx 2 --ny 0 --nsteps 10 --niter 10 --stiffness 1000.0 --cpu 1 --anchor-mode left 
python run_vbd_cloth.py --nx 2 --ny 0 --nsteps 10 --niter 10 --stiffness 100.0  --cpu 1 --anchor-mode left 
python run_vbd_cloth.py --nx 2 --ny 0 --nsteps 10 --niter 10 --stiffness 10.0   --cpu 1 --anchor-mode left 

python run_vbd_cloth.py --nx 3 --ny 0 --nsteps 10 --niter 10 --stiffness 1000.0 --cpu 1 --anchor-mode both 
python run_vbd_cloth.py --nx 3 --ny 0 --nsteps 10 --niter 10 --stiffness 100.0  --cpu 1 --anchor-mode both 
python run_vbd_cloth.py --nx 3 --ny 0 --nsteps 10 --niter 10 --stiffness 10.0   --cpu 1 --anchor-mode both 


'''

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="Run a simple Vertex Block Descent step on a cloth-like truss grid." )
    parser.add_argument("--nx",            type=int,   default=1, help="Number of cells in the X direction (grid has nx+1 points).")
    parser.add_argument("--ny",            type=int,   default=0, help="Number of cells in the Y direction (grid has ny+1 points).")
    parser.add_argument("--dt",            type=float, default=0.05, help="Simulation time step (seconds).")
    parser.add_argument("--mass",          type=float, default=1.0, help="Mass assigned to interior grid vertices.")
    parser.add_argument("--anchor-mass",   type=float, default=1.0e+6, help="Mass assigned to fixed corner vertices (handled by Truss.build_grid_2d).")
    parser.add_argument("--edge",          type=float, default=1.0, help="Rest length of grid edges (meters).")
    parser.add_argument("--stiffness",     type=float, default=1000000.0, help="Spring stiffness for axis-aligned edges.")
    parser.add_argument("--diag-stiffness",type=float, default=0.0, help="Spring stiffness for diagonal edges (if >0).")
    parser.add_argument("--gravity",       type=float, nargs=3, default=(0.0, -9.81, 0.0), help="Gravity vector (m/s^2).")
    parser.add_argument("--extra-accel",   type=float, nargs=3, default=(0.0, 0.0, 0.0), help="Additional acceleration applied to all vertices (m/s^2).")
    parser.add_argument("--nsteps",        type=int,   default=10, help="Number of simulation time steps.")
    parser.add_argument("--niter",         type=int,   default=10, help="Number of VBD solver iterations per step.")
    parser.add_argument("--det-eps",       type=float, default=1e-6, help="Determinant threshold for skipping ill-conditioned Hessians.")
    parser.add_argument("--solver",        type=str,   default="vbd", choices=ALL_SOLVER_CHOICES,help=f"Solver name (shared choice). Available: {', '.join(ALL_SOLVER_CHOICES)}.")
    parser.add_argument("--verbose",       type=int,   default=1, help="Print iteration progress from the solver.")
    parser.add_argument("--no-pin",        type=int,   default=0, help="Do not pin the default corner vertices; treat all vertices as free.")
    parser.add_argument("--no-plot",       type=int,   default=0, help="Skip plotting the initial and final grids.")
    parser.add_argument("--savefig",       type=str,   default="run_vbd_cloth.png", help="Path to save the comparison plot instead of displaying it.")
    parser.add_argument("--nloc",          type=int,   default=32, help="OpenCL local work size for new GPU solver.")
    parser.add_argument("--device",        type=int,   default=0, help="OpenCL device index for new GPU solver.")
    parser.add_argument("--cpu",           type=int,   default=0, help="Also run CPU solver for comparison.")
    parser.add_argument("--gpu",           type=int,   default=-1, help="Run GPU solver (defaults to opposite of --cpu when -1).")
    parser.add_argument("--chain",         type=int,   default=0, help="Treat the grid as a 1D chain (forces ny=0).")
    parser.add_argument("--anchor-mode",   type=str,   default="left", choices=["none", "left", "right", "both"], help="Override endpoint anchoring.")
    parser.add_argument("--track",         type=str,   default="all", help="Comma-separated vertex indices to plot trajectories for (CPU solver only).")
    parser.add_argument("--verb",          type=int,   default=0, help="Diagnostic verbosity: 0=off, 1=per step, 2=per solver iteration.")
    parser.add_argument("--test-coloring", type=int,   default=0, help="Run graph coloring test and plotting before simulation.")
    parser.add_argument("--color-seed",    type=int,   default=-1, help="Seed for graph coloring RNG (negative for random).")
    args = parser.parse_args()

    truss = Truss()

    verb_level = max(0, int(args.verb))
    truss_solver_module.set_verbosity(verb_level)
    truss_solver_ocl_mod.set_verbosity(verb_level)

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

    coloring_enabled = bool(args.test_coloring)
    coloring_seed = None if args.color_seed < 0 else int(args.color_seed)
    vertex_colors = np.array([], dtype=int)
    color_partitions = []
    if coloring_enabled:
        vertex_colors, color_partitions = truss.color_graph(seed=coloring_seed)
        truss.verify_graph_coloring(vertex_colors)
        n_colors = 0 if vertex_colors.size == 0 else int(vertex_colors.max()) + 1
        print(f"Graph coloring succeeded: {n_colors} colors")
        print(f"  partition sizes: {[len(part) for part in color_partitions]}")
        if not bool(args.no_plot) and vertex_colors.size:
            fig_color, ax_color = plt.subplots(figsize=(6, 6))
            plot_graph_coloring(truss, vertex_colors, color_partitions, ax=ax_color)
            plt.tight_layout()
        plt.show()
        exit()

    anchor_mode = args.anchor_mode.lower()
    if anchor_mode not in {"none", "left", "right", "both"}:
        raise ValueError(f"unsupported anchor mode: {anchor_mode}")

    if anchor_mode == "none":
        new_fixed = set()
    elif anchor_mode == "left":
        new_fixed = {0}
    elif anchor_mode == "right":
        new_fixed = {args.nx}
    else:  # both
        new_fixed = {0, args.nx}

    truss.masses[:] = args.mass
    truss.fixed = new_fixed
    if truss.fixed:
        idx = np.array(sorted(truss.fixed), dtype=int)
        truss.masses[idx] = args.anchor_mass

    no_pin = bool(args.no_pin)
    fixed_points = None if no_pin else list(truss.fixed)

    gravity = np.asarray(args.gravity, dtype=np.float32)
    extra_accel = np.asarray(args.extra_accel, dtype=np.float32)
    total_accel = gravity + extra_accel

    initial_points = truss.points.copy()

    track_indices = None
    if args.track:
        if args.track.strip().lower() == "all":
            track_indices = list(range(len(truss.points)))
        else:
            items = [tok.strip() for tok in args.track.split(',') if tok.strip()]
            if items:
                track_indices = [int(tok) for tok in items]
                for idx in track_indices:
                    if idx < 0 or idx >= len(truss.points):
                        raise ValueError(f"track vertex index {idx} out of range (0..{len(truss.points)-1})")

    run_cpu = bool(args.cpu)
    if int(args.gpu) >= 0:
        run_gpu = bool(args.gpu)
    else:
        run_gpu = not run_cpu

    solver_name = args.solver
    if run_cpu and solver_name not in CPU_SOLVER_CHOICES:
        raise ValueError(f"CPU solver '{solver_name}' not available; choices: {sorted(CPU_SOLVER_CHOICES)}")
    if run_gpu:
        gpu_solver_set = GPU_SOLVER_CHOICES
        desired_gpu_solver = solver_name
        if desired_gpu_solver not in gpu_solver_set:
            raise ValueError(f"GPU solver '{solver_name}' not available for new OpenCL path; GPU provides: {sorted(gpu_solver_set)}")

    cpu_positions = cpu_velocities = None
    gpu_positions = gpu_velocities = None
    cpu_traj = gpu_traj = None

    if run_cpu:
        solver_callback = get_solver(args.solver)
        solver_config = {'verbose': args.verbose}
        if args.solver in {"vbd", "vbd_serial"}:
            solver_config.update({'niter': args.niter, 'det_eps': args.det_eps})
        elif args.solver in {"momentum", "jacobi_diff", "jacobi_diff_cpu", "gs_diff", "gs_diff_cpu",
                             "jacobi_fly", "jacobi_fly_cpu", "gs_fly", "gs_fly_cpu"}:
            solver_config['niter'] = args.niter
        cpu_solver = TrussSolver(
            truss,
            dt=args.dt,
            gravity=total_accel.astype(np.float64),
            solver=solver_callback,
            solver_config=solver_config,
            fixed_points=fixed_points,
            track_indices=track_indices,
            verbose=0,
        )
        cpu_positions, cpu_velocities, cpu_traj = cpu_solver.run(args.nsteps)

    if run_gpu:
        if run_cpu:
            truss.points = initial_points.copy()
        gpu_solver_name = desired_gpu_solver
        gpu_solver_callback = get_solver_ocl(gpu_solver_name)
        gpu_config = {"verbose": args.verbose, "niter": args.niter, "det_eps": args.det_eps}
        gpu_solver = TrussSolverOCLGPU(
            truss,
            dt=args.dt,
            gravity=total_accel.astype(np.float64),
            fixed_points=fixed_points,
            nloc=int(args.nloc),
            device_index=int(args.device),
        )
        gpu_positions, gpu_velocities, gpu_traj = gpu_solver.run(
            args.nsteps,
            solver_callback=gpu_solver_callback,
            solver_config=gpu_config,
            track_indices=track_indices,
            verbose=bool(args.verbose),
        )

    if run_cpu and gpu_positions is not None:
        diff_pos = gpu_positions - cpu_positions
        diff_vel = None
        if gpu_velocities is not None and cpu_velocities is not None:
            diff_vel = gpu_velocities - cpu_velocities
        print("Comparison CPU vs GPU:")
        print("  |pos| max =", np.linalg.norm(diff_pos, axis=1).max())
        print("  pos min =", diff_pos.min(axis=0))
        print("  pos max =", diff_pos.max(axis=0))
        if diff_vel is not None:
            print("  |vel| max =", np.linalg.norm(diff_vel, axis=1).max())
            print("  vel min =", diff_vel.min(axis=0))
            print("  vel max =", diff_vel.max(axis=0))

    positions = cpu_positions if run_cpu else gpu_positions
    velocities = cpu_velocities if run_cpu else gpu_velocities
    traj_to_plot = None
    if track_indices is not None:
        if run_cpu and cpu_traj is not None:
            traj_to_plot = cpu_traj
        elif gpu_traj is not None:
            traj_to_plot = gpu_traj

    if fixed_points is not None and len(fixed_points) > 0:
        positions[fixed_points] = initial_points[fixed_points]
        if velocities is not None:
            velocities[fixed_points] = 0.0

    min_y = float(np.min(positions[:, 1]))
    max_drop = float(np.max(truss.points[:, 1] - positions[:, 1]))
    delta = positions - initial_points
    delta_min = delta.min(axis=0)
    delta_max = delta.max(axis=0)
    print(f"  delta position min: {delta_min}")
    print(f"  delta position max: {delta_max}")

    print(f"Simulation finished: {args.nsteps} steps, {args.niter} solver iterations/step")
    print(f"  lowest y-position: {min_y:.6f} m")
    print(f"  max vertical drop: {max_drop:.6f} m")
    print(f"  pinned vertices: {('none' if fixed_points is None else fixed_points)}")
    if velocities is not None:
        vel_min = velocities.min(axis=0)
        vel_max = velocities.max(axis=0)
        print(f"  velocity min: {vel_min}")
        print(f"  velocity max: {vel_max}")

    if not bool(args.no_plot):
        fig, ax = plt.subplots(figsize=(8, 8))
        plot_truss(initial_points, truss.bonds, ax=ax, edge_color='tab:blue', point_color='tab:blue', label='initial')
        plot_truss(positions, truss.bonds, ax=ax, edge_color='tab:red', point_color='tab:red', label='final')

        if traj_to_plot is not None and track_indices is not None:
            for i, vid in enumerate(track_indices):
                label = f"track {vid}" if i == 0 else None
                ax.plot(traj_to_plot[:, i, 0], traj_to_plot[:, i, 1], linestyle='-', marker='o', markersize=2, color='tab:green', label=label)

        ax.legend()
        ax.set_title("Cloth grid before/after VBD")
        plt.tight_layout()
        if args.savefig:
            plt.savefig(args.savefig)
        plt.show()

