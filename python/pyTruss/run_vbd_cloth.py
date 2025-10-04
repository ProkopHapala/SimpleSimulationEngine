import argparse
import numpy as np
import matplotlib.pyplot as plt

from truss import Truss, solve_vbd_numpy
from truss_ocl import TrussOpenCLSolver
from plot_utils import plot_truss


'''
python run_vbd_cloth.py --nx 3 --ny 0 --niter 10 --serial 1 --cpu 1
'''

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="Run a simple Vertex Block Descent step on a cloth-like truss grid." )
    parser.add_argument("--nx",            type=int,   default=2, help="Number of cells in the X direction (grid has nx+1 points).")
    parser.add_argument("--ny",            type=int,   default=0, help="Number of cells in the Y direction (grid has ny+1 points).")
    parser.add_argument("--dt",            type=float, default=0.1, help="Simulation time step (seconds).")
    parser.add_argument("--mass",          type=float, default=1.0, help="Mass assigned to interior grid vertices.")
    parser.add_argument("--anchor-mass",   type=float, default=1.0e+6, help="Mass assigned to fixed corner vertices (handled by Truss.build_grid_2d).")
    parser.add_argument("--edge",          type=float, default=1.0, help="Rest length of grid edges (meters).")
    parser.add_argument("--stiffness",     type=float, default=1000000.0, help="Spring stiffness for axis-aligned edges.")
    parser.add_argument("--diag-stiffness",type=float, default=0.0, help="Spring stiffness for diagonal edges (if >0).")
    parser.add_argument("--gravity",       type=float, nargs=3, default=(0.0, -9.81, 0.0), help="Gravity vector (m/s^2).")
    parser.add_argument("--extra-accel",   type=float, nargs=3, default=(0.0, 0.0, 0.0), help="Additional acceleration applied to all vertices (m/s^2).")
    parser.add_argument("--niter",         type=int,   default=10, help="Number of VBD iterations to run.")
    parser.add_argument("--det-eps",       type=float, default=1e-6, help="Determinant threshold for skipping ill-conditioned Hessians.")
    parser.add_argument("--verbose",       type=int,   default=1, help="Print iteration progress from the solver.")
    parser.add_argument("--no-pin",        type=int,   default=0, help="Do not pin the default corner vertices; treat all vertices as free.")
    parser.add_argument("--no-plot",       type=int,   default=0, help="Skip plotting the initial and final grids.")
    parser.add_argument("--savefig",       type=str,   default="run_vbd_cloth.png", help="Path to save the comparison plot instead of displaying it.")
    parser.add_argument("--serial",        type=int,   default=1, help="Use OpenCL serial VBD kernel (simpler, for debugging).")
    parser.add_argument("--cpu",           type=int,   default=1, help="Run NumPy reference VBD solver.")
    parser.add_argument("--compare",       type=int,   default=0, help="Run both CPU and GPU solvers and report differences.")
    parser.add_argument("--chain",         type=int,   default=0, help="Treat the grid as a 1D chain (forces ny=0).")
    parser.add_argument("--anchor-mode",   type=str,   default="both", choices=["none", "left", "right", "both"], help="Override endpoint anchoring.")
    parser.add_argument("--track",         type=str,   default="all", help="Comma-separated vertex indices to plot trajectories for (CPU solver only).")
    args = parser.parse_args()

    solver = TrussOpenCLSolver()
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

    run_cpu = bool(args.cpu or args.compare)
    run_gpu = not args.compare or True

    cpu_positions = cpu_velocities = None
    gpu_positions = gpu_velocities = None
    cpu_traj = [] if run_cpu and track_indices is not None else None

    if track_indices is not None and not run_cpu:
        print("Warning: --track requires --cpu or --compare; trajectories will not be recorded for GPU-only runs.")

    if run_cpu:
        cpu_positions, cpu_velocities = solve_vbd_numpy(
            truss,
            dt=args.dt,
            gravity=total_accel.astype(np.float64),
            velocities=np.zeros_like(initial_points, dtype=np.float64),
            fixed_points=fixed_points,
            niter=args.niter,
            det_eps=args.det_eps,
            verbose=args.verbose,
            track_indices=track_indices,
            trajectory=cpu_traj,
        )

    if bool(args.compare) or not run_cpu:
        if bool(args.serial):
            gpu_positions, gpu_velocities = solver.solve_vbd_serial(
                truss,
                dt=args.dt,
                gravity=total_accel,
                fixed_points=fixed_points,
                niter=args.niter,
                det_eps=args.det_eps,
                bPrint=bool(args.verbose),
            )
        else:
            gpu_positions, gpu_velocities = solver.solve_vbd(
                truss,
                dt=args.dt,
                gravity=total_accel,
                fixed_points=fixed_points,
                niter=args.niter,
                det_eps=args.det_eps,
                bPrint=bool(args.verbose),
            )

    if bool(args.compare) and run_cpu and gpu_positions is not None:
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

    print(f"VBD finished after {args.niter} iterations")
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

        if cpu_traj and len(cpu_traj) > 0 and track_indices is not None:
            traj_arr = np.stack(cpu_traj, axis=0)  # [nstep, ntrack, 3]
            for i, vid in enumerate(track_indices):
                label = f"track {vid} (CPU)" if i == 0 else None
                ax.plot(traj_arr[:, i, 0], traj_arr[:, i, 1], linestyle='--', marker='o', markersize=2, color='tab:green', label=label)

        ax.legend()
        ax.set_title("Cloth grid before/after VBD")
        plt.tight_layout()
        if args.savefig:
            plt.savefig(args.savefig)
        plt.show()

