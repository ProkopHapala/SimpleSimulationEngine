import argparse
import numpy as np
import matplotlib.pyplot as plt

from truss import Truss
from truss_ocl import TrussOpenCLSolver
from plot_utils import plot_truss

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="Run a simple Vertex Block Descent step on a cloth-like truss grid." )
    parser.add_argument("--nx",            type=int,   default=10, help="Number of cells in the X direction (grid has nx+1 points).")
    parser.add_argument("--ny",            type=int,   default=10, help="Number of cells in the Y direction (grid has ny+1 points).")
    parser.add_argument("--dt",            type=float, default=0.033, help="Simulation time step (seconds).")
    parser.add_argument("--mass",          type=float, default=1.0, help="Mass assigned to interior grid vertices.")
    parser.add_argument("--anchor-mass",   type=float, default=100.0, help="Mass assigned to fixed corner vertices (handled by Truss.build_grid_2d).")
    parser.add_argument("--edge",          type=float, default=1.0, help="Rest length of grid edges (meters).")
    parser.add_argument("--stiffness",     type=float, default=200.0, help="Spring stiffness for axis-aligned edges.")
    parser.add_argument("--diag-stiffness",type=float, default=0.0, help="Spring stiffness for diagonal edges (if >0).")
    parser.add_argument("--gravity",       type=float, nargs=3, default=(0.0, -9.81, 0.0), help="Gravity vector (m/s^2).")
    parser.add_argument("--extra-accel",   type=float, nargs=3, default=(0.0, 0.0, 0.0), help="Additional acceleration applied to all vertices (m/s^2).")
    parser.add_argument("--niter",         type=int,   default=100, help="Number of VBD iterations to run.")
    parser.add_argument("--det-eps",       type=float, default=1e-6, help="Determinant threshold for skipping ill-conditioned Hessians.")
    parser.add_argument("--verbose",       type=int,   default=1, help="Print iteration progress from the solver.")
    parser.add_argument("--no-pin",        type=int,   default=0, help="Do not pin the default corner vertices; treat all vertices as free.")
    parser.add_argument("--no-plot",       type=int,   default=0, help="Skip plotting the initial and final grids.")
    parser.add_argument("--savefig",       type=str,   default="run_vbd_cloth.png", help="Path to save the comparison plot instead of displaying it.")
    args = parser.parse_args()

    solver = TrussOpenCLSolver()
    truss = Truss()

    truss.build_grid_2d(
        nx=args.nx,
        ny=args.ny,
        m=args.mass,
        m_end=args.anchor_mass,
        l=args.edge,
        k=args.stiffness,
        k_diag=args.diag_stiffness,
    )

    no_pin = bool(args.no_pin)
    fixed_points = None if no_pin else list(truss.fixed)

    gravity = np.asarray(args.gravity, dtype=np.float32)
    extra_accel = np.asarray(args.extra_accel, dtype=np.float32)
    total_accel = gravity + extra_accel

    initial_points = truss.points.copy()

    positions, velocities = solver.solve_vbd(
        truss,
        dt=args.dt,
        gravity=total_accel,
        fixed_points=fixed_points,
        niter=args.niter,
        det_eps=args.det_eps,
        bPrint=bool(args.verbose),
    )

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
        ax.legend()
        ax.set_title("Cloth grid before/after VBD")
        plt.tight_layout()
        if args.savefig:
            plt.savefig(args.savefig, dpi=150)
        plt.show()

