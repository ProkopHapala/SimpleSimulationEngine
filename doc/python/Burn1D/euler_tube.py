import argparse
import numpy as np
import matplotlib.pyplot as plt
from euler_tube_solver import CompressibleSolver
from euler_tube_plotting import setup_figure, set_fixed_limits, update_plot, save_snapshot

def main():
    parser = argparse.ArgumentParser(description="1D Euler Solver (Well-Balanced)")
    parser.add_argument("--cells", type=int, default=100, help="Number of cells")
    parser.add_argument("--length", type=float, default=2.0, help="Tube length (m)")
    parser.add_argument("--dt", type=float, default=None, help="Fixed time step")
    parser.add_argument("--prop_force", type=float, default=10.0, help="Propeller force (N)")
    parser.add_argument("--pulse", type=int, default=0, help="Enable pulse mode")
    parser.add_argument("--period_on", type=float, default=0.005, help="Pulse ON time")
    parser.add_argument("--period_off", type=float, default=0.05, help="Pulse OFF time")
    parser.add_argument("--u0", type=float, default=1.0, help="Initial velocity inside tube at t=0 (m/s)")
    parser.add_argument("--u_amb", type=float, default=10.0, help="Ambient/free-stream velocity at boundaries (m/s)")
    parser.add_argument("--verbosity", type=int, default=1, help="Verbosity level (0 = off)")
    parser.add_argument("--bc_mode", type=str, choices=["ambient", "original"], default="ambient")
    parser.add_argument("--autoscale", type=int, default=1, help="Dynamic autoscaling of plots")
    parser.add_argument("--n-steps", type=int, default=500, help="Total steps (for --show 0)")
    parser.add_argument("--show", type=int, default=1, help="1=interactive animation, 0=save PNG")
    args = parser.parse_args()

    geom = [
        (0.0, 0.1, 0.0),
        (0.3 * args.length, 0.1, 0.0),
        (0.5 * args.length, 0.05, 0.0),
        (0.7 * args.length, 0.1, 0.0),
        (args.length, 0.1, 0.0)
    ]

    solver = CompressibleSolver(nx=args.cells, length=args.length, geom_segments=geom,
                                u0=args.u0, u_amb=args.u_amb, p_init=101325.0, rho_init=1.225,
                                bc_mode=args.bc_mode)
    rho0, u0, p0, T0 = solver.get_primitive()

    prop_config = {
        'force': args.prop_force,
        'x_start': 0.1 * args.length,
        'x_end': 0.2 * args.length,
        'pulse_mode': args.pulse,
        'period_on': args.period_on,
        'period_off': args.period_off
    }

    if args.show:
        fig, ax_geom, ax_p, ax_u, ax_rho, ax_hist, ax_hist2, line_p, line_u, line_rho, line_mass, line_eng = setup_figure(solver, prop_config, args)
        if not args.autoscale:
            set_fixed_limits(ax_p, ax_rho, ax_u, p0, rho0, u0)
        history_t, history_mass, history_energy = [], [], []
        frame = 0
        try:
            while True:
                dt_sim = solver.step(dt_fixed=args.dt, prop_params=prop_config)
                if args.u0 == 0 and args.prop_force == 0 and frame % 50 == 0:
                    v_max = np.max(np.abs(solver.get_primitive()[1]))
                    if v_max > 1e-2:
                        print(f"WARNING: Drift detected! Max Vel = {v_max:.5f} m/s")
                if args.verbosity > 0:
                    rho_v, u_v, _, _ = solver.get_primitive()
                    total_mom = np.sum(solver.U[1, :] * solver.dx)
                    total_E = np.sum(solver.U[2, :] * solver.dx)
                    total_ke = np.sum(0.5 * rho_v * u_v**2 * solver.A * solver.dx)
                    print(f"t={solver.t:.6f}s  mom={total_mom:.6e} kg·m/s  E={total_E:.6e} J  KE={total_ke:.6e} J")
                if frame % 5 == 0:
                    update_plot(solver, prop_config, args, fig, ax_geom, ax_p, ax_u, ax_rho,
                                ax_hist, ax_hist2, line_p, line_u, line_rho, line_mass, line_eng,
                                history_t, history_mass, history_energy, dt_sim)
                    if not plt.fignum_exists(fig.number):
                        break
                frame += 1
        except KeyboardInterrupt:
            print("\nSimulation stopped by user.")
        finally:
            solver.print_profiles("final")
    else:
        for _ in range(args.n_steps):
            solver.step(dt_fixed=args.dt, prop_params=prop_config)
        solver.print_profiles("final")
        save_snapshot(solver)

if __name__ == "__main__":
    main()
