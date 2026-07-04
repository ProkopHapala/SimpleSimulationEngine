import argparse
import matplotlib.pyplot as plt
from lagrangian_tube_solver import WaveTubeSolver
from wave_tube_plotting import setup_figure, update_plot, save_snapshot

def main():
    parser = argparse.ArgumentParser(description="1D Lagrangian wave propagation in variable geometry")
    parser.add_argument("--geo", choices=["nozzle", "cone"], default="nozzle", help="Geometry type")
    parser.add_argument("--length", type=float, default=2.0, help="Tube length (m)")
    parser.add_argument("--n-cells", type=int, default=1000, help="Number of elements")
    parser.add_argument("--kick-vel", type=float, default=50.0, help="Strength of the kick (m/s)")
    parser.add_argument("--dt", type=float, default=2e-6, help="Time step (s)")
    parser.add_argument("--steps-per-frame", type=int, default=50, help="Sim steps per frame")
    parser.add_argument("--n-steps", type=int, default=500, help="Total steps (for --show 0)")
    parser.add_argument("--show", type=int, default=1, help="1=interactive animation, 0=save PNG")
    args = parser.parse_args()

    solver = WaveTubeSolver(geo=args.geo, length=args.length, n_cells=args.n_cells,
                            kick_vel=args.kick_vel, dt=args.dt)

    if args.show:
        fig, ax_geo, ax_data, ax_p, ax_hist, ax_mom_r, lc_lines, line_v, line_p, line_mom, line_E = setup_figure(solver, args)
        solver.step()
        E0 = solver.hist_E[0]
        ax_mom_r.set_ylim(E0 * 0.9995, E0 * 1.0005)
        import matplotlib.animation as animation
        def update(frame):
            for _ in range(args.steps_per_frame):
                p = solver.step()
            update_plot(solver, p, lc_lines, line_v, line_p, line_mom, line_E, ax_hist, ax_mom_r)
            return lc_lines, line_v, line_p, line_mom, line_E
        ani = animation.FuncAnimation(fig, update, interval=20, blit=False)
        plt.show()
    else:
        for _ in range(args.n_steps):
            solver.step()
        save_snapshot(solver)

if __name__ == "__main__":
    main()
