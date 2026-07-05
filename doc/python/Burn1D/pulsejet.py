import argparse
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from lagrangian_tube_solver import PulsejetSolver
from pulsejet_plotting import (setup_figure_engine, update_plot_engine,
                               setup_figure_inertia, update_plot_inertia, save_snapshot)

def main():
    parser = argparse.ArgumentParser(description="1D Pulsejet / Inertia test solver")
    parser.add_argument("--mode", choices=["ENGINE", "INERTIA"], default="ENGINE", help="Engine mode")
    parser.add_argument("--cfl", type=float, default=0.4, help="CFL number")
    parser.add_argument("--dt-max", type=float, default=5e-5, help="Max timestep [s]")
    parser.add_argument("--steps-per-frame", type=int, default=20, help="Physics steps per frame")
    parser.add_argument("--target-mass", type=float, default=200e-5, help="Target mass per element [kg]")
    parser.add_argument("--friction", type=float, default=0.1, help="Wall friction coefficient")
    parser.add_argument("--viscosity", type=float, default=0.1, help="Artificial viscosity")
    parser.add_argument("--thermal-cond", type=float, default=0.1, help="Thermal mixing coefficient")
    parser.add_argument("--species-diff", type=float, default=0.1, help="Species mixing coefficient")
    parser.add_argument("--relax-time", type=float, default=0.003, help="Pre-injection relaxation time [s]")
    parser.add_argument("--inject-period", type=float, default=0.020, help="Injection period [s]")
    parser.add_argument("--inject-start", type=float, default=0.005, help="First injection start [s]")
    parser.add_argument("--inject-duration", type=float, default=0.005, help="Injection duration [s]")
    parser.add_argument("--inject-point", type=float, default=0.0, help="Injection point (abs x in m)")
    parser.add_argument("--fuel-mass", type=float, default=3e-5, help="Fuel mass per cycle [kg]")
    parser.add_argument("--ignite-delay", type=float, default=0.001, help="Ignition delay [s]")
    parser.add_argument("--burn-time", type=float, default=0.02, help="Burn time (INERTIA mode) [s]")
    parser.add_argument("--burn-mult", type=float, default=4.0, help="Burn multiplier (INERTIA mode)")
    parser.add_argument("--jet-damping", type=float, default=5.0, help="Jet damping coefficient for valveless mechanism")
    parser.add_argument("--n-steps", type=int, default=500, help="Total steps (for --show 0)")
    parser.add_argument("--show", type=int, default=1, help="1=interactive animation, 0=save PNG")
    args = parser.parse_args()

    solver = PulsejetSolver(
        mode=args.mode, target_mass=args.target_mass, cfl=args.cfl, dt_max=args.dt_max,
        friction=args.friction, viscosity=args.viscosity, thermal_cond=args.thermal_cond,
        species_diff=args.species_diff, relax_time=args.relax_time,
        inject_period=args.inject_period, inject_start=args.inject_start,
        inject_duration=args.inject_duration, inject_point=args.inject_point,
        fuel_mass=args.fuel_mass, ignite_delay=args.ignite_delay,
        burn_time=args.burn_time, burn_mult=args.burn_mult, jet_damping=args.jet_damping)

    if args.show:
        paused = [False]
        if args.mode == "ENGINE":
            fig, ax_geo, ax_spec, ax_thermo, ax_vel, ax_mom, lc_elements, l_fuel, l_n2, l_o2, l_prod, line_p, line_t, line_rho, line_v, line_mom, ax_thrust, line_thrust = setup_figure_engine(solver, args)
            def update(frame):
                if not paused[0]:
                    for _ in range(args.steps_per_frame):
                        solver.step()
                update_plot_engine(solver, args, fig, ax_geo, ax_spec, ax_thermo, ax_vel, ax_mom,
                                   lc_elements, l_fuel, l_n2, l_o2, l_prod, line_p, line_t, line_rho, line_v, line_mom,
                                   ax_thrust, line_thrust)
                return lc_elements, line_p, line_v, line_mom
        else:
            fig, ax_vis, ax_vel, ax_p, ax_hist, lc, line_v, line_p, line_mom = setup_figure_inertia(solver)
            def update(frame):
                if not paused[0]:
                    for _ in range(args.steps_per_frame):
                        solver.step()
                update_plot_inertia(solver, fig, ax_vis, ax_vel, ax_p, ax_hist, lc, line_v, line_p, line_mom)
                return lc, line_v, line_p, line_mom

        def on_key(event):
            if event.key == 'i':
                solver.manual_inject()
                print(f"[INJECT] {solver.fuel_mass:.2e} kg at x={solver.inject_point:.3f}m  t={solver.sim_time:.5f}s")
            elif event.key == 's':
                solver.manual_ignite()
                print(f"[SPARK] Ignited chamber  t={solver.sim_time:.5f}s")
            elif event.key == 'a':
                solver.auto_inject = not solver.auto_inject
                print(f"[AUTO] {'ON' if solver.auto_inject else 'OFF'}")
            elif event.key == ' ':
                paused[0] = not paused[0]
                print(f"[{'PAUSED' if paused[0] else 'RESUMED'}]")
            elif event.key == 'r':
                solver.manual_inject()
                solver.manual_ignite()
                print(f"[RAM] Inject+Ignite  t={solver.sim_time:.5f}s")
            elif event.key == '?':
                print("\n--- Pulsejet Controls ---")
                print("  i     = inject fuel")
                print("  s     = spark (ignite)")
                print("  r     = inject + ignite (ram)")
                print("  a     = toggle auto injection/ignition")
                print("  space = pause/resume")
                print("  ?     = this help")
                print("------------------------\n")

        fig.canvas.mpl_connect('key_press_event', on_key)
        print("Interactive mode. Keys: i=inject  s=spark  r=ram  a=auto  space=pause  ?=help")
        crashed = [False]
        def safe_update(frame):
            if crashed[0]:
                return lc_elements if args.mode == "ENGINE" else lc
            try:
                return update(frame)
            except Exception as e:
                crashed[0] = True
                print(f"\n*** SIMULATION CRASHED at t={solver.sim_time:.6f}s ***")
                print(f"    {type(e).__name__}: {e}")
                print("    Close the window to exit.")
                ani.event_source.stop()
                return lc_elements if args.mode == "ENGINE" else lc
        ani = animation.FuncAnimation(fig, safe_update, interval=10, blit=False, cache_frame_data=False)
        plt.show()
    else:
        for _ in range(args.n_steps):
            solver.step()
        save_snapshot(solver)

if __name__ == "__main__":
    main()
