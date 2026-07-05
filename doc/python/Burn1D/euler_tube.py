import argparse
import numpy as np
import matplotlib.pyplot as plt
from euler_tube_solver import CompressibleSolver
from euler_tube_plotting import setup_figure, set_fixed_limits, update_plot, save_snapshot

def main():
    parser = argparse.ArgumentParser(description="1D Euler Solver (Well-Balanced)")
    parser.add_argument("--cells", type=int, default=400, help="Number of cells")
    parser.add_argument("--length", type=float, default=5.0, help="Tube length (m)")
    parser.add_argument("--dt", type=float, default=None, help="Fixed time step")
    parser.add_argument("--prop_force", type=float, default=0.0, help="Propeller force (N)")
    parser.add_argument("--pulse", type=int, default=0, help="Enable pulse mode")
    parser.add_argument("--period_on", type=float, default=0.005, help="Pulse ON time")
    parser.add_argument("--period_off", type=float, default=0.05, help="Pulse OFF time")
    parser.add_argument("--u0", type=float, default=0.0, help="Initial velocity inside tube at t=0 (m/s)")
    parser.add_argument("--u_amb", type=float, default=0.0, help="Ambient/free-stream velocity at boundaries (m/s)")
    parser.add_argument("--verbosity", type=int, default=1, help="Verbosity level (0 = off)")
    parser.add_argument("--bc_mode", type=str, choices=["nonreflecting", "ambient", "original"], default="nonreflecting")
    parser.add_argument("--autoscale", type=int, default=1, help="Dynamic autoscaling of plots")
    parser.add_argument("--combustion", type=int, default=0, help="Enable combustion mode (species transport)")
    parser.add_argument("--jet-damping", type=float, default=0.0, help="Jet damping coefficient for valveless mechanism")
    parser.add_argument("--flux", type=str, choices=["hllc", "rusanov"], default="hllc", help="Flux scheme")
    parser.add_argument("--muscl", type=int, default=1, help="Enable MUSCL 2nd-order reconstruction")
    parser.add_argument("--inject-point", type=float, default=None, help="Fuel injection point (m)")
    parser.add_argument("--inject-period", type=float, default=0.020, help="Injection period (s)")
    parser.add_argument("--inject-duration", type=float, default=0.005, help="Injection duration (s)")
    parser.add_argument("--fuel-mass", type=float, default=3e-5, help="Fuel mass per cycle (kg)")
    parser.add_argument("--ignite-delay", type=float, default=0.001, help="Ignition delay (s)")
    parser.add_argument("--engine-start", type=float, default=0.8, help="Engine start x (m)")
    parser.add_argument("--engine-end", type=float, default=5.0, help="Engine end x (m)")
    parser.add_argument("--n-steps", type=int, default=500, help="Total steps (for --show 0)")
    parser.add_argument("--show", type=int, default=1, help="1=interactive animation, 0=save PNG")
    args = parser.parse_args()

    # Very smooth pulsejet: gentle asymmetric bump + continuously expanding exhaust
    d_res = 0.12       # left reservoir diameter
    d_tube = 0.07      # intake/throat diameter
    d_chamber = 0.16   # chamber diameter (2.3x tube — pronounced bump!)
    d_exit = 0.12      # exhaust exit diameter (gently expanding from throat)
    x_fine = np.arange(0, args.length + 0.005, 0.005)
    d_fine = np.full_like(x_fine, d_res)
    for i, xx in enumerate(x_fine):
        if xx < 0.40:
            d_fine[i] = d_res
        elif xx < 0.80:  # gentle contraction over 0.4m: reservoir → intake
            t = (xx - 0.40) / 0.40
            s = 0.5 - 0.5 * np.cos(np.pi * t)
            d_fine[i] = d_res * (1-s) + d_tube * s
        elif xx < 1.00:  # intake (flat 0.2m)
            d_fine[i] = d_tube
        elif xx < 1.40:  # gentle expansion over 0.4m: intake → chamber
            t = (xx - 1.00) / 0.40
            s = 0.5 - 0.5 * np.cos(np.pi * t)
            d_fine[i] = d_tube * (1-s) + d_chamber * s
        elif xx < 1.60:  # chamber (flat 0.2m)
            d_fine[i] = d_chamber
        elif xx < 2.00:  # gentle contraction over 0.4m: chamber → throat
            t = (xx - 1.60) / 0.40
            s = 0.5 - 0.5 * np.cos(np.pi * t)
            d_fine[i] = d_chamber * (1-s) + d_tube * s
        else:  # long expanding exhaust cone: throat → exit (3.0m long)
            t = (xx - 2.00) / 3.00
            d_fine[i] = d_tube * (1-t) + d_exit * t
    geom = [(x, d, 0.0) for x, d in zip(x_fine, d_fine)]
    if args.inject_point is None:
        args.inject_point = 1.50  # Chamber center

    solver = CompressibleSolver(nx=args.cells, length=args.length, geom_segments=geom,
                                u0=args.u0, u_amb=args.u_amb, p_init=101325.0, rho_init=1.225,
                                bc_mode=args.bc_mode, combustion=bool(args.combustion),
                                jet_damping=args.jet_damping, engine_start=args.engine_start,
                                engine_end=args.engine_end, inject_point=args.inject_point,
                                inject_period=args.inject_period, inject_duration=args.inject_duration,
                                fuel_mass=args.fuel_mass, ignite_delay=args.ignite_delay,
                                muscl=bool(args.muscl), flux_type=args.flux)
    rho0, u0, p0, T0, _ = solver.get_primitive()

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
                    rho_v, u_v, _, _, _ = solver.get_primitive()
                    total_mom = np.sum(solver.U[1, :] * solver.dx)
                    total_E = np.sum(solver.U[2, :] * solver.dx)
                    total_ke = np.sum(0.5 * rho_v * u_v**2 * solver.A * solver.dx)
                    thrust = solver.hist_thrust[-1] if solver.hist_thrust else 0.0
                    print(f"t={solver.t:.6f}s  mom={total_mom:.6e} kg·m/s  E={total_E:.6e} J  KE={total_ke:.6e} J  thrust={thrust:.3f} N")
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
