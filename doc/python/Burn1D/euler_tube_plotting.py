import numpy as np
import matplotlib.pyplot as plt
from tube_plotting import draw_tube_outline

def setup_figure(solver, prop_config, args):
    plt.ion()
    fig = plt.figure(figsize=(12, 8))
    gs = fig.add_gridspec(3, 2)
    ax_geom = fig.add_subplot(gs[0, 0])
    ax_p = fig.add_subplot(gs[1, 0], sharex=ax_geom)
    ax_u = fig.add_subplot(gs[2, 0], sharex=ax_geom)
    ax_rho = fig.add_subplot(gs[1, 1], sharex=ax_geom)
    ax_hist = fig.add_subplot(gs[2, 1])
    r = np.sqrt(solver.A / np.pi)
    ax_geom.plot(solver.x, r, 'k-', lw=2)
    ax_geom.plot(solver.x, -r, 'k-', lw=2)
    ax_geom.fill_between(solver.x, r, -r, color='gray', alpha=0.3)
    ax_geom.set_ylabel("Radius (m)")
    ax_geom.set_title("Tube Geometry")
    if args.prop_force != 0:
        ax_geom.axvspan(prop_config['x_start'], prop_config['x_end'], color='red', alpha=0.3, label='Propeller')
        ax_geom.legend()
    line_p, = ax_p.plot([], [], 'b-')
    line_u, = ax_u.plot([], [], 'g-')
    line_rho, = ax_rho.plot([], [], 'r-')
    line_mass, = ax_hist.plot([], [], 'k-', label='Total Mass')
    ax_hist2 = ax_hist.twinx()
    line_eng, = ax_hist2.plot([], [], 'r--', label='Total Energy')
    ax_p.set_ylabel("Pressure (Pa)")
    ax_u.set_ylabel("Velocity (m/s)")
    ax_rho.set_ylabel("Density (kg/m3)")
    ax_u.set_xlabel("Position (m)")
    ax_hist.set_xlabel("Time (s)")
    ax_hist.set_ylabel("Mass (kg)")
    ax_hist2.set_ylabel("Energy (J)")
    ax_p.axhline(solver.p_amb, color='gray', linestyle=':', linewidth=1, label='Ambient p')
    ax_rho.axhline(solver.rho_amb, color='gray', linestyle=':', linewidth=1, label='Ambient rho')
    ax_u.axhline(0.0, color='gray', linestyle=':', linewidth=1, label='u=0')
    lines = [line_mass, line_eng]
    labels = [l.get_label() for l in lines]
    ax_hist.legend(lines, labels, loc='center right')
    return fig, ax_geom, ax_p, ax_u, ax_rho, ax_hist, ax_hist2, line_p, line_u, line_rho, line_mass, line_eng

def set_fixed_limits(ax_p, ax_rho, ax_u, p0, rho0, u0):
    def padded_range(arr, frac=0.05, min_span=1e-3):
        a_min = float(np.min(arr)); a_max = float(np.max(arr))
        span = a_max - a_min
        pad = max(frac * max(abs(a_min), abs(a_max), 1.0), min_span)
        if span < min_span: a_min -= pad; a_max += pad
        else: a_min -= pad; a_max += pad
        return a_min, a_max
    p_lim = padded_range(p0, frac=0.01, min_span=100.0)
    rho_lim = padded_range(rho0, frac=0.01, min_span=0.1)
    u_span = max(np.max(np.abs(u0)), 0.5)
    u_lim = (-1.2 * u_span, 1.2 * u_span)
    ax_p.set_ylim(*p_lim)
    ax_rho.set_ylim(*rho_lim)
    ax_u.set_ylim(*u_lim)

def update_plot(solver, prop_config, args, fig, ax_geom, ax_p, ax_u, ax_rho, ax_hist, ax_hist2,
                line_p, line_u, line_rho, line_mass, line_eng, history_t, history_mass, history_energy, dt_sim):
    rho, u, p, T, Y_fuel = solver.get_primitive()
    line_p.set_data(solver.x, p)
    line_u.set_data(solver.x, u)
    line_rho.set_data(solver.x, rho)
    if args.autoscale:
        ax_p.set_ylim(np.min(p) * 0.95, np.max(p) * 1.05)
        ax_u.set_ylim(np.min(u) - 10, np.max(u) + 10)
        ax_rho.set_ylim(np.min(rho) * 0.95, np.max(rho) * 1.05)
    history_t.append(solver.t)
    history_mass.append(np.sum(solver.U[0, :] * solver.dx))
    history_energy.append(np.sum(solver.U[2, :] * solver.dx))
    if len(history_t) > 200:
        history_t.pop(0); history_mass.pop(0); history_energy.pop(0)
    line_mass.set_data(history_t, history_mass)
    line_eng.set_data(history_t, history_energy)
    ax_hist.set_xlim(min(history_t), max(history_t) + dt_sim)
    ax_hist.set_ylim(min(history_mass) * 0.999, max(history_mass) * 1.001)
    ax_hist2.set_ylim(min(history_energy) * 0.95, max(history_energy) * 1.05)
    ax_geom.set_title(f"t={solver.t:.4f}s | Propeller: {'ON' if (solver.t % (prop_config['period_on'] + prop_config['period_off']) < prop_config['period_on']) and args.pulse else 'OFF'}")
    plt.pause(0.001)

def save_snapshot(solver, filename='euler_tube.png'):
    rho, u, p, T, Y_fuel = solver.get_primitive()
    fig, axes = plt.subplots(3, 1, figsize=(10, 10))
    r = np.sqrt(solver.A / np.pi)
    axes[0].plot(solver.x, r, 'k-', lw=2)
    axes[0].plot(solver.x, -r, 'k-', lw=2)
    axes[0].fill_between(solver.x, r, -r, color='gray', alpha=0.3)
    axes[0].set_ylabel("Radius (m)")
    axes[0].set_title(f"Tube Geometry | t={solver.t:.4f}s")
    axes[1].plot(solver.x, p, 'b-', label='Pressure')
    axes[1].plot(solver.x, u, 'g-', label='Velocity')
    axes[1].set_ylabel("p [Pa], u [m/s]")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    axes[2].plot(solver.x, rho, 'r-', label='Density')
    axes[2].set_ylabel("Density (kg/m3)")
    axes[2].set_xlabel("Position (m)")
    axes[2].grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    print(f"Saved {filename}")
    plt.close(fig)
