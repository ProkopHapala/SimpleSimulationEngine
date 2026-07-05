import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from tube_plotting import draw_tube_outline, draw_tube_outline_interp, make_element_segments

def setup_figure_engine(solver, args):
    fig = plt.figure(figsize=(14, 12))
    gs = fig.add_gridspec(5, 2)
    ax_geo = fig.add_subplot(gs[0, :])
    ax_geo.set_title(f"1D Pulsejet: Composition & Combustion (t=0.00s)")
    ax_geo.set_xlim(-0.2, solver.length + 0.2)
    ax_geo.set_ylim(-0.15, 0.15)
    ax_geo.set_aspect('equal')
    gx = np.linspace(0, solver.length, 300)
    gr = np.array([solver.get_radius(x) for x in gx])
    ax_geo.plot(gx, gr, 'k', lw=2)
    ax_geo.plot(gx, -gr, 'k', lw=2)
    lc_elements = LineCollection([], linewidths=0.5)
    ax_geo.add_collection(lc_elements)
    ax_spec = fig.add_subplot(gs[1, :])
    ax_spec.set_title("Chemical Composition (Mass Fraction)")
    ax_spec.set_xlim(-0.2, solver.length + 0.2)
    ax_spec.set_ylim(0, 1.0)
    ax_spec.grid(True, alpha=0.3)
    l_fuel, = ax_spec.plot([], [], 'orange', label='Fuel', lw=2)
    l_n2, = ax_spec.plot([], [], 'gray', label='N2', lw=1.5)
    l_o2, = ax_spec.plot([], [], 'b', label='O2', lw=1.5)
    l_prod, = ax_spec.plot([], [], 'g', label='Products (CO2+H2O)', lw=1.5)
    ax_spec.legend(loc='upper right')
    ax_thermo = fig.add_subplot(gs[2, :])
    ax_thermo.set_xlim(-0.2, solver.length + 0.2)
    ax_thermo.set_title("Thermodynamics")
    line_p, = ax_thermo.plot([], [], 'r-', lw=1, label='Pressure (atm)')
    line_t, = ax_thermo.plot([], [], 'm-', lw=1, label='Temp (x300K)')
    line_rho, = ax_thermo.plot([], [], 'b-', lw=1, label='Density (kg/m3)')
    ax_thermo.legend(loc='upper right')
    ax_thermo.grid(True, alpha=0.3)
    ax_vel = fig.add_subplot(gs[3, :])
    ax_vel.set_xlim(-0.2, solver.length + 0.2)
    ax_vel.set_ylim(-400, 400)
    ax_vel.set_ylabel("Velocity (m/s)")
    line_v, = ax_vel.plot([], [], 'k-', lw=1)
    ax_vel.grid(True, alpha=0.3)
    ax_mom = fig.add_subplot(gs[4, :])
    ax_mom.set_title("Total Momentum & Thrust")
    ax_mom.set_xlabel("Time (s)")
    line_mom, = ax_mom.plot([], [], 'k-', lw=1.5, label='Momentum')
    ax_thrust = ax_mom.twinx()
    line_thrust, = ax_thrust.plot([], [], 'r-', lw=1, label='Thrust (N)')
    ax_mom.legend(loc='upper left')
    ax_thrust.legend(loc='upper right')
    ax_mom.grid(True)
    plt.tight_layout()
    return fig, ax_geo, ax_spec, ax_thermo, ax_vel, ax_mom, lc_elements, l_fuel, l_n2, l_o2, l_prod, line_p, line_t, line_rho, line_v, line_mom, ax_thrust, line_thrust

def update_plot_engine(solver, args, fig, ax_geo, ax_spec, ax_thermo, ax_vel, ax_mom,
                       lc_elements, l_fuel, l_n2, l_o2, l_prod, line_p, line_t, line_rho, line_v, line_mom,
                       ax_thrust=None, line_thrust=None):
    centers, rho, p, t, Y, vols, m_t = solver.get_state()
    if solver.sim_time >= solver.inject_start:
        phase = (solver.sim_time - solver.inject_start) % solver.inject_period
        t_to_inj = 0.0 if phase < solver.inject_duration else (solver.inject_period - phase)
        t_to_ign = max(0.0, solver.inject_duration + solver.ignite_delay - phase)
    else:
        t_to_inj = solver.inject_start - solver.sim_time
        t_to_ign = t_to_inj + solver.inject_duration + solver.ignite_delay
    segs = []
    colors = []
    norm_t = np.clip((t - 300) / 400, 0, 1)
    radii = solver.get_radius(solver.nodes_x)
    for i in range(len(solver.nodes_x)):
        x = solver.nodes_x[i]
        r = radii[i]
        segs.append([[x, -r], [x, r]])
        if i < len(t): val = norm_t[i]
        elif i > 0: val = norm_t[i - 1]
        else: val = 0
        colors.append((val, 0, 1 - val, 1))
    lc_elements.set_segments(segs)
    lc_elements.set_color(colors)
    ax_geo.set_title(f"t={solver.sim_time * 1000:.1f}ms | next_inj={t_to_inj * 1000:.1f}ms | next_ign={t_to_ign * 1000:.1f}ms")
    if Y is not None:
        l_fuel.set_data(centers, Y[:, 2])
        l_o2.set_data(centers, Y[:, 1])
        l_prod.set_data(centers, Y[:, 3] + Y[:, 4])
    line_p.set_data(centers, p / P_ATM)
    line_t.set_data(centers, t / 300.0)
    line_rho.set_data(centers, rho)
    p_norm = p / P_ATM
    t_norm = t / 300.0
    all_vals = np.concatenate([p_norm, t_norm, rho])
    all_vals = all_vals[np.isfinite(all_vals)]
    if len(all_vals) > 0:
        y_min, y_max = all_vals.min(), all_vals.max()
        pad = 0.1 * max(1e-3, y_max - y_min)
        ax_thermo.set_ylim(y_min - pad, y_max + pad)
    v_elem = 0.5 * (solver.nodes_v[:-1] + solver.nodes_v[1:])
    line_v.set_data(centers, v_elem)
    if solver.hist_time:
        line_mom.set_data(solver.hist_time, solver.hist_momentum)
        ax_mom.set_xlim(0, max(solver.hist_time[-1], 0.05))
        vals = np.array(solver.hist_momentum)
        vals = vals[np.isfinite(vals)]
        if len(vals) > 5:
            ax_mom.set_ylim(vals.min(), vals.max())
    if ax_thrust is not None and line_thrust is not None and solver.hist_time:
        line_thrust.set_data(solver.hist_time, solver.hist_thrust)
        thrust_vals = np.array(solver.hist_thrust)
        thrust_vals = thrust_vals[np.isfinite(thrust_vals)]
        if len(thrust_vals) > 5:
            ax_thrust.set_ylim(thrust_vals.min(), thrust_vals.max())

def setup_figure_inertia(solver):
    fig = plt.figure(figsize=(12, 10))
    gs = fig.add_gridspec(4, 1)
    ax_vis = fig.add_subplot(gs[0])
    ax_vis.set_aspect('equal')
    ax_vis.set_xlim(-0.2, solver.length + 0.2)
    ax_vis.set_ylim(-0.15, 0.15)
    gx = np.linspace(-0.5, solver.length + 0.5, 300)
    gr = np.array([solver.get_radius(x) for x in gx])
    ax_vis.plot(gx, gr, 'k', lw=1)
    ax_vis.plot(gx, -gr, 'k', lw=1)
    lc = LineCollection([], linewidths=0.5)
    ax_vis.add_collection(lc)
    ax_vel = fig.add_subplot(gs[1])
    ax_vel.set_ylabel("Velocity (m/s)")
    line_v, = ax_vel.plot([], [], 'g.-', lw=0.5, markersize=2)
    ax_vel.grid(True)
    ax_p = fig.add_subplot(gs[2])
    ax_p.set_ylabel("Pressure (atm)")
    line_p, = ax_p.plot([], [], 'r-', lw=1)
    ax_p.axhline(1.0, color='gray', ls=':')
    ax_p.grid(True)
    ax_hist = fig.add_subplot(gs[3])
    ax_hist.set_ylabel("Total Momentum (kg m/s)")
    line_mom, = ax_hist.plot([], [], 'k-', label='Momentum')
    ax_hist.legend()
    ax_hist.grid(True)
    plt.tight_layout()
    return fig, ax_vis, ax_vel, ax_p, ax_hist, lc, line_v, line_p, line_mom

def update_plot_inertia(solver, fig, ax_vis, ax_vel, ax_p, ax_hist, lc, line_v, line_p, line_mom):
    centers, rho, p, t, Y, vols, m_t = solver.get_state()
    radii = solver.get_radius(solver.nodes_x)
    segs = []
    for i in range(len(solver.nodes_x)):
        segs.append([[solver.nodes_x[i], -radii[i]], [solver.nodes_x[i], radii[i]]])
    lc.set_segments(segs)
    line_v.set_data(centers, 0.5 * (solver.nodes_v[:-1] + solver.nodes_v[1:]))
    if solver.mode == "INERTIA":
        ax_vel.set_ylim(0, 300)
    line_p.set_data(centers, p / P_ATM)
    if solver.hist_time:
        line_mom.set_data(solver.hist_time, solver.hist_momentum)
        ax_hist.set_xlim(0, max(solver.hist_time[-1], 0.1))
        ax_hist.set_ylim(min(solver.hist_momentum) * 1.1, max(solver.hist_momentum) * 1.1)

def save_snapshot(solver, filename='pulsejet.png'):
    centers, rho, p, t, Y, vols, m_t = solver.get_state()
    fig, axes = plt.subplots(3, 1, figsize=(10, 10))
    gx = np.linspace(0, solver.length, 300)
    gr = np.array([solver.get_radius(x) for x in gx])
    axes[0].plot(gx, gr, 'k', lw=2)
    axes[0].plot(gx, -gr, 'k', lw=2)
    radii = solver.get_radius(solver.nodes_x)
    for i in range(0, len(solver.nodes_x), max(1, len(solver.nodes_x) // 100)):
        axes[0].plot([solver.nodes_x[i], solver.nodes_x[i]], [-radii[i], radii[i]], 'b', lw=0.5)
    axes[0].set_xlim(-0.2, solver.length + 0.2)
    axes[0].set_ylim(-0.15, 0.15)
    axes[0].set_aspect('equal')
    axes[0].set_title(f"Pulsejet [{solver.mode}] | t={solver.sim_time * 1000:.1f}ms")
    axes[1].plot(centers, p / P_ATM, 'r-', label='Pressure (atm)')
    axes[1].plot(centers, 0.5 * (solver.nodes_v[:-1] + solver.nodes_v[1:]), 'g-', label='Velocity')
    axes[1].set_ylabel("p [atm], v [m/s]")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    if solver.hist_time:
        axes[2].plot(solver.hist_time, solver.hist_momentum, 'k-', label='Momentum')
        axes[2].set_xlabel("Time (s)")
        axes[2].set_ylabel("Momentum")
        axes[2].grid(True)
    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    print(f"Saved {filename}")
    plt.close(fig)

from tube_common import P_ATM
