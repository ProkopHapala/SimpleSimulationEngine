import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from tube_plotting import draw_tube_outline, make_element_segments, draw_element_lines

def setup_figure(solver, args):
    fig, (ax_geo, ax_data, ax_hist) = plt.subplots(3, 1, figsize=(10, 12))
    ax_geo.set_title(f"Kick Test: Geometry = {args.geo.upper()}")
    ax_geo.set_aspect('equal')
    draw_tube_outline(ax_geo, solver.get_radius, solver.length)
    lc_lines = LineCollection([], lw=0.5, colors='b')
    ax_geo.add_collection(lc_lines)
    ax_geo.set_xlim(0, solver.length)
    ax_geo.set_ylim(-0.15, 0.15)
    line_v, = ax_data.plot([], [], 'g-', label='Velocity')
    ax_data.set_ylabel("Velocity (m/s)")
    ax_data.set_ylim(-args.kick_vel * 1.5, args.kick_vel * 1.5)
    ax_data.grid(True)
    ax_data.legend(loc='upper left')
    ax_p = ax_data.twinx()
    line_p, = ax_p.plot([], [], 'r-', label='Pressure')
    ax_p.set_ylabel("Pressure (atm)", color='r')
    ax_p.set_ylim(0.8, 1.2)
    ax_p.legend(loc='upper right')
    line_mom, = ax_hist.plot([], [], 'b-', alpha=0.5, label='Momentum (Gas)')
    ax_mom_r = ax_hist.twinx()
    line_E, = ax_mom_r.plot([], [], 'k-', lw=2, label='Total Energy (MUST BE FLAT)')
    ax_hist.set_xlabel("Time (s)")
    ax_hist.set_ylabel("Momentum")
    ax_mom_r.set_ylabel("Energy")
    ax_hist.grid(True)
    ax_hist.legend(loc='upper left')
    ax_mom_r.legend(loc='upper right')
    return fig, ax_geo, ax_data, ax_p, ax_hist, ax_mom_r, lc_lines, line_v, line_p, line_mom, line_E

def update_plot(solver, p, lc_lines, line_v, line_p, line_mom, line_E, ax_hist, ax_mom_r):
    centers = solver.get_centers()
    radii = solver.get_radius(solver.nodes_x)
    segs = []
    step_vis = max(1, len(solver.nodes_x) // 200)
    for i in range(0, len(solver.nodes_x), step_vis):
        segs.append([[solver.nodes_x[i], -radii[i]], [solver.nodes_x[i], radii[i]]])
    lc_lines.set_segments(segs)
    line_v.set_data(centers, solver.get_vel_elements())
    line_p.set_data(centers, p / 101325.0)
    if len(solver.hist_t) > 0:
        line_mom.set_data(solver.hist_t, solver.hist_mom)
        line_E.set_data(solver.hist_t, solver.hist_E)
        ax_hist.set_xlim(0, max(solver.hist_t[-1], 0.01))
        m_max = max(np.max(np.abs(solver.hist_mom)), 1e-3)
        ax_hist.set_ylim(-m_max * 1.2, m_max * 1.2)

def save_snapshot(solver, filename='wave_tube.png'):
    fig, (ax_geo, ax_data, ax_hist) = plt.subplots(3, 1, figsize=(10, 12))
    draw_tube_outline(ax_geo, solver.get_radius, solver.length)
    draw_element_lines(ax_geo, solver.nodes_x, solver.get_radius, max(1, len(solver.nodes_x) // 200))
    ax_geo.set_xlim(0, solver.length)
    ax_geo.set_ylim(-0.15, 0.15)
    ax_geo.set_aspect('equal')
    ax_geo.set_title(f"Wave Tube: {solver.geo} | t={solver.sim_t:.4e}s")
    centers = solver.get_centers()
    ax_data.plot(centers, solver.get_vel_elements(), 'g-', label='Velocity')
    ax_data.set_ylabel("Velocity (m/s)")
    ax_data.grid(True)
    ax_data.legend()
    p = 101325.0 * (solver.cell_mass / ((np.diff(solver.nodes_x) / 3.0) * (
        solver.get_radius(solver.nodes_x[:-1])**2 * np.pi + solver.get_radius(solver.nodes_x[1:])**2 * np.pi +
        np.sqrt(solver.get_radius(solver.nodes_x[:-1])**2 * np.pi * solver.get_radius(solver.nodes_x[1:])**2 * np.pi))) / 1.225)**1.4
    ax_p = ax_data.twinx()
    ax_p.plot(centers, p / 101325.0, 'r-', label='Pressure')
    ax_p.set_ylabel("Pressure (atm)", color='r')
    if solver.hist_t:
        ax_hist.plot(solver.hist_t, solver.hist_mom, 'b-', alpha=0.5, label='Momentum')
        ax_hist.set_xlabel("Time (s)")
        ax_hist.set_ylabel("Momentum")
        ax_hist.grid(True)
        ax_hist.legend()
    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    print(f"Saved {filename}")
    plt.close(fig)
