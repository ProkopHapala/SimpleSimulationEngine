import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

def draw_tube_outline(ax, get_radius, length, n_samples=500, lw=2, color='k'):
    """Draw upper and lower tube wall curves."""
    gx = np.linspace(0, length, n_samples)
    gr = np.array([get_radius(x) for x in gx])
    ax.plot(gx, gr, color, lw=lw)
    ax.plot(gx, -gr, color, lw=lw)
    return gx, gr

def draw_tube_outline_interp(ax, get_radius, x_min, x_max, n_samples=500, lw=2, color='k'):
    """Draw tube outline for arbitrary x range (for solvers with negative coords)."""
    gx = np.linspace(x_min, x_max, n_samples)
    gr = np.array([get_radius(x) for x in gx])
    ax.plot(gx, gr, color, lw=lw)
    ax.plot(gx, -gr, color, lw=lw)
    return gx, gr

def make_element_segments(nodes_x, get_radius, step=1):
    """Return list of [[x,-r],[x,r]] segments for LineCollection."""
    radii = get_radius(nodes_x)
    segs = []
    for i in range(0, len(nodes_x), step):
        segs.append([[nodes_x[i], -radii[i]], [nodes_x[i], radii[i]]])
    return segs

def draw_element_lines(ax, nodes_x, get_radius, step=1, color='b', lw=0.5):
    """Draw vertical element lines at node positions."""
    radii = get_radius(nodes_x)
    for i in range(0, len(nodes_x), step):
        ax.plot([nodes_x[i], nodes_x[i]], [-radii[i], radii[i]], color, lw=lw)

def save_snapshot_generic(filename, solver, profile_fn, history_fn=None, title_prefix=""):
    """Generic snapshot: geometry + profiles + optional history.
    profile_fn(ax_data) plots velocity/pressure/etc on given axis.
    history_fn(ax_hist) plots conservation history on given axis."""
    fig, axes = plt.subplots(3, 1, figsize=(10, 10))
    # Geometry
    gx = np.linspace(0, solver.length, 300)
    gr = np.array([solver.get_radius(x) for x in gx])
    axes[0].plot(gx, gr, 'k', lw=2)
    axes[0].plot(gx, -gr, 'k', lw=2)
    radii = solver.get_radius(solver.nodes_x)
    step = max(1, len(solver.nodes_x) // 100)
    for i in range(0, len(solver.nodes_x), step):
        axes[0].plot([solver.nodes_x[i], solver.nodes_x[i]], [-radii[i], radii[i]], 'b', lw=0.5)
    axes[0].set_xlim(-0.2, solver.length + 0.2)
    axes[0].set_ylim(-0.15, 0.15)
    axes[0].set_aspect('equal')
    t_str = f"t={solver.sim_t:.4e}s" if hasattr(solver, 'sim_t') else f"t={solver.sim_time*1000:.1f}ms"
    axes[0].set_title(f"{title_prefix} | {t_str}")
    # Profiles
    profile_fn(axes[1])
    # History
    if history_fn is not None:
        history_fn(axes[2])
    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    print(f"Saved {filename}")
    plt.close(fig)
