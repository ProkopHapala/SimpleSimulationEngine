# === AUTO-DOC BEGIN ===
"""
@brief Matplotlib visualization utilities for truss structures and graph coloring.

**plot_truss** renders points and bonds as a LineCollection with optional per-bond stiffness
coloring and per-node RGBA colors. **plot_graph_coloring** overlays graph coloring partitions
with annotated color IDs. Used by `run_vbd_cloth.py`, `run_vbd_cloth_new.py`, and
`test_graph_coloring.py`. Design: all functions accept an optional `ax` parameter for
composable subplot layouts.
"""
# === AUTO-DOC END ===

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.collections import LineCollection

def plot_truss( points, bonds, ax=None, edge_color='k', edge_alpha=1.0, point_color='b', point_size=20, color_by_stiffness=False, cmap='viridis',  show_colorbar=True, ks=None, label=None, node_colors=None, margin=0.2):
    """
    Plot the truss efficiently using LineCollection.
    
    Args:
        points: Nx3 array of point coordinates
        bonds: List of (i,j) tuples defining connections
        ax: matplotlib axis to plot on (creates new if None)
        edge_color: color for edges (ignored if color_by_stiffness=True)
        edge_alpha: transparency for edges
        point_color: color for points (ignored if node_colors is provided)
        point_size: size of points
        color_by_stiffness: if True, color edges by their stiffness
        cmap: colormap to use when color_by_stiffness=True
        show_colorbar: whether to show colorbar when color_by_stiffness=True
        ks: array of stiffness values for each bond (required if color_by_stiffness=True)
        label: label for legend
        node_colors: optional per-node RGBA colors, overrides point_color if provided
    """
    #if ax is None: _, ax = plt.subplots(figsize=(10, 10))
    if ax is None: ax = plt.gca()
    segments = []
    for i, j in bonds: segments.append([points[i, :2], points[j, :2]])
    lc = LineCollection(segments, linewidths=1)
    if color_by_stiffness:
        if ks is None:
            raise ValueError("ks (stiffness values) must be provided when color_by_stiffness=True")
        lc.set_array(np.array(ks))
        lc.set_cmap(cmap)
        if show_colorbar: plt.colorbar(lc, ax=ax, label='Stiffness')
    else:
        lc.set_color(edge_color)
    lc.set_alpha(edge_alpha)
    ax.add_collection(lc)
    
    if point_size is not None:
        color_data = node_colors if node_colors is not None else point_color
        ax.scatter(points[:, 0], points[:, 1], c=color_data, s=point_size, zorder=2, label=label)
    ax.set_aspect('equal')
    ax.grid(True)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    #ax.set_xlim(points[:, 0].min()-margin, points[:, 0].max()+margin)
    #ax.set_ylim(points[:, 1].min()-margin, points[:, 1].max()+margin)
    #ax.autoscale()
    return ax

def plot_graph_coloring(truss, vertex_colors, partitions, ax=None, cmap_name="tab20"):
    """Plot truss with vertices colored by graph coloring partition."""
    if len(vertex_colors) == 0:
        return ax
    if ax is None:
        ax = plt.gca()
    n_colors = int(vertex_colors.max()) + 1
    palette = cm.get_cmap(cmap_name, n_colors)(np.linspace(0.0, 1.0, n_colors))
    node_colors = palette[vertex_colors.astype(int)]
    ax = plot_truss(truss.points, truss.bonds, ax=ax, edge_color='0.3', edge_alpha=0.5, point_size=60, node_colors=node_colors)
    for color_id, verts in enumerate(partitions):
        p = truss.points[verts[0], :2]
        ax.text(p[0], p[1], str(color_id), color='k', fontsize=10, ha='center', va='center')
    ax.set_title("Graph coloring (color id annotated)")
    return ax