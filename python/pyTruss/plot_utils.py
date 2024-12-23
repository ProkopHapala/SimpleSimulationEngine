import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

def plot_truss( points, bonds, ax=None, edge_color='k', edge_alpha=1.0, point_color='b', point_size=20, color_by_stiffness=False, cmap='viridis',  show_colorbar=True, ks=None, label=None):
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
        node_colors: optional list of colors for nodes, overrides point_color if provided
    """
    #if ax is None: _, ax = plt.subplots(figsize=(10, 10))
    if ax is None: ax = plt.gca()
    
    # Prepare line segments for LineCollection
    segments = []
    for i, j in bonds:
        segments.append([points[i, :2], points[j, :2]])
    
    # Create LineCollection
    lc = LineCollection(segments, linewidths=1)
    
    if color_by_stiffness:
        if ks is None:
            raise ValueError("ks (stiffness values) must be provided when color_by_stiffness=True")
        # Normalize colors by stiffness
        lc.set_array(np.array(ks))
        lc.set_cmap(cmap)
        if show_colorbar:
            plt.colorbar(lc, ax=ax, label='Stiffness')
    else:
        lc.set_color(edge_color)
    
    lc.set_alpha(edge_alpha)
    ax.add_collection(lc)
    
    if point_size is not None:
        ax.scatter(points[:, 0], points[:, 1], c=point_color, s=point_size, zorder=2, label=label)
    
    # Set equal aspect ratio and grid
    ax.set_aspect('equal')
    ax.grid(True)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    
    # Update limits
    ax.autoscale()
    
    return ax