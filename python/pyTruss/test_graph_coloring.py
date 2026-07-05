# === AUTO-DOC BEGIN ===
"""
@brief Test script for graph coloring on 2D grid trusses.

Builds a grid truss, constructs a point-neighbor adjacency list, runs `color_graph` from
`sparse.py`, and visualizes the coloring using `plot_utils.plot_truss`. The coloring
partitions vertices into independent sets — used by parallel Gauss-Seidel solvers to
update all vertices of one color simultaneously without data races. Includes a local
`build_point_neighbor_list` helper (point-indexed, distinct from `sparse.build_neighbor_list`
which is bond-indexed).
"""
# === AUTO-DOC END ===

import numpy as np
import matplotlib.pyplot as plt
import argparse
from truss import Truss
from sparse import color_graph
import plot_utils as pu

def build_point_neighbor_list(bonds, n_points):
    """Build list of neighboring points for each point"""
    neighbs = [[] for _ in range(n_points)]
    for i_, j_ in bonds:
        neighbs[i_].append(j_)
        neighbs[j_].append(i_)
    return neighbs

def test_truss_coloring(nx=5, ny=5):
    """Test the graph coloring algorithm on a 2D grid truss
    Args:
        nx (int): Number of cells in x direction
        ny (int): Number of cells in y direction
    """
    # Create a truss system
    truss = Truss()
    truss.build_grid_2d(nx=nx, ny=ny, m=1.0, m_end=1000.0, l=1.0, k=10000.0, k_diag=1000.0)
    
    # Get neighbor list for coloring
    n_points = len(truss.points)
    neighs = build_point_neighbor_list(truss.bonds, n_points)
    
    # Color the graph
    colors, color_groups = color_graph(neighs)
    
    # Create a figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))
    
    # Plot original truss
    pu.plot_truss(truss.points, truss.bonds, ax=ax1, point_color='blue')
    ax1.set_title('Original Truss')

    colors_ = [ 'rgbcmykw'[i] for i in colors ]
    
    # Plot the colored truss using plot_truss
    pu.plot_truss(truss.points, truss.bonds, ax=ax2, edge_color='gray', edge_alpha=0.5,   point_color=colors_, point_size=100)
    
    # Add legend manually
    from matplotlib.lines import Line2D
    legend_elements = [Line2D([0], [0], marker='o', color='w',    markerfacecolor='rgbcmykw'[i],   markersize=10, label=f'Color {i}') for i in range(len(color_groups))]
    ax2.legend(handles=legend_elements)    
    ax2.set_title(f'Colored Truss ({len(color_groups)} colors)')

    plt.tight_layout()
    return fig

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Test graph coloring on a 2D grid truss.")
    parser.add_argument("--nx", type=int, default=5)
    parser.add_argument("--ny", type=int, default=5)
    parser.add_argument("--savefig", type=str, default="", help="Path to save the plot.")
    parser.add_argument("--noshow", action="store_true", help="Skip plt.show() for headless execution.")
    args = parser.parse_args()
    fig = test_truss_coloring(nx=args.nx, ny=args.ny)
    if args.savefig:
        fig.savefig(args.savefig, dpi=150)
    if not args.noshow: plt.show()
