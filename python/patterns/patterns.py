import numpy as np

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
import math

def generate_star_segments(n_fold):
    """
    Generates the line segments for an n-fold star centered at the origin,
    consisting of radial lines from the center.

    Args:
        n_fold: The number of radial lines (and points) in the star.

    Returns:
        A list of tuples, where each tuple represents a line segment
        [(x1, y1), (x2, y2)] of the star, originating from (0,0).
    """
    outer_radius = 1.0
    angles = np.linspace(0, 2 * np.pi, n_fold, endpoint=False)

    segments = []
    center = (0, 0)
    for angle in angles:
        x = outer_radius * np.cos(angle)
        y = outer_radius * np.sin(angle)
        segments.append([center, (x, y)])

    return segments

def generate_triangular_grid(rows, cols, spacing=1.0):
    nodes = []
    row_height = spacing * np.sqrt(3) / 2
    for r in range(rows):
        y = r * row_height
        start_x = (r % 2) * spacing / 2  # Stagger rows
        for c in range(cols):
            x = start_x + c * spacing
            nodes.append((x, y,0.0))
    return nodes

def generate_square_grid(rows, cols, spacing=1.0):
    nodes = []
    for r in range(rows):
        y = r * spacing
        for c in range(cols):
            x = c * spacing
            nodes.append((x, y,0.0))
    return nodes

def generate_hexagonal_grid(rows, cols, spacing=1.0):
    nodes = []
    row_height = spacing * np.sqrt(3) / 2
    for r in range(rows):
        y       = r * row_height
        start_x = (r % 2) * spacing / 2  # Stagger rows
        for c in range(cols):
            x = start_x + c * spacing
            nodes.append((x, y,0.0))
            nodes.append((x + spacing/2, y + row_height*(1./3.),np.pi))
    return nodes

def place_and_rotate_stars(n_fold, grid_type, rotation_angle, grid_nodes, star_segments):
    """
    Places and rotates star segments at each grid node.

    Args:
        n_fold: The number of points in the star (not directly used in transformation, but part of definition).
        grid_type: The type of grid (3, 4, or 6) (not directly used in transformation, but part of definition).
        rotation_angle: The rotation angle in degrees.
        grid_nodes: A list of (x, y, angle) coordinates for the grid nodes.
        star_segments: A list of line segments [(x1, y1), (x2, y2)] for a single star centered at the origin.

    Returns:
        A list of transformed line segments for all stars placed on the grid.
    """
    angle_rad = math.radians(rotation_angle)
    transformed_segments = []

    for node in grid_nodes:
        node_x, node_y, a = node
        a += angle_rad
        for segment in star_segments:
            transformed_segment = []
            for point in segment:
                x, y = point
                # Apply rotation around the origin
                x_rotated = x * math.cos(a) - y * math.sin(a)
                y_rotated = x * math.sin(a) + y * math.cos(a)
                # Translate to the grid node
                x_transformed = x_rotated + node_x
                y_transformed = y_rotated + node_y
                transformed_segment.append((x_transformed, y_transformed))
            transformed_segments.append(transformed_segment)

    return transformed_segments

def plot_chiral_pattern(n_fold, grid_type, rotation_angle, grid_rows, grid_cols, spacing):
    """
    Generates and plots a single chiral pattern for a given combination
    of star type, grid type, and rotation angle.

    Args:
        n_fold: The number of radial lines in the star.
        grid_type: The type of grid (3 for triangular, 4 for square, or 6 for hexagonal).
        rotation_angle: The rotation angle in degrees.
        grid_rows: The number of rows in the grid.
        grid_cols: The number of columns in the grid.
        spacing: The spacing between grid nodes.
    """
    # Generate the star segments
    star_segments = generate_star_segments(n_fold)

    # Determine the grid type and generate the grid nodes
    if grid_type == 3:
        grid_nodes = generate_triangular_grid(grid_rows, grid_cols, spacing)
        grid_name = "Triangular"
    elif grid_type == 4:
        grid_nodes = generate_square_grid(grid_rows, grid_cols, spacing)
        grid_name = "Square"
    elif grid_type == 6:
        grid_nodes = generate_hexagonal_grid(grid_rows, grid_cols, spacing)
        grid_name = "Hexagonal"
    else:
        print(f"Warning: Invalid grid_type {grid_type}. Skipping this pattern.")
        return

    # Place and rotate the stars on the grid
    transformed_segments = place_and_rotate_stars(n_fold, grid_type, rotation_angle, grid_nodes, star_segments)

    # Create a new figure and axes
    fig, ax = plt.subplots(figsize=(8, 8))

    # Create a LineCollection object
    line_collection = LineCollection(transformed_segments)

    # Add the LineCollection to the axes
    ax.add_collection(line_collection)

    # Set the aspect of the plot to be equal
    ax.set_aspect('equal', adjustable='box')

    # Automatically adjust the plot limits
    ax.autoscale_view()

    # Set title and labels
    ax.set_title(f'{n_fold}-fold star on {grid_name} grid, rotated {rotation_angle}°')
    ax.set_xlabel("X-coordinate")
    ax.set_ylabel("Y-coordinate")

    # Display the plot
    plt.show()


# -----------------------------------------------------------------------------
# General dictionary-driven pattern API
# Lattice defined by (a, b, gamma_deg). Shapes placed by fractional coords in unit cell.
# Each shape: { name, pos:[u,v], angle_deg, arms, lengths:[...], angles_deg:[...] }
# lengths has m elements, angles_deg has m-1 turning angles between consecutive segments.

def lattice_vectors(lattice):
    a = lattice.get('a', 1.0); b = lattice.get('b', 1.0); g = math.radians(lattice.get('gamma', 60.0))
    a_vec = (a, 0.0)
    b_vec = (b*math.cos(g), b*math.sin(g))
    return a_vec, b_vec

def lattice_nodes(rows, cols, a_vec, b_vec):
    nodes = []
    for i in range(rows):
        for j in range(cols):
            x = i*a_vec[0] + j*b_vec[0]
            y = i*a_vec[1] + j*b_vec[1]
            nodes.append((x, y))
    return nodes

def polyline_points(lengths, angles_deg):
    # Start at origin, initial heading along +x; angles are incremental turns between segments
    pts = [(0.0, 0.0)]
    ang = 0.0
    for i, L in enumerate(lengths):
        x = pts[-1][0] + L*math.cos(ang)
        y = pts[-1][1] + L*math.sin(ang)
        pts.append((x, y))
        if i < len(angles_deg): ang += math.radians(angles_deg[i])
    return pts

def segments_from_points(pts):
    return [[pts[k], pts[k+1]] for k in range(len(pts)-1)]

def rot(p, a):
    return (p[0]*math.cos(a) - p[1]*math.sin(a), p[0]*math.sin(a) + p[1]*math.cos(a))

def transform_segments(segs, a, tx, ty):
    out = []
    for s in segs:
        p0 = rot(s[0], a); p1 = rot(s[1], a)
        out.append([(p0[0]+tx, p0[1]+ty), (p1[0]+tx, p1[1]+ty)])
    return out

def shape_base_segments(arms, lengths, angles_deg):
    base = segments_from_points(polyline_points(lengths, angles_deg))
    segs = []
    for k in range(arms):
        ak = 2*math.pi*k/arms
        segs += transform_segments(base, ak, 0.0, 0.0)
    return segs

def plot_pattern_from_dict(cfg, rows, cols, lw=2.0):
    """
    cfg: {
      'lattice': { 'a':float, 'b':float, 'gamma':float },
      'shapes': [ { 'name':str, 'pos':[u,v], 'angle':float, 'arms':int, 'lengths':[...], 'angles':[...]} ]
    }
    rows, cols: number of cells to tile in a and b directions
    """
    a_vec, b_vec = lattice_vectors(cfg['lattice'])
    nodes        = lattice_nodes(rows, cols, a_vec, b_vec)

    all_segs = []
    for sh in cfg.get('shapes', []):
        arms = int(sh['arms'])
        lengths = list(sh['lengths'])
        angles  = list(sh.get('angles', sh.get('angles', [])))
        pos     = sh.get('pos', sh.get('position', [0.0, 0.0]))
        ang0    = math.radians(sh.get('angle', sh.get('angle', 0.0)))
        color   = sh.get('color', 'black')

        base = shape_base_segments(arms, lengths, angles)
        shape_segs = []
        for X, Y in nodes:
            tx = X + pos[0]*a_vec[0] + pos[1]*b_vec[0]
            ty = Y + pos[0]*a_vec[1] + pos[1]*b_vec[1]
            shape_segs += transform_segments(base, ang0, tx, ty)
        all_segs.append((shape_segs, color))

    fig, ax = plt.subplots(figsize=(8, 8))
    for segs, color in all_segs:
        ax.add_collection(LineCollection(segs, colors=color, linewidths=lw))
    ax.set_aspect('equal', adjustable='box')
    ax.autoscale_view()
    ax.set_title('General lattice pattern')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    

if __name__ == "__main__":
    # how to run:
    # python -u -m python.patterns.patterns

    # Example usage of the plot_chiral_pattern function with potentially corrected hexagonal grid
    # pattern_parameters_updated = [
    #     (3, 6, -5, 5, 5, 2.0 ),  # 3-fold star on hexagonal grid, rotated 30 deg
    #     (4, 4, 45/2., 10, 10, 1.4),  # 4-fold star on square grid, rotated 0 deg
    #     (6, 3, 15, 10, 10, 1.5),  # 6-fold star on triangular grid, rotated 0 deg
    # ]
    # for params in pattern_parameters_updated:
    #     plot_chiral_pattern(*params)
    # plt.show()



    cfg = {
        'lattice': { 'a': 1.5, 'b': 1.5, 'gamma': 60.0 },  # triangular lattice
        'shapes': [

            # Corrugation 1
            # { 'name':'star3', 'color':'k', 'pos':[0.0,  0.0 ],  'angle':   0.0, 'arms': 3, 'lengths':[0.5,0.5], 'angles':[-60] },
            # { 'name':'star4', 'color':'r', 'pos':[1/3., 1/3.],  'angle':  60.0, 'arms': 3, 'lengths':[0.5,0.5], 'angles':[-120] },

            # Corrugation 2
            # { 'name':'star3', 'color':'k',   'pos':[0.0,  0.0 ],  'angle':   0.0, 'arms': 3, 'lengths':[0.5], 'angles':[] },
            # { 'name':'star4', 'color':'r', 'pos':[1/3., 1/3.],  'angle':  60.0, 'arms': 3, 'lengths':[0.5,0.5], 'angles':[-60] },

            # Corrugation 3
            { 'name':'star3', 'color':'k', 'pos':[0.0,  0.0 ],  'angle':  0.0, 'arms': 3, 'lengths':[0.5        ], 'angles':[]      },
            { 'name':'star4', 'color':'k', 'pos':[1/3., 1/3.],  'angle': 60.0, 'arms': 3, 'lengths':[0.5,0.5    ], 'angles':[-120]  },
            { 'name':'star4', 'color':'r', 'pos':[2/3., 2/3.],  'angle': 30.0, 'arms': 3, 'lengths':[0.5,0.3,0.5], 'angles':[60,30] },
        ]
    }

    # cfg = {
    #     'lattice': { 'a': 1.5, 'b': 1.5, 'gamma': 90.0 },  # triangular lattice
    #     'shapes': [
    #         # Corrugation 1
    #         # { 'name':'star3', 'color':'k', 'pos':[0.0,  0.0 ],  'angle':   0.0, 'arms': 4, 'lengths':[0.5,0.5], 'angles':[-90] },
    #         # { 'name':'star4', 'color':'r', 'pos':[1/2., 1/2.],  'angle':  90.0, 'arms': 4, 'lengths':[0.5,0.5], 'angles':[-90] },
        
    #         # Corrugation 2
    #         # { 'name':'star3', 'color':'k', 'pos':[0.0,  0.0 ],  'angle':   0.0, 'arms': 4, 'lengths':[0.5,0.5], 'angles':[-90] },
    #         # { 'name':'star4', 'color':'r', 'pos':[1/2., 1/2.],  'angle':  90.0, 'arms': 4, 'lengths':[0.5,0.5*np.sqrt(0.5)], 'angles':[-45] },
        
    #         # Corrugation 3
    #         # { 'name':'star3', 'color':'k', 'pos':[0.0,  0.0 ],  'angle':   0.0, 'arms': 4, 'lengths':[0.5,0.5], 'angles':[-90] },
    #         # { 'name':'star4', 'color':'r', 'pos':[1/2., 1/2.],  'angle':  90.0, 'arms': 4, 'lengths':[1.0,0.5], 'angles':[90] },

    #         # Corrugation 4
    #         # { 'name':'star3', 'color':'k', 'pos':[0.0,  0.0 ],  'angle':   0.0, 'arms': 2, 'lengths':[0.5,1.0,0.5,0.25], 'angles':[-90,-90,-90] },
    #         # { 'name':'star4', 'color':'r', 'pos':[1/2., 1/2.],  'angle':  90.0, 'arms': 2, 'lengths':[0.5,1.0,0.5], 'angles':[-90,-90] },
            
            
    #         # { 'name':'star4', 'color':'r', 'pos':[1/2., 1/2.],  'angle':  90.0, 'arms': 2, 'lengths':[0.5,1.0], 'angles':[-90] },
    #         # { 'name':'star4', 'color':'r', 'pos':[0.0, -1/2  ],  'angle':  0.0, 'arms': 2, 'lengths':[0.25,0.5], 'angles':[90] },
    #     ]
    # }



    plot_pattern_from_dict(cfg, rows=5, cols=5, lw=2.0)

    plt.savefig('pattern.svg')
    plt.show()





    # Example usage (for testing)
    # triangular_nodes = generate_triangular_grid(5, 5)
    # square_nodes = generate_square_grid(5, 5)
    # hexagonal_nodes = generate_hexagonal_grid(5, 5)
    # print("Triangular grid nodes:", triangular_nodes[:5])
    # print("Square grid nodes:", square_nodes[:5])
    # print("Hexagonal grid nodes:", hexagonal_nodes[:5])