

# From:
# https://aistudio.google.com/app/prompts?state=%7B%22ids%22:%5B%221YkBMWN1skW7OUihPhY9y7T5wZjjDDZIF%22%5D,%22action%22:%22open%22,%22userId%22:%22100958146796876347936%22,%22resourceKeys%22:%7B%7D%7D&usp=sharing


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

def triangulate_strip(pts_a, pts_b, mode='shortest'):
    """
    Triangulates the strip between two rows of points.
    mode='parametric': Syncs based on index percentage (0..1).
    mode='shortest': Greedy choice of the shortest diagonal edge.
    """
    count_a = len(pts_a)
    count_b = len(pts_b)
    faces = []
    idx_a = 0
    idx_b = 0
    
    while idx_a < count_a - 1 or idx_b < count_b - 1:
        # Determine which vertex to advance
        choice = None
        
        # 1. Check strict boundaries
        if idx_a == count_a - 1:
            choice = 'B' # A is finished, must advance B
        elif idx_b == count_b - 1:
            choice = 'A' # B is finished, must advance A
        else:
            # 2. Both are available, decide based on Mode
            if mode == 'parametric':
                u_next_a = (idx_a + 1) / (count_a - 1) 
                u_next_b = (idx_b + 1) / (count_b - 1) 
                # Advance the one that is 'lagging behind'
                choice = 'A' if u_next_a <= u_next_b else 'B'
                
            elif mode == 'shortest':
                # Candidate A: creates diagonal (A_next, B_curr)
                # Candidate B: creates diagonal (A_curr, B_next)
                va_curr = pts_a[idx_a]
                vb_curr = pts_b[idx_b]
                va_next = pts_a[idx_a+1]
                vb_next = pts_b[idx_b+1]
                
                # We want the diagonal that splits the quad to be the SHORTER one
                dist_diag_if_A = np.linalg.norm(va_next - vb_curr)
                dist_diag_if_B = np.linalg.norm(va_curr - vb_next)
                
                if dist_diag_if_A < dist_diag_if_B:
                    choice = 'A' # Pick A to create edge (A_next, B_curr)
                elif dist_diag_if_B < dist_diag_if_A:
                    choice = 'B' # Pick B to create edge (A_curr, B_next)
                else:
                    # Tie-break: use parametric
                    u_next_a = (idx_a + 1) / (count_a - 1)
                    u_next_b = (idx_b + 1) / (count_b - 1)
                    choice = 'A' if u_next_a <= u_next_b else 'B'

        # 3. Create Triangle
        if choice == 'A':
            # Triangle (A_curr, B_curr, A_next)
            faces.append([idx_a, count_a + idx_b, idx_a + 1]) 
            idx_a += 1
        else:
            # Triangle (A_curr, B_curr, B_next)
            faces.append([idx_a, count_a + idx_b, count_a + idx_b + 1])
            idx_b += 1
            
    return np.array(faces)

def generate_mesh_variations(nx_top, nx_bottom, ny_rows, mode='shortest'):
    # Interpolate vertex counts
    row_counts = np.linspace(nx_top, nx_bottom, ny_rows).round().astype(int)
    
    all_vertices = []
    row_start_indices = [0]
    
    # Generate Geometry
    for i in range(ny_rows):
        v_progress = i / (ny_rows - 1)
        current_y = 1.0 - v_progress
        # Trapezoid shape
        top_width = 1.0; bottom_width = 2.0 
        current_width = top_width + (bottom_width - top_width) * v_progress
        
        count = row_counts[i]
        xs = np.linspace(-current_width/2, current_width/2, count)
        ys = np.full(count, current_y)
        all_vertices.append(np.column_stack((xs, ys)))
        
        if i < ny_rows - 1:
            row_start_indices.append(row_start_indices[-1] + count)
            
    all_vertices_flat = np.vstack(all_vertices)
    
    # Triangulate Rows
    all_faces = []
    for r in range(ny_rows - 1):
        pts_a = all_vertices[r]
        pts_b = all_vertices[r+1]
        
        # Get relative faces
        local_faces = triangulate_strip(pts_a, pts_b, mode=mode)
        
        # Adjust to global indices
        offset_a = row_start_indices[r]
        offset_b = row_start_indices[r+1]
        threshold = len(pts_a)
        
        for face in local_faces:
            global_face = []
            for idx in face:
                if idx < threshold:
                    global_face.append(offset_a + idx)
                else:
                    global_face.append(offset_b + (idx - threshold))
            all_faces.append(global_face)
            
    return all_vertices_flat, np.array(all_faces), row_counts

# Plotting Comparison
nx1, nx2, ny = 8, 20, 8

fig, axs = plt.subplots(1, 3, figsize=(18, 6))

modes = [('shortest', "Shortest Diagonal"), ('parametric', "Parametric")]

for ax, (m, title) in zip(axs, modes):
    v, f, row_counts = generate_mesh_variations(nx1, nx2, ny, mode=m)
    
    col = PolyCollection(v[f], facecolors='cyan', edgecolors='black', alpha=0.5, linewidths=0.5)
    ax.add_collection(col)
    ax.scatter(v[:,0], v[:,1], s=5, c='red')
    ax.set_aspect('equal')
    ax.autoscale()
    ax.set_title(title)

plt.show()