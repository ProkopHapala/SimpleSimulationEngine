import numpy as np

def generate_christmas_tree(
    n_levels=4,
    level_spacing=2.0,
    n_branches_per_level=6,
    branch_length=5.0,
    branch_width=1.0,
    tilt_angle_deg=-15.0,
    radial_segments=3,
    chord_segments=2
):
    """
    Generate a 'Christmas-tree' spacecraft geometry.
    Returns:
        vertices: (N, 3) float array of vertex positions.
        faces: list of lists, where each inner list contains vertex indices for a polygon.
    """
    vertices = []
    faces = []
    
    # Simple central backbone (a tube or just a line, let's just make the branches for now)
    # Actually, let's represent the backbone as a thin vertical quad or just ignore it if it's too thin.
    
    tilt_rad = np.radians(tilt_angle_deg)
    
    for level in range(n_levels):
        z_base = level * level_spacing
        
        # Scale length to make it look like a tree
        current_length = branch_length * (1.0 - 0.15 * level)
        current_width = branch_width * (1.0 - 0.1 * level)
        
        for branch in range(n_branches_per_level):
            angle = branch * (2 * np.pi / n_branches_per_level)
            
            # Direction of the branch
            dir_x = np.cos(angle)
            dir_y = np.sin(angle)
            
            # Create a simple fan / rectangular strip for the branch
            # We will create a grid of vertices
            
            v_start_idx = len(vertices)
            
            for r in range(radial_segments + 1):
                r_frac = r / radial_segments
                r_dist = current_length * r_frac
                
                # Apply tilt
                z = z_base + r_dist * np.sin(tilt_rad)
                xy_dist = r_dist * np.cos(tilt_rad)
                
                center_x = xy_dist * dir_x
                center_y = xy_dist * dir_y
                
                # Width vector (perpendicular to direction in xy plane)
                w_x = -dir_y
                w_y = dir_x
                
                for c in range(chord_segments + 1):
                    c_frac = c / chord_segments - 0.5 # -0.5 to 0.5
                    w_dist = current_width * c_frac
                    
                    x = center_x + w_x * w_dist
                    y = center_y + w_y * w_dist
                    
                    vertices.append([x, y, z])
            
            # Create faces
            for r in range(radial_segments):
                for c in range(chord_segments):
                    i0 = v_start_idx + r * (chord_segments + 1) + c
                    i1 = i0 + 1
                    i2 = i0 + (chord_segments + 1)
                    i3 = i2 + 1
                    
                    # two triangles per quad
                    faces.append([i0, i1, i3])
                    faces.append([i0, i3, i2])

    return np.array(vertices, dtype=np.float64), faces

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    
    verts, faces = generate_christmas_tree()
    print(f"Generated {len(verts)} vertices, {len(faces)} faces.")
    
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    polys = [[verts[idx] for idx in f] for f in faces]
    mesh = Poly3DCollection(polys, alpha=0.5, facecolor='cyan', edgecolor='k', linewidths=0.5)
    ax.add_collection3d(mesh)
    
    all_pts = verts
    mins = all_pts.min(axis=0)
    maxs = all_pts.max(axis=0)
    span = max(maxs - mins)
    mid = 0.5 * (mins + maxs)
    
    ax.set_xlim(mid[0] - span/2, mid[0] + span/2)
    ax.set_ylim(mid[1] - span/2, mid[1] + span/2)
    ax.set_zlim(mid[2] - span/2, mid[2] + span/2)
    
    plt.title("Spacecraft Geometry Preview")
    plt.savefig('geom_preview.png')
    print("Saved geom_preview.png")
