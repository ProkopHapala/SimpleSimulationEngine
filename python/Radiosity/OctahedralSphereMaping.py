import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

# --- User Configuration ---
point_3d = np.array([0.6, -0.75, 0.5])

n_segments_per_edge = 3
# --------------------------

# ===== STAGE 1: PROJECTION =====
def map_3d_to_uv(p):
    """Projects a 3D point onto the octahedral UV map."""
    p = p / np.linalg.norm(p)
    uv = p[:2] / (np.abs(p[0]) + np.abs(p[1]) + np.abs(p[2]))
    if p[2] < 0:
        u, v = uv
        uv = np.array([(1 - np.abs(v)) * np.sign(u), (1 - np.abs(u)) * np.sign(v)])
    return uv

# ===== STAGE 2 & 3: LOCALIZATION AND INTERPOLATION =====
def get_barycentric_coords(point, v0, v1, v2):
    """Calculates barycentric coordinates of a point in a 2D triangle."""
    # Using a robust matrix inversion method
    mat = np.array([[v1[0] - v0[0], v2[0] - v0[0]],
                    [v1[1] - v0[1], v2[1] - v0[1]]])
    try:
        w1w2 = np.linalg.solve(mat, point - v0)
        w1, w2 = w1w2
        w0 = 1.0 - w1 - w2
        return np.array([w0, w1, w2])
    except np.linalg.LinAlgError:
        return np.array([np.nan, np.nan, np.nan])

def find_interpolation_data(uv_point, n_points):
    """
    Finds the 3 grid nodes and barycentric coefficients for a given UV point.
    This version uses the corrected logic based on the visualized diagonal.
    """
    n_quads = n_points - 1
    # Convert UV coordinate from [-1, 1] to grid index space [0, n_quads]
    u_grid = (uv_point[0] + 1.0) / 2.0 * n_quads
    v_grid = (uv_point[1] + 1.0) / 2.0 * n_quads

    col_idx = min(int(np.floor(u_grid)), n_quads - 1)
    row_idx = min(int(np.floor(v_grid)), n_quads - 1)

    # Local coordinates within the quad (from 0 to 1)
    u_local = u_grid - col_idx
    v_local = v_grid - row_idx

    # Define the four vertices of the quad in (row, col) format
    v_bl = (row_idx, col_idx)          # Bottom-left
    v_br = (row_idx, col_idx + 1)      # Bottom-right
    v_tl = (row_idx + 1, col_idx)      # Top-left
    v_tr = (row_idx + 1, col_idx + 1)  # Top-right

    mid_point_idx = n_quads / 2.0
    v_indices = []

    # Determine which diagonal is used for this quad (this logic matches the visual plot)
    in_top_right_or_bottom_left_quadrant = \
        (col_idx >= mid_point_idx and row_idx >= mid_point_idx) or \
        (col_idx < mid_point_idx and row_idx < mid_point_idx)

    # --- THE CORRECTED LOGIC BLOCK ---
    if in_top_right_or_bottom_left_quadrant:
        # Diagonal is '\' (from top-left to bottom-right).
        # The two vertices on the diagonal are v_tl and v_br.
        diag_v1, diag_v2 = v_tl, v_br
        # The line equation is u_local + v_local = 1.
        # If u_local + v_local > 1, the point is "above" the line, closer to v_tr.
        if (u_local + v_local) > 1.0:
            off_diag_v = v_tr
        else:
            off_diag_v = v_bl
        v_indices = [diag_v1, diag_v2, off_diag_v]
    else:
        # Diagonal is '/' (from bottom-left to top-right).
        # The two vertices on the diagonal are v_bl and v_tr.
        diag_v1, diag_v2 = v_bl, v_tr
        # The line equation is v_local = u_local.
        # If v_local > u_local, the point is "above" the line, closer to v_tl.
        if v_local > u_local:
            off_diag_v = v_tl
        else:
            off_diag_v = v_br
        v_indices = [diag_v1, diag_v2, off_diag_v]

    # The rest of the function proceeds as before
    u_coords = np.linspace(-1.0, 1.0, n_points)
    v_coords = np.linspace(-1.0, 1.0, n_points)
    v0_uv = np.array([u_coords[v_indices[0][1]], v_coords[v_indices[0][0]]])
    v1_uv = np.array([u_coords[v_indices[1][1]], v_coords[v_indices[1][0]]])
    v2_uv = np.array([u_coords[v_indices[2][1]], v_coords[v_indices[2][0]]])
    
    triangle_verts = np.array([v0_uv, v1_uv, v2_uv])
    coeffs = get_barycentric_coords(uv_point, v0_uv, v1_uv, v2_uv)
    
    return v_indices, coeffs, triangle_verts

# =======================================================
# --- Main Execution & Output ---
# =======================================================
uv_point = map_3d_to_uv(point_3d)
n_points = (2 * n_segments_per_edge) + 1
indices, coeffs, triangle_verts = find_interpolation_data(uv_point, n_points)

print("--- Octahedral Map Interpolation (Corrected Logic) ---")
print(f"Input 3D Point: {tuple(np.round(point_3d, 3))}")
print(f"Mapped UV Coordinate: ({uv_point[0]:.4f}, {uv_point[1]:.4f})\n")
print("Interpolation Triangle Nodes (row, col):")
print(f"  - Node 1 Index: {indices[0]}, Coefficient: {coeffs[0]:.4f}")
print(f"  - Node 2 Index: {indices[1]}, Coefficient: {coeffs[1]:.4f}")
print(f"  - Node 3 Index: {indices[2]}, Coefficient: {coeffs[2]:.4f}")
print(f"Sum of Coefficients: {np.sum(coeffs):.4f}")

# =======================================================
# --- Visualization ---
# =======================================================
fig, ax = plt.subplots(figsize=(14, 14))
n_quads = n_points - 1
u_coords = np.linspace(-1.0, 1.0, n_points)
v_coords = np.linspace(-1.0, 1.0, n_points)

# Plot grid and diagonals (same as before)
ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.5, alpha=0.5)
ax.set_xticks(u_coords, minor=True); ax.set_yticks(v_coords, minor=True)
mid_point_idx = n_quads / 2.0
for i in range(n_quads):
    for j in range(n_quads):
        u_start, u_end = u_coords[j], u_coords[j+1]
        v_start, v_end = v_coords[i], v_coords[i+1]
        in_tr_or_bl = (j >= mid_point_idx and i >= mid_point_idx) or (j < mid_point_idx and i < mid_point_idx)
        if in_tr_or_bl: ax.plot([u_start, u_end], [v_end, v_start], color='cyan', linewidth=2)
        else: ax.plot([u_start, u_end], [v_start, v_end], color='cyan', linewidth=2)

# Plot octahedron face boundaries
ax.plot([0, 1], [1, 0], 'k-', [1, 0], [0, -1], 'k-', [0, -1], [-1, 0], 'k-', [-1, 0], [0, 1], 'k-', lw=3)
ax.plot([0, 0], [-1, 1], 'k-', [-1, 1], [0, 0], 'k-', lw=3)

# Highlight interpolation triangle and point
triangle_patch = patches.Polygon(triangle_verts, closed=True, facecolor='yellow', alpha=0.7, edgecolor='orange', linewidth=2)
ax.add_patch(triangle_patch)
ax.plot(uv_point[0], uv_point[1], '*', color='green', markersize=20, markeredgecolor='black', label='Target Point')

# Plot Main Octahedron Vertices and Labels
octa_vertices = {
    0: ((0, 0, 1), "North"), 1: ((1, 0, 0), ""), 2: ((0, 1, 0), ""),
    3: ((-1, 0, 0), ""), 4: ((0, -1, 0), ""), 5: ((0, 0, -1), "South"),
}
uv_map = {
    0: (0, 0), 1: (1, 0), 2: (0, 1), 3: (-1, 0), 4: (0, -1),
    5: [(-1, -1), (1, -1), (-1, 1), (1, 1)]
}
label_style = dict(boxstyle='round,pad=0.3', fc='yellow', ec='black', lw=1, alpha=0.9)
for idx in range(6):
    coords_3d, name = octa_vertices[idx]
    label_text = f"{idx}:{coords_3d} - {name} Pole" if name else f"{idx}:{coords_3d}"
    if idx != 5:
        u, v = uv_map[idx]
        ax.plot(u, v, 'o', markersize=12, color='orange' if idx==0 else 'red')
        ax.text(u, v, label_text, fontsize=12, ha='center', va='center', bbox=label_style, weight='bold')
    else:
        for u, v in uv_map[idx]:
            ax.plot(u, v, 'ro', markersize=12)
            ax.text(u, v, label_text, fontsize=12, ha='center', va='center', bbox=label_style, weight='bold')

# Final plot styling
ax.set_aspect('equal', adjustable='box')
ax.set_xlim(-1.1, 1.1); ax.set_ylim(-1.1, 1.1)
ax.set_xlabel("U Coordinate", fontsize=14); ax.set_ylabel("V Coordinate", fontsize=14)
ax.set_title("Octahedral Map Interpolation (Corrected Logic)", fontsize=16)
ax.legend()
plt.show()