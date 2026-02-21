import numpy as np
import pyopencl as cl
import pyopencl.cltypes
import matplotlib.pyplot as plt
import sys
import os

# Ensure the paths are correct for imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'Radiosity')))
try:
    from PolygonRasterization import HexRasterizerNP
    from Radiosity3D import assemble_elements, compute_view_factors_3d, solve_radiosity_system_3d, calculate_temperatures_3d, rasterize_face_to_elems
except ImportError:
    print("Warning: Could not import Radiosity modules.")

from SpacecraftGeometry import generate_christmas_tree

def compute_occlusion(elements, V, F):
    ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)
    
    n_points = len(elements['center'])
    centers = elements['center'].astype(np.float32)
    
    # Assume face IDs are contiguous and correspond to triangles for simplicity in this test
    # We will triangulate faces to ensure they are triangles
    tris = []
    tids = []
    for fi, f in enumerate(F):
        if len(f) < 3: continue
        v0 = V[f[0]]
        for k in range(1, len(f)-1):
            v1 = V[f[k]]
            v2 = V[f[k+1]]
            tris.append([v0, v1, v2])
            tids.append(fi)
    tris_xyz = np.array(tris, dtype=np.float32)
    tri_face_ids = np.array(tids, dtype=np.float32)
    n_tris = len(tri_face_ids)
    
    # Pack points
    pts_packed = np.zeros((n_points, 4), dtype=np.float32)
    pts_packed[:, :3] = centers
    pts_packed[:, 3] = np.arange(n_points, dtype=np.float32) # Simplification: assigning unique ID
    
    # Pack tris
    tri_seq = np.zeros((n_tris*3, 4), dtype=np.float32)
    tri_seq[0::3, :3] = tris_xyz[:, 0, :]
    tri_seq[1::3, :3] = tris_xyz[:, 1, :]
    tri_seq[2::3, :3] = tris_xyz[:, 2, :]
    tri_seq[:, 3] = tri_face_ids.repeat(3)
    
    occ_matrix = np.zeros(n_points * n_points, dtype=np.float32)
    
    mf = cl.mem_flags
    buf_pts = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=pts_packed)
    buf_tris = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=tri_seq)
    buf_occ = cl.Buffer(ctx, mf.WRITE_ONLY, size=occ_matrix.nbytes)
    
    with open('cl/SpacecraftSimulation.cl', 'r') as f:
        kernel_source = f.read()
    
    prg = cl.Program(ctx, kernel_source).build()
    
    prg.compute_occlusion(queue, (n_points, n_points), None,
        buf_pts, buf_tris, buf_occ,
        np.int32(n_points), np.int32(n_tris))
    
    queue.finish()
    cl.enqueue_copy(queue, occ_matrix, buf_occ)
    
    return occ_matrix.reshape((n_points, n_points))

def run_radiosity_pipeline(V, F, hex_size=0.5, P_hot=100.0, hot_face_idx=0):
    print("--- Running Radiosity Pipeline ---")
    
    # Heating per face
    nF = len(F)
    face_heat = np.zeros(nF, dtype=np.float64)
    if 0 <= hot_face_idx < nF:
        face_heat[hot_face_idx] = P_hot

    all_elems = []
    for fi, face in enumerate(F):
        pts = V[np.array(face, dtype=int)]
        elems = rasterize_face_to_elems(pts, face_id=fi, hex_size=hex_size, area_keep_ratio=0.5, verbosity=0)
        all_elems.extend(elems)
        
    elements = assemble_elements(all_elems, eps_front=0.8, eps_back=0.8, face_heat_W_per_m2=face_heat)
    
    # Solve
    Fgeom = compute_view_factors_3d(elements)
    
    print("Computing occlusion on GPU...")
    occ = compute_occlusion(elements, V, F)
    Fgeom *= (1.0 - occ)
    print(f"Occlusion complete. Blocked pairs: {occ.sum()} / {occ.size}")
    
    B = solve_radiosity_system_3d(elements)
    T4 = calculate_temperatures_3d(B, elements)
    
    print(f"Radiosity complete. Max T4: {T4.max():.2f}")
    
    # Plotting
    try:
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        centers = elements['center']
        sc = ax.scatter(centers[:, 0], centers[:, 1], centers[:, 2], c=T4, cmap='inferno', s=16)
        plt.colorbar(sc, label='T^4')
        plt.title('Radiosity Heat Map')
        plt.savefig('radiosity_heatmap.png')
        print("Saved radiosity_heatmap.png")
    except Exception as e:
        print("Plotting failed:", e)

if __name__ == '__main__':
    V, F = generate_christmas_tree()
    run_radiosity_pipeline(V, F)
