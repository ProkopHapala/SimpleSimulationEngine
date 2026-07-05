import numpy as np
import pyopencl as cl
import pyopencl.cltypes
import matplotlib.pyplot as plt
import sys
import os
import math

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

def _triangulate(V, F):
    tris = []; tids = []
    for fi, f in enumerate(F):
        if len(f) < 3: continue
        v0 = V[f[0]]
        for k in range(1, len(f)-1):
            v1 = V[f[k]]; v2 = V[f[k+1]]
            tris.append([v0, v1, v2]); tids.append(fi)
    return np.array(tris, dtype=np.float32), np.array(tids, dtype=np.float32)

def _build_clusters(tris_xyz, tri_face_ids, cluster_size=64):
    n_tris = len(tri_face_ids)
    centroids = tris_xyz.mean(axis=1)  # (n_tris, 3)
    # Spatial sort: z primary, x secondary
    sort_key = centroids[:, 2] * 1e6 + centroids[:, 0]
    sort_idx = np.argsort(sort_key)
    sorted_tris = tris_xyz[sort_idx]
    sorted_tids = tri_face_ids[sort_idx]
    n_clusters = (n_tris + cluster_size - 1) // cluster_size
    cluster_spheres = np.zeros((n_clusters, 4), dtype=np.float32)
    cluster_ranges = np.zeros((n_clusters, 2), dtype=np.int32)
    for ci in range(n_clusters):
        start = ci * cluster_size
        end = min(start + cluster_size, n_tris)
        chunk_pts = sorted_tris[start:end].reshape(-1, 3)
        cmin = chunk_pts.min(axis=0); cmax = chunk_pts.max(axis=0)
        center = (cmin + cmax) / 2.0
        radius = np.max(np.linalg.norm(chunk_pts - center, axis=1)) + 1e-4
        cluster_spheres[ci] = [center[0], center[1], center[2], radius]
        cluster_ranges[ci] = [start, end - start]
    # Pack triangles as sequential float4s (A,B,C) with w=face_id
    tri_seq = np.zeros((n_tris * 3, 4), dtype=np.float32)
    tri_seq[0::3, :3] = sorted_tris[:, 0, :]
    tri_seq[1::3, :3] = sorted_tris[:, 1, :]
    tri_seq[2::3, :3] = sorted_tris[:, 2, :]
    tri_seq[0::3, 3] = sorted_tids
    tri_seq[1::3, 3] = sorted_tids
    tri_seq[2::3, 3] = sorted_tids
    return tri_seq, cluster_spheres, cluster_ranges, n_clusters

def compute_occlusion_tiled(elements, V, F, cluster_size=32, verbosity=0):
    assert cluster_size <= 32, "SpacecraftSimulation.cl CLUSTER_CAPACITY is 32"
    n_points = len(elements['center'])
    centers = elements['center'].astype(np.float32)
    tris_xyz, tri_face_ids = _triangulate(V, F)
    n_tris = len(tri_face_ids)
    if n_tris == 0:
        return np.zeros((n_points, n_points), dtype=np.float32)
    tri_seq, cluster_spheres, cluster_ranges, n_clusters = _build_clusters(tris_xyz, tri_face_ids, cluster_size)
    if verbosity > 0:
        print(f"  Tiled occlusion: {n_points} points, {n_tris} tris, {n_clusters} clusters (size={cluster_size})")

    pts_packed = np.zeros((n_points, 4), dtype=np.float32)
    pts_packed[:, :3] = centers
    pts_packed[:, 3] = np.arange(n_points, dtype=np.float32)

    occ_matrix = np.zeros(n_points * n_points, dtype=np.float32)

    ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)
    mf = cl.mem_flags
    buf_pts = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=pts_packed)
    buf_tris = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=tri_seq)
    buf_spheres = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=cluster_spheres)
    buf_ranges = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=cluster_ranges)
    buf_occ = cl.Buffer(ctx, mf.WRITE_ONLY, size=occ_matrix.nbytes)

    cl_path = os.path.join(os.path.dirname(__file__), 'cl', 'SpacecraftSimulation.cl')
    with open(cl_path, 'r') as f:
        kernel_source = f.read()
    prg = cl.Program(ctx, kernel_source).build()

    TILE = 8
    gw = (math.ceil(n_points / TILE) * TILE, math.ceil(n_points / TILE) * TILE)
    lw = (TILE, TILE)

    prg.compute_occlusion_tiled(queue, gw, lw,
        buf_pts, buf_tris, buf_spheres, buf_ranges, buf_occ,
        np.int32(n_points), np.int32(n_clusters))

    queue.finish()
    cl.enqueue_copy(queue, occ_matrix, buf_occ)
    return occ_matrix.reshape((n_points, n_points))

def _pack_element_arrays(elements):
    centers = elements['center'].astype(np.float32)
    normals = elements['normal'].astype(np.float32)
    areas = elements['area'].astype(np.float32)
    face_ids = elements.get('face_id', np.arange(len(centers), dtype=np.int32)).astype(np.float32)
    pts_packed = np.zeros((len(centers), 4), dtype=np.float32)
    norm_area = np.zeros((len(centers), 4), dtype=np.float32)
    pts_packed[:, :3] = centers
    pts_packed[:, 3] = face_ids
    norm_area[:, :3] = normals
    norm_area[:, 3] = areas
    return pts_packed, norm_area

def compute_sparse_channels(elements, V, F, kmax=32, min_weight=1e-6, max_dist=0.0, cluster_size=32, row_chunk=256, verbosity=0):
    assert kmax <= 32, "compute_sparse_channels kernel supports kmax <= 32"
    assert cluster_size <= 32, "SpacecraftSimulation.cl CLUSTER_CAPACITY is 32"
    n_points = len(elements['center'])
    tris_xyz, tri_face_ids = _triangulate(V, F)
    n_tris = len(tri_face_ids)
    idx = -np.ones((n_points, kmax), dtype=np.int32)
    weights = np.zeros((n_points, kmax), dtype=np.float32)
    if n_points == 0 or n_tris == 0:
        return idx, weights
    tri_seq, cluster_spheres, cluster_ranges, n_clusters = _build_clusters(tris_xyz, tri_face_ids, cluster_size)
    pts_packed, norm_area = _pack_element_arrays(elements)
    if verbosity > 0:
        print(f"  Sparse channels: {n_points} points, {n_tris} tris, {n_clusters} clusters, kmax={kmax}, row_chunk={row_chunk}")

    ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)
    mf = cl.mem_flags
    buf_pts = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=pts_packed)
    buf_norm_area = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=norm_area)
    buf_tris = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=tri_seq)
    buf_spheres = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=cluster_spheres)
    buf_ranges = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=cluster_ranges)
    chunk_idx = -np.ones((row_chunk, kmax), dtype=np.int32)
    chunk_w = np.zeros((row_chunk, kmax), dtype=np.float32)
    buf_idx = cl.Buffer(ctx, mf.WRITE_ONLY, size=chunk_idx.nbytes)
    buf_w = cl.Buffer(ctx, mf.WRITE_ONLY, size=chunk_w.nbytes)

    cl_path = os.path.join(os.path.dirname(__file__), 'cl', 'SpacecraftSimulation.cl')
    with open(cl_path, 'r') as f:
        kernel_source = f.read()
    prg = cl.Program(ctx, kernel_source).build()

    lsz = 64
    k_sparse = prg.compute_sparse_channels
    for row0 in range(0, n_points, row_chunk):
        n_rows = min(row_chunk, n_points - row0)
        gw = (math.ceil(n_rows / lsz) * lsz,)
        k_sparse(queue, gw, (lsz,),
            buf_pts, buf_norm_area, buf_tris, buf_spheres, buf_ranges, buf_idx, buf_w,
            np.int32(n_points), np.int32(n_clusters), np.int32(row0), np.int32(n_rows), np.int32(kmax), np.float32(min_weight), np.float32(max_dist))
        queue.finish()
        cl.enqueue_copy(queue, chunk_idx, buf_idx)
        cl.enqueue_copy(queue, chunk_w, buf_w)
        idx[row0:row0+n_rows, :] = chunk_idx[:n_rows, :]
        weights[row0:row0+n_rows, :] = chunk_w[:n_rows, :]
        if verbosity > 1:
            print(f"    rows {row0}:{row0+n_rows}, nnz={np.count_nonzero(chunk_idx[:n_rows] >= 0)}")
    return idx, weights

def solve_radiosity_sparse(elements, idx, weights, n_iter=80, damping=1.0):
    area = elements['area'] + 1e-12
    P_density = elements['heat_in'] / area
    rho_front = elements['rho_front']
    rho_back = elements['rho_back']
    B = P_density.copy()
    for _ in range(n_iter):
        Bnew = P_density.copy()
        for i in range(len(B)):
            for s in range(idx.shape[1]):
                j = idx[i, s]
                if j < 0: continue
                w = weights[i, s]
                rho = rho_front[i] if w >= 0.0 else rho_back[i]
                Bnew[i] += rho * abs(w) * B[j]
        B = (1.0 - damping) * B + damping * Bnew
    return B

def calculate_temperatures_sparse(B, elements, idx, weights):
    area = elements['area'] + 1e-12
    P_density = elements['heat_in'] / area
    H_front = np.zeros_like(B)
    H_back = np.zeros_like(B)
    for i in range(len(B)):
        for s in range(idx.shape[1]):
            j = idx[i, s]
            if j < 0: continue
            w = weights[i, s]
            if w >= 0.0:
                H_front[i] += abs(w) * B[j]
            else:
                H_back[i] += abs(w) * B[j]
    num = P_density + elements['eps_front'] * H_front + elements['eps_back'] * H_back
    den = elements['eps_front'] + elements['eps_back'] + 1e-12
    return num / den

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
    idx, weights = compute_sparse_channels(elements, V, F, kmax=32, min_weight=1e-6, max_dist=0.0, cluster_size=32, verbosity=1)
    print(f"Sparse channels complete. nnz={np.count_nonzero(idx >= 0)} / {idx.size}")
    
    B = solve_radiosity_sparse(elements, idx, weights, n_iter=80, damping=1.0)
    T4 = calculate_temperatures_sparse(B, elements, idx, weights)
    
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
