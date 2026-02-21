import numpy as np
import pyopencl as cl
import pyopencl.array
import matplotlib.pyplot as plt
import sys
import os

from SpacecraftGeometry import generate_christmas_tree

def triangulate_faces(vertices, faces):
    tris = []
    tids = []
    for fi, f in enumerate(faces):
        if len(f) < 3: continue
        v0 = vertices[f[0]]
        for k in range(1, len(f)-1):
            v1 = vertices[f[k]]
            v2 = vertices[f[k+1]]
            tris.append([v0, v1, v2])
            tids.append(fi)
    if len(tris) == 0:
        return np.zeros((0,3,3), np.float32), np.zeros((0,), np.int32)
    return np.asarray(tris, dtype=np.float32), np.asarray(tids, dtype=np.int32)

def run_attenuation_pipeline(V, F):
    print("--- Running Attenuation Pipeline ---")
    
    tris_xyz, tri_face_ids = triangulate_faces(V, F)
    ntris = len(tri_face_ids)
    
    print(f"Total Triangles: {ntris}")
    
    # 1. Setup OpenCL
    ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)
    
    # Pack triangles into float4
    tri_seq = np.zeros((ntris*3, 4), dtype=np.float32)
    tri_seq[0::3, :3] = tris_xyz[:, 0, :]
    tri_seq[1::3, :3] = tris_xyz[:, 1, :]
    tri_seq[2::3, :3] = tris_xyz[:, 2, :]
    tri_seq[:, 3] = tri_face_ids.repeat(3)
    
    # Material properties
    mu_vals = np.ones(ntris, dtype=np.float32) * 5.0 # Random high mu
    
    # Detector settings
    WIDTH, HEIGHT = 400, 300
    SOURCE_POS = np.array([0.0, -20.0, 10.0], dtype=np.float32)
    DET_ORIGIN = np.array([-10.0, 20.0, 15.0], dtype=np.float32) # Top Left
    DET_U      = np.array([20.0/WIDTH, 0.0, 0.0], dtype=np.float32) # Right
    DET_V      = np.array([0.0, 0.0, -15.0/HEIGHT], dtype=np.float32) # Down
    
    image_out = np.zeros(WIDTH * HEIGHT, dtype=np.float32)
    
    mf = cl.mem_flags
    buf_tris = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=tri_seq)
    buf_mu = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=mu_vals)
    buf_img = cl.Buffer(ctx, mf.WRITE_ONLY, size=image_out.nbytes)
    
    # Compile
    kernel_source = ""
    with open('cl/SpacecraftSimulation.cl', 'r') as f:
        kernel_source = f.read()
    
    prg = cl.Program(ctx, kernel_source).build()
    
    gw = (WIDTH, HEIGHT)
    lw = None
    
    prg.compute_attenuation(queue, gw, lw,
        buf_tris, buf_mu, buf_img,
        np.int32(ntris),
        cl.array.vec.make_float3(*SOURCE_POS),
        cl.array.vec.make_float3(*DET_ORIGIN),
        cl.array.vec.make_float3(*DET_U),
        cl.array.vec.make_float3(*DET_V),
        np.int32(WIDTH), np.int32(HEIGHT)
    )
    
    queue.finish()
    cl.enqueue_copy(queue, image_out, buf_img)
    
    img_2d = image_out.reshape((HEIGHT, WIDTH))
    
    print(f"Attenuation complete. Transmission range: [{img_2d.min():.4f}, {img_2d.max():.4f}]")
    
    plt.figure(figsize=(10, 8))
    plt.imshow(img_2d, cmap='gray', vmin=0, vmax=1)
    plt.colorbar(label="Transmission")
    plt.title("Ionizing Radiation Attenuation")
    plt.savefig("attenuation_preview.png")
    print("Saved attenuation_preview.png")

if __name__ == '__main__':
    V, F = generate_christmas_tree()
    run_attenuation_pipeline(V, F)
