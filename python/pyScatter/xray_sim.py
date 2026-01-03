import argparse
from pathlib import Path
import numpy as np
import pyopencl as cl
import pyopencl.array
import time
import math

parser = argparse.ArgumentParser(description="Monte-Carlo X-ray simulation")
parser.add_argument("--width",         type=int,   default=800,  help="Detector width")
parser.add_argument("--height",        type=int,   default=600,  help="Detector height")
parser.add_argument("--tubes",         type=int,   default=200,  help="Number of cylinder primitives")
parser.add_argument("--tris",          type=int,   default=200,  help="Number of triangle primitives")
parser.add_argument("--cluster-size",  type=int,   default=64,   help="Max primitives per cluster (batch size)")
parser.add_argument("--clusters",      type=int,   default=None, help="Force approximate number of clusters; overrides cluster-size")
parser.add_argument("--source-radius", type=float, default=2.0,  help="Source spot radius")
parser.add_argument("--seed",          type=int,   default=42,   help="Random seed")
args = parser.parse_args()

# --- 1. CONFIGURATION ---
WIDTH, HEIGHT = args.width, args.height
N_TUBES, N_TRIS = args.tubes, args.tris
CLUSTER_SIZE = args.cluster_size # Max primitives per cluster (matches CACHE_CAPACITY in kernel)

# Physics
SOURCE_POS = np.array([-50.0, 0.0, 0.0], dtype=np.float32)
DET_ORIGIN = np.array([50.0, -20.0, 20.0], dtype=np.float32) # Top Left
DET_U      = np.array([0.0, 40.0/WIDTH, 0.0], dtype=np.float32) # Right
DET_V      = np.array([0.0, 0.0, -40.0/HEIGHT], dtype=np.float32) # Down
SOURCE_RADIUS = np.float32(args.source_radius) # Controls blur size

ctx = cl.create_some_context()
queue = cl.CommandQueue(ctx)

# --- 2. GEOMETRY GENERATION ---
print("Generating Geometry...")
np.random.seed(args.seed)

# Structure to hold raw primitive data before packing
# type: 0=Tube, 1=Tri
# Tube: p1(start), p2(end), p3(ignored), R, Mu
# Tri:  p1(v0), p2(v1), p3(v2), Thickness, Mu
prims_raw = []

# Generate Random Tubes (Slender cylinders)
for _ in range(N_TUBES):
    # Random position in a box (-10,-10,-10) to (10,10,10)
    center = (np.random.rand(3) - 0.5) * 20.0
    axis   = (np.random.rand(3) - 0.5)
    axis   = axis / np.linalg.norm(axis)
    length = np.random.rand() * 5.0 + 1.0
    p1 = center - axis * length * 0.5
    p2 = center + axis * length * 0.5
    
    radius = 0.1 + np.random.rand() * 0.2
    mu = 0.5 + np.random.rand() * 2.0 # Attenuation coeff
    
    prims_raw.append({
        'type': 0.0,
        'p1': p1, 'p2': p2, 'p3': np.array([0,0,0]),
        'param': radius, 'mu': mu,
        'center': center # Used for clustering
    })

# Generate Random Triangles (Plates)
for _ in range(N_TRIS):
    center = (np.random.rand(3) - 0.5) * 15.0 # Slightly tighter cluster
    # Create random triangle near center
    v0 = center + (np.random.rand(3)-0.5)*3.0
    v1 = center + (np.random.rand(3)-0.5)*3.0
    v2 = center + (np.random.rand(3)-0.5)*3.0
    
    thickness = 0.05
    mu = 5.0 # High attenuation (metal)
    
    prims_raw.append({
        'type': 1.0,
        'p1': v0, 'p2': v1, 'p3': v2,
        'param': thickness, 'mu': mu,
        'center': (v0+v1+v2)/3.0
    })

# --- 3. CLUSTERING & PACKING ---
print("Clustering & Packing...")

# Simple Spatial Sort for mock clustering (Sort by Z then X)
# A real engine would use BVH or Octree construction here.
prims_raw.sort(key=lambda p: (p['center'][2], p['center'][0]))

# Arrays for GPU
gpu_prims = []    # The float4 stream
gpu_spheres = []  # Cluster bounds
gpu_ranges = []   # {start, count}

current_prim_idx = 0

# Chunk into clusters
num_prims = len(prims_raw)
effective_cluster_size = max(1, math.ceil(num_prims / args.clusters)) if args.clusters else CLUSTER_SIZE
for i in range(0, num_prims, effective_cluster_size):
    chunk = prims_raw[i : min(i + effective_cluster_size, num_prims)]
    
    # Calculate Bounding Sphere for Chunk
    centers = np.array([p['center'] for p in chunk])
    min_pt = np.min(centers, axis=0)
    max_pt = np.max(centers, axis=0)
    b_center = (min_pt + max_pt) / 2.0
    # Radius covers all primitive centers + max primitive size approx
    dists = np.linalg.norm(centers - b_center, axis=1)
    b_radius = np.max(dists) + 3.0 # +3.0 padding for object extent
    
    # Pack Primitives
    for p in chunk:
        # Row 1: {x,y,z, type}
        gpu_prims.append([p['p1'][0], p['p1'][1], p['p1'][2], p['type']])
        # Row 2: {x,y,z, param} (param = R or Thickness)
        gpu_prims.append([p['p2'][0], p['p2'][1], p['p2'][2], p['param']])
        # Row 3: {x,y,z, mu}
        gpu_prims.append([p['p3'][0], p['p3'][1], p['p3'][2], p['mu']])
        
    # Meta Data
    gpu_spheres.append([b_center[0], b_center[1], b_center[2], b_radius])
    gpu_ranges.append([current_prim_idx, len(chunk)])
    
    current_prim_idx += len(chunk)

# Convert to numpy format for OpenCL
# NOTE: float4 in OpenCL is 16 bytes. 
np_prims = np.array(gpu_prims, dtype=np.float32).flatten()
np_spheres = np.array(gpu_spheres, dtype=np.float32).flatten()
np_ranges = np.array(gpu_ranges, dtype=np.int32).flatten() # int2

# --- 4. GPU SETUP ---
print(f"Total Prims: {num_prims}, Clusters: {len(gpu_ranges)}")

mf = cl.mem_flags
buf_prims   = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=np_prims)
buf_spheres = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=np_spheres)
buf_ranges  = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=np_ranges)
buf_out     = cl.Buffer(ctx, mf.WRITE_ONLY, size=WIDTH * HEIGHT * 4) # float


kernel_path = Path(__file__).parent / "cl" / "cl" / "Xray.cl"
kernel_source = kernel_path.read_text()
# Compile Kernel (Assuming kernel_source is defined above)
prg = cl.Program(ctx, kernel_source).build()

# --- 5. EXECUTION ---
print("Running Kernel...")
start_time = time.time()

# Global work size must be multiple of local size (16,16)
gw = (math.ceil(WIDTH/16)*16, math.ceil(HEIGHT/16)*16)
lw = (16, 16)

prg.xray_simulation(queue, gw, lw,
    buf_prims, buf_spheres, buf_ranges, buf_out,
    np.int32(len(gpu_ranges)),
    cl.array.vec.make_float3(*SOURCE_POS),
    np.float32(SOURCE_RADIUS),
    cl.array.vec.make_float3(*DET_ORIGIN),
    cl.array.vec.make_float3(*DET_U),
    cl.array.vec.make_float3(*DET_V),
    np.int32(WIDTH)
)

queue.finish()
print(f"Execution time: {time.time() - start_time:.4f} s")

# --- 6. VISUALIZATION ---
result = np.empty(WIDTH * HEIGHT, dtype=np.float32)
cl.enqueue_copy(queue, result, buf_out)
result = result.reshape((HEIGHT, WIDTH))




try:
    import matplotlib.pyplot as plt
    plt.figure(figsize=(10, 8))
    plt.imshow(result, cmap='gray', vmin=0, vmax=1)
    plt.title("Monte Carlo Attenuation (X-Ray Projection)")
    plt.colorbar(label="Transmission")
    plt.show()
except ImportError:
    print("Matplotlib not installed. Saving to raw file 'result.raw'")
    result.tofile("result.raw")



## RUN AS
# export PYOPENCL_CTX='1'; python xray_sim.py