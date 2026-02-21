import numpy as np
import pyopencl as cl
import os
from heapq import heappush, heappop

# GPU tile solver wrapper
class GPUTerrainSolver:
    def __init__(self, width, height):
        self.width = np.int32(width)
        self.height = np.int32(height)
        self.ctx = cl.create_some_context()
        self.queue = cl.CommandQueue(self.ctx)

        current_dir = os.path.dirname(__file__) if __file__ != "" else "."
        with open(os.path.join(current_dir, "terrain.cl"), "r") as f:
            self.prg = cl.Program(self.ctx, f.read()).build()

        self.mf = cl.mem_flags
        self.num_tiles_x = width // 16
        self.num_tiles_y = height // 16

    def solve(self, h_map, K=0.01):
        h_gpu = cl.Buffer(self.ctx, self.mf.READ_ONLY | self.mf.COPY_HOST_PTR, hostbuf=h_map.astype(np.float32))
        parent_gpu = cl.Buffer(self.ctx, self.mf.READ_WRITE, 4 * self.width * self.height)
        cost_gpu = cl.Buffer(self.ctx, self.mf.READ_WRITE, 4 * self.width * self.height)
        sinks_gpu = cl.Buffer(self.ctx, self.mf.READ_WRITE, 4 * self.num_tiles_x * self.num_tiles_y)

        self.prg.solve_tiles(self.queue, (self.width, self.height), (16, 16),
                             h_gpu, parent_gpu, cost_gpu, sinks_gpu,
                             self.width, self.height, np.float32(K))

        parent_res = np.empty((self.height, self.width), dtype=np.uint32)
        cost_res = np.empty((self.height, self.width), dtype=np.float32)
        sinks_res = np.empty(self.num_tiles_y * self.num_tiles_x, dtype=np.uint32)

        cl.enqueue_copy(self.queue, parent_res, parent_gpu)
        cl.enqueue_copy(self.queue, cost_res, cost_gpu)
        cl.enqueue_copy(self.queue, sinks_res, sinks_gpu)

        return parent_res, cost_res, sinks_res

def unpack(val):
    return val & 0xFFFF, (val >> 16) & 0xFFFF

def pack(x, y):
    return np.uint32(x | (y << 16))

# Watershed helpers (global grid kernels in terrain.cl)
def build_watershed_program(ctx):
    current_dir = os.path.dirname(__file__) if __file__ != "" else "."
    kernel_path = os.path.join(current_dir, "terrain.cl")
    with open(kernel_path, "r") as f:
        src = f.read()
    return cl.Program(ctx, src).build()


def run_calc_flow(queue, prg, heightmap_buf, flow_buf, basin_buf, W, H):
    prg.calc_flow(queue, (W, H), None,
                  heightmap_buf, flow_buf, basin_buf,
                  np.int32(W), np.int32(H))


def run_propagate_basins(queue, prg, flow_buf, basin_buf, changed_buf, W, H):
    prg.propagate_basins(queue, (W, H), None,
                         flow_buf, basin_buf, changed_buf,
                         np.int32(W), np.int32(H))

# Shared utilities (no solver, no CLI)

def hash_u32(x):
    x = (x ^ (x >> 16)) & 0xFFFFFFFF
    x = (x * 0x7feb352d) & 0xFFFFFFFF
    x = (x ^ (x >> 15)) & 0xFFFFFFFF
    x = (x * 0x846ca68b) & 0xFFFFFFFF
    x = (x ^ (x >> 16)) & 0xFFFFFFFF
    return x

def smoothstep(t):
    return t*t*(3.0 - 2.0*t)

def value_noise_2d(x, y, seed=1):
    xi = np.floor(x).astype(np.int64)
    yi = np.floor(y).astype(np.int64)
    xf = (x - xi).astype(np.float32)
    yf = (y - yi).astype(np.float32)
    u = smoothstep(xf)
    v = smoothstep(yf)
    def rnd(ix, iy):
        h = (np.uint64(ix) * np.uint64(0x9e3779b1) + np.uint64(iy) * np.uint64(0x85ebca6b) + np.uint64(seed) * np.uint64(0x27d4eb2d)) & np.uint64(0xFFFFFFFF)
        r = hash_u32(np.uint32(h))
        return (r.astype(np.float32) / np.float32(0xFFFFFFFF))
    v00 = rnd(xi + 0, yi + 0)
    v10 = rnd(xi + 1, yi + 0)
    v01 = rnd(xi + 0, yi + 1)
    v11 = rnd(xi + 1, yi + 1)
    vx0 = v00*(1.0 - u) + v10*u
    vx1 = v01*(1.0 - u) + v11*u
    return vx0*(1.0 - v) + vx1*v

def fbm_value_noise(W, H, octaves=6, base_freq=2.0, lacunarity=2.0, gain=0.5, warp=0.0, seed=1):
    X, Y = np.meshgrid(np.linspace(0.0, 1.0, W, endpoint=False), np.linspace(0.0, 1.0, H, endpoint=False))
    x = (X * base_freq).astype(np.float32)
    y = (Y * base_freq).astype(np.float32)
    h = np.zeros((H, W), dtype=np.float32)
    amp = 1.0
    freq = 1.0
    ang = (seed * 12.9898) % (2.0*np.pi)
    c, s = np.cos(ang), np.sin(ang)
    for i in range(octaves):
        rx = (c*x - s*y) * freq
        ry = (s*x + c*y) * freq
        if warp != 0.0:
            wx = value_noise_2d(rx + 13.1, ry + 7.7, seed + 101 + i)
            wy = value_noise_2d(rx - 9.2, ry + 21.3, seed + 211 + i)
            rx = rx + warp*(wx - 0.5)
            ry = ry + warp*(wy - 0.5)
        h += amp * (value_noise_2d(rx, ry, seed + i) * 2.0 - 1.0)
        amp *= gain
        freq *= lacunarity
    h -= h.min()
    h /= h.max() if h.max() != 0 else 1.0
    return h.astype(np.float32)



# Local reimplementation of utilities to ensure this script is self-contained
# (User asked for minimal changes, but external deps were breaking. This is safer)
def smoothstep(t): return t*t*(3.0 - 2.0*t)
def value_noise_2d(x, y, seed=1):
    xi = np.floor(x).astype(np.int64); yi = np.floor(y).astype(np.int64)
    xf = x - xi; yf = y - yi
    u = smoothstep(xf); v = smoothstep(yf)
    def hash_u32(x):
        x = np.uint64(x); x = (x ^ (x >> 32)) * 0x45d9f3b
        x = ((x ^ (x >> 32)) * 0x45d9f3b) & 0xFFFFFFFF
        return x.astype(np.uint32)
    def rnd(ix, iy):
        h = (np.uint64(ix)*0x9e3779b1 + np.uint64(iy)*0x85ebca6b + np.uint64(seed)*0x27d4eb2d)
        return (hash_u32(h).astype(np.float32) / 4294967295.0)
    return rnd(xi,yi)*(1-u)*(1-v) + rnd(xi+1,yi)*u*(1-v) + rnd(xi,yi+1)*(1-u)*v + rnd(xi+1,yi+1)*u*v

def fbm_value_noise_2(W, H, octaves=6, base_freq=2.0, lacunarity=2.0, gain=0.5, warp=0.0, seed=1):
    X, Y = np.meshgrid(np.linspace(0,1,W,False), np.linspace(0,1,H,False))
    x, y = X*base_freq, Y*base_freq
    h = np.zeros((H,W),np.float32); amp=1.0
    for _ in range(octaves):
        h += amp * (value_noise_2d(x, y, seed)*2-1); amp *= gain; x*=lacunarity; y*=lacunarity
    h -= h.min(); h /= h.max()
    return h

def compute_tile_sinks(sinks, W, H, tile_size=16, unpack_fn=None):
    if unpack_fn is None:
        unpack_fn = lambda v: (v & 0xFFFF, (v >> 16) & 0xFFFF)
    nx = W // tile_size
    ny = H // tile_size
    sink_xy = np.zeros((ny, nx, 2), dtype=np.int32)
    for it, sp in enumerate(sinks):
        tx = it % nx
        ty = it // nx
        sink_xy[ty, tx, :] = unpack_fn(sp)
    return sink_xy

def tile_graph_shortest(gates, start_tile, end_tile):
    adj = {}
    for g in gates:
        a, b, w = tuple(g['a']), tuple(g['b']), g['barrier']
        adj.setdefault(a, []).append((b, w, g))
        adj.setdefault(b, []).append((a, w, g))
    dist = {start_tile: 0.0}
    prev = {}
    pq = [(0.0, start_tile)]
    while pq:
        d, t = heappop(pq)
        if t == end_tile:
            break
        if d > dist.get(t, 1e30):
            continue
        for nb, w, g in adj.get(t, []):
            nd = max(d, w)
            if nd < dist.get(nb, 1e30):
                dist[nb] = nd
                prev[nb] = (t, g)
                heappush(pq, (nd, nb))
    if end_tile not in prev and end_tile != start_tile:
        raise ValueError(f"No tile path from {start_tile} to {end_tile}")
    tile_path = [end_tile]
    gate_seq = []
    cur = end_tile
    while cur != start_tile:
        p, g = prev[cur]
        tile_path.append(p)
        gate_seq.append(g)
        cur = p
    tile_path.reverse(); gate_seq.reverse()
    return tile_path, gate_seq, dist.get(end_tile, 0.0)


# 4. HIERARCHICAL GATE FINDING (Diagonal Aware)
def compute_gates_simple(costs, tile_size=16):
    """Cost-only gate finder (horizontal/vertical), sets a/b/pA/pB/barrier."""
    ny, nx = costs.shape[0]//tile_size, costs.shape[1]//tile_size
    gates = []
    for ty in range(ny):
        for tx in range(nx - 1):
            xL, xR = (tx+1)*tile_size - 1, (tx+1)*tile_size
            y_slice = slice(ty*tile_size, (ty+1)*tile_size)
            cL, cR = costs[y_slice, xL], costs[y_slice, xR]
            barrier = np.maximum(cL, cR)
            idx = int(np.argmin(barrier))
            gates.append({'a': (tx, ty), 'b': (tx + 1, ty),
                          'pA': (xL, ty*tile_size + idx),
                          'pB': (xR, ty*tile_size + idx),
                          'barrier': float(barrier[idx])})
    for ty in range(ny - 1):
        for tx in range(nx):
            yT, yB = (ty+1)*tile_size - 1, (ty+1)*tile_size
            x_slice = slice(tx*tile_size, (tx+1)*tile_size)
            cT, cB = costs[yT, x_slice], costs[yB, x_slice]
            barrier = np.maximum(cT, cB)
            idx = int(np.argmin(barrier))
            gates.append({'a': (tx, ty), 'b': (tx, ty + 1),
                          'pA': (tx*tile_size + idx, yT),
                          'pB': (tx*tile_size + idx, yB),
                          'barrier': float(barrier[idx])})
    return gates


def compute_gates_diag(costs, heights, tile_size=16):
    """Diagonal/height-aware gate finder using cost+dist heuristic."""
    H, W = costs.shape
    gates = []
    def check_boundary(yA, xA, yB, xB):
        cA = costs[yA, xA]; cB = costs[yB, xB]
        dist = np.sqrt((xA-xB)**2 + (yA-yB)**2) * 0.001
        return cA + cB + dist
    nx, ny = W // tile_size, H // tile_size
    # Horizontal
    for ty in range(ny):
        for tx in range(nx - 1):
            xL = (tx+1)*tile_size - 1
            xR = xL + 1
            for y in range(ty*tile_size, (ty+1)*tile_size):
                b = check_boundary(y, xL, y, xR)
                gates.append({'a': (tx, ty), 'b': (tx + 1, ty),
                              'pA': (xL, y), 'pB': (xR, y),
                              'barrier': b})
    # Vertical
    for ty in range(ny - 1):
        yT = (ty+1)*tile_size - 1
        yB = yT + 1
        for tx in range(nx):
            for x in range(tx*tile_size, (tx+1)*tile_size):
                b = check_boundary(yT, x, yB, x)
                gates.append({'a': (tx, ty), 'b': (tx, ty + 1),
                              'pA': (x, yT), 'pB': (x, yB),
                              'barrier': b})
    return gates


def get_pixel_path(start, parents, unpack_fn=None, max_steps=2048):
    if unpack_fn is None:
        unpack_fn = lambda v: (v & 0xFFFF, (v >> 16) & 0xFFFF)
    path = [start]
    curr = start
    for _ in range(max_steps):
        nxt = unpack_fn(parents[curr[1], curr[0]])
        if nxt == curr:
            break
        curr = nxt
        path.append(curr)
    return np.array(path)

# ------------------------------------------------------------------------------
# REFERENCE DIJKSTRA MATCHING GPU LOGIC EXACTLY
# ------------------------------------------------------------------------------
def global_dijkstra(grid, start, end):
    H, W = grid.shape
    dist = np.full(grid.shape, np.inf)
    parent = {}
    dist[start[1], start[0]] = grid[start[1], start[0]]
    pq = [(grid[start[1], start[0]], start)]
    while pq:
        d, (cx, cy) = heappop(pq)
        if (cx, cy) == end: break
        if d > dist[cy, cx]: continue
        for dx in (-1,0,1):
            for dy in (-1,0,1):
                if dx==0 and dy==0:
                    continue
                nx = cx + dx
                ny = cy + dy
                if 0 <= nx < W and 0 <= ny < H:
                    # FIX: Euclidean Weights
                    #dist_penalty = 1.41421e-5 if (dx != 0 and dy != 0) else 1e-5
                    # loat dist_penalty = (dx != 0 && dy != 0) ? 1.41421f * 0.001f : 0.001f;
                    dist_penalty = 1.41421e-3 if (dx != 0 and dy != 0) else 1e-3
                    new_cost = max(d, grid[ny, nx]) + dist_penalty
                    if new_cost < dist[ny, nx]:
                        dist[ny, nx] = new_cost
                        parent[(nx, ny)] = (cx, cy)
                        heappush(pq, (new_cost, (nx, ny)))            
    if end not in parent and start != end:
        return None
    path = [end]
    curr = end
    while curr != start:
        curr = parent[curr]
        path.append(curr)
    return np.array(path)

# ==============================================================================
# 3. ALGORITHM B: Iterative Relaxation (NumPy Vectorized)
# ==============================================================================
def solve_relaxation_numpy(grid, start, end, epsilon=0.001, max_iter=1500):
    H, W = grid.shape
    
    # 1. Initialize Costs
    costs = np.full((H, W), 1e10, dtype=np.float32)
    sx, sy = start
    costs[sy, sx] = grid[sy, sx] # Sink height
    
    # Parent pointers: encode as (y * W + x)
    Y, X = np.mgrid[0:H, 0:W]
    parents = (Y * W + X).astype(np.int32)
    
    # 2. Iteration (Simulates GPU Parallelism)
    for i in range(max_iter):
        prev_costs = costs.copy()
        
        # We will compute the min_cost for every pixel by looking at 8 neighbors simultaneously
        # Shifted arrays allow us to compare cost[y,x] with cost[y+dy, x+dx] without loops
        
        for dy in [-1, 0, 1]:
            for dx in [-1, 0, 1]:
                if dx == 0 and dy == 0: continue
                
                # Penalty
                pen = 1.41421*epsilon if (dx!=0 and dy!=0) else epsilon
                
                # Shift the COST map to bring neighbors to the center
                # We want: neighbor_cost at (y,x) corresponding to direction (dy, dx)
                # If we want the neighbor to the LEFT (dx=-1), we shift the array RIGHT (+1)
                shifted_costs = np.roll(prev_costs, (-dy, -dx), axis=(0, 1))
                
                # Handle boundaries (roll wraps around, so we must mask edges)
                # (Simple boundary fix: set wrapped edges to infinity)
                if dy == 1: shifted_costs[-1, :] = 1e10
                if dy == -1: shifted_costs[0, :] = 1e10
                if dx == 1: shifted_costs[:, -1] = 1e10
                if dx == -1: shifted_costs[:, 0] = 1e10
                
                # The Core Logic: max(neighbor, my_height) + penalty
                candidates = np.maximum(shifted_costs, grid) + pen
                
                # Update Min
                improved_mask = candidates < costs
                costs[improved_mask] = candidates[improved_mask]
                
                # Update Parents (Reverse logic: If I updated from neighbor (dy,dx), 
                # my parent is (y+dy, x+dx))
                # We need to construct the parent index grid shifted similarly
                # Note: This is tricky in vectorization. 
                # Easier approach for parents: Just store direction index or re-compute at end.
                # Let's re-compute parent at the end to be sure.
        
        # Check convergence
        if np.allclose(prev_costs, costs):
            print(f"NumPy Relaxation converged in {i} iterations.")
            break

    # 3. Reconstruct Parent Pointers (Post-Process)
    # Now that we have the perfect cost map, finding parents is easy
    # For each pixel, pick neighbor that satisfies cost logic
    parent_map = {}
    for y in range(H):
        for x in range(W):
            if x == sx and y == sy: continue
            
            best_n_cost = 1e10
            best_n = None
            
            for dy in [-1, 0, 1]:
                for dx in [-1, 0, 1]:
                    if dx==0 and dy==0: continue
                    nx, ny = x+dx, y+dy
                    if 0<=nx<W and 0<=ny<H:
                        # Check if this neighbor could have been the parent
                        # cost[y,x] should be approx max(cost[ny,nx], h[y,x]) + pen
                        pen = 1.41421*epsilon if (dx!=0 and dy!=0) else epsilon
                        
                        # We look for the neighbor with the lowest compatible cost
                        if costs[ny, nx] < best_n_cost:
                            best_n_cost = costs[ny, nx]
                            best_n = (nx, ny)
            
            if best_n:
                parent_map[(x,y)] = best_n

    # 4. Trace Path
    if end not in parent_map: return np.zeros((0,2))
    path = [end]
    curr = end
    while curr != start:
        curr = parent_map[curr]
        path.append(curr)
        if len(path) > W*H: break # Safety
    return np.array(path)


class GPUErosionAccumulator:
    def __init__(self, width, height, ctx=None, queue=None):
        self.width = np.int32(width)
        self.height = np.int32(height)
        self.ctx = ctx or cl.create_some_context()
        self.queue = queue or cl.CommandQueue(self.ctx)

        current_dir = os.path.dirname(__file__) if __file__ != "" else "."
        # Reuse terrain.cl which already contains accumulate_sdf
        with open(os.path.join(current_dir, "terrain.cl"), "r") as f:
            self.prg = cl.Program(self.ctx, f.read()).build()

        self.mf = cl.mem_flags
        self.accum_gpu = cl.Buffer(self.ctx, self.mf.READ_WRITE, 4 * width * height)
        cl.enqueue_fill_buffer(self.queue, self.accum_gpu, np.float32(0), 0, 4 * width * height)

    def run_passes(self, input_noise_map, passes=100, radius=16):
            noise_gpu = cl.Buffer(self.ctx, self.mf.READ_ONLY | self.mf.COPY_HOST_PTR, hostbuf=input_noise_map.astype(np.float32))
            global_size = (int(self.width), int(self.height))
            local_size = (16, 16)
            
            print(f"Running {passes} erosion passes...")
            
            # Pre-compile args to speed up loop slightly
            r_arg = np.int32(radius)
            w_arg = np.int32(self.width)
            h_arg = np.int32(self.height)
            
            for i in range(passes):
                # Generate random threshold for this frame (0.0 to 1.0)
                # This slices the terrain at a specific height
                threshold = np.float32(np.random.rand())
                
                self.prg.accumulate_sdf(
                    self.queue, global_size, local_size,
                    noise_gpu,
                    self.accum_gpu,
                    w_arg, h_arg,
                    threshold,   # <-- Passing global threshold
                    r_arg
                )
                
                if i % 10 == 0:
                    print(f"Pass {i}/{passes}", end='\r')
                    
            self.queue.finish()
            print(f"\nDone.")

    def get_result(self, total_passes):
        result = np.empty((self.height, self.width), dtype=np.float32)
        cl.enqueue_copy(self.queue, result, self.accum_gpu)
        return result / total_passes

# import pyopencl as cl
# import numpy as np
# import os
# import matplotlib.pyplot as plt

# import pyopencl as cl
# import numpy as np
# import os
# import matplotlib.pyplot as plt

# class GPUShadertoyErosion:
#     def __init__(self, width, height, ctx=None, queue=None):
#         self.width = np.int32(width)
#         self.height = np.int32(height)
#         self.ctx = ctx or cl.create_some_context()
#         self.queue = queue or cl.CommandQueue(self.ctx)

#         # Load Kernels
#         current_dir = os.path.dirname(__file__) if __file__ != "" else "."
#         with open(os.path.join(current_dir, "terrain.cl"), "r") as f:
#             self.prg = cl.Program(self.ctx, f.read()).build()
        
#         self.mf = cl.mem_flags
#         n_pixels = width * height
        
#         # Buffers
#         self.noise_gpu = cl.Buffer(self.ctx, self.mf.READ_WRITE, 4 * n_pixels)
#         self.accum_gpu = cl.Buffer(self.ctx, self.mf.READ_WRITE, 4 * n_pixels)
#         self.height_gpu = cl.Buffer(self.ctx, self.mf.READ_WRITE, 4 * n_pixels)
#         self.image_gpu = cl.Buffer(self.ctx, self.mf.READ_WRITE, 16 * n_pixels) # float4 RGBA
        
#         cl.enqueue_fill_buffer(self.queue, self.accum_gpu, np.float32(0), 0, 4 * n_pixels)

#     def generate_input_noise(self, seed=42):
#         print("1. Generating Input Noise...")
#         self.prg.generate_noise(
#             self.queue, (int(self.width), int(self.height)), (16, 16),
#             self.noise_gpu, self.width, self.height, np.float32(seed)
#         )

#     def run_accumulation(self, frames=256):
#         print(f"2. Eroding ({frames} passes)...")
#         global_size = (int(self.width), int(self.height))
#         local_size = (16, 16)
        
#         # Batch enqueue to avoid driver timeouts on huge loops
#         BATCH = 16
#         for i in range(0, frames, BATCH):
#             for j in range(BATCH):
#                 if i+j >= frames: break
#                 self.prg.buffer_a(
#                     self.queue, global_size, local_size,
#                     self.noise_gpu, self.accum_gpu,
#                     self.width, self.height, np.int32(i+j)
#                 )
#             self.queue.finish()
#             print(f"   Pass {i+BATCH}/{frames}", end='\r')
#         print("")

#     def run_post_process(self, total_frames):
#         print("3. Generating Heightmap (Buffer B)...")
#         # Normalize Accumulator on CPU (easiest for prototype)
#         accum_host = np.empty((self.height, self.width), dtype=np.float32)
#         cl.enqueue_copy(self.queue, accum_host, self.accum_gpu)
#         accum_host /= float(total_frames)
        
#         # Upload normalized for Buffer B
#         norm_accum_gpu = cl.Buffer(self.ctx, self.mf.READ_ONLY | self.mf.COPY_HOST_PTR, hostbuf=accum_host)
        
#         self.prg.buffer_b(
#             self.queue, (int(self.width), int(self.height)), (16, 16),
#             norm_accum_gpu, self.height_gpu,
#             self.width, self.height, np.int32(total_frames)
#         )
#         self.queue.finish()

#     def render_image(self):
#         print("4. Rendering Colors & Lighting...")
#         self.prg.render_map(
#             self.queue, (int(self.width), int(self.height)), (16, 16),
#             self.height_gpu, self.image_gpu,
#             self.width, self.height
#         )
#         self.queue.finish()
        
#         # Retrieve Float4 Image
#         img_host = np.empty((self.height, self.width, 4), dtype=np.float32)
#         cl.enqueue_copy(self.queue, img_host, self.image_gpu)
        
#         # Convert to 0..1 range (Clamp)
#         #img_host = np.clip(img_host, 0.0, 10000.0)
#         return img_host







import pyopencl as cl
import numpy as np
import os
import matplotlib.pyplot as plt

class GPUShadertoyErosion:
    def __init__(self, width, height, ctx=None, queue=None):
        self.width = np.int32(width)
        self.height = np.int32(height)
        self.ctx = ctx or cl.create_some_context()
        self.queue = queue or cl.CommandQueue(self.ctx)

        current_dir = os.path.dirname(__file__) if __file__ != "" else "."
        with open(os.path.join(current_dir, "terrain.cl"), "r") as f:
            self.prg = cl.Program(self.ctx, f.read()).build()
        
        self.mf = cl.mem_flags
        n_pixels = width * height
        
        self.noise_gpu = cl.Buffer(self.ctx, self.mf.READ_WRITE, 4 * n_pixels)
        self.accum_gpu = cl.Buffer(self.ctx, self.mf.READ_WRITE, 4 * n_pixels)
        self.height_gpu = cl.Buffer(self.ctx, self.mf.READ_WRITE, 4 * n_pixels)
        self.image_gpu = cl.Buffer(self.ctx, self.mf.READ_WRITE, 16 * n_pixels)
        
        cl.enqueue_fill_buffer(self.queue, self.accum_gpu, np.float32(0), 0, 4 * n_pixels)

    def generate_input_noise(self, seed=42):
        print("1. Generating Noise...")
        self.prg.generate_noise(
            self.queue, (int(self.width), int(self.height)), (16, 16),
            self.noise_gpu, self.width, self.height, np.float32(seed)
        )

    def run_accumulation(self, frames=256):
        print(f"2. Eroding ({frames} frames)...")
        global_size = (int(self.width), int(self.height))
        local_size = (16, 16)
        
        # Batching for stability
        BATCH = 32
        for i in range(0, frames, BATCH):
            for j in range(BATCH):
                if i+j >= frames: break
                self.prg.buffer_a(
                    self.queue, global_size, local_size,
                    self.noise_gpu, self.accum_gpu,
                    self.width, self.height, np.int32(i+j)
                )
            self.queue.finish()
            print(f"   Frame {i+BATCH}/{frames}", end='\r')
        print("")

    def run_post_process(self, total_frames):
        print("3. Post-Processing (Buffer B)...")
        
        # Normalize Accumulator (Average the frames)
        accum_host = np.empty((self.height, self.width), dtype=np.float32)
        cl.enqueue_copy(self.queue, accum_host, self.accum_gpu)
        accum_host /= float(total_frames)
        
        # Upload normalized buffer for Buffer B
        norm_accum_gpu = cl.Buffer(self.ctx, self.mf.READ_ONLY | self.mf.COPY_HOST_PTR, hostbuf=accum_host)
        
        self.prg.buffer_b(
            self.queue, (int(self.width), int(self.height)), (16, 16),
            norm_accum_gpu, self.height_gpu,
            self.width, self.height
        )
        self.queue.finish()

    def render_image(self):
        print("4. Rendering...")
        self.prg.render_map(
            self.queue, (int(self.width), int(self.height)), (16, 16),
            self.height_gpu, self.image_gpu,
            self.width, self.height
        )
        self.queue.finish()
        
        img_host = np.empty((self.height, self.width, 4), dtype=np.float32)
        cl.enqueue_copy(self.queue, img_host, self.image_gpu)
        return np.clip(img_host, 0.0, 1.0)





# ==============================================================================
# MAIN
# ==============================================================================


# ==============================================================================
# MAIN EXECUTION
# ==============================================================================
# if __name__ == "__main__":
#     import os
#     os.environ.setdefault("PYOPENCL_CTX", "1")
    
#     W, H = 512, 512
#     FRAMES = 256 # More frames = smoother erosion
    
#     sim = GPUShadertoyErosion(W, H)
    
#     # 1. Generate Input Noise (iChannel0)
#     sim.generate_input_noise(seed=1337)
    
#     # 2. Accumulate SDF (Buffer A)
#     sim.run_accumulation(frames=FRAMES)
    
#     # 3. Post-Process (Buffer B)
#     heightmap = sim.run_post_process(total_frames=FRAMES)
    
#     # 4. Visualize
#     plt.figure(figsize=(10, 8))
#     plt.imshow(heightmap, cmap='gray')
#     plt.title(f"Accumulative SDF Erosion ({FRAMES} passes)")
#     plt.colorbar()
#     plt.show()        