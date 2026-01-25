import numpy as np
import pyopencl as cl
import matplotlib.pyplot as plt
from heapq import heappush, heappop
import argparse, importlib.util, os

# ==============================================================================
ap = argparse.ArgumentParser(description="Verify tile solver with fbm terrain")
ap.add_argument('--seed', type=int, default=1, help='random seed for fbm')
ap.add_argument('--octaves', type=int, default=6)
ap.add_argument('--base-freq', type=float, default=2.0)
ap.add_argument('--lacunarity', type=float, default=2.0)
ap.add_argument('--gain', type=float, default=0.5)
ap.add_argument('--warp', type=float, default=0.4)
ap.add_argument('--width', type=int, default=32)
ap.add_argument('--height', type=int, default=32)
args = ap.parse_args()

_tt_path = os.path.join(os.path.dirname(__file__), "test_terrain.py")
_spec = importlib.util.spec_from_file_location("tt_fbm", _tt_path)
tt_fbm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(tt_fbm)
fbm_value_noise = tt_fbm.fbm_value_noise

# 1. TERRAIN GENERATION (32x32)
# fbm terrain
W, H = args.width, args.height
terrain = fbm_value_noise(W, H, octaves=args.octaves, base_freq=args.base_freq, lacunarity=args.lacunarity, gain=args.gain, warp=args.warp, seed=args.seed)


# ==============================================================================
# 2. REFERENCE: GLOBAL DIJKSTRA (CPU)
# ==============================================================================
def global_dijkstra(grid, start, end):
    dist = np.full(grid.shape, np.inf)
    parent = {}
    dist[start[1], start[0]] = grid[start[1], start[0]]
    pq = [(grid[start[1], start[0]], start)]
    
    while pq:
        d, (cx, cy) = heappop(pq)
        if (cx, cy) == end: break
        if d > dist[cy, cx]: continue
        
        for dx, dy in [(-1,0), (1,0), (0,-1), (0,1), (-1,-1), (1,-1), (-1,1), (1,1)]:
            nx, ny = cx + dx, cy + dy
            if 0 <= nx < W and 0 <= ny < H:
                # Flooding Cost: max(current_path_cost, next_height)
                new_cost = max(d, grid[ny, nx]) 
                if new_cost < dist[ny, nx]:
                    dist[ny, nx] = new_cost
                    parent[(nx, ny)] = (cx, cy)
                    heappush(pq, (new_cost, (nx, ny)))
                    
    path = [end]
    while path[-1] != start:
        path.append(parent[path[-1]])
    return np.array(path)

# choose start/end
tile_size = 16
if (W % tile_size)==0 and (H % tile_size)==0:
    # derive sinks from GPU result later; temporary placeholders
    ref_path = None
else:
    ref_path = global_dijkstra(terrain, (W//4, H//4), (3*W//4, 3*H//4))

# ==============================================================================
# 3. GPU TILE SOLVER
# ==============================================================================
cl_code = """
uint pack(int x, int y) { return (uint)x | ((uint)y << 16); }

__kernel void solve_tile(
    __global const float* h, __global float* cost, __global uint* parent, __global uint* sinks) 
{
    int lx = get_local_id(0); int ly = get_local_id(1);
    int gx = get_group_id(0) * 16 + lx; int gy = get_group_id(1) * 16 + ly;
    int tid = get_group_id(1) * 2 + get_group_id(0);

    __local float lc[16][16];
    __local uint lp[16][16];
    
    // 1. Find Local Min (Sink)
    float my_h = h[gy * 32 + gx];
    // Quadratic potential to force center-preference on flats
    lc[ly][lx] = my_h + 0.005f * ((lx-8)*(lx-8) + (ly-8)*(ly-8)); 
    barrier(CLK_LOCAL_MEM_FENCE);

    __local int2 sink;
    if(lx==0 && ly==0) {
        float mv = 1e10f;
        for(int j=0; j<16; j++) for(int i=0; i<16; i++) {
            if(lc[j][i] < mv) { mv = lc[j][i]; sink = (int2)(i,j); }
        }
        sinks[tid] = pack(get_group_id(0)*16 + sink.x, get_group_id(1)*16 + sink.y);
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // 2. Intra-Tile Dijkstra (Relaxation)
    lc[ly][lx] = (lx == sink.x && ly == sink.y) ? my_h : 1e10f;
    lp[ly][lx] = pack(gx, gy);
    barrier(CLK_LOCAL_MEM_FENCE);

    for(int it=0; it<32; it++) {
        float cur_c = lc[ly][lx];
        uint cur_p = lp[ly][lx];
        for(int dy=-1; dy<=1; dy++) for(int dx=-1; dx<=1; dx++) {
            if(dx==0 && dy==0) continue;
            int nx = lx+dx; int ny = ly+dy;
            if(nx>=0 && nx<16 && ny>=0 && ny<16) {
                float cand = max(lc[ny][nx], my_h);
                if(cand < cur_c) {
                    cur_c = cand;
                    cur_p = pack(get_group_id(0)*16 + nx, get_group_id(1)*16 + ny);
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        lc[ly][lx] = cur_c;
        lp[ly][lx] = cur_p;
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    cost[gy * 32 + gx] = lc[ly][lx];
    parent[gy * 32 + gx] = lp[ly][lx];
}
"""

ctx = cl.create_some_context()
queue = cl.CommandQueue(ctx)
prg = cl.Program(ctx, cl_code).build()

h_g = cl.Buffer(ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=terrain)
c_g = cl.Buffer(ctx, cl.mem_flags.WRITE_ONLY, 4 * W * H)
p_g = cl.Buffer(ctx, cl.mem_flags.WRITE_ONLY, 4 * W * H)
s_g = cl.Buffer(ctx, cl.mem_flags.WRITE_ONLY, 4 * 4)

prg.solve_tile(queue, (W,H), (16,16), h_g, c_g, p_g, s_g)

cost_map = np.empty((H, W), dtype=np.float32)
parent_map = np.empty((H, W), dtype=np.uint32)
sinks_packed = np.empty(4, dtype=np.uint32)
cl.enqueue_copy(queue, cost_map, c_g)
cl.enqueue_copy(queue, parent_map, p_g)
cl.enqueue_copy(queue, sinks_packed, s_g)

# ==============================================================================
# 4. HIERARCHICAL GATE FINDING (The "Bridge")
# ==============================================================================
# For 2x2 tiles: examine boundary between bottom-left and bottom-right tiles
if W==32 and H==32:
    cost_left_bot = cost_map[16:32, 15]
    cost_right_bot = cost_map[16:32, 16]
    barrier = np.maximum(cost_left_bot, cost_right_bot)
    best_y_local = np.argmin(barrier)
    best_y_global = best_y_local + 16
    print(f"Hierarchical GPU Found Gate at Y: {best_y_global}")
else:
    best_y_global = H//2
    print(f"Gate scan not implemented for size {W}x{H}, using midline")

# ==============================================================================
# 5. PATH RECONSTRUCTION (Gate -> Sink)
# ==============================================================================
def unpack(val): return val & 0xFFFF, (val >> 16) & 0xFFFF

def reconstruct_gpu_path(start_x, start_y):
    path = [(start_x, start_y)]
    curr = (start_x, start_y)
    for _ in range(32):
        nxt = unpack(parent_map[curr[1], curr[0]])
        if nxt == curr: break
        curr = nxt
        path.append(curr)
    return np.array(path)

# Path from Gate Left to Sink Left
p_left = reconstruct_gpu_path(15 if W==32 else tile_size-1, best_y_global)
# Path from Gate Right to Sink Right
p_right = reconstruct_gpu_path(16 if W==32 else tile_size, best_y_global)

# ==============================================================================
# 6. VISUALIZATION
# ==============================================================================
fig, axs = plt.subplots(1, 2, figsize=(12, 6))

img0 = axs[0].imshow(terrain, cmap='terrain')
axs[0].contour(terrain, levels=10, colors='k', linewidths=0.5, alpha=0.5)
axs[0].axvline(15.5, color='red', linestyle='--', label='Tile Boundary')
axs[0].axhline(15.5, color='red', linestyle='--')
axs[0].axhline(best_y_global, color='k', linestyle=':', alpha=0.6, label=f'Pass y={best_y_global}')
if ref_path is not None:
    axs[0].plot(ref_path[:, 0], ref_path[:, 1], 'w-', linewidth=4, alpha=0.5, label='Global Reference Path')
axs[0].plot(p_left[:, 0], p_left[:, 1], 'r-', linewidth=2, label='Tile A Path (Sink -> Gate)')
axs[0].plot(p_right[:, 0], p_right[:, 1], 'b-', linewidth=2, label='Tile B Path (Gate -> Sink)')
axs[0].plot([15, 16], [best_y_global, best_y_global], 'k-', linewidth=3, label='Found Gate')
sinks = [unpack(s) for s in sinks_packed]
sx, sy = zip(*sinks)
axs[0].scatter(sx, sy, c='yellow', s=100, edgecolors='k', zorder=5, label='Tile Sinks')
axs[0].set_title(f"Terrain (seed={args.seed})")
fig.colorbar(img0, ax=axs[0], fraction=0.046, pad=0.04, label='height')
axs[0].legend()

# Cost map
img1 = axs[1].imshow(cost_map, cmap='inferno')
axs[1].axvline(15.5, color='red', linestyle='--')
axs[1].axhline(15.5, color='red', linestyle='--')
axs[1].axhline(best_y_global, color='k', linestyle=':', alpha=0.6)
axs[1].set_title("GPU intra-tile cost map")
fig.colorbar(img1, ax=axs[1], fraction=0.046, pad=0.04, label='cost')

plt.tight_layout()
plt.show()