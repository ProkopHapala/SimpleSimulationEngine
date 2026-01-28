import numpy as np
import pyopencl as cl
import matplotlib.pyplot as plt
import argparse, importlib.util, os, heapq
os.environ.setdefault("PYOPENCL_CTX", "1")

# Robustly import terrain_utils from the sibling directory or similar
# Assuming the user has the 'terrain.py' (wrapper) and 'terrain.cl' in the same folder as this script
# We will use the local files if possible
import sys
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)

try:
    from terrain import GPUTerrainSolver, unpack, pack
except ImportError:
    # Fallback to the user's specific path if local import fails
    _utils_path = os.path.join(current_dir, "../terrain_ocl_2/terrain.py")
    if os.path.exists(_utils_path):
        _spec = importlib.util.spec_from_file_location("terrain_utils", _utils_path)
        terrain_utils = importlib.util.module_from_spec(_spec)
        _spec.loader.exec_module(terrain_utils)
        GPUTerrainSolver = terrain_utils.GPUTerrainSolver
        unpack = terrain_utils.unpack
        pack = terrain_utils.pack
    else:
        raise ImportError("Could not find terrain.py in local folder or ../terrain_ocl_2/")

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

def fbm_value_noise(W, H, octaves=6, base_freq=2.0, lacunarity=2.0, gain=0.5, warp=0.0, seed=1):
    X, Y = np.meshgrid(np.linspace(0,1,W,False), np.linspace(0,1,H,False))
    x, y = X*base_freq, Y*base_freq
    h = np.zeros((H,W),np.float32); amp=1.0
    for _ in range(octaves):
        h += amp * (value_noise_2d(x, y, seed)*2-1); amp *= gain; x*=lacunarity; y*=lacunarity
    h -= h.min(); h /= h.max()
    return h

# Shared tile graph minimax path helper (imported to avoid duplication)
_tg_spec = importlib.util.spec_from_file_location("terrain_utils_tile", os.path.join(current_dir, "../terrain_ocl_2/terrain.py"))
_tg_mod = importlib.util.module_from_spec(_tg_spec)
_tg_spec.loader.exec_module(_tg_mod)
tile_graph_shortest = _tg_mod.tile_graph_shortest
compute_gates = _tg_mod.compute_gates
solve_relaxation_numpy = getattr(_tg_mod, "solve_relaxation_numpy", None)

def compute_tile_sinks(sinks_packed, W, H, tile_size=16, unpack_fn=unpack):
    nx, ny = W//tile_size, H//tile_size
    grid = np.zeros((ny, nx, 2), dtype=np.int32)
    for i, s in enumerate(sinks_packed):
        grid[i//nx, i%nx] = unpack_fn(s)
    return grid

def global_dijkstra(grid, start, end):
    H, W = grid.shape
    dist = np.full_like(grid, np.inf, dtype=np.float32)
    parent = {}
    dist[start[1], start[0]] = grid[start[1], start[0]]
    pq = [(float(grid[start[1], start[0]]), start)]
    while pq:
        d, (cx, cy) = heapq.heappop(pq)
        if (cx, cy) == end: break
        if d > dist[cy, cx]: continue
        for dx in (-1,0,1):
            for dy in (-1,0,1):
                if dx==0 and dy==0: continue
                nx = cx + dx; ny = cy + dy
                if 0 <= nx < W and 0 <= ny < H:
                    cand = max(d, grid[ny, nx])
                    if cand < dist[ny, nx]:
                        dist[ny, nx] = cand
                        parent[(nx, ny)] = (cx, cy)
                        heapq.heappush(pq, (cand, (nx, ny)))
    path = [end]
    while path[-1] != start:
        path.append(parent[path[-1]])
    path.reverse()
    return np.array(path), dist[end[1], end[0]]

def get_pixel_path(start, parents, unpack_fn, max_steps=2000):
    path = [start]; curr = start
    for _ in range(max_steps):
        nxt = unpack_fn(parents[curr[1], curr[0]])
        if nxt == curr: break
        curr = nxt; path.append(curr)
    return np.array(path)

# ------------------------------------------------------------------------------

def log_ref(name, path, cost, verb=0):
    pts = [(int(p[0]), int(p[1])) for p in path]
    print(f"{name}: cost={cost:.6f}, start={pts[0]}, end={pts[-1]}, len={len(pts)}")
    if verb >= 2:
        if verb >= 3 and len(pts) > 50:
            sample = pts[:10] + ['...'] + pts[-10:]
            print(f"{name} full (truncated): {sample}")
        else:
            print(f"{name} full: {pts}")

ap = argparse.ArgumentParser(description="Verify tile solver with fbm terrain")
ap.add_argument('--seed',       type=int,   default=9, help='random seed for fbm')
ap.add_argument('--octaves',    type=int,   default=6)
ap.add_argument('--base-freq',  type=float, default=2.0)
ap.add_argument('--lacunarity', type=float, default=2.0)
ap.add_argument('--gain',       type=float, default=0.5)
ap.add_argument('--warp',       type=float, default=0.4)
ap.add_argument('--width',      type=int,   default=32)
ap.add_argument('--height',     type=int,   default=32)
ap.add_argument('--start-tile', type=str,   default='0,1', help='start tile as tx,ty')
ap.add_argument('--end-tile',   type=str,   default='1,1', help='end tile as tx,ty')
ap.add_argument('--verbosity',  type=int,   default=2, help='0: quiet, 1: endpoints, 2: full paths')
ap.add_argument('--ref',        type=str,   default='dijkstra', choices=['dijkstra','relax'], help='CPU reference: dijkstra or relax')
args = ap.parse_args()

# 1. TERRAIN GENERATION (32x32)
W, H = args.width, args.height
terrain = fbm_value_noise(W, H, octaves=args.octaves, base_freq=args.base_freq, lacunarity=args.lacunarity, gain=args.gain, warp=args.warp, seed=args.seed)

# 3. GPU TILE SOLVER
solver = GPUTerrainSolver(W, H)
parent_map, cost_map, sinks_packed = solver.solve(terrain, K=np.float32(0.005))
tile_size = 16
sink_xy = compute_tile_sinks(sinks_packed, W, H, tile_size=tile_size, unpack_fn=unpack)

def parse_tile(arg, nx, ny, name):
    tx, ty = [int(s) for s in arg.split(',')]
    if not (0 <= tx < nx and 0 <= ty < ny):
        raise ValueError(f"{name} out of range: ({tx},{ty}) not in [0,{nx-1}]x[0,{ny-1}]")
    return tx, ty

nx = W // tile_size
ny = H // tile_size
start_tile = parse_tile(args.start_tile, nx, ny, "start-tile")
end_tile = parse_tile(args.end_tile, nx, ny, "end-tile")
start_sink = tuple(sink_xy[start_tile[1], start_tile[0]])
end_sink = tuple(sink_xy[end_tile[1], end_tile[0]])
if args.ref == 'relax':
    if solve_relaxation_numpy is None:
        raise ValueError("solve_relaxation_numpy not available in terrain_ocl_2/terrain.py")
    ref_path = solve_relaxation_numpy(terrain, start_sink, end_sink, epsilon=0.001, max_iter=500)
    ref_cost = max(terrain[y, x] for x, y in ref_path) if len(ref_path) > 0 else np.inf
else:
    ref_path, ref_cost = global_dijkstra(terrain, start_sink, end_sink)

# 4. HIERARCHICAL GATE FINDING
gates = compute_gates(cost_map, terrain, tile_size=16)
# Normalize gates to include tile ids if missing
for g in gates:
    if 'a' not in g or 'b' not in g:
        g['a'] = (g['pA'][0] // tile_size, g['pA'][1] // tile_size)
        g['b'] = (g['pB'][0] // tile_size, g['pB'][1] // tile_size)
tile_path, gate_seq, tile_path_cost = tile_graph_shortest(gates, start_tile, end_tile)
best_gate = gate_seq[0]
last_gate = gate_seq[-1]
# Orient gates to start/end tiles
if best_gate['b'] == start_tile:
    best_gate = {'a': best_gate['b'], 'b': best_gate['a'], 'pA': best_gate['pB'], 'pB': best_gate['pA'], 'barrier': best_gate['barrier']}
if last_gate['b'] != end_tile:
    last_gate = {'a': last_gate['b'], 'b': last_gate['a'], 'pA': last_gate['pB'], 'pB': last_gate['pA'], 'barrier': last_gate['barrier']}
best_y_global = best_gate['pA'][1]
print(f"Tile path: {tile_path} (cost={tile_path_cost:.6f})")
print(f"Hierarchical GPU First Gate at Y: {best_y_global}")
log_ref("CPU reference", ref_path, ref_cost, args.verbosity)

# 5. PATH RECONSTRUCTION
def log_path(name, path, verb=0):
    pts = [(int(p[0]), int(p[1])) for p in path]
    print(f"{name}: start={pts[0]}, end={pts[-1]}, len={len(pts)}")
    if verb >= 2:
        if verb >= 3 and len(pts) > 50:
            sample = pts[:10] + ['...'] + pts[-10:]
            print(f"{name} full (truncated): {sample}")
        else:
            print(f"{name} full: {pts}")

def log_ref(name, path, cost, verb=0):
    pts = [(int(p[0]), int(p[1])) for p in path]
    print(f"{name}: cost={cost:.6f}, start={pts[0]}, end={pts[-1]}, len={len(pts)}")
    if verb >= 2:
        if verb >= 3 and len(pts) > 50:
            sample = pts[:10] + ['...'] + pts[-10:]
            print(f"{name} full (truncated): {sample}")
        else:
            print(f"{name} full: {pts}")

p_left = get_pixel_path(best_gate['pA'], parent_map, unpack_fn=unpack)
p_right = get_pixel_path(last_gate['pB'], parent_map, unpack_fn=unpack)
log_path("p_left", p_left, args.verbosity)
log_path("p_right", p_right, args.verbosity)

# 6. VISUALIZATION
fig, axs = plt.subplots(1, 2, figsize=(12, 6))

img0 = axs[0].imshow(terrain, cmap='terrain')
axs[0].contour(terrain, levels=10, colors='k', linewidths=0.5, alpha=0.5)
axs[0].axvline(15.5, color='red', linestyle='--', label='Tile Boundary')
axs[0].axhline(15.5, color='red', linestyle='--')
#axs[0].axhline(best_y_global, color='k', linestyle=':', alpha=0.6, label=f'Pass y={best_y_global}')
axs[0].plot(ref_path[:, 0], ref_path[:, 1], 'w-', linewidth=3, alpha=0.9, label='CPU ref path')
axs[0].plot(p_left[:, 0], p_left[:, 1], 'r-', linewidth=2, label='Start tile path (Sink->Gate)')
axs[0].plot(p_right[:, 0], p_right[:, 1], 'b-', linewidth=2, label='End tile path (Gate->Sink)')
axs[0].plot([best_gate['pA'][0], best_gate['pB'][0]], [best_gate['pA'][1], best_gate['pB'][1]], 'k-', linewidth=3, label='First Gate')
if last_gate != best_gate:
    axs[0].plot([last_gate['pA'][0], last_gate['pB'][0]], [last_gate['pA'][1], last_gate['pB'][1]], 'k--', linewidth=2, label='Last Gate')
sinks = [unpack(s) for s in sinks_packed]
sx, sy = zip(*sinks)
axs[0].scatter(sx, sy, c='yellow', s=100, edgecolors='k', zorder=5, label='Tile Sinks')
for idx, (x, y) in enumerate(sinks):
    axs[0].text(x, y, str(idx), color='k', fontsize=8, ha='center', va='center', zorder=6,
                bbox=dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.7))
axs[0].set_title(f"Terrain (seed={args.seed})")
fig.colorbar(img0, ax=axs[0], fraction=0.046, pad=0.04, label='height')
axs[0].legend()

img1 = axs[1].imshow(cost_map, cmap='inferno')
axs[1].axvline(15.5, color='red', linestyle='--')
axs[1].axhline(15.5, color='red', linestyle='--')
axs[1].axhline(best_y_global, color='k', linestyle=':', alpha=0.6)
axs[1].set_title("GPU intra-tile cost map")
fig.colorbar(img1, ax=axs[1], fraction=0.046, pad=0.04, label='cost')

plt.tight_layout()
plt.show()