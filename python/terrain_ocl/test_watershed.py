import numpy as np
import pyopencl as cl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import argparse, os, sys
os.environ.setdefault("PYOPENCL_CTX", "1")

from terrain import (
    fbm_value_noise_2,
    build_watershed_program,
    run_calc_flow,
    run_propagate_basins,
)

def main():
    ap = argparse.ArgumentParser(description="Verify global watershed (steepest descent + pointer jumping)")
    ap.add_argument('--seed',       type=int,   default=1)
    ap.add_argument('--octaves',    type=int,   default=6)
    ap.add_argument('--base-freq',  type=float, default=2.0)
    ap.add_argument('--lacunarity', type=float, default=2.0)
    ap.add_argument('--gain',       type=float, default=0.5)
    ap.add_argument('--warp',       type=float, default=0.4)
    ap.add_argument('--width',      type=int,   default=64)
    ap.add_argument('--height',     type=int,   default=64)
    ap.add_argument('--max-conns',  type=int,   default=50, help='plot top-N lowest gates')
    args = ap.parse_args()

    W, H = args.width, args.height
    terrain = fbm_value_noise_2(W, H, octaves=args.octaves, base_freq=args.base_freq,
                                lacunarity=args.lacunarity, gain=args.gain,
                                warp=args.warp, seed=args.seed)

    ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)
    prg = build_watershed_program(ctx)

    # Buffers
    mf = cl.mem_flags
    h_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=terrain)
    flow_buf = cl.Buffer(ctx, mf.READ_WRITE, 4 * W * H)
    basin_buf = cl.Buffer(ctx, mf.READ_WRITE, 4 * W * H)
    changed_buf = cl.Buffer(ctx, mf.READ_WRITE, 4)

    # 1. Calculate Flow
    run_calc_flow(queue, prg, h_buf, flow_buf, basin_buf, W, H)

    # 2. Propagate Basins (Pointer Jumping)
    for i in range(W*2): # Safety limit
        changed = np.zeros(1, dtype=np.int32)
        cl.enqueue_copy(queue, changed_buf, changed)
        run_propagate_basins(queue, prg, flow_buf, basin_buf, changed_buf, W, H)
        cl.enqueue_copy(queue, changed, changed_buf)
        if changed[0] == 0:
            print(f"Basins converged in {i} iterations.")
            break

    # 3. Retrieve Data
    flow_map = np.empty((H, W), dtype=np.uint32)
    basin_map = np.empty((H, W), dtype=np.int32)
    cl.enqueue_copy(queue, flow_map, flow_buf)
    cl.enqueue_copy(queue, basin_map, basin_buf)

    # 4. Find Ridge Gates (CPU Side for now - easy to vectorize)
    # Horizontal diff
    b_curr = basin_map[:, :-1]
    b_right = basin_map[:, 1:]
    diff_mask = (b_curr != b_right)
    
    # Candidate Gates (Coords)
    y_coords, x_coords = np.where(diff_mask)
    # The 'Barrier Height' is max(h[y,x], h[y,x+1])
    # Actually, physically, the gate is at the 'Pass', which is the Saddle.
    # The water flows A -> Gate -> B? No.
    # Water flows Gate -> A and Gate -> B.
    # So the gate height IS the barrier.
    
    # Store edges: (Basin1, Basin2) -> (MinHeight, GateX, GateY)
    adj = {} 
    
    for y, x in zip(y_coords, x_coords):
        id1 = b_curr[y, x]
        id2 = b_right[y, x]
        if id1 > id2: id1, id2 = id2, id1 # Canonize
        
        h_val = max(terrain[y,x], terrain[y,x+1])
        
        key = (id1, id2)
        if key not in adj or h_val < adj[key][0]:
            adj[key] = (h_val, x, y, 'H') # 'H' for horizontal boundary

    # Do Vertical diff similarly... (omitted for brevity)

    print(f"Found {len(adj)} connections between basins.")

    # 5. Visualization
    plt.figure(figsize=(12, 6))
    plt.subplot(1, 2, 1)
    unique_ids = np.unique(basin_map)
    colors = np.zeros((unique_ids.size, 3))
    for i, bid in enumerate(unique_ids):
        # hash basin id and seed into a uint64 to avoid overflow
        h = (np.uint64(bid) ^ (np.uint64(args.seed) * np.uint64(0x9e3779b97f4a7c15))) & np.uint64((1 << 64) - 1)
        rng_i = np.random.default_rng(int(h))
        colors[i] = rng_i.random(3)
    cmap = mcolors.ListedColormap(colors)
    # Map basin IDs to [0, n) indices
    id_to_idx = {bid: i for i, bid in enumerate(unique_ids)}
    mapped = np.vectorize(id_to_idx.get)(basin_map)
    plt.imshow(mapped, cmap=cmap, interpolation='nearest')
    plt.title(f"Natural Watershed Basins (N={len(unique_ids)})")
    
    # Plot Sinks
    unique_sinks = np.unique(basin_map)
    for sid in unique_sinks:
        sy, sx = divmod(sid, W)
        plt.plot(sx, sy, 'ko', markersize=2)

    # Plot Gates and Paths
    plt.subplot(1, 2, 2)
    plt.imshow(terrain, cmap='gray')
    
    def trace_down(sx, sy):
        path = []
        curr = (sx, sy)
        for _ in range(1000):
            path.append(curr)
            packed = flow_map[curr[1], curr[0]]
            nx, ny = packed & 0xFFFF, (packed >> 16) & 0xFFFF
            if (nx, ny) == curr: break
            curr = (nx, ny)
        return np.array(path)

    # Draw Top 50 Connections (Lowest Gates)
    sorted_edges = sorted(adj.values(), key=lambda x: x[0])
    
    for _, gx, gy, type in sorted_edges[:50]:
        plt.plot(gx, gy, 'r*', markersize=5)
        
        # Trace down to Basin A
        p1 = trace_down(gx, gy)
        plt.plot(p1[:,0], p1[:,1], 'c-', linewidth=0.5, alpha=0.7)
        
        # Trace down to Basin B (Neighbor pixel)
        nx, ny = (gx+1, gy) if type == 'H' else (gx, gy+1)
        p2 = trace_down(nx, ny)
        plt.plot(p2[:,0], p2[:,1], 'c-', linewidth=0.5, alpha=0.7)

    plt.title("Ridge Gates & Downstream Flows")
    plt.show()

if __name__ == "__main__":
    main()