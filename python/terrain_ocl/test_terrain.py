import argparse
import numpy as np
import matplotlib.pyplot as plt
from terrain import GPUTerrainSolver, unpack


def smoothstep(t):
    return t*t*(3.0 - 2.0*t)


def hash_u32(x):
    x = (x ^ (x >> 16)) & 0xFFFFFFFF
    x = (x * 0x7feb352d) & 0xFFFFFFFF
    x = (x ^ (x >> 15)) & 0xFFFFFFFF
    x = (x * 0x846ca68b) & 0xFFFFFFFF
    x = (x ^ (x >> 16)) & 0xFFFFFFFF
    return x


def value_noise_2d(x, y, seed=1):
    xi = np.floor(x).astype(np.int64)
    yi = np.floor(y).astype(np.int64)
    xf = (x - xi).astype(np.float32)
    yf = (y - yi).astype(np.float32)
    u = smoothstep(xf)
    v = smoothstep(yf)
    def rnd(ix, iy):
        h = (np.uint32(ix) * np.uint32(0x9e3779b1) + np.uint32(iy) * np.uint32(0x85ebca6b) + np.uint32(seed) * np.uint32(0x27d4eb2d)) & np.uint32(0xFFFFFFFF)
        r = hash_u32(h)
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


def get_pixel_path(start_coord, parents, max_steps=2000):
    path = [start_coord]
    curr = start_coord
    for _ in range(max_steps):
        packed = parents[curr[1], curr[0]]
        nxt = unpack(packed)
        path.append(nxt)
        if nxt == curr:
            break
        curr = nxt
    return np.array(path)


def compute_tile_sinks(sinks, W, H, tile_size=16):
    nx = W // tile_size
    ny = H // tile_size
    sink_xy = np.zeros((ny, nx, 2), dtype=np.int32)
    for it, sp in enumerate(sinks):
        tx = it % nx
        ty = it // nx
        sink_xy[ty, tx, :] = unpack(sp)
    return sink_xy


def compute_gates(costs, tile_size=16):
    H, W = costs.shape
    nx = W // tile_size
    ny = H // tile_size
    gates = []
    for ty in range(ny):
        y0 = ty * tile_size
        ys = slice(y0, y0 + tile_size)
        for tx in range(nx - 1):
            xL = (tx + 1) * tile_size - 1
            xR = xL + 1
            cL = costs[ys, xL]
            cR = costs[ys, xR]
            b = np.maximum(cL, cR)
            i = int(b.argmin())
            gates.append({'a': (tx, ty), 'b': (tx + 1, ty), 'pA': (xL, y0 + i), 'pB': (xR, y0 + i), 'barrier': float(b[i])})
    for ty in range(ny - 1):
        yT = (ty + 1) * tile_size - 1
        yB = yT + 1
        for tx in range(nx):
            x0 = tx * tile_size
            xs = slice(x0, x0 + tile_size)
            cT = costs[yT, xs]
            cB = costs[yB, xs]
            b = np.maximum(cT, cB)
            i = int(b.argmin())
            gates.append({'a': (tx, ty), 'b': (tx, ty + 1), 'pA': (x0 + i, yT), 'pB': (x0 + i, yB), 'barrier': float(b[i])})
    return gates


def main():
    ap = argparse.ArgumentParser(description="GPU tile flooding demo with hierarchical noise")
    ap.add_argument('--width', type=int, default=256)
    ap.add_argument('--height', type=int, default=256)
    ap.add_argument('--octaves', type=int, default=7)
    ap.add_argument('--base-freq', type=float, default=2.0)
    ap.add_argument('--lacunarity', type=float, default=2.0)
    ap.add_argument('--gain', type=float, default=0.5)
    ap.add_argument('--warp', type=float, default=0.45)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--K', type=float, default=0.005, help='quadratic potential weight')
    ap.add_argument('--tile-size', type=int, default=16)
    ap.add_argument('--draw-flow', action='store_true')
    ap.add_argument('--flow-stride', type=int, default=8, help='pixel spacing for drawn flow paths (when --draw-flow)')
    ap.add_argument('--max-flow-paths', type=int, default=20000)
    ap.add_argument('--edge-alpha', type=float, default=0.6)
    ap.add_argument('--edge-lw', type=float, default=0.8)
    args = ap.parse_args()

    W, H = args.width, args.height
    if (W % args.tile_size) != 0 or (H % args.tile_size) != 0:
        raise SystemExit('width and height must be multiples of tile-size')
    h_map = fbm_value_noise(W, H, octaves=args.octaves, base_freq=args.base_freq, lacunarity=args.lacunarity, gain=args.gain, warp=args.warp, seed=args.seed)

    solver = GPUTerrainSolver(W, H)
    parents, costs, sinks = solver.solve(h_map, K=np.float32(args.K))
    sink_xy = compute_tile_sinks(sinks, W, H, tile_size=args.tile_size)
    gates = compute_gates(costs, tile_size=args.tile_size)

    plt.figure(figsize=(16, 7))

    plt.subplot(1, 2, 1)
    plt.imshow(costs, cmap='terrain')
    plt.title("Flooding Cost (Local Basins)")
    for x in range(args.tile_size, W, args.tile_size):
        plt.axvline(x - 0.5, color='k', alpha=0.15, linewidth=0.5)
    for y in range(args.tile_size, H, args.tile_size):
        plt.axhline(y - 0.5, color='k', alpha=0.15, linewidth=0.5)

    plt.subplot(1, 2, 2)
    plt.imshow(h_map, cmap='gray')
    for x in range(args.tile_size, W, args.tile_size):
        plt.axvline(x - 0.5, color='w', alpha=0.08, linewidth=0.6)
    for y in range(args.tile_size, H, args.tile_size):
        plt.axhline(y - 0.5, color='w', alpha=0.08, linewidth=0.6)

    if args.draw_flow:
        n_draw = 0
        for tx in range(0, W, args.flow_stride):
            for ty in range(0, H, args.flow_stride):
                p = get_pixel_path((tx, ty), parents)
                plt.plot(p[:, 0], p[:, 1], color='cyan', alpha=0.10, linewidth=0.5)
                n_draw += 1
                if n_draw >= args.max_flow_paths:
                    break
            if n_draw >= args.max_flow_paths:
                break

    for ty in range(sink_xy.shape[0]):
        for tx in range(sink_xy.shape[1]):
            sx, sy = sink_xy[ty, tx]
            plt.plot(sx, sy, 'ro', markersize=2.5)

    for g in gates:
        xA, yA = g['pA']
        xB, yB = g['pB']
        plt.plot(xA, yA, 'ys', markersize=3.0, markeredgecolor='k', markeredgewidth=0.4)
        plt.plot(xB, yB, 'ys', markersize=3.0, markeredgecolor='k', markeredgewidth=0.4)
        pA = get_pixel_path((xA, yA), parents)
        pB = get_pixel_path((xB, yB), parents)
        plt.plot(pA[:, 0], pA[:, 1], color='orange', alpha=args.edge_alpha, linewidth=args.edge_lw)
        plt.plot(pB[:, 0], pB[:, 1], color='orange', alpha=args.edge_alpha, linewidth=args.edge_lw)
        sAx, sAy = sink_xy[g['a'][1], g['a'][0]]
        sBx, sBy = sink_xy[g['b'][1], g['b'][0]]
        plt.plot([sAx, xA, xB, sBx], [sAy, yA, yB, sBy], color='lime', alpha=0.35, linewidth=0.8)

    plt.title("Tile sink connectivity via lowest-barrier gates (yellow); gate->sink paths (orange); schematic links (green)")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()