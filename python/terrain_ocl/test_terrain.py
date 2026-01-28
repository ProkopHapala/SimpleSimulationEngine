import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
os.environ.setdefault("PYOPENCL_CTX", "1")
from terrain import GPUTerrainSolver, unpack
import importlib.util

_utils_path = os.path.join(os.path.dirname(__file__), "../terrain_ocl_2/terrain.py")
_spec = importlib.util.spec_from_file_location("terrain_utils", _utils_path)
terrain_utils = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(terrain_utils)
fbm_value_noise = terrain_utils.fbm_value_noise
compute_gates = terrain_utils.compute_gates
compute_tile_sinks = terrain_utils.compute_tile_sinks
get_pixel_path = terrain_utils.get_pixel_path


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