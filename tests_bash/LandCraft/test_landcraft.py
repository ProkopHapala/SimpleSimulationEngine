#!/usr/bin/env python3
import sys
import os
import numpy as np
import argparse
import time
sys.path.append('../../python')
print("DEBUG == 3")
from LandCraft import landcraft as lc
print("DEBUG == 4")
# ---- Helpers ----

def get_buffers(nx, ny):
    """Wrap C buffers as NumPy arrays (views). Shapes are (ny,nx)."""
    lc.init_buffers()
    ground = lc.getBuff('ground', nx*ny).reshape((ny, nx))
    water  = lc.getBuff('water',  nx*ny).reshape((ny, nx))
    return ground, water

def save_pgm(path, A):
    """Save 2D array as 8-bit grayscale PGM (no matplotlib)."""
    H, W = A.shape
    amin = float(np.nanmin(A))
    amax = float(np.nanmax(A))
    rng  = (amax - amin) if (amax > amin) else 1.0
    B = ((A - amin) / rng * 255.0).clip(0,255).astype(np.uint8)
    with open(path, 'wb') as f:
        f.write(f"P5 {W} {H} 255\n".encode('ascii'))
        f.write(B.tobytes(order='C'))


def python_set_terrain(ground, kind='gauss', seed=123, max_height=500.0):
    ny, nx = ground.shape
    if kind == 'flat':
        ground[:] = 0.0
    elif kind == 'rand':
        rng = np.random.default_rng(seed)
        ground[:] = rng.standard_normal(size=ground.shape)
        ground *= (max_height / 8.0)
    elif kind == 'gauss':
        xs = (np.arange(nx) - 0.5*nx) / nx
        ys = (np.arange(ny) - 0.5*ny) / ny
        X, Y = np.meshgrid(xs, ys)
        R2 = X*X + Y*Y
        ground[:] = max_height * np.exp(-8.0 * R2)
    elif kind == 'sine':
        xs = np.arange(nx) / nx
        ys = np.arange(ny) / ny
        X, Y = np.meshgrid(xs, ys)
        ground[:] = max_height * (0.3*np.sin(6.28*X) + 0.2*np.sin(4.71*Y))
    else:
        raise ValueError(f"Unknown terrain kind '{kind}'")


def basin_fill(outflow_xy=(0,0), n_iter=50, verbose=1):
    """Simple basin filling using water relaxation with a designated outflow cell.
    - Sets outflow at (ix,iy), then applies multiple lc.relax_all() iterations.
    """
    ix, iy = outflow_xy
    lc.set_outflow_at(ix, iy)
    for it in range(n_iter):
        lc.relax_all()
        if verbose and (it % 10 == 0 or it == n_iter-1):
            print(f"relax_all iter {it+1}/{n_iter}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='LandCraft ctypes smoke test')
    ap.add_argument('--nx',              type=int, default=256)
    ap.add_argument('--ny',              type=int, default=256)
    ap.add_argument('--data-folder',     type=str, default=None, help='Optional base data folder (e.g., tests_bash/apps/data)')
    ap.add_argument('--terrain',         type=str, default='cpp_bisec', choices=['none','cpp','cpp_bisec','gauss','rand','sine','flat'], help='Terrain source: none, cpp (generate_terrain), cpp_bisec (bisec+optional erosion), or python models')
    ap.add_argument('--seed',            type=int, default=16464)
    ap.add_argument('--max-height',      type=float, default=500.0)
    # Neighborhood
    ap.add_argument('--neighbors',       type=int, default=6, choices=[4,6,8], help='Neighborhood topology for hydraulics (4/6/8)')
    # Erosion (used when terrain=cpp_bisec)
    ap.add_argument('--erosion-iters',   type=int, default=400, help='Number of erosion iterations (cpp_bisec only; 0=skip)')
    ap.add_argument('--erosion-drops',   type=int, default=400)
    ap.add_argument('--erosion-steps',   type=int, default=500)
    ap.add_argument('--erosion-margin',  type=int, default=25)
    ap.add_argument('--erosion-min',     type=float, default=0.10)
    ap.add_argument('--erosion-max',     type=float, default=0.15)
    ap.add_argument('--erosion-prob',    type=float, default=0.90)
    # Hydro / rivers
    ap.add_argument('--outflow-x',       type=int, default=0)
    ap.add_argument('--outflow-y',       type=int, default=0)
    ap.add_argument('--relax-iters',     type=int, default=50)
    ap.add_argument('--force-relax-on-load', action='store_true', help='If set, run basin fill even when terrain was loaded from cache')
    ap.add_argument('--rivers-min-flow', type=float, default=50.0)
    ap.add_argument('--plot',            type=int, default=2, help='Show/save matplotlib plots of ground/water')
    # Persistence
    ap.add_argument('--load-prefix',     type=str, default="./", help='If set, load ground.bin and water.bin from this directory/prefix before generation')
    ap.add_argument('--save-prefix',     type=str, default="./", help='If set, save ground.bin and water.bin to this directory/prefix at the end')
    args = ap.parse_args()


    lc.world_init(args.data_folder)
    lc.map_init(args.nx, args.ny)
    # Set neighborhood topology early
    lc.set_neighbors(args.neighbors)
 
    # 2) Init buffers and wrap as NumPy arrays
    ground, water = get_buffers(args.nx, args.ny)
    print(f"Buffers ready: ground.shape={ground.shape}, water.shape={water.shape}")

    # 3) Load/generate terrain with cache guard
    loaded = False
    will_save = False
    load_prefix = args.load_prefix
    save_prefix = args.save_prefix if args.save_prefix is not None else args.load_prefix
    if load_prefix is not None:
        gpath = os.path.join(load_prefix, 'ground.bin')
        wpath = os.path.join(load_prefix, 'water.bin')
        if os.path.isfile(gpath) and os.path.isfile(wpath):
            rc = lc.load(gpath, wpath)
            print(f"Loaded terrain: rc={rc} from {gpath}, {wpath}")
            loaded = (rc == 0)
            if not loaded:
                print(f"WARNING: cache present but load failed (rc={rc}); will regenerate.")
                will_save = True
        else:
            print(f"Cache not found in {load_prefix}; will generate and save.")
            will_save = True

    if (not loaded) and args.terrain == 'cpp':
        lc.generate_terrain(args.seed, args.max_height)
        print('Terrain generated by C++ (noise+erosion).')
    elif (not loaded) and args.terrain == 'cpp_bisec':
        # Random bisect noise via C++; optional erosion; then scale/initialize water in Python
        lc.make_terrain_bisec(args.seed)
        if args.erosion_iters > 0:
            lc.droplet_erosion(args.erosion_iters, args.erosion_drops, args.erosion_steps, args.erosion_margin, args.erosion_min, args.erosion_max, args.erosion_prob)
        # Normalize/scale like generate_terrain: water = ground*H; ground *= H
        ground *= args.max_height
        water[:] = ground
        print('Terrain generated by bisect noise (cpp_bisec) with optional erosion and scaled in Python.')
    elif (not loaded) and args.terrain in ('gauss','rand','sine','flat'):
        python_set_terrain(ground, kind=args.terrain, seed=args.seed, max_height=args.max_height)
        # Set initial water level equal to ground for stability, similar to C++ generate
        water[:] = ground
        print(f"Terrain set in Python: kind={args.terrain}")
    elif not loaded:
        print('Terrain left unchanged (none).')
  
    # 5) Basin filling (set outflow + relax)
    if loaded and not args.force_relax_on_load:
        print("Loaded terrain; skipping basin fill (use --force-relax-on-load to run it).")
    else:
        basin_fill(outflow_xy=(args.outflow_x, args.outflow_y), n_iter=args.relax_iters, verbose=1)

    # 6) Save if requested or if we promised to save cache
    if (save_prefix is not None) and (will_save or args.save_prefix is not None):
        os.makedirs(save_prefix, exist_ok=True)
        gpath = os.path.join(save_prefix, 'ground.bin')
        wpath = os.path.join(save_prefix, 'water.bin')
        rc = lc.save(gpath, wpath)
        print(f"Saved terrain: rc={rc} to {gpath}, {wpath}")

    # Optional: compute rivers after rain gathering
    wmax = lc.gather_rain(100.0)
    n_riv = lc.find_all_rivers(args.rivers_min_flow)
    print(f"Rivers found: {n_riv} (wmax={wmax:.3f})")


    #min,max of ground and water
    print("ground min,max=", np.min(ground), np.max(ground))
    print("water  min,max=", np.min(water), np.max(water))

    #print('ground[0,:10]=', np.array2string(ground[:,:], precision=3))
    #print('water [0,:10]=', np.array2string(water [:,:], precision=3))

    # # Write images without matplotlib (ASan-safe)
    # save_pgm('ground.pgm', ground)
    # save_pgm('water.pgm',  water)
    
    # Optional plotting
    # mode 1 (default): save images via matplotlib.image (no pyplot, avoids ft2font)
    # mode 2: interactive pyplot (titles/colorbars). Use without ASan preload.
    if args.plot == 1:
        import matplotlib
        matplotlib.use('Agg')  # non-interactive backend
        from matplotlib import image as mpimg, cm
        mpimg.imsave('ground.png', ground, cmap=cm.gist_earth, origin='lower')
        mpimg.imsave('water.png',  water,  cmap=cm.magma, origin='lower')
    elif args.plot == 2:
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(1, 2, figsize=(10, 4))
        im0 = axs[0].imshow(ground, origin='lower', cmap='gist_earth', interpolation='nearest')
        axs[0].set_title('Ground height')
        fig.colorbar(im0, ax=axs[0], fraction=0.046, pad=0.04)
        im1 = axs[1].imshow(water, origin='lower', cmap='viridis', interpolation='nearest')
        axs[1].set_title('Water level')
        fig.colorbar(im1, ax=axs[1], fraction=0.046, pad=0.04)
        plt.tight_layout()
        plt.savefig('landcraft_buffers.png', dpi=150)
        plt.show()

    # Show a few samples for sanity
    # print('ground[0,:10]=', np.array2string(ground[:,:], precision=3))
    # print('water [0,:10]=', np.array2string(water [:,:], precision=3))
    # print("DEBUG 9")


