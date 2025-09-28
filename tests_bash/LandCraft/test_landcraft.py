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

def plot_hydro(ground, water, moveCost=None, seeds=None, mark_xy=None):
    import matplotlib.pyplot as plt
    depth = np.maximum(water - ground, 0.0)
    ncols = 4 if moveCost is not None else 3
    plt.figure(figsize=(4*ncols,4))
    plt.subplot(1,ncols,1); plt.title('ground'); plt.imshow(ground, cmap='terrain');
    if moveCost is not None:
        plt.subplot(1,ncols,2); plt.title('moveCost'); plt.imshow(moveCost, cmap='viridis');
        k=3
    else:
        k=2
    plt.subplot(1,ncols,k); plt.title('water'); plt.imshow(water, cmap='Blues')
    plt.subplot(1,ncols,k+1); plt.title('depth (water-ground)'); plt.imshow(depth, cmap='magma')
    # overlays
    if seeds is not None:
        H, W = ground.shape
        ys, xs = np.divmod(seeds, W)
        if moveCost is not None:
            plt.subplot(1,ncols,2); plt.scatter(xs, ys, s=6, c='r', marker='o')
        plt.subplot(1,ncols,k);   plt.scatter(xs, ys, s=6, c='r', marker='o')
        plt.subplot(1,ncols,k+1); plt.scatter(xs, ys, s=6, c='r', marker='o')
    if mark_xy is not None:
        mx, my = mark_xy
        if moveCost is not None:
            plt.subplot(1,ncols,2); plt.scatter([mx],[my], s=40, c='y', marker='x')
        plt.subplot(1,ncols,k);   plt.scatter([mx],[my], s=40, c='y', marker='x')
        plt.subplot(1,ncols,k+1); plt.scatter([mx],[my], s=40, c='y', marker='x')
    plt.tight_layout(); 
    #plt.show()


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
    ap.add_argument('--nx',              type=int, default=512)
    ap.add_argument('--ny',              type=int, default=512)
    ap.add_argument('--data-folder',     type=str, default=None, help='Optional base data folder (e.g., tests_bash/apps/data)')
    ap.add_argument('--terrain',         type=str, default='cpp_bisec', choices=['none','cpp','cpp_bisec','gauss','rand','sine','flat'], help='Terrain source: none, cpp (generate_terrain), cpp_bisec (bisec+optional erosion), or python models')
    ap.add_argument('--seed',            type=int, default=16464)
    ap.add_argument('--max-height',      type=float, default=500.0)
    # Neighborhood
    ap.add_argument('--neighbors',       type=int, default=6, choices=[4,6,8], help='Neighborhood topology for hydraulics (4/6/8)')
    # Erosion (used when terrain=cpp_bisec)
    ap.add_argument('--erosion-iters',   type=int, default=3200, help='Number of erosion iterations (cpp_bisec only; 0=skip)')
    ap.add_argument('--erosion-drops',   type=int, default=800)
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
    # Basin filling options
    ap.add_argument('--basin',            type=str, default='contour', choices=['none','bellman','dijkstra','contour'], help='Basin filling algorithm to run')
    ap.add_argument('--basin-iters',      type=int, default=100, help='Max iterations for Bellman-Ford')
    ap.add_argument('--basin-apply',      type=int, default=1, choices=[0,1], help='Apply spill levels to water after computation')
    ap.add_argument('--basin-level-cap',  type=float, default=0.0, help='Level cap for contour flood (<=0 disables)')
    ap.add_argument('--drain-x',          type=int, default=-1, help='Optional explicit drainage seed x (used as seed)')
    ap.add_argument('--drain-y',          type=int, default=-1, help='Optional explicit drainage seed y (used as seed)')
    ap.add_argument('--drain-level',      type=float, default=None, help='Optional initial water level at the drain seed for Bellman/Dijkstra')
    ap.add_argument('--force-rivers',     action='store_true', help='Force river computation even when basin mode is used')
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

    # =========== Droplet erosion ===========

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
        fig, axs = plt.subplots(1, 2, figsize=(20,10))
        im0 = axs[0].imshow(ground, origin='lower', cmap='gist_earth', interpolation='nearest')
        axs[0].set_title('Ground height')
        fig.colorbar(im0, ax=axs[0], fraction=0.046, pad=0.04)
        im1 = axs[1].imshow(water, origin='lower', cmap='viridis', interpolation='nearest')
        axs[1].set_title('Water level')
        fig.colorbar(im1, ax=axs[1], fraction=0.046, pad=0.04)
        plt.tight_layout()
        plt.savefig('landcraft_buffers.png', dpi=150)
        #plt.show()

    # =========== Basin filling ===========

    # 5) Basin filling tests (priority-flood/Bellman/contour) or relax fallback
    def boundary_seeds(nx, ny):
        ids = []
        # top/bottom rows
        for ix in range(nx): ids.append(ix); ids.append((ny-1)*nx + ix)
        # left/right cols
        for iy in range(ny): ids.append(iy*nx); ids.append(iy*nx + (nx-1))
        return np.array(ids, dtype=np.int32)

    if args.basin != 'none':
        apply = bool(args.basin_apply)
        seeds = boundary_seeds(args.nx, args.ny)
        # Optional explicit drain seed
        drain_seed = None
        if args.drain_x >= 0 and args.drain_y >= 0:
            if 0 <= args.drain_x < args.nx and 0 <= args.drain_y < args.ny:
                drain_seed = args.drain_y*args.nx + args.drain_x
                seeds = np.concatenate([seeds, np.array([drain_seed], dtype=np.int32)])
                # For Bellman/Dijkstra, initial level uses water[s] if present; set it if requested
                if args.drain_level is not None:
                    water[args.drain_y, args.drain_x] = max(ground[args.drain_y, args.drain_x], float(args.drain_level))
            else:
                print(f"WARNING: drain seed ({args.drain_x},{args.drain_y}) out of bounds; ignored.")
        t0 = time.time()
        if args.basin == 'bellman':
            niter = lc.basin_bellman_boundary(args.basin_iters, apply)
            print(f"Basin Bellman boundary: iters={niter} apply={apply}")
        elif args.basin == 'dijkstra':
            lc.basin_dijkstra_boundary(apply)
            print(f"Basin Dijkstra boundary: apply={apply}")
        elif args.basin == 'contour':
            lc.basin_contour_begin_flood(seeds, args.basin_level_cap)
            steps=0
            while True:
                added = lc.basin_contour_flood_step(args.basin_level_cap)
                steps+=1
                if added==0: break
            lc.basin_contour_begin_drain()
            dsteps=0
            while True:
                dadded = lc.basin_contour_drain_step()
                dsteps+=1
                if dadded==0: break
            lc.basin_contour_finish(apply)
            print(f"Basin Contour two-wave: flood_steps={steps} drain_steps={dsteps} apply={apply}")
        dt = time.time()-t0
        # Diagnostics on moveCost
        moveCost = lc.getBuff('moveCost', args.nx*args.ny).reshape((args.ny,args.nx))
        print(f"moveCost min,max= {float(np.nanmin(moveCost))} {float(np.nanmax(moveCost))}  (dt={dt:.3f}s)")
        # Optional quick plot of seeds/drain overlay including depth
        if args.plot == 2:
            mark_xy = (args.drain_x, args.drain_y) if drain_seed is not None else None
            plot_hydro(ground, water, moveCost=moveCost, seeds=seeds, mark_xy=mark_xy)
        # Skip relax if we already performed basin filling
    else:
        if loaded and not args.force_relax_on_load:
            print("Loaded terrain; skipping basin fill (use --force-relax-on-load to run it).")
        else:
            basin_fill(outflow_xy=(args.outflow_x, args.outflow_y), n_iter=args.relax_iters, verbose=1)
            if args.plot == 2:
                # Show relax-based outflow seed and depth as well
                mark_xy = (args.outflow_x, args.outflow_y)
                plot_hydro(ground, water, moveCost=None, seeds=None, mark_xy=mark_xy)

    # 6) Save if requested or if we promised to save cache
    if (save_prefix is not None) and (will_save or args.save_prefix is not None):
        os.makedirs(save_prefix, exist_ok=True)
        gpath = os.path.join(save_prefix, 'ground.bin')
        wpath = os.path.join(save_prefix, 'water.bin')
        rc = lc.save(gpath, wpath)
        print(f"Saved terrain: rc={rc} to {gpath}, {wpath}")

    # Optional: compute rivers after rain gathering (disabled by default to avoid heavy recursion)
    if args.force_rivers:
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
    


    # Show a few samples for sanity
    # print('ground[0,:10]=', np.array2string(ground[:,:], precision=3))
    # print('water [0,:10]=', np.array2string(water [:,:], precision=3))
    # print("DEBUG 9")

    plt.show()


