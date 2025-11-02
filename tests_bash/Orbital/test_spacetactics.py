#!/usr/bin/env python3
import argparse
import os
import sys
from pathlib import Path

import numpy as np

# Ensure python/ is on path for module imports
THIS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = THIS_DIR.parent.parent
sys.path.append(str(PROJECT_ROOT / "python"))

from OrbitalWar import spacetactics as st  # noqa: E402

# Match thrust timeline used in C++ SpaceTactics app
THRUST_TIMES = np.array([-10.0, 35000.0, 40000.0, 60000.0, 65000.0, 100000.0], dtype=np.double)
THRUST_VALUES = np.array(
    [
        [0.000, 0.0, 0.0],
        [0.001, 0.0, 0.0],
        [-0.200, 0.0, 0.0],
        [-0.100, 0.0, 0.0],
        [-0.001, 0.0, 0.0],
        [0.000, 0.0, 0.0],
    ],
    dtype=np.double,
)

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="SpaceTactics ctypes smoke test")
    ap.add_argument("--n-steps", type=int, default=500)
    ap.add_argument("--dt", type=float, default=200.0)
    ap.add_argument("--save-npz", type=str, default=None, help="Optional output file for trajectories")
    args = ap.parse_args()

    # Enable verbose logging from the native layer so we can inspect buffer registration.
    st.set_debug(verbosity=0, idebug=0)

    info = st.setup_jupiter_moons_skirmish()

    n_steps = args.n_steps
    dt = args.dt

    print("py.DEBUG 0")

    st.allocate_trjs(n_steps)
    # Ensure the integrator step matches the thrust spline interpolation
    st.predict_trjs(0, dt)
    st.set_ship_thrust(3, THRUST_TIMES, THRUST_VALUES)
    st.predict_trjs(n_steps, dt)

    print("py.DEBUG 1")

    st.init_buffers()
    traj_len, traj_dt = st.get_trj_metadata()
    print(f"Trajectory buffers: len={traj_len}, dt={traj_dt}")

    ships = info["ships"]
    planets = info["planets"]

    print("py.DEBUG 2")

    ship_trajs = {}
    for idx, name in enumerate(ships):
        buf = st.get_buffer(f"ship[{idx}]:{name}.trjPos")
        if buf is None:
            raise RuntimeError(f"Missing buffer for {name}")
        ship_trajs[name] = buf
        print(f"Ship {name} trajectory shape: {buf.shape}, first point: {buf[0]}")

    print("py.DEBUG 3")    

    jupiter_traj = st.get_buffer("planet[1]:Jupiter.trjPos")
    if jupiter_traj is None:
        raise RuntimeError("Jupiter trajectory buffer missing")
    print(f"Jupiter trajectory sample (first 3 points):\n{jupiter_traj[:3]}")

    print("py.DEBUG 4")

    # Sanity checks duplicating LandCraft pattern: verify zero-copy by slicing
    for name, buf in ship_trajs.items():
        assert buf.shape[1] == 3, f"Unexpected stride for {name}"

    print("py.DEBUG 5")

    if args.save_npz:
        np.savez(args.save_npz, **{name: np.array(buf) for name, buf in ship_trajs.items()}, Jupiter=jupiter_traj)
        print(f"Saved trajectories to {args.save_npz}")

    print("ALL DONE")