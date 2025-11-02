#!/usr/bin/env python3
import argparse
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
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


def select_names(all_names: list[str], include: list[str] | None, exclude: list[str] | None) -> list[str]:
    chosen = list(all_names) if include is None else [name for name in include if name in all_names]
    if exclude:
        exclusions = set(exclude)
        chosen = [name for name in chosen if name not in exclusions]
    return chosen


def fetch_traj_map(prefix: str, names: list[str], *, label: str) -> dict[str, np.ndarray]:
    trajs: dict[str, np.ndarray] = {}
    for idx, name in enumerate(names):
        buf = st.get_buffer(f"{prefix}[{idx}]:{name}.trjPos")
        if buf is None:
            raise RuntimeError(f"Missing buffer for {label.lower()} {name}")
        trajs[name] = buf
        print(f"{label} {name} trajectory shape: {buf.shape}, first point: {buf[0]}")
    return trajs


def compute_origin(frame: str, planet_trajs: dict[str, np.ndarray]) -> np.ndarray | None:
    frame = frame.lower()
    if frame == "jupiter":
        origin = planet_trajs.get("Jupiter")
    elif frame == "sun":
        origin = planet_trajs.get("Sun")
    else:
        origin = None
    if origin is None and frame != "absolute":
        raise RuntimeError(f"Origin trajectory for frame '{frame}' not available")
    return None if origin is None else np.asarray(origin)


def relative_traj(traj: np.ndarray, origin: np.ndarray | None) -> np.ndarray:
    arr = np.asarray(traj)
    if origin is None:
        return arr
    if origin.shape != arr.shape:
        raise ValueError("Trajectory length mismatch with origin frame")
    return arr - origin


def print_bounds(label: str, traj: np.ndarray) -> None:
    xmin, xmax = traj[:, 0].min(), traj[:, 0].max()
    ymin, ymax = traj[:, 1].min(), traj[:, 1].max()
    print(f"{label:>12} | X [{xmin: .3e}, {xmax: .3e}]  Y [{ymin: .3e}, {ymax: .3e}]")


def report_bounds(
    planet_trajs: dict[str, np.ndarray],
    planet_names: list[str],
    ship_trajs: dict[str, np.ndarray],
    ship_names: list[str],
    origin: np.ndarray | None,
    frame_label: str,
) -> None:
    print("\nTrajectory bounds (absolute frame):")
    abs_points: list[np.ndarray] = []
    for name in planet_names:
        arr = np.asarray(planet_trajs[name])
        print_bounds(f"Planet {name}", arr)
        abs_points.append(arr)
    for name in ship_names:
        arr = np.asarray(ship_trajs[name])
        print_bounds(f"Ship {name}", arr)
        abs_points.append(arr)
    if abs_points:
        print_bounds("Overall", np.vstack(abs_points))

    print("\nTrajectory bounds (frame =", frame_label, "):")
    rel_points: list[np.ndarray] = []
    if origin is not None:
        print_bounds("Origin", relative_traj(origin, origin))
    for name in planet_names:
        rel = relative_traj(planet_trajs[name], origin)
        print_bounds(f"Planet {name}", rel)
        rel_points.append(rel)
    for name in ship_names:
        rel = relative_traj(ship_trajs[name], origin)
        print_bounds(f"Ship {name}", rel)
        rel_points.append(rel)
    if rel_points:
        print_bounds("Overall", np.vstack(rel_points))


def plot_trajectories(
    planet_trajs: dict[str, np.ndarray],
    planet_names: list[str],
    ship_trajs: dict[str, np.ndarray],
    ship_names: list[str],
    origin: np.ndarray | None,
    frame_label: str,
) -> None:
    fig, ax = plt.subplots(figsize=(8, 8))
    for name in planet_names:
        rel = relative_traj(planet_trajs[name], origin)
        ax.plot(rel[:, 0], rel[:, 1], label=f"Planet {name}")
    for name in ship_names:
        rel = relative_traj(ship_trajs[name], origin)
        ax.plot(rel[:, 0], rel[:, 1], linestyle="--", label=f"Ship {name}")
    ax.set_title(f"SpaceTactics trajectories (XY plane, {frame_label} frame)")
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(loc="best")
    ax.grid(True)
    plt.show()

def _normalize_many(value: str | list[str] | None, *, all_token: str | None = None, default: list[str] | None = None) -> list[str] | None:
        if value is None:
            return default
        if isinstance(value, list):
            return value if value else default
        token = value.strip()
        if all_token is not None and token.lower() == all_token.lower():
            return None
        if not token:
            return [] if default is None else default
        return [token]

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="SpaceTactics ctypes smoke test")
    ap.add_argument("--n-steps",  type=int, default=500)
    ap.add_argument("--dt",       type=float, default=200.0)
    ap.add_argument("--save-npz", type=str, default=None, help="Optional output file for trajectories")
    ap.add_argument("--frame",      default="absolute", choices=["absolute", "jupiter", "sun"], help="Frame of reference for plotting trajectories")
    ap.add_argument("--planets",      nargs="+", default="all", help="Planet names to plot (default: all)")
    ap.add_argument("--omit-planets", nargs="+", default="Sun", help="Planet names to skip when plotting (default: Sun)")
    ap.add_argument("--ships",        nargs="+", default="all", help="Ship names to plot (default: all)")
    ap.add_argument("--omit-ships",   nargs="+", default="", help="Ship names to skip when plotting")
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

    ship_trajs = fetch_traj_map("ship", ships, label="Ship")
    planet_trajs = fetch_traj_map("planet", planets, label="Planet")

    print("py.DEBUG 3")    

    jupiter_traj = planet_trajs.get("Jupiter")
    if jupiter_traj is None:
        raise RuntimeError("Jupiter trajectory buffer missing")
    print(f"Jupiter trajectory sample (first 3 points):\n{jupiter_traj[:3]}")

    print("py.DEBUG 4")

    # Sanity checks duplicating LandCraft pattern: verify zero-copy by slicing
    for name, buf in ship_trajs.items():
        assert buf.shape[1] == 3, f"Unexpected stride for {name}"

    include_planets = _normalize_many(args.planets, all_token="all")
    omit_planets = _normalize_many(args.omit_planets, default=["Sun"])
    include_ships = _normalize_many(args.ships, all_token="all")
    omit_ships = _normalize_many(args.omit_ships, default=[])

    planet_names = select_names(planets, include_planets, omit_planets)
    ship_names = select_names(ships, include_ships, omit_ships)

    # Determine frame of reference
    origin_traj = compute_origin(args.frame, planet_trajs)
    frame_label = args.frame.capitalize() if args.frame.lower() != "absolute" else "Absolute"

    # Plot and report diagnostics
    plot_trajectories(planet_trajs, planet_names, ship_trajs, ship_names, origin_traj, frame_label)
    report_bounds(planet_trajs, planet_names, ship_trajs, ship_names, origin_traj, frame_label)

    if args.save_npz:
        np.savez(args.save_npz, **{name: np.array(buf) for name, buf in ship_trajs.items()}, Jupiter=planet_trajs["Jupiter"])
        print(f"Saved trajectories to {args.save_npz}")

    print("ALL DONE")