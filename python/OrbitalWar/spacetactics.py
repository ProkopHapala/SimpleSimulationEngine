import numpy as np
import ctypes
from ctypes import c_int, c_double, c_char_p, c_void_p
import os, sys

# Append the parent directory to sys.path to allow imports from pyMeta
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from pyMeta import cpp_utils

# --- Library Loading
# Adjust the build path and library name as needed
cpp_utils.BUILD_PATH = os.path.normpath(
    os.path.join(cpp_utils.PACKAGE_PATH, "../../cpp/Build/libs/OrbitalWar")
)
lib = cpp_utils.loadLib(
    "SpaceTacticsLib"
)  # Assumes the library is named libSpaceTacticsLib.so or similar

# --- ctypes Helpers
c_double_p = ctypes.POINTER(c_double)

array1d = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags="CONTIGUOUS")


def _wrap_ptr(name: str, length: int, stride: int = 1, shape=None) -> np.ndarray | None:
    """Turn a raw pointer from the C layer into a NumPy array view.

    Args:
        name: Buffer name registered by the C layer (e.g., 'ship[0]:A.trjPos').
        length: Expected number of doubles in the buffer.
        stride: Logical stride for reshaping. Defaults to 1 (flat array).
        shape: Optional target shape; if omitted, inferred from length/stride.
    """
    ptr = lib.sw_getBuff(name.encode("utf-8"))
    if not ptr:
        return None
    size = lib.sw_getBuffSize(name.encode("utf-8"))
    if size <= 0:
        return None
    if length > 0 and size < length:
        raise ValueError(f"Buffer {name} smaller than expected: {size} < {length}")
    if shape is None:
        if stride > 1:
            n = size // stride
            shape = (n, stride)
        else:
            shape = (size,)
    arr = np.ctypeslib.as_array(
        ptr, shape=(shape[0] * shape[1] if len(shape) == 2 else shape[0],)
    )
    return arr.reshape(shape)


def _wrap_ptr(name: str, length: int, stride: int = 1, shape=None) -> np.ndarray | None:
    """Turn a raw pointer from the C layer into a NumPy array view.

    Args:
        name: Buffer name registered by the C layer (e.g., 'ship[0]:A.trjPos').
        length: Expected number of doubles in the buffer.
        stride: Logical stride for reshaping. Defaults to 1 (flat array).
        shape: Optional target shape; if omitted, inferred from length/stride.
    """
    ptr = lib.sw_getBuff(name.encode("utf-8"))
    if not ptr:
        return None
    size = lib.sw_getBuffSize(name.encode("utf-8"))
    if size <= 0:
        return None
    if length > 0 and size < length:
        raise ValueError(f"Buffer {name} smaller than expected: {size} < {length}")
    if shape is None:
        if stride > 1:
            n = size // stride
            shape = (n, stride)
        else:
            shape = (size,)
    arr = np.ctypeslib.as_array(
        ptr, shape=(shape[0] * shape[1] if len(shape) == 2 else shape[0],)
    )
    return arr.reshape(shape)


# --- C API Function Definitions

# World Management
lib.sw_init.argtypes = []
lib.sw_init.restype = None

def init():
    """Initializes the simulation world."""
    lib.sw_init()
    lib.sw_init_buffers()
    lib.sw_init_buffers()

lib.sw_clear.argtypes = []
lib.sw_clear.restype = None
def clear():
    """Clears all state from the simulation world."""
    lib.sw_clear()
    lib.sw_init_buffers()
    lib.sw_init_buffers()

lib.sw_set_debug.argtypes = [c_int, c_int]
lib.sw_set_debug.restype = None
def set_debug(verbosity: int = -1, idebug: int = -1):
    """Sets verbosity/idebug flags in the native library.

    Args:
        verbosity: Non-negative to set, or -1 to leave unchanged.
        idebug: Non-negative to set, or -1 to leave unchanged.
    """
    lib.sw_set_debug(verbosity, idebug)


lib.sw_reserve_bodies.argtypes = [c_int, c_int]
lib.sw_reserve_bodies.restype = None
def reserve_bodies(n_planets: int, n_ships: int):
    """Pre-allocates memory for planets and ships to avoid vector reallocations."""
    lib.sw_reserve_bodies(n_planets, n_ships)


# Object Management
lib.sw_add_planet.argtypes = [c_char_p, c_double, c_double, c_double_p, c_double_p]
lib.sw_add_planet.restype = c_int
def add_planet(name: str, mass: float, radius: float, pos: tuple, vel: tuple) -> int:
    """Adds a planet to the simulation."""
    pos_arr = (c_double * 3)(*pos)
    vel_arr = (c_double * 3)(*vel)
    return lib.sw_add_planet(name.encode("utf-8"), mass, radius, pos_arr, vel_arr)


lib.sw_add_ship.argtypes = [c_char_p, c_double, c_double, c_double_p, c_double_p]
lib.sw_add_ship.restype = c_int
def add_ship(name: str, mass: float, radius: float, pos: tuple, vel: tuple) -> int:
    """Adds a ship to the simulation."""
    pos_arr = (c_double * 3)(*pos)
    vel_arr = (c_double * 3)(*vel)
    return lib.sw_add_ship(name.encode("utf-8"), mass, radius, pos_arr, vel_arr)


lib.sw_get_n_planets.argtypes = []
lib.sw_get_n_planets.restype = c_int
def get_n_planets() -> int:
    """Returns the number of planets in the simulation."""
    return lib.sw_get_n_planets()


lib.sw_get_n_ships.argtypes = []
lib.sw_get_n_ships.restype = c_int
def get_n_ships() -> int:
    """Returns the number of ships in the simulation."""
    return lib.sw_get_n_ships()


# Simulation
lib.sw_allocate_trjs.argtypes = [c_int]
lib.sw_allocate_trjs.restype = None
def allocate_trjs(n: int):
    """Allocates memory for trajectories with n points."""
    lib.sw_allocate_trjs(n)


lib.sw_predict_trjs.argtypes = [c_int, c_double]
lib.sw_predict_trjs.restype = None
def predict_trjs(n: int, dt: float):
    """Runs the simulation to predict trajectories."""
    lib.sw_predict_trjs(n, dt)


lib.sw_set_ship_thrust.argtypes = [c_int, c_int, array1d, array1d]
lib.sw_set_ship_thrust.restype = None
def set_ship_thrust(ship_idx: int, ts: np.ndarray, thrusts: np.ndarray):
    """Sets the thrust profile for a ship from a series of time-stamped thrust vectors."""
    n_points = len(ts)
    # Ensure thrusts is a flat array of doubles (x1,y1,z1, x2,y2,z2, ...)
    thrusts_flat = np.ascontiguousarray(thrusts.flatten(), dtype=np.double)
    ts_arr = np.ascontiguousarray(ts, dtype=np.double)
    lib.sw_set_ship_thrust(ship_idx, n_points, ts_arr, thrusts_flat)


# Data Access
lib.sw_get_trj_pos.argtypes = [c_int, c_int, array1d, c_int]
lib.sw_get_trj_pos.restype = c_int
def get_trj_pos(body_type: int, body_idx: int, max_points: int) -> np.ndarray:
    """Gets the trajectory positions for a body.

    Args:
        body_type (int): 0 for planet, 1 for ship.
        body_idx (int): The index of the body.
        max_points (int): The maximum number of points to retrieve.

    Returns:
        np.ndarray: A (n, 3) numpy array of positions, or an empty array if not found.
    """
    out_pos = np.zeros((max_points, 3), dtype=np.double)
    n_copied = lib.sw_get_trj_pos(body_type, body_idx, out_pos.flatten(), max_points)
    return out_pos[:n_copied]


# Buffer registry helpers -------------------------------------------------


lib.sw_getBuffSize.argtypes = [c_char_p]
lib.sw_getBuffSize.restype = c_int


lib.sw_inertial_transform.argtypes = [c_double_p, c_double_p]
lib.sw_inertial_transform.restype = None

lib.sw_init_buffers.argtypes = []
lib.sw_init_buffers.restype = None
def init_buffers():
    lib.sw_init_buffers()


lib.sw_getBuff.argtypes = [c_char_p]
lib.sw_getBuff.restype = ctypes.POINTER(c_double)
def get_buffer(name: str, components: int = 3) -> np.ndarray | None:
    n = lib.sw_get_trj_len()
    return _wrap_ptr(name, n * components, stride=components, shape=(n, components))


lib.sw_get_trj_dt.argtypes = []
lib.sw_get_trj_dt.restype = c_double
lib.sw_get_trj_len.argtypes = []
lib.sw_get_trj_len.restype = c_int
def get_trj_metadata() -> tuple[int, float]:
    return lib.sw_get_trj_len(), lib.sw_get_trj_dt()


def inertial_transform(pos_shift: np.ndarray | None, vel_shift: np.ndarray | None):
    pos_ptr = (
        None
        if pos_shift is None
        else np.ascontiguousarray(pos_shift, dtype=np.double).ctypes.data_as(c_double_p)
    )
    vel_ptr = (
        None
        if vel_shift is None
        else np.ascontiguousarray(vel_shift, dtype=np.double).ctypes.data_as(c_double_p)
    )
    lib.sw_inertial_transform(pos_ptr, vel_ptr)


# Scenario setup -----------------------------------------------------------


def setup_jupiter_moons_skirmish():
    """Replicates the 3v3 Jupiter moon scenario from the C++ SpaceTactics constructor."""

    init()

    # Reserve memory to prevent vector reallocation and the subsequent shallow-copy bug
    reserve_bodies(6, 6)

    vJup = 13.1e3
    rJup = 7.784120e11

    # Planets (Sun + Jupiter system)
    add_planet("Sun", 1.99e30, 696.0e6, (0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
    add_planet("Jupiter", 1.90e27, 71.4920e6, (0.0, rJup, 0.0), (vJup, 0.0, 0.0))
    add_planet("Io", 8.93190e23, 1.83e6, (4.21700e8, rJup, 0.0), (vJup, 17330.0, 0.0))
    add_planet("Europa", 4.8e23, 1.5608e6, (6.71034e8, rJup, 0.0), (vJup, 13740.0, 0.0))
    add_planet(
        "Ganymede", 1.48190e24, 2.6312e6, (1.07041e9, rJup, 0.0), (vJup, 10880.0, 0.0)
    )
    add_planet(
        "Calisto", 1.07590e24, 2.4103e6, (1.88271e9, rJup, 0.0), (vJup, 8204.0, 0.0)
    )

    # Ships for faction A (relative to Ganymede)
    ganymede_pos = np.array([1.07041e9, rJup, 0.0])
    ganymede_vel = np.array([vJup, 10880.0, 0.0])
    reference_pos = ganymede_pos
    reference_vel = ganymede_vel
    offsets = np.array([[-4e6, 0.0, 0.0], [-3e6, 0.0, 0.0], [-5e6, 0.0, 0.0]])
    vel_offsets = np.array([[0.0, -5.6e3, 0.0]] * 3)
    ship_names = []
    for i, (pos_off, vel_off) in enumerate(zip(offsets, vel_offsets), start=1):
        add_ship(
            f"ShipA{i}",
            1.0e4,
            2.5e3,
            tuple(reference_pos + pos_off),
            tuple(reference_vel + vel_off),
        )
        ship_names.append(f"ShipA{i}")

    # Ships for faction B (relative to Jupiter)
    jupiter_pos = np.array([0.0, rJup, 0.0])
    jupiter_vel = np.array([vJup, 0.0, 0.0])
    reference_pos = jupiter_pos
    reference_vel = jupiter_vel
    offsets = np.array(
        [
            [0.0, -4.5e8, 0.0],
            [0.0, -4.6e8, 0.0],
            [0.0, -4.7e8, 0.0],
        ]
    )
    vel_offsets = np.array([[20.6e3, 0.0, 0.0]] * 3)
    for i, (pos_off, vel_off) in enumerate(zip(offsets, vel_offsets), start=1):
        add_ship(
            f"ShipB{i}",
            1.0e4,
            2.5e3,
            tuple(reference_pos + pos_off),
            tuple(reference_vel + vel_off),
        )
        ship_names.append(f"ShipB{i}")

    # Center everything on Jupiter (matches C++ inertial transform)
    inertial_transform(-jupiter_pos, -jupiter_vel)

    planet_names = ["Sun", "Jupiter", "Io", "Europa", "Ganymede", "Calisto"]

    return {
        "n_planets": get_n_planets(),
        "n_ships": get_n_ships(),
        "planets": planet_names,
        "ships": ship_names,
    }


def get_ship_trj_buffer(
    index: int, name: str, kind: str = "trjPos"
) -> np.ndarray | None:
    return get_buffer(f"ship[{index}]:{name}.{kind}")


def get_planet_trj_buffer(
    index: int, name: str, kind: str = "trjPos"
) -> np.ndarray | None:
    return get_buffer(f"planet[{index}]:{name}.{kind}")


lib.sw_get_trj_pos.argtypes = [c_int, c_int, array1d, c_int]
lib.sw_get_trj_pos.restype = c_int


def get_trj_pos(body_type: int, body_idx: int, max_points: int) -> np.ndarray:
    """Gets the trajectory positions for a body.

    Args:
        body_type (int): 0 for planet, 1 for ship.
        body_idx (int): The index of the body.
        max_points (int): The maximum number of points to retrieve.

    Returns:
        np.ndarray: A (n, 3) numpy array of positions, or an empty array if not found.
    """
    out_pos = np.zeros((max_points, 3), dtype=np.double)
    n_copied = lib.sw_get_trj_pos(body_type, body_idx, out_pos.flatten(), max_points)
    return out_pos[:n_copied]


# Buffer registry helpers -------------------------------------------------


lib.sw_getBuffSize.argtypes = [c_char_p]
lib.sw_getBuffSize.restype = c_int

lib.sw_get_trj_dt.argtypes = []
lib.sw_get_trj_dt.restype = c_double

lib.sw_inertial_transform.argtypes = [c_double_p, c_double_p]
lib.sw_inertial_transform.restype = None

lib.sw_init_buffers.argtypes = []
lib.sw_init_buffers.restype = None


def init_buffers():
    lib.sw_init_buffers()


lib.sw_getBuff.argtypes = [c_char_p]
lib.sw_getBuff.restype = ctypes.POINTER(c_double)


def get_buffer(name: str, components: int = 3) -> np.ndarray | None:
    n = lib.sw_get_trj_len()
    return _wrap_ptr(name, n * components, stride=components, shape=(n, components))


lib.sw_get_trj_len.argtypes = []
lib.sw_get_trj_len.restype = c_int


def get_trj_metadata() -> tuple[int, float]:
    return lib.sw_get_trj_len(), lib.sw_get_trj_dt()


def inertial_transform(pos_shift: np.ndarray | None, vel_shift: np.ndarray | None):
    pos_ptr = (
        None
        if pos_shift is None
        else np.ascontiguousarray(pos_shift, dtype=np.double).ctypes.data_as(c_double_p)
    )
    vel_ptr = (
        None
        if vel_shift is None
        else np.ascontiguousarray(vel_shift, dtype=np.double).ctypes.data_as(c_double_p)
    )
    lib.sw_inertial_transform(pos_ptr, vel_ptr)


# Scenario setup -----------------------------------------------------------


def setup_jupiter_moons_skirmish():
    """Replicates the 3v3 Jupiter moon scenario from the C++ SpaceTactics constructor."""

    init()

    # Reserve memory to prevent vector reallocation and the subsequent shallow-copy bug
    reserve_bodies(6, 6)

    vJup = 13.1e3
    rJup = 7.784120e11

    # Planets (Sun + Jupiter system)
    add_planet("Sun", 1.99e30, 696.0e6, (0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
    add_planet("Jupiter", 1.90e27, 71.4920e6, (0.0, rJup, 0.0), (vJup, 0.0, 0.0))
    add_planet("Io", 8.93190e23, 1.83e6, (4.21700e8, rJup, 0.0), (vJup, 17330.0, 0.0))
    add_planet("Europa", 4.8e23, 1.5608e6, (6.71034e8, rJup, 0.0), (vJup, 13740.0, 0.0))
    add_planet(
        "Ganymede", 1.48190e24, 2.6312e6, (1.07041e9, rJup, 0.0), (vJup, 10880.0, 0.0)
    )
    add_planet(
        "Calisto", 1.07590e24, 2.4103e6, (1.88271e9, rJup, 0.0), (vJup, 8204.0, 0.0)
    )

    # Ships for faction A (relative to Ganymede)
    ganymede_pos = np.array([1.07041e9, rJup, 0.0])
    ganymede_vel = np.array([vJup, 10880.0, 0.0])
    reference_pos = ganymede_pos
    reference_vel = ganymede_vel
    offsets = np.array([[-4e6, 0.0, 0.0], [-3e6, 0.0, 0.0], [-5e6, 0.0, 0.0]])
    vel_offsets = np.array([[0.0, -5.6e3, 0.0]] * 3)
    ship_names = []
    for i, (pos_off, vel_off) in enumerate(zip(offsets, vel_offsets), start=1):
        add_ship(
            f"ShipA{i}",
            1.0e4,
            2.5e3,
            tuple(reference_pos + pos_off),
            tuple(reference_vel + vel_off),
        )
        ship_names.append(f"ShipA{i}")

    # Ships for faction B (relative to Jupiter)
    jupiter_pos = np.array([0.0, rJup, 0.0])
    jupiter_vel = np.array([vJup, 0.0, 0.0])
    reference_pos = jupiter_pos
    reference_vel = jupiter_vel
    offsets = np.array(
        [
            [0.0, -4.5e8, 0.0],
            [0.0, -4.6e8, 0.0],
            [0.0, -4.7e8, 0.0],
        ]
    )
    vel_offsets = np.array([[20.6e3, 0.0, 0.0]] * 3)
    for i, (pos_off, vel_off) in enumerate(zip(offsets, vel_offsets), start=1):
        add_ship(
            f"ShipB{i}",
            1.0e4,
            2.5e3,
            tuple(reference_pos + pos_off),
            tuple(reference_vel + vel_off),
        )
        ship_names.append(f"ShipB{i}")

    # Center everything on Jupiter (matches C++ inertial transform)
    inertial_transform(-jupiter_pos, -jupiter_vel)

    planet_names = ["Sun", "Jupiter", "Io", "Europa", "Ganymede", "Calisto"]

    return {
        "n_planets": get_n_planets(),
        "n_ships": get_n_ships(),
        "planets": planet_names,
        "ships": ship_names,
    }


def get_ship_trj_buffer(
    index: int, name: str, kind: str = "trjPos"
) -> np.ndarray | None:
    return get_buffer(f"ship[{index}]:{name}.{kind}")


def get_planet_trj_buffer(
    index: int, name: str, kind: str = "trjPos"
) -> np.ndarray | None:
    return get_buffer(f"planet[{index}]:{name}.{kind}")
