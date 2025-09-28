import numpy as np
import ctypes
from ctypes import c_int, c_double, c_bool, c_float, c_char_p, c_void_p, c_uint16
import os, sys

# Reuse cpp_utils loader pattern like in pyMolecular/eFF.py
sys.path.append('../')
from pyMeta import cpp_utils

# Configure build path and load shared library
# Adjust if your build outputs a different name/path
# Use the CombatModels build dir; LandCraftLib target lives there per CMake
cpp_utils.BUILD_PATH = os.path.normpath(cpp_utils.PACKAGE_PATH + '../../../cpp/Build/libs/CombatModels')
lib = cpp_utils.loadLib('LandCraftLib')

# ctypes helpers
c_double_p = ctypes.POINTER(c_double)
c_int_p    = ctypes.POINTER(c_int)
c_uint16_p = ctypes.POINTER(c_uint16)

array1i  = np.ctypeslib.ndpointer(dtype=np.int32,   ndim=1, flags='CONTIGUOUS')
array1u16= np.ctypeslib.ndpointer(dtype=np.uint16,  ndim=1, flags='CONTIGUOUS')
array1d  = np.ctypeslib.ndpointer(dtype=np.double,  ndim=1, flags='CONTIGUOUS')


def _np_as(arr, atype):
    if arr is None:
        return None
    return arr.ctypes.data_as(atype)

# ------------- World -------------
lib.lc_world_init.argtypes = [c_char_p]
lib.lc_world_init.restype  = None
def world_init(data_folder: str | None):
    # Pass NULL when data_folder is None so C++ can decide default behavior
    return lib.lc_world_init(None if data_folder is None else data_folder.encode('utf-8'))

# ------------- Buffers -------------
lib.lc_init_buffers.argtypes = []
lib.lc_init_buffers.restype  = None
def init_buffers():
    return lib.lc_init_buffers()

lib.lc_getBuff.argtypes = [c_char_p]
lib.lc_getBuff.restype  = c_double_p
def getBuff(name: str, shape):
    ptr = lib.lc_getBuff(name.encode('utf-8'))
    if not isinstance(shape, tuple):
        shape = (shape,)
    return np.ctypeslib.as_array(ptr, shape=shape)

lib.lc_getIBuff.argtypes = [c_char_p]
lib.lc_getIBuff.restype  = c_int_p
def getIBuff(name: str, shape):
    ptr = lib.lc_getIBuff(name.encode('utf-8'))
    if not isinstance(shape, tuple):
        shape = (shape,)
    return np.ctypeslib.as_array(ptr, shape=shape)

# ------------- Map & Terrain -------------
lib.lc_map_init.argtypes = [c_int, c_int]
lib.lc_map_init.restype  = None
def map_init(nx: int, ny: int):
    return lib.lc_map_init(nx, ny)

lib.lc_generate_terrain.argtypes = [c_int, c_double]
lib.lc_generate_terrain.restype  = None
def generate_terrain(seed: int, max_height: float):
    return lib.lc_generate_terrain(seed, max_height)

# Flexible terrain helpers
lib.lc_make_terrain_bisec.argtypes = [c_int]
lib.lc_make_terrain_bisec.restype  = None
def make_terrain_bisec(seed: int = 16464):
    return lib.lc_make_terrain_bisec(seed)

lib.lc_droplet_erosion.argtypes = [c_int, c_int, c_int, c_int, c_double, c_double, c_double]
lib.lc_droplet_erosion.restype  = None
def droplet_erosion(niter: int, nDrops: int, nStepMax: int, margin: int, erodeMin: float, erodeMax: float, erodeProb: float):
    return lib.lc_droplet_erosion(niter, nDrops, nStepMax, margin, erodeMin, erodeMax, erodeProb)

lib.lc_set_neighbors.argtypes = [c_int]
lib.lc_set_neighbors.restype  = None
def set_neighbors(n: int):
    """Set neighborhood topology: 4 (von Neumann), 6 (hex-like), or 8 (Moore)."""
    return lib.lc_set_neighbors(n)

lib.lc_save.argtypes = [c_char_p, c_char_p]
lib.lc_save.restype  = c_int
def save(ground_path: str, water_path: str) -> int:
    return lib.lc_save(ground_path.encode('utf-8'), water_path.encode('utf-8'))

# Note: generate_terrain_bisec is not exported by C layer; use make_terrain_bisec + normalize via generate_terrain if needed.
lib.lc_load.argtypes = [c_char_p, c_char_p]
lib.lc_load.restype  = c_int
def load(ground_path: str, water_path: str) -> int:
    return lib.lc_load(ground_path.encode('utf-8'), water_path.encode('utf-8'))

# ------------- Hydraulics 2D -------------
lib.lc_relax_all.argtypes = []
lib.lc_relax_all.restype  = None
def relax_all():
    return lib.lc_relax_all()

lib.lc_relax_hex.argtypes = [c_int, c_int]
lib.lc_relax_hex.restype  = None
def relax_hex(ix: int, iy: int):
    return lib.lc_relax_hex(ix, iy)

lib.lc_set_outflow_at.argtypes = [c_int, c_int]
lib.lc_set_outflow_at.restype  = None
def set_outflow_at(ix: int, iy: int):
    return lib.lc_set_outflow_at(ix, iy)

lib.lc_set_inflow_at.argtypes = [c_int, c_int, c_double]
lib.lc_set_inflow_at.restype  = None
def set_inflow_at(ix: int, iy: int, delta: float):
    return lib.lc_set_inflow_at(ix, iy, delta)

lib.lc_gather_rain.argtypes = [c_double]
lib.lc_gather_rain.restype  = c_double
def gather_rain(min_sink_flow: float) -> float:
    return lib.lc_gather_rain(min_sink_flow)

lib.lc_find_all_rivers.argtypes = [c_double]
lib.lc_find_all_rivers.restype  = c_int
def find_all_rivers(min_flow: float) -> int:
    return lib.lc_find_all_rivers(min_flow)

# River data
lib.lc_rivers_count.argtypes = []
lib.lc_rivers_count.restype  = c_int
def rivers_count() -> int:
    return lib.lc_rivers_count()

lib.lc_river_length.argtypes = [c_int]
lib.lc_river_length.restype  = c_int
def river_length(river_id: int) -> int:
    return lib.lc_river_length(river_id)

lib.lc_river_get_path.argtypes = [c_int, array1i, c_int]
lib.lc_river_get_path.restype  = c_int
def river_get_path(river_id: int, out_idx: np.ndarray) -> int:
    return lib.lc_river_get_path(river_id, out_idx, out_idx.size)

lib.lc_river_get_flow.argtypes = [c_int, array1d, c_int]
lib.lc_river_get_flow.restype  = c_int
def river_get_flow(river_id: int, out_flow: np.ndarray) -> int:
    return lib.lc_river_get_flow(river_id, out_flow, out_flow.size)

# Droplet
lib.lc_trace_droplet.argtypes = [c_int, c_int, array1i, c_int]
lib.lc_trace_droplet.restype  = c_int
def trace_droplet(ix: int, iy: int, out_idx: np.ndarray) -> int:
    return lib.lc_trace_droplet(ix, iy, out_idx, out_idx.size)

# ------------- Roads -------------
lib.lc_road_build_straight.argtypes = [c_int, c_int, c_int, c_int]
lib.lc_road_build_straight.restype  = c_int
def road_build_straight(ax: int, ay: int, bx: int, by: int) -> int:
    return lib.lc_road_build_straight(ax, ay, bx, by)

lib.lc_road_length.argtypes = [c_int]
lib.lc_road_length.restype  = c_int
def road_length(road_id: int) -> int:
    return lib.lc_road_length(road_id)

lib.lc_road_get_path_xy.argtypes = [c_int, array1u16, c_int]
lib.lc_road_get_path_xy.restype  = c_int
def road_get_path_xy(road_id: int, out_xy_pairs: np.ndarray) -> int:
    # out_xy_pairs should be uint16 of length 2*N (ia,ib pairs)
    return lib.lc_road_get_path_xy(road_id, out_xy_pairs, out_xy_pairs.size)

lib.lc_road_profile_heights.argtypes = [c_int, array1d, array1d, c_int]
lib.lc_road_profile_heights.restype  = c_int
def road_profile_heights(road_id: int, out_ground: np.ndarray, out_water: np.ndarray) -> int:
    nmax = min(out_ground.size, out_water.size)
    return lib.lc_road_profile_heights(road_id, out_ground, out_water, nmax)

lib.lc_roads_clear.argtypes = []
lib.lc_roads_clear.restype  = None
def roads_clear():
    return lib.lc_roads_clear()

# ------------- Vehicles -------------
lib.lc_vehicle_type_create_default.argtypes = []
lib.lc_vehicle_type_create_default.restype  = c_int
def vehicle_type_create_default() -> int:
    return lib.lc_vehicle_type_create_default()

lib.lc_vehicle_spawn_on_road.argtypes = [c_int, c_int]
lib.lc_vehicle_spawn_on_road.restype  = c_int
def vehicle_spawn_on_road(road_id: int, type_id: int) -> int:
    return lib.lc_vehicle_spawn_on_road(road_id, type_id)

lib.lc_vehicle_step_all.argtypes = [c_double]
lib.lc_vehicle_step_all.restype  = None
def vehicle_step_all(dt: float):
    return lib.lc_vehicle_step_all(dt)

lib.lc_vehicle_get_state.argtypes = [c_int, c_int_p, c_int_p, c_int_p]
lib.lc_vehicle_get_state.restype  = c_int
def vehicle_get_state(vid: int):
    ipath = c_int()
    idir  = c_int()
    onWay = c_int()
    ret = lib.lc_vehicle_get_state(vid, ctypes.byref(ipath), ctypes.byref(idir), ctypes.byref(onWay))
    return ret, ipath.value, idir.value, (onWay.value != 0)

# ------------- Economy -------------
lib.lc_econ_load_technologies.argtypes = [c_char_p]
lib.lc_econ_load_technologies.restype  = c_int
def econ_load_technologies(fname: str) -> int:
    return lib.lc_econ_load_technologies(fname.encode('utf-8'))

lib.lc_econ_get_tech_count.argtypes = []
lib.lc_econ_get_tech_count.restype  = c_int
def econ_get_tech_count() -> int:
    return lib.lc_econ_get_tech_count()

lib.lc_econ_get_tech_name.argtypes = [c_int, ctypes.c_char_p, c_int]
lib.lc_econ_get_tech_name.restype  = c_int
def econ_get_tech_name(tech_id: int) -> str:
    buf = ctypes.create_string_buffer(256)
    n = lib.lc_econ_get_tech_name(tech_id, buf, ctypes.sizeof(buf))
    if n < 0: return ''
    return buf.value.decode('utf-8')

lib.lc_factory_create.argtypes = []
lib.lc_factory_create.restype  = c_int
def factory_create() -> int:
    return lib.lc_factory_create()

lib.lc_factory_set_technology.argtypes = [c_int, c_int]
lib.lc_factory_set_technology.restype  = c_int
def factory_set_technology(fid: int, tech_id: int) -> int:
    return lib.lc_factory_set_technology(fid, tech_id)

lib.lc_factory_set_stock.argtypes = [c_int, c_char_p, c_double]
lib.lc_factory_set_stock.restype  = c_int
def factory_set_stock(fid: int, commodity: str, amount: float) -> int:
    return lib.lc_factory_set_stock(fid, commodity.encode('utf-8'), amount)

lib.lc_factory_get_stock.argtypes = [c_int, c_char_p, c_double_p]
lib.lc_factory_get_stock.restype  = c_int
def factory_get_stock(fid: int, commodity: str):
    val = c_double()
    ret = lib.lc_factory_get_stock(fid, commodity.encode('utf-8'), ctypes.byref(val))
    return ret, val.value

lib.lc_factory_produce.argtypes = [c_int, c_double]
lib.lc_factory_produce.restype  = c_double
def factory_produce(fid: int, N: float) -> float:
    return lib.lc_factory_produce(fid, N)

# ------------- PathFinder -------------
lib.lc_pf_bind_to_map.argtypes = []
lib.lc_pf_bind_to_map.restype  = c_int
def pf_bind_to_map() -> int:
    return lib.lc_pf_bind_to_map()

lib.lc_pf_set_cost_params.argtypes = [c_double, c_double, c_double]
lib.lc_pf_set_cost_params.restype  = None
def pf_set_cost_params(ch2: float, chminus: float, chplus: float):
    return lib.lc_pf_set_cost_params(ch2, chminus, chplus)

lib.lc_pf_clear_centers.argtypes = []
lib.lc_pf_clear_centers.restype  = None
def pf_clear_centers():
    return lib.lc_pf_clear_centers()

lib.lc_pf_add_center.argtypes = [c_int, c_int]
lib.lc_pf_add_center.restype  = None
def pf_add_center(ix: int, iy: int):
    return lib.lc_pf_add_center(ix, iy)

lib.lc_pf_prepare.argtypes = []
lib.lc_pf_prepare.restype  = None
def pf_prepare():
    return lib.lc_pf_prepare()

lib.lc_pf_step.argtypes = []
lib.lc_pf_step.restype  = c_int
def pf_step() -> int:
    return lib.lc_pf_step()

lib.lc_pf_find_connections.argtypes = []
lib.lc_pf_find_connections.restype  = c_int
def pf_find_connections() -> int:
    return lib.lc_pf_find_connections()

lib.lc_pf_make_paths.argtypes = []
lib.lc_pf_make_paths.restype  = c_int
def pf_make_paths() -> int:
    return lib.lc_pf_make_paths()

lib.lc_pf_get_num_paths.argtypes = []
lib.lc_pf_get_num_paths.restype  = c_int
def pf_get_num_paths() -> int:
    return lib.lc_pf_get_num_paths()

lib.lc_pf_get_path_length.argtypes = [c_int]
lib.lc_pf_get_path_length.restype  = c_int
def pf_get_path_length(path_id: int) -> int:
    return lib.lc_pf_get_path_length(path_id)

lib.lc_pf_get_path.argtypes = [c_int, array1i, c_int]
lib.lc_pf_get_path.restype  = c_int
def pf_get_path(path_id: int, out_idx: np.ndarray) -> int:
    return lib.lc_pf_get_path(path_id, out_idx, out_idx.size)

