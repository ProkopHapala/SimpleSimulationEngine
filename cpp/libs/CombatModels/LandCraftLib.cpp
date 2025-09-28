#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <unordered_map>
#include <string>

#include "fastmath.h"
#include "Vec2.h"
#include "IO_utils.h"

#include "Grid2DAlgs.h"
#include "TerrainHydraulics.h"
#include "PathFinder.h"
#include "../apps/LandCraft/Roads.h"
#include "../apps/LandCraft/Economy.h"
#include "../apps/LandCraft/LandCraftWorld.h"

// Minimal, test-focused C API for LandCraft subsystems.
// Pattern follows eFF_lib.cpp and LandCombatLib.cpp: extern "C", init_buffers/getBuff exposure.

// Dedicated C-layer instance of the world for ctypes users
LandCraftWorld W;
// C-layer readiness flag
// Interface-only buffer maps (for NumPy/ctypes)
static std::unordered_map<std::string, double*> g_buffers;
static std::unordered_map<std::string, int*>    g_ibuffers;

// -------------------- Helpers --------------------

static void lc_init_buffers_impl(){
    g_buffers["ground"] = W.hydro.ground;
    g_buffers["water" ] = W.hydro.water;
}

static void lc_clear_roads(){ W.roadsClear(); }


int* hi=0;


// -------------------- API --------------------
extern "C"{

// ---- World ----
void lc_world_init(const char* dataFolder){ 
    // no buffering stdout and stderr for printf not using iostream pure C not C++ std
    //not this shit: std::cout.rdbuf();std::cerr.rdbuf();
    setbuf(stdout, NULL);
    setbuf(stderr, NULL);

    //printf("hi[10] %i", hi[10]);

    W.init(dataFolder); 
}

// ---- Buffers / Pointers ----
void lc_init_buffers(){ lc_init_buffers_impl(); }

double* lc_getBuff(const char* name){ auto it=g_buffers.find(name); return (it==g_buffers.end())?nullptr:it->second; }
int*    lc_getIBuff(const char* name){ auto it=g_ibuffers.find(name); return (it==g_ibuffers.end())?nullptr:it->second; }

// ---- Map & Terrain ----
// Initialize grid and allocate ground/water arrays
void lc_map_init(int nx, int ny){ W.allocMap(Vec2i{nx,ny}); lc_init_buffers_impl(); }
// Procedural terrain generation (simple: bisecNoise + scale); seed is optional
void lc_generate_terrain(unsigned int seed, double maxHeight){ W.generateTerrain(seed,maxHeight); }

// Flexible terrain helpers
void lc_make_terrain_bisec(int seed){ W.makeTerrainBisec(seed); }
void lc_droplet_erosion(int niter, int nDrops, int nStepMax, int margin, double erodeMin, double erodeMax, double erodeProb){
    W.droplerErosion(niter,nDrops,nStepMax,margin,erodeMin,erodeMax,erodeProb);
}
void lc_set_neighbors(int n){ W.setNeighbors(n); }

// Save/Load terrain arrays (raw binary)
int lc_save(const char* ground_path, const char* water_path){ return W.save(ground_path, water_path); }

int lc_load(const char* ground_path, const char* water_path){ return W.load(ground_path, water_path); }

// ---- Hydraulics 2D ----
void lc_relax_all(){ W.relaxAll(); }
void lc_relax_hex(int ix, int iy){ W.relaxHex(ix,iy); }

void lc_set_outflow_at(int ix, int iy){ W.setOutflowAt(ix,iy); }

void lc_set_inflow_at(int ix, int iy, double delta){ W.setInflowAt(ix,iy,delta); }

double lc_gather_rain(double minSinkFlow){ return W.gatherRain(minSinkFlow); }
int     lc_find_all_rivers(double minFlow){ return W.findAllRivers(minFlow); }

// River data retrieval
int lc_rivers_count(){ return W.riverCount(); }
int lc_river_length(int river_id){ return W.riverLength(river_id); }
int lc_river_get_path(int river_id, int* out_idx, int n){ return W.riverGetPath(river_id,out_idx,n); }
int lc_river_get_flow(int river_id, double* out_flow, int n){ return W.riverGetFlow(river_id,out_flow,n); }

// Droplet tracing
int lc_trace_droplet(int ix, int iy, int* out_idx, int max_len){ return W.traceDroplet(ix,iy,out_idx,max_len); }

// ---- Roads ----
int lc_road_build_straight(int ax, int ay, int bx, int by){ return W.makeRoadStraight({ax,ay},{bx,by}); }

int lc_road_length(int road_id){ return W.roadLength(road_id); }

int lc_road_get_path_xy(int road_id, uint16_t* out_xy_pairs, int max_pairs){ return W.roadGetPathXY(road_id,out_xy_pairs,max_pairs); }

int lc_road_profile_heights(int road_id, double* out_ground, double* out_water, int nmax){ return W.roadProfileHeights(road_id,out_ground,out_water,nmax); }

void lc_roads_clear(){ lc_clear_roads(); }

// ---- Vehicles ----
int lc_vehicle_type_create_default(){ return W.vehicleTypeCreateDefault(); }

int lc_vehicle_spawn_on_road(int road_id, int type_id){ return W.vehicleSpawnOnRoad(road_id,type_id); }

void lc_vehicle_step_all(double dt){ W.vehicleStepAll(dt); }

int lc_vehicle_get_state(int vid, int* ipath, int* idir, int* onWay){ return W.vehicleGetState(vid,ipath,idir,onWay); }

// ---- Economy ----
int lc_econ_load_technologies(const char* fname){ return W.econLoadTechnologies(fname); }

int lc_econ_get_tech_count(){ return W.econGetTechCount(); }

int lc_econ_get_tech_name(int tech_id, char* out, int out_len){ return W.econGetTechName(tech_id,out,out_len); }

int lc_factory_create(){ return W.factoryCreate(); }

int lc_factory_set_technology(int fid, int tech_id){ return W.factorySetTechnology(fid,tech_id); }

int lc_factory_set_stock(int fid, const char* commodity, double amount){ return W.factorySetStock(fid,commodity,amount); }

int lc_factory_get_stock(int fid, const char* commodity, double* amount){ return W.factoryGetStock(fid,commodity,amount); }

double lc_factory_produce(int fid, double N){ return W.factoryProduce(fid,N); }

// ---- PathFinder ----
// Bind to current map and prepare internal buffers
int lc_pf_bind_to_map(){ return W.BindPathFinderToMap(); }

void lc_pf_set_cost_params(double ch2, double chminus, double chplus){ W.pf.setCostParams(ch2,chminus,chplus); }

void lc_pf_clear_centers(){ W.pf.clearCenters(); }
void lc_pf_add_center(int ix, int iy){ W.pf.addCenter(ix,iy); }

void lc_pf_prepare(){ W.pf.pepare(); }

int lc_pf_step(){ return W.pf.stepDelta(); }

int lc_pf_find_connections(){ W.pf.findConnections(); return (int)W.pf.pass.size(); }
int lc_pf_make_paths(){ W.pf.makePaths(); return (int)W.pf.paths.size(); }

int lc_pf_get_num_paths(){ return W.pf.getNumPaths(); }
int lc_pf_get_path_length(int path_id){ return W.pf.getPathLength(path_id); }
int lc_pf_get_path(int path_id, int* out_idx, int maxn){ return W.pf.getPath(path_id,out_idx,maxn); }

} // extern "C"
