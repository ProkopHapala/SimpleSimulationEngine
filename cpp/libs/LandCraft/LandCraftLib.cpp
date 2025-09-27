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

// Minimal, test-focused C API for LandCraft subsystems.
// Pattern follows eFF_lib.cpp and LandCombatLib.cpp: extern "C", init_buffers/getBuff exposure.

// -------------------- Global State --------------------
static HydraulicGrid2D g_hydro;
static bool            g_hydro_ready = false;

static std::vector<River*>* g_rivers = nullptr;  // alias to g_hydro.rivers

static RoadBuilder g_rb;
static std::vector<Road*> g_roads;               // simple registry

static PathFinder  g_pf;                         // not yet wired; placeholder

static std::unordered_map<std::string, double*> g_buffers;
static std::unordered_map<std::string, int*>    g_ibuffers;

// Vehicles
static std::vector<RoadVehicleType*> g_vehicleTypes;
static std::vector<RoadVehicle*>     g_vehicles;

// Economy
static std::vector<Technology*>                 g_techs;
static std::unordered_map<std::string,int>      g_techIndex;
static std::vector<Factory*>                    g_factories;

// -------------------- Helpers --------------------
static inline bool in_bounds(int ix, int iy){ return (ix>=0) && (ix<g_hydro.n.x) && (iy>=0) && (iy<g_hydro.n.y); }
static inline int  ip2i_safe(int ix, int iy){ return iy*g_hydro.n.x + ix; }

static void lc_init_buffers(){
    g_buffers["ground"] = g_hydro.ground;
    g_buffers["water" ] = g_hydro.water;
}

static void lc_clear_roads(){
    for(Road* r : g_roads){ delete r; }
    g_roads.clear();
}

// -------------------- API --------------------
extern "C"{

// ---- Buffers / Pointers ----
void lc_init_buffers(){ lc_init_buffers(); }

double* lc_getBuff(const char* name){ auto it=g_buffers.find(name); return (it==g_buffers.end())?nullptr:it->second; }
int*    lc_getIBuff(const char* name){ auto it=g_ibuffers.find(name); return (it==g_ibuffers.end())?nullptr:it->second; }

// ---- Map & Terrain ----
// Initialize grid and allocate ground/water arrays
void lc_map_init(int nx, int ny){
    g_hydro.allocate({nx,ny});
    g_hydro.initNeighs_6(false);
    g_hydro.allocate_outflow();
    g_hydro_ready = true;
    lc_init_buffers();
}

// Procedural terrain generation (simple: bisecNoise + scale); seed is optional
void lc_generate_terrain(unsigned int seed, double maxHeight){
    if(!g_hydro_ready) return;
    srand(seed);
    // Use bisecNoise on ground heights in place; assumes ntot is 2^k x 2^k or similar
    // Fallback: simply call bisecNoise with a fixed npow estimating from nx
    int npow = 7; // reasonable default; caller can re-run if needed
    g_hydro.ground[0] = 0.2;
    bisecNoise(npow, g_hydro.ground, -1.0/256, 1.0/256);
    // Optionally small erosion batches for variability
    for(int j=0; j<64; j++){
        int isz = 16;
        int ix0 = rand()%(g_hydro.n.x - isz);
        int iy0 = rand()%(g_hydro.n.y - isz);
        g_hydro.errodeDroples(200, 300, +0.05, 0.10, 0.8, {ix0,iy0}, {ix0+isz,iy0+isz});
    }
    for(int i=0;i<g_hydro.ntot;i++){ g_hydro.water[i]=g_hydro.ground[i]*maxHeight; g_hydro.ground[i]*=maxHeight; }
}

// Save/Load terrain arrays (raw binary)
int lc_save(const char* ground_path, const char* water_path){
    if(!g_hydro_ready) return -1;
    int er1 = saveBin(ground_path, sizeof(double)*g_hydro.ntot, (char*)g_hydro.ground);
    int er2 = saveBin(water_path , sizeof(double)*g_hydro.ntot, (char*)g_hydro.water );
    return (er1||er2)?-1:0;
}

int lc_load(const char* ground_path, const char* water_path){
    if(!g_hydro_ready) return -1;
    int er1 = loadBin(ground_path, sizeof(double)*g_hydro.ntot, (char*)g_hydro.ground, false);
    int er2 = loadBin(water_path , sizeof(double)*g_hydro.ntot, (char*)g_hydro.water , false);
    return (er1||er2)?-1:0;
}

// ---- Hydraulics 2D ----
void lc_relax_all(){ if(g_hydro_ready) g_hydro.relaxWater(); }
void lc_relax_hex(int ix, int iy){ if(g_hydro_ready && in_bounds(ix,iy)) g_hydro.relaxWater({ix,iy}); }

void lc_set_outflow_at(int ix, int iy){
    if(!g_hydro_ready || !in_bounds(ix,iy)) return;
    for(int i=0;i<g_hydro.ntot;i++){ g_hydro.known[i]=false; }
    int ihex = ip2i_safe(ix,iy);
    g_hydro.contour2[0]=ihex; g_hydro.nContour=1; g_hydro.isOutflow=true;
    g_hydro.water[ihex] = g_hydro.ground[ihex];
}

void lc_set_inflow_at(int ix, int iy, double delta){
    if(!g_hydro_ready || !in_bounds(ix,iy)) return;
    for(int i=0;i<g_hydro.ntot;i++){ g_hydro.known[i]=false; }
    int ihex = ip2i_safe(ix,iy);
    g_hydro.contour2[0]=ihex; g_hydro.nContour=1; g_hydro.isOutflow=false;
    g_hydro.water[ihex] = g_hydro.ground[ihex] + delta;
}

double lc_gather_rain(double minSinkFlow){ if(!g_hydro_ready) return 0.0; return g_hydro.gatherRain(minSinkFlow); }
int     lc_find_all_rivers(double minFlow){ if(!g_hydro_ready) return 0; return g_hydro.findAllRivers(minFlow); }

// River data retrieval
int lc_rivers_count(){ if(!g_hydro_ready) return 0; return (int)g_hydro.rivers.size(); }
int lc_river_length(int river_id){ if(!g_hydro_ready) return -1; if(river_id<0||river_id>=(int)g_hydro.rivers.size()) return -1; return (int)g_hydro.rivers[river_id]->path.size(); }
int lc_river_get_path(int river_id, int* out_idx, int n){
    if(!g_hydro_ready) return -1; if(river_id<0||river_id>=(int)g_hydro.rivers.size()) return -1;
    River* r = g_hydro.rivers[river_id];
    int m = (int)r->path.size(); if(n<m) m=n;
    for(int i=0;i<m;i++){ out_idx[i]=r->path[i]; }
    return m;
}
int lc_river_get_flow(int river_id, double* out_flow, int n){
    if(!g_hydro_ready) return -1; if(river_id<0||river_id>=(int)g_hydro.rivers.size()) return -1;
    River* r = g_hydro.rivers[river_id];
    int m = (int)r->flow.size(); if(n<m) m=n;
    for(int i=0;i<m;i++){ out_flow[i]=r->flow[i]; }
    return m;
}

// Droplet tracing
int lc_trace_droplet(int ix, int iy, int* out_idx, int max_len){
    if(!g_hydro_ready || !in_bounds(ix,iy)) return -1;
    return g_hydro.traceDroplet({ix,iy}, max_len, out_idx);
}

// ---- Roads ----
int lc_road_build_straight(int ax, int ay, int bx, int by){
    if(!g_hydro_ready || !in_bounds(ax,ay) || !in_bounds(bx,by)) return -1;
    Road* road  = new Road();
    g_roads.push_back(road);
    g_rb.road = road;
    g_rb.path.clear();
    g_rb.path.push_back({(uint16_t)ax,(uint16_t)ay,0.0});
    g_rb.pushStright({bx,by});
    // fill heights from terrain
    for(auto& p : g_rb.path){ int i = g_hydro.ip2i({(int)p.ia,(int)p.ib}); p.height = g_hydro.ground[i]; }
    g_rb.writeIt();
    return (int)g_roads.size()-1;
}

int lc_road_length(int road_id){ if(road_id<0 || road_id>=(int)g_roads.size()) return -1; return g_roads[road_id]->n; }

int lc_road_get_path_xy(int road_id, uint16_t* out_xy_pairs, int max_pairs){
    if(road_id<0 || road_id>=(int)g_roads.size()) return -1;
    Road* r = g_roads[road_id];
    int n = r->n; if(max_pairs<n) n=max_pairs;
    for(int i=0;i<n;i++){ out_xy_pairs[2*i+0]=r->path[i].ia; out_xy_pairs[2*i+1]=r->path[i].ib; }
    return n;
}

int lc_road_profile_heights(int road_id, double* out_ground, double* out_water, int nmax){
    if(!g_hydro_ready) return -1;
    if(road_id<0 || road_id>=(int)g_roads.size()) return -1;
    Road* r = g_roads[road_id];
    int n = r->n; if(nmax<n) n=nmax;
    for(int i=0;i<n;i++){
        int ix = r->path[i].ia; int iy = r->path[i].ib;
        int idx = g_hydro.ip2i({ix,iy});
        out_ground[i] = g_hydro.ground[idx];
        out_water [i] = g_hydro.water [idx];
    }
    return n;
}

void lc_roads_clear(){ lc_clear_roads(); }

// ---- Vehicles ----
int lc_vehicle_type_create_default(){
    RoadVehicleType* t = new RoadVehicleType();
    g_vehicleTypes.push_back(t);
    return (int)g_vehicleTypes.size()-1;
}

int lc_vehicle_spawn_on_road(int road_id, int type_id){
    if(road_id<0 || road_id>=(int)g_roads.size()) return -1;
    if(type_id<0 || type_id>=(int)g_vehicleTypes.size()) return -1;
    RoadVehicle* v = new RoadVehicle();
    v->road = g_roads[road_id];
    v->type = g_vehicleTypes[type_id];
    v->ipath = 0;
    v->idir = 1;
    v->t_rest = 0.0;
    v->onWay = true;
    g_vehicles.push_back(v);
    return (int)g_vehicles.size()-1;
}

void lc_vehicle_step_all(double dt){ for(auto* v: g_vehicles){ v->move(dt); } }

int lc_vehicle_get_state(int vid, int* ipath, int* idir, int* onWay){
    if(vid<0 || vid>=(int)g_vehicles.size()) return -1;
    RoadVehicle* v = g_vehicles[vid];
    if(ipath) *ipath = v->ipath;
    if(idir ) *idir  = v->idir;
    if(onWay) *onWay = v->onWay ? 1:0;
    return 0;
}

// ---- Economy ----
int lc_econ_load_technologies(const char* fname){
    FILE* pFile = fopen(fname, "r");
    if(!pFile) return -1;
    int n=0;
    while(true){
        Technology* tech = new Technology();
        if( tech->fromFile(pFile) ){
            g_techIndex[tech->name] = (int)g_techs.size();
            g_techs.push_back(tech);
            n++;
        }else{ delete tech; break; }
    }
    fclose(pFile);
    return n;
}

int lc_econ_get_tech_count(){ return (int)g_techs.size(); }

int lc_econ_get_tech_name(int tech_id, char* out, int out_len){
    if(tech_id<0 || tech_id>=(int)g_techs.size() || out_len<=0) return -1;
    const std::string& s = g_techs[tech_id]->name;
    int n = (int)std::min((int)s.size(), out_len-1);
    memcpy(out, s.c_str(), n); out[n]='\0';
    return n;
}

int lc_factory_create(){ Factory* f = new Factory(); g_factories.push_back(f); return (int)g_factories.size()-1; }

int lc_factory_set_technology(int fid, int tech_id){
    if(fid<0||fid>=(int)g_factories.size()) return -1;
    if(tech_id<0||tech_id>=(int)g_techs.size()) return -2;
    g_factories[fid]->setTechnology(g_techs[tech_id]);
    return 0;
}

int lc_factory_set_stock(int fid, const char* commodity, double amount){
    if(fid<0||fid>=(int)g_factories.size()) return -1;
    g_factories[fid]->stored[ std::string(commodity) ] = amount;
    return 0;
}

int lc_factory_get_stock(int fid, const char* commodity, double* amount){
    if(fid<0||fid>=(int)g_factories.size()) return -1;
    auto& m = g_factories[fid]->stored;
    auto it = m.find( std::string(commodity) );
    if(it==m.end()) return -2;
    if(amount) *amount = it->second;
    return 0;
}

double lc_factory_produce(int fid, double N){
    if(fid<0||fid>=(int)g_factories.size()) return 0.0;
    return g_factories[fid]->produce(N);
}

// ---- PathFinder ----
// Bind to current map and prepare internal buffers
int lc_pf_bind_to_map(){
    if(!g_hydro_ready) return -1;
    g_pf.set(g_hydro);               // copy grid dims and neighbor topology
    g_pf.bind(g_hydro.ground, nullptr); // default: cost from terrain height only
    g_pf.allocate();
    return 0;
}

void lc_pf_set_cost_params(double ch2, double chminus, double chplus){
    g_pf.ch2 = ch2; g_pf.chminus = chminus; g_pf.chplus = chplus;
}

void lc_pf_clear_centers(){ g_pf.centers.clear(); }
void lc_pf_add_center(int ix, int iy){ g_pf.centers.push_back({ix,iy}); }

void lc_pf_prepare(){ g_pf.pepare(); }

int lc_pf_step(){
    int n0 = g_pf.nContour;
    g_pf.path_step();
    return g_pf.nContour - n0; // frontier delta
}

int lc_pf_find_connections(){ g_pf.findConnections(); return (int)g_pf.pass.size(); }
int lc_pf_make_paths(){ g_pf.makePaths(); return (int)g_pf.paths.size(); }

int lc_pf_get_num_paths(){ return (int)g_pf.paths.size(); }
int lc_pf_get_path_length(int path_id){ if(path_id<0||path_id>=(int)g_pf.paths.size()) return -1; return (int)g_pf.paths[path_id]->path.size(); }
int lc_pf_get_path(int path_id, int* out_idx, int maxn){
    if(path_id<0||path_id>=(int)g_pf.paths.size()) return -1;
    auto* w = g_pf.paths[path_id];
    int m = (int)w->path.size(); if(maxn<m) m=maxn;
    for(int i=0;i<m;i++){ out_idx[i]=w->path[i]; }
    return m;
}

} // extern "C"
