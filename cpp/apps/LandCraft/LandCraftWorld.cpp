#include "LandCraftWorld.h"
#include "IO_utils.h"
#include "testUtils.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <string>
#include <cstring>
#include <cstdint>
#include <unistd.h> // getcwd

// NOTE: This file moves non-GUI initialization from the app into the shared world.
// TODO: Make all these steps config-driven (similar to LTWorld), e.g. parse
//       dataPath + "/world.ini" or a Lua script that lists initialization steps.
// ---------- Economy ----------
void LandCraftWorld::loadTechnologies(const char* fname){
    printf("LandCraftWorld::loadTechnologies(%s)\n", fname);
    const char* path = fname;
    FILE* pFile = fopen(path, "r" );
    if (!pFile && dataPath.size()>0){
        std::string full = dataPath + "/" + fname;
        path = full.c_str();
        pFile = fopen(path, "r");
    }
    if (!pFile){ printf("[LandCraftWorld] Unable to open tech file: %s\n", path); return; }
    while(true){
        Technology* tech = new Technology();
        if( tech->fromFile(pFile) ){
            tech->print();
            techIndex[tech->name] = (int)techs.size();
            techs.push_back(tech);
        }else{ delete tech; break; }
    }
    fclose(pFile);
}

bool LandCraftWorld::loadConfig(const char* fname){
    printf("LandCraftWorld::loadConfig(%s)\n", fname);
    // If caller passed a relative path with a directory (e.g., "data/world.ini"), adopt that directory as dataPath.
    if(fname){
        const char* slash = strrchr(fname, '/');
        if(slash){
            std::string dir = std::string(fname, slash - fname);
            if(dir.size()>0 && dir[0] != '/'){
                char cwd[1024]; if(getcwd(cwd, sizeof(cwd))){ dataPath = std::string(cwd) + "/" + dir; } else { dataPath = dir; }
            }else{
                dataPath = dir;
            }
        }
    }
    // Resolve path without hardcoding: use the path as given
    const char* path = fname;
    std::string full;
    FILE* f = nullptr;
    f = fopen(path, "r");
    if(!f){
        printf("LandCraftWorld::loadConfig(%s): ERROR: config not found: %s\n", fname, path);
        char cwd[1024];
        if(getcwd(cwd, sizeof(cwd))){ printf("Current working directory: %s\n", cwd); }
        else{ perror("getcwd"); }
        exit(0);
    }
    // state
    int mapSize = -1; double mapStep = 50.0; bool haveMap=false;
    double rivers_min_flow = -1.0;

    char buf[2048];
    while(fgets(buf, sizeof(buf), f)){
        // strip comments
        char* p = buf;
        // trim leading spaces
        while(*p==' '||*p=='\t') p++;
        if(*p=='#' || *p==';' || *p=='\n' || *p=='\0') continue;
        // read key
        char key[64]; int nconsumed=0;
        if(sscanf(p, "%63s %n", key, &nconsumed) < 1) continue;
        p += nconsumed;

        // Strict, compact command list (one line = one command)
        if( strcmp(key, "techs")==0 ){
            char rel[512]; if(sscanf(p, "%511s", rel)==1) loadTechnologies(rel);
        }else if( strcmp(key, "terrain")==0 ){
            int sz; double step; if(sscanf(p, "%d %lf", &sz, &step)==2){
                makeMapCached(sz,step,false);
                haveMap=true; mapSize=sz; mapStep=step;
            }
        }else if( strcmp(key, "rivers")==0 ){
            if(sscanf(p, "%lf", &rivers_min_flow)==1){ double wmax = hydro.gatherRain(100.0); (void)wmax; hydro.findAllRivers(rivers_min_flow); }
        }else if( strcmp(key, "straight_road")==0 ){
            int x0,y0,x1,y1; if(sscanf(p, "%d %d %d %d", &x0,&y0,&x1,&y1)==4){ makeRoadStraight({x0,y0},{x1,y1}); }
        }else if( strcmp(key, "vehicles")==0 ){
            int n=1; char type[64]={0}; int nread = sscanf(p, "%d %63s", &n, type);
            if(n<1) n=1; for(int i=0;i<n;i++) makeVehicles();
        }else if( strcmp(key, "load_terrain")==0 ){
            if(!ready){ printf("[LandCraftWorld] load_terrain requested before terrain allocation; ignoring.\n"); }
            else{ if(!loadTerrainDefault()){ printf("[LandCraftWorld] load_terrain failed; keeping current.\n"); } }
        }else if( strcmp(key, "save_terrain")==0 ){
            if(!ready){ printf("[LandCraftWorld] save_terrain requested before terrain allocation; ignoring.\n"); }
            else{ if(!saveTerrainDefault()){ printf("[LandCraftWorld] save_terrain failed.\n"); } }
        }
    }
    if(!haveMap && mapSize>0){ allocMap(mapSize,mapStep); haveMap=true; }
    fclose(f);
    return true;
}

void LandCraftWorld::allocMap(int sz, double step ){
    printf("LandCraftWorld::makeMap(%d, %f)\n", sz, step);
    // World does not own the GUI ruler; we only size and initialize hydraulics here.
    hydro.allocate( {sz,sz} );
    hydro.initNeighs_6(false);
    hydro.allocate_outflow();
    ready = true;
    // buffer maps are owned by the C API layer; nothing to do here
}

void LandCraftWorld::makeMapCached(int sz, double step, bool newMap){
    allocMap(sz,step);
    bool loaded=false;
    if(!newMap){ loaded = loadTerrainDefault(); }
    if(!loaded){
        printf("[LandCraftWorld] Cached terrain missing/corrupt -> generating new terrain and saving cache.\n");
        generateTerrain();
        saveTerrainDefault();
    }
}

void LandCraftWorld::generateTerrain(){
    printf("LandCraftWorld::generateTerrain()\n");
    // Simple synthetic terrain
    srand(16464);
    hydro.ground[0]=0.2;
    bisecNoise( 7, hydro.ground, -1.0/256, 1.0/256 );
    hydro.initNeighs_6(false);
    for( int j=0; j<500; j++ ){
        int isz = 25;
        int ix0 = rand()%(hydro.n.x-isz);
        int iy0 = rand()%(hydro.n.y-isz);
        hydro.errodeDroples( 400, 500, +0.1, 0.15, 0.9, {ix0, iy0}, {ix0+isz, iy0+isz} );
    }
    for(int i=0; i<hydro.ntot; i++){ hydro.ground[i] *= maxHeight; hydro.water[i] = hydro.ground[i]; }
}

void LandCraftWorld::makeRivers(){
    printf("LandCraftWorld::makeRivers()\n");
    double wmax = hydro.gatherRain( 100.0 ); (void)wmax;
    hydro.findAllRivers( 50.0 );
    std::sort(hydro.rivers.begin(), hydro.rivers.end(), [](River* a, River* b){ return a->path.size() > b->path.size(); } );
}

// ---------- Transport ----------
int LandCraftWorld::makeRoadStraight(Vec2i p1, Vec2i p2){
    printf("LandCraftWorld::makeRoadStraight(%d, %d, %d, %d)\n", p1.x, p1.y, p2.x, p2.y);
    Road* road = new Road();
    roads.push_back( road );
    RoadBuilder rb; rb.road = road; rb.path.clear();
    rb.path.push_back( (RoadTile){(uint16_t)p1.x,(uint16_t)p1.y,0.0} );
    rb.pushStright( p2 );
    // Sample heights if needed later; currently not used
    rb.writeIt();
    return (int)roads.size()-1;
}

// int LandCraftWorld::roadBuildStraight(int ax, int ay, int bx, int by){
//     if(!ready) return -1;
//     Road* road = new Road();
//     roads.push_back(road);
//     rb.road = road; rb.path.clear();
//     rb.path.push_back({(uint16_t)ax,(uint16_t)ay,0.0});
//     rb.pushStright({bx,by});
//     for(auto& p: rb.path){ int i = hydro.ip2i({(int)p.ia,(int)p.ib}); p.height = hydro.ground[i]; }
//     rb.writeIt();
//     return (int)roads.size()-1;
// }

int LandCraftWorld::roadLength(int road_id) const{
    if(road_id<0 || road_id>=(int)roads.size()) return -1; return roads[road_id]->n;
}

int LandCraftWorld::roadGetPathXY(int road_id, uint16_t* out_xy_pairs, int max_pairs) const{
    if(road_id<0 || road_id>=(int)roads.size()) return -1;
    Road* r = roads[road_id];
    int n = r->n; if(max_pairs<n) n=max_pairs;
    for(int i=0;i<n;i++){ out_xy_pairs[2*i+0]=r->path[i].ia; out_xy_pairs[2*i+1]=r->path[i].ib; }
    return n;
}

int LandCraftWorld::roadProfileHeights(int road_id, double* out_ground, double* out_water, int nmax) const{
    if(!ready) return -1;
    if(road_id<0 || road_id>=(int)roads.size()) return -1;
    Road* r = roads[road_id];
    int n = r->n; if(nmax<n) n=nmax;
    for(int i=0;i<n;i++){
        int ix = r->path[i].ia; int iy = r->path[i].ib;
        int idx = hydro.ip2i({ix,iy});
        out_ground[i] = hydro.ground[idx];
        out_water [i] = hydro.water [idx];
    }
    return n;
}

void LandCraftWorld::roadsClear(){
    for(Road* r: roads){ delete r; }
    roads.clear();
}

void LandCraftWorld::generateTerrain(unsigned int seed, double maxH){
    if(!ready) return;
    srand(seed);
    int npow = 7;
    hydro.ground[0]=0.2;
    bisecNoise( npow, hydro.ground, -1.0/256, 1.0/256 );
    for(int j=0;j<64;j++){
        int isz=16; int ix0=rand()%(hydro.n.x-isz); int iy0=rand()%(hydro.n.y-isz);
        hydro.errodeDroples(200,300,+0.05,0.10,0.8, {ix0,iy0}, {ix0+isz,iy0+isz});
    }
    for(int i=0;i<hydro.ntot;i++){ hydro.water[i]=hydro.ground[i]*maxH; hydro.ground[i]*=maxH; }
}

int LandCraftWorld::save(const char* ground_path, const char* water_path){
    if(!ready) return -1;
    int er1 = saveBin(ground_path, sizeof(double)*hydro.ntot, (char*)hydro.ground);
    int er2 = saveBin(water_path , sizeof(double)*hydro.ntot, (char*)hydro.water );
    return (er1||er2)?-1:0;
}

int LandCraftWorld::load(const char* ground_path, const char* water_path){
    if(!ready) return -1;
    int er1 = loadBin(ground_path, sizeof(double)*hydro.ntot, (char*)hydro.ground, false);
    int er2 = loadBin(water_path , sizeof(double)*hydro.ntot, (char*)hydro.water , false);
    return (er1||er2)?-1:0;
}

// Resolve default terrain paths using dataPath if set, otherwise CWD data/
static inline std::string lc_join(const std::string& a, const char* b){ return a + "/" + b; }

bool LandCraftWorld::loadTerrainDefault(){
    std::string base = (dataPath.size()>0)? dataPath : std::string("data");
    std::string g = lc_join(base, "ground.bin");
    std::string w = lc_join(base, "water.bin");
    int er = load(g.c_str(), w.c_str());
    return (er==0);
}

bool LandCraftWorld::saveTerrainDefault(){
    std::string base = (dataPath.size()>0)? dataPath : std::string("data");
    std::string g = lc_join(base, "ground.bin");
    std::string w = lc_join(base, "water.bin");
    int er = save(g.c_str(), w.c_str());
    return (er==0);
}

void LandCraftWorld::setOutflowAt(int ix, int iy){
    if(!ready) return;
    for(int i=0;i<hydro.ntot;i++){ hydro.known[i]=false; }
    int ihex = iy*hydro.n.x + ix;
    hydro.contour2[0]=ihex; hydro.nContour=1; hydro.isOutflow=true;
    hydro.water[ihex] = hydro.ground[ihex];
}

void LandCraftWorld::setInflowAt(int ix, int iy, double delta){
    if(!ready) return;
    for(int i=0;i<hydro.ntot;i++){ hydro.known[i]=false; }
    int ihex = iy*hydro.n.x + ix;
    hydro.contour2[0]=ihex; hydro.nContour=1; hydro.isOutflow=false;
    hydro.water[ihex] = hydro.ground[ihex] + delta;
}

int LandCraftWorld::riverLength(int river_id) const{
    if(!ready) return -1;
    if(river_id<0 || river_id>=(int)hydro.rivers.size()) return -1;
    return (int)hydro.rivers[river_id]->path.size();
}

int LandCraftWorld::riverGetPath(int river_id, int* out_idx, int n) const{
    if(!ready) return -1;
    if(river_id<0 || river_id>=(int)hydro.rivers.size()) return -1;
    River* r = hydro.rivers[river_id];
    int m = (int)r->path.size(); if(n<m) m=n;
    for(int i=0;i<m;i++){ out_idx[i]=r->path[i]; }
    return m;
}

int LandCraftWorld::riverGetFlow(int river_id, double* out_flow, int n) const{
    if(!ready) return -1;
    if(river_id<0 || river_id>=(int)hydro.rivers.size()) return -1;
    River* r = hydro.rivers[river_id];
    int m = (int)r->flow.size(); if(n<m) m=n;
    for(int i=0;i<m;i++){ out_flow[i]=r->flow[i]; }
    return m;
}
int LandCraftWorld::traceDroplet(int ix, int iy, int* out_idx, int max_len){
    if(!ready) return -1;
    return hydro.traceDroplet({ix,iy}, max_len, out_idx);
}

int LandCraftWorld::makeVehicles(){
    printf("LandCraftWorld::makeVehicles()\n");
    if(roads.empty()) return -1;
    RoadVehicleType* vtype = new RoadVehicleType();
    vehicleTypes.push_back(vtype);
    RoadVehicle* v = new RoadVehicle();
    v->road = roads[0];
    v->type = vtype;
    vehicles.push_back( v );
    return (int)vehicles.size()-1;
}

// ---------- Vehicles API ----------
int LandCraftWorld::vehicleTypeCreateDefault(){
    printf("LandCraftWorld::vehicleTypeCreateDefault()\n");
    RoadVehicleType* t = new RoadVehicleType();
    vehicleTypes.push_back(t);
    return (int)vehicleTypes.size()-1;
}

int LandCraftWorld::vehicleSpawnOnRoad(int road_id, int type_id){
    printf("LandCraftWorld::vehicleSpawnOnRoad(%d, %d)\n", road_id, type_id);
    if(road_id<0 || road_id>=(int)roads.size()) return -1;
    if(type_id<0 || type_id>=(int)vehicleTypes.size()) return -2;
    RoadVehicle* v = new RoadVehicle();
    v->road = roads[road_id];
    v->type = vehicleTypes[type_id];
    v->ipath = 0; v->idir = 1; v->t_rest = 0.0; v->onWay = true;
    vehicles.push_back(v);
    return (int)vehicles.size()-1;
}

void LandCraftWorld::vehicleStepAll(double dt){ for(auto* v: vehicles){ v->move(dt); } }

int LandCraftWorld::vehicleGetState(int vid, int* ipath, int* idir, int* onWay) const{
    if(vid<0 || vid>=(int)vehicles.size()) return -1;
    RoadVehicle* v = vehicles[vid];
    if(ipath) *ipath = v->ipath;
    if(idir ) *idir  = v->idir;
    if(onWay) *onWay = v->onWay ? 1:0;
    return 0;
}

// ---------- Economy / Factories ----------
int LandCraftWorld::econLoadTechnologies(const char* fname){
    int n0 = (int)techs.size();
    loadTechnologies(fname);
    return (int)techs.size() - n0;
}

int LandCraftWorld::econGetTechName(int tech_id, char* out, int out_len) const{
    if(tech_id<0 || tech_id>=(int)techs.size() || out_len<=0) return -1;
    const std::string& s = techs[tech_id]->name;
    int n = (int)std::min((int)s.size(), out_len-1);
    memcpy(out, s.c_str(), n); out[n]='\0';
    return n;
}

int LandCraftWorld::factoryCreate(){ Factory* f = new Factory(); factories.push_back(f); return (int)factories.size()-1; }

int LandCraftWorld::factorySetTechnology(int fid, int tech_id){
    if(fid<0||fid>=(int)factories.size()) return -1;
    if(tech_id<0||tech_id>=(int)techs.size()) return -2;
    factories[fid]->setTechnology(techs[tech_id]);
    return 0;
}

int LandCraftWorld::factorySetStock(int fid, const char* commodity, double amount){
    if(fid<0||fid>=(int)factories.size()) return -1;
    factories[fid]->stored[ std::string(commodity) ] = amount;
    return 0;
}

int LandCraftWorld::factoryGetStock(int fid, const char* commodity, double* amount) const{
    if(fid<0||fid>=(int)factories.size()) return -1;
    auto& m = factories[fid]->stored;
    auto it = m.find( std::string(commodity) );
    if(it==m.end()) return -2;
    if(amount) *amount = it->second;
    return 0;
}

double LandCraftWorld::factoryProduce(int fid, double N){
    if(fid<0||fid>=(int)factories.size()) return 0.0;
    return factories[fid]->produce(N);
}

// ---------- Pathfinder ----------
int LandCraftWorld::BindPathFinderToMap(){
    if(!ready) return -1;
    pf.set(hydro);
    pf.bind(hydro.ground, nullptr);
    pf.allocate();
    return 0;
}

// ---------- Pure simulation steps (no rendering) ----------
void LandCraftWorld::hydroRelaxStep(){
    // TODO: optionally return (wsum_before, wsum_after, E_before, E_after) for diagnostics
    double wsum_before=0.0, E_before=0.0;
    hydro.sumWater(wsum_before, E_before);
    hydro.relaxWater();
    for(int iy=0; iy<hydro.n.y; iy++){
        hydro.relaxWaterRasterX( iy, 0, hydro.n.x, 1.0 );
        hydro.relaxWaterRasterY( iy, 0, hydro.n.x, 1.0 );
    }
    double wsum_after=0.0, E_after=0.0;
    hydro.sumWater(wsum_after, E_after);
    // Note: caller may check conservation; we keep it pure here.
}

void LandCraftWorld::vehiclesStep(double dt){
    for(RoadVehicle* veh : vehicles){
        // Prefer deterministic advancement; keep path-following logic encapsulated in RoadVehicle
        veh->move(dt);
        if(!veh->onWay) veh->depart();
    }
}
