#ifndef LandCraftWorld_h
#define LandCraftWorld_h

#include <string>
#include <vector>
#include <unordered_map>

#include "TerrainHydraulics.h"
#include "PathFinder.h"
#include "Roads.h"
#include "Economy.h"

// Centralized game state for LandCraft. Holds terrain/hydraulics, rivers, roads, vehicles, economy, and pathfinder.
// A single instance (W) is shared by the app and the C API for consistent testing and debugging.
class LandCraftWorld{
public:
    // Core simulation subsystems
    HydraulicGrid2D hydro;
    std::vector<River*>* rivers = nullptr; // alias to hydro.rivers
    double maxHeight = 500.0; // terrain height scale

    // World state / integration flags
    bool ready = false; // set true when map allocated

    // // Named buffer exposure (for C API, tools, tests)
    // std::unordered_map<std::string,double*> buffers;
    // std::unordered_map<std::string,int*>    ibuffers;

    // Transport
    std::vector<Road*>           roads;
    std::vector<RoadVehicleType*> vehicleTypes;
    std::vector<RoadVehicle*>     vehicles;

    // Reusable builders/helpers
    RoadBuilder rb;

    // Economy
    std::vector<Technology*>                techs;
    std::unordered_map<std::string,int>     techIndex;
    std::vector<Factory*>                   factories;

    // Pathfinding
    PathFinder pf;

    // Data path
    std::string dataPath;

    // Initialize world (load tech tables etc.). Data folder can be used to load files.
    inline void init(const char* dataFolder){
        if(dataFolder) dataPath = dataFolder; else dataPath.clear();
        rivers = nullptr;
        // Note: defer heavy loads to explicit calls from app/script (e.g., econ_load)
    }

    // ===== Initialization API (to be called by app/script) =====
    // TODO: Make these config-driven similar to LTWorld (see apps/LandTactics/LTWorld.*)
    //       e.g., parse dataPath + "/world.ini" or Lua for steps and assets.

    void loadTechnologies(const char* fname);

    // Terrain/Hydraulics
    void makeMap(int sz, double step, bool newMap);
    void makeMapCached(int sz, double step, bool newMap);
    void generateTerrain();
    void mapInit(int nx, int ny){ makeMap(nx, (double)0.0, false); /* step handled by GUI; keep grid only */ }
    void generateTerrain(unsigned int seed, double maxHeight);
    int  save(const char* ground_path, const char* water_path);
    int  load(const char* ground_path, const char* water_path);
    bool loadTerrainDefault();
    bool saveTerrainDefault();
    void relaxAll(){ if(ready) hydro.relaxWater(); }
    void relaxHex(int ix, int iy){ if( ready && in_bounds(ix,iy) ) hydro.relaxWater({ix,iy}); }
    void setOutflowAt(int ix, int iy);
    void setInflowAt(int ix, int iy, double delta);
    double gatherRain(double minSinkFlow){ return hydro.gatherRain(minSinkFlow); }
    int    findAllRivers(double minFlow){ return hydro.findAllRivers(minFlow); }
    int    riverCount() const { return (int)hydro.rivers.size(); }
    int    riverLength(int river_id) const;
    int    riverGetPath(int river_id, int* out_idx, int n) const;
    int    riverGetFlow(int river_id, double* out_flow, int n) const;
    int    traceDroplet(int ix, int iy, int* out_idx, int max_len);

    void makeRivers();

    // ---------- Transport ----------
    int  makeRoadStraight(Vec2i p1, Vec2i p2);
    int  roadLength(int road_id) const;
    int  roadGetPathXY(int road_id, uint16_t* out_xy_pairs, int max_pairs) const;
    int  roadProfileHeights(int road_id, double* out_ground, double* out_water, int nmax) const;
    void roadsClear();
    int  makeVehicles();

    // ---------- Vehicles: world-facing API ----------
    int  vehicleTypeCreateDefault();
    int  vehicleSpawnOnRoad(int road_id, int type_id);
    void vehicleStepAll(double dt);
    int  vehicleGetState(int vid, int* ipath, int* idir, int* onWay) const;

    // ---------- Economy / Factories ----------
    int  econLoadTechnologies(const char* fname);
    int  econGetTechCount() const { return (int)techs.size(); }
    int  econGetTechName(int tech_id, char* out, int out_len) const;
    int  factoryCreate();
    int  factorySetTechnology(int fid, int tech_id);
    int  factorySetStock(int fid, const char* commodity, double amount);
    int  factoryGetStock(int fid, const char* commodity, double* amount) const;
    double factoryProduce(int fid, double N);


    inline bool in_bounds(int ix, int iy){ return (ix>=0) && (ix<hydro.n.x) && (iy>=0) && (iy<hydro.n.y); }
    inline int  ip2i_safe(int ix, int iy){ return iy*hydro.n.x + ix; }

    // ---------- Pathfinder ----------
    int  BindPathFinderToMap();
    // void pfSetCostParams(double ch2, double chminus, double chplus);      // TODO/DEBUG deprecated
    // void pfClearCenters();                                                // TODO/DEBUG deprecated
    // void pfAddCenter(int ix, int iy);                                     // TODO/DEBUG deprecated
    // void pfPrepare();                                                     // TODO/DEBUG deprecated
    // int  pfStep();                                                        // TODO/DEBUG deprecated
    // int  pfFindConnections();                                             // TODO/DEBUG deprecated
    // int  pfMakePaths();                                                   // TODO/DEBUG deprecated
    // int  pfGetNumPaths() const;                                           // TODO/DEBUG deprecated
    // int  pfGetPathLength(int path_id) const;                              // TODO/DEBUG deprecated
    // int  pfGetPath(int path_id, int* out_idx, int maxn) const;            // TODO/DEBUG deprecated

    // ===== Pure simulation steps (no rendering) =====
    // TODO: expose stats outputs (E_before/after, wsum_before/after) via return struct if needed by GUI
    void hydroRelaxStep();
    void vehiclesStep(double dt); // alias used by app; calls vehicleStepAll

    // ===== Config-driven initialization =====
    // Returns true on success, false if file missing or parse error. Supported keys (space-separated):
    //   techs.file <path>
    //   map.size <int>
    //   map.step <double>
    //   terrain.useCache <0|1>
    //   terrain.generate <0|1>
    //   rivers.min_flow <double>
    //   roads.straight <x0> <y0> <x1> <y1>
    //   vehicles.default <0|1>
    bool loadConfig(const char* fname);
};



#endif // LandCraftWorld_h
