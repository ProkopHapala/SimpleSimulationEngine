#ifndef LandCraftLua_h
#define LandCraftLua_h

#include "LuaHelpers.h"
#include "../../libs/LandCraft/LandCraftLib.h"

// Skeleton Lua bindings for LandCraft subsystems.
// Modeled on the spacecraft Lua pattern in `cpp/common/Orbital/EditSpaceCraft.h`.
// These functions call the C API from LandCraftLib so they can be used in-game or in tools.

namespace LandCraftLua {

// --- Helpers ---
inline void pushIntArray(lua_State* L, const int* a, int n){
    lua_createtable(L, n, 0);
    for(int i=0;i<n;i++){ lua_pushinteger(L, (lua_Integer)a[i]); lua_rawseti(L, -2, i+1); }
}

inline void pushU16Pairs(lua_State* L, const unsigned short* a, int n_pairs){
    // returns array of {x,y} tables
    lua_createtable(L, n_pairs, 0);
    for(int i=0;i<n_pairs;i++){
        lua_createtable(L, 2, 0);
        lua_pushinteger(L, (lua_Integer)a[2*i+0]); lua_rawseti(L, -2, 1);
        lua_pushinteger(L, (lua_Integer)a[2*i+1]); lua_rawseti(L, -2, 2);
        lua_rawseti(L, -2, i+1);
    }
}

inline void pushDoubleArray(lua_State* L, const double* a, int n){
    lua_createtable(L, n, 0);
    for(int i=0;i<n;i++){ lua_pushnumber(L, a[i]); lua_rawseti(L, -2, i+1); }
}

// Example: l_GenerateTerrain{ seed:int }
inline int l_GenerateTerrain(lua_State* L){
    unsigned int seed = (unsigned int)Lua::getIntField(L, "seed");
    double maxH = Lua::getDoubleField(L, "maxHeight", 500.0);
    lc_generate_terrain(seed, maxH);
    return 0;
}

// Example: l_RelaxAll{ }
inline int l_RelaxAll(lua_State* L){
    (void)L;
    lc_relax_all();
    return 0;
}

// Example: l_RelaxHex{ ix:int, iy:int }
inline int l_RelaxHex(lua_State* L){
    int ix = Lua::getIntField(L, "ix");
    int iy = Lua::getIntField(L, "iy");
    lc_relax_hex(ix,iy);
    return 0;
}

// Example: l_GatherRain{ minFlow:double } -> wmax:double
inline int l_GatherRain(lua_State* L){
    double minFlow = Lua::getDoubleField(L, "minFlow");
    double wmax = lc_gather_rain(minFlow);
    lua_pushnumber(L, wmax);
    return 1;
}

// Example: l_FindRivers{ minFlow:double } -> count:int
inline int l_FindRivers(lua_State* L){
    double minFlow = Lua::getDoubleField(L, "minFlow");
    int n = lc_find_all_rivers(minFlow);
    lua_pushinteger(L, n);
    return 1;
}

// Example: l_RiverPath{ id:int } -> indices:table
inline int l_RiverPath(lua_State* L){
    int id = Lua::getIntField(L, "id");
    int len = lc_river_length(id);
    if(len<=0){ lua_newtable(L); return 1; }
    std::vector<int> tmp(len);
    int m = lc_river_get_path(id, tmp.data(), len);
    pushIntArray(L, tmp.data(), m);
    return 1;
}

// Example: l_TraceDroplet{ ix:int, iy:int, max:int } -> indices:table
inline int l_TraceDroplet(lua_State* L){
    int ix = Lua::getIntField(L, "ix");
    int iy = Lua::getIntField(L, "iy");
    int maxn = Lua::getIntField(L, "max", 256);
    std::vector<int> tmp(maxn);
    int n = lc_trace_droplet(ix,iy,tmp.data(), maxn);
    if(n<0) n=0; pushIntArray(L, tmp.data(), n);
    return 1;
}

// Roads
inline int l_BuildRoad(lua_State* L){
    int ax = Lua::getIntField(L, "ax");
    int ay = Lua::getIntField(L, "ay");
    int bx = Lua::getIntField(L, "bx");
    int by = Lua::getIntField(L, "by");
    int id = lc_road_build_straight(ax,ay,bx,by);
    lua_pushinteger(L, id);
    return 1;
}

inline int l_RoadProfile(lua_State* L){
    int id = Lua::getIntField(L, "id");
    int n  = lc_road_length(id);
    if(n<=0){ lua_newtable(L); lua_newtable(L); return 2; }
    std::vector<double> g(n), w(n);
    int m = lc_road_profile_heights(id, g.data(), w.data(), n);
    if(m<0) m=0; pushDoubleArray(L, g.data(), m); pushDoubleArray(L, w.data(), m);
    return 2;
}

// Registration helper
inline void register_api(lua_State* L){
    lua_register(L, "generate_terrain", l_GenerateTerrain);
    lua_register(L, "relax_all",       l_RelaxAll);
    lua_register(L, "relax_hex",       l_RelaxHex);
    lua_register(L, "gather_rain",     l_GatherRain);
    lua_register(L, "find_rivers",     l_FindRivers);
    lua_register(L, "river_path",      l_RiverPath);
    lua_register(L, "trace_droplet",   l_TraceDroplet);
    lua_register(L, "build_road",      l_BuildRoad);
    lua_register(L, "road_profile",    l_RoadProfile);

    // ---- Hydraulics inflow/outflow/save/load ----
    lua_register(L, "inflow", [](lua_State* L)->int{
        int ix = Lua::getIntField(L, "ix"); int iy = Lua::getIntField(L, "iy"); double d = Lua::getDoubleField(L, "delta", 10.0); lc_set_inflow_at(ix,iy,d); return 0; });
    lua_register(L, "outflow", [](lua_State* L)->int{
        int ix = Lua::getIntField(L, "ix"); int iy = Lua::getIntField(L, "iy"); lc_set_outflow_at(ix,iy); return 0; });
    lua_register(L, "save_terrain", [](lua_State* L)->int{
        const char* g = Lua::getStringField(L, "ground"); const char* w = Lua::getStringField(L, "water"); int er = lc_save(g,w); lua_pushinteger(L, er); return 1; });
    lua_register(L, "load_terrain", [](lua_State* L)->int{
        const char* g = Lua::getStringField(L, "ground"); const char* w = Lua::getStringField(L, "water"); int er = lc_load(g,w); lua_pushinteger(L, er); return 1; });

    // ---- Vehicles ----
    lua_register(L, "vehicle_type_default", [](lua_State* L)->int{ int id = lc_vehicle_type_create_default(); lua_pushinteger(L,id); return 1; });
    lua_register(L, "vehicle_spawn", [](lua_State* L)->int{ int road = Lua::getIntField(L, "road"); int type = Lua::getIntField(L, "type", 0); int id = lc_vehicle_spawn_on_road(road,type); lua_pushinteger(L,id); return 1; });
    lua_register(L, "vehicle_step_all", [](lua_State* L)->int{ double dt = Lua::getDoubleField(L, "dt", 1.0); lc_vehicle_step_all(dt); return 0; });
    lua_register(L, "vehicle_state", [](lua_State* L)->int{ int vid = Lua::getIntField(L, "id"); int ip=0,dir=0,on=0; int er = lc_vehicle_get_state(vid,&ip,&dir,&on); lua_pushinteger(L, er); lua_pushinteger(L, ip); lua_pushinteger(L, dir); lua_pushboolean(L, on!=0); return 4; });

    // ---- Economy ----
    lua_register(L, "econ_load", [](lua_State* L)->int{ const char* f = Lua::getStringField(L, "file"); int n = lc_econ_load_technologies(f); lua_pushinteger(L, n); return 1; });
    lua_register(L, "factory_new", [](lua_State* L)->int{ int id = lc_factory_create(); lua_pushinteger(L,id); return 1; });
    lua_register(L, "factory_set_tech", [](lua_State* L)->int{ int fid = Lua::getIntField(L, "id"); int tid = Lua::getIntField(L, "tech"); int er = lc_factory_set_technology(fid,tid); lua_pushinteger(L, er); return 1; });
    lua_register(L, "factory_set", [](lua_State* L)->int{ int fid = Lua::getIntField(L, "id"); const char* c = Lua::getStringField(L, "commodity"); double a = Lua::getDoubleField(L, "amount"); int er = lc_factory_set_stock(fid,c,a); lua_pushinteger(L, er); return 1; });
    lua_register(L, "factory_get", [](lua_State* L)->int{ int fid = Lua::getIntField(L, "id"); const char* c = Lua::getStringField(L, "commodity"); double a=0; int er = lc_factory_get_stock(fid,c,&a); lua_pushinteger(L, er); lua_pushnumber(L,a); return 2; });
    lua_register(L, "factory_produce", [](lua_State* L)->int{ int fid = Lua::getIntField(L, "id"); double N = Lua::getDoubleField(L, "N"); double m = lc_factory_produce(fid,N); lua_pushnumber(L,m); return 1; });

    // ---- PathFinder ----
    lua_register(L, "pf_bind", [](lua_State* L)->int{ (void)L; int er = lc_pf_bind_to_map(); lua_pushinteger(L, er); return 1; });
    lua_register(L, "pf_params", [](lua_State* L)->int{ double ch2 = Lua::getDoubleField(L, "ch2", 1.0); double chm = Lua::getDoubleField(L, "chminus", 0.0); double chp = Lua::getDoubleField(L, "chplus", 0.0); lc_pf_set_cost_params(ch2,chm,chp); return 0; });
    lua_register(L, "pf_clear_centers", [](lua_State* L)->int{ (void)L; lc_pf_clear_centers(); return 0; });
    lua_register(L, "pf_add_center", [](lua_State* L)->int{ int ix = Lua::getIntField(L, "ix"); int iy = Lua::getIntField(L, "iy"); lc_pf_add_center(ix,iy); return 0; });
    lua_register(L, "pf_prepare", [](lua_State* L)->int{ (void)L; lc_pf_prepare(); return 0; });
    lua_register(L, "pf_step", [](lua_State* L)->int{ (void)L; int d = lc_pf_step(); lua_pushinteger(L,d); return 1; });
    lua_register(L, "pf_find_connections", [](lua_State* L)->int{ (void)L; int n = lc_pf_find_connections(); lua_pushinteger(L,n); return 1; });
    lua_register(L, "pf_make_paths", [](lua_State* L)->int{ (void)L; int n = lc_pf_make_paths(); lua_pushinteger(L,n); return 1; });
    lua_register(L, "pf_paths", [](lua_State* L)->int{ int np = lc_pf_get_num_paths(); lua_createtable(L, np, 0); for(int i=0;i<np;i++){ int len = lc_pf_get_path_length(i); std::vector<int> tmp(len); int m = lc_pf_get_path(i, tmp.data(), len); pushIntArray(L, tmp.data(), m); lua_rawseti(L, -2, i+1);} return 1; });
}

} // namespace LandCraftLua

#endif
