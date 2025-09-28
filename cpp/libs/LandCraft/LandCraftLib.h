#ifndef LandCraftLib_h
#define LandCraftLib_h


extern "C" {

// Buffers
void     lc_init_buffers();
double*  lc_getBuff(const char* name);
int*     lc_getIBuff(const char* name);

// World
void lc_world_init(const char* dataFolder);

// Map & Terrain
void lc_map_init(int nx, int ny);
void lc_generate_terrain(unsigned int seed, double maxHeight);
int  lc_save(const char* ground_path, const char* water_path);
int  lc_load(const char* ground_path, const char* water_path);

// Hydraulics 2D
void   lc_relax_all();
void   lc_relax_hex(int ix, int iy);
void   lc_set_outflow_at(int ix, int iy);
void   lc_set_inflow_at(int ix, int iy, double delta);
double lc_gather_rain(double minSinkFlow);
int    lc_find_all_rivers(double minFlow);

// Rivers data
int lc_rivers_count();
int lc_river_length(int river_id);
int lc_river_get_path(int river_id, int* out_idx, int n);
int lc_river_get_flow(int river_id, double* out_flow, int n);

// Droplet tracing
int lc_trace_droplet(int ix, int iy, int* out_idx, int max_len);

// Roads
int  lc_road_build_straight(int ax, int ay, int bx, int by);
int  lc_road_length(int road_id);
int  lc_road_get_path_xy(int road_id, unsigned short* out_xy_pairs, int max_pairs);
int  lc_road_profile_heights(int road_id, double* out_ground, double* out_water, int nmax);
void lc_roads_clear();

// Vehicles
int  lc_vehicle_type_create_default();
int  lc_vehicle_spawn_on_road(int road_id, int type_id);
void lc_vehicle_step_all(double dt);
int  lc_vehicle_get_state(int vid, int* ipath, int* idir, int* onWay);

// Economy
int    lc_econ_load_technologies(const char* fname);
int    lc_econ_get_tech_count();
int    lc_econ_get_tech_name(int tech_id, char* out, int out_len);
int    lc_factory_create();
int    lc_factory_set_technology(int fid, int tech_id);
int    lc_factory_set_stock(int fid, const char* commodity, double amount);
int    lc_factory_get_stock(int fid, const char* commodity, double* amount);
double lc_factory_produce(int fid, double N);

// PathFinder
int  lc_pf_bind_to_map();
void lc_pf_set_cost_params(double ch2, double chminus, double chplus);
void lc_pf_clear_centers();
void lc_pf_add_center(int ix, int iy);
void lc_pf_prepare();
int  lc_pf_step();
int  lc_pf_find_connections();
int  lc_pf_make_paths();
int  lc_pf_get_num_paths();
int  lc_pf_get_path_length(int path_id);
int  lc_pf_get_path(int path_id, int* out_idx, int maxn);



} // extern "C"

#endif // LandCraftLib_h