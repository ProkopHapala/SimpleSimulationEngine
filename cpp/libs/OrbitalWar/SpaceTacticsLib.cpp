#include <vector>
#include <string>
#include <unordered_map>
#include <sstream>
#include <cstring>
#include <cstdio>

#if __has_include(<sanitizer/asan_interface.h>)
#include <sanitizer/asan_interface.h>
#define HAVE_ASAN_INTERFACE 1
#endif

// Core simulation logic is in header files
#include "SpaceWorld.h"
#include "spline_hermite.h"
#include "globals.h"

// Singleton instance of the simulation world
static SpaceWorld* W = nullptr;

// Buffer registry mirroring LandCraft pattern
static std::unordered_map<std::string, double*> g_buffers;
static std::unordered_map<std::string, int>     g_buffer_sizes;

static void sw_clear_buffers(){
    g_buffers.clear();
    g_buffer_sizes.clear();
}

static void sw_register_body_buffers(const char* prefix, const SpaceBody& body){
    if(!W) return;
    const int n = W->trj_n;
    if(n <= 0 || !prefix) return;
    const size_t prefix_len = strlen(prefix);
    if(verbosity>2){
        printf("    sw_register_body_buffers enter prefix_ptr=%p text='%s' len=%zu\n", (const void*)prefix, prefix, prefix_len);
#ifdef HAVE_ASAN_INTERFACE
        printf("      ASAN prefix poisoned? %d\n", __asan_address_is_poisoned(prefix));
        printf("      ASAN body name poisoned? %d\n", __asan_address_is_poisoned(body.name.c_str()));
#endif
    }
    // Trajectory positions (Vec3d => double[3])
    if(body.trjPos){
        std::string key(prefix);
        key += ".trjPos";
        if(verbosity>2){
            printf("    key=%s body=%p name_ptr=%p size=%zu\n", key.c_str(), (const void*)&body, (const void*)body.name.c_str(), body.name.size());
#ifdef HAVE_ASAN_INTERFACE
            printf("      ASAN key poisoned? %d\n", __asan_address_is_poisoned(key.c_str()));
#endif
        }
        g_buffers[key] = (double*)body.trjPos;
        g_buffer_sizes[key] = n * 3;
    }
    if(body.trjThrust){
        std::string key(prefix);
        key += ".trjThrust";
        if(verbosity>2){
            printf("    key=%s body=%p name_ptr=%p size=%zu\n", key.c_str(), (const void*)&body, (const void*)body.name.c_str(), body.name.size());
        }
        g_buffers[key] = (double*)body.trjThrust;
        g_buffer_sizes[key] = n * 3;
    }
    if(verbosity>2){
        printf("    sw_register_body_buffers exit prefix_ptr=%p text='%s'\n", (const void*)prefix, prefix);
#ifdef HAVE_ASAN_INTERFACE
        printf("      ASAN prefix poisoned (exit)? %d\n", __asan_address_is_poisoned(prefix));
#endif
    }
}

static void sw_register_buffers(){
    if(!W) return;
    if(verbosity>0){
        printf("sw_register_buffers: trj_n=%d planets=%zu ships=%zu\n", W->trj_n, W->planets.size(), W->ships.size());
    }
    sw_clear_buffers();
    const int n_planets = W->planets.size();
    const int n_ships   = W->ships.size();
    for(int i=0; i<n_planets; ++i){
        const SpaceBody& body = W->planets[i];
        char prefix[128];
        snprintf(prefix, sizeof(prefix), "planet[%d]:%s", i, body.name.c_str());
        if(verbosity>1){
            printf("  registering %s trjPos=%p trjThrust=%p prefix_ptr=%p size=%zu\n", body.name.c_str(), (void*)body.trjPos, (void*)body.trjThrust, (const void*)prefix, strlen(prefix));
        }
        sw_register_body_buffers(prefix, body);
        if(verbosity>2){
            printf("  prefix post planet[%d] buf_ptr=%p text='%s' len=%zu\n", i, (const void*)prefix, prefix, strlen(prefix));
#ifdef HAVE_ASAN_INTERFACE
            printf("    ASAN prefix poisoned (post)? %d\n", __asan_address_is_poisoned(prefix));
#endif
        }
    }
    for(int i=0; i<n_ships; ++i){
        const SpaceBody& body = W->ships[i];
        char prefix[128];
        snprintf(prefix, sizeof(prefix), "ship[%d]:%s", i, body.name.c_str());
        if(verbosity>1){
            printf("  registering %s trjPos=%p trjThrust=%p prefix_ptr=%p size=%zu\n", body.name.c_str(), (void*)body.trjPos, (void*)body.trjThrust, (const void*)prefix, strlen(prefix));
        }
        sw_register_body_buffers(prefix, body);
        if(verbosity>2){
            printf("  prefix post ship[%d] buf_ptr=%p text='%s' len=%zu\n", i, (const void*)prefix, prefix, strlen(prefix));
#ifdef HAVE_ASAN_INTERFACE
            printf("    ASAN prefix poisoned (post)? %d\n", __asan_address_is_poisoned(prefix));
#endif
        }
    }
}

// C-style API for Python ctypes
extern "C" {
// ---- World Management

void sw_init() {
    if (!W) {
        W = new SpaceWorld();
    }
    sw_clear_buffers();
}

void sw_clear() {
    if (W) {
        delete W;
        W = new SpaceWorld(); // Re-initialize to a clean state
    }
    sw_clear_buffers();
}

void sw_set_debug(int verbosity_, int idebug_){
    if(verbosity_>=0) verbosity = verbosity_;
    if(idebug_>=0)    idebug    = idebug_;
    if(verbosity>0){
        printf("sw_set_debug: verbosity=%d idebug=%d\n", verbosity, idebug);
    }
}

// ---- Object Management

int sw_add_planet(const char* name, double mass, double radius, double* pos, double* vel) {
    if (!W) return -1;
    W->addPlanet(std::string(name), mass, radius, *(Vec3d*)pos, *(Vec3d*)vel);
    return W->planets.size() - 1;
}

int sw_add_ship(const char* name, double mass, double radius, double* pos, double* vel) {
    if (!W) return -1;
    W->addShip(std::string(name), mass, radius, *(Vec3d*)pos, *(Vec3d*)vel);
    return W->ships.size() - 1;
}

int sw_get_n_planets() {
    return W ? W->planets.size() : 0;
}

int sw_get_n_ships() { return W ? W->ships.size() : 0; }



void sw_inertial_transform(double* pos_shift, double* vel_shift){
    if(!W) return;
    Vec3d pos = pos_shift ? *(Vec3d*)pos_shift : Vec3dZero;
    Vec3d vel = vel_shift ? *(Vec3d*)vel_shift : Vec3dZero;
    W->intertialTansform(pos, vel);
}

void sw_reserve_bodies(int n_planets, int n_ships){
    if (!W) return;
    if (n_planets > 0) {
        W->planets.reserve(n_planets);
    }
    if (n_ships > 0) {
        W->ships.reserve(n_ships);
    }
}

// ---- Simulation

void sw_allocate_trjs(int n) {
    if (!W) return;
    W->allocateTrjs(n);
    W->allocateODE();
    if(verbosity>0){
        printf("sw_allocate_trjs: n=%d planets=%zu ships=%zu masses=%p\n", n, W->planets.size(), W->ships.size(), (void*)W->masses);
        if(verbosity>1){
            for(size_t i=0;i<W->planets.size();++i){
                const SpaceBody& b = W->planets[i];
                printf("    planet[%zu]=%s body=%p name_ptr=%p trjPos=%p trjThrust=%p\n", i, b.name.c_str(), (const void*)&b, (const void*)b.name.c_str(), (void*)b.trjPos, (void*)b.trjThrust);
            }
            for(size_t i=0;i<W->ships.size();++i){
                const SpaceBody& b = W->ships[i];
                printf("    ship[%zu]=%s body=%p name_ptr=%p trjPos=%p trjThrust=%p\n", i, b.name.c_str(), (const void*)&b, (const void*)b.name.c_str(), (void*)b.trjPos, (void*)b.trjThrust);
            }
        }
    }
    sw_register_buffers();
}

void sw_predict_trjs(int n, double dt) {
    if (!W) return;
    W->predictTrjs(n, dt);
    if(verbosity>1){
        printf("sw_predict_trjs: n=%d dt=%g trj_n=%d masses=%p\n", n, dt, W->trj_n, (void*)W->masses);
    }
    sw_register_buffers();
}

void sw_set_ship_thrust(int ship_idx, int n_points, double* ts, double* thrusts) {
    if (!W || ship_idx < 0 || ship_idx >= W->ships.size()) return;
    SpaceBody& ship = W->ships[ship_idx];
    if (!ship.trjThrust || W->trj_n <= 0) return; // Must be allocated first

    // nonUni2spline expects Vec3d*, so we cast the double*
    nonUni2spline(0.0, W->trj_dt, n_points, ts, (Vec3d*)thrusts, W->trj_n, ship.trjThrust);
}


// ---- Data Access

int sw_get_trj_pos(int body_type, int body_idx, double* out_pos, int max_points) {
    if (!W) return 0;

    SpaceBody* body = nullptr;
    if (body_type == 0) { // Planet
        if (body_idx < 0 || body_idx >= W->planets.size()) return 0;
        body = &W->planets[body_idx];
    } else if (body_type == 1) { // Ship
        if (body_idx < 0 || body_idx >= W->ships.size()) return 0;
        body = &W->ships[body_idx];
    } else {
        return 0;
    }

    if (!body->trjPos) return 0;

    int n_to_copy = W->trj_n;
    if (n_to_copy > max_points) {
        n_to_copy = max_points;
    }

    memcpy(out_pos, body->trjPos, n_to_copy * sizeof(Vec3d));

    return n_to_copy;
}

void sw_init_buffers();
double* sw_getBuff(const char* name);
int  sw_getBuffSize(const char* name);
int  sw_get_trj_len();
double sw_get_trj_dt();

void sw_init_buffers(){
    sw_register_buffers();
}

double* sw_getBuff(const char* name){
    if(!name) return nullptr;
    auto it = g_buffers.find(name);
    return (it==g_buffers.end())?nullptr:it->second;
}

int sw_getBuffSize(const char* name){
    if(!name) return 0;
    auto it = g_buffer_sizes.find(name);
    return (it==g_buffer_sizes.end())?0:it->second;
}

int sw_get_trj_len(){
    return W ? W->trj_n : 0;
}

double sw_get_trj_dt(){
    return W ? W->trj_dt : 0.0;
}

} // extern "C"