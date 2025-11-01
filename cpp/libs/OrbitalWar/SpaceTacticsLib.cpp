#include <vector>
#include <string>
#include <unordered_map>
#include <sstream>
#include <cstring>

// Core simulation logic is in header files
#include "SpaceWorld.h"
#include "spline_hermite.h"

// Singleton instance of the simulation world
static SpaceWorld* W = nullptr;

// Buffer registry mirroring LandCraft pattern
static std::unordered_map<std::string, double*> g_buffers;
static std::unordered_map<std::string, int>     g_buffer_sizes;

static void sw_clear_buffers(){
    g_buffers.clear();
    g_buffer_sizes.clear();
}

static void sw_register_body_buffers(const std::string& prefix, const SpaceBody& body){
    if(!W) return;
    const int n = W->trj_n;
    if(n <= 0) return;
    // Trajectory positions (Vec3d => double[3])
    if(body.trjPos){
        std::string key = prefix + ".trjPos";
        g_buffers[key] = (double*)body.trjPos;
        g_buffer_sizes[key] = n * 3;
    }
    if(body.trjThrust){
        std::string key = prefix + ".trjThrust";
        g_buffers[key] = (double*)body.trjThrust;
        g_buffer_sizes[key] = n * 3;
    }
}

static void sw_register_buffers(){
    if(!W) return;
    sw_clear_buffers();
    const int n_planets = W->planets.size();
    const int n_ships   = W->ships.size();
    for(int i=0; i<n_planets; ++i){
        const SpaceBody& body = W->planets[i];
        std::ostringstream ss;
        ss << "planet[" << i << "]:" << body.name;
        sw_register_body_buffers(ss.str(), body);
    }
    for(int i=0; i<n_ships; ++i){
        const SpaceBody& body = W->ships[i];
        std::ostringstream ss;
        ss << "ship[" << i << "]:" << body.name;
        sw_register_body_buffers(ss.str(), body);
    }
}

// C-style API for Python ctypes
extern "C" {

// ---- Function Declarations

void sw_init();
void sw_clear();
int  sw_add_planet(const char* name, double mass, double radius, double* pos, double* vel);
int  sw_add_ship(const char* name, double mass, double radius, double* pos, double* vel);
int  sw_get_n_planets();
int  sw_get_n_ships();
void sw_allocate_trjs(int n);
void sw_predict_trjs(int n, double dt);
void sw_set_ship_thrust(int ship_idx, int n_points, double* ts, double* thrusts);
int  sw_get_trj_pos(int body_type, int body_idx, double* out_pos, int max_points);
void sw_init_buffers();
double* sw_getBuff(const char* name);
int  sw_getBuffSize(const char* name);
int  sw_get_trj_len();
double sw_get_trj_dt();
void sw_inertial_transform(double* pos_shift, double* vel_shift);

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

int sw_get_n_ships() {
    return W ? W->ships.size() : 0;
}

// ---- Simulation

void sw_allocate_trjs(int n) {
    if (!W) return;
    W->allocateTrjs(n);
    W->allocateODE();
    sw_register_buffers();
}

void sw_predict_trjs(int n, double dt) {
    if (!W) return;
    W->predictTrjs(n, dt);
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

void sw_inertial_transform(double* pos_shift, double* vel_shift){
    if(!W) return;
    Vec3d pos = pos_shift ? *(Vec3d*)pos_shift : Vec3dZero;
    Vec3d vel = vel_shift ? *(Vec3d*)vel_shift : Vec3dZero;
    W->intertialTansform(pos, vel);
}

} // extern "C"