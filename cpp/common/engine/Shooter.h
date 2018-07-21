

#ifndef Shooter_h
#define Shooter_h

#include <vector>
//#include <unordered_set>
//#include <algorithm>

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "Body.h"

//#include "NBodyWorld2D.h"       // TODO: We would like to add this later

#include "Interfaces.h"

#include "Object3D.h"
#include "Terrain25D.h"
#include "Warrior3D.h"
#include "Warrior25D.h"
#include "Projectile3D.h"

// class NonInertWorld : public NBodyWorld2D {       //  TODO:  We would like to add this later
class Shooter {
	public:

	int debug_shit=19;

    int perFrame  = 10;
    double dt     = 0.005d;
    double inv_dt = 1/dt;
    double time   = 0.0d;

    double ground_height = 0.0d; // only if there is no terrain

    double airDrag     = -0.05d;
    double gravity     = -9.81;
    double objR        = 10.0d;

    Vec3d wind_speed,watter_speed;

    double projLifetime = 10.0; //
    int nMaxBurst       = 10;   // when there is more shots in burst, new burst is created
    int nBurstReserve   = 10;   // optimization - we avoid re-allocation of Burst3d::shots when new shots are added

    int shotsCount=0,warriorCount=0;

    // -- types
    std::vector<ProjectileType*> projectile_types;
    std::vector<ObjectType*>     objectTypes;
    std::vector<VehicleType*>    vehicleTypes;

    // -- Objects
	std::vector<Warrior3D*>      warriors;     // DEPRECATED : use Vehicle3d instead
	std::vector<Warrior25D*>     warriors25D;  // TODO : is it worth it? used just for optimization in some Tactics games ?
	std::vector<Projectile3D*>   projectiles;  // Reserve some size? see http://stackoverflow.com/questions/11457571/how-to-set-initial-size-of-stl-vector
	//std::vector<Object3D*>       objects;      // DEPRECATED : use Object3d instead
	std::vector<AnyControler*>   controlers;

	std::vector<Object3d*>  objects;
	std::vector<Vehicle3d*> vehicles;
	std::vector<Burst3d*>   bursts;

	Terrain25D * terrain = NULL;

    // -- Temporaries
    //std::vector<Vec3d> tmpPos; // We probably don't need this if we use

    // ==== function declarations

    void burstObjectCollision( Object3d& obj, Burst3d& burst );

    virtual void update_warriors3D();
    virtual void update_warriors25D();
    virtual void update_projectiles3D();
    virtual void update_bursts();

    virtual void update_world( );
    virtual void init_world  ( );

    //virtual void drawEnvironment();

    Projectile3D* fireProjectile( Warrior3D * w, double speed, int kind );
    Warrior3D*    makeWarrior   ( const Vec3d& pos, const Vec3d& dir, const Vec3d& up, int kind );
    int registrWarrior( Warrior3D * w );

};

#endif
