

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

#include "Object3D.h"
#include "Terrain25D.h"
#include "Warrior3D.h"
#include "Warrior25D.h"
#include "Projectile3D.h"

// class NonInertWorld : public NBodyWorld2D {       //  TODO:  We would like to add this later
class Shooter {
	public:

	int debug_shit=19;

    int perFrame = 10;
    double dt    = 0.005d;

    double restitution = -0.8d;
    double airDrag     = -0.05d;
    double landDrag    = -0.5d;
    double gravity     = -9.81;
    double objR        = 10.0d;
    double ground_height = 0.0;

    Vec3d wind_speed,watter_speed;

    double projLifetime = 10.0;

    int shotsCount=0,warriorCount=0;

	std::vector<Warrior3D*>     warriors;
	std::vector<Warrior25D*>    warriors25D;
	std::vector<Projectile3D*>  projectiles;  // see http://stackoverflow.com/questions/11457571/how-to-set-initial-size-of-stl-vector
	std::vector<Object3D*>      objects;

	Terrain25D * terrain = NULL;

    // ==== function declarations

    virtual void update_warriors3D();
    virtual void update_warriors25D();
    virtual void update_projectiles3D();

    virtual void update_world( );
    virtual void init_world  ( );

    //virtual void drawEnvironment();

    Projectile3D* fireProjectile( Warrior3D * w, double speed, int kind );
    Warrior3D*    makeWarrior   ( const Vec3d& pos, const Vec3d& dir, const Vec3d& up, int kind );
    int registrWarrior( Warrior3D * w );

};

#endif
