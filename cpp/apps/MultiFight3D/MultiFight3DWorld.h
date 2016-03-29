

#ifndef MultiFight3DWorld_h
#define MultiFight3DWorld_h

#include <vector>
//#include <unordered_set>
//#include <algorithm>

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "Body.h"

//#include "NBodyWorld2D.h"       // TODO: We would like to add this later

#include "Warrior3D.h"
#include "Projectile3D.h"

// class NonInertWorld : public NBodyWorld2D {       //  TODO:  We would like to add this later
class MultiFight3DWorld {
	public:

    int perFrame = 10;
    double dt    = 0.005d;

    int defaultWarriorShape, defaultObjectShape;

	std::vector<Warrior3D*>    warriors;
	std::vector<Projectile3D*> projectiles;  // see http://stackoverflow.com/questions/11457571/how-to-set-initial-size-of-stl-vector

    std::vector<KinematicBody*> objects;

    // ==== function declarations

    virtual void update_world( );
    virtual void init_world  ( );

    //virtual void drawEnvironment();

    //void makeWarrior   ( const Vec2d& pos, double angle, char * filename, int shape );
    //void fireProjectile( Warrior2D * w );

    // ==== inline functions

    inline void addEnviroForces( const Vec3d& pos, const Vec3d& vel, Vec3d& fout, bool landed ){
    }

    inline bool collideWithWorld( const Vec3d& pos, Vec3d& vel, Vec3d& normal ){
        return false;
    }



};

#endif
