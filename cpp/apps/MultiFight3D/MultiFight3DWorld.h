

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

    double restitution = -0.8d;
    double airDrag     = -0.05d;
    double landDrag    = -0.5d;
    double objR        = 2.0d;

    double reload_time = 0.4;

    int defaultWarriorShape, defaultObjectShape, defaultObjectHitShape, defaultProjectileShape;

	std::vector<Warrior3D*>    warriors;
	std::vector<Projectile3D*> projectiles;  // see http://stackoverflow.com/questions/11457571/how-to-set-initial-size-of-stl-vector

    std::vector<KinematicBody*> objects;

    // ==== function declarations

    virtual void update_world( );
    virtual void init_world  ( );

    //virtual void drawEnvironment();

    Warrior3D* makeWarrior( const Vec3d& pos, const Vec3d& dir, const Vec3d& Up, int shape );
    void fireProjectile( Warrior3D * w );

    // ==== inline functions

    inline void addEnviroForces( const Vec3d& pos, const Vec3d& vel, Vec3d& fout, bool landed ){
    }

    inline bool collideWithObject(  const Vec3d& objpos, double objR, const Vec3d& pos, Vec3d& vel, Vec3d& normal ){
        Vec3d d; d.set_sub( pos, objpos );
        double r2 = d.norm2( );
        if( r2 < sq(objR) ){
            double r     = sqrt(r2);
            double inv_r = 1/( r + 1.0e-8 );
            normal.set_mul( d, inv_r );
            //printf( "normal (%3.3f,%3.3f)\n", normal.x, normal.y );
            //if( r > (objR+0.5) ){
                double vnormal = normal.dot( vel );
                //printf( "vel (%3.3f,%3.3f,%3.3f) normal (%3.3f,%3.3f,%3.3f) vnormal %3.3f \n" );
                if ( vnormal < 0 ){
                    vel.add_mul( normal, (restitution-1)*vnormal );
                    //vel.add_mul( normal, -vnormal );
                    //vel.mul( restitution );
                }
            //}
            //if( norm.dot_perp( pos ) > 1.0e-4 ) normal.set_mul( pos, pos.norm() ); // we set normal only if it changed more than 1e-4 to save computation of sqrt
            return true;
        }
        return false;
    }

};

#endif
