
#ifndef Warrior25D_h
#define Warrior25D_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "Body2D.h"
#include "Mesh.h"
#include "Projectile3D.h"

class Warrior25D : public RigidBody2D {
	public:

	int id=-1,kind=-1,shape=0;
    bool    landed  = false;
    bool    trigger = false;

    double power;
    Vec3d  span;
    Mesh * mesh;

    double  hground  = 0.5;
    double  throttle = 0.0;

    Mat3d rot3d;
    Vec3d pos3d;

    // ==== functions

    inline void sync3D(){
        //rotmat.set( {rot.x,rot.y,0.0}, {rot.y,-rot.x,0.0}, {1.0,0.0,0.0} );
        //pos3d .set(  pos.x,pos.y,0                                       );
        rot3d.set( {rot.x,0,rot.y}, {0.0,1.0,0.0}, {rot.y,0.0,-rot.x} );
        pos3d .set(  pos.x,0,pos.y  );
    }

    virtual void fromString( char const* str );
    virtual void update( double dt, const Vec3d& wind_speed, const Vec3d& watter_speed ) = 0;

    virtual int shoot( std::vector<Projectile3D*>& projectiles ){};

};

#endif
