
#ifndef Warrior25D_h
#define Warrior25D_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "Body2D.h"
#include "Mesh.h"

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

    // ==== functions

    virtual void fromString( char const* str );
    virtual void update( double dt, const Vec3d& wind_speed ) = 0;

};

#endif
