#ifndef Projectile3D_h
#define Projectile3D_h

#include "fastmath.h"
#include "Vec3.h"
#include "Body.h"

class Projectile3D : public PointBody {
	public:
    double time;

    virtual void move( double dt ){
        time += dt;
        move_PointBody(dt);
    };

};

#endif
