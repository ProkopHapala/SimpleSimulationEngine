#ifndef Projectile2D_h
#define Projectile2D_h

#include "fastmath.h"
#include "Vec2.h"
#include "Body2D.h"

class Projectile2D : public PointBody2D {
	public:
    double time;

    virtual ~Projectile2D(){};
    virtual void move( double dt ){
        time += dt;
        move_PointBody2D(dt);
    }

};

#endif
