
#ifndef Warrior2D_h
#define Warrior2D_h

#include "fastmath.h"
#include "Vec2.h"
#include "Body2D.h"

class Warrior2D : public RigidBody2D {
	public:
    bool  landed = false;
    Vec2d surf;

    void tryJump(){
        printf( "Warrior2D::tryJump : %i (%3.3f,%3.3f) \n", landed, pos.x, pos.y );
        if( landed ){
            vel.add_mul( surf, 5.0 );
        }
    }

};

#endif
