
#ifndef Warrior3D_h
#define Warrior3D_h

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "Body.h"

class Warrior3D : public RigidBody { public:
	int id=-1,kind=-1,shape=0;
    bool    landed  = false;
    bool    trigger = false;
    Vec3d   surf;
    Quat4d  look;
    Vec3d   gun_rot;
    double  until_reaload = 0;
    double  hground       = 0.5;

    void setPose( const Vec3d& pos_, const Vec3d& dir, const Vec3d& up ){
        //w->kind = kind; w->id = warriorCount; warriorCount++;
        initOne();
        pos.set           ( pos_    );
        rotMat.a.set      ( dir     );
        rotMat.b.set      ( up      );
        rotMat.c.set_cross( dir, up );
        qrot.fromMatrix   ( rotMat );
        printf( "pos (%g,%g,%g) qrot (%g,%g,%g,%g)\n", pos.x, pos.x, pos.x, qrot.x,qrot.y,qrot.z,qrot.w );
    }

};

#endif
