
#ifndef Warrior3D_h
#define Warrior3D_h

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "Body.h"
#include "Terrain25D.h" // if we don't use functions from Terrain25D we do not need to link it against Terrain25D.o or Terrain25D.cpp

class Warrior3D : public RigidBody { public:
	int id=-1,kind=-1,shape=0;
    bool    landed  = false;
    bool    trigger = false;
    Vec3d   surf;
    Quat4d  look;
    Vec3d   gun_rot = (Vec3d){1.0,0.0,0.0};
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

    virtual void interact( Terrain25D * terrain ){
        Vec2d dv;
        double h  = terrain->eval( {pos.x,pos.z}, dv );
        double dh = pos.y - hground - h;
        if( dh  < 0.0 ){
            //w->force.add( {0.0,gravity,0.0} );
            //w->force.add( {dv.x, dh*(-1-0.5*w->vel.y), dv.y} );
            force.add( { -dv.x, dh*(-100+2.5*vel.y), -dv.y } );
            force.add( {vel.x*terrain->drag,0.0,vel.z*terrain->drag} );
            //printf( " dv (%3.3f,%3.3f) (%3.3f,%3.3f)  \n", dv.x, dv.y, w->pos.x,w->pos.y );
            //w->force.add( {0, dh*(-100+2.5*w->vel.y), 0} );

        }
    }

};

#endif
