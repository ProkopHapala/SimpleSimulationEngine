#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cstdint>

#include "AeroSurf.h"
#include "AeroCraft.h"

AeroCraft craft1;

extern "C"{

void loadFromFile( char* fname ){
    craft1.fromFile(fname);
}

void setPose( double* pos, double* vel, double* rot ){
    craft1.pos    = *(Vec3d*)pos;
    craft1.vel    = *(Vec3d*)vel;
    craft1.rotMat = *(Mat3d*)rot;
    //printf( "pos (%g,%g,%g) \n", craft1.pos.x, craft1.pos.y, craft1.pos.z );
    //printf( "vel (%g,%g,%g) \n", craft1.vel.x, craft1.vel.y, craft1.vel.z );
    //printf( "rot\n"); craft1.rotMat.print();
}

void setTilt( int iwing, double tilt ){
    craft1.panels[iwing].lrot.rotate( tilt, craft1.panels[iwing].lrot.a );
}

void fly( int n, int nsub, double dt, double* pos_, double* vel_, double* rot_ ){
    Vec3d* pos=(Vec3d*)pos_;
    Vec3d* vel=(Vec3d*)vel_;
    Mat3d* rot=(Mat3d*)rot_;
    //printf( "n %i nsub %i ntot \n", n, nsub, n*nsub );
    for(int i=0; i<n; i++){
        //printf( "fly %i\n", i );
        for(int j=0; j<nsub; j++){
            //printf( "fly i %i j %i \n ", i, j );
            craft1.clean_temp();
            craft1.applyAeroForces( Vec3dZero );
            craft1.force.add( {0.0,-9.81,0.0} );
            craft1.move(dt);
            //printf( "fly i %i j %i pos (%g,%g,%g)  \n ", i, j, craft1.pos.x, craft1.pos.y, craft1.pos.z );
        }
        pos[i] = craft1.pos;
        vel[i] = craft1.vel;
        rot[i] = craft1.rotMat;
    }
}

}
