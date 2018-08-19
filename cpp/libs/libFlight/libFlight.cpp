#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cstdint>

/*
#include "arrayAlgs.h"
#include "AeroSurf.h"
#include "AeroCraft.h"
#include "geom3D.h"

#include "ShooterCommon.h"
#include "Projectile3D.h"
*/

#include "libFlight.h"

FlightWorld world;

// ==================== EXPORTED FUNCTIONS ===================

extern "C"{

void* getWorldPointer(){ return &world; };

void loadFromFile( char* fname ){


    int nint = 10
    //int nint = 1
    , nleft = nint;

    //int arr[nint] = {0,1,2,3,-4,5,-6,7,-8,9,10,11,-12,13};
    //int arr[nint] = {-1,-2,-3,-4,-5,-6,-7,-8,-9,-10};
    int arr[nint] = {1,2,-3,-4,5,6,7,8,-9,-10};
    //int arr[nint] = {-1};

    nleft = prune( nint, arr , [](int i){  return i<0; } );
    printf( "nint %i nleft %i \n", nint, nleft );
    for(int i=0; i<nint; i++){ if(i==nleft)printf(" | "); printf( " %i", arr[i] ); } printf("\n");
    //exit(0);


    world.craft1   .fromFile(fname);
    world.craft_bak.fromFile(fname);

    world.shotType.mass    = 42e-3;   // kg
    world.shotType.caliber = 12.7e-3; // m
    world.shotType.updateAux(0.33);
    world.shotType.balisticCoef = 0.0;
    printf( "projectile ballistic coef: %g ", world.shotType.balisticCoef );


    world.controlsState[iPitch    ].setLimits( 0.5, -0.2,0.2 );
    world.controlsState[iYaw      ].setLimits( 0.5, -0.2,0.2 );
    world.controlsState[iRoll     ].setLimits( 0.5, -0.2,0.2 );
    world.controlsState[iThrottle ].setLimits( 0.5,  0.0,1.0 );
    for(int i=0; i<world.nControls; i++) world.controlsState[i].x = 0.0;
    world.controlsState[iThrottle ].x = 1.0;

}

void setPose( double* pos, double* vel, double* rot ){
    world.craft1.pos    = *(Vec3d*)pos;
    world.craft1.vel    = *(Vec3d*)vel;
    world.craft1.rotMat = *(Mat3d*)rot;
    //printf( "pos (%g,%g,%g) \n", craft1.pos.x, craft1.pos.y, craft1.pos.z );
    //printf( "vel (%g,%g,%g) \n", craft1.vel.x, craft1.vel.y, craft1.vel.z );
    //printf( "rot\n"); craft1.rotMat.print();
}

void setTilt( int iwing, double tilt ){
    //world.craft1.panels[iwing].lrot.rotate( tilt, world->craft1.panels[iwing].lrot.a );
    world.setWing( iwing, tilt );
}

void setTargets( int n, double* targets_){
    world.nTargets = n;
    world.targets  = (Sphere3d*)targets_;
}

void setWing( int i, double angle ){ world.setWing( i, angle ); };
void fly    ( int n, int nsub, double dt, Vec3d* pos, Vec3d* vel, Mat3d* rot ){ world.fly( n, nsub, dt, (Vec3d*)pos, (Vec3d*)vel, (Mat3d*)rot ); };
void flyAndShootTargets( int nsub, double dt, double* controlBuff, double* stateBuff, double* targetShot ){ world.flyAndShootTargets( nsub, dt, controlBuff, stateBuff, targetShot ); };

} // extern "C"
