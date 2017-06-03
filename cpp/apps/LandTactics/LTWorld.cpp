
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"

#include "LTFaction.h"
#include "LTWorld.h" // THE HEADER

#define DEBUG_PLOT_INTERACTION( pa, pb, R, G, B ) if(interacts){ \
    glColor3f( R, G, B ); \
    Draw2D::drawLine_d( pa->pos, pb->pos ); };

void LTWorld::update( ){
    for( int i=0; i<per_frame; i++  ){
        simulationStep( dt );
    }
}

void LTWorld::simulationStep( double dt ){
    for( LTUnit* f : units ){
        f->update( dt );
    }
};


int LTWorld::getUnitAt( const Vec2d& p, LTFaction * faction ){
    int imin=-1;
    double r2min = 1e+300;
    //for( Unit * u : units ){ // this does not seem to work
    for( int i=0; i<units.size(); i++ ){
        LTUnit * u = units[i];
        if( u->faction == faction ) continue;
        double r2 = p.dist2( u->pos );
        if( (r2<sq(u->radius))&&(r2 < r2min) ){ r2min=r2; imin=i; }
    }
    //printf( " imin %i r2min %f \n", imin, r2min );
    return imin;
};


void LTWorld::init(){
    printf( " LTWorld::init() \n" );
    evalAuxSimParams();

    ruler.setSize(128,128);
    ruler.setStep(50);

    map_center = (Vec2d){ruler.na*0.75*ruler.step,ruler.nb*0.5*ruler.step};

    ground = new double[ruler.ntot];
    hydraulics.setSize(ruler.na,ruler.nb);
    hydraulics.ground = ground;

    //world.hydraulics.allocate( 512, 512 );

    hydraulics.genTerrainNoise( 8, 2.0, 1.0,  0.5, 0.8, 45454, {100.0,100.0} );
    for( int j=0; j<500; j++ ){
        int isz = 25;
        int ix0 = rand()%(hydraulics.nx-isz);
        int iy0 = rand()%(hydraulics.ny-isz);
        hydraulics.errodeDroples( 200, 100, 0.02, 0.15, 0.5, ix0, iy0, ix0+isz, iy0+isz );
    }
    for(int i=0; i<ruler.ntot; i++){ ground[i] *= maxHeight; };
    //for(int i=0; i<ruler.ntot; i++){ ground[i] = randf(0.0,1.0); };
    //for(int ib=0; ib<ruler.nb; ib++){  for(int ia=0; ia<ruler.na; ia++){  ground[ib*ruler.na+ia] = ia/(float)ruler.na;  } };

    //soldierTypes.push_back( SoldierType(){"pikemen",1.0d,0.25d,1.0d} );
    //soldierTypes.push_back( {"pikemen",1.0d,0.25d,1.0d, 1.0, 1.0 } );
    //unitTypes.push_back( UnitType() );

    LTUnit * u;

    LTFaction* fac1 = new LTFaction( "RedArmy" , {1.0f,0.25f,0.0f} );
    factions.push_back( fac1 );
    u = new LTUnit( &unitTypes[0], fac1, map_center+(Vec2d){-50.0,-30.0} ); fac1->units.push_back(u); units.push_back(u);
    u = new LTUnit( &unitTypes[0], fac1, map_center+(Vec2d){ 0.0,-30.0} ); fac1->units.push_back(u); units.push_back(u);
    u = new LTUnit( &unitTypes[0], fac1, map_center+(Vec2d){+50.0,-30.0} ); fac1->units.push_back(u); units.push_back(u);

    LTFaction* fac2 = new LTFaction( "BlueArmy", {0.0f,0.5f, 1.0f} );
    factions.push_back( fac2 );
    u = new LTUnit( &unitTypes[0], fac2, map_center+(Vec2d){-50.0, 30.0} ); fac2->units.push_back(u); units.push_back(u);
    u = new LTUnit( &unitTypes[0], fac2, map_center+(Vec2d){ 00.0, 30.0} ); fac2->units.push_back(u); units.push_back(u);
    u = new LTUnit( &unitTypes[0], fac2, map_center+(Vec2d){+50.0, 30.0} ); fac2->units.push_back(u); units.push_back(u);


};






