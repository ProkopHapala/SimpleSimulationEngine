
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"

#include "Faction.h"
#include "TacWorld.h" // THE HEADER

#define DEBUG_PLOT_INTERACTION( pa, pb, R, G, B ) if(interacts){ \
    glColor3f( R, G, B ); \
    Draw2D::drawLine_d( pa->pos, pb->pos ); };

void TacWorld::update( ){
    for( int i=0; i<per_frame; i++  ){
        simulationStep( dt );
    }
}

void TacWorld::simulationStep( double dt ){

    for( Unit* f : units ){
        f->update( dt );
    }
};


int TacWorld::getUnitAt( const Vec2d& p, Faction * faction ){
    int imin=-1;
    double r2min = 1e+300;
    //for( Unit * u : units ){ // this does not seem to work
    for( int i=0; i<units.size(); i++ ){
        Unit * u = units[i];
        if( u->faction == faction ) continue;
        double r2 = p.dist2( u->pos );
        if( (r2<sq(u->radius))&&(r2 < r2min) ){ r2min=r2; imin=i; }
    }
    //printf( " imin %i r2min %f \n", imin, r2min );
    return imin;
};


void TacWorld::init(){
    printf( " TacWorld::init() \n" );
    evalAuxSimParams();

    terrain.init( 100, 100, 30.0 );
    terrain.x0 = -0.5 * terrain.nx * terrain.step;
    terrain.y0 = -0.5 * terrain.ny * terrain.step;
    terrain.allocate( );
    terrain.generateRandom( 0.0, 1.0 );

    //soldierTypes.push_back( SoldierType(){"pikemen",1.0d,0.25d,1.0d} );
    //soldierTypes.push_back( {"pikemen",1.0d,0.25d,1.0d, 1.0, 1.0 } );
    //unitTypes.push_back( UnitType() );

    Unit * u;

    Faction* fac1 = new Faction( "RedArmy" , {1.0f,0.25f,0.0f} );
    factions.push_back( fac1 );
    u = new Unit( &unitTypes[0], fac1, {-5.0,-3.0} ); fac1->units.push_back(u); units.push_back(u);
    u = new Unit( &unitTypes[0], fac1, { 0.0,-3.0} ); fac1->units.push_back(u); units.push_back(u);
    u = new Unit( &unitTypes[0], fac1, {+5.0,-3.0} ); fac1->units.push_back(u); units.push_back(u);

    Faction* fac2 = new Faction( "BlueArmy", {0.0f,0.5f, 1.0f} );
    factions.push_back( fac2 );
    u = new Unit( &unitTypes[0], fac2, {-5.0, 3.0} ); fac2->units.push_back(u); units.push_back(u);
    u = new Unit( &unitTypes[0], fac2, { 0.0, 3.0} ); fac2->units.push_back(u); units.push_back(u);
    u = new Unit( &unitTypes[0], fac2, {+5.0, 3.0} ); fac2->units.push_back(u); units.push_back(u);


};






