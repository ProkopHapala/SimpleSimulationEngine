
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"

#include "FormationWorld.h" // THE HEADER

#define DEBUG_PLOT_INTERACTION( pa, pb, R, G, B ) if(interacts){ \
    glColor3f( R, G, B ); \
    Draw2D::drawLine_d( pa->pos, pb->pos ); };

void FormationWorld::update( ){
    for( int i=0; i<per_frame; i++  ){
        //n_moves=0; n_interactions=0;
        simulationStep( dt );
        //simulationStep_semiBruteForce( dt );
        //simulationStep_BruteForce( dt );
        //printf( " ==== DONE sub_step %i  v2max %f f2max %f n_moves %i n_inter %i \n", i, v2max, f2max, n_moves, n_interactions );
    }
}

void FormationWorld::simulationStep( double dt ){

    for( Formation* f : formations ){
       // if( f != NULL ){
            f->clean_temp();
            f->applyWillForce( );
            f->interactInside( );
            f->update_bbox( );
       // }
    }

    for( int i=0; i<formations.size(); i++ ){
        Formation * fi = formations[i];
        //if( fi != NULL ){
            for( int j=0; j<i; j++ ){
                fi->interact( formations[j] );
            }
        //}
    }

    for( Formation* f : formations ){
        //if( f != NULL ){
            //f->moveBy( {0.01, 0.01 } );
            f->update( dt );
        //}
    }
};


void FormationWorld::refreshFormations( ){
    formations.clear();
    for( Faction* fa : factions ){
        for( Formation* fm : fa->formations ){
            formations.push_back( fm );
        }
    }
}

void FormationWorld::init(){
    printf( " FormationWorld::init() \n" );
    evalAuxSimParams();

    terrain.init( 100, 100, 5.0 );
    terrain.x0 = -0.5 * terrain.nx * terrain.step;
    terrain.y0 = -0.5 * terrain.ny * terrain.step;
    terrain.allocate( );
    terrain.generateRandom( 0.0, 1.0 );

    //soldierTypes.push_back( SoldierType(){"pikemen",1.0d,0.25d,1.0d} );
    soldierTypes.push_back( {"pikemen",1.0d,0.25d,1.0d} );

    Faction* fac1 = new Faction( "RedArmy" , {1.0f,0.25f,0.0f} );
    factions.push_back( fac1 );
    fac1->initFaction( 4, 4, 16, soldierTypes, {-20.0,3.0}, {+20.0,3.0}, 1.0 );
    //fac1->battleLines[0]->setTargetLine( {-10.0,3.0}, {+10.0,3.0} );

    Faction* fac2 = new Faction( "BlueArmy", {0.0f,0.5f, 1.0f} );
    factions.push_back( fac2 );
    fac2->initFaction( 4, 4, 16, soldierTypes, {+20.0,-3.0}, {-20.0,-3.0}, 1.0 );
    //fac2->battleLines[0]->setTargetLine( {+10.0,-3.0}, {-10.0,-3.0} );


    refreshFormations( );

/*
    formations.reserve( 16 );
    SoldierType * pikemen  = new SoldierType();
    Formation * formation1 = new Formation( 4, 4, pikemen );
    formation1->setEnds( {-2.0,-1.0}, {3.0,2.0}, 2.0 );
    formation1->deploySoldiers();
    formations.push_back( formation1 );
*/

};


/*

void NBodyWorld::assembleForces( ULONG i ){
    // BE WARE : particle->force should be cleaned before we start
    // onside step
    Particle2D*  buf_i_[256];
    Particle2D** buf_i = &buf_i_[0];
    UINT ni = map.HashMap<Particle2D>::getBucketObjects( i, buf_i );
    for(int ii=0; ii<ni; ii++ ){
        Particle2D* pi = buf_i[ii];
        for(int jj=0; jj<ii; jj++ ){
            Particle2D* pj = buf_i[jj];
            Vec2d fout;
            double qq = pi->charge * pj->charge;
            bool interacts = pairwiseForce( pi->pos, pj->pos, qq, fout );
            pi->force.add( fout );
            pj->force.sub( fout );
            n_interactions++;
            DEBUG_PLOT_INTERACTION( pi, pj, 0.1f, 0.9f, 0.1f )
        }
    }
    // offside part
    UHALF ix,iy;
    map.unfoldBucketInt( i, ix, iy );
    assembleForces_offside( i, map.getBucketInt( ix-1, iy-1 ), ni, buf_i );
    assembleForces_offside( i, map.getBucketInt( ix  , iy-1 ), ni, buf_i );
    assembleForces_offside( i, map.getBucketInt( ix+1, iy-1 ), ni, buf_i );
    assembleForces_offside( i, map.getBucketInt( ix-1, iy   ), ni, buf_i );
    //         onside                          ix   iy
    assembleForces_offside( i, map.getBucketInt( ix+1, iy   ), ni, buf_i );
    assembleForces_offside( i, map.getBucketInt( ix-1, iy+1 ), ni, buf_i );
    assembleForces_offside( i, map.getBucketInt( ix  , iy+1 ), ni, buf_i );
    assembleForces_offside( i, map.getBucketInt( ix+1, iy+1 ), ni, buf_i );
};

void NBodyWorld::assembleForces_offside( ULONG i, ULONG j, UINT ni, Particle2D** buf_i ){
    // BE WARE : particle->force should be cleaned before we start
    //printf( " assembleForces_offside === %i %i %i \n", i, j, ni );
    //if( activeCellsNeighbors.find(j) != activeCellsNeighbors.end() ){
    //if( i < j ){ // this will ensure that we do not double-count // WARRNIG : We miss situation when i>j and j is not active cell
        Particle2D* buf_j[256];
        UINT nj = map.HashMap<Particle2D>::getBucketObjects( j, buf_j );
        for(int ii=0; ii<ni; ii++ ){
            Particle2D* pi = buf_i[ii];
            for(int jj=0; jj<nj; jj++ ){
                Particle2D* pj = buf_j[jj];
                Vec2d fout;
                double qq = pi->charge * pj->charge;
                bool interacts = pairwiseForce( pi->pos, pj->pos, qq, fout );
                fout.mul(0.5d);      // if double counting ( not i<j condition )
                pi->force.add( fout );
                pj->force.sub( fout );
                //printf( " %i %i   %i %i  (%3.3f,%3.3f)(%3.3f,%3.3f)\n",   i, j,  ii, jj,  pi->pos.x,pi->pos.y,  pj->pos.x,pj->pos.y );
                n_interactions++;
                DEBUG_PLOT_INTERACTION( pi, pj, 0.9f, 0.1f, 0.9f )
            }
        }
    //}
};


*/



