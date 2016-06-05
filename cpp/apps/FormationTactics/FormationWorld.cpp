
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

    nSoldiers = 0;
    for( Formation* f : formations ){
       // if( f != NULL ){
            f->eliminateInvalids( );
            f->clean_temp();
            //f->applyWillForce( );
            //f->interactInside( ); // this is now in formationInteractions
            f->update_bbox( );
            nSoldiers += f->nCapable;
       // }
    }

    nSoldierInteractions = formationInteractions( );
    //nSoldierInteractions = formationInteractions_buff( );

    for( Formation* f : formations ){
        //if( f != NULL ){
            //f->moveBy( {0.01, 0.01 } );
            f->update( dt );
        //}
    }
};

int FormationWorld::formationInteractions( ){
    int nInteractions = 0;
    for( int i=0; i<formations.size(); i++ ){
    Formation * fi = formations[i];
    //if( fi != NULL ){
        //for( int j=0; j<i; j++ ){ // This will be more complicated to resolve symetrically
        for( int j=0; j<formations.size(); j++ ){
            if( i==j ){
                nInteractions += fi->interactInside( );
            }else{
                nInteractions += fi->interact( formations[j] );
            }
        }
    //}
    }
    return nInteractions;
}

int FormationWorld::formationInteractions_buff( ){ // this is version of interactions accelerated by collision buffer
    int nInteractions = 0;
    for( int i=0; i<formations.size(); i++ ){
        Formation * fi = formations[i];
        double r = RmaxInteract; // FICME - this should be specific for formation pair
        colruler.setup(  {fi->bbox.x0-3*r,fi->bbox.y0-3*r}, {2*r+0.1,2*r+0.1} );
        //colruler.setup(  {fi->bbox.x0-3*r,fi->bbox.y0-3*r}, {3*r,3*r} );
        double rf = r/colruler.step.x;
        //printf( " (%3.3f,%3.3f) (%3.3f,%3.3f) (%3.3f,%3.3f)\n",  fi->bbox.x0, fi->bbox.y0, fi->bbox.x1, fi->bbox.y1, colruler.step.x, colruler.step.y );
        //Draw2D::drawGrid(  fi->bbox.x0, fi->bbox.y0, fi->bbox.x1, fi->bbox.y1, colruler.step.x, colruler.step.y  );
        colbuf.clear( );

        for( int k=0; k<fi->nCapable; k++){ // put soldiers to collision buffer
            Soldier * s = fi->soldiers + k;
            Vec2d dipos; Vec2i ipos;
            colruler.pos2index( s->pos, dipos, ipos );
            colbuf.insert( s, ipos, dipos, rf );
            //int h = colbuf.xy2i( ipos.x, ipos.y );
            //Draw2D::color_of_hash( h+1000 );
            //Draw2D::drawCircle_d( s->pos, 0.5, 8, true );
        }

/*
        glBegin(GL_POINTS);
        for( int k=0; k<1000; k++){ // put soldiers to collision buffer
            Vec2d pos;
            pos.set( randf(fi->bbox.x0, fi->bbox.x1), randf(fi->bbox.y0, fi->bbox.y1) );
            Vec2d dipos; Vec2i ipos;
            colruler.pos2index( pos, dipos, ipos );
            //colbuf.insert( s, ipos, dipos, r );
            int h = colbuf.xy2i( ipos.x, ipos.y );
            Draw2D::color_of_hash( h+1000 );
            glVertex3f( (float)pos.x, (float)pos.y, 0.0 );
        }
        glEnd();
        glColor3f(0.0f,1.0f,0.0f);
        //break;
*/

/*
        // just debug
        printf( " ==== i,k %i %i  \n", i, colbuf.count );
        for( int icell=0; icell<colbuf.count; icell++ ){
            int ix,iy;
            colbuf.i2xy( icell, ix, iy );
            if( colbuf.counts[icell] > 0 ) printf( " (%i,%i,%i) ", ix,iy, colbuf.counts[icell] );
            // printf( " >>> %i (%i,%i) %i \n", i,ix,ix, colbuf.counts[icell] );
            //for( int im=0; im<colbuf.counts[icell]; im++ ){
            //    Soldier * si = colbuf.get( icell, im );
            //    printf( " %i (%3.3f,%3.3f) \n", im, si->pos.x,si->pos.y );
            //}
        }
        printf( "\n" );
*/
        // PROBLEM : may interact several-times with the same neighbor

        for( int j=0; j<formations.size(); j++ ){ // cast other formations against that buffer
            if(i==j){
                for( int kj=0; kj<fi->nCapable; kj++ ){
                    Soldier * sj = fi->soldiers + kj;
                    int ix  = colruler.x2i( sj->pos.x );
                    int iy  = colruler.y2i( sj->pos.y );
                    int ixy = colbuf.xy2i(ix,iy);
                    for( int im=0; im<colbuf.counts[ixy]; im++ ){
                        Soldier * si = colbuf.get( ixy, im );
                        //if( (si == NULL) ){ printf( "error ixy %i im %i is NULL \n", ixy, im ); exit(0); }
                        nInteractions++;
                        //if ( si < sj ) si->friend_interaction( sj );
                        if( si != sj ) sj->friend_interaction( si );
                        //if( si != sj ){  si->friend_interaction( sj ); Draw2D::drawLine_d ( si->pos, sj->pos );  }
                        //if( si != sj ){  if( si->friend_interaction( sj ) ) Draw2D::drawLine_d ( si->pos, sj->pos );  }
                    }
                }
            }else{
                Formation * fj = formations[j];
                if ( fi->bbox.notOverlaps( fj->bbox ) ) continue;
                //printf( " interact formations %i %i \n", i, j );
                bool enemy = ( fi->faction != fj->faction );
                for( int kj=0; kj<fj->nCapable; kj++ ){
                    Soldier * sj = fj->soldiers + kj;
                    int ix  = colruler.x2i( sj->pos.x );   if ( (ix<0) || (ix>=colbuf.NX) ) continue;
                    int iy  = colruler.y2i( sj->pos.y );   if ( (iy<0) || (iy>=colbuf.NY) ) continue;
                    //if( (ix<0) || (ix>=colbuf.NX) ){ printf( "error ix %i is out of range \n", ix, 0, colbuf.NX ); return; }
                    //if( (iy<0) || (iy>=colbuf.NY) ){ printf( "error ix %i is out of range \n", iy, 0, colbuf.NY ); return; }
                    int ixy = colbuf.xy2i(ix,iy);
                    for( int im=0; im<colbuf.counts[ixy]; im++ ){
                        //Soldier * si = colbuf.buff[ colbuf.im2i( ixy, im ) ];
                        Soldier * si = colbuf.get( ixy, im );
                        if( (si == NULL) ){ printf( "error ixy %i im %i is NULL \n", ixy, im ); exit(0); }
                        nInteractions++;
                        if( enemy ) { si->enemy_interaction ( sj, fi->melee );  }
                        else        { si->friend_interaction( sj );             }
                    }
                }
            }
        }
        //exit(0);
        //return nInteractions;
    }
    return nInteractions;
}


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

    terrain.init( 100, 100, 30.0 );
    terrain.x0 = -0.5 * terrain.nx * terrain.step;
    terrain.y0 = -0.5 * terrain.ny * terrain.step;
    terrain.allocate( );
    terrain.generateRandom( 0.0, 1.0 );

    //soldierTypes.push_back( SoldierType(){"pikemen",1.0d,0.25d,1.0d} );
    //soldierTypes.push_back( {"pikemen",1.0d,0.25d,1.0d, 1.0, 1.0 } );
    soldierTypes.push_back( SoldierType() );

    Faction* fac1 = new Faction( "RedArmy" , {1.0f,0.25f,0.0f} );
    factions.push_back( fac1 );
    fac1->initFaction( 4, 4, 16, soldierTypes, {+20.0,3.0}, {-20.0,3.0}, 1.0 );
    //fac1->battleLines[0]->setTargetLine( {-10.0,3.0}, {+10.0,3.0} );

    Faction* fac2 = new Faction( "BlueArmy", {0.0f,0.5f, 1.0f} );
    factions.push_back( fac2 );
    fac2->initFaction( 4, 4, 16, soldierTypes, {-20.0,-3.0}, {+20.0,-3.0}, 1.0 );
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



