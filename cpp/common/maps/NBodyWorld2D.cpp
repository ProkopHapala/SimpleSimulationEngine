
#include "NBodyWorld2D.h" // THE HEADER


/*
    TODO:
    there is some problem with reinserting :
    it is probably appearing when hash==0

 output obtained as:
  >> ./NBodyWorld > OUT
  >> grep 020-th OUT
 insering 020-th particle  (-8.517,-4.435) to (32762,32764) 2147254266
 reinsert: 020-th 2147254266=(32762,32764) -> 2147319802=(32762,32765)  121
 reinsert: 020-th 2147319802=(32762,32765) -> 2147254266=(32762,32764)  199
 reinsert: 020-th 2147254266=(32762,32764) -> 2147254267=(32763,32764)  0
 !!! cannot remove !!! : 020-th 2147254267=(32763,32764)

*/


/* Algorithm Pseudocode:
// force assembling step
clean set "ActiveParticles"
for all cells "i" from "ActiveCells" set:
    for all neighboring cells "j":
        for all particles "a" form cell "i":
            for all particles "b" from cell "j":
                - evaluate force between particles "a","b" and accumulate it in "a" and "b"
                - store the particles into temporary set "ActiveParticles"
// movement step
clean set "ActiveCells"
for all particles "a" from "ActiveParticles" set
    if force > forceCutoff:
        move "a"
        put cell "i" of particle "a" into set "ActiveCells"
*/

//#include <SDL2/SDL.h>
//#include <SDL2/SDL_opengl.h>
//#include "Draw2D.h"
//#define DEBUG_PLOT_INTERACTION( pa, pb, R, G, B ) if(interacts){ \
//    glColor3f( R, G, B ); \
//    Draw2D::drawLine_d( pa->pos, pb->pos ); };

void NBodyWorld::update( ){
    n_moves_frame=0; n_interactions_frame=0;
    for( int i=0; i<per_frame; i++  ){
        n_moves=0; n_interactions=0;
        simulationStep( dt );
        //simulationStep_semiBruteForce( dt );
        //simulationStep_BruteForce( dt );
       // printf( " ==== DONE sub_step %i  v2max %f f2max %f n_moves %i n_inter %i \n", i, v2max, f2max, n_moves, n_interactions );
       n_moves_frame        += n_moves;
       n_interactions_frame += n_interactions;
    }
}

void NBodyWorld::init( int power, double map_step, int nParticles_, Particle2D* particles_ ){

    evalAuxSimParams();

    nParticles = nParticles_;
    particles  = particles_;
    activeParticles = new Particle2D*[ nParticles ];

    //int power = 16; int nside = 300;
    //nParticles = (2*nside+1)*(2*nside+1);
    //nParticles = 4*nside*nside;

    //map.init( 4.0d, power );

    map.init( map_step, power );
	printf( "map: %i %i %i %i \n", map.power, map.mask, map.capacity, map.filled );

    for( int i=0; i<nParticles; i++ ){
        //map.insertNoTest( &(points[i]), points[i].x, points[i].y  );
        map.insertIfNew( &(particles[i]), particles[i].pos.x, particles[i].pos.y  );
        printf( " insering %03i-th particle  (%3.3f,%3.3f) to (%i,%i) %i \n", i,  particles[i].pos.x,  particles[i].pos.y, map.getIx( particles[i].pos.x), map.getIy( particles[i].pos.y), map.getBucket( particles[i].pos.x, particles[i].pos.y) );
    }

    //nParticles = makeParticleGrid( Particle2D * particles,  nside, nside, 1.0, 1.0, 0.4, 0.4 );

	printf( "map: %i %i %i %i \n", map.power, map.mask, map.capacity, map.filled );

    ULONG icell_old = 0;
    for( int i=0; i<nParticles; i++ ){ activateAroundParticle( &particles[i], icell_old ); };
    printf( "number of active cells %i\n", activeCells.size() );

};

bool NBodyWorld::interact( Particle2D * pi, Particle2D * pj ){
    Vec2d fout;
    double qq = pi->charge * pj->charge;
    bool interacts = pairwiseForce( pi->pos, pj->pos, qq, fout );
    //printf( "praticles  %i %i force %f %f \n", i, j, fout.x, fout.y );
    pi->force.add( fout );
    pj->force.sub( fout );
    return interacts;
};

// ==========================
//    Algorithm subroutines
// ==========================

void NBodyWorld::simulationStep( double dt ){
    // find particles in active cells, clean forces and put the particles to activeParticles buffer
    activeCellsNeighbors.clear();
    nActiveParticles = 0;
    for( ULONG icell : activeCells ){
        UHALF ix,iy;
        map.unfoldBucketInt( icell, ix, iy );
        activateCell( map.getBucketInt( ix-1, iy-1 ) );
        activateCell( map.getBucketInt( ix  , iy-1 ) );
        activateCell( map.getBucketInt( ix+1, iy-1 ) );
        activateCell( map.getBucketInt( ix-1, iy   ) );
        activateCell( icell                          );
        activateCell( map.getBucketInt( ix+1, iy   ) );
        activateCell( map.getBucketInt( ix-1, iy+1 ) );
        activateCell( map.getBucketInt( ix  , iy+1 ) );
        activateCell( map.getBucketInt( ix+1, iy+1 ) );
    }
    // evaluate pairwise forces
    //for( ULONG icell : activeCells ){ assembleForces( icell ); } // performance intensive step
    for( ULONG icell : activeCellsNeighbors ){ assembleForces( icell ); }

    if( picked != NULL ){
        Vec2d fstring;
        stringForce( picked->pos, anchor, anchorStiffness, fstring );
        picked->force.add( fstring );
    }

    // move the active particles and update activeCells accordingly
    // CONSIDERATION :  would be better to iterate over activeParticles or over activeCellsNighbors ?
    bool picked_done = false;
    activeCells.clear();
    ULONG icell_old = 0;
    v2max=0; f2max=0;
    for( int i=0; i<nActiveParticles; i++ ){
        Particle2D* pi = activeParticles[i];
        if( pi == picked ) picked_done = true;

        bool active = moveParticle( pi );
        if( active )  activateAroundParticle( pi, icell_old );
    }

    if( ( !picked_done ) && ( picked != NULL ) ){
        moveParticle( picked );
        activateAroundParticle( picked, icell_old );
    }

};

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
            bool interacts = interact( pi, pj );
            /*
            Vec2d fout;
            double qq = pi->charge * pj->charge;
            bool interacts = pairwiseForce( pi->pos, pj->pos, qq, fout );
            pi->force.add( fout );
            pj->force.sub( fout );
            */
            n_interactions++;
            //DEBUG_PLOT_INTERACTION( pi, pj, 0.1f, 0.9f, 0.1f )
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
                bool interacts = interact( pi, pj );
                /*
                Vec2d fout;
                double qq = pi->charge * pj->charge;
                bool interacts = pairwiseForce( pi->pos, pj->pos, qq, fout );
                fout.mul(0.5d);      // if double counting ( not i<j condition )
                pi->force.add( fout );
                pj->force.sub( fout );
                */
                //printf( " %i %i   %i %i  (%3.3f,%3.3f)(%3.3f,%3.3f)\n",   i, j,  ii, jj,  pi->pos.x,pi->pos.y,  pj->pos.x,pj->pos.y );
                n_interactions++;
                //DEBUG_PLOT_INTERACTION( pi, pj, 0.9f, 0.1f, 0.9f )
            }
        }
    //}
};

bool NBodyWorld::moveParticle( Particle2D* pi ){

    double v2 = pi->vel  .norm2(); v2max = (v2>v2max) ? v2 : v2max;
    double f2 = pi->force.norm2(); f2max = (f2>f2max) ? f2 : f2max;
    //bool forceOff = ( pi->force.norm2() < f2conv );
    //bool velOff   = ( pi->vel.norm2()   < v2conv );


    // Freezing condition
    //if( v2 < v2conv ) pi->vel  .set( 0.0d, 0.0d );
    //if( f2 < f2conv ) pi->force.set( 0.0d, 0.0d );
    if( f2 > f2conv ) pi->stepsConverged = 0;
    if( v2 < v2conv ){
        if( pi->stepsConverged > 10 ) return false;
        pi->stepsConverged++;
    };

    pi->vel.mul( damp );
    ULONG old_index = map.getBucket( pi->pos.x, pi->pos.y );
    pi->move_PointBody2D( dt );
    ULONG new_index = map.getBucket( pi->pos.x, pi->pos.y );
    n_moves++;

    //double v2 = pi->vel.  norm2();   v2max = (v2>v2max) ? v2 : v2max;
    //double f2 = pi->force.norm2();   f2max = (f2>f2max) ? f2 : f2max;

    if( old_index != new_index ){
        bool removed = map.HashMap<Particle2D>::tryRemove  ( pi, old_index );
        if( removed ){
            int iinsert = map.HashMap<Particle2D>::insertIfNew( pi, new_index );
            if( iinsert < 0 ){
                printf( " cannot insert! inconsistent HashMap! \n" );
                exit(0);
            }
        }else{
            printf( " cannot remove! inconsistent HashMap! \n" );
            exit(0);
        }
    }

    return true;
}

void NBodyWorld::activateCell( ULONG i ){
    // check if this is first time we visit this cell
    if( activeCellsNeighbors.find(i) != activeCellsNeighbors.end() ) return;
    // find all particles, put them to ActiveParticles, and clean the forces
    Particle2D* buf_i[256];
    UINT ni = map.HashMap<Particle2D>::getBucketObjects( i, &buf_i[0] );
    for(int ii=0; ii<ni; ii++ ){
        Particle2D* pi = buf_i[ii];
        pi->force.set( 0.0, 0.0 );
        activeParticles[nActiveParticles] = pi;
        nActiveParticles++;
    }
    // store cell into the visited cells
    activeCellsNeighbors.insert( i );
}

void NBodyWorld::activateAroundParticle( Particle2D* pi, ULONG& icell_old ){
    // CONSIDERATION : we can optimize here ... in activeParticles are particles from one cell grouped => we can check if icell changed from previous
    ULONG icell = map.getBucket( pi->pos.x, pi->pos.y );
    if( icell != icell_old ){
        if( activeCells.find(icell) == activeCells.end() ){
            activeCells.insert(icell);


            UHALF ix,iy,ix_,iy_;
            map.unfoldBucketInt( icell, ix_, iy_ );
            ix = map.getIx( pi->pos.x );
            iy = map.getIy( pi->pos.y );
            //printf( "activate cell %i=(%i,%i)=(%i,%i) pi (%3.3f,%3.3f) %i-th  \n", icell, ix_, iy_, ix,iy, pi->pos.x, pi->pos.y, pi-particles );
        }
        icell_old = icell;
    }
}

// ==========================
//    Just for debuging
// ==========================

void NBodyWorld::moveParticleDebug( Particle2D* pi, int i ){
    ULONG old_index = map.getBucket( pi->pos.x, pi->pos.y );
    pi->vel.mul( damp );
    pi->move_PointBody2D( dt );
    ULONG new_index = map.getBucket( pi->pos.x, pi->pos.y );
    // reinsert if particle moved to other cell
    if( old_index != new_index ){
        bool removed = map.HashMap<Particle2D>::tryRemove  ( pi, old_index );
        if( removed ){
            UHALF ixn,iyn, ixo,iyo;
            int iinsert = map.HashMap<Particle2D>::insertIfNew( pi, new_index );
            map.unfoldBucketInt( old_index, ixo, iyo );
            map.unfoldBucketInt( new_index, ixn, iyn );
            if( iinsert >= 0 ){
                printf( " reinsert: %03i-th %i=(%i,%i) -> %i=(%i,%i)  %i\n", i, old_index, ixo,iyo, new_index, ixn,iyn, iinsert );
            }else{
                printf( "!!! cannot insert !!! : %03i-th %i=(%i,%i) -> %i=(%i,%i) %i\n", i, old_index, ixo,iyo, new_index, ixn,iyn, iinsert  );
                exit(0);
            }
            //printf( " map after reinsert %i -> %i   filled %i capacity %i\n", old_index, new_index, map.filled, map.capacity );
        }else{
            UHALF ixn,iyn, ixo,iyo;
            map.unfoldBucketInt( old_index, ixo, iyo );
            //printf( " cannot remove! %03i-th (%3.3f,%3.3f) (%i,%i) bucket %i (%i,%i) \n", i, ox,oy, oix, oiy, old_index, ix, iy );
            printf( "!!! cannot remove !!! : %03i-th %i=(%i,%i) \n", i, old_index, ixo,iyo );
            printf( " map.fields[134]:  %i %i %i \n", map.fields[0].bucket, map.fields[0].object, pi );
            exit(0);
        }
    }
}

void NBodyWorld::simulationStep_BruteForce( double dt ){

    for( int i=0; i<nParticles; i++ ){ particles[i].force.set( 0.0, 0.0 ); }
    for( int i=0; i<nParticles; i++ ){
        Particle2D* pi = particles+i;
        for(int j=0; j<i; j++ ){
            Particle2D* pj = particles+j;
            bool interacts = interact( pi, pj );
            /*
            Vec2d fout;
            double qq = pi->charge * pj->charge;
            bool interacts = pairwiseForce( pi->pos, pj->pos, qq, fout );
            //printf( "praticles  %i %i force %f %f \n", i, j, fout.x, fout.y );
            pi->force.add( fout );
            pj->force.sub( fout );
            */
            n_interactions++;
            //DEBUG_PLOT_INTERACTION( pi, pj, 0.1f, 0.9f, 0.1f )
        }
    }

    if( picked != NULL ){
        Vec2d fstring;
        stringForce( picked->pos, anchor, anchorStiffness, fstring );
        picked->force.add( fstring );
    }

    activeCells.clear();
    ULONG icell_old = 0;
    v2max=0; f2max=0;
    for( int i=0; i<nParticles; i++ ){
        Particle2D* pi = particles+i;
        bool active = moveParticle( pi );
        //moveParticleDebug( pi, i );
        if( active ) activateAroundParticle( pi, icell_old );
    }

};


void NBodyWorld::simulationStep_semiBruteForce( double dt ){
    for( int i=0; i<nParticles; i++ ){
        particles[i].force.set( 0.0, 0.0 );
    }


    for( ULONG icell : activeCells ){ assembleForces( icell ); }


/*
    for( int ix=-10; ix<10; ix++ ){
        for( int iy=-10; iy<10; iy++ ){
            ULONG icell = map.getBucketInt( ix+MAP_OFFSET, iy+MAP_OFFSET );
            assembleForces( icell );
        }
    }
*/

    //exit(0);

    if( picked != NULL ){
        Vec2d fstring;
        stringForce( picked->pos, anchor, anchorStiffness, fstring );
        picked->force.add( fstring );
    }

    activeCells.clear();
    printf( "activeCells size : %i \n", activeCells.size() );
    ULONG icell_old = 0;
    v2max=0; f2max=0;
    for( int i=0; i<nParticles; i++ ){
        Particle2D* pi = particles+i;
        bool active = moveParticle( pi );
        //moveParticleDebug( pi, i );
        if( active ) activateAroundParticle( pi, icell_old );
    }
    printf( "activeCells size - : %i \n", activeCells.size() );
    checkHashMapConsistency( );

};


void NBodyWorld::checkHashMapConsistency( ){
    Particle2D* temp[256];
    for( int i=0; i<nParticles; i++ ){
        Particle2D* pi = particles+i;
        ULONG icell = map.getBucket( pi->pos.x, pi->pos.y );
        UINT nfound = map.HashMap<Particle2D>::getBucketObjects( icell, temp );
        bool found = false;
        for( int j=0; j<nfound; j++ ){
            if( temp[j] == pi ){ found = true; break; }
        }
        if( !found ){
            UHALF ix,iy;
            map.unfoldBucketInt( icell, ix, iy );
            printf( "!!! cannot find %i-th (%3.3f,%3.3f) particle in cell %i=(%i,%i)\n", i, pi->pos.x, pi->pos.y, icell, ix, iy );
            int ifield = map.HashMap<Particle2D>::findBruteForce( pi );
            if(ifield>=0){
                printf( "bruteForce found %i-th particle in field %i | object %i bucket %i n %i \n", i, ifield, map.fields[ifield].object, map.fields[ifield].bucket, map.fields[ifield].n  );
                //map.DEBUG_LEVEL = 3;
                int ifield_ = map.HashMap<Particle2D>::find( pi, icell );
                printf( " find returned %i \n", ifield_  );
                nfound = map.HashMap<Particle2D>::getBucketObjects( icell, temp );
                printf( " nfound %i \n", nfound );
                for( int j=0; j<nfound; j++ ){
                    printf( " temp[i] %i=(3.3f,3.3f) pi %i=(3.3f,3.3f) ", temp[j], temp[j]->pos.x, temp[j]->pos.y,   pi, pi->pos.x, pi->pos.y );
                }
            }
            UINT hash = map.mask&hashFunc( icell );
            int nhash = map.HashMap<Particle2D>::checkHashConsisent( hash );
            if( nhash == 0) printf( "hash %i (n=%i) checked fine \n", hash, map.fields[hash].n );
            exit(0);
        }
    }

}




