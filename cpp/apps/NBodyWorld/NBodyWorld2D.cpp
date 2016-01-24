
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

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "drawMath2D.h"

#define DEBUG_PLOT_INTERACTION( pa, pb, R, G, B ) \
    glColor3f( R, G, B ); \
    Draw2D::drawLine_d( pa->pos, pb->pos );

void NBodyWorld::update( ){
    for( int i=0; i<per_frame; i++  ){
        //simulationStep( dt );
        simulationStep_semiBruteForce( dt );
        //simulationStep_BruteForce( dt );
    }
}

void NBodyWorld::moveParticle( Particle2D* pi ){
    ULONG old_index = map.getBucket( pi->pos.x, pi->pos.y );
    pi->vel.mul( damp );
    pi->move_PointBody2D( dt );
    ULONG new_index = map.getBucket( pi->pos.x, pi->pos.y );
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
}

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
    for(int i=0; i<nParticles; i++ ){
        Particle2D* pi = particles+i;
        for(int j=0; j<i; j++ ){
            Particle2D* pj = particles+j;
            Vec2d fout;
            double qq = pi->charge * pj->charge;
            pairwiseForce( pi->pos, pj->pos, qq, fout );
            //printf( "praticles  %i %i force %f %f \n", i, j, fout.x, fout.y );
            pi->force.add( fout );
            pj->force.sub( fout );
        }
    }

    if( picked != NULL ){
        Vec2d fstring;
        stringForce( picked->pos, anchor, anchorStiffness, fstring );
        picked->force.add( fstring );
    }

    for( int i=0; i<nParticles; i++ ){
        Particle2D* pi = particles+i;
        moveParticle( pi );
        //moveParticleDebug( pi, i );
    }

};


void NBodyWorld::simulationStep_semiBruteForce( double dt ){
    activeCells.clear();
    for( int i=0; i<nParticles; i++ ){
        particles[i].force.set( 0.0, 0.0 );
        ULONG icell = map.getBucket( particles[i].pos.x, particles[i].pos.y );
        if( activeCells.find( icell ) == activeCells.end() ){
            activeCells.insert(icell);
        }
    }
    for( ULONG icell : activeCells ){ assembleForces( icell ); }
    //exit(0);

    if( picked != NULL ){
        Vec2d fstring;
        stringForce( picked->pos, anchor, anchorStiffness, fstring );
        picked->force.add( fstring );
    }

    for( int i=0; i<nParticles; i++ ){
        Particle2D* pi = particles+i;
        moveParticle( pi );
        //moveParticleDebug( pi, i );
    }

};


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
        activateCell( icell                       );
        activateCell( map.getBucketInt( ix+1, iy   ) );
        activateCell( map.getBucketInt( ix-1, iy+1 ) );
        activateCell( map.getBucketInt( ix  , iy+1 ) );
        activateCell( map.getBucketInt( ix+1, iy+1 ) );
    }
    // evaluate pairwise forces
    for( ULONG icell : activeCells ){ assembleForces( icell ); } // performance intensive step
    // move the active particles and update activeCells accordingly
    // CONSIDERATION :  would be better to iterate over activeParticles or over activeCellsNighbors ?
    activeCells.clear();
    ULONG icell_old = 0;
    for( int i=0; i<nActiveParticles; i++ ){
        Particle2D* pi = activeParticles[i];
        pi->move_PointBody2D( dt );
        activateAroundParticle( pi, icell_old );
    }
};

void NBodyWorld::activateCell( ULONG i ){
    // check if this is first time we visit this cell
    if( activeCellsNeighbors.find(i) != activeCellsNeighbors.end() ) return;
    // find all particles, put them to ActiveParticles, and clean the forces
    Particle2D* buf_i[256];
    UINT ni = map.HashMap<Particle2D>::getBucketObjects( i, &buf_i[0] );
    for(int ii=0; ii<ni; ii++ ){
        Particle2D* pi = buf_i[ii];
        pi->force.set( 0.0d, 0.0d );
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
        if( activeCells.find(icell) == activeCells.end() ){ activeCells.insert(icell); }
        icell_old = icell;
    }
}

void NBodyWorld::assembleForces( ULONG i ){
    // BE WARE : particle->force should be cleaned before we start
    // onside step
    Particle2D* buf_i[256];
    UINT ni = map.HashMap<Particle2D>::getBucketObjects( i, &buf_i[0] );
    for(int ii=0; ii<ni; ii++ ){
        Particle2D* pi = buf_i[ii];
        for(int jj=0; jj<ii; jj++ ){
            Particle2D* pj = buf_i[jj];
            Vec2d fout;
            double qq = pi->charge * pj->charge;
            pairwiseForce( pi->pos, pj->pos, qq, fout );
            pi->force.add( fout );
            pj->force.sub( fout );
            DEBUG_PLOT_INTERACTION( pi, pj, 0.1f, 0.9f, 0.1f )
        }
    }
    // offside part
    UHALF ix,iy;
    map.unfoldBucketInt( i, ix, iy );
    assembleForces_offside( i, map.getBucketInt( ix-1, iy-1 ), ni, &buf_i[0] );
    assembleForces_offside( i, map.getBucketInt( ix  , iy-1 ), ni, &buf_i[0] );
    assembleForces_offside( i, map.getBucketInt( ix+1, iy-1 ), ni, &buf_i[0] );
    assembleForces_offside( i, map.getBucketInt( ix-1, iy   ), ni, &buf_i[0] );
    //         onside                          ix   iy
    assembleForces_offside( i, map.getBucketInt( ix+1, iy   ), ni, &buf_i[0] );
    assembleForces_offside( i, map.getBucketInt( ix-1, iy+1 ), ni, &buf_i[0] );
    assembleForces_offside( i, map.getBucketInt( ix  , iy+1 ), ni, &buf_i[0] );
    assembleForces_offside( i, map.getBucketInt( ix+1, iy+1 ), ni, &buf_i[0] );
};

void NBodyWorld::assembleForces_offside( ULONG i, ULONG j, UINT ni, Particle2D** buf_i ){
    // BE WARE : particle->force should be cleaned before we start
    //printf( " assembleForces_offside === %i %i %i \n", i, j, ni );
    if( i < j ){ // this will ensure that we do not double-count
        Particle2D* buf_j[256];
        UINT nj = map.HashMap<Particle2D>::getBucketObjects( j, buf_j );
        for(int ii=0; ii<ni; ii++ ){
            Particle2D* pi = buf_i[ii];
            for(int jj=0; jj<nj; jj++ ){
                Particle2D* pj = buf_j[jj];
                Vec2d fout;
                double qq = pi->charge * pj->charge;
                pairwiseForce( pi->pos, pj->pos, qq, fout );
                pi->force.add( fout );
                pj->force.sub( fout );
                //printf( " %i %i   %i %i  (%3.3f,%3.3f)(%3.3f,%3.3f)\n",   i, j,  ii, jj,  pi->pos.x,pi->pos.y,  pj->pos.x,pj->pos.y );
                DEBUG_PLOT_INTERACTION( pi, pj, 0.9f, 0.1f, 0.9f )
            }
        }
    }
};

void NBodyWorld::init(){
    evalAuxSimParams();

    activeParticles = new Particle2D*[ 1<<20 ];

    int power = 8; int nside = 5;
    //int power = 11; int nside = 20;
    //int power = 16; int nside = 300;
    nParticles = (2*nside+1)*(2*nside+1);
    //nParticles = 4*nside*nside;
	particles  = new Particle2D[nParticles];
    map.init( 2.0f, power );
	printf( "map: %i %i %i %i \n", map.power, map.mask, map.capacity, map.filled );
	int i = 0;
	for( int iy=-nside; iy<nside; iy++ ){
		for( int ix=-nside; ix<nside; ix++ ){
			particles[i].charge = ((int)randf(0.0, 1.9999999))*2 -1;
			particles[i].vel.set( 0.0, 0.0 );
			particles[i].setMass( 1.0 );
			particles[i].pos.set(
                ( ix + randf(0.2,0.8) ) * map.step,
                ( iy + randf(0.2,0.8) ) * map.step );
                //( ix + 0.5 ) * map.step,
                //( iy + 0.5 ) * map.step );
			//map.insertNoTest( &(points[i]), points[i].x, points[i].y  );
			map.insertIfNew( &(particles[i]), particles[i].pos.x, particles[i].pos.y  );
			printf( " insering %03i-th particle  (%3.3f,%3.3f) to (%i,%i) %i \n", i,  particles[i].pos.x,  particles[i].pos.y, map.getIx( particles[i].pos.x), map.getIy( particles[i].pos.y), map.getBucket( particles[i].pos.x, particles[i].pos.y) );
            i++;
		};
	};
	nParticles = i;
	printf( "map: %i %i %i %i \n", map.power, map.mask, map.capacity, map.filled );


    ULONG icell_old = 0;
    for( int i=0; i<nParticles; i++ ){ activateAroundParticle( &particles[i], icell_old ); };
    printf( "number of active cells %i\n", activeCells.size() );


};


