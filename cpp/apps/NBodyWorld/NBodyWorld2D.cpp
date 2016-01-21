
#include "NBodyWorld2D.h" // THE HEADER


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


void NBodyWorld::simulationStep( double dt ){
    // find particles in active cells, clean forces and put the particles to activeParticles buffer
    for( ULONG icell : activeCells ){
        UHALF ix,iy;
        map.unfoldBucket( icell, ix, iy );
        activateParticles( map.getBucket( ix-1, iy-1 ) );
        activateParticles( map.getBucket( ix  , iy-1 ) );
        activateParticles( map.getBucket( ix+1, iy-1 ) );
        activateParticles( map.getBucket( ix-1, iy   ) );
        activateParticles( icell                           );
        activateParticles( map.getBucket( ix+1, iy   ) );
        activateParticles( map.getBucket( ix-1, iy+1 ) );
        activateParticles( map.getBucket( ix  , iy+1 ) );
        activateParticles( map.getBucket( ix+1, iy+1 ) );
    }
    // evaluate pairwise forces
    for( ULONG icell : activeCells ){ assembleForces( icell ); } // performance intensive step
    // move the active particles and update activeCells accordingly
    // CONSIDERATION :  would be better to iterate over activeParticles or over activeCellsNighbors ?
    ULONG icell_old = 0;
    for( int i=0; i<nActiveParticles; i++ ){
        Particle2D* pi = activeParticles[i];
        pi->move_PointBody2D( dt );
        // CONSIDERATION : we can optimize here ... in activeParticles are particles from one cell grouped => we can check if icell changed from previous
        ULONG icell = map.getBucket( pi->pos.x, pi->pos.y );
        if( icell != icell_old ){
            if( activeCells.find(icell) == activeCells.end() ){ activeCells.insert(icell); }
            icell_old = icell;
        }
    }

};

void NBodyWorld::activateParticles( ULONG i ){
    // check if this is first time we visit this cell
    if( activeCellsNighbors.find(i) != activeCellsNighbors.end() ) return;
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
    activeCellsNighbors.insert( i );
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
            shortRangeForce( pi->pos, pj->pos, fout );
            pi->force.add( fout );
            pj->force.sub( fout );
        }
    }
    // offside part
    UHALF ix,iy;
    map.unfoldBucket( i, ix, iy );
    assembleForces_offside( i, map.getBucket( ix-1, iy-1 ), ni, &buf_i[0] );
    assembleForces_offside( i, map.getBucket( ix  , iy-1 ), ni, &buf_i[0] );
    assembleForces_offside( i, map.getBucket( ix+1, iy-1 ), ni, &buf_i[0] );
    assembleForces_offside( i, map.getBucket( ix-1, iy   ), ni, &buf_i[0] );
    //         onside                          ix   iy
    assembleForces_offside( i, map.getBucket( ix+1, iy   ), ni, &buf_i[0] );
    assembleForces_offside( i, map.getBucket( ix-1, iy+1 ), ni, &buf_i[0] );
    assembleForces_offside( i, map.getBucket( ix  , iy+1 ), ni, &buf_i[0] );
    assembleForces_offside( i, map.getBucket( ix+1, iy+1 ), ni, &buf_i[0] );
};

void NBodyWorld::assembleForces_offside( ULONG i, ULONG j, UINT ni, Particle2D** buf_i ){
    // BE WARE : particle->force should be cleaned before we start
    if( i < j ){ // this will ensure that we do not double-count
        Particle2D* buf_j[256];
        UINT nj = map.HashMap<Particle2D>::getBucketObjects( j, buf_j );
        for(int ii=0; ii<ni; ii++ ){
            Particle2D* pi = buf_i[ii];
            for(int jj=0; jj<nj; jj++ ){
                Particle2D* pj = buf_j[jj];
                Vec2d fout;
                shortRangeForce( pi->pos, pj->pos, fout );
                pi->force.add( fout );
                pj->force.sub( fout );
            }
        }
    }
};

void NBodyWorld::init(){

};


