
#include "MoleculeWorld2D.h" // THE HEADER

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"

#define DEBUG_PLOT_INTERACTION( pa, pb, R, G, B ) if(interacts){ \
    glColor3f( R, G, B ); \
    Draw2D::drawLine_d( pa->pos, pb->pos ); };




void interMolForce(
	int na, double * aRs, double * aQs, Vec2d * aPos, Vec2d * aFs,
	int nb, double * bRs, double * bQs, Vec2d * bPos, Vec2d * bFs
){
	//printf( "---------------------\n" );
	for(int ia=0; ia<na; ia++){
		for(int ib=0; ib<nb; ib++){
			Vec2d f;
			//dR.set_sub( aRs[ia], bRs[ib] );
			double qq = aQs[ia]*bQs[ib];
			double R  = aRs[ia]+bRs[ib];
			pairwiseForce( aPos[ia], bPos[ib], qq, R*R, f );
			aFs[ia].add( f ); bFs[ib].sub( f );
			//aFs[ia].sub( f ); bFs[ib].add( f );
			//printf( " %i %i   %f %f    %f %f %f \n", ia, ib, C6, C12,    f.x, f.y, f.z );
		}
	}
}

void forceTransform( int n, Vec2d * poss, Vec2d * fs,  Vec2d center,  Vec2d& f, double& tq ){
    for( int i = 0; i < n; i++ ){
        Vec2d d; d.set_sub( poss[i], center );
        tq += d.cross( fs[i] );
        f.add( fs[i] );
    }
}


/*
void MoleculeWorld::optStep( ){
// one step of the relaxation

	// initialize everything associated with forces to zero
	for( int i = 0; i < nmols; i++ ){
		fpos[i].set( 0.0 );			// set all "forces positions" to zero
		frot[i].set( 0.0, 0.0, 0.0, 0.0 );	// set all "forces rotations" to zero
		rot[i].normalize();			// keep quaternion normalized, otherwise unstable !!!
	}

	Vec3d vec;

	// calculate forces
	forcesMolecules();
	if( tip != NULL ){
		forcesTip( vec );
	}
	pixelDataListItem = vec.z;

	// in frot[i] left only its part which is perpendicular to rot[i]
	for( int i = 0; i < nmols; i++ ){
		double qfq = rot[i].dot( frot[i] );	// projection of frot[i] onto rot[i]
		frot[i].add_mul( rot[i], -qfq );	// from frot[i] subtract its part parallel to rot[i]
	}

	// update positions and velocities of molecules
	optimizer->move();
}
*/



void MoleculeWorld::forcesMolecules( ){
	for( int i = 0; i < nMols; i++ ){
        int itype = mols[i];
		MoleculeType2D* moli = &molTypes[itype];				                // in moli store i-th molecule
		int npi = moli->natoms;						                        // npi is a number of atoms in i-th molecule
		moli->rigidTransform( rot[i], pos[i], Tps_i );
		//cleanPointForce( npi, fs_i );					                    // initialize fs_i to zeros
		VecN::set( 2*npi, 0.0, (double*)fs_i );
		for( int j = 0; j < i; j++ ){
		    int jtype = mols[i];
			MoleculeType2D* molj = &molTypes[jtype];				            // in molj store j-th molecule
			int npj = molj->natoms;						                        // npj is a number of atoms in j-th molecule
			molj->rigidTransform( rot[j], pos[j], Tps_i );	                    // to Tps_j store system coordinates of j-th molecule's atom positions
			//cleanPointForce( npj, fs_j );
			VecN::set( 2*npj, 0.0, (double*)fs_j );				                    // initialize fs_j to zeros
			interMolForce(
				npi, moli->Rs, moli->Qs, Tps_i, fs_i,
				npj, molj->Rs, molj->Qs, Tps_j, fs_j
			);
			if( notFixed[j] ) forceTransform( npj, Tps_j, fs_j,  pos[j],  fpos[j], frot[j] );
		}
		if( notFixed[i] ) forceTransform( npi, Tps_i, fs_i,  pos[i], fpos[i], frot[i] );
	}
}

