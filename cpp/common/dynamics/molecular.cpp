
#include "molecular.h"

int findBonds( int natoms, Atom * atoms, AtomType * atomTypes, Bond * bonds ){
	int nbonds = 0;
	for( int i=0; i<n; i++ ){
		Atom* atomi = Atom + i;
		Vec3d posi  = *(atomi->pos);
		double rbi  =   atomi->type->rBondMax;
		for( int j=0; j<i; j++ ){     // TODO : this could be accelerated by spatial table (Array,HashMap)
			Atom* atomj = Atom + j;
			Vec d;
			d.set_sub( *(Atom[i].pos), posi );
			double r2   = d.norm2();
			double rbij = rbi + atomj->type->rBondMax;
			// We may also check for rBondMin, to know if geometry is valid
			if( r2 < sq( rbij ) ){
				double r = sqrt( r2 );
				Bond bij = bonds + nbonds;
				int iwi = atomi->bondPossible( dr );
				if( iwi >= 0 ){
					int iwj = atomj->bondPossible( dr );
					if( iwj >= 0 ){
						atomi->bonds[iwi]->broken = true;   
						atomj->bonds[iwj]->broken = true;   
						atomi->bonds[iwi] = bij;
						atomj->bonds[iwj] = bij;
						bij->a=atomi;
						bij->b=btomj;
						bij->r=r;
						bij->dir.set_mul( 1/r );
						nbonds++;
					}
				}
			}
		}
	}
	return nbonds;
}

bool check_disturb( int iatom, int natoms, Atom * atoms, AtomType * atomTypes ){
	// TODO : this could be accelerated by spatial table (Array,HashMap)
	Atom* atomi = Atom + i;
	Vec3d posi  = *(atomi->pos);
	double rbi  =   atomi->type->rBondMin;
	double Etot;
	for( int j=0; j<natoms; j++ ){
		Atom* atomj = Atom + j;
		Vec d;
		d.set_sub( *(Atom[i].pos), posi );
		double r2   = d.norm2();
		double rbij = rbi + atomj->type->rBondMin;
		if( r2 < sq( rbij ) ){ // check if some atoms are not too close
			return false;
		}
		// TODO : we can also compute other leading terms of binding energy here
	}
	//if( Etot > EtotMax ) return false;
	return true;
}

int disturb( int ndisturbs, int retry, double Rdist, int natoms, Atom * atoms, AtomType * atomTypes ){
	int nsuccess = 0;
	for( int i=0; i<ndisturbs; i++ ){
		int iatom = rand() % natoms;
		for( int j=0; j<nretry; j++ ){
			Vec3d d; d.set( randf( -Rdist, Rdist ), randf( -Rdist, Rdist ), randf( -Rdist, Rdist ) );
			Vec3d old_pos; old_pos.set( *atoms[i].pos );   
			if( check_disturb( iatom, natoms, atoms, atomTypes ); ){
				nsuccess++;
			}else{
				atoms[i].pos->set( old_pos );
			}
		}
	}
	return nsuccess;
}


