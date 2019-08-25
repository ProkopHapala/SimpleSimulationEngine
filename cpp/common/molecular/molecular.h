
#ifndef molecular_h
#define molecular_h

#include "fastmath.h"
#include "Vec3.h"
//#include "Mat3.h"
//#include "quaternion.h"

#define max_bonds 4

class Atom;
class Bond;

class AtomType{
	public:
	int    nbonds;
	double RiiLJ;     // vdW optimal radius
	double EiiLJ;     // vdW energy minimum 
	double rBondMin;  // 
	double rBondMax;  // 
	double rBondOpt;  // 
	double kBond;     // stiff
	double kAngle;
}

class Atom{
	public:
	double    charge;
	AtomType* type;
	Vec3d*    pos;  // this is in independent array (?)
	Vec3d*    force;
 	Bond*     bonds[ max_bonds ];

	inline int bondPossible( double dr ){
		int    iworst  =  -1;
		double dr_worst = 0.0;
		for( int i=0; i<max_bonds; i++ ){
			Bond* bi = bonds[i];
			if( bi !=NULL ){
				double rMaxi = bi->a->type->rBondMax + bi->a->type->rBondMax ; 
				double dri   = Bonds[b]->r - rMaxi;
				if( dri > dr_worst ){ iworst = i; dr_worst = dri; }
			}else{
				//bonds[i] = b;
				return i;
			}
		}
		if( dr_worst > dr ){
			//bonds[iworst]->broken = true;
			//bonds[iworst] = b;
			return iworst;
		}
		return -1;
	}

	inline void addBendForces( ){
		for( int i=max_bonds-1; i>=0; i++ ){
			Bond* a = bonds[i];
			if( a != NULL ){
				for( int j=0; j<max_bonds; i++ ){
				applyBendForce( a, bonds[j] );
			}
		}
	}

	inline void applyBendForce( Bond* a, Bond* b ){
		Vec2d da,db, ftmp;
		bool asg = a->a != this;
		bool bsg = b->a != this;
		da.set( a->dir ); if( asg ){ da.mul(-1); }  /// should be possible to optimize it without "asg" and "bsg" ?
		db.set( b->dir ); if( bsg ){ db.mul(-1); }
		double dot_ab  = da.dot( db );
		double renorm = type->kAngle * ( 1 - dot_ab ); 		
		// force on a
		ftmp.set_mul( a, dot_ab ); ftmp.sub( b ); ftmp.mul( renorm );
		force->sub( ftmp ); 
		if( asg ){	a->b->force->add( ftmp ); }else{ a->a->force->add( ftmp ); }
		// force on b
		ftmp.set_mul( b, dot_ab ); ftmp.sub( a ); ftmp.mul( renorm );
		force->sub( ftmp ); 
		if( bsg ){ b->b->force->sub( ftmp ); }else{ b->a->force->add( ftmp ); }
	}

}

class Bond{
	public:
	Vec3d  dir;
	double r;
	Atom*  a;
	Atom*  b;
	bool broken = false;

	double rOpt;
	double k;

	inline void update(){
		dir.set_sub( *(b->pos), *(a->pos) );
		r = d.norm();
		dir.mul( 1/r );
	}

	inline void applyForce(){
		Vec2d f;
		f.set_mul( dir, k * ( rOpt - r ) );
		a->force->add( f );
		b->force->sub( f );
	}

}


#endif 
