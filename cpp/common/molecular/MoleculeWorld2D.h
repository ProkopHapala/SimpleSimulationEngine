
#ifndef MoleculeWorld2D_h
#define MoleculeWorld2D_h

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"
#include "HashMap2D.h"
#include "Body2D.h"

#include "VecN.h"
#include  "DynamicOpt.h"

// ==== Fucntion declarations

void interMolForce(
	int na, double * aRs, double * aQs, Vec2d * aPos, Vec2d * aFs,
	int nb, double * bRs, double * bQs, Vec2d * bPos, Vec2d * bFs
);

void forceTransform( int n, Vec2d * poss, Vec2d * fs,  Vec2d center,  Vec2d& f, double& tq );


// ==== Inline functions

inline void stringForce( const Vec2d& pa, const Vec2d& pb, double k, Vec2d& Fout ){
    Vec2d d;
    d.set_sub( pb, pa );
    Fout.set( d.x*k, d.y*k );
};

/*
inline bool pairwiseForce( const Vec2d& pa, const Vec2d& pb, double qq, double R2, Vec2d& Fout ){
    const double r2max = 4.0d;
    Vec2d d;
    //d.set_sub( pb, pa );
    d.set_sub( pa, pb );
    double r2 = d.norm2();
    if( r2 < r2max ){
        double mr2 = r2max-r2;
        double fr  = ( R2/(r2+0.01) + qq )*mr2*mr2;
        Fout.set( d.x*fr, d.y*fr );
        return true;
    }else{
        Fout.set( 0.0, 0.0 );
        return false;
    }
};
*/


inline bool pairwiseForce( const Vec2d& pa, const Vec2d& pb, double qq, double R2, Vec2d& Fout ){

    Vec2d d;
    //d.set_sub( pb, pa );
    d.set_sub( pa, pb );
    double r2 = d.norm2();

    double s2  = R2/(r2+1e-8);
    double s6  = s2*s2*s2;
    double fr  = s6*s6 - s6 +qq*s2;
    Fout.set( d.x*fr, d.y*fr );

    return true;

    //printf( " (%3.3f,%3.3f) (%3.3f,%3.3f) \n", pa.x ,pa.y,  pb.x ,pb.y );
    //printf( " (%3.3f,%3.3f) (%3.3f,%3.3f) (%3.3f,%3.3f)  %3.3f %3.3f %3.3f  (%3.3f,%3.3f) \n", pa.x ,pa.y,  pb.x ,pb.y, d.x, d.y,  r2, s2, fr, Fout.x, Fout.y );

    /*
    const double r2max = 4.0d*R2;
    if( r2 < r2max ){
        double mr2 = r2max-r2;
        double fr  = ( R2/(r2+0.01) + qq )*mr2*mr2;
        Fout.set( d.x*fr, d.y*fr );
        return true;
    }else{
        Fout.set( 0.0, 0.0 );
        return false;
    }
    */

}



class MoleculeType2D{
    public:
    int       natoms;
    Vec2d   * poss = NULL;
    double  * Qs   = NULL;
    double  * Rs   = NULL;
    //int     * atomTypes = NULL;

    void rigidTransform( double angle, Vec2d pos0, Vec2d * poss_ ){
        Vec2d rot; rot.fromAngle( angle );
        //printf( " angle %3.3f rot (%3.3f,%3.3f) \n", angle, rot.x, rot.y );
        //printf(" Tpos: ");
        for(int i=0; i<natoms; i++){
            //printf( " atom %i \n", i );
            Vec2d pi; pi.set( poss[i] );
            pi.mul_cmplx(rot);
            pi.add(pos0);
            poss_[i].set(pi);
            //printf("(%3.3f,%3.3f)", pi.x, pi.y );
        }
        //printf("\n");

    };

    inline void allocate( int natoms_ ){
        natoms = natoms_;
        poss   = new Vec2d [natoms];
        Qs     = new double[natoms];
        Rs     = new double[natoms];
    }
};

class MoleculeWorld{
	public:

    int              nMolTypes;
    MoleculeType2D * molTypes;

	int      nMols;
	int    * mols     = NULL;
	bool   * notFixed = NULL;

    Vec2d  *pos = NULL, *vpos = NULL, *fpos = NULL;	// positions of molecules, velocities and forces
	double *rot = NULL, *vrot = NULL, *frot = NULL;	// rotations of molecules, velocities and forces

    int nMaxAtoms;
    Vec2d  *Tps_i = NULL, *Tps_j = NULL;
	Vec2d  *fs_i  = NULL, *fs_j  = NULL;

    DynamicOpt optimizer;

    // ==== function declaration

	void forcesMolecules( );

	inline void allocate( int nMolTypes_, int nMols_, int nMaxAtoms_ ){

	    nMolTypes = nMolTypes_; molTypes = new MoleculeType2D[nMolTypes];

        nMols = nMols_;
        mols = new int[nMols]; notFixed = new bool[nMols];
        optimizer.realloc( nMols * 3 );
        pos = (Vec2d*)optimizer.pos; vpos = (Vec2d*)optimizer.vel; fpos = (Vec2d*)optimizer.force;	// positions of molecules, velocities and forces
        rot = optimizer.pos+2*nMols; vrot = optimizer.force+2*nMols; frot = optimizer.force+2*nMols;	// rotations of molecules, velocities and forces

        nMaxAtoms = nMaxAtoms_;
        Tps_i = new Vec2d[nMaxAtoms]; Tps_j = new Vec2d[nMaxAtoms]; fs_i = new Vec2d[nMaxAtoms]; fs_j = new Vec2d[nMaxAtoms];
	};

};

#endif
