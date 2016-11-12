
#ifndef forcefield_h
#define forcefield_h

#include <math.h>
#include "fastmath.h"
#include "Vec3.h"

// =============== Forcefield setup

void makeLJparams( int ntypes, double * Rs, double * Es, double*& C6s, double*& C12s ){
	int n2 = ntypes*ntypes;
	C6s  = new double[n2];
	C12s = new double[n2];
	for( int i=0; i<ntypes; i++ ){
		for( int j=0; j<=i; j++ ){
			double r   = Rs[i] + Rs[j];
			double e   = sqrt( Es[i] * Es[j]);
			double r6  = r*r*r; r6*=r6;
			double c6  = -2*e*r6;
			double c12 =    e*r6*r6;
			C6s [ i*ntypes + j ] = c6;
			C12s[ i*ntypes + j ] = c12;
			C6s [ j*ntypes + i ] = c6;
			C12s[ j*ntypes + i ] = c12;
		}
	}
}

// =============== inter-Atomic Force-Filed

inline void forceLJ( const Vec3d& dR, double c6, double c12, Vec3d& f ){
	double ir2  = 1/dR.norm2( );
	double ir6  = ir2*ir2*ir2;
	double ir12 = ir6*ir6;
	f.set_mul( dR , ( 6*ir6*c6 + 12*ir12*c12 ) * ir2   );
}

inline void forceLJE( const Vec3d& dR, double c6, double c12, double qq, Vec3d& f ){
	const double kcoulomb   = 14.3996448915;   //   (e/Angstrome) /(4*pi*epsilon0)  =  ( (1.60217657e-19/1e-10)^2)/(4*pi*8.85418782e-12)
	double ir2  = 1/dR.norm2( );
	double ir6  = ir2*ir2*ir2;
	double ir12 = ir6*ir6;
	double fr = 6*ir6*c6 + 12*ir12*c12;
	//double fr = 12*ir12*c12;
	//double fr = 0;
	if( fabs(qq) > 1e-8 ){ fr += kcoulomb*qq*sqrt(ir2); };
	//fr += kcoulomb*qq*sqrt(ir2);
	f.set_mul( dR ,  fr * ir2 );

	//printf( " %f %f %f %f %f %f \n", qq, sqrt(1/ir2), sqrt(ir2), ir2*sqrt(ir2),  kcoulomb*qq*ir2*sqrt(ir2)*dR.x,    f.x );
}



inline double forceLJE_2( const Vec3d& dR, double rmin, double Emin, double qq, Vec3d& f ){
	const double kcoulomb   = 14.3996448915;   //   (e/Angstrome) /(4*pi*epsilon0)  =  ( (1.60217657e-19/1e-10)^2)/(4*pi*8.85418782e-12)
	double ir2  = 1/dR.norm2( );
	double is2  = (rmin*rmin)/ir2;
	double is6  = is2*is2*is2;
	double vdW   = -2*Emin*is6;
	double pauli =    Emin*is6*is6;
	double elec  = 0;
	if( fabs(qq) > 1e-8 ){ 	elec = kcoulomb*qq*sqrt(ir2);	};
	f.set_mul( dR ,  ( 6*vdW + 12*pauli + elec  ) * -ir2 );
	double  E =          vdW +    pauli + elec;
	return E;
	//printf( " %f %f %f %f %f %f \n", qq, sqrt(1/ir2), sqrt(ir2), ir2*sqrt(ir2),  kcoulomb*qq*ir2*sqrt(ir2)*dR.x,    f.x );
}


// =============== molecule-surface plane

inline void forceMolSurf(
	double z0, double zMin, double 	Emin, const Vec3d& hat,
	int n,  Vec3d * Rs, Vec3d * Fs
){
	for(int i=0; i<n; i++){
		Vec3d dR,f;
		double  z   = hat.dot( Rs[i] ) - z0;
		double iz   = 1/z;
		double ir   = zMin*iz;
		double ir3  = ir*ir*ir;
		double v    = Emin*( ir3*ir3 - 2*ir3 );
		Fs[i].add_mul( hat, v * iz );
		//Fs[i].add_mul( hat, -z*0.1  );
	}
}

// =============== inter-Molecular Force-Filed

inline int interMolForceLJE(
	int na, int * atypes, double * Qas, Vec3d * aRs, Vec3d * aFs,
	int nb, int * btypes, double * Qbs, Vec3d * bRs, Vec3d * bFs,
	int ntypes, double * C6s, double * C12s
){
	//printf( "---------------------\n" );
	//Vec3d fsum; fsum.set(0.0);
	for(int ia=0; ia<na; ia++){
		int atyp  = atypes[ia];
		int ityp0 = ntypes*atyp;
		double qa = Qas[ia];
		for(int ib=0; ib<nb; ib++){
			int btyp   = btypes[ib];
			int ityp   = ityp0 + btyp;
			double C6  = C6s [ ityp ];
			double C12 = C12s[ ityp ];
			Vec3d dR,f;
			dR.set_sub( aRs[ia], bRs[ib] );
			//printf( " %i %i %i %i \n", ia, ib, atyp, btyp );
			forceLJE( dR, C6, C12, qa*Qbs[ib], f );
			aFs[ia].add( f ); bFs[ib].sub( f );
			//aFs[ia].sub( f ); bFs[ib].add( f );

			//printf( " %i %i %i %i   %f   %f   %f %f %f \n", ia, ib, atyp, btyp,    qa*Qbs[ib], f.dot(dR),  f.x, f.y, f.z );
			//printf( " %i %i   %f %f    %f %f %f \n", ia, ib, C6, C12,    f.x, f.y, f.z );
			//fsum.add( f );
		}
	}
	return na*nb;
	//printf( "fsum %f %f %f \n", fsum.x, fsum.y, fsum.z );
}

inline void interMolForceLJ(
	int na, int * atypes, Vec3d * aRs, Vec3d * aFs,
	int nb, int * btypes, Vec3d * bRs, Vec3d * bFs,
	int ntypes, double * C6s, double * C12s
){
	//printf( "---------------------\n" );
	for(int ia=0; ia<na; ia++){
		int atyp  = atypes[ia];
		int ityp0 = ntypes*atyp;
		for(int ib=0; ib<nb; ib++){
			int btyp  = btypes[ib];
			int ityp  = ityp0 + btyp;
			double C6  = C6s[ ityp ];
			double C12 = C12s[ ityp ];
			Vec3d dR,f;
			dR.set_sub( aRs[ia], bRs[ib] );
			forceLJ( dR, C6, C12, f );
			aFs[ia].add( f ); bFs[ib].sub( f );
			//aFs[ia].sub( f ); bFs[ib].add( f );
			//printf( " %i %i   %f %f    %f %f %f \n", ia, ib, C6, C12,    f.x, f.y, f.z );
		}
	}
}

#endif


