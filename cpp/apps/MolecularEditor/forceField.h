
#ifndef forcefield_h
#define forcefield_h

#include <math.h>
#include "fastmath.h"
#include "Vec3.h"


namespace ForceField{
	static int      ntypes  = 0;
	static double * C6s     = NULL;
    static double * C12s    = NULL;
    static double   Rcut,Rcut2;
    static double * FcutCoulomb = NULL;
    static double * FcutLJ      = NULL;

    inline void forceLJE( const Vec3d& dR, double c6, double c12, double qq, double fcut, Vec3d& f ){
        //printf("forceLJE %g %g   (%3.3f,%3.3f,%3.3f)\n", c6, c12, dR.x, dR.y, dR.z);
        const double kcoulomb   = 14.3996448915;   //   (e/Angstrome) /(4*pi*epsilon0)  =  ( (1.60217657e-19/1e-10)^2)/(4*pi*8.85418782e-12)
        double ir2  = 1/dR.norm2( );
        double ir6  = ir2*ir2*ir2;
        double ir12 = ir6*ir6;
        double fr = 6*ir6*c6 + 12*ir12*c12;
        //double fr = 12*ir12*c12;
        //double fr = 0;
        //if( fabs(qq) > 1e-8 ){ fr += kcoulomb*qq*sqrt(ir2); };
        //fr += fcut; // NOTE : fcut is not true force - it is divided by |r| in order to save one sqrt()
        double ir = sqrt(ir2);
        fr += kcoulomb*qq*ir;
        //fr += kcoulomb*qq*sqrt(ir2);
        f.set_mul( dR , fr * ir2  + fcut*ir );
        //printf( " %f %f %f %f %f %f \n", qq, sqrt(1/ir2), sqrt(ir2), ir2*sqrt(ir2),  kcoulomb*qq*ir2*sqrt(ir2)*dR.x,    f.x );
    };

    void init( double Rcut_, int ntypes_, double * RLJs, double * ELJs ){
        ntypes = ntypes_;
        Rcut   = Rcut_; Rcut2 = Rcut*Rcut;
        int n2 = ntypes*ntypes;
        C6s          = new double[n2];
        C12s         = new double[n2];
        FcutCoulomb  = new double[n2];
        FcutLJ       = new double[n2];
        for( int i=0; i<ntypes; i++ ){
            for( int j=0; j<=i; j++ ){
                double r   = RLJs[i] + RLJs[j];
                double e   = sqrt( ELJs[i] * ELJs[j]);
                double r6  = r*r*r; r6*=r6;
                double c6  = -2*e*r6;
                double c12 =    e*r6*r6;

                int ii = i*ntypes + j;
                int jj = j*ntypes + i;

                C6s [ii] = c6;  C6s [jj] = c6;
                C12s[ii] = c12; C12s[jj] = c12;

                Vec3d dR,f; dR.set(0.0,0.0,Rcut);
                forceLJE( dR, c6, c12, 0.0, 0.0, f ); double fcutLJ      = f.z;
                forceLJE( dR, c6, c12, 1.0, 0.0, f ); double fcutCoulomb = f.z - fcutLJ;

                FcutLJ     [ii] = fcutLJ;      FcutLJ     [jj] = fcutLJ;
                FcutCoulomb[ii] = fcutCoulomb; FcutCoulomb[jj] = fcutCoulomb;

            }
        }
    }

};


// =============== Forcefield setup

/*
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
*/

// =============== inter-Atomic Force-Filed

inline Vec3d radialSpringForce( const Vec3d& dp, double k, double l0){

    Vec3d f;
    //printf( "r2 %g  l0 %g\n", sqrt(r2), l0 );
    double r = dp.norm();
    f.set_mul( dp, k*(l0-r)/r );

    /*
    double r2 = dp.norm2();
    if( r2 > l0*l0 ){
        double r = sqrt(r2);
        f.set_mul( dp, k*(r-l0)/r );
    }else{
        f.set(0.0);
    }
    */
    return f;
}


inline Vec3d radialBondForce( const Vec3d& dp, double k, double l0, double dlmax){

    Vec3d f;
    //printf( "r2 %g  l0 %g\n", sqrt(r2), l0 );
    double r  = dp.norm();
    double dr = r-l0;
    if( dr > dlmax ) dr=dlmax;
    f.set_mul( dp, -k*dr/r );
    return f;
}

inline void forceLJ( const Vec3d& dR, double c6, double c12, Vec3d& f ){
	double ir2  = 1/dR.norm2( );
	double ir6  = ir2*ir2*ir2;
	double ir12 = ir6*ir6;
	f.set_mul( dR , ( 6*ir6*c6 + 12*ir12*c12 ) * ir2   );
}

inline void forceLJE( const Vec3d& dR, double c6, double c12, double qq, Vec3d& f ){
    //printf("forceLJE %g %g   (%3.3f,%3.3f,%3.3f)\n", c6, c12, dR.x, dR.y, dR.z);
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
	int nb, int * btypes, double * Qbs, Vec3d * bRs, Vec3d * bFs
){
	//printf( "---------------------\n" );
	//Vec3d fsum; fsum.set(0.0);
	for(int ia=0; ia<na; ia++){
		int atyp  = atypes[ia];
		int ityp0 = ForceField::ntypes*atyp;
		double qa = Qas[ia];
		for(int ib=0; ib<nb; ib++){
            Vec3d dR;
            dR.set_sub( aRs[ia], bRs[ib] );

            double r2 = dR.norm2();
            if( r2>ForceField::Rcut2 ) continue;

			int btyp   = btypes[ib];
			int ityp   = ityp0 + btyp;
			double C6  = ForceField::C6s [ ityp ];
			double C12 = ForceField::C12s[ ityp ];
			//printf( "interMolForceLJE %i %i  %i %i  %i %i  %g %g \n", ia, ib, atyp, btyp, ityp, ntypes,   C6, C12 );
			Vec3d f;

			//printf( " %i %i %i %i \n", ia, ib, atyp, btyp );

			//forceLJE( dR, C6, C12, qa*Qbs[ib], f );

			double qq   = qa*Qbs[ib];
			double fcut = ForceField::FcutLJ[ityp] + ForceField::FcutCoulomb[ityp]*qq;
			//ForceField::forceLJE( dR, C6, C12, qq, -fcut, f );
			ForceField::forceLJE( dR, C6, C12, qq, 0.0, f );

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


