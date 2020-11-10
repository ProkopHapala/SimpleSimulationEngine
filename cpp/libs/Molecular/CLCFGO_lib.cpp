
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <vector>
#include <unordered_map>
#include <string>

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
//#include "VecN.h"

#include "integration.h"
#include "AOIntegrals.h"
//#include "AOrotations.h"

#include "Grid.h"
#include "CLCFGO.h"
#include "CLCFGO_tests.h"


// ============ Global Variables

CLCFGO solver;

std::unordered_map<std::string, double*>  buffers;
std::unordered_map<std::string, int*>     ibuffers;

extern "C"{
// ========= Grid initialization

void init( int natom_, int nOrb_, int perOrb_, int natypes_  ){
    solver.realloc( natom_, nOrb_, perOrb_, natypes_ );
    solver.setDefaultValues();
    /*
    int natom =0; ///< number of atoms (nuclei, ions)
    int perOrb=0; //!< Brief number of spherical functions per orbital
    int nOrb  =0; //!< Brief number of single-electron orbitals in system
    // this is evaluated automaticaly
    int nBas  =0; ///< number of basis functions
    int nqOrb =0; ///< number of charges (axuliary density elements) per orbital
    int nQtot =0; ///< total number of charge elements

    // atoms (ions)
    Vec3d*  apos   =0;  ///< positioon of atoms
    Vec3d*  aforce =0;  ///< positioon of atoms
    Vec3d*  aQs    =0;  ///< charge of atom
    int*    atype  =0;  ///< type of atom (in particular IKinetic pseudo-potential)

    // orbitals
    Vec3d*  opos =0;   ///< store positions for the whole orbital
    Vec3d*  odip =0;   ///< Axuliary array to store dipoles for the whole orbital
    double* oEs  =0;   ///< orbital energies
    double* oQs  =0;   ///< total charge in orbital before renormalization (just DEBUG?)
    int*    onq  =0;   ///< number of axuliary density functions per orbital

    // --- Wave-function components for each orbital
    Vec3d*  epos  =0; ///< position of spherical function for expansion of orbitals
    double* esize =0;
    double* ecoef =0;  ///< expansion coefficient of expansion of given orbital
    // --- Forces acting on wave-functions components
    Vec3d*  efpos  =0; ///<   force acting on position of orbitals
    double* efsize =0; ///<   force acting on combination coefficnet of orbitals
    double* efcoef =0; ///<  force acting on size of gaussians

    // --- Auxuliary electron density expansion basis functions
    Vec3d * rhoP  =0; ///< position of density axuliary functio
    double* rhoQ  =0; ///< temporary array to store density projection on pair overlap functions
    double* rhoS  =0;
    // --- Forces acting on auxuliary density basis functions
    Vec3d * rhofP =0; ///< position of density axuliary functio
    double* rhofQ =0; ///< temporary array to store density projection on pair overlap functions
    double* rhofS =0;
    double* rhoEQ =0; /// coulomb energy
    */
    
    // atoms (ions)
    buffers.insert( { "apos",   (double*)solver.apos   } );
    buffers.insert( { "aforce", (double*)solver.aforce } );
    buffers.insert( { "aQs",    (double*)solver.aQs    } );
    ibuffers.insert( { "atype",          solver.atype  } );
    // orbitals
    buffers.insert( { "opos", (double*)solver.opos  } );
    buffers.insert( { "odip", (double*)solver.odip  } );
    buffers.insert( { "oEs",           solver.oEs   } );
    buffers.insert( { "oQs",           solver.oQs   } );
    buffers.insert( { "onq",  (double*)solver.onq   } );
    // --- Wave-function components for each orbital
    buffers.insert( { "epos", (double*)solver.epos   } );
    buffers.insert( { "esize",         solver.esize  } );
    buffers.insert( { "ecoef",         solver.ecoef  } );
    // --- Forces acting on wave-functions components
    buffers.insert( { "efpos", (double*)solver.efpos  } );
    buffers.insert( { "efsize",         solver.efsize } );
    buffers.insert( { "efcoef",         solver.efcoef } );
    // --- Auxuliary electron density expansion basis functions
    buffers.insert( { "rhoP", (double*)solver.rhoP } );
    buffers.insert( { "rhoQ",          solver.rhoQ } );
    buffers.insert( { "rhoS",          solver.rhoS } );
    // --- Forces acting on auxuliary density basis functions
    buffers.insert( { "rhofP", (double*)solver.rhofP } );
    buffers.insert( { "rhofQ",          solver.rhofQ } );
    buffers.insert( { "rhofS",          solver.rhofS } );
    //buffers.insert( { "rhoEQ",          solver.rhoEQ } );

}

/*
int*    getTypes (){ return (int*)   ff.types;  }
double* getPoss  (){ return (double*)ff.poss;   }
double* getQrots (){ return (double*)ff.qrots;  }
double* getHbonds(){ return (double*)ff.hbonds; }
double* getEbonds(){ return (double*)ff.ebonds; }
double* getBondCaps(){ return (double*)ff.bondCaps; }
*/

double* getBuff(const char* name){ 
    auto got = buffers.find( name );
    if( got==buffers.end() ){ return 0;        }
    else                    { return got->second; }
}

void setBuff(const char* name, double* buff){ 
    buffers[name] = buff;
    //auto got = buffers.find( name );
    //if( got==buffers.end() ){ return null;        }
    //else                    { return got->second; }
}

int* getIBuff(const char* name){ 
    auto got = ibuffers.find( name );
    if( got == ibuffers.end() ){ return 0;        }
    else                    { return got->second; }
}

void setIBuff(const char* name, int* buff){ 
    ibuffers[name] = buff;
    //auto got = buffers.find( name );
    //if( got==buffers.end() ){ return null;        }
    //else                    { return got->second; }
}

#define NEWBUFF(name,N)   double* name = new double[N]; buffers.insert( {#name, name} );

void testDerivsP_Coulomb_model( int n, double x0, double dx ){
    //initTestElectrons( );
    for(int i=0; i<solver.nBas; i++){ printf( "epos[%i] (%g,%g,%g) \n", i, solver.epos[i].x,solver.epos[i].y,solver.epos[i].z); }
    solver.toRho(0,1, 0);   
    solver.toRho(2,3, 1);
    NEWBUFF(l_xs,n)
    NEWBUFF(l_Q,n)
    NEWBUFF(l_dQ_ana,n)
    NEWBUFF(l_dQ_num,n)
    NEWBUFF(l_r,n)
    NEWBUFF(l_E,n)
    NEWBUFF(l_Fana,n)
    NEWBUFF(l_Fnum,n)
    solver.bEvalKinetic = false;
    for(int i=0; i<n; i++){
        solver.DEBUG_iter=i;
        solver.cleanForces();
        double x = x0 + i*dx;
        l_xs[i] = x;
        solver.epos[0].x=x;   
        solver.toRho  (0,1, 0);                            // w0*w1 -> q0
        solver.toRho  (2,3, 1);                            // w2*w3 -> q1
        double E_  = solver.CoublombElement(0,1);          // E = C( q0, q1 )
        double E   = E_ * solver.rhoQ[0] * solver.rhoQ[1]; // E(q0,q1) * q0 * q1
        solver.fromRho( 0,1, 0 );                   // w0,w1 <- q0
        l_Q     [i] = solver.rhoQ[0];
        l_dQ_ana[i] = solver.DEBUG_dQdp.x;
        if(i>1) l_dQ_num[i-1]  = (l_Q[i] - l_Q[i-2])/(2*dx);
        l_r[i]       = (solver.rhoP[0] - solver.rhoP[1]).norm();
        l_Fana[i]    = solver.efpos[0].x;                             // This is used when fromRho() is  modified
        l_E[i]       = E; //func( line_E->xs[i], line_Fana->ys[i] );
        if(i>1) l_Fnum[i-1] = (l_E[i] - l_E[i-2])/(2*dx);
    }
    //double* buff = buffers["l_xs"];
    //for(int i=0; i<n; i++){ printf("%i : %g \n", i, buff[i]); }
}

void testDerivsS_Coulomb_model( int n, double x0, double dx ){
    //initTestElectrons( );
    for(int i=0; i<solver.nBas; i++){ printf( "epos[%i] (%g,%g,%g) \n", i, solver.epos[i].x,solver.epos[i].y,solver.epos[i].z); }
    solver.toRho(0,1, 0);   
    solver.toRho(2,3, 1);
    NEWBUFF(l_xs,n)
    NEWBUFF(l_Q,n)
    NEWBUFF(l_dQ_ana,n)
    NEWBUFF(l_dQ_num,n)
    NEWBUFF(l_r,n)
    NEWBUFF(l_E,n)
    NEWBUFF(l_Fana,n)
    NEWBUFF(l_Fnum,n)
    solver.bEvalKinetic = false;
    for(int i=0; i<n; i++){
        solver.DEBUG_iter=i;
        //printf("testDerivsS_Coulomb_model[%i] \n", i );
        solver.cleanForces();
        double x = x0 + i*dx;
        l_xs[i] = x;
        solver.esize[0]=x;   
        solver.toRho  (0,1, 0);                            // w0*w1 -> q0
        solver.toRho  (2,3, 1);                            // w2*w3 -> q1
        double E_  = solver.CoublombElement(0,1);          // E = C( q0, q1 )
        double E   = E_ * solver.rhoQ[0] * solver.rhoQ[1]; // E(q0,q1) * q0 * q1
        solver.fromRho( 0,1, 0 );                   // w0,w1 <- q0
        l_Q     [i] = solver.rhoQ[0];
        l_dQ_ana[i] = solver.DEBUG_dQdp.x;
        if(i>1) l_dQ_num[i-1]  = (l_Q[i] - l_Q[i-2])/(2*dx);
        l_r[i]       = (solver.rhoP[0] - solver.rhoP[1]).norm();
        l_Fana[i]    =  solver.efsize[0];                             // This is used when fromRho() is  modified
        l_E[i]       =  E; //func( line_E->xs[i], line_Fana->ys[i] );
        if(i>1) l_Fnum[i-1] = (l_E[i] - l_E[i-2])/(2*dx);
        //return;
    }
    //double* buff = buffers["l_xs"];
    //for(int i=0; i<n; i++){ printf("%i : %g \n", i, buff[i]); }
}

void testDerivsP_Total( int n, double x0, double dx ){
    for(int i=0; i<solver.nBas; i++){ printf( "epos[%i] (%g,%g,%g) s %g c %g \n", i, solver.epos[i].x,solver.epos[i].y,solver.epos[i].z, solver.esize[i], solver.ecoef[i] ); }
    //NEWBUFF(l_r,n)
    NEWBUFF(l_xs,n)
    NEWBUFF(l_Q,n)
    NEWBUFF(l_E,n)
    NEWBUFF(l_Fana,n)
    NEWBUFF(l_Fnum,n)
    solver.bEvalKinetic = false;
    for(int i=0; i<n; i++){
        solver.DEBUG_iter=i;
        double x = x0 + i*dx;
        l_xs[i] = x;
        solver.epos[0].x=x;
        //printf( ">>> testDerivs_Total [%i] \n", i );
        solver.cleanForces();
        //solver.projectOrb( 0, dip, false );
        //solver.projectOrb( 1, dip, false );
        //double E  = solver.CoulombOrbPair( 0, 1 );
        double E  = solver.eval();
        double xc = solver.rhoP[2].x; //line_Qc->ys[i] = xc;
        l_Fana[i]= solver.efpos[0].x;
        l_E[i]   = E;
        l_Q[i]   = solver.oQs[0];
        //printf( "<<< testDerivs_Total i[%i] x %g E %g \n", i, x, E );
        if(i>1)l_Fnum[i-1] = (l_E[i] - l_E[i-2])/(2*dx);
        //return;
    }
    //for(int i=0;i<n;i++){ line_E->ys[i]   -= line_E->ys[n-1]; };
}

void testDerivsS_Total( int n, double x0, double dx ){
    for(int i=0; i<solver.nBas; i++){ printf( "epos[%i] (%g,%g,%g) s %g c %g \n", i, solver.epos[i].x,solver.epos[i].y,solver.epos[i].z, solver.esize[i], solver.ecoef[i] ); }
    //NEWBUFF(l_r,n)
    NEWBUFF(l_xs,n)
    NEWBUFF(l_Q,n)
    NEWBUFF(l_E,n)
    NEWBUFF(l_Fana,n)
    NEWBUFF(l_Fnum,n)
    solver.bEvalKinetic = false;
    for(int i=0; i<n; i++){
        solver.DEBUG_iter=i;
        double x = x0 + i*dx;
        l_xs[i] = x;
        solver.esize[0]=x;
        //printf( ">>> testDerivs_Total [%i] \n", i );
        solver.cleanForces();
        double E  = solver.eval();
        //double xc = solver.rhoP[2].x; //line_Qc->ys[i] = xc;
        l_Fana[i]= solver.efsize[0];
        //printf( "i %i x %g efsize[0] %g \n", i, x, solver.efsize[0] );
        l_E[i]   = E;
        l_Q[i]   = solver.oQs[0];
        //printf( "<<< testDerivs_Total i[%i] x %g E %g \n", i, x, E );
        if(i>1)l_Fnum[i-1] = (l_E[i] - l_E[i-2])/(2*dx);
        //return;
    }
    //for(int i=0;i<n;i++){ line_E->ys[i]   -= line_E->ys[n-1]; };
}


void testDerivsTotal( int n, double* xs, double* Es, double* Fs, int what ){
    solver.bEvalKinetic = false;
    return testDerivsTotal( solver, n, xs, Es, Fs, what );
}

void testDerivsCoulombModel( int n, double* xs, double* Es, double* Fs, int what  ){
    return testDerivsCoulombModel( solver, n, xs, Es, Fs, what  );
}

} // extern "C"
