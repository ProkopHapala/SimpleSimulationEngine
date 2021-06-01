
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

#include "testUtils.h"

#include "integration.h"
#include "AOIntegrals.h"
//#include "AOrotations.h"


Vec3d DEBUG_dQdp;
int DEBUG_iter     = 0;
int DEBUG_log_iter = 0;
int i_DEBUG=0;

#include "Grid.h"
#include "GaussianBasis.h"
#include "CLCFGO.h"
#include "CLCFGO_tests.h"


// ============ Global Variables

CLCFGO solver;

std::unordered_map<std::string, double*>  buffers;
std::unordered_map<std::string, int*>     ibuffers;

extern "C"{
// ========= Grid initialization

void makeDefaultBuffers(){
    // atoms (ions)
    buffers["apos"  ]= (double*)solver.apos ;
    buffers["aforce"]= (double*)solver.aforce ;
    buffers["aQs"   ]=            solver.aQs  ;
    buffers["aQsize"]=          solver.aQsize ;
    buffers["aPsize"]=          solver.aPsize ;
    buffers["aPcoef"]=          solver.aPcoef ;
    ibuffers["atype"]=          solver.atype  ;
    // orbitals
    buffers["opos"]= (double*)solver.opos  ;
    buffers["odip"]= (double*)solver.odip  ;
    buffers["oEs" ]=           solver.oEs  ;
    buffers["oQs" ]=           solver.oQs  ;
    buffers["onq" ]=  (double*)solver.onq  ;
    ibuffers["ospin"]=  solver.ospin   ;
    // --- Wave-function components for each orbital
    buffers["epos" ]= (double*)solver.epos  ;
    buffers["esize"]=         solver.esize  ;
    buffers["ecoef"]=         solver.ecoef  ;
    // --- Forces acting on wave-functions components
    buffers["efpos" ]= (double*)solver.efpos ;
    buffers["efsize"]=         solver.efsize ;
    buffers["efcoef"]=         solver.efcoef ;
    // --- Forces acting on wave-functions components
    buffers["enfpos" ]= (double*)solver.enfpos ;
    buffers["enfsize"]=          solver.enfsize ;
    buffers["enfcoef"]=          solver.enfcoef ;
    // --- Auxuliary electron density expansion basis functions
    buffers["rhoP"]= (double*)solver.rhoP ;
    buffers["rhoQ"]=          solver.rhoQ ;
    buffers["rhoS"]=          solver.rhoS ;
    // --- Forces acting on auxuliary density basis functions
    buffers["rhofP"]= (double*)solver.rhofP ;
    buffers["rhofQ"]=          solver.rhofQ ;
    buffers["rhofS"]=          solver.rhofS ;
    //buffers.insert( { "rhoEQ",          solver.rhoEQ } );
}

bool loadFromFile( char const* filename, bool bCheck ){
    bool b = solver.loadFromFile( filename, bCheck );
    //solver.setDefaultValues();
    makeDefaultBuffers();
    return b;
}

void init( int natom_, int nOrb_, int perOrb_, int natypes_ ){
    solver.realloc( natom_, nOrb_, perOrb_, natypes_ );
    solver.setDefaultValues();
    makeDefaultBuffers();
}

double  eval(){ return solver.eval(); };
double* getEnergyPointer(){ return &solver.Ek; }
int*    getDimPointer   (){ return &solver.natypes; }

void printSetup            (){ solver.printSetup();                          }
void printAtomsAndElectrons(){ solver.printAtoms(); solver.printElectrons(); }

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

void atomsPotAtPoints(         int n, double* ps, double* out, double s, double Q ){ solver.atomsPotAtPoints( n, (Vec3d*)ps, out, s, Q ); };
void orbAtPoints     ( int io, int n, double* ps, double* out ){ solver.orbAtPoints     ( io, n, (Vec3d*)ps, out ); };
void rhoAtPoints     ( int io, int n, double* ps, double* out ){ solver.rhoAtPoints     ( io, n, (Vec3d*)ps, out ); };
void hartreeAtPoints ( int io, int n, double* ps, double* out ){ solver.hartreeAtPoints ( io, n, (Vec3d*)ps, out ); };

double test_Poisson( int io, double Rmax, double gStep, double * line_rho, double* line_rho_, bool bPrint, bool bSave, bool useWf ){
    return test_Poisson( solver, io, Rmax, gStep, line_rho, line_rho_, bPrint, bSave, useWf );
}

double test_OrbInteraction( int iMODE, int io, int jo, int nint, double dx, double Rmax, double gStep, double * line_Ek, double* line_Ek_g, double * line_f1, double* line_f2, int bPrint, bool bSave ){
    return test_OrbInteraction( solver, iMODE, io, jo, nint, dx, Rmax, gStep, line_Ek, line_Ek_g, line_f1, line_f2, bPrint, bSave );
}

void test_GaussIntegral_ST( int iMODE, int n, double sj, double* sis, double* rs, double* E, double* fr, double* fs ){
    for(int i=0; i<n; i++){
        double si = sis[i];
        double r_ = rs [i];
        double r2 = r_*r_;
        _Gauss_sij_aux( si, sj ) 
        if       ( iMODE==1 ){
            _Gauss_overlap( r2, si, sj );
            E[i]=S; fr[i]=dS_dr*r_; fs[i]=dS_dsi;   
        }else if ( iMODE==2 ) {
            _Gauss_tau    ( r2, si, sj );
            E[i]=tau; fr[i]=dTau_dr*r_; fs[i]=dTau_dsi;
        }else if ( iMODE==3 ) { 
            _Gauss_overlap( r2, si, sj );
            _Gauss_tau    ( r2, si, sj );
            _Gauss_kinetic( r2, si, sj );
            E[i]=T; fr[i]=dT_dr*r_; fs[i]=dT_dsi;
        }
    }
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
        DEBUG_iter=i;
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
        l_dQ_ana[i] = DEBUG_dQdp.x;
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
        DEBUG_iter=i;
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
        l_dQ_ana[i] = DEBUG_dQdp.x;
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
        DEBUG_iter=i;
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
        DEBUG_iter=i;
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
    /*
    solver.bNormalize     = 0;
    solver.bEvalKinetic   = 0;
    solver.bEvalCoulomb   = 0;
    solver.bEvalExchange  = 0;
    solver.bEvalPauli     = 0;
    solver.iPauliModel    = 1;
    solver.bEvalAA        = 0;
    solver.bEvalAE        = 1;
    solver.bEvalAECoulomb = 1;
    solver.bEvalAEPauli   = 0;
    */
    return testDerivsTotal( solver, n, xs, Es, Fs, what );
}

void testDerivsCoulombModel( int n, double* xs, double* Es, double* Fs, int what  ){
    return testDerivsCoulombModel( solver, n, xs, Es, Fs, what  );
}

void setSwitches(bool bNormalize, bool bEvalKinetic, bool bEvalCoulomb, bool  bEvalExchange, bool  bEvalPauli, int iPauliModel,  bool bEvalAA, bool  bEvalAE, bool  bEvalAECoulomb, bool  bEvalAEPauli ){
    solver.bNormalize     = bNormalize;
    solver.bEvalKinetic   = bEvalKinetic;
    solver.bEvalCoulomb   = bEvalCoulomb;
    solver.bEvalExchange  = bEvalExchange;
    solver.bEvalPauli     = bEvalPauli;
    solver.iPauliModel    = iPauliModel;
    solver.bEvalAA        = bEvalAA;
    solver.bEvalAE        = bEvalAE;
    solver.bEvalAECoulomb = bEvalAECoulomb;
    solver.bEvalAEPauli   = bEvalAEPauli;
}

void setSwitches_(int bNormalize, int bNormForce, int bEvalKinetic, int bEvalCoulomb, int  bEvalExchange, int  bEvalPauli, int bEvalAA, int bEvalAE, int bEvalAECoulomb, int bEvalAEPauli ){
#define _setbool(b,i) { if(i>0){b=true;}else if(i<0){b=false;} }
    _setbool( solver.bNormalize     , bNormalize     );
    _setbool( solver.bNormForce     , bNormForce     );
    _setbool( solver.bEvalKinetic   , bEvalKinetic   );
    _setbool( solver.bEvalCoulomb   , bEvalCoulomb   );
    _setbool( solver.bEvalExchange  , bEvalExchange  );
    _setbool( solver.bEvalPauli     , bEvalPauli     );
    _setbool( solver.bEvalAA        , bEvalAA        );
    _setbool( solver.bEvalAE        , bEvalAE        );
    _setbool( solver.bEvalAECoulomb , bEvalAECoulomb );
    _setbool( solver.bEvalAEPauli   , bEvalAEPauli   );
#undef _setbool
}

void setPauli( int iPauli ){ solver.iPauliModel = iPauli; }

} // extern "C"
