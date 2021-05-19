
#ifndef CLCFGO_tests_h
#define CLCFGO_tests_h

#include "integration.h"
#include "AOIntegrals.h"
//#include "AOrotations.h"
#include "Grid.h"
#include "GridUtils.h"
#include "CLCFGO.h"
#include "VecN.h"


#include "Plot2D.h"

#define iTEST_POS_DERIV  0
#define iTEST_SIZE_DERIV 1

void testDerivsCoulombModel( CLCFGO& solver, int n, double* xs, double* Es, double* Fs, int what  ){
    //initTestElectrons( );
    solver.toRho(0,1, 0);
    solver.toRho(2,3, 1);
    for(int i=0; i<n; i++){
        DEBUG_iter=i;
        solver.clearAuxDens();
        solver.cleanForces();
        switch(what){
            case iTEST_POS_DERIV : solver.epos [0].x=xs[i];       break;
            case iTEST_SIZE_DERIV: solver.esize[0]  =fabs(xs[i]); break;
        }
        solver.toRho  (0,1, 0);                          // w0*w1 -> q0
        solver.toRho  (2,3, 1);                          // w2*w3 -> q1
        double E_  = solver.CoublombElement(0,1);        // E = C( q0, q1 )
        double E   = E_ * solver.rhoQ[0]*solver.rhoQ[1]; // E(q0,q1) * q0 * q1
        solver.fromRho( 0,1, 0 );                   // w0,w1 <- q0
        Es[i] = E;
        switch(what){
            case iTEST_POS_DERIV : Fs[i]=solver.efpos [0].x; break;
            case iTEST_SIZE_DERIV: Fs[i]=solver.efsize[0];   break;
        }
        //l_Q     [i] = solver.rhoQ[0];
        //l_dQ_ana[i] = solver.DEBUG_dQdp.x;
    }
}

void testDerivsTotal( CLCFGO& solver, int n, double* xs, double* Es, double* Fs, int what ){
    solver.reportOrbitals();
    for(int i=0; i<n; i++){
        //printf( "===== testDerivsTotal[%i]\n", i  );
        DEBUG_iter=i;
        //if(DEBUG_iter==DEBUG_log_iter){ printf("before switch(what) \n"); solver.reportOrbitals(); }
        switch(what){
            case iTEST_POS_DERIV : solver.epos [0].x=xs[i];        break;
            case iTEST_SIZE_DERIV: solver.esize[0]  =fabs(xs[i]);  break;
        }
        //if(DEBUG_iter==DEBUG_log_iter){ printf("before cleanForces() \n"); solver.reportOrbitals(); }
        solver.cleanForces();
        //if(DEBUG_iter==DEBUG_log_iter){ printf("after cleanForces() \n"); solver.reportOrbitals(); }
        double E  = solver.eval();
        //double xc = solver.rhoP[2].x; //This may be in aux struct
        Es[i] = E;
        switch(what){
            case iTEST_POS_DERIV : Fs[i]=solver.efpos [0].x; break;
            case iTEST_SIZE_DERIV: Fs[i]=solver.efsize[0];   break;
        }
    }
}

void testDerivsTotal( int n, double x0, double dx, CLCFGO& solver, Plot2D& plot, int iline0=-1 ){
    DataLine2D *line_E,*line_Fnum,*line_Fana;
    double x_bak=0;
    if(iline0<0){
        //printf( "n  %i dx %g  \n", n , dx );
        line_E    = new DataLine2D( n, x0, dx, 0xFFFF0000, "E"     ); plot.add(line_E    );
        line_Fnum = new DataLine2D( n, x0, dx, 0xFF0080FF, "F_num" ); plot.add(line_Fnum );
        line_Fana = new DataLine2D( n, x0, dx, 0xFF0000FF, "F_ana" ); plot.add(line_Fana );
    }else{
        //printf( "testDerivsTotal update \n" );
        line_E    = plot.lines[iline0+0];
        line_Fnum = plot.lines[iline0+1];
        line_Fana = plot.lines[iline0+2];
        n=line_E->n;
    }
    x_bak = solver.epos[0].x;
    //double DEBUG_sum = 0.0;
    for(int i=0; i<n; i++){
        //solver.cleanForces(); // this is already in solver.eval()
        solver.epos[0].x  =  line_E->xs[i];
        line_E   ->ys[i]  =  solver.eval();
        //DEBUG_sum        +=  line_E->ys[i];
        line_Fana->ys[i]  = -solver.efpos[0].x;
        if(i>1)line_Fnum->ys[i-1] = (line_E->ys[i] - line_E->ys[i-2])/(2*dx);
    }
    solver.epos[0].x = x_bak;
    //printf( "DEBUG_sum %g size %g \n", DEBUG_sum, solver.esize[0] );
}

void plotOrb( CLCFGO& solver, DataLine2D *line, int io, Vec3d p0, Vec3d dp, float sc=1.0 ){
    Vec3d ps[line->n];
    for(int i=0; i<line->n; i++){  ps[i]=p0+dp*line->xs[i]; }
    solver.orbAtPoints( io, line->n, ps, line->ys );
    for(int i=0; i<line->n; i++){  line->ys[i]*=sc; }
}

void plotAtomsPot( CLCFGO& solver, DataLine2D *line, Vec3d p0, Vec3d dp, float sc=1.0, double s=0.0 ){
    Vec3d ps[line->n];
    for(int i=0; i<line->n; i++){  ps[i]=p0+dp*line->xs[i]; }
    //solver.orbAtPoints( io, line->n, ps, line->ys );
    solver.atomsPotAtPoints( line->n, ps, line->ys, s, 1.0 );
    for(int i=0; i<line->n; i++){  line->ys[i]*=sc; }
}



template<typename Func>
void testDerivs( int n, double x0, double dx, CLCFGO& solver, Plot2D& plot1, Func func ){
    // ======= Test Orbital Wavefunction Overlap
    //printf( "n  %i dx %g  \n", n , dx );
    DataLine2D* line_E    = new DataLine2D( n, x0, dx, 0xFFFF0000, "E"     ); plot1.add(line_E    );
    DataLine2D* line_Fnum = new DataLine2D( n, x0, dx, 0xFF0080FF, "F_num" ); plot1.add(line_Fnum );
    DataLine2D* line_Fana = new DataLine2D( n, x0, dx, 0xFF0000FF, "F_ana" ); plot1.add(line_Fana );
    for(int i=0; i<n; i++){
        solver.cleanForces();
        //solver.epos[2].x = line_E->xs[i];
        //line_E   ->ys[i] = solver.evalOverlap( 0, 1 );
        //line_Fana->ys[i] = solver.efpos[0].x;
        line_E->ys[i] =  func( line_E->xs[i], line_Fana->ys[i] );
        if(i>1)line_Fnum->ys[i-1] = (line_E->ys[i] - line_E->ys[i-2])/(2*dx);
    }
}


void testDerivs_Coulomb( int n, double x0, double dx, CLCFGO& solver, Plot2D& plot1 ){
    // ======= Test Orbital Wavefunction Overlap
    //printf( "n  %i dx %g  \n", n , dx );
    DataLine2D* line_E    = new DataLine2D( n, x0, dx, 0xFFFF0000, "E"     ); plot1.add(line_E    );
    DataLine2D* line_Fnum = new DataLine2D( n, x0, dx, 0xFF0080FF, "F_num" ); plot1.add(line_Fnum );
    DataLine2D* line_Fana = new DataLine2D( n, x0, dx, 0xFF0000FF, "F_ana" ); plot1.add(line_Fana );

    DataLine2D* line_Qc   = new DataLine2D( n, x0, dx, 0xFFFF80FF, "Qc" ); plot1.add(line_Qc );
    DataLine2D* line_Xc   = new DataLine2D( n, x0, dx, 0xFFFF00FF, "Xc" ); plot1.add(line_Xc );
    solver.bNormalize = false;
    for(int i=0; i<n; i++){
        solver.cleanForces();
        double x = x0 + i*dx;
        solver.epos[0].x=x;
        Vec3d dip;
        solver.projectOrb( 0, dip );
        solver.projectOrb( 1, dip );
        //solver.evalElectrostatICoulomb();

        double xc = solver.rhoP[2].x; line_Qc->ys[i] = xc;

        //printf( "line_Qc->ys[i] %g %g | %g \n", line_Qc->ys[i], line_Qc->xs[i],     ( line_Qc->ys[i]-line_Qc->ys[0] ) / ( line_Qc->xs[i] - line_Qc->xs[0] )   );

        //double qc = solver.rhoQ[0];   line_Xc->ys[i] = qc;
        double E  = solver.CoulombOrbPair( 0, 1 );
        solver.assembleOrbForces_fromRho(0);
        //f=solver.efpos[0].x * M_PI * 2;
        //line_Fana->ys[i]= solver.efpos[0].x * M_PI * 2;
        line_Fana->ys[i]= solver.efpos[0].x * 7;
        line_E->ys[i]   = E; //func( line_E->xs[i], line_Fana->ys[i] );
        if(i>1)line_Fnum->ys[i-1] = (line_E->ys[i] - line_E->ys[i-2])/(2*dx);
    }
}



void initTestElectrons( CLCFGO& solver ){
    /*
    {auto& _=solver;
        _.ecoef[0] =   1.0;
        _.ecoef[1] =   1.0;
        _.ecoef[2] =  -0.5;
        _.ecoef[3] =   1.0;
        _.epos [0] = (Vec3d){ 0.0, 0.0,0.0};
        _.epos [1] = (Vec3d){ 0.0, 0.0,0.0};
        _.epos [2] = (Vec3d){-3.0,-0.1,0.0};
        _.epos [3] = (Vec3d){-3.0,+1.5,0.0};
        //_.ecoef[3] = +0.3;
    }
    */
    {auto& _=solver;
        _.ecoef[0] =   1.0;
        _.ecoef[1] =   1.0;
        _.ecoef[2] =   1.0;
        _.ecoef[3] =   1.0;

        _.esize[0] =   1.0;
        _.esize[1] =   1.0;
        _.esize[2] =   1.0;
        _.esize[3] =   1.0;

        _.epos [0] = (Vec3d){ 0.0, 0.0,0.0};
        _.epos [1] = (Vec3d){ 0.0, 0.0,0.0};
        _.epos [2] = (Vec3d){-1.5, 0.0,0.0};
        _.epos [3] = (Vec3d){ 0.0, 0.0,0.0};
        //_.ecoef[3] = +0.3;
    }
}


/*


### Coulomb Force Derivatives
        evalElectrostatICoulomb();
        for(int i=0; i<nOrb; i++) assembleOrbForces(i);
   Problem is derivative of Q
     qi = Sab(ra-rb) * Ca * Cb
     qj = Scd(rc-rd) * Cc * Cd                   ... where Ca,Cb,Cc,Cd are expansion coeficients in given basis and S is overlap integral for given pair of basis functions
   EQij(r) = qi(ra,rb) * qj(rc,rd) * Kij(ri,rj)  ... where qi,qj are charges in some overlap cloud and Kij is Coulomb Matrix kernel between the two clouds
   Force calculated by derivatives as:
     FQ_xa = dEQ/dxa =  (qi*qj) * (dKij/ri)/(dri/dxa) + ( Kij*qj ) * (dqi/dxa)
*/


void testDerivs_Total( int n, double x0, double dx, CLCFGO& solver, Plot2D& plot1 ){
    DataLine2D* line_E    = new DataLine2D( n, x0, dx, 0xFFFF0000, "E"     ); plot1.add(line_E    );
    DataLine2D* line_Fnum = new DataLine2D( n, x0, dx, 0xFF0080FF, "F_num" ); plot1.add(line_Fnum );
    DataLine2D* line_Fana = new DataLine2D( n, x0, dx, 0xFF0000FF, "F_ana" ); plot1.add(line_Fana );
    DataLine2D* line_Q    = new DataLine2D( n, x0, dx, 0xFF00FF00, "Q" );     plot1.add(line_Q    );
    initTestElectrons( solver );
    for(int i=0; i<n; i++){
        solver.cleanForces();
        double x = x0 + i*dx;
        solver.epos[0].x=x;
        //solver.projectOrb( 0, dip, false );
        //solver.projectOrb( 1, dip, false );
        double xc = solver.rhoP[2].x; //line_Qc->ys[i] = xc;
        //double E  = solver.CoulombOrbPair( 0, 1 );
        double E  = solver.eval();

        line_Fana->ys[i]= solver.efpos[0].x;
        line_E->ys[i]   = E;
        line_Q->ys[i]   = solver.oQs[0];


        printf( "testDerivs_Total i[%i] x %g E %g \n", i, x, E );
        if(i>1)line_Fnum->ys[i-1] = (line_E->ys[i] - line_E->ys[i-2])/(2*dx);
    }
    //for(int i=0;i<n;i++){ line_E->ys[i]   -= line_E->ys[n-1]; };
}


void testDerivs_Coulomb_model( int n, double x0, double dx, CLCFGO& solver, Plot2D& plot1 ){
    // ======= Test Orbital Wavefunction Overlap
    //printf( "n  %i dx %g  \n", n , dx );

    //DataLine2D* line_Eq    = new DataLine2D( n, x0, dx, 0xFFFF8080, "Eq"     ); //plot1.add(line_Eq    );
    //DataLine2D* line_Fqnum = new DataLine2D( n, x0, dx, 0xFF008077, "Fq_num" ); //plot1.add(line_Fqnum ); // orange
    //DataLine2D* line_Fqana = new DataLine2D( n, x0, dx, 0xFF000077, "Fq_ana" ); //plot1.add(line_Fqana ); // red     // efpos[0].x;
    DataLine2D* line_E     = new DataLine2D( n, x0, dx, 0xFFFF0000, "E"     ); plot1.add(line_E    );
    DataLine2D* line_Fnum  = new DataLine2D( n, x0, dx, 0xFF0077FF, "F_num" ); plot1.add(line_Fnum ); // orange
    DataLine2D* line_Fana  = new DataLine2D( n, x0, dx, 0xFF0000FF, "F_ana" ); plot1.add(line_Fana ); // red     // efpos[0].x;

    //DataLine2D* line_Frho = new DataLine2D( n, x0, dx, 0xFF00FFFF, "F_rho" ); plot1.add(line_Frho ); // yellow  // rhofP[0].x

    DataLine2D* line_dQi_num  = new DataLine2D( n, x0, dx, 0xFF080FF00, "dQi_num" ); //plot1.add(line_dQi_num );
    DataLine2D* line_dQi_ana  = new DataLine2D( n, x0, dx, 0xFF0FF8000, "dQi_ana" ); //plot1.add(line_dQi_ana );

    DataLine2D* line_Qi   = new DataLine2D( n, x0, dx, 0xFF008000, "Qi" ); plot1.add(line_Qi );
    //DataLine2D* line_Si   = new DataLine2D( n, x0, dx, 0xFFFF00FF, "Si" ); //plot1.add(line_Si );
    DataLine2D* line_Pi   = new DataLine2D( n, x0, dx, 0xFFFFFFFF, "Pi" ); plot1.add(line_Pi );
    //DataLine2D* line_Qj   = new DataLine2D( n, x0, dx, 0xFF00F000, "Qj" ); plot1.add(line_Qj );

    //DataLine2D* line_S1   = new DataLine2D( n, x0, dx, 0xFF8800FF, "S1" ); //plot1.add(line_S1 );
    //DataLine2D* line_S2   = new DataLine2D( n, x0, dx, 0xFFFF0088, "S2" ); //plot1.add(line_S2 );

    //DataLine2D* line_Qc   = new DataLine2D( n, x0, dx, 0xFFFF80FF, "Qc" ); plot1.add(line_Qc );
    //DataLine2D* line_Xc   = new DataLine2D( n, x0, dx, 0xFFFF00FF, "Xc" ); plot1.add(line_Xc );


    initTestElectrons( solver );

    solver.toRho(0,1,0);
    solver.toRho(2,3,1);

    //int ie=0,je=1;
    for(int i=0; i<n; i++){
        //printf( "-------------- %i | testDerivs_Coulomb_model \n", i );
        solver.cleanForces();
        double x = x0 + i*dx + 0.01;
        solver.epos[0].x=x;    // set position of  wf basis function   xhi[0]
        //solver.rhoP[0].x=x;      // set position of  density blob        rho[0]

        // project wave-function (molecular orbitals) to auxiliary density basis
        solver.toRho  (0,1,0);
        solver.toRho  (2,3,1);
        // evaluate coulomb interaction between density basis functions
        double E_  = solver.CoublombElement(0,1);
        //line_Eq->ys[i] = E_;
        double E   = E_ * solver.rhoQ[0];
        double aij;
        solver.fromRho( 0,1,0 );
        //solver.fromRho( 0,1,0,   aij );

        //E /= solver.rhoQ[0];

        //line_Eq->ys[i] = E /( -2* solver.rhoQ[0] );
        //line_S1->ys[i] = solver.esize[0];
        //line_S2->ys[i] = solver.esize[1];

        line_Qi     ->ys[i] = solver.rhoQ[0];
        //line_dQi_ana->ys[i] = solver.DEBUG_dQdp.x;
        //line_dQi_ana->ys[i] = dQ * -x * 0.5;
        if(i>1){
            line_dQi_num->ys[i-1]  = (line_Qi->ys[i] - line_Qi->ys[i-2])/(2*dx);
        //    printf( "[%i] /%g   num %g  ana %g \n", i, line_dQi_num->ys[i-1]/line_dQi_ana->ys[i-1], line_dQi_num->ys[i-1], line_dQi_ana->ys[i-1] );
        }
        //line_Si->ys[i] = solver.rhoS[0];
        //line_Pi->ys[i] = solver.rhoP[0].x;
        //line_Pi->ys[i] = solver.rhoP[0].x - solver.rhoP[1].x;
        line_Pi->ys[i] = (solver.rhoP[0] - solver.rhoP[1]).norm();


        //double xc = solver.rhoP[0].x;
        // x-axis according to density blob
        //line_E   ->xs[i] = xc;
        //line_Fnum->xs[i] = xc;
        //line_Fana->xs[i] = xc;
        //line_Fana->ys[i] = solver.rhofP[0].x;

        //line_Frho->ys[i] = solver.rhofP[0].x;

        //printf( "[%i]", i );

        //solver.fromRho(2,3,1);
        //line_Fqana->ys[i] = 0.5*solver.efpos[0].x/line_Qi->ys[i];  // This is derivative of force for Q = const. (  dQ/dx =0 )
        //line_Fana->ys[i]  = 0.5*solver.efpos[0].x + E_*line_dQi_ana->ys[i];    // This is used with fromRho() not modified
        //line_Fana->ys[i]  = 0.5*solver.efpos[0].x + E_*dQdp.x;                   // This is used with fromRho() not modified
        line_Fana->ys[i]  = solver.efpos[0].x;                             // This is used when fromRho() is  modified
        //printf( "testDerivs_Coulomb_model E_ %g dQdx %g \n", E_, dQdp.x );

        // TODO :  we should first make work this "Fana" with fromRho()
        //         Why is there coef 0.5 ?


        line_E->ys[i]   = E; //func( line_E->xs[i], line_Fana->ys[i] );
        if(i>1) line_Fnum ->ys[i-1] = (line_E ->ys[i] - line_E ->ys[i-2])/(2*dx);
        //if(i>1) line_Fqnum->ys[i-1] = (line_Eq->ys[i] - line_Eq->ys[i-2])/(2*dx);

        //if(i>1) line_dQi->ys[i-1]  = (line_Qi->ys[i] - line_Qi->ys[i-2])/(2*dx);

        //   d(E*Q)/dx = (dE/dx) * Q   +   E * (dQ/dx);
        //if(i>1) line_Fana->ys[i]   =  line_Fana->ys[i-1]*line_Qi->ys[i-1]   +  line_E->ys[i-1]*line_dQi->ys[i-1];

    }
}



void testDerivs_Coulomb_model_S( int n, double x0, double dx, CLCFGO& solver, Plot2D& plot1 ){

    //DataLine2D* line_Eq    = new DataLine2D( n, x0, dx, 0xFFFF8080,  "Eq"     ); //plot1.add(line_Eq    );
    //DataLine2D* line_Fqnum = new DataLine2D( n, x0, dx, 0xFF008077, "Fq_num" ); //plot1.add(line_Fqnum ); // orange
    //DataLine2D* line_Fqana = new DataLine2D( n, x0, dx, 0xFF000077, "Fq_ana" ); //plot1.add(line_Fqana ); // red     // efpos[0].x;

    //DataLine2D* line_E     = new DataLine2D( n, x0, dx, 0xFFFF0000, "E"     ); plot1.add(line_E    );
    //DataLine2D* line_Fnum  = new DataLine2D( n, x0, dx, 0xFF0077FF, "F_num" ); plot1.add(line_Fnum ); // orange
    //DataLine2D* line_Fana  = new DataLine2D( n, x0, dx, 0xFF0000FF, "F_ana" ); plot1.add(line_Fana ); // red     // efpos[0].x;

    //DataLine2D* line_dQi_num = new DataLine2D( n, x0, dx, 0xFFFF8888, "dQi_num" ); plot1.add(line_dQi_num );
    //DataLine2D* line_dQi_ana = new DataLine2D( n, x0, dx, 0xFFFF0000, "dQi_ana" ); plot1.add(line_dQi_ana );

    DataLine2D* line_Qi      = new DataLine2D( n, x0, dx, 0xFF008800, "Qi" ); plot1.add(line_Qi );
    DataLine2D* line_aij     = new DataLine2D( n, x0, dx, 0xFF888800, "aij" );   //plot1.add(line_aij );
    DataLine2D* line_Qaij    = new DataLine2D( n, x0, dx, 0xFF008888, "Q/aij" ); //plot1.add(line_Qaij );

    solver.toRho(0,1,0);
    solver.toRho(2,3,1);

    //int ie=0,je=1;
    for(int i=0; i<n; i++){
        solver.cleanForces();
        double x = x0 + i*dx;
        solver.esize[0]=x;        // set size of wf basis function  xhi[0]

        solver.toRho(0,1,0);
        solver.toRho(2,3,1);
        double E_q = solver.CoublombElement(0,1);
        double Q   = solver.rhoQ[0];
        double E   = E_q * Q;

        double aij;
        //double dCsi, dCsj; Vec3d dQdp;
        //solver.fromRho(0,1,0, aij, dCsi, dCsj, dQdp);
        solver.fromRho(0,1,0 );
        //Vec3d dQdp =  solver.fromRho(0,1,0);

        line_Qi     ->ys[i] = Q;
        line_aij    ->ys[i] = aij;
        line_Qaij   ->ys[i] = Q/aij;
        //printf( "Q/aij %g  \n", Q/aij );

        //line_dQi_ana->ys[i] = dCsi*(Q/aij); // ToDo : this works but there should be better way (using only cij, without overlap sij)

        //solver.fromRho(2,3,1);
        //line_Fqana->ys[i] = 0.5*solver.efsize[0]/line_Qi->ys[i];            // This is derivative of force for Q = const. (  dQ/dx =0 )
        //line_Fana ->ys[i] = 0.5*solver.efsize[0] + E_q*line_dQi_ana->ys[i];

        //line_Eq->ys[i] = E_q;
        //line_E->ys[i] = E;
        //if(i>1) line_Fnum-> ys[i-1] = ( line_E ->ys[i] - line_E ->ys[i-2] )/(2*dx);
        //if(i>1) line_Fqnum->ys[i-1] = ( line_Eq->ys[i] - line_Eq->ys[i-2] )/(2*dx);
        //if(i>1){
        //    line_dQi_num->ys[i-1] = ( line_Qi->ys[i] - line_Qi->ys[i-2] )/(2*dx);
        //    printf( "[%i] /%g  dQ_ana  %g   dQ_num %g \n", i, line_dQi_num->ys[i-1]/line_dQi_ana->ys[i],  line_dQi_ana->ys[i],  line_dQi_num->ys[i-1] );
        //}


    }
}



/*
// ===================================================
///        test   Wave Function Overlap
// ===================================================
void test_WfOverlap( CLCFGO& solver, Plot2D& plot1 ){
    // ======= Test Orbital Wavefunction Overlap
    printf( " test_WfOverlap 1 \n" );
    int    nint = 10;
    double Lmax = 5.0;
    DataLine2D* line_ISgrid = new DataLine2D( nint, 0, Lmax/nint, 0xFFFF0000, "IS_grid" ); plot1.add(line_ISgrid );
    DataLine2D* line_ISana  = new DataLine2D( 100,  0, 0.1      , 0xFFFF8000, "IS_ana"  );  plot1.add(line_ISana  );
    {
    DEBUG_saveFile1="temp/wf0.xsf";
    DEBUG_saveFile2="temp/wf1.xsf";
    auto func1 = [&](GridShape& grid, double* f, double x ){                     solver.orb2grid( 0, grid, f ); };
    auto func2 = [&](GridShape& grid, double* f, double x ){ solver.epos[2].x=x; solver.orb2grid( 1, grid, f ); };
    gridNumIntegral( nint, 0.2, 6.0, Lmax, line_ISgrid->ys, func1, func2, true );
    }
    printf( " test_WfOverlap 2 \n" );
    for(int i=0; i<line_ISana->n; i++){
        solver.epos[2].x=line_ISana->xs[i];
        line_ISana->ys[i] = solver.evalOverlap( 0, 1 );
    }
}
*/

// ===================================================
///        test   Wave Function Overlap
// ===================================================

void test_Kinetic( CLCFGO& solver, Plot2D& plot1 ){
    // ======= Test Orbital Wavefunction Overlap
    int    nint = 20;
    double Lmax = 8.0;
    DataLine2D* line_ITgrid = new DataLine2D( nint, 0, Lmax/nint, 0xFFFF0000, "IT_grid" ); plot1.add(line_ITgrid );
    DataLine2D* line_ITana  = new DataLine2D( 100,  0, 0.1      , 0xFFFF8000, "IT_ana"  ); plot1.add(line_ITana  );
    {
    DEBUG_saveFile1="temp/wf0.xsf";
    DEBUG_saveFile2="temp/Lwf1.xsf";
    auto func1 = [&](GridShape& grid, double* f, double x ){ solver.epos[0].x=x; solver.orb2grid( 0, grid, f );      };
    auto func2 = [&](GridShape& grid, double* f, double x ){
        double* tmp = new double[grid.n.totprod()];
        solver.epos[0].x=x;
        solver.orb2grid( 0, grid, tmp );
        grid.Laplace   ( tmp, f );
        delete [] tmp;
    };
    gridNumIntegral( nint, 0.2, 8.0, Lmax, line_ITgrid->ys, func1, func2, true );
    }
    Vec3d dip;
    solver.bNormalize = false;
    for(int i=0; i<line_ITana->n; i++){
        solver.epos[0].x=line_ITana->xs[i];
        line_ITana->ys[i] = solver.projectOrb( 0, dip );
    }
    printf( "Ek[0] ana %g num %g / %g \n", line_ITana->ys[0], line_ITgrid->ys[0], line_ITana->ys[0]/line_ITgrid->ys[0] );
    printf( "KineticIntegral(0) Grid %g Ana %g ratio %g /%g \n", line_ITgrid->ys[0], line_ITana->ys[0],  line_ITgrid->ys[0]/line_ITana->ys[0],  line_ITana->ys[0]/line_ITgrid->ys[0]  );
}

// ===================================================
///        test   Density Projection
// ===================================================

void test_ProjectDensity( CLCFGO& solver, Plot2D& plot1 ){
    GridShape grid;
    int mpow = grid.init( 6.0, 0.1, true );
    int ng = grid.n.totprod();
    std::vector<double> grho1(ng,0.0);
    std::vector<double> gorb1(ng,0.0);
    std::vector<double> dgrho(ng,0.0);
    for(int i=0; i<solver.nBas; i++){ solver.ecoef[i]=0; solver.epos[i]=Vec3dZero;  }
    solver.esize[0] = 1.2; solver.ecoef[0] =  1.0; solver.epos[0]=Vec3dZero;
    //solver.esize[1] = 0.7; solver.ecoef[1] = -0.7; solver.epos[1]=(Vec3d){1.0,0.0,0.0};
    Vec3d dip;
    //solver.bNormalize = false;
    //solver.projectOrb( 0, dip );
    solver.projectOrb( 0, dip );
    solver.orb2grid  ( 0, grid, gorb1.data() );
    //solver.rho2grid  ( 0, grid, grho1.data() );
    solver.hartree2grid  ( 0, grid, grho1.data() );
    double sumRho = 0;
    double sumWf2 = 0;
    for(int i=0;i<ng;i++){
        double rho = grho1[i];
        double wf  = gorb1[i];
        wf=wf*wf;
        gorb1[i]=wf;
        dgrho[i]=rho-wf;
        sumRho += rho;
        sumWf2 += wf;
    }
    double dV = grid.voxelVolume();
    printf      ( "DEBUG |rho| %g |wf^2| %g \n", sumRho*dV, sumWf2*dV );
    float DEBUG_sc = 0.1;
    int ixy = grid.n.x*grid.n.y*grid.n.z/2 + grid.n.x*grid.n.y/2;
    DataLine2D* line_wf2  = new DataLine2D( grid.n.x, 0, 0.1      , 0xFF0080FF, "wf2"  ); plot1.add(line_wf2  );
    DataLine2D* line_rho  = new DataLine2D( grid.n.x, 0, 0.1      , 0xFFFF00FF, "rho"  ); plot1.add(line_rho  );
    for(int i=0; i<grid.n.x; i++){ line_wf2->ys[i]=gorb1[ixy+i]*DEBUG_sc; line_rho->ys[i]=grho1[ixy+i]*DEBUG_sc; }
    grid.saveXSF( "temp/rho_orb.xsf", gorb1.data() );
    grid.saveXSF( "temp/rho_aux.xsf", grho1.data() );
    grid.saveXSF( "temp/dgrho.xsf"  , dgrho.data() );
}

/*
// ===================================================
///        test   Density Overlap
// ===================================================
void test_DensityOverlap( CLCFGO& solver, Plot2D& plot1 ){
     // ======= Test Orbital Density Overlap
    int    nint = 10;
    double Lmax = 5.0;
    DataLine2D* line_IrhoGrid = new DataLine2D( nint, 0, Lmax/nint, 0xFF0000FF, "Irho_grid" ); plot1.add(line_IrhoGrid );
    DataLine2D* line_IrhoAna  = new DataLine2D( 100, 0, 0.1      , 0xFF0080FF, "Irho_ana"  ); plot1.add(line_IrhoAna  );
    //DataLine2D* line_IrhoWf   = new DataLine2D( 100, 0, 0.1      , 0xFFFF00FF, "Irho_Wf"   ); plot1.add(line_IrhoWf  );
    //for(int i=0; i<line_ISana->n; i++){ line_IrhoWf ->ys[i] = sq(line_ISana->ys[i]); }
    for(int i=0; i<line_IrhoAna->n; i++){
        solver.epos[2].x=line_IrhoAna->xs[i];
        solver.projectOrbs( true );
        line_IrhoAna->ys[i] = solver.DensOverlapOrbPair( 0, 1 );
    }
    { // ---- Test On Grid
    DEBUG_saveFile1="temp/rho0.xsf";
    DEBUG_saveFile2="temp/rho1.xsf";
    auto func1 = [&](GridShape& grid, double* f, double x ){
        //solver.epos[0].x=x;
        //solver.projectOrbs();
        solver.orb2grid( 0, grid, f );
        int ntot=grid.n.totprod();
        for(int i=0; i<ntot; i++){ double fi=f[i]; f[i]=fi*fi; }
    };
    auto func2 = [&](GridShape& grid, double* f, double x ){
        solver.epos[2].x=x;
        //solver.projectOrbs();
        int ntot=grid.n.totprod();
        solver.orb2grid( 1, grid, f );
        for(int i=0; i<ntot; i++){ double fi=f[i]; f[i]=fi*fi; }
    };
    gridNumIntegral( nint, 0.2, 6.0, Lmax, line_IrhoGrid->ys, func1, func2, true );
    printf( "DEBUG rho*rho grid %g | | %g \n", line_IrhoGrid->ys[0], line_IrhoAna->ys[0], Gauss::norm3Ds(1) );
    }
}
*/

// =========================================================================
///        test   Electrostatics   ( density * Hartree-Potential overlap )
// =========================================================================

void test_ElectroStatics( CLCFGO& solver, Plot2D& plot1 ){
    // ======= Test Orbital Density Overlap
    int    nint = 60;
    double Lmax = 6.0;
    DataLine2D* line_IElGrid = new DataLine2D( nint, 0, Lmax/nint, 0xFF0000FF, "IEl_grid" ); plot1.add(line_IElGrid );
    DataLine2D* line_IElAna  = new DataLine2D( 60, 0, 0.1       , 0xFF0080FF, "IEl_ana"  ); plot1.add(line_IElAna  );
    //DataLine2D* line_IrhoWf   = new DataLine2D( 100, 0, 0.1      , 0xFFFF00FF, "Irho_Wf"   ); plot1.add(line_IrhoWf  );

    solver.epos[0]=Vec3dZero;
    solver.epos[1]=Vec3dZero;

    for(int i=0; i<line_IElAna->n; i++){
        solver.epos[0].x=line_IElAna->xs[i];
        solver.projectOrbs( true );
        //line_IrhoAna->ys[i] = solver.DensOverlapOrbPair( 0, 1 );
        line_IElAna->ys[i] = solver.CoulombOrbPair( 0, 1 );
    }
    { // ---- Test On Grid
    DEBUG_saveFile1="temp/VH0.xsf";
    DEBUG_saveFile2="temp/rho1.xsf";
    auto func1 = [&](GridShape& grid, double* f, double x ){     //    ---- Potencial of Orb_1
        int ntot=grid.n.totprod();
        solver.epos[0].x=0;
        solver.projectOrbs( true );
        //solver.orb2grid( 0, grid, f );
        //solver.rho2grid( 0, grid, f );
        solver.hartree2grid( 0, grid, f );
        //for(int i=0; i<ntot; i++){ double fi=f[i]; f[i]=fi*fi; }
    };
    auto func2 = [&](GridShape& grid, double* f, double x ){     //    ---- Density of Orb_2
        int ntot=grid.n.totprod();
        solver.epos[1].x=x;
        //solver.projectOrbs();
        solver.orb2grid( 1, grid, f );
        for(int i=0; i<ntot; i++){ double fi=f[i]; f[i]=fi*fi; }
    };
    gridNumIntegral( nint, 0.2, 6.0, Lmax, line_IElGrid->ys, func1, func2, true );
    for(int i=0; i<line_IElAna->n; i++){ printf( "%i %g %g  %g %g %g \n", i, line_IElAna->xs[i], line_IElAna->ys[i], line_IElGrid->xs[i], line_IElGrid->ys[i],  line_IElGrid->ys[i]*2 ); };
    printf( "DEBUG rho*rho grid %g | | %g \n", line_IElGrid->ys[0], line_IElAna->ys[0], Gauss::norm3Ds(1) );
    }
}



// =========================================================================
///        test   Electrostatics   ( density * Hartree-Potential overlap )
// =========================================================================

void test_ElectroStaticsBrute( CLCFGO& solver, Plot2D& plot1 ){
    // ======= Test Orbital Density Overlap
    int    nint = 10;
    double Lmax = 6.0;
    DataLine2D* line_IElGrid = new DataLine2D( nint, 0, Lmax/nint, 0xFF0000FF, "IEl_grid" ); plot1.add(line_IElGrid );
    DataLine2D* line_IElAna  = new DataLine2D( 100, 0, 0.1       , 0xFF0080FF, "IEl_ana"  ); plot1.add(line_IElAna  );
    //DataLine2D* line_IrhoWf   = new DataLine2D( 100, 0, 0.1      , 0xFFFF00FF, "Irho_Wf"   ); plot1.add(line_IrhoWf  );

    // --- define grid
    GridShape gsh;
    double Rmax  = 4.0;
    double gStep = 0.15;
    int ng = (2*Rmax)/gStep;
    gsh.cell = (Mat3d){ (2*Rmax),0.0,0.0,  0.0,(2*Rmax),0.0,  0.0,0.0,(2*Rmax) };
    gsh.n    = {ng,ng,ng};
    gsh.pos0 = (Vec3d){-Rmax,-Rmax,-Rmax};
    gsh.updateCell();
    double  dV = gsh.voxelVolume();

    solver.epos[0]=Vec3dZero;
    solver.epos[1]=Vec3dZero;

    double * grid1 = new double[ gsh.n.totprod() ];
    double * grid2 = new double[ gsh.n.totprod() ];
    solver.orb2grid( 0, gsh, grid1 );
    solver.orb2grid( 0, gsh, grid2 );

    double Q1=0,Q2=0;
    for(int i=0; i<gsh.n.totprod(); i++){
        double dq;
        dq = grid1[i]; dq*=dq; grid1[i]=dq; Q1+=dq;
        dq = grid2[i]; dq*=dq; grid2[i]=dq; Q2+=dq;
    }

    printf( "gsh n%i(%i,%i,%i) Giga-ops %g dV %g Q1 %g Q2 %g const_El_eVA %g \n", gsh.n.totprod(), gsh.n.x, gsh.n.y, gsh.n.z, sq( (double)gsh.n.totprod() )*1e-9, dV, Q1, Q2, const_El_eVA );

    for(int i=0; i<line_IElAna->n; i++){
        double r = line_IElAna->xs[i];
        //printf( "test_ElectroStaticsBrute()[%i] r=%g ", i, r );
        solver.epos[0].x=r;

        // --- Analytical
        solver.projectOrbs( true );
        line_IElAna->ys[i] = solver.CoulombOrbPair( 0, 1 );

        // Numerical
        //line_IElGrid->ys[i] = coulombGrid_brute( gsh.n, gsh.n, Vec3dZero, Vec3d{r,0,0}, gsh.dCell, gsh.dCell, grid1, grid2, 0, gsh.dCell.a.norm() )*(dV*dV) * const_El_eVA;
        line_IElGrid->ys[i]=0;

        //printf( " Ana %g Grid %g \n", line_IElAna->ys[i], line_IElGrid->ys[i] );

        printf( "%i %g %g %g \n", i, r, line_IElAna->ys[i], line_IElGrid->ys[i] );
    }

    delete [] grid1;
    delete [] grid2;
}



#endif
