
#ifndef  Multipoles_h
#define  Multipoles_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "fastmath.h"
#include "Vec3.h"

// https://en.wikipedia.org/wiki/Multipole_expansion
// https://en.wikipedia.org/wiki/Axial_multipole_moments
// https://en.wikipedia.org/wiki/Spherical_multipole_moments
// http://physics.stackexchange.com/questions/106564/clarification-of-multipole-expansion-for-a-point-charge


// fast multipole method
// https://github.com/davidson16807/fast-multipole-method/blob/master/fast-multipole-method.js
// https://github.com/barbagroup/pyfmm/tree/master/pyFMM

inline double Eelec( const Vec3d& dR, double qq ){
    double r = dR.norm();
    return qq/r;
}

double evalEelectrostatic( const Vec3d& p, double Q, int n, Vec3d * ps, double * Qs ){
    double E = 0;
    for( int i=0; i<n; i++){
        E += Eelec( ps[i]-p, Q*Qs[i] );
    }
    return E;
}

inline void getMultipole( const Vec3d& dR, double Q, int order, double * coefs ){
    coefs[0] += Q;
    if(order<1) return;
    coefs[1] += dR.x*Q;
    coefs[2] += dR.y*Q;
    coefs[3] += dR.z*Q;
    if(order<2) return;
    coefs[4] += dR.x*dR.x*Q;
    coefs[5] += dR.x*dR.y*Q;
    coefs[6] += dR.y*dR.y*Q;
    coefs[7] += dR.y*dR.z*Q;
    coefs[8] += dR.z*dR.z*Q;
    coefs[9] += dR.z*dR.x*Q;
}

void getMultiPole( const Vec3d& center, int n, Vec3d * ps, double * Qs, int order, double * coefs ){
    for( int i=0; i<10; i++ ) coefs[i]=0;
    for( int i=0; i<n; i++){
        getMultipole( ps[i]-center, Qs[i], order, coefs );
        /*
        Vec3d dR = ps[i] - center;
        double Q = Qs[i];
        coefs[0] += Q;
        if(order<1) continue;
        coefs[1] += dR.x*Q;
        coefs[2] += dR.y*Q;
        coefs[3] += dR.z*Q;
        if(order<2) continue;
        coefs[4] += dR.x*dR.x*Q;
        coefs[5] += dR.x*dR.y*Q;
        coefs[6] += dR.y*dR.y*Q;
        coefs[7] += dR.y*dR.z*Q;
        coefs[8] += dR.z*dR.z*Q;
        coefs[9] += dR.z*dR.x*Q;
        */
    }
}

double Emultipole( const Vec3d& dR, int order, double * coefs ){
    //double r   = dR.norm();
    //double ir  = 1 / r;
    //double ir2 = ir*ir;
    double ir2 = 1/dR.norm2();
    double E   = coefs[0];
    if( order>0 ) E += ir2    *( coefs[1]*dR.x + coefs[2]*dR.y + coefs[3]*dR.z );
    if( order>1 ) E += ir2*ir2*((coefs[4]*dR.x + coefs[5]*dR.y)*dR.x +
                                (coefs[6]*dR.y + coefs[7]*dR.z)*dR.y +
                                (coefs[8]*dR.z + coefs[9]*dR.x)*dR.z );
    return sqrt(ir2)*E;
}

#endif

