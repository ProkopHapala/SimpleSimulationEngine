
#ifndef SuperSonic2D_h
#define SuperSonic2D_h

/*

https://en.wikipedia.org/wiki/Oblique_shock
http://www.i-asem.org/publication_conf/anbre15/T2I.5.AS502_1358F1.pdf
http://www.tfd.chalmers.se/~nian/courses/compflow/notes/TME085_L05.pdf

*  tangential velocity component does not change across the shock
It is always possible to convert an oblique shock into a normal shock by a Galilean transformation.

tg(theta) = 2*cotg(beta) ( (M^2)*sin(beta)^2 -1 ) / ( M^2 (kapa + cos(2*beta) ) +2 )

sin(theta)/cos(theta) =    2*(cos(beta)/sin(beta))   ( (M^2)*sin(beta)^2 -1 ) / ( M^2 (kapa + cos(2*beta) ) +2 )

(st/ct)  = 2* (cb/sb)  * ( (M^2)*sb^2 -1 ) / ( M^2 (k + cb^2 - sb^2) ) +2 )
( M^2 (k + cb^2 - sb^2) ) + 2 )*(st/ct) =  2* (cb/sb)  ( (M^2)*sb^2 -1 )
(st/ct) = T

T*M2*k + 2*T +   T*M2*(cb^2 - sb^2) = 2*(cb/sb)  ( M2*sb^2 -1 )

 T*M2*k*s + 2*T*s  + T*M2*(c*c - s*s)*s - 2*M2*s*s*c + 2*c = 0
(T*M2*k   + 2*T)*s + T*M2*(c*c - s*s)*s - 2*M2*s*s*c + 2*c = 0

NASA Oblique Shock calculator
https://www.grc.nasa.gov/WWW/K-12/airplane/oblique.html
https://web.archive.org/web/20121021100737/http://www.aerostudents.com/files/aerodynamicsC/obliqueShockWaves.pdf

https://www.youtube.com/watch?v=JvrWmg8m5is

*/

#include "fastmath.h"
#include "Vec2.h"
//#include "Vec3.h"

#include "approximation.h"


namespace ObliqueShock{

inline double tgTheta( double c, double s, double M, double k ){
    double MM   = M*M;
    double cotg = c/s;
    double C    = MM*(k+1) + 2;
    return -cotg*(1 + (2-C)/( C - 2*MM*s*s ) );
}

inline double pressureRatio( double s, double M, double k ){
    /// See: https://en.wikipedia.org/wiki/Oblique_shock#Maximum_deflection_angle
    double M2s2 = M*s;    M2s2 *= M2s2;
    return 1 + (2*k/(k+1))*( M2s2 - 1 );
}

inline double densityRatio( double s, double M, double k ){
    /// See: https://en.wikipedia.org/wiki/Oblique_shock#Maximum_deflection_angle
    double M2s2 = M*s;    M2s2 *= M2s2;
    return (k+1)*M2s2/( (k-1)*M2s2 + 2);
}

inline double temperatureRatio( double s, double M, double k ){
    /// See: https://en.wikipedia.org/wiki/Oblique_shock#Maximum_deflection_angle
    double M2s2 = M*s;    M2s2 *= M2s2;
    double k1 = k+1;
    double pRatio   = 1 + (2*k/k1)*( M2s2 - 1 );
    double rhoRatio = k1*M2s2/( (k-1)*M2s2 + 2);
    return pRatio*rhoRatio;
    //return
}

double evalM2( double s, double sbt, double M, double k ){
    // sbt .... sin(beta+theta)  = s(a-b) = sa*cb - sb*ca
    /// See: https://en.wikipedia.org/wiki/Oblique_shock#Maximum_deflection_angle
    double M2s2 = M*s;  M2s2 *= M2s2;
    double kmod = (k-1)*0.5;
    return sqrt( ( 1 + kmod*M2s2 )/( k*M2s2 - kmod ) )/sbt;
}

void makeInterpTable( double sinThetaMin, double M, double k, int m, double* poly ){
    const int nmax=100;
    double sbs[nmax]; // sin(beta)
    double sts[nmax]; // sin(theta)
    // tg = s/c = s/sqrt( 1-s*s )
    double dsb = M_PI/nmax;
    double ot=-100;
    int i0=-1,i1;
    for(int i=0; i<nmax; i++ ){
        double s   = i*dsb;
        sbs[i]=s;
        double t   = tgTheta( sqrt(1-(s*s)), s, M, k );
        if(t<ot){i1=i; break; }
        double t2  = t*t;
        double st = sqrt(t2/(1+t2)); // tan(theta) -> sin(theta)
        if(t<0)st=-st;
        sts[i]=st;
        if((i0<0)&&(t>sinThetaMin)){ i0=i; };
    }
    Approx::polyFit( i1-i0, m, sts+i0, sbs+i0, poly );
}


class Solver{

    double Min;    // incident mach number
    Vec2d  normal; // direction of shock

    double sinThetaMin;
    double dM,invdM;
    int nM=0, npoly=0;
    double* table=0;

    void makeTables( int npoly_, int nM_, double dM_, double kappa, double sinThetaMin_ ){
        int npoly    = npoly_;
        sinThetaMin  = sinThetaMin_;
        dM           = dM_;
        invdM        = 1/dM;
        table        = new double[npoly*nM];
        for( int i=0; i<nM; i++ ){
            double M = 1 + i*dM;
            makeInterpTable( sinThetaMin, M, kappa, npoly, table+i*npoly );
        }
    }

    double getSinBeta( double sinTheta, double M ){
        double kM = (M-1)*invdM;
        int    iM = (int)kM;
        double cM = kM-iM;
        double y0 = Approx::pval( sinTheta, npoly, table+iM   );
        double y1 = Approx::pval( sinTheta, npoly, table+iM+1 );
        return y0*(1-cM) + y1*cM;
    }

    ~Solver(){ delete [] table; };

};

}


struct SScell{
    Vec2d  v; // velocity x,y
    double p; // pressure
};

class SuperSonic2D{ public:

    Vec2i ns;
    double** lines;
    double*  line;
    Vec2d* v;


    void eval(){
        // push neighboring samples asside to fullfill equation of state and continuity equation
        //pV = nRT
        for(int i=0; i<ns.y;  i++){
        }
    }



};

#endif

