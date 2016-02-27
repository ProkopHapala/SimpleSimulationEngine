
#ifndef buoyancy_h
#define buoyancy_h

#include "fastmath.h"
#include "Vec2.h"
//#include "VecN.h"
#include "geom2D.h"
#include "Convex2d.h"
//#include "PolyLinear1d.h"

/*
 centre of mass:   yc(x) = 0.5*(yl(x) + yr(x))
 volume:           yv(y) =      yl(x) - yr(x)
 V = Integral_dc yv(x)
 M = Integral_dx yc(x)*yv(x)
 yc(x) = 0.5* ( x*(yc-oyc)/(x-ox) + oyc ) = 0.5*(x*ac + oyc)   where   ac = (yc-oyc)/(x-ox)   =   ( yl + yr - oyl - oyr )/dx
 yv(y) = x*(yv-oyv)/(x-ox) + oyv          = x*av + oyc         where   av = (yv-oyv)/(x-ox)   =   ( yl - yr - oyl + oyr )/dx
 V =       x    *oyv  +  x^2/2 *       av
 M = 0.5*( x*oyc*oyv  +  x^2/2 * ( oyc*av + ac*oyv )  +  x^3/3 * ac*av  )
*/

double integrate_moment( int n, double * xs, double * yLs, double * yRs, double displacement, double& watterline ){
    double Vsum = 0.0d;
    double Msum = 0.0d;
    double ox   = xs [0];
    double oyl  = yLs[0];
    double oyr  = yRs[0];
    double oyc  = oyl + oyr;
    double oyv  = oyl - oyr;
    for( int i=1; i<n; i++ ){
        double x  = xs [i];
        double yl = yLs[i];
        double yr = yRs[i];
        double dx = x  - ox;
        double yc = yl + yr;
        double yv = yl - yr;
        double ac = yc - oyc;
        double av = yv - oyv;
        double dV = dx*( oyv + av*0.5 );
        double Vrest = displacement - Vsum;
        if( dV < Vrest ) {
            //Msum += dx*( oyc*oyv + dx*( 0.5d*( oyc*av/dx + ac*oyv/dx ) + dx*0.33333333333d*ac*av/dx/dx ) );
            Msum += dx*( oyc*oyv + 0.5d*( oyc*av + ac*oyv ) + 0.33333333333d*ac*av );
            Vsum += dV;
            //printf( "|_| %i %f  (%3.3f,%3.3f) \n", i, Vsum, ox, x );
            ox=x; oyc=yc; oyv=yv;
        }else{
            //printf( " integrate_moment %i %f %f \n", i, Vrest, dV );
            if( dx < 1e-8 ){
                watterline = ox;
                return Msum * 0.5d;
            };
            double idx = 1/dx;
            ac *= idx;
            av *= idx;
            if( fabs( av ) < 1e-8 ){
                dx = Vrest / oyv;
                //printf( " linear dx %f   %f %f %f \n", dx,   Vrest, av, oyv );
                Msum += dx*( oyc*oyv + dx*0.5*ac*oyv );
                //watterline = ox + dx;
                //return Msum * 0.5d;
            }else{
                double x1,x2;
                dx =  quadratic_roots( 0.5d*av, oyv, -Vrest, x1, x2 );
                //dx = _max( x1, x2 );
                if( av > 0 ){ dx = _max( x1, x2 ); } else { dx = _min( x1, x2 ); };
                //printf( " quadrt dx %f   %f %f %f \n", dx,  Vrest, av, oyv );
                //printf( " integrate_moment quadratic_roots %f %f %f %f %f \n", av, oyv, Vrest, x1, x2  );
                Msum += dx*( oyc*oyv + dx*( 0.5*( oyc*av + ac*oyv ) + dx*0.33333333333*ac*av ) );
                //watterline = ox + dx;
                //return Msum * 0.5d;
            }
            //Msum += dx*( oyc*oyv + dx*( 0.5*( oyc*av + ac*oyv ) + dx*0.33333333333*ac*av ) );
            watterline = ox + dx;
            return Msum * 0.5d;
        }
    }
    return NAN;
}

double buoy_moment_2D( const Convex2d& hull, const Vec2d& dir, const Vec2d& cog, double displacement ){
    double yLs [hull.n];
    double yRs [hull.n];
    double ys  [hull.n];
    double xs  [hull.n];
	hull.projectToLine( dir, xs, yLs, yRs );
	double watterline;
	double moment = integrate_moment( hull.n, xs, yLs, yRs, displacement, watterline );
	double ycog = dir.dot_perp( cog );
    return moment/displacement - ycog;
}

#endif





