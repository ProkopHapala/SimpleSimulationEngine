
#ifndef buoyancy_h
#define buoyancy_h

#include "fastmath.h"
#include "Vec2.h"
#include "VecN.h"
#include "geom2D.h"
#include "Convex2d.h"
#include "PolyLinear1d.h"


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
            Msum  = dx*( oyc*oyv + 0.5d*( oyc*av + ac*oyv ) + 0.33333333333*ac*av );
            Vsum += dV;
            //printf( "|_| %i %f  (%3.3f,%3.3f) \n", i, Vsum, ox, x );
            ox=x; oyc=yc; oyv=yv;
        }else{
            //printf( " integrate_moment %i %f %f \n", i, Vrest, dV );
            double idx = 1/dx;
            ac *= idx;
            av *= idx;
            double x1,x2;
            dx =  quadratic_roots( 0.5d*av, oyv, -Vrest, x1, x2 );
            dx = _max( x1, x2 );
            //printf( " integrate_moment quadratic_roots %f %f %f %f %f \n", 0.5d*av, oyv, -Vrest, x1, x2  );
            Msum += dx*( oyc*oyv + dx*( ( 0.5d*oyc*av + ac*oyv ) + dx*0.33333333333*ac*av ) );
            watterline = ox + dx;
            return Msum * 0.5d;
        }
    }
    return NAN;
}

/*
// FIXME : more effective would be compute watterline together with moment
double integrate_moment( int n, double xs, double yLs, double yRs, double displacement ){
    double Vsum = 0.0d;
    double Msum = 0.0d;
    double ox =xs [i];
    double oyl=yLs[i];
    double oyr=yRs[i];
    for( int i=1; i<n; i++ ){
        double x   = xs [i];
        double yl  = yLs[i];
        double yr  = yRs[i];
        double dx  = x  - ox;
        double bc =           oyl + oyr;
        double bv =           oyl - oyr;
        double ac = yl + yr - oyl - oyr;
        double av = yl - yr - oyl + oyr;
        double dV = dx*( bv + av );
        double Vrest = displacement - Vsum;
        if( dV < Vrest ) {
            Msum  = dx*( bc*bv + 0.5d*( bc*av + ac*bv ) + 0.33333333333*ac*av );
            Vsum += dV;
            //printf( "|_| %i %f  (%3.3f,%3.3f) \n", i, Isum, ox, x );
            ox=x; oyl=yl; oyr=yr;
        }else{
            double idx = 1/dx;
            ac *= idx;
            av *= idx;
            double x1,x2;
            dx =  quadratic_roots( 0.5d*av, bv, -Vrest, x1, x2 );
            dx = _max( x1, x2 );
            Msum += dx*( bc*bv + dx*( ( 0.5d*bc*av + ac*bv ) + dx*0.33333333333*ac*av ) );
            return Msum * 0.5d;
        }
    }
    return NAN;
}
*/


double buoy_moment_2D( const Convex2d& hull, const Vec2d& dir, const Vec2d& cog, double displacement ){
    double yLs [hull.n];
    double yRs [hull.n];
    double ys  [hull.n];
    double xs  [hull.n];
	hull.projectToLine( dir, xs, yLs, yRs );
	double watterline;
	double moment = integrate_moment( hull.n, xs, yLs, yRs, displacement, watterline );
	//VecN::sub( np, yLs, yRs, ys );
	//PolyLinear1d pline( np, xs, ys );
	//double x_watterline = pline.x_of_integral( displacement );
	//double xcog = dir.dot_perp( cog );
	//pline.detach();
    return moment;
}

#endif





