#ifndef  functions_h
#define  functions_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

//#include "fastmath.h"


int N_MANDELBROT = 8;

// ==== Types

typedef double (*Function1d)( double x);
typedef double (*DiffFunc1d)( double x, double& dx );

typedef double (*Function2d)( double x, double y);
typedef double (*DiffFunc2d)( double x, double y, double& dfdx, double& dfdy );

typedef double (*FunctionNd)( int n, double * xs );
typedef double (*DiffFuncNd)( int n, double * xs, double * dfs );

// ==== functions 1D

inline double sigmoideAbs( double x, double& dfdx ){
	double D  = 1/( 1 + fabs(x) );
	dfdx      = D*D;
	return  x*D;
}

inline double sigmoideSqrt( double x, double& dfdx ){
	double D2 = 1/(1 + x*x);
	double D  = sqrt( D2 );
	dfdx      = D*D2;
	return  x*D;
}

inline double lorenz( double x, double& dfdx ){
	double f = 1/(1 + x*x);
	dfdx     = 2*x*f*f;
	return  f;
}

inline double x1period( double x, double& dfdx ){
	int ix = (int)(x+1000.0d)-1000;
	double dx = x-ix;
	if( ix&0x1 ){
		dfdx =    2.0d;
		return 2.0d*dx-1.0d;
	}else{
		dfdx =   -2.0d;
		return  1.0d-2.0d*dx;
	}
}

inline double x2period( double x, double& dfdx ){
	int ix = (int)(x+1000.0)-1000;
	double dx = x-ix-0.5d;
	if( ix&0x1 ){
		dfdx =  -8*dx;
		return 1-dx*dx*4;
	}else{
		dfdx =   8*dx;
		return dx*dx*4-1;
	}
}

// ==== functions 2D

/*
double warp_x2period( double x, double y, double& dfdx, double& dfdy ){
	double cwarp = 0.5;
	for(int i=0; i<5; i++){
		double dsx,dsy;
		double sx = x2period( x, dsx );
		double sy = x2period( y, dsy );
 		x = x + cwarp*sx;
		y = y + cwarp*sy;
	}
}

*/

inline double lorenz( double x, double y, double& dfdx, double& dfdy ){
	double D = 1/(1+x*x + y*y);
	double D2 = D*D*2;
	dfdx = -D2*x;
	dfdy = -D2*y;
	return D;
}

inline double mlorenz( double x, double y, double& dfdx, double& dfdy ){
	double f = 1 - lorenz( x, y, dfdx, dfdy );
	dfdx = -dfdx; dfdy = -dfdy;
	return f;
}

inline double harmonic( double x, double y ){
	return x*x + y*y;
}

inline double harmonic( double x, double y, double& dfdx, double& dfdy ){
	dfdx = x;
	dfdy = y;
	return x*x + y*y;
}

inline double rosenbrok( double x, double y ){
	double f = (x*x - y);
	return f*f + x*x*0.1;
}

inline double sinValey( double x, double y ){
	double f = ( sin(x)  - y);
	return f*f + x*x*0.1;
}

inline double cosValey( double x, double y ){
	double f = ( cos(4*x) - y);
	return (0.2*f*f + x*x*0.05)*0.2;
}

inline double spiral( double x, double y ){
	double  phi = atan2 (y,x);
	double  r2  = x*x + y*y;
	double  r   = sqrt(r2);
	return  (1+sin( 2*M_PI*r + phi ))*0.1     + 0.05*r2;
}

double mandelbort( double cX, double cY ){
	double x = 0, y = 0;
	int i=0;
	for (i=0; i<N_MANDELBROT; i++){
		double x_ =  x * x - y * y + cX;
        y         = 2.0 * x * y    + cY;
        x = x_;
	}
	double r = sqrt(x*x + y*y);
	//return 0.5*(sin(r)+1.0);
	return tanh(r);
}

#endif

