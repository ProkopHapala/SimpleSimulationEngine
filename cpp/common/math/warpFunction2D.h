
#ifndef  warpFunction2D_h
#define  warpFunction2D_h

//#include "fastmath.h"
#include "functions.h"

int NWRAP = 4;

/*
typedef double (*Function2d)( double x, double y );
typedef double (*DiffFunc2d)( double x, double y, double& dfdx, double& dfdy );
typedef void   (*Warp2d    )(  double& Tx,   double& Ty,   double& dTxdx, double& dTxdy, double& dTydx, double& dTydy );


// ===================== 1D Functions primitives

double x1period( double x, double& dfdx ){
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

double x2period( double x, double& dfdx ){
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

// ===================== 2D Functions primitives

double harmonic( double x, double y, double& dfdx, double& dfdy ){
	dfdx = x;
	dfdy = y;
	return x*x + y*y;
}

double lorenz( double x, double y, double& dfdx, double& dfdy ){
	double D = 1/(1+x*x + y*y);
	double D2 = D*D*2;
	dfdx = -D2*x;
	dfdy = -D2*y;
	return D;
}

double mlorenz( double x, double y, double& dfdx, double& dfdy ){
	double f = 1 - lorenz( x, y, dfdx, dfdy );
	dfdx = -dfdx; dfdy = -dfdy;
	return f;
}

*/

// ===================== Warping


#define _WARP_DERIVATIVS_2D_ \
	double dFxdx = dFxdTx * dTxdx  +  dFxdTy * dTydx;	\
	double dFxdy = dFxdTx * dTxdy  +  dFxdTy * dTydy;	\
	double dFydx = dFydTx * dTxdx  +  dFydTy * dTydx;	\
	double dFydy = dFydTx * dTxdy  +  dFydTy * dTydy;	\
	dTxdx = dFxdx;  dTxdy = dFxdy;  					\
	dTydx = dFydx;  dTydy = dFydy;						\

void warp_sin( double& x, double& y,   double& dTxdx, double& dTxdy, double& dTydx, double& dTydy ){
	double strength=0.8;
	double sx = sin(x);
	double sy = sin(y);
	double dFxdTx = 1;						double dFxdTy = cos(y)*strength;
	double dFydTx = cos(x)*strength;		double dFydTy = 1;
	x += sy*strength;
	y += sx*strength;
	_WARP_DERIVATIVS_2D_;
}

void warp_x2period( double& x, double& y,   double& dTxdx, double& dTxdy, double& dTydx, double& dTydy ){
	double strength=0.5;
	double dFxdTx = 1, dFxdTy;
	double dFydTx, dFydTy = 1;
	double sy = x2period(y, dFxdTy );   dFxdTy*=strength;
	double sx = x2period(x, dFydTx );	 dFydTx*=strength;
	x += sy*strength;
	y += sx*strength;
	_WARP_DERIVATIVS_2D_;
}

void warp_x1period( double& x, double& y,   double& dTxdx, double& dTxdy, double& dTydx, double& dTydy ){
	double strength=0.5;
	double dFxdTx = 1, dFxdTy;
	double dFydTx, dFydTy = 1;
	double sy = x1period(y, dFxdTy );   dFxdTy*=strength;
	double sx = x1period(x, dFydTx );	 dFydTx*=strength;
	x += sy*strength;
	y += sx*strength;
	_WARP_DERIVATIVS_2D_;
}

double warped_function( double x, double y, double& dfdx, double& dfdy ){

	double fscale  = 2.0;
	double ifscale = 1.0/fscale;
	x *= fscale;	y *= fscale;
	double dTxdx = 1,dTxdy = 0;
	double dTydx = 0,dTydy = 1;
	for(int i=0; i<NWRAP; i++){
		//warp_sin( x, y,   dTxdx, dTxdy, dTydx, dTydy );
		warp_x2period( x, y,   dTxdx, dTxdy, dTydx, dTydy );
		//warp_x1period( x, y,   dTxdx, dTxdy, dTydx, dTydy );
		//x+=0.1-y*0.25;
		//y+=0.1+x*0.25;
		x *= 1.1;
	}
	x *= ifscale;	y *= ifscale;
	double dfdTx, dfdTy;
	//double f = harmonic( x, y, dfdTx, dfdTy );
	double f = mlorenz( x, y, dfdTx, dfdTy );
	dfdx = ( dfdTx * dTxdx + dfdTy * dTydx )*ifscale;
	dfdy = ( dfdTx * dTxdy + dfdTy * dTydy )*ifscale;

	//double f = x2period( x, dfdx )*0.5 + 0.5;
	//double f = x1period( x, dfdx )*0.5 + 0.5;

	return f;
}

double warped_function_noDiff( double x, double y ){
	double dfdx,dfdy;
	//return harmonic( x, y, dfdx, dfdy );
	return warped_function( x, y, dfdx, dfdy );
}

#endif

