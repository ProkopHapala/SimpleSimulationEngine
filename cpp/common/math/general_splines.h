
#ifndef  radial_splines_h
#define  radial_splines_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include <Vec3.h>

inline double splineR2_0246( double r2, double c0, double c1, double c2, double c3 ){
	return c0 + r2*(c1+r2*(c2+r2*c3));
}

inline double get_splineR2_0246( double r2, double * coefs ){
	double ir2 = 1/r2;
	double  r4 = r2*r2;
	int i0     = (ir2*invstep)<<2;
	return coefs[o0] + coefs[o0+1]*r2 + coefs[o0+2]*r4 +  coefs[o0+3]*ir2;
}

#endif



