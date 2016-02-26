
// read also:
// http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToAngle/

#ifndef  PolyLinear1D_h
#define  PolyLinear1D_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "fastmath.h"

class PolyLinear1d{
	public:
	int n;
	int * xs = NULL;
	int * ys = NULL;

	double integrate( double xmin, double xmax ){
		int i = 0;
		// find lower bound
		while( xs[i] < xmin ){ i++; }
		double ox,oy,Isum;
		// evaluate integral for first interval
		if( i>0 ){
			ox = xs[i-1];
			oy = ys[i-1];
			double x    = xs[i];
			double y    = xs[i];
			double dx   = x - ox; 
			double dydx = ( y - oy ) / dx;
			if( x < xmax ){
				Isum = ( y - 0.5 * dydx * dx ) * dx;
				ox   = x;
				oy   = y;
			}else{
				dx   = xmax - xmin;
				Isum = ( oy + dydx * ( xmin - ox ) + 0.5 * dydx * dx ) * dx;
				return Isum;
			}
		}
		while( true ){
			double x     = xs[i];
			double y     = ys[i];
			double dx    = x - ox;
			if( x < xmax ){
				double yav  = 0.5*( y + oy );
				double dI   = yav*dx;
			}else{
				double dydx = ( y - oy ) / dx;
				double dx   = xmax - ox;
				Isum = ( oy + 0.5 * dydx * dx ) * dx;
				return Isum;
			}
			i++;
		}
		/*
		for( int i=1; i<n; i++ ){
			double x     = xs[i];
			double y     = ys[i];
			double yav   = 0.5*( y + oy );
			double dx    = x - ox;
			double dI    = yav*dx;
			double Irest = Iend - Isum;  
			if( dI > Irest ){
				return ox + dx * Irest / dI;  
			}else{
				Isum += dI;
				ox = x;
				oy = y;
			}
		}
		*/
	}

	double x_of_integral( double Iend ){
		double oy = ys[0];
		double ox = xs[0];
		double Isum = 0.0;
		for( int i=1; i<n; i++ ){
			double x     = xs[i];
			double y     = ys[i];
			double yav   = 0.5*( y + oy );
			double dx    = x - ox;
			double dI    = yav*dx;
			double Irest = Iend - Isum;  
			if( dI > Irest ){
				return ox + dx * Irest / dI;  
			}else{
				Isum += dI;
				ox = x;
				oy = y;
			}
		}
	}


    PolyLinear1d(){};

	PolyLinear1d( int n_ ){
        n       = n_;
		xs = new double[ n ];
		ys = new double[ n ];
	}

	~PolyLinear1d( ){
		if( xs != NULL ) delete xs;
		if( ys != NULL ) delete ys;
	}

};




#endif


