
// read also:
// http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToAngle/

#ifndef  PolyLinear1d_h
#define  PolyLinear1d_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "fastmath.h"

class PolyLinear1d{
	public:
	int n;
	double * xs = NULL;
	double * ys = NULL;

    double integrate_ibounds( int imin, int imax ){
        double Isum = 0.0d;
        double ox  = xs[imin];
        double oy  = ys[imin];
        for( int i = imin+1; i<imax; i++ ){
            double x  = xs[i];
            double y  = ys[i];
            Isum     += ( y + oy ) * ( x - ox );
            //printf( " %i %f %f %f \n ", i, Isum, x , y );
            ox        = x;
            oy        = y;
        }
        return Isum * 0.5;
    }

	double integrate( double xmin, double xmax ){
		int i = 0;
		// find lower bound
		while( xs[i] < xmin ){ i++; }
		double ox,oy,x,y,Isum=0.0d;
		//printf( " lower %i %f %f \n", i, xs[i], xmin );
		// evaluate integral for first interval
		if( i>0 ){ // this will happen only if xmin is inside polyline range
            ox = xs[i-1]; oy = ys[i-1];
			x  = xs[i  ];  y = ys[i  ];
			double dydx = ( y - oy ) / ( x - ox );
			if( x < xmax ){
                double dx = x - xmin;
				Isum += ( y - 0.5 * dydx * dx ) * dx;
				//printf( " _| %i %f (%3.3f,%3.3f) \n", i, i, Isum,  xmin, x  );
				ox = x; oy = y;
                i++;
			}else{
				double dx = xmax - xmin;
				Isum += ( oy + dydx * ( xmin - ox ) + 0.5 * dydx * dx ) * dx;
				//printf( " _  %i %f (%3.3f,%3.3f) \n", i, Isum,  xmin, xmax );
				return Isum;
			}
		}else{
            ox = xs[i];
            oy = ys[i];
            i++;
		}
		while( true ){
			x  = xs[i]; y = ys[i];
			if( x > xmax ) break;
            Isum += 0.5*( y + oy )*( x - ox );
            //printf( "|_| %i %f  (%3.3f,%3.3f) \n", i, Isum, ox, x );
            ox = x; oy = y;
			i++;
		}
        double dydx = ( y - oy ) / ( x - ox );
        double dx   = xmax - ox;
        Isum += ( oy + 0.5 * dydx * dx ) * dx;
        //printf( "|_  %i %f  (%3.3f,%3.3f) \n", i, Isum,  ox, xmax );
		return Isum;
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
			if( dI < Irest ){
				Isum += dI;
				//printf( "|_|  %i %f  (%3.3f,%3.3f) \n", i, Isum,  ox, x );
				ox = x;
				oy = y;
			}else{
			    //printf( " x_of_integral %i %f %f \n", i, Irest, dI );
                double dxdy  = ( y - oy ) / dx;
                if( fabs(dxdy) < 1e-8 ){
                    dx = Irest / oy;
                    //return ox + Irest / oy;
                }else{
                    double x1,x2;
                    quadratic_roots( 0.5d*dxdy, oy, -Irest, x1, x2 );
                    if( dxdy > 0 ){ dx = _max( x1, x2 ); } else { dx = _min( x1, x2 ); };
                    //printf( " x_of_integral quadratic_roots %f %f %f %f %f \n", dxdy, oy, Irest, x1, x2  );
                }
                //printf( " x1 %f x2 %f   %f %f \n", x1, x2, ox, ox + _max( x1, x2 ) );
				return ox + dx;
			}
		}
		return NAN;
	}

    void detach(){ xs = NULL; ys = NULL; };

    PolyLinear1d(){};
    PolyLinear1d( int n_, double * xs_, double * ys_ ){ n=n_; xs=xs_; ys=ys_; };

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


