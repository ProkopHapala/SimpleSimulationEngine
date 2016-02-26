
#ifndef  Convex2d_h
#define  Convex2d_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"

#include "arrayAlgs.h" // FIXME move this to .cpp 

/////////////////////////
//   CLASS :   Convex2d
//////////////////////////

class Convex2d{
	public:
	int id;
	int n;              // number of lines or corners ( it is the same )
	Vec2d  * corners = NULL;
	Line2d * lines   = NULL;

	// ======= function declarations

    bool cut       ( const Line2d& line, Convex2d& left, Convex2d& right );
	void fromPoints( int np, Vec2d * points, int * permut, double * vals );
	void fromPoints( int np, Vec2d * points );
	void make_corners( );

	// ======= inline functions

	inline void update_lines( ){
		int nm1 = n - 1;
		for ( int i=0; i<nm1; i++ ){
			//lines[i].set( corners[i], corners[i+1] );
			lines[i].set( corners[i+1], corners[i] );
		}
		//lines[nm1].set( corners[nm1], corners[0] );
		lines[nm1].set( corners[0], corners[nm1] );
	}

	inline bool pointIn( const Vec2d& p ) const{
		return pointInConvexPolygon<double>( p, (double*)lines, n );
	}

    inline void boundingBox( Rect2d * rect ) const{
        //double a,b;
        double xmin,xmax,ymin,ymax;
        xmin=xmax=corners[0].x;
        ymin=ymax=corners[0].y;
        for ( int i=0; i<n; i++ ){
            double x = corners[0].x;
            double y = corners[0].y;
            if( x<xmin ) xmin = x;
            if( x>xmax ) xmax = x;
            if( y<ymin ) ymin = y;
            if( y>ymax ) ymax = y;
        }
        rect->set( xmin, xmax, ymin, ymax );
    }

    void reallocate( int n_ ){
        n = n_;
        if( corners != NULL ) delete corners;   corners = new Vec2d [ n ];
        if( lines   != NULL ) delete lines;     lines   = new Line2d[ n ];
    }

    Convex2d(){};

	Convex2d( int n_ ){
        n       = n_;
		corners = new Vec2d [ n ];
		lines   = new Line2d[ n ];
	}

	~Convex2d( ){
		delete corners;
		delete lines;
	}


// ======== temp

	void projectToLine( const Vec2d& dir, double * xs, double * ys ){   // FIXME move this to .cpp 
		// -- order corners by pos along direction
		int    permut[ n ];
		for( int i=0; i<n; i++ ){ 
			xs [i] = dir.dot     ( corners[i] ); 
		}		
		quickSort<double>( xs, permut, 0, n );
		// -- initialize left and right bonder from bottom point "ibot" to top point "itop" 
		int ibot  = permut[ 0   ];
		int itop  = permut[ n-1 ];
		int index = 0;
		ys [ 0 ]  = 0.0d;
		Vec2d oleft,oright,pleft,pright;
		oleft .set( xs[ibot], dir.dot_perp( corners[ibot] ) );
		oright.set( oleft );
		int ileft  = ibot;  _circ_inc( ileft,  n );
		int iright = ibot;  _circ_dec( iright, n );
		pleft .set( xs[ileft],  dir.dot_perp( corners[ileft ] ) );
		pright.set( xs[iright], dir.dot_perp( corners[iright] ) );
		// -- iterate over left and right border resolving points acording to its order along the direction
		do {
			index++;
			if( pleft.x < pright.x ){ // left is closer 
				double yright = oright.y +  ( pleft.x - oleft.x ) * ( pright.y - oright.y ) / ( pright.x - oright.x );
				ys[ index ]   = pleft.y - yright;
				oleft.set( pleft );
				_circ_inc( ileft,  n );
				pleft .set( xs[ileft],  dir.dot_perp( corners[ileft ] ) );
				// FIXME : we should take care when it come to end ? probably not then ileft=itop
			}else{
				double yleft = oleft.y +  ( pright.x - oright.x ) * ( pleft.y - oleft.y ) / ( pleft.x - oleft.x );
				ys[ index ]  = yleft - pright.y;
				oright.set( pright );
				_circ_dec( iright,  n );
				pright .set( xs[iright],  dir.dot_perp( corners[iright ] ) );
			}
			if( index >= (n-1) ){ break; printf( " loop should end %i %i %i %i \n", index, n, ileft, iright ); } // FIXME DEBUG just to prevent infinite loop
		} while( !( ( itop == ileft ) && ( itop == iright ) ) );
	}

};


#endif


