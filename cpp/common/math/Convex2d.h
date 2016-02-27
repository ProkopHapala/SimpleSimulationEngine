
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
	void projectToLine( const Vec2d& dir, double * xs, double * yLs, double * yRs ) const ;

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


};


#endif


