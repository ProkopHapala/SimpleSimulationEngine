#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cstdint>

#include "lineSearch.h"

// area of polygon without triangulation
// http://stackoverflow.com/questions/451426/how-do-i-calculate-the-area-of-a-2d-polygon
// http://geomalgorithms.com/a01-_area.html#area2D_Triangle%28%29




double getHullVolume2D( int ntri, int * trinagles, Vec2d * points, Line2d line ){
	for( int i=0; i<; i++ ){
		
	}	
}

double getHullMoment2D( int ntri, int * trinagles, Vec2d * points, Line2d line ){

}

void double getHullMoment( int ntri, int * trinagles, Vec2d * points, double angle, double volume ){
	double angle = angles[i];
	double a = sin( angle );
	double b = cos( angle );
	c = lineSearch_Brent( // this will find the watterline particular ship displacement
		[=](double x){ Line2d line; line.set( a,b,c ); return getHullVolume2D( ntri, trinagles, points, line ); },
		double x0, double a, double b, double tol 
	);
	Line2d line; line.set( a,b,c );
	return getHullMoment2D( ntri, trinagles, points, line );
}

void staticStability2D( 
	double volume,
	int    nTri,   int * trinagles, Vec2d * points, 
	int nAngles,   double * angles, double * moments, 
){
	double c = 0.0d;
	for( int i=0; i<nAngles; i++ ){
		moments[i] = getHullMoment( ntri, trinagles, points, angles[i], volume );
	}	
}

extern "C"{



}
