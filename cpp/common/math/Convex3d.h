
#ifndef  Convex3d_h
#define  Convex3d_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"

#include "arrayAlgs.h" // FIXME move this to .cpp


/*

Comvex Hull to Mesh pseudocode

- for each if N planes
    - evaluate all (N-1) lines resulting from intersection with each other
    - trim lines by other planes

Convex Hull pseudocode (slow):




*/




class Convex2d{ public:
    std::vector<Plane3D> planes;
    std::vector<Vec2i>   edges;
    std::vector<Vec3d>   verts;



void toMesh(){
    Line3d lines[N]; // some lines may be empty
    for(int i=0; i<N; i++){
        Plane3d& pi = planes[i];
        for(int j=i+1; j<N; j++){
            LineSegment3d lij; lij.fromLine( pi, planes[j] );
            int k1,k2;
            for( int k=i+1; k<N; k++ ){
                if( k==j ) continue; // cannot cut by self
                int side = lij.trim( plane[k] ); // side can be 1,-1,0 depending on which side of line segment was trimmed
                if      ( side > 0 ){ k1=k; }
                else if ( side < 0 ){ k2=k; };
            }
            pushVert( i,j,k );
        }
    }
}




/*

void fromPoints(int n, Vec3d* points){
    // inspierd by http://mindthenerd.blogspot.com/2012/05/fastest-convex-hull-algorithm-ever.html
    std::vector<>;
    double xmin,xmax,ymin,ymax,zmin,zmax;
    int    ixmin,ixmax,iymin,iymax,izmin,izmax;
    for(int i=0; i<n; i++){
        if( xmin > points[i].x ){ xmin = points[i].x; ixmin=i; };
        if( xmax < points[i].x ){ xmax = points[i].x; ixmax=i; };
        if( ymin > points[i].y ){ ymin = points[i].y; iymin=i; };
        if( ymax < points[i].y ){ ymax = points[i].y; iymax=i; };
        if( zmin > points[i].z ){ zmin = points[i].x; izmin=i; };
        if( zmax < points[i].z ){ zmax = points[i].x; izmax=i; };
    }
    Plane3D f000,f001,f010,f011,f100,f101,f110,f111;
    f000.fromPoints( points[ixmin], points[iymin], points[izmin] );
    f001.fromPoints( points[ixmax], points[iymin], points[izmin] );
    f010.fromPoints( points[ixmin], points[iymax], points[izmin] );
    f011.fromPoints( points[ixmax], points[iymax], points[izmin] );
    f100.fromPoints( points[ixmin], points[iymin], points[izmax] );
    f101.fromPoints( points[ixmax], points[iymin], points[izmax] );
    f110.fromPoints( points[ixmin], points[iymax], points[izmax] );
    f111.fromPoints( points[ixmax], points[iymax], points[izmax] );
    // this can be perhaps optimized so that we dont have to allocate 8 vectors;
    std::vector<Vec3d> ps000,ps001,ps010,ps011,ps100,ps101,ps110,ps111;
    double c000,c001,c010,c011,c100,c101,c110,c111;
    int    ixmin,ixmax,iymin,iymax,izmin,izmax;
    for(int i=0; i<n; i++){
        double cl
        c = p000.dist( points[i]; if( 0 < c ) ) p000.push( points[i] );
        if( 0 < p001.dist( points[i] ) ) p001.push( points[i] );
        if( 0 < p010.dist( points[i] ) ) p010.push( points[i] );
        if( 0 < p011.dist( points[i] ) ) p011.push( points[i] );
        if( 0 < p100.dist( points[i] ) ) p100.push( points[i] );
        if( 0 < p101.dist( points[i] ) ) p101.push( points[i] );
        if( 0 < p110.dist( points[i] ) ) p110.push( points[i] );
        if( 0 < p111.dist( points[i] ) ) p111.push( points[i] );
    }
}
*/


};


#endif


