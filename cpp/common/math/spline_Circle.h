
#ifndef  Spline_Circle_h
#define  Spline_Circle_h

/*
# Circle Spline is a Line composed of circulr segments with conrners
*/

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "geom2D.h"

/*
struct CircPoint{
    Vec2d  p;
    double r;
}
*/

class CircleSpline{ public:

    std::vector<Circle2d> CPs;
    Ray2d* rays=0;
    Arc2d* arcs=0;

    void realloc_side(){
        int n=CPs.size();
        //printf( "arcs %i %i \n", (long)arcs, (long)rays );
        //_realloc( rays, n-1 );
        //_realloc( arcs, n-2 );
        _realloc( rays, n ); // closed
        _realloc( arcs, n );
    }

    inline int pickNode(const Vec2d& p){
        int    ipick   =-1;
        double dist=0;
        for(int i=0; i<CPs.size(); i++){
            double r = CPs[i].pointDist(p)/ sq(CPs[i].r);
            if(r<dist){ dist=r; ipick=i; }
        }
        return ipick;
    }

    inline bool getLine(int i, Vec2d& p0, Vec2d& p1, int isign ){
        // R1    = sin(a) * ( L + K )
        // R2    = sin(a) * ( K )
        // R1-R2 = sin(a) * L
        //if((i<0)||(i>=(CPs.size()-1)))return false;
        int n = CPs.size();
        Circle2d c1 = CPs[ wrap_index_fast(i  ,n) ];
        Circle2d c2 = CPs[ wrap_index_fast(i+1,n) ];
        bool bswap=c1.r<c2.r;
        if(bswap){ _swap( c1,c2 ); isign=-isign; }
        Vec2d  d  = c2.p0 - c1.p0;
        double il = 1/d.norm();
        double sa = (c1.r - c2.r)*il;
        double ca = sqrt(1-sa*sa)*isign;
        Vec2d hat; hat.set_mul_cmplx( d, {sa*il,ca*il} );
        p0.set_add_mul(c1.p0,hat,c1.r);
        p1.set_add_mul(c2.p0,hat,c2.r);
        //printf( "p0(%g,%g) p1(%g,%g)\n",  );
        if(bswap){ _swap( p0,p1 ); }
        return bswap;
    }

    //void getArc(int i, Vec2d& d0, double& da ){}

    void evalSide( int isign ){
        if( !arcs ){ realloc_side(); };
        int n=CPs.size();
        // lines
        for(int i=0; i<n; i++){
            Ray2d& ray = rays[ i ];
            Vec2d p0,p1;
            getLine( i, p0, p1, isign );
            ray.fromPoints(p0,p1);
            //printf( "ray (%g,%g) (%g,%g) l %g \n", ray.p0.x,ray.p0.y,  ray.dir.x,ray.dir.y,  ray.l );
        }
        // points
        for(int i=0; i<n; i++){
            Arc2d& arc = arcs[i];
            int im=(i>0)?(i-1):(n-1);
            Vec2d op;
            rays[im].getEnd(op);
            //glColor3f(1.0,0.0,0.0); Draw2D::drawPointCross_d( op, 0.1 ); Draw2D::drawLine_d( op, rays[im].p0 );
            //glColor3f(0.0,0.0,1.0); Draw2D::drawPointCross_d( rays[i].p0, 0.1 );
            arc.fromCenter2points( CPs[i].p0, op, rays[i].p0 );
        }
    }

    ~CircleSpline(){
        _dealloc(arcs);
        _dealloc(rays);
    }

};

#endif



