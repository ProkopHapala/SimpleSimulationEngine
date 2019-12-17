
// Look also here
// common/maps/  Ruler2DFast
// common/maps/  CubeGridRuler

#ifndef  SquareRuler_h
#define  SquareRuler_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "GridIndex2D.h"

#include "spline_hermite.h"

class SquareRuler : public GridIndex2D { public:
    double    step = 1.0d, invStep=1.0d;
    Vec2d     pos  = (Vec2d){0.0,0.0};

    // ray tracing
    Vec3i  ray_i,ray_si;
    double ray_t;
    int maxRayIter = 100;

    // --- inline functions

    inline int    x2i( double x  ) const { return (int)( invStep*(x - pos.x ) ); }
    inline int    y2i( double y  ) const { return (int)( invStep*(y - pos.y ) ); }
    inline double i2x( double ix ) const { return      (step*ix)    + pos.x; }
    inline double i2y( double iy ) const { return      (step*iy)    + pos.y; }

    //inline int dx ( double ix, double x ) const { return (x-pos0.x) - (ix*step.x);  }
    //inline int dxi( double ix, double x ) const { return invStep.x*(x-pos0.x) - ix; }

    inline int x2id( double x, double& dix ) const { double x_ = invStep*(x-pos.x); int ix = (int)x_; dix = x_ - ix; return ix; }
    inline int y2id( double y, double& diy ) const { double y_ = invStep*(y-pos.y); int iy = (int)y_; diy = y_ - iy; return iy; }

    inline int x2idx( double x, double& dx ) const { double x_ = invStep*(x-pos.x); int ix = (int)x_; dx = step*(x - ix); return ix; }
    inline int y2idy( double y, double& dy ) const { double y_ = invStep*(y-pos.y); int iy = (int)y_; dy = step*(y - iy); return iy; }

    inline void pos2index( const Vec2d& pos, Vec2i& ipos, Vec2d& dipos ) const {
	    ipos.x = x2id( pos.x, dipos.x );
	    ipos.y = y2id( pos.y, dipos.y );
    }

    inline void index2pos( const Vec2i& ipos, const Vec2d& dipos, Vec2d& pos ) const {
        pos.x = i2x( ipos.x + dipos.x );
        pos.y = i2y( ipos.y + dipos.y );
    }

    inline void setup( const Vec2d& pos_, double step_ ){
        pos  .set( pos_ );
        step = step_, invStep=1/step;
    }

    inline double getValue_cubic( Vec2d p, const double * hs )const { // TODO somthing related to HermiteSpline
        Vec2i ind; Vec2d dind;
        pos2index( p, ind, dind );
        double h00 = fetchWraped( {ind.a-1,ind.b-1}, hs );
        double h01 = fetchWraped( {ind.a  ,ind.b-1}, hs );
        double h02 = fetchWraped( {ind.a+1,ind.b-1}, hs );
        double h03 = fetchWraped( {ind.a+2,ind.b-1}, hs );
        double h10 = fetchWraped( {ind.a-1,ind.b  }, hs );
        double h11 = fetchWraped( {ind.a  ,ind.b  }, hs );
        double h12 = fetchWraped( {ind.a+1,ind.b  }, hs );
        double h13 = fetchWraped( {ind.a+2,ind.b  }, hs );
        double h20 = fetchWraped( {ind.a-1,ind.b+1}, hs );
        double h21 = fetchWraped( {ind.a  ,ind.b+1}, hs );
        double h22 = fetchWraped( {ind.a+1,ind.b+1}, hs );
        double h23 = fetchWraped( {ind.a+2,ind.b+1}, hs );
        double h30 = fetchWraped( {ind.a-1,ind.b+2}, hs );
        double h31 = fetchWraped( {ind.a  ,ind.b+2}, hs );
        double h32 = fetchWraped( {ind.a+1,ind.b+2}, hs );
        double h33 = fetchWraped( {ind.a+2,ind.b+2}, hs );
        return Spline_Hermite::val2D( dind.x, dind.y, h00,h01,h02,h03, h00,h01,h02,h03, h00,h01,h02,h03, h00,h01,h02,h03 );
    }

    inline double getDeriv( Vec2d p, Vec2d& dF, double * hs )const { // TODO somthing related to HermiteSpline
        Vec2i ind; Vec2d dind;
        pos2index( p, ind, dind );
        double h00 = fetchWraped( {ind.a-1,ind.b-1}, hs );
        double h01 = fetchWraped( {ind.a  ,ind.b-1}, hs );
        double h02 = fetchWraped( {ind.a+1,ind.b-1}, hs );
        double h03 = fetchWraped( {ind.a+2,ind.b-1}, hs );
        double h10 = fetchWraped( {ind.a-1,ind.b  }, hs );
        double h11 = fetchWraped( {ind.a  ,ind.b  }, hs );
        double h12 = fetchWraped( {ind.a+1,ind.b  }, hs );
        double h13 = fetchWraped( {ind.a+2,ind.b  }, hs );
        double h20 = fetchWraped( {ind.a-1,ind.b+1}, hs );
        double h21 = fetchWraped( {ind.a  ,ind.b+1}, hs );
        double h22 = fetchWraped( {ind.a+1,ind.b+1}, hs );
        double h23 = fetchWraped( {ind.a+2,ind.b+1}, hs );
        double h30 = fetchWraped( {ind.a-1,ind.b+2}, hs );
        double h31 = fetchWraped( {ind.a  ,ind.b+2}, hs );
        double h32 = fetchWraped( {ind.a+1,ind.b+2}, hs );
        double h33 = fetchWraped( {ind.a+2,ind.b+2}, hs );
        return Spline_Hermite::dval2D( dind.x, dind.y, dF.x, dF.y, h00,h01,h02,h03, h00,h01,h02,h03, h00,h01,h02,h03, h00,h01,h02,h03 );
    }

    inline void setStep( double step_ ){ step = step_; invStep = 1.0d/step; };

    void rayStart( Vec2d ray0, Vec2d hray ){}
    void rayStep(){}

    void sampleLine( int n, Vec2d p, Vec2d dp, const double * vmap, double* vals ){
        for(int i=0; i<n; i++){ vals[i] = getValue_cubic( p, vmap ); p.add(dp); }
    }

    void rayCut( Vec2d ray0, Vec2d hray, const double * vmap, double * ts, double* vals, double tmax ){}
    void rayHorizonts( Vec2d p, Vec2d hray, const double * vmap, double rayh0, int ndhs, double * dhs, double * ts, double tmax ){}

};

#endif

