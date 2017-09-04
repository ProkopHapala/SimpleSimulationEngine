
#ifndef  Ruler2DFast_h
#define  Ruler2DFast_h

#include "Vec2.h"

const static int Ruler2D_nEdges = 4;
const static int Ruler2D_nVerts = 4;

// ===== class CubicRuler

class Ruler2DFast { public:
    // --- variables
    Vec2d pos0    = (Vec2d){0.0,0.0};
	Vec2d step    = (Vec2d){1.0,1.0};
	Vec2d invStep = (Vec2d){1.0,1.0};

    Vec2i n       = (Vec2i){0,0};
	int   ntot    = 0;

	// --- inline functions

    inline void  setN(Vec2i n_)       { n=n_; ntot=n.x*n.y;                }
    inline Vec2i i2ip(int   i ) const { return {i%n.x,i/n.x};   } // https://stackoverflow.com/questions/7070346/c-best-way-to-get-integer-division-and-remainder
    inline int   ip2i(Vec2i ip) const { return (n.x*ip.y+ip.x); }

	inline int    x2i( double x  ) const { return (int)( invStep.x*(x - pos0.x ) ); }
	inline int    y2i( double y  ) const { return (int)( invStep.y*(y - pos0.y ) ); }
	inline double i2x( double ix ) const { return      (step.x*ix)   + pos0.x; }
	inline double i2y( double iy ) const { return      (step.y*iy)   + pos0.y; }

	//inline int dx ( double ix, double x ) const { return (x-pos0.x) - (ix*step.x);  }
	//inline int dxi( double ix, double x ) const { return invStep.x*(x-pos0.x) - ix; }

	inline int x2id( double x, double& dix ) const { double x_ = invStep.x*(x-pos0.x); int ix = (int)x_; dix = x_ - ix; return ix; }
	inline int y2id( double y, double& diy ) const { double y_ = invStep.y*(y-pos0.y); int iy = (int)y_; diy = y_ - iy; return iy; }

	inline int x2idx( double x, double& dx ) const { double x_ = invStep.x*(x-pos0.x); int ix = (int)x_; dx = step.x*(x - ix); return ix; }
	inline int y2idy( double y, double& dy ) const { double y_ = invStep.y*(y-pos0.y); int iy = (int)y_; dy = step.y*(y - iy); return iy; }

	inline void pos2index( const Vec2d& pos, Vec2d& dipos, Vec2i& ipos ) const {
		ipos.x = x2id( pos.x, dipos.x );
		ipos.y = y2id( pos.y, dipos.y );
	}

    inline void index2pos( const Vec2i& ipos, const Vec2d& dipos, Vec2d& pos ) const {
        pos.x = i2x( ipos.x + dipos.x );
        pos.y = i2y( ipos.y + dipos.y );
    }

    inline void setup( const Vec2d& pos0_, const Vec2d& step_ ){
        pos0  .set( pos0_ );
        step  .set( step_  );
		invStep.set_inv( step );
    }


    int getOverlapingTiles( Vec2d pos, double r, int * results ){
        Vec2d dpos;
        Vec2i ipos;
        pos2index( pos, dpos, ipos );
        dpos.mul( step );
        results[0] = ip2i(ipos); int nret=1;
        int dix=0,diy=0;
        double dr2   = 0;
        double mr    = 1-r;
        if     (  dpos.x < r  ){ results[nret]=ip2i( {ipos.x-1  , ipos.y    }); nret++; dix=-1; dr2 += sq(  dpos.x); }
        else if(  dpos.x > mr ){ results[nret]=ip2i( {ipos.x+1  , ipos.y    }); nret++; dix=+1; dr2 += sq(1-dpos.x); }
        if     (  dpos.y < r  ){ results[nret]=ip2i( {ipos.x    , ipos.y-1  }); nret++; diy=-1; dr2 += sq(  dpos.y); }
        else if(  dpos.y > mr ){ results[nret]=ip2i( {ipos.x    , ipos.y+1  }); nret++; diy=+1; dr2 += sq(1-dpos.y); }
        if     ( dr2 < (r*r)  ){ results[nret]=ip2i( {ipos.x+dix, ipos.y+diy}); }
        //if( (dix!=0)&&(diy!=0) ){ insert( o, ipos.x+dix, ipos.y+diy ); }
        //printf( " %1.3f %1.3f  (%1.3f,%1.3f) (%i,%i) %1.3f \n", r, mr, dpos.x,dpos.y, dix, diy, dr2 );
        return nret;
    }

};

#endif

