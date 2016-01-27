
#ifndef  Map2D_h
#define  Map2D_h

#include <vector>

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"

// ==============================================
// ========== CLASS STUB : GridMap2D_stub
// ==============================================

class Map2D {
    public:
    double        x0=0,y0=0;
	double        step, invStep;
	int           nx, ny, nxy;

	inline int    getIx  ( double  x ){ return (int)( invStep * (x-x0) ); };
	inline int    getIy  ( double  y ){ return (int)( invStep * (y-y0) ); };
    inline double getIx_f( double  x ){ return      ( invStep * (x-x0) ); };
	inline double getIy_f( double  y ){ return      ( invStep * (y-y0) ); };
	inline double getX   ( int    ix ){ return (step * ix)+x0; };
	inline double getY   ( int    iy ){ return (step * iy)+y0; };

	inline int    getIndexI( int ix,   int iy   ){ return nx*iy + ix;                     };
	inline int    getIndex ( double x, double y ){ return getIndex( getIx(x), getIy(y) ); };

    inline void init( int nx_, int ny_, double step_ ){
        step    = step_;
		invStep = 1/step;
		nx = nx_; ny=ny_;
		nxy = nx*ny;
    }

};

#endif

