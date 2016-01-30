
#ifndef  TerrainCubic_h
#define  TerrainCubic_h

//#include "fastmath.h"
//#include "Vec2.h"
#include "Map2D.h"
#include "spline_hermite.h"

//#include "TiledView.h"

// ==============================================
// ========== CLASS STUB : GridMap2D_stub
// ==============================================

class TerrainCubic : public Map2D {
    public:
	double * heights;

    // ==== functions
    // evaluation
    double getVal( double x, double y );
    // drawing
	int    renderRect( double x0, double y0, double x1, double y1, int nx );

	// ==== inline functions

    inline void allocate( ){ heights = new double[nxy]; }

    inline void generateRandom( double vmin, double vmax ){
        for( int i=0; i<nxy; i++ ){ heights[i] = randf( vmin, vmax );  }
    }

};

#endif

