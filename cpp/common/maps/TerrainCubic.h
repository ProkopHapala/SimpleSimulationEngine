
#ifndef  TerrainCubic_h
#define  TerrainCubic_h

#include "fastmath.h"
#include "Vec2.h"
#include "Map2D.h"
#include "spline_hermite.h"

//#include "TiledView.h"

// ==============================================
// ========== CLASS STUB : GridMap2D_stub
// ==============================================

class TerrainCubic : public Map2D {
    public:

	double * heights;

	double getVal( double x, double y ){
		double f_ix = getIx_f( x );
		double f_iy = getIx_f( y );
		int ix    = (int)f_ix;
		int iy    = (int)f_iy;
		double tx = f_ix - ix;
		double ty = f_iy - iy;

/*
		int i1  = getIndex( ix-1, iy );
		int i0  = i1 - nx;
		int i2  = i1 + nx;
		int i3  = i2 + nx;
		printf( "%i %i %f %f %f %f\n", ix, iy, tx, ty, x, y );
		return Spline_Hermite::val2D( tx, ty,
			heights[i0], heights[i0+1], heights[i0+2], heights[i0+3],
			heights[i1], heights[i1+1], heights[i1+2], heights[i1+3],
			heights[i2], heights[i2+1], heights[i2+2], heights[i2+3],
			heights[i3], heights[i3+1], heights[i3+2], heights[i3+3]
		);
*/

		double * h0 = heights + getIndexI( ix-1, iy );
		double * h1 = h0 + nx;
		double * h2 = h1 + nx;
		double * h3 = h2 + nx;
		return Spline_Hermite::val2D( tx, ty, h0, h1, h2, h3 );

	};

    int renderRect( double x0, double y0, double x1, double y1, int nx ){
		//int ix0 = getIx( x0 );  int iy0 = getIy( y0 );
		//int ix1 = getIx( x1 );  int iy1 = getIy( y1 );
		printf( " TerrainCubic::renderRect %f %f %f %f %i \n", x0, y0, x1, y1, nx );
		int ny       = nx/0.86602540378;
		float dx = (x1-x0)/nx;
        float dy = (y1-y0)/ny;
        float dxhalf = 0.5f*dx;
		int nverts=0;
		float ylo,yhi;
		for( int iy = 0; iy<ny; iy++ ){
            if( iy & 1 ){ ylo = y0 + iy*dy; yhi = ylo+dy; }else{  yhi = y0 + iy*dy; ylo = yhi+dy;  }
            float x     = x0;
            glBegin( GL_TRIANGLE_STRIP );
			for( int ix = 0; ix<=nx; ix++ ){
                float val;
                val = (float) getVal( x, yhi ); glColor3f( val, val, val ); glVertex3f( x, yhi, 0 ); x+=dxhalf;
                val = (float) getVal( x, ylo ); glColor3f( val, val, val ); glVertex3f( x, ylo, 0 ); x+=dxhalf;
				nverts+=2;
			}
            glEnd();
		}
		return nverts;
    }

    void allocate( ){ heights = new double[nxy]; }

    void generateRandom( double vmin, double vmax ){
        for( int i=0; i<nxy; i++ ){ heights[i] = randf( vmin, vmax );  }
    }

/*
    int renderTile( double x0, double y0, double x1, double y1, int nx ){
            renderRect( x0, y0, x1, y1, nx );
        return( ilist );
    }
*/

};

#endif

