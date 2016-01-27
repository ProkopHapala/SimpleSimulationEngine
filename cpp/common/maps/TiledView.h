
#ifndef  TiledView_h
#define  TiledView_h

#include <vector>

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"

// ==============================================
// ========== CLASS STUB : GridMap2D_stub
// ==============================================

class TiledView : public Map2D {
    public:
    int * tiles;
    int itx0,ity0;
    int perScreen;

    virtual int tileToList( float x0, float y0, float x1, float y1 )=0;

/*
    int update( float x0, float y0, float x1, float x1 ){
        float szx = x1-x0;
        float szy = y1-y0;
        if      ( ( szx > 1.5*step ) || ( szy > 1.5*step ) ){
            reconstruct( x0, y0, x1, x1 );
        }else if( ( szx < 0.7*step ) || ( szy < 0.7*step ) ){
            reconstruct( x0, y0, x1, x1 );
        }else{
            int ix0 = getIx( x0 );  int iy0 = getIy( y0 );
            int ix1 = getIx( x1 );  int iy1 = getIy( y1 );
        }
    }
*/

    void reconstruct( float x0, float y0, float x1, float y1 ){
        float szx = x1-x0;
        float szy = y1-y0;
        step = fmax( szx, szy )/perScreen; invStep = 1/step;
        int isx = (nx-perScreen)/2;
        int isy = (ny-perScreen)/2;
        itx0 = getIx( x0 ) - isx;
        ity0 = getIy( y0 ) - isy;
        for( int iy=0; iy<=ny; iy++ ){
            for( int ix=0; ix<=nx; ix++ ){
                int i = getIndexI( ix, iy );
                double x = getX(ix - isx);
                double y = getY(iy - isy);
                tiles[i] = tileToList( x, y, x+step, y+step );
            }
        }
    };

    void render( float x0, float y0, float x1, float y1){
        int ix0 = getIx( x0 )-itx0;  int iy0 = getIy( y0 );
        int ix1 = getIx( x1 )-itx0;  int iy1 = getIy( y1 );
        for( int iy=iy0; iy<=iy1; iy++ ){
            if( ( iy > 0 ) && ( iy<ny ) ){
                for( int ix=ix0; ix<=ix1; ix++ ){
                    if( ( ix > 0 ) && ( ix<nx ) ){
                        int i = getIndexI( ix, iy );
                        if( tiles[i]!=0 ) glCallList( tiles[i] );
                    }
                }
            }
        }
    };

};

#endif

