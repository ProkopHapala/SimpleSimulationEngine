
#ifndef  TiledView_h
#define  TiledView_h

#include <vector>

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"
#include "Map2D.h"

// ==============================================
// ========== CLASS STUB : GridMap2D_stub
// ==============================================

class TiledView : public Map2D {
    public:
    int * tiles;
    int itx0,ity0;
    int isx,isy;
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
        printf( " TiledView ::reconstruct nx, ny, step ,invStep %i %i %f %f\n", nx, ny, step ,invStep );
        float szx = x1-x0;
        float szy = y1-y0;
        step = fmax( szx, szy )/perScreen; invStep = 1/step;
        printf( "--TiledView ::reconstruct nx, ny, perScreen, szx, szy, step ,invStep %i %i %i %f %f %f %f\n", nx, ny, perScreen, szx, szy, step ,invStep );
        isx     = (nx-perScreen)/2;
        isy     = (ny-perScreen)/2;
        itx0    = getIx( x0 ) - isx;
        ity0    = getIy( y0 ) - isy;
        printf( "--TiledView ::reconstruct isx, isy, itx0 ,ity0 %i %i %i %i\n", isx, isy, itx0 ,ity0 );
        for( int iy=0; iy<ny; iy++ ){
            for( int ix=0; ix<nx; ix++ ){
                int i    = getIndexI( ix, iy );
                double x = getX(ix - isx);
                double y = getY(iy - isy);
                printf( " %i %i %i %i %f %f\n", ix, iy, i, nxy, x ,y );
                //if ( (ix&1)^(iy&1) )
                tiles[i] = tileToList( x, y, x+step, y+step );
            }
        }
    };

    void render( float x0, float y0, float x1, float y1){
        int ix0 = getIx( x0 )-itx0+1;  int iy0 = getIy( y0 )-ity0+1;
        int ix1 = getIx( x1 )-itx0+1;  int iy1 = getIy( y1 )-ity0+1;
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
        glColor3f( 0.2f, 0.9f, 0.2f ); Draw2D::drawRectangle( x0, y0, x1, y1, false );
    };

    void init(  int nx_, int ny_, int perScreen_ ){
        perScreen = perScreen_;
        setnxy( nx_, ny_ );
        tiles = new int[ nxy ];
    }

    ~TiledView(){
        if( tiles!= NULL ) delete tiles;
        for( int i=0; i<nxy; i++ ){ tiles[i] = 0; }
    }

};

#endif

