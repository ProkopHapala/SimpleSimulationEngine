
#ifndef  TiledView_h
#define  TiledView_h

#include "Map2D.h"

// ==============================================
// ========== CLASS STUB : GridMap2D_stub
// ==============================================

class TiledView : public Map2D {
    public:
    int * tiles;
    // ==== virtual functions
    virtual int tileToList( float xmin, float ymin, float xmax, float ymax )=0;
     // ==== functions
    void renderAll   ( float xmin, float ymin, float xmax, float ymax );
    bool checkRender ( float xmin, float ymin, float xmax, float ymax );
    void draw_raw    ( float xmin, float ymin, float xmax, float ymax );
    void draw        ( float xmin, float ymin, float xmax, float ymax );
    void shiftRender ( int dix, int diy );
    void printIndexes( );

    // ==== inline functions

    void init(  int nx_, int ny_ ){
        setnxy( nx_, ny_ );
        tiles = new int[ nxy ];
    }

    ~TiledView(){
        if( tiles!= NULL ) delete tiles;
        for( int i=0; i<nxy; i++ ){ tiles[i] = 0; }
    }

};

#endif

