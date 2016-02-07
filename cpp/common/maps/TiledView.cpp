
#include <stdlib.h>
#include <stdio.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "drawMath2D.h"

#include "TiledView.h" // THE HEADER

void TiledView::printIndexes( ){
    for( int iy=0; iy<ny; iy++ ){
        for( int ix=0; ix<nx; ix++ ){
            int i   = getIndexI( ix, iy );
            printf( "%i ", tiles[i] );
        }
        printf( "\n" );
    }
}

void TiledView::draw( float xmin, float ymin, float xmax, float ymax ){
    checkRender( xmin, ymin, xmax, ymax );
    draw_raw   ( xmin, ymin, xmax, ymax );
    //printf( "TiledView::draw %3.3f %3.3f %3.3f %3.3f  \n", xmin, ymin, xmax, ymax );
    //glColor3f( 0.2f, 0.9f, 0.2f ); Draw2D::drawRectangle( xmin, ymin, xmax, ymax, false );
};


void TiledView::renderAll( float xmin, float ymin, float xmax, float ymax ){
    //float szx = xmax-xmin; float stepX = 2*szx/nx;
    //float szy = ymax-ymin; float stepY = 2*szy/ny;
    float szx = xmax-xmin; float stepX = szx/(nx-2);
    float szy = ymax-ymin; float stepY = szy/(ny-2);
    setStep( fmax( stepX, stepY ) );
    x0 = xmin-0.5*step; y0 = ymin-0.5*step;
    printf( " TiledView ::reconstruct nx, ny, step ,invStep %i %i %3.3f %3.3f  osz,osy %3.3f %3.3f\n", nx, ny, step ,invStep, step*nx, step*ny );
    for( int iy=0; iy<ny; iy++ ){
        for( int ix=0; ix<nx; ix++ ){
            int i    = getIndexI( ix, iy );
            double x = getX(ix);
            double y = getY(iy);
            //printf( " %i %i %i %i %f %f\n", ix, iy, i, nxy, x ,y );
            //if ( (ix&1)^(iy&1) )
            tiles[i] = tileToList( x, y, x+step, y+step );
        }
    }
};

void TiledView::shiftRender( int dix, int diy ){
    int tiles_[nxy];
    for( int i=0; i<nxy; i++ ){ tiles_[i]=tiles[i]; }
    printf( "TiledView::shiftRender %i %i \n", dix, diy );
    //printIndexes( );
    x0 += dix*step;
    y0 += diy*step;
    int ndeleted  = 0;
    int nrendered = 0;
    for( int iy=0; iy<ny; iy++ ){
        int iy__ = iy - diy;
        bool toDelete = ( iy__ < 0 ) || ( iy__ >= ny );
        int iy_  = iy + diy;
        bool toRender = ( iy_  < 0 ) || ( iy_  >= ny );
        for( int ix=0; ix<nx; ix++ ){
            int i   = getIndexI( ix, iy );
            // delete tiles outside the box
            int ix__ = ix - dix;
            if ( toDelete || ( ix__ < 0 ) || ( ix__ >= nx ) ){
                glDeleteLists( tiles_[i], 1 ); ndeleted++;
            }
            // fill in proper tiles
            int ix_ = ix + dix;
            if ( toRender || ( ix_ < 0 ) || ( ix_ >= nx ) ){
                // render tiles newly in box
                double x = getX(ix);
                double y = getY(iy);
                tiles[i] = tileToList( x, y, x+step, y+step ); nrendered++;
                //printf( "%i \n", tiles[i] );
            } else{
                // relocate useful old tiles
                int i_   = getIndexI( ix_, iy_ );
                tiles[i] = tiles_[i_];
            }
        }
    }
    //printf( "deleted %i rendered %i\n", ndeleted, nrendered );
    //printIndexes( );
    //exit(0);
}

bool TiledView::checkRender( float xmin, float ymin, float xmax, float ymax ){
    float oszx = nx*step;
    float oszy = ny*step;
    float szx  = xmax-xmin;
    float szy  = ymax-ymin;
    if       ( (szx>(oszx-step) ) || ( szy>(oszy-step) ) ){
        //printf( "oszx,oszy, szx,szy   %3.3f %3.3f   %3.3f %3.3f \n", oszx, oszy, szx,szy );
        //printf( " TiledView::reconstruct Zoom out \n" );
        renderAll( xmin, ymin, xmax, ymax );
        return true;
    }else if ( (szx<(0.3*oszx)) && (szy<(0.3*oszy)) ){
        //printf( "oszx,oszy, szx,szy   %3.3f %3.3f   %3.3f %3.3f \n", oszx, oszy, szx,szy );
        //printf( " TiledView::reconstruct Zoom in \n" );
        renderAll( xmin, ymin, xmax, ymax );
        return true;
    }else {
        bool shift = false;
        int dix=0,diy=0;
        if      (xmin<x0       ){ dix = (int)((xmin-x0     )*invStep)-1; shift=true; }
        else if (xmax>(x0+oszx)){ dix = (int)((xmax-x0-oszx)*invStep)+1; shift=true; }
        if      (ymin<y0       ){ diy = (int)((ymin-y0     )*invStep)-1; shift=true; }
        else if (ymax>(y0+oszy)){ diy = (int)((ymax-y0-oszy)*invStep)+1; shift=true; }
        if (shift){
            printf( " %3.3f %3.3f    %3.3f %3.3f | %3.3f %3.3f    %3.3f %3.3f   \n", xmin, xmax, x0, x0+oszx,     ymin, ymax, y0, y0+oszy );
            shiftRender( dix, diy );
        }
        return true;
    }
    return false;
}

void TiledView::draw_raw( float xmin, float ymin, float xmax, float ymax ){
    int ix0 = getIx( xmin );  int iy0 = getIy( ymin );
    int ix1 = getIx( xmax );  int iy1 = getIy( ymax );
    for( int iy=iy0; iy<=iy1; iy++ ){
        if( ( iy >= 0 ) && ( iy<ny ) ){
            for( int ix=ix0; ix<=ix1; ix++ ){
                if( ( ix >= 0 ) && ( ix<nx ) ){
                    int i = getIndexI( ix, iy );
                    if( tiles[i]!=0 ) glCallList( tiles[i] );
                }
            }
        }
    }
    //printf( "TiledView::draw %3.3f %3.3f %3.3f %3.3f  \n", xmin, ymin, xmax, ymax );
    glColor3f( 0.2f, 0.9f, 0.2f ); Draw2D::drawRectangle( xmin, ymin, xmax, ymax, false );
};

