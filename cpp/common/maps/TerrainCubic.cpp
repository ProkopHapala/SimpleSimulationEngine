

#include <stdlib.h>
#include <stdio.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "TerrainCubic.h" // THE HEADER

double TerrainCubic::getVal( double x, double y ){
    double f_ix = getIx_f( x );
    double f_iy = getIx_f( y );
    int ix      = (int)f_ix;
    int iy      = (int)f_iy;
    double tx   = f_ix - ix;
    double ty   = f_iy - iy;
    double * h0 = heights + getIndexI( ix-1, iy );
    double * h1 = h0 + nx;
    double * h2 = h1 + nx;
    double * h3 = h2 + nx;
    return Spline_Hermite::val2D( tx, ty, h0, h1, h2, h3 );
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
};

int TerrainCubic::rayLine( Vec2d hdir, Vec2d p0, double hg0, double dr, double rmax, int ntg, double * tgs, Vec3d * poss ){
    double h0 = getVal( p0.x, p0.y )+hg0;
    double r  = 0.0d;
    for(int itg=0; itg<ntg; itg++){
        double dh = tgs[itg]*dr;
        double h  = h0 + dh*r;
        double hg = h  - getVal( p0.x+hdir.x*r, p0.y+hdir.y*r );
        while( r<rmax ){
            double r_  = r+dr;
            double h_  = h+dh;
            double hg_ = h_ - getVal( p0.x+hdir.x*r_, p0.y+hdir.y*r_ );
            if( hg_<0 ){
                r += dr*hg/(hg+hg_);
                poss[itg] = {p0.x+hdir.x*r,p0.y+hdir.y*r, h0 + dh*r };
                //printf( "%i %g %g (%g,%g,%g)\n", itg, dh, r, poss[itg].z, poss[itg].x, poss[itg].z );
                break;
            }else{
                r=r_;  h=h_;
            }
        }
    }
}

int TerrainCubic::renderRect( double x0, double y0, double x1, double y1, int nx ){
    //int ix0 = getIx( x0 );  int iy0 = getIy( y0 );
    //int ix1 = getIx( x1 );  int iy1 = getIy( y1 );
    //printf( " TerrainCubic::renderRect %f %f %f %f %i \n", x0, y0, x1, y1, nx );
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
