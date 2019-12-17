
#include <stdlib.h>
#include <stdio.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "fastmath.h"

#include "TerrainRBF.h" // THE HEADER

double TerrainRBF::getVal( double x, double y ){
    //ULONG i = getBucket(x,y);
    UINT nfound = getBucketObjects( x, y, buff );
    //UINT nfound = HashMap2D<RBF2D>::getBucketObjects( x, y, buff );


    //printf( " nfound %i at point (%3.3f,%3.3f) \n", nfound, x,y );

    double val = 0;
    for( int i=0; i<nfound; i++ ){
        val += buff[i]->getVal( {x,y} );
    }
    return val;
};

int TerrainRBF::renderRect    ( double xmin, double ymin, double xmax, double ymax, int nx ){
    //int ix0 = getIx( x0 );  int iy0 = getIy( y0 );
    //int ix1 = getIx( x1 );  int iy1 = getIy( y1 );
    //printf( " TerrainCubic::renderRect %f %f %f %f %i \n", x0, y0, x1, y1, nx );
    int ny       = nx/0.86602540378;
    float dx     = (xmax-xmin)/nx;
    float dy     = (ymax-ymin)/ny;
    float dxhalf = 0.5f*dx;
    int nverts=0;
    float ylo,yhi;
    for( int iy = 0; iy<ny; iy++ ){
        if( iy & 1 ){ ylo = ymin + iy*dy; yhi = ylo+dy; }else{  yhi = ymin + iy*dy; ylo = yhi+dy;  }
        float x     = xmin;
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

int TerrainRBF::insertRBF( RBF2D * rbf ){
    int ninserts=0;
    double xmin = rbf->pos.x - rbf->rmax;   UHALF ixmin = getIx( xmin );
    double xmax = rbf->pos.x + rbf->rmax;   UHALF ixmax = getIx( xmax );
    double ymin = rbf->pos.y - rbf->rmax;   UHALF iymin = getIy( ymin );
    double ymax = rbf->pos.y + rbf->rmax;   UHALF iymax = getIy( ymax );
    for( UHALF iy=iymin; iy<=iymax; iy++ ){
        for( UHALF ix=ixmin; ix<=ixmax; ix++ ){
            //rgf->getVal( const Vec2d& where ); // TODO : checking for rect / circle overlap would narrow the selection
            if( 0<= insertIfNewInt( rbf, ix, iy ) ) ninserts++;
        }
    }
    return ninserts;
}

int TerrainRBF::removeRBF( RBF2D * rbf ){
    int ninserts=0;
    double xmin = rbf->pos.x - rbf->rmax;   UHALF ixmin = getIx( xmin );
    double xmax = rbf->pos.x + rbf->rmax;   UHALF ixmax = getIx( xmax );
    double ymin = rbf->pos.y - rbf->rmax;   UHALF iymin = getIy( ymin );
    double ymax = rbf->pos.y + rbf->rmax;   UHALF iymax = getIy( ymax );
    for( UHALF iy=iymin; iy<=iymax; iy++ ){
        for( UHALF ix=ixmin; ix<=ixmax; ix++ ){
            //rgf->getVal( const Vec2d& where ); // TODO : checking for rect / circle overlap would narrow the selection
            if ( tryRemoveInt( rbf, ix, iy ) ) ninserts++;
        }
    }
    return ninserts;
}

int TerrainRBF::insertRBFs( int n_, RBF2D * rbfs_ ){
    nrbfs  = n_;
    rbfs   = rbfs_;
    int ninserts = 0;
    for( int i=0; i<nrbfs; i++ ){ ninserts += insertRBF( rbfs+i ); }
    return ninserts;
}

/*
void makeCliff( RBF2D * rbf ){
    double rf = randf(0.0,1.0);
    if( rf > 0.5 ){
        if( rf > 0.9 ){
            rbf->ncliffs = 2;
            rbf->cliffs  = new Line2d[2];
        }else{
            rbf->ncliffs = 2;
            rbf->cliffs  = new Line2d[2];
        }
        for(int i=0; i<rbf->ncliffs; i++){
            rbf->cliffs[i].dir.set( randf( -0.5, 0.5 ), randf( -0.5, 0.5 ) );
            //rbf->dit.set( randf( -0.5, 0.5 ), randf( -0.5, 0.5 ) );
        }
    }
}
*/

void TerrainRBF::generateRandom( double xmin, double ymin, double xmax, double ymax, double rmin, double rmax,  double vmin, double vmax, int n_ ){
    RBF2D* rbfs_   = new RBF2D[n_];
    for( int i=0; i<n_; i++ ){
        //printf( "rbfs[%i]\n", i );
        rbfs_[i].pos.set( randf( xmin, xmax ), randf( ymin, ymax ) );
        //makeCliff( rbfs_+i );
        rbfs_[i].rmax   = randf( rmin, rmax );
        rbfs_[i].height = randf( vmin, vmax );
        //printf( "inserting rbfs[%i]\n to hashmap", i );
        //insertRBF( &rbfs[i] );
    }
    insertRBFs( n_, rbfs_ );
}


