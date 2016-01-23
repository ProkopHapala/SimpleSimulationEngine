
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "drawMath2D.h"

#include "AppSDL2OGL.h"
#include "testUtils.h"

#include "NBodyWorld2D.h"

class NBodyWorldApp : public AppSDL2OGL {
	public:
    NBodyWorld world;

	virtual void draw   ();
	virtual void drawHUD();
	NBodyWorldApp( int& id, int WIDTH_, int HEIGHT_ );

};

NBodyWorldApp::NBodyWorldApp( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {
    world.init();
}

void NBodyWorldApp::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    world.update();


/*
    glColor3f( 0.2f, 0.2f, 0.2f );
    for( int i=0; i<world.nParticles; i++ ){
        Particle2D* pi = &(world.particles[i]);
        ULONG icell = world.map.getBucket( pi->pos.x, pi->pos.y );
        double x,y;
        world.map.unfoldBucket( icell, x, y );
        Draw2D::drawRectangle( (float)x, (float)y, (float)(x+world.map.step), (float)(y+world.map.step), false );
    }
*/

/*
    glColor3f( 0.2f, 0.2f, 0.2f );
    for( ULONG icell : world.activeCells ){
        double x,y;
        world.map.unfoldBucket( icell, x, y );
        Draw2D::drawRectangle( (float)x, (float)y, (float)(x+world.map.step), (float)(y+world.map.step), false );
    }
*/

    glColor3f( 0.9f, 0.9f, 0.9f );
    Particle2D* screenObjects[65536];
    //float camXmin_ =-1; float camXmax_ =+1;
    //float camYmin_ =-1; float camYmax_ =+1;
    float camXmin_ =camXmin; float camXmax_ =camXmax;
    float camYmin_ =camYmin; float camYmax_ =camYmax;
/*
    Draw2D::drawRectangle( camXmin_, camYmin_, camXmax_, camYmax_, false );
	UINT nfound = world.map.getObjectsInRect( camXmin_, camYmin_, camXmax_, camYmax_, &(screenObjects[0]) );
	//glBegin(GL_POINTS);
	for( int i=0; i<nfound; i++ ){
        Particle2D* p = screenObjects[i];
		//glVertex3f( (float) p->pos.x, (float)p->pos.y, 1.0f );
		Draw2D::drawCircle_d( p->pos, 0.25, 8, true );
	}
	//glEnd();
	printf( "nfound %i filled %i \n", nfound, world.map.filled );
*/
    UHALF ix0 = world.map.getIx( camXmin_ );  UHALF iy0 = world.map.getIy( camYmin_ );
    UHALF ix1 = world.map.getIx( camXmax_ );  UHALF iy1 = world.map.getIy( camYmax_ );
    UINT nfound = 0;
    int ncells  = 0;
    for( UHALF iy = iy0; iy<=iy1; iy++ ){
        for( UHALF ix = ix0; ix<=ix1; ix++ ){
            UINT nfound = world.map.getBucketObjectsInt( ix, iy, screenObjects );
            if( nfound > 0 ){
                double x = world.map.getX(ix);
                double y = world.map.getY(iy);
                Draw2D::drawRectangle( (float)x, (float)y, (float)(x+world.map.step), (float)(y+world.map.step), false );
                for( int i=0; i<nfound; i++ ){
                    Particle2D* p = screenObjects[i];
                    Draw2D::drawCircle_d( p->pos, 0.25, 8, true );
                }
            }
            //printf( " ix %i iy %i  \n", ix, iy, ni );
        }
    }
    printf( " ======== frame %i DONE ( map.filled %i )\n", frameCount, world.map.filled );
	//STOP = true;


};

void NBodyWorldApp::drawHUD(){

}

// ===================== MAIN

NBodyWorldApp * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	thisApp = new NBodyWorldApp( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















