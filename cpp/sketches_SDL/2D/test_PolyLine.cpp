
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "AppSDL2OGL.h"
#include "testUtils.h"

#include "SoftPolyLine2D.h"


Vec2d stringFroce(Vec2d p, Vec2d p0, double l0, double k){
    Vec2d d  = p0-p;
    double l2 = d.norm2();
    if( l2>(l0*l0) ){
        double l = sqrt(l2);
        d.mul( k*(l-l0)/l );
    }else{
        d.set(0.0);
    }
    return d;
}



class TestAppPolyLine : public AppSDL2OGL { public:
    SoftPolyLine2D pline1;
    int ipick=-1;

	virtual void draw   ();
	virtual void drawHUD();
	virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );

	TestAppPolyLine( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppPolyLine::TestAppPolyLine( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

    pline1.init( 5, {-10,0.0}, {10.0, 0.0}, 1.0, 1.0, 2.0 );

    pline1.ps[2].y += 4.0;
    printf("=====\n");
    pline1.setCurrentRef();
    //for(int i=0; i<pline1.n; i++){ pline1.ps[i].y+=randf(-3.0,3.0); }

    //pline1.ang0s[0].fromAngle( -0.5 );

}

void TestAppPolyLine::draw(){
    long tstart = getCPUticks();
    glClearColor( 0.9f, 0.9f, 0.9f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );


	pline1.cleanForce();
	pline1.updateAux();
	pline1.radialForces();
	pline1.angularForces();
	if(ipick>=0){ pline1.fs[ipick].add( stringFroce(pline1.ps[ipick], {mouse_begin_x,mouse_begin_y}, 0.0, 1.0) );  };
	//pline1.fs[0].set(0.0);
	//pline1.fs[1].set(0.0);
	//pline1.fs[2].set(0.0);
	pline1.move( 0.1, 0.02 );
	//exit(0);


	pline1.updateAux();
	Vec2d pm = (Vec2d){mouse_begin_x,mouse_begin_y};
	Vec2d dp = pline1.dpmin( pm );
	glColor3f(0.0f,0.8f,0.0f); Draw2D::drawVecInPos_d( dp*(10.0/dp.norm2()), pm );
	//printf( "dp (%f,%f)\n", dp.x, dp.y );

	glColor3f(0.0f,1.0f,0.0f);
	for(int iy=-20; iy<20; iy++){
		for(int ix=-50; ix<50; ix++){
            Vec2d pm = (Vec2d){ix*0.4,iy*0.5};
            Vec2d dp = pline1.dpmin( pm );
            float r = dp.norm();
            if(r<2.0){
                Draw2D::drawVecInPos_d( dp*(0.3/dp.norm()), pm );
            }
        }
	}

	//for(int i=0; i<pline1.n; i++){ printf( "%i f=(%f,%f) \n", i, pline1.fs[i].x, pline1.fs[i].y ); }
	glColor3f(0.0f,0.0f,0.0f); Draw2D::drawLines( pline1.n, pline1.ps );
	glColor3f(1.0f,0.0f,0.0f); Draw2D::drawLines( pline1.n, pline1.ps, pline1.fs, 10.0 );

};

void TestAppPolyLine::mouseHandling( ){
    Uint32 buttons = SDL_GetMouseState( &mouseX, &mouseY );
    defaultMouseHandling( mouseX, mouseY );
}

void TestAppPolyLine::eventHandling ( const SDL_Event& event  ){
    //printf( "TestAppPolyLine::eventHandling() \n" );
    switch( event.type ){
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    //pickParticle( world.picked );
                    ipick = pline1.nearestPoint( {mouse_begin_x,mouse_begin_y} );
                break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    //world.picked = NULL;
                    ipick = -1;
                    break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}

void TestAppPolyLine::drawHUD(){

}

// ===================== MAIN

TestAppPolyLine * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	thisApp = new TestAppPolyLine( junk , 800, 600 );
	thisApp->zoom = 30;
	thisApp->loop( 1000000 );
	return 0;
}
















