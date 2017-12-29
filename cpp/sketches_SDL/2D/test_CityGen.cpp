
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw.h"
#include "Draw2D.h"
#include "AppSDL2OGL.h"

#include "CityGeneration.h"
//#include "TerrainCubic.h"
//#include "TiledView.h"

void drawQuadNode( QuadNode& qnd, std::vector<Vec2d>& ps ){
    glBegin(GL_LINE_LOOP);
        Vec2d p;
        p=ps[qnd.i00]; glVertex3f(p.x,p.y,0.0); //printf( " %f %f \n", p.x, p.y );
        p=ps[qnd.i01]; glVertex3f(p.x,p.y,0.0); //printf( " %f %f \n", p.x, p.y );
        p=ps[qnd.i11]; glVertex3f(p.x,p.y,0.0); //printf( " %f %f \n", p.x, p.y );
        p=ps[qnd.i10]; glVertex3f(p.x,p.y,0.0); //printf( " %f %f \n", p.x, p.y );
    glEnd();
    //exit(0);
}

void drawQuadRecur( QuadNode& qnd, std::vector<Vec2d>& ps ){
    if( qnd.subs ){
        int ntot = qnd.n.totprod();
        for(int i=0; i<ntot; i++){
            drawQuadRecur( qnd.subs[i], ps );
        };
    }else{
        //Draw::color_of_hash( rand() );
        drawQuadNode( qnd, ps );
    };
}

class TestAppCityGen : public AppSDL2OGL {
	public:
	int viewList;

	std::vector<Vec2d> ps;

	QuadNode rootQuad;

	// ---- function declarations

	virtual void draw   ();
	virtual void drawHUD();

	TestAppCityGen( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppCityGen::TestAppCityGen( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {
    double sz = 4.0;
    double rndsz = 1.0;
    srand(15454587);
    ps.push_back( {-sz+randf(-rndsz,rndsz), -sz+randf(-rndsz,rndsz) } );
    ps.push_back( {-sz+randf(-rndsz,rndsz),  sz+randf(-rndsz,rndsz) } );
    ps.push_back( { sz+randf(-rndsz,rndsz), -sz+randf(-rndsz,rndsz) } );
    ps.push_back( { sz+randf(-rndsz,rndsz),  sz+randf(-rndsz,rndsz) } );
    rootQuad.setCorners(0,1,2,3);

    //rootQuad.inset( {0.1,0.1},{0.9,0.9}, ps );
    //rootQuad.split( {4,3}, {0.5,0.5} );
    //rootQuad.makeSubPoints( ps );
    //rootQuad.makeSubQuads();

    rootQuad.splitRecursive( 3, (Vec2i){2,2}, (Vec2i){3,3}, (Vec2d){0.3,0.3}, 0.1, ps );
}

void TestAppCityGen::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

	//drawQuadNode( rootQuad, ps );

	srand(15454587);
	glColor3f(0.0f,0.0f,0.0f);
    drawQuadRecur( rootQuad, ps );

    glColor3f(1.0f,1.0f,1.0f);
	glBegin(GL_POINTS);
	for( Vec2d& p : ps ){
        glVertex3f( p.x, p.y, 0.0 );
	}
	glEnd();
};

void TestAppCityGen::drawHUD(){}



// ===================== MAIN

TestAppCityGen * testApp;

int main(int argc, char *argv[]){

	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestAppCityGen( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}
















