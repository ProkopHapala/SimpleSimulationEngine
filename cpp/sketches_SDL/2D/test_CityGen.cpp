/// @file @brief This program demonstrates procedural 2D city generation. It creates a city layout, including road networks, building plots, and potentially simple building representations, allowing the user to explore different urban patterns generated algorithmically.


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
#include "CityGeneration2.h"
//#include "TerrainCubic.h"
//#include "TiledView.h"




template<typename Func>
class BinSearch{ public:
    Func f;
    double tol = 1e-8;

    BinSearch(Func f_):f(f_){};
    double exec( double xmin, double xmax ){
        double x = (xmin+xmax)*0.5;
        //printf( "to f() \n" );
        double y = f(x);
        if( fabs(y)<tol) return x;
        if( y > 0 ){ return exec(xmin,x); }
        else       { return exec(x,xmax); }
    }
};



void drawQuad( Quad2d& q){
    glBegin(GL_LINE_LOOP);
        glVertex3f(q.p00.x,q.p00.y, 10.0);
        glVertex3f(q.p01.x,q.p01.y, 10.0);
        glVertex3f(q.p11.x,q.p11.y, 10.0);
        glVertex3f(q.p10.x,q.p10.y, 10.0);
    glEnd();
}


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

	std::vector<Quad2d> quads;

	QuadNode rootQuad;

	QuadNode2 q2;

	// ---- function declarations

	virtual void draw   ();
	virtual void drawHUD();

	TestAppCityGen( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppCityGen::TestAppCityGen( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {
    double sz = 4.0;
    double rndsz = 1.0;
    srand(1544579);
    ps.push_back( {-sz+randf(-rndsz,rndsz), -sz+randf(-rndsz,rndsz) } );
    ps.push_back( {-sz+randf(-rndsz,rndsz),  sz+randf(-rndsz,rndsz) } );
    ps.push_back( { sz+randf(-rndsz,rndsz), -sz+randf(-rndsz,rndsz) } );
    ps.push_back( { sz+randf(-rndsz,rndsz),  sz+randf(-rndsz,rndsz) } );
    rootQuad.setCorners(0,1,2,3);

    //rootQuad.inset( {0.1,0.1},{0.9,0.9}, ps );
    //rootQuad.split( {4,3}, {0.5,0.5} );
    //rootQuad.makeSubPoints( ps );
    //rootQuad.makeSubQuads();

    //rootQuad.splitRecursive( 3, (Vec2i){2,2}, (Vec2i){3,3}, (Vec2d){0.3,0.3}, 0.1, ps );

    //splitOpen( q2, 4, 0b111, ps );

    /*
    int iter=0;
    auto f = [&](double x){ iter++; double y=x*x*x; printf( "i,x,y %i %f %f \n", iter, x, y ); return y;  };
    BinSearch<decltype(f)> binsearch(f);
    binsearch.exec(-15.1,2.5);
    printf( "final iter %i \n", iter );
    */

    int added=0;
    auto fcond = [&](const Quad2d& q){ return true; };
    auto fleaf = [&](const Quad2d& q){ printf("in fleaf \n" ); added++; quads.push_back(q); };
    QuadSpliter<decltype(fcond),decltype(fleaf)> quadspliter( fcond, fleaf );
    //QuadSpliter<auto,auto> quadspliter( fcond, fleaf );

    Quad2d q; q.set( ps[0],ps[1],ps[2],ps[3] );
    quadspliter.splitBinRec( q, 4, 0);
    printf( "added = %i\n", added );


}

void TestAppCityGen::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

	//drawQuadNode( rootQuad, ps );

	//drawQuadNode2Rec( q2, ps );

	for( Quad2d& q : quads ){ drawQuad( q); }

    /*
	srand(15454587);
	glColor3f(0.0f,0.0f,0.0f);
    drawQuadRecur( rootQuad, ps );

    glColor3f(1.0f,1.0f,1.0f);
	glBegin(GL_POINTS);
	for( Vec2d& p : ps ){
        glVertex3f( p.x, p.y, 0.0 );
	}
	glEnd();
	*/
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
















