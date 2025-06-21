/// @file @brief This demo explores the conversion between triangular and rectangular data representations in 2D, using `Tris2Rect.h`. It visualizes how data from a triangular mesh can be mapped onto a regular rectangular grid, or vice-versa, which is useful for various data processing tasks.

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "geom2D.h"

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw.h"
#include "Draw2D.h"
#include "AppSDL2OGL.h"
#include "Ruler2DFast.h"

//#include "TerrainRBF.h"
#include "TiledView.h"

// ======================  TestApp
// note
// GridMap2D contains function to include LineSegment and Trinagle into rectangular map


void rectGrid2TriMesh( Vec2i n, Vec2d* points, std::vector<Vec3i>& tris ){
    int i=0;
    for( int iy=1; iy<n.y; iy++ ){
        for( int ix=1; ix<n.x; ix++ ){
            Vec2d d;
            int i00=i;     int i01=i+1;
            int i10=i+n.x; int i11=i+n.x+1;
            d.set_sub( points[i00], points[i11] ); double rr1 = d.norm2();
            d.set_sub( points[i10], points[i01] ); double rr2 = d.norm2();
            if( rr1 < rr2 ){          //  \     //
                tris.push_back({i00,i10,i11});
                tris.push_back({i00,i11,i01});
                //printf( "%i %i (%i,%i,%i) (%i,%i,%i)\n", iy,ix,  i00,i10,i11, i00,i11,i01 );
            }else{                    //  /    //
                tris.push_back({i00,i10,i01});
                tris.push_back({i11,i01,i10});
            }
            i++;
        }
        i++;
    };
}



class TestAppTris2Rect : public AppSDL2OGL, public TiledView {
	public:
	int npoints;
	//RBF2D * rbfs;
	//TerrainRBF terrain;
	int viewList;

	bool clicked = false;

	Vec2i   Npoints;
	Vec2d * points  = NULL;
	double * heights = NULL;
	std::vector<Vec3i> tris;

	// ---- function declarations

	virtual void draw   ();
	virtual void drawHUD();
    virtual void eventHandling( const SDL_Event& event );
	virtual int tileToList( float x0, float y0, float x1, float y1 );

	TestAppTris2Rect( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppTris2Rect::TestAppTris2Rect( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

    Npoints = {16,16};
    points  = new Vec2d[Npoints.x*Npoints.y];
    heights = new double[Npoints.x*Npoints.y];

    int i=0;
    Vec2d p;
    for( int iy=0; iy<Npoints.y; iy++ ){
        for( int ix=0; ix<Npoints.x; ix++ ){
            p.x=ix*1.0+randf(-0.4,0.4);
            p.y=iy*1.0+randf(-0.4,0.4);
            points[i] = p;
            heights[i] = randf(0.0,1.0);
            i++;
        }
    }

    rectGrid2TriMesh( Npoints, points, tris );

    //viewList = tileToList( -10, -10, 10, 10  );

    //TiledView::init( 6, 6 ); // initialize parent
    //x0 -= 1000;
    //y0 -= 1000;
    //tiles    = new int[ nxy ];

    //TiledView::renderAll( -10, -10, 10, 10 );
    //exit(0);

}

void TestAppTris2Rect::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

    /*
	for( Vec3i tri : tris ){
        Draw::color_of_hash( (tri.a+54)*(tri.b+87)*(tri.c+59) );
        Draw2D::drawTriangle_d( points[tri.a], points[tri.b], points[tri.c] );
	}
    */


    glBegin(GL_TRIANGLES);
    double c; Vec2d p;
    for( Vec3i tri : tris ){
        c = heights[tri.a]; glColor3f(c,c,c); p=points[tri.a]; glVertex3f( (float)p.x, (float)p.y, 100.0);
        c = heights[tri.b]; glColor3f(c,c,c); p=points[tri.b]; glVertex3f( (float)p.x, (float)p.y, 100.0);
        c = heights[tri.c]; glColor3f(c,c,c); p=points[tri.c]; glVertex3f( (float)p.x, (float)p.y, 100.0);
	}
	glEnd();


    //terrain.renderRect( 3, 3, 8, 8, 30 );
    //if( viewList != 0 ) glCallList( viewList );

    //exit(0);

    //TiledView::render(  3, 3, 8, 8 );

    /*
    float szview = 4.0f;
    TiledView::draw(  mouse_begin_x-szview, mouse_begin_y-szview, mouse_begin_x+szview, mouse_begin_y+szview );
    glColor3f(0.2f,0.2f,0.9f); Draw2D::drawRectangle( -10, -10, 10, 10, false );
    glColor3f(0.9f,0.2f,0.9f); Draw2D::drawRectangle( x0, y0, x0+nx*step, y0+ny*step, false );
    //exit(0);
    */


};


int TestAppTris2Rect::tileToList( float x0, float y0, float x1, float y1 ){
    int ilist=glGenLists(1);
    glNewList( ilist, GL_COMPILE );
        //terrain.renderRect( x0, y0, x1, y1, 9 );
        //glColor3f(0.9f,0.2f,0.2f); Draw2D::drawRectangle( x0+0.1, y0+0.1, x1-0.1, y1-0.1, false );
    glEndList();
    return ilist;
}

void TestAppTris2Rect::drawHUD(){
}


void TestAppTris2Rect::eventHandling( const SDL_Event& event ){
    switch( event.type ){
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    clicked = true;
                break;
            }
            /*
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    world.picked = NULL;
                    break;
            }
            break;
            */
    };
    AppSDL2OGL::eventHandling( event );
};

// ===================== MAIN

TestAppTris2Rect * testApp;

int main(int argc, char *argv[]){




	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestAppTris2Rect( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}
