
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"
#include "drawMath2D.h"
#include "AppSDL2OGL.h"

#include "SimplexGrid.h"

// ======================  TestApp

using MySimplexField = SimplexField < bool, bool >;
using MySimplexGrid  = SimplexGrid  < MySimplexField >;

class TestAppSimplexGrid : public AppSDL2OGL{
	public:
    MySimplexGrid grid;
    int shape;

	// ---- function declarations

	virtual void draw   ();
	virtual void drawHUD();
    virtual void eventHandling( const SDL_Event& event );
	//virtual int tileToList( float x0, float y0, float x1, float y1 );

	TestAppSimplexGrid( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppSimplexGrid::TestAppSimplexGrid( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

    grid.init( 1.0, 8 );

    MySimplexField* p;
    p = new MySimplexField(); p->lo = true; p->hi = true; grid.insertNoTest( p, 1.15454, 2.15454 );
    p = new MySimplexField(); p->lo = true; p->hi = false; grid.insertNoTest( p, 3.15454, 2.15454 );
    p = new MySimplexField(); p->lo = false; p->hi = true; grid.insertNoTest( p, 1.15454, 5.15454 );


    shape=glGenLists(1);
	glNewList( shape, GL_COMPILE );
	glBegin   ( GL_POINTS   );
        for( int i=0; i<10000; i++ ){
            double x = randf( 0,5 );
            double y = randf( 0,5 );
            double da,db;
            UHALF   ia,ib;
            //bool s = simplexIndex( x+100, y+100, ia,ib, da, db );
            bool s = grid.simplexIndex( x+100, y+100, ia,ib, da, db );
            //printf( " (%3.3f,%3.3f) (%3.3f,%3.3f) (%i,%i) \n", x, y, da, db, ia, ib );
            s = false;
            if( s ){
                glColor3f( da, db, 1-da-db );
            }else{
                //int h = rand_hash( (ia*359) ^ (ib*353)  );
                int h = (ia*920419823) ^ (ib*941083981);
                glColor3f( (h&0x0000FF)/255.0f, (h&0x00FF00)/65535.0f, (h&0xFF0000)/16777216.0f );
            }
            glVertex3f( (float)x,(float)y, 0.0f );
        }
    glEnd();
	glEndList();

}

void TestAppSimplexGrid::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

    glColor3f( 0.1f,0.1f,0.1f );
    glBegin( GL_TRIANGLES );
	for( int i=0; i<grid.capacity; i++ ){
        MySimplexField * field = grid.fields[ i ].object;
        //printf( " %i : object %i bucket %i \n", i, field, grid.fields[ i ].bucket );
        if ( field != NULL ){
            UHALF  ia,ib; grid.unfoldBucketInt( grid.fields[ i ].bucket, ia, ib );
            double x,y;   grid.nodePoint( ia, ib, x, y );
            if( field->lo ){
                glVertex3f( (float)x,       (float)y,                 0.0f );
                glVertex3f( (float)(x+1),   (float)y,                 0.0f );
                glVertex3f( (float)(x+0.5), (float)(y+0.86602540378), 0.0f );
            };
            if( field->hi ){
                glVertex3f( (float)(x+0.5 ), (float)(y+0.86602540378),  0.0f );
                glVertex3f( (float)(x+1.5 ), (float)(y+0.86602540378),  0.0f );
                glVertex3f( (float)(x+1.0 ), (float) y               ,  0.0f );
            };
            //printf( " %i %i %3.3f %3.3f  %i %i \n", ia, ib, x, y, field->lo, field->hi );
        }
	}
	glEnd();

    //exit(0);
	//glCallList( shape );

};


/*
int TestAppSimplexGrid::tileToList( float x0, float y0, float x1, float y1 ){
    int ilist=glGenLists(1);
    glNewList( ilist, GL_COMPILE );
        terrain.renderRect( x0, y0, x1, y1, 9 );
        glColor3f(0.9f,0.2f,0.2f); Draw2D::drawRectangle( x0+0.1, y0+0.1, x1-0.1, y1-0.1, false );
    glEndList();
    return ilist;
}
*/

void TestAppSimplexGrid::drawHUD(){
}


void TestAppSimplexGrid::eventHandling( const SDL_Event& event ){
    switch( event.type ){
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    MySimplexField* p = new MySimplexField();
                    p->lo = true;
                    p->hi = true;
                    grid.insertIfNew( p, mouse_begin_x, mouse_begin_y );
                break;
            }
            break;
        /*
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

TestAppSimplexGrid * testApp;

int main(int argc, char *argv[]){

	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestAppSimplexGrid( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}
















