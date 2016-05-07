
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "AppSDL2OGL.h"


#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"


#include "SimplexGrid.h"

// ======================  TestApp

using MySimplexField = SimplexField < bool, bool >;
using MySimplexGrid  = SimplexGrid  < MySimplexField >;

class TestAppSimplexGrid : public AppSDL2OGL{
	public:
    MySimplexGrid grid;
    int shape;

    bool mouse_left = false;
    bool mouse_right = false;

    int nhits;
    Vec2d p1,p2;
    Vec2d hits[1024];
    int   boundaries[1024];

	// ---- function declarations

	virtual void draw   ();
	virtual void drawHUD();
    virtual void eventHandling( const SDL_Event& event );
	//virtual int tileToList( float x0, float y0, float x1, float y1 );

    void paintSimplex( double x, double y );
    void eraseSimplex( double x, double y );

    void renderMapContent ( );
	void renderSimplex    ( int ia, int ib, bool s, float step );
	//void drawSimplexGrid( int n, float step );

	TestAppSimplexGrid( int& id, int WIDTH_, int HEIGHT_ );

};

void TestAppSimplexGrid::paintSimplex( double x, double y ){
	double da,db; UHALF ia,ib;
	bool  s      = grid.simplexIndex        ( mouse_begin_x, mouse_begin_y, ia,ib, da, db );
    ULONG bucket = grid.getBucketInt        ( ia, ib );
    int index  = grid.getFirstBucketIndex ( bucket );
    MySimplexField* p;
    if( index >= 0 ){
        //printf( " modifying existing tile, index %i \n", index );
        p = grid.fields[ index ].object;
        if( s ){ p->hi=true; }else{ p->lo=true; };
    }else{
        printf( " inserted (%i,%i) index %i hash %i bucket %i \n", ia, ib, index, grid.mask & hashFunc( bucket ), bucket );
        //printf( " inserting new tile \n" );
        p = new MySimplexField();
        if( s ){ p->lo=false; p->hi=true; }else{ p->lo=true; p->hi=false; };
        grid.HashMap<MySimplexField>::insertNoTest( p, bucket );
    }
};

void TestAppSimplexGrid::eraseSimplex( double x, double y ){
	double da,db; UHALF ia,ib;
	bool  s      = grid.simplexIndex        ( mouse_begin_x, mouse_begin_y, ia,ib, da, db );
    ULONG bucket = grid.getBucketInt        ( ia, ib );
    int index  = grid.getFirstBucketIndex   ( bucket );
    MySimplexField* p;
    if( index >= 0 ){
        //printf( " errasing tile, index %i \n", index );
        p = grid.fields[ index ].object;
        if( s ){ p->hi=false; }else{ p->lo=false; };
        if( !( p->hi || p->lo ) ){
            printf( " removed (%i,%i) index %i hash %i bucket %i \n", ia, ib, index, grid.mask & hashFunc( bucket ), bucket );
            //printf( " empty field %i removed from HashMap \n", index );
            grid.HashMap<MySimplexField>::remove( index, grid.mask & hashFunc( bucket ) );
            delete p;
        }
    };
};

void TestAppSimplexGrid::renderSimplex( int ia, int ib, bool s, float step ){
    double x,y;
    grid.nodePoint( ia, ib, x, y );
    int h   = (ia*920419823) ^ (ib*941083981);
    glColor3f( (h&0x0000FF)/255.0f, (h&0x00FF00)/65535.0f, (h&0xFF0000)/16777216.0f );
    Draw2D::drawSimplex( (float)x, (float)y, s, step );
}

void TestAppSimplexGrid::renderMapContent( ){
    glColor3f( 0.1f,0.1f,0.1f );
    glBegin( GL_TRIANGLES );
	for( int i=0; i<grid.capacity; i++ ){
        MySimplexField * field = grid.fields[ i ].object;
        //printf( " %i : object %i bucket %i \n", i, field, grid.fields[ i ].bucket );
        if ( field != NULL ){
            UHALF  ia,ib; grid.unfoldBucketInt( grid.fields[ i ].bucket, ia, ib );
            double x,y;   grid.nodePoint( ia, ib, x, y );
            if( field->lo ){
                glVertex3f( (float)(x              ), (float)y,                           0.0f );
                glVertex3f( (float)(x+    grid.step), (float)y,                           0.0f );
                glVertex3f( (float)(x+0.5*grid.step), (float)(y+0.86602540378*grid.step), 0.0f );
            };
            if( field->hi ){
                glVertex3f( (float)(x+0.5*grid.step ), (float)(y+0.86602540378*grid.step), 0.0f );
                glVertex3f( (float)(x+1.5*grid.step ), (float)(y+0.86602540378*grid.step), 0.0f );
                glVertex3f( (float)(x+1.0*grid.step ), (float) y                         , 0.0f );
            };
            //printf( " %i %i %3.3f %3.3f  %i %i \n", ia, ib, x, y, field->lo, field->hi );
        }
	}
	glEnd();
}

TestAppSimplexGrid::TestAppSimplexGrid( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

    grid.init( 1.7, 8 );

    shape=glGenLists(1);
	glNewList( shape, GL_COMPILE );

	glBegin   ( GL_POINTS   );
        for( int i=0; i<10000; i++ ){
            double x = randf( -5,5 );
            double y = randf( -5,5 );
            double da,db;
            UHALF   ia,ib;
            //bool s = simplexIndex( x+100, y+100, ia,ib, da, db );
            bool s = grid.simplexIndex( x, y, ia,ib, da, db );
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

    glColor3f(0.1f,0.1f,0.1f); Draw2D::drawSimplexGrid( 10, (float)grid.step );

    Vec2d dirHat;
    p1.set(0.9,0.6);
    p2.set(8.2,8.3);
    dirHat.set_sub( p2, p1 ); dirHat.normalize();
    nhits = grid.raster_line( dirHat, p1, p2, hits, boundaries );

	glEndList();

}

void TestAppSimplexGrid::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );


    if      ( mouse_left  ){ paintSimplex(mouse_begin_x,mouse_begin_y); }
    else if ( mouse_right ){ eraseSimplex(mouse_begin_x,mouse_begin_y); }

    if( mouse_left || mouse_right ){ printf( " frame %i grid.filled %i \n", frameCount, grid.filled ); };

    renderMapContent( );

	glCallList( shape );

	double da,db; UHALF ia,ib;
	bool s = grid.simplexIndex( mouse_begin_x, mouse_begin_y, ia,ib, da, db );
	renderSimplex( ia, ib, s, grid.step );
	glColor3f( 0.8f,0.8f,0.8f ); Draw2D::drawPointCross( {mouse_begin_x, mouse_begin_y}, 0.2f );

    glColor3f( 0.8, 0.0, 0.8 );
    Draw2D::drawLine_d( p1, p2 );
	for(int i=0; i<nhits; i++){
	    switch( boundaries[i] ){
            case 0: glColor3f( 1.0,0.0,0.0 ); break;
            case 1: glColor3f( 0.0,1.0,0.0 ); break;
            case 2: glColor3f( 0.0,0.0,1.0 ); break;
	    }
        Draw2D::drawPointCross_d( hits[i], 0.1 );
	}

};

void TestAppSimplexGrid::drawHUD(){
}


void TestAppSimplexGrid::eventHandling( const SDL_Event& event ){
    switch( event.type ){
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //paintSimplex( mouse_begin_x, mouse_begin_y );
                    mouse_left  = true;
                    break;
                case SDL_BUTTON_RIGHT:
                    mouse_right = true;
                    //eraseSimplex( mouse_begin_x, mouse_begin_y );
                    break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    mouse_left = false;
                    break;
                case SDL_BUTTON_RIGHT:
                    mouse_right = false;
                    break;
            }
            break;
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
















