
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw.h"
#include "Draw2D.h"
#include "AppSDL2OGL.h"


#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"


#include "SimplexGrid.h"

#include "SimplexRuler.h"


// ======================  TestApp

using MySimplexField = SimplexField < bool, bool >;
using MySimplexGrid  = SimplexGrid  < MySimplexField >;

class TestAppSimplexGrid : public AppSDL2OGL{
	public:
    MySimplexGrid grid;

    SimplexRuler ruler;

    int shape;

    bool mouse_left = false;
    bool mouse_right = false;

    int nhits;
    Vec2d p1,p2;
    Vec2d hits      [1024];
    int   boundaries[1024];
    UHALF edges     [1024];


    Vec2d hray = (Vec2d){0.0,1.0};
    Vec2d ray0 = (Vec2d){0.6,0.6};

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

    ruler.MAP_OFFSET = 1000;
    grid.init( 1.0, 8 );

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
    nhits = grid.raster_line( dirHat, p1, p2, hits, boundaries, edges );

	glEndList();

	zoom = 5;

}

void TestAppSimplexGrid::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

    //return;

    if      ( mouse_left  ){ paintSimplex(mouse_begin_x,mouse_begin_y); }
    else if ( mouse_right ){ eraseSimplex(mouse_begin_x,mouse_begin_y); }

    if( mouse_left || mouse_right ){ printf( " frame %i grid.filled %i \n", frameCount, grid.filled ); };

    renderMapContent( );

	glCallList( shape );

	double da,db; UHALF ia,ib;

	//bool s = grid.simplexIndex( mouse_begin_x, mouse_begin_y, ia,ib, da, db );
    //renderSimplex( ia, ib, s, grid.step );

	Vec2i ind; Vec2d dind; Vec2d p;
	bool s = 1& ruler.simplexIndex( {mouse_begin_x, mouse_begin_y}, ind, dind );
	//printf( " (%i,%i) \n", ind.x, ind.y );
    ruler.nodePoint( ind, p );
    Draw2D::drawSimplex( p.x, p.y, s, ruler.step );

	glColor3f( 0.8f,0.8f,0.8f ); Draw2D::drawPointCross( {mouse_begin_x, mouse_begin_y}, 0.2f );

	/*
    glColor3f( 0.8, 0.0, 0.8 );
    Draw2D::drawLine_d( p1, p2 );
	for(int i=0; i<nhits; i++){
	    switch( boundaries[i] ){
            case 0: glColor3f( 1.0,0.0,0.0 ); break;
            case 1: glColor3f( 0.0,1.0,0.0 ); break;
            case 2: glColor3f( 0.0,0.0,1.0 ); break;
	    }
        Draw2D::drawPointCross_d( hits[i], 0.1 );

        Vec2d nd1,nd2;
        int ii = i<<2;
        grid.nodePoint( edges[ii]  , edges[ii+1], nd1.x, nd1.y );
        grid.nodePoint( edges[ii+2], edges[ii+3], nd2.x, nd2.y );
        Draw2D::drawLine_d( nd1, nd2 );
	}
    */

    printf("===============\n");
    //Vec2d hray = (Vec2d){ 0.0, 1.0};
    //Vec2d ray0 = (Vec2d){ 1.9, 0.6};
    Vec2d ray0 = (Vec2d){ mouse_begin_x, mouse_begin_y};
    hray.normalize();
    ruler.rayStart( ray0, hray );
    for(int i=0; i<10; i++){
        int edgeKind = ruler.rayStep();
        Draw::setRGB(0xFF<<(8*edgeKind));
        //glColor3f(0.0f,1.0f,0.0f);
        Draw2D::drawPointCross_d( ray0 + hray * ruler.ray_t, 0.1 );
        Vec2d p1,p2;
        //ruler.nodePoint ( {ruler.ray_i.a+kind2edge[edgeKind][0], ruler.ray_i.b+kind2edge[edgeKind][1]}, p1 );
        //ruler.nodePoint ( {ruler.ray_i.a+kind2edge[edgeKind][2], ruler.ray_i.b+kind2edge[edgeKind][3]}, p2 );
        //printf( "%i (%i,%i)(%i,%i)\n", edgeKind, kind2edge[edgeKind][0].x, kind2edge[edgeKind][0].y, kind2edge[edgeKind][1].x, kind2edge[edgeKind][1].y );
        Vec2i ip1,ip2; ip1 = ruler.ray_i.xy() + kind2edge[edgeKind][0]; ip2 = ruler.ray_i.xy() + kind2edge[edgeKind][1];

        Draw::color_of_hash(edgeKind*154);
        ruler.nodePoint ( ip1, p1 );
        ruler.nodePoint ( ip2, p2 );
        //printf( "%i (%i,%i)(%i,%i) (%g,%g)(%g,%g) \n", edgeKind, ip1.x, ip1.y, ip2.x, ip2.y, p1.x, p1.y, p2.x, p2.y );
        Draw2D::drawLine_d( p1, p2 );
        //printf("%f \n", ruler.ray_t );
    }
    glColor3f(1.0f,1.0f,1.0f);
    Draw2D::drawPointCross_d( ray0, 0.3 );
    Draw2D::drawLine_d      ( ray0, ray0+hray*ruler.ray_t );

    //STOP = true;
};

void TestAppSimplexGrid::drawHUD(){
}


void TestAppSimplexGrid::eventHandling( const SDL_Event& event ){
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case   SDLK_r: hray=(Vec2d){randf(-5.0,5.0),randf(-5.0,5.0)}; ray0=(Vec2d){randf(-5.0,5.0),randf(-5.0,5.0)}; break;
            }
            break;
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
















