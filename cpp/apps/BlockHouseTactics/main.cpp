
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"
#include "Solids.h"

#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"

#include "BlockHouseWorld.h"

class TestAppBlockBuilder : public AppSDL2OGL_3D {
	public:

	BlockHouseWorld world;

    int ix=127,iy=127,iz=127;
    int cursorShape;
	int cursorSide=0;
    int cursorWallType=1;

	// ---- function declarations
	virtual void draw   ();
    //virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );

	void drawWall( const Vec3d& pos, int iSide, int shape ){
        glPushMatrix();
        Mat3d rotmat = world.rotations[ iSide ];
        rotmat.a.mul( world.scaling.x );
        rotmat.b.mul( world.scaling.y );
        rotmat.c.mul( world.scaling.z );
        float glMat[16];
        Draw3D::toGLMat( pos, rotmat, glMat );
        glMultMatrixf( glMat );
        glCallList( shape );
        glPopMatrix();
	}

	void drawBlock( const Block& block ){
	    Vec3d pos;
	    world.index2pos( {block.ix,block.iy,block.iz}, pos );
	    //printf( "block %i %i %i (%3.3f,%3.3f,%3.3f)", block.ix,block.iy,block.iz, pos.x, pos.y, pos.z );
		for(int iSide=0; iSide<6; iSide++){
            int type = block.sides[iSide];
            if(type<nMaxTypes){
                int shape = world.wallTypes[ type ].shape; // this seems to be a bit long dereferencing
                //printf( " side %i type %i shape %i \n", iSide, block.sides[iSide], shape );
                drawWall( pos, iSide, shape );
            }
		}
	}

	void drawTruss( bool DEBUG ){
	    if( ( world.truss.bonds != NULL )&&( world.truss.points != NULL )  ){
            glBegin( GL_LINES );
            for( int i=0; i<world.truss.nbonds; i++ ){
                Bond& bond = world.truss.bonds[i];
                Vec3d& pi  =  world.truss.points[bond.i];
                Vec3d& pj  =  world.truss.points[bond.j];
                glVertex3f( (float)pi.x, (float)pi.y, (float)pi.z );
                glVertex3f( (float)pj.x, (float)pj.y, (float)pj.z );
                if( DEBUG ){
                    printf( " %i  %i %i   (%3.3f,%3.3f,%3.3f)  (%3.3f,%3.3f,%3.3f)\n",  i,   bond.i, bond.j,  pi.x, pi.y, pi.z,  pj.x, pj.y, pj.z  );
                }
            }
            glEnd();
	    }
	}

	TestAppBlockBuilder( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppBlockBuilder::TestAppBlockBuilder( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    camPos.set( 0,0,-10 );

	world.init();

    cursorShape = glGenLists(1);
    glNewList( cursorShape , GL_COMPILE );
        glPushMatrix();
        glScalef( 0.5f, 0.5f, 0.5f );
        glColor3f(0.01f,0.01f,0.01f); Draw3D::drawLines    ( Solids::Cube_nedges, Solids::Cube_edges, Solids::Cube_verts );
        glPopMatrix();
    glEndList();

    world.wallTypes[world.nTypes].shape = glGenLists(1); // this is a bit strange, later we should initialize wallTypes more properly
    glNewList( world.wallTypes[world.nTypes].shape , GL_COMPILE );
        glBegin( GL_LINE_LOOP );
        glColor3f (  0.3f, 0.3f, 0.3f );
        glVertex3f(  0.5f,  0.5f, 0.51f);        glVertex3f(  0.5f, -0.5f, 0.51f);        glVertex3f( -0.5f, -0.5f, 0.51f);        glVertex3f( -0.5f,  0.5f, 0.51f);
        glEnd  ();
        glBegin( GL_QUADS );
		glColor3f (  0.2f, 0.2f, 0.2f );
        glNormal3f(  0.0f, 0.0f, 1.0f );
        glVertex3f(  0.5f,  0.5f, 0.5f);        glVertex3f(  0.5f, -0.5f, 0.5f);        glVertex3f( -0.5f, -0.5f, 0.5f);        glVertex3f( -0.5f,  0.5f, 0.5f);
        glEnd();
    glEndList();
    world.nTypes++;

    world.wallTypes[world.nTypes].shape = glGenLists(1);
    glNewList( world.wallTypes[world.nTypes].shape , GL_COMPILE );
        glBegin( GL_LINE_LOOP );
        glColor3f (  0.3f, 0.3f, 0.3f );
        glVertex3f(  0.5f,  0.5f, 0.51f);        glVertex3f(  0.5f, -0.5f, 0.51f);        glVertex3f( -0.5f, -0.5f, 0.51f);        glVertex3f( -0.5f,  0.5f, 0.51f);
        glEnd  ();
        glBegin( GL_QUADS );
		glColor3f (  0.8f, 0.2f, 0.2f );
        glNormal3f(  0.0f, 0.0f, 1.0f );
        glVertex3f(  0.5f,  0.5f, 0.5f);        glVertex3f(  0.5f, -0.5f, 0.5f);        glVertex3f( -0.5f, -0.5f, 0.5f);        glVertex3f( -0.5f,  0.5f, 0.5f);
        glEnd();
    glEndList();
    world.nTypes++;

    world.wallTypes[world.nTypes].shape = glGenLists(1);
    glNewList( world.wallTypes[world.nTypes].shape , GL_COMPILE );
        glBegin( GL_LINE_LOOP );
        glColor3f (  0.3f, 0.3f, 0.3f );
        glVertex3f(  0.5f,  0.5f, 0.51f);        glVertex3f(  0.5f, -0.5f, 0.51f);        glVertex3f( -0.5f, -0.5f, 0.51f);        glVertex3f( -0.5f,  0.5f, 0.51f);
        glEnd  ();
        glBegin( GL_QUADS );
		glColor3f (  0.2f, 0.8f, 0.2f );
        glNormal3f(  0.0f, 0.0f, 1.0f );
        glVertex3f(  0.5f,  0.5f, 0.5f);        glVertex3f(  0.5f, -0.5f, 0.5f);        glVertex3f( -0.5f, -0.5f, 0.5f);        glVertex3f( -0.5f,  0.5f, 0.5f);
        glEnd();
    glEndList();
    world.nTypes++;

    world.wallTypes[world.nTypes].shape = glGenLists(1);
    glNewList( world.wallTypes[world.nTypes].shape , GL_COMPILE );
        glBegin( GL_LINE_LOOP );
        glColor3f (  0.3f, 0.3f, 0.3f );
        glVertex3f(  0.5f,  0.5f, 0.5f);        glVertex3f(  0.5f, -0.5f, 0.5f);        glVertex3f( -0.5f, -0.5f, 0.5f);        glVertex3f( -0.5f,  0.5f, 0.5f);
        glEnd  ();
        glBegin( GL_QUADS );
		glColor3f (  0.2f, 0.2f, 0.8f );
        glNormal3f(  0.0f, 0.0f, 1.0f );
        glVertex3f(  0.5f,  0.5f, 0.5f);        glVertex3f(  0.5f, -0.5f, 0.5f);        glVertex3f( -0.5f, -0.5f, 0.5f);        glVertex3f( -0.5f,  0.5f, 0.5f);
        glEnd();
    glEndList();
    world.nTypes++;

};


void TestAppBlockBuilder::draw   (){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glEnable(GL_DEPTH_TEST);

    //glEnable     ( GL_LIGHTING         );

    //printf( " ==== frame %i \n", frameCount );
    //printf( " perspective %i first_person %i \n", perspective, first_person );
    //printf( " rot %3.3f %3.3f %3.3f \n", camMat.a.x, camMat.a.y, camMat.a.z );
    //printf( " rot %3.3f %3.3f %3.3f \n", camMat.b.x, camMat.b.y, camMat.b.z );
    //printf( " rot %3.3f %3.3f %3.3f \n", camMat.c.x, camMat.c.y, camMat.c.z );
    //printf( " %i %i %i   side %i type %i shape %i \n", ix,iy,iz, cursorSide, cursorWallType, world.wallTypes[ cursorWallType ].shape );


    Vec3d pos;
	world.index2pos( {ix, iy, iz}, pos );
    drawWall( pos, cursorSide, world.wallTypes[ cursorWallType ].shape  );

    //printf( "nBlocks %i nMaxBlocks %i \n", world.nBlocks, world.nMaxBlocks );
    for( int i=0; i<world.nBlocks; i++ ){
        if( !world.blocks[i].isEmpty() ){
            drawBlock( world.blocks[i] );
        }
    }

    glDisable ( GL_LIGHTING );

    glColor3f(0.5f,0.0f,0.0f); Draw3D::drawScale( {   0.0,   0.0,0.0}, {ix-127,0.0   ,   0.0}, {0.0,1.0,0.0}, 1.0, 0.1, 0.1 );
    glColor3f(0.0f,0.5f,0.0f); Draw3D::drawScale( {ix-127,   0.0,0.0}, {ix-127,iy-127,   0.0}, {1.0,0.0,0.0}, 1.0, 0.1, 0.1 );
    glColor3f(0.0f,0.0f,0.5f); Draw3D::drawScale( {ix-127,iy-127,0.0}, {ix-127,iy-127,iz-127}, {1.0,0.0,0.0}, 1.0, 0.1, 0.1 );
	Draw3D::drawAxis ( 3.0f );

	glDisable(GL_DEPTH_TEST);

	drawTruss( false );

	drawWall( pos, cursorSide, cursorShape );

};


void TestAppBlockBuilder::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    int rot;
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_a:  ix ++; if( ix >= world.nMax.x ) ix = 255; break;
                case SDLK_d:  ix --; if( ix <  0          ) ix = 0;   break;
                case SDLK_w:  iy ++; if( iy >= world.nMax.y ) iy = 255; break;
                case SDLK_s:  iy --; if( iy <  0          ) iy = 0;   break;
                case SDLK_q:  iz ++; if( iz >= world.nMax.z ) iz = 255; break;
                case SDLK_e:  iz --; if( iz <  0          ) iz = 0;   break;

                case SDLK_r:  cursorSide++; if( cursorSide >= 6 ) cursorSide = 0; break;
                case SDLK_f:  cursorSide--; if( cursorSide <  0 ) cursorSide = 5; break;
                case SDLK_LEFTBRACKET:  cursorWallType ++; if( cursorWallType>=world.nTypes ) cursorWallType = 0;            break;
                case SDLK_RIGHTBRACKET: cursorWallType --; if( cursorWallType< 0            ) cursorWallType = world.nTypes-1;    break;
                case SDLK_RETURN:    world.changeBlock( ix, iy, iz, cursorSide, cursorWallType ); break;
                case SDLK_BACKSPACE: world.eraseBlock ( ix, iy, iz ); break;
                case SDLK_u:  qCamera.setOne();  break;

                case SDLK_t: world.blocks2truss( ); drawTruss( true ); break;

                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            break;

        /*
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    pickParticle( world.picked );
                break;
            }
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
    AppSDL2OGL_3D::eventHandling( event );
}



// ===================== MAIN

TestAppBlockBuilder * testApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestAppBlockBuilder( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}
















