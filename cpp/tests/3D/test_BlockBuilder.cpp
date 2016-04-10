
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
#include "Solids.h"

#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"

namespace BlockWorld {

    static Vec3d scaling;
    static Vec3d pos0;
    static Mat3d rotations[6];

    void drawBlock( Uint16 shape, Uint8 orientation, Uint8 ix, Uint8 iy, Uint8 iz ){
        Mat3d rotmat = rotations[ orientation & 0b00011111 ];
        rotmat.a.mul( scaling.x );
        rotmat.b.mul( scaling.y );
        rotmat.c.mul( scaling.z );
        if (orientation & 0b10000000) rotmat.a.mul( -1.0d );
        if (orientation & 0b01000000) rotmat.b.mul( -1.0d );
        if (orientation & 0b00100000) rotmat.c.mul( -1.0d );
        Vec3d pos; pos.set( ix*scaling.x+pos0.x, iy*scaling.y+pos0.y, iz*scaling.z+pos0.z );
        float glMat[16];
        glPushMatrix ( );
        /*
        printf( " pos (%3.3f,%3.3f,%3.3f) \n", pos.x, pos.y, pos.z );
        printf( " a   (%3.3f,%3.3f,%3.3f) \n", rotmat.a.x, rotmat.a.x, rotmat.a.x );
        printf( " b   (%3.3f,%3.3f,%3.3f) \n", rotmat.b.x, rotmat.b.x, rotmat.b.x );
        printf( " c   (%3.3f,%3.3f,%3.3f) \n", rotmat.c.x, rotmat.c.x, rotmat.c.x );
        */
        Draw3D::toGLMat( pos, rotmat, glMat );
        glMultMatrixf( glMat );
        glCallList  ( shape );
        glPopMatrix ( );
    }

    void buildRotations( const Mat3d& rot ){
        rotations[0].set( rot.a, rot.b, rot.c );
        rotations[1].set( rot.a, rot.c, rot.b );
        rotations[2].set( rot.b, rot.a, rot.c );
        rotations[3].set( rot.b, rot.c, rot.a );
        rotations[4].set( rot.c, rot.a, rot.b );
        rotations[5].set( rot.c, rot.b, rot.a );
    }

    void setupBlockWorld( const Vec3d& pos0_, const Vec3d& scaling_, const Mat3d& rot ){
        pos0   .set( pos0_    );
        scaling.set( scaling_ );
        buildRotations( rot );
        /*
        printf( " pos      (%3.3f,%3.3f,%3.3f) \n", pos0.x, pos0.y, pos0.z );
        //printf( " scaling_ (%3.3f,%3.3f,%3.3f) \n", scaling_.x, scaling_.y, scaling_.z );
        printf( " scaling  (%3.3f,%3.3f,%3.3f) \n", scaling.x,  scaling.y, scaling.z );
        printf( " a        (%3.3f,%3.3f,%3.3f) \n", rotations[0].a.x, rotations[0].a.y, rotations[0].a.z );
        printf( " b        (%3.3f,%3.3f,%3.3f) \n", rotations[0].b.x, rotations[0].b.y, rotations[0].b.z );
        printf( " c        (%3.3f,%3.3f,%3.3f) \n", rotations[0].c.x, rotations[0].c.y, rotations[0].c.z );
        */
    }

};


class Block{
    public:
    Uint16 shape;
    Uint8  orientation;
    Uint8  ix;
    Uint8  iy;
    Uint8  iz;
};


// ======================  TestApp

class TestAppBlockBuilder : public AppSDL2OGL_3D {
	public:

    const static int nMaxBlocks = 1024;
    int nBlocks = 0;
    int iBlock  = 0;
	Block blocks[nMaxBlocks];

	const static int nMaxShapes = 256;
	int nShapes = 0;
	int iShape  = 0;
    int shapes[nMaxShapes];

    Uint8 ix=127,iy=127,iz=127,orientation=0;
    int cursorShape;

	// ---- function declarations
	virtual void draw   ();
    //virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );

	int findBlock        ( Uint8 ix, Uint8 iy, Uint8 iz );
	int changeBlock      ( Uint8 ix, Uint8 iy, Uint8 iz, Uint8 orientation, Uint16 shape );
	int defragmentBlocks ( );

	TestAppBlockBuilder( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppBlockBuilder::TestAppBlockBuilder( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

	float l0_pos  [] = {0.1f,0.2f,1.0f,0}; glLightfv    ( GL_LIGHT0, GL_POSITION, l0_pos  );
	float l0_amb  [] = {0.0f,0.0f,0.0f,0}; glLightfv    ( GL_LIGHT0, GL_AMBIENT,  l0_amb  );
	float l0_diff [] = {0.4f,0.5f,0.4f,0}; glLightfv    ( GL_LIGHT0, GL_DIFFUSE,  l0_diff );
	float l0_spec [] = {0.0f,0.0f,0.0f,0}; glLightfv    ( GL_LIGHT0, GL_SPECULAR, l0_spec );
	glEnable     ( GL_LIGHT0           );

    float l1_pos  [] = {+1.0f,1.0f,1.0f,0}; glLightfv    ( GL_LIGHT1, GL_POSITION, l1_pos  );
	float l1_amb  [] = {0.0f,0.0f,0.0f,0}; glLightfv    ( GL_LIGHT1, GL_AMBIENT,  l1_amb  );
	float l1_diff [] = {0.2f,0.0f,0.0f,0}; glLightfv    ( GL_LIGHT1, GL_DIFFUSE,  l1_diff );
	float l1_spec [] = {0.0f,0.0f,0.0f,0}; glLightfv    ( GL_LIGHT1, GL_SPECULAR, l1_spec );
	glEnable     ( GL_LIGHT1           );

    float l2_pos  [] = {-1.0f,0.0f,1.0f,0}; glLightfv   ( GL_LIGHT2, GL_POSITION, l2_pos  );
	float l2_amb  [] = {0.0f,0.0f,0.0f,0}; glLightfv    ( GL_LIGHT2, GL_AMBIENT,  l2_amb  );
	float l2_diff [] = {0.0f,0.0f,0.2f,0}; glLightfv  ( GL_LIGHT2, GL_DIFFUSE,  l2_diff );
	float l2_spec [] = {0.0f,0.0f,0.0f,0}; glLightfv    ( GL_LIGHT2, GL_SPECULAR, l2_spec );
	glEnable     ( GL_LIGHT2           );

    glEnable     ( GL_LIGHTING         );

    Mat3d rotMat;    rotMat.set( {1.0d,0.0d,0.0d}, {0.0d,1.0d,0.0d}, {0.0d,0.0d,1.0d} );
    BlockWorld::setupBlockWorld( {-127.0d,-127.0d,-127.0d}, {1.0d,1.0d,1.0d}, rotMat           );

    cursorShape = glGenLists(1);
    glNewList( cursorShape , GL_COMPILE );
        glPushMatrix();
        glScalef( 0.5f, 0.5f, 0.5f );
        glColor3f(0.01f,0.01f,0.01f); Draw3D::drawLines    ( Solids::Cube_nedges, Solids::Cube_edges, Solids::Cube_verts );
        Draw3D::drawAxis( 0.5f );
        glPopMatrix();
    glEndList();

/*
    shapes[nShapes] = glGenLists(1);
    glNewList( shapes[nShapes] , GL_COMPILE );
    glDisable ( GL_LIGHTING ); Draw3D::drawAxis( 0.5f );
    glEndList();
    nShapes++;
*/

    shapes[nShapes] = glGenLists(1);
    glNewList( shapes[nShapes] , GL_COMPILE );
        //glPushMatrix();
        //glScalef( 0.5f, 0.5f, 0.5f );
        glEnable ( GL_LIGHTING );
        glColor3f( 0.9f, 0.9f, 0.9f );
        //Draw3D::drawConeFan        ( 16, 0.5f, {0.0f,0.0f,-0.5f}, {0.0f,0.0f,0.5f} );
        Draw3D::drawCone        ( 4, M_PI/4.0f, 9*M_PI/4.0f, (float)sqrt(0.5f), 0.01f, {0.0f,0.0f,0.5f}, {0.0f,0.0f,0.0f}, false );
        glBegin( GL_QUADS );
        glNormal3f( 0.0f, 0.0f, 1.0f );
        glVertex3f(  0.5f,  0.5f, 0.5f);
        glVertex3f(  0.5f, -0.5f, 0.5f);
        glVertex3f( -0.5f, -0.5f, 0.5f);
        glVertex3f( -0.5f,  0.5f, 0.5f);
        glEnd();
        //Draw3D::drawCylinderStrip  ( 16, 0.5f, 0.01f, {0.0f,0.0f,-0.5f}, {0.0f,0.0f,0.5f} );
        //glPopMatrix();
    glEndList();
    nShapes++;


    shapes[nShapes] = glGenLists(1);
    glNewList( shapes[nShapes] , GL_COMPILE );
        //glPushMatrix();
        //glScalef( 0.5f, 0.5f, 0.5f );
        glEnable ( GL_LIGHTING );
        glColor3f( 0.9f, 0.9f, 0.9f );
        //Draw3D::drawConeFan        ( 16, 0.5f, {0.0f,0.0f,-0.5f}, {0.0f,0.0f,0.5f} );
        Draw3D::drawCone        ( 16, 0.0f, M_PI/2.0f, 0.5f, 0.1f, {0.0f,0.0f,-0.5f}, {0.0f,0.0f,0.5f}, false );
        //Draw3D::drawCylinderStrip  ( 16, 0.5f, 0.01f, {0.0f,0.0f,-0.5f}, {0.0f,0.0f,0.5f} );
        //glPopMatrix();
    glEndList();
    nShapes++;

    shapes[nShapes] = glGenLists(1);
    glNewList( shapes[nShapes] , GL_COMPILE );
        //glPushMatrix();
        //glScalef( 0.5f, 0.5f, 0.5f );
        glEnable ( GL_LIGHTING );
        glColor3f( 0.9f, 0.9f, 0.9f );
        //Draw3D::drawConeFan        ( 16, 0.5f, {0.0f,0.0f,-0.5f}, {0.0f,0.0f,0.5f} );
        Draw3D::drawCone        ( 5, 0.0f, 2*M_PI, 0.5f, 0.1f, {0.0f,0.0f,-0.5f}, {0.0f,0.0f,0.5f}, false );
        //Draw3D::drawCylinderStrip  ( 16, 0.5f, 0.01f, {0.0f,0.0f,-0.5f}, {0.0f,0.0f,0.5f} );
        //glPopMatrix();
    glEndList();
    nShapes++;

    shapes[nShapes] = glGenLists(1);
    glNewList( shapes[nShapes] , GL_COMPILE );
        glPushMatrix();
        //glScalef( 0.5f, 0.5f, 0.5f );
        glTranslatef( -0.5f, 0.5f, 0.0f );
        glEnable ( GL_LIGHTING );
        glColor3f( 0.9f, 0.9f, 0.9f );
        //Draw3D::drawConeFan        ( 16, 0.5f, {0.0f,0.0f,-0.5f}, {0.0f,0.0f,0.5f} );
        Draw3D::drawCone        ( 3, 0.0f, M_PI/2, 1.0f, 1.0f, {0.0f,0.0f,-0.5f}, {0.0f,0.0f,0.5f}, false );
        //Draw3D::drawCylinderStrip  ( 16, 0.5f, 0.01f, {0.0f,0.0f,-0.5f}, {0.0f,0.0f,0.5f} );
        glPopMatrix();
    glEndList();
    nShapes++;

    shapes[nShapes] = glGenLists(1);
    glNewList( shapes[nShapes] , GL_COMPILE );
        //glPushMatrix();
        //glScalef( 0.5f, 0.5f, 0.5f );
        glEnable ( GL_LIGHTING );
        glColor3f( 0.9f, 0.9f, 0.9f );
        Draw3D::drawCylinderStrip  ( 16, 0.5f, 0.5f, {0.0f,0.0f,-0.5f}, {0.0f,0.0f,0.5f} );
        //glPopMatrix();
    glEndList();
    nShapes++;

    shapes[nShapes] = glGenLists(1);
    glNewList( shapes[nShapes] , GL_COMPILE );
        glPushMatrix();
        glScalef( 0.5f, 0.5f, 0.5f );
        glEnable ( GL_LIGHTING );
        glColor3f( 0.9f, 0.9f, 0.9f );
        Draw3D::drawPolygons( Solids::Cube_nfaces,        Solids::Cube_ngons,        Solids::Cube_faces,        Solids::Cube_verts        );
        glPopMatrix();
    glEndList();
    nShapes++;

    shapes[nShapes] = glGenLists(1);
    glNewList( shapes[nShapes] , GL_COMPILE );
        glPushMatrix();
        glScalef( 0.5f, 0.5f, 0.5f );
        glEnable ( GL_LIGHTING );
        glColor3f( 0.9f, 0.9f, 0.9f );
        Draw3D::drawPolygons( Solids::Tetrahedron_nfaces, Solids::Tetrahedron_ngons, Solids::Tetrahedron_faces, Solids::Tetrahedron_verts );
        glPopMatrix();
    glEndList();
    nShapes++;

    shapes[nShapes] = glGenLists(1);
    glNewList( shapes[nShapes] , GL_COMPILE );
        glPushMatrix();
        glScalef( 0.5f, 0.5f, 0.5f );
        glEnable ( GL_LIGHTING );
        glColor3f( 0.9f, 0.9f, 0.9f );
        Draw3D::drawPolygons( Solids::Octahedron_nfaces,  Solids::Octahedron_ngons,  Solids::Octahedron_faces,  Solids::Octahedron_verts  );
        glPopMatrix();
    glEndList();
    nShapes++;

    shapes[nShapes] = glGenLists(1);
    glNewList( shapes[nShapes] , GL_COMPILE );
        glPushMatrix();
        glScalef( 0.5f, 0.5f, 0.5f );
        glEnable ( GL_LIGHTING );
        glColor3f( 0.9f, 0.9f, 0.9f );
        Draw3D::drawPolygons( Solids::RhombicDodecahedron_nfaces,        Solids::RhombicDodecahedron_ngons,        Solids::RhombicDodecahedron_faces,        Solids::RhombicDodecahedron_verts        );
        glPopMatrix();
    glEndList();
    nShapes++;

    shapes[nShapes] = glGenLists(1);
    glNewList( shapes[nShapes] , GL_COMPILE );
        glPushMatrix();
        glScalef( 0.5f, 0.5f, 0.5f );
        glEnable ( GL_LIGHTING );
        glColor3f( 0.9f, 0.9f, 0.9f );
        Draw3D::drawPolygons( Solids::Icosahedron_nfaces,        Solids::Icosahedron_ngons,        Solids::Icosahedron_faces,        Solids::Icosahedron_verts        );
        glPopMatrix();
    glEndList();
    nShapes++;

/*
    blocks[nBlocks].shape       = shapes[iShape];
    blocks[nBlocks].orientation = 0;
    blocks[nBlocks].ix          = 2;
    blocks[nBlocks].iy          = 3;
    blocks[nBlocks].iz          = 4;
    nBlocks++;
*/

};


void TestAppBlockBuilder::draw   (){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glEnable(GL_DEPTH_TEST);

    BlockWorld::drawBlock( shapes[iShape], orientation, ix, iy, iz );

    for( int i=0; i<nBlocks; i++ ){
        if( blocks[i].shape != 0 ){
            //printf( " drawing shape %i orientation %i pos (%i,%i,%i) \n", blocks[i].shape, blocks[i].orientation, blocks[i].ix, blocks[i].iy, blocks[i].iz );
            BlockWorld::drawBlock( blocks[i].shape, blocks[i].orientation, blocks[i].ix, blocks[i].iy, blocks[i].iz );
        }
    }

	glDisable ( GL_LIGHTING );
	Draw3D::drawAxis ( 3.0f );

	glDisable(GL_DEPTH_TEST);
    BlockWorld::drawBlock( cursorShape   , orientation, ix, iy, iz );

};


int TestAppBlockBuilder::findBlock ( Uint8 ix, Uint8 iy, Uint8 iz ){
    for( int i=0; i<nBlocks; i++ ){
        if( ( blocks[i].ix == ix ) && ( blocks[i].iy == iy ) && ( blocks[i].iz == iz ) ){
            return i;
        };
    }
    return -1;
};

int TestAppBlockBuilder::changeBlock( Uint8 ix, Uint8 iy, Uint8 iz, Uint8 orientation, Uint16 shape ){
    int i = findBlock ( ix, iy, iz );
    if( i < 0 ){
        if( nBlocks >= nMaxBlocks ) return -1;
        i=nBlocks;
        nBlocks++;
        blocks[i].ix = ix;
        blocks[i].iy = iy;
        blocks[i].iz = iz;
    };
    blocks[i].orientation = orientation;
    blocks[i].shape = shape;
    return i;
};

int TestAppBlockBuilder::defragmentBlocks ( ){
    int nFound = 0;
    int nEmpty = 0;
    int i      = nBlocks;
    int j      = 0;
    while( i>j ){
        if( blocks[i].shape == 0 ){ nEmpty++; i--; continue; }
        int k;
        for( k = j; k<i; k++ ){
            if( blocks[k].shape == 0 ){
                blocks[k] = blocks[i];
                i--;
                nEmpty++;
                break;
            }
            nFound++;
        }
        nFound++;
        j = k+1;
    }
    nBlocks = nFound;
    return nEmpty;
}

void TestAppBlockBuilder::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    int rot;
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_a:  ix ++; if( ix >= nMaxShapes ) ix = 255; break;
                case SDLK_d:  ix --; if( ix <  0          ) ix = 0;   break;
                case SDLK_w:  iy ++; if( iy >= nMaxShapes ) iy = 255; break;
                case SDLK_s:  iy --; if( iy <  0          ) iy = 0;   break;
                case SDLK_q:  iz ++; if( iz >= nMaxShapes ) iz = 255; break;
                case SDLK_e:  iz --; if( iz <  0          ) iz = 0;   break;
                case SDLK_x:  orientation ^= 0b10000000; break;
                case SDLK_y:  orientation ^= 0b01000000; break;
                case SDLK_z:  orientation ^= 0b00100000; break;
                case SDLK_r:  rot = orientation & 0b00011111; rot++; if(rot>= 6 ) rot=0; orientation = rot | ( orientation & 0b11100000 );  break;
                case SDLK_f:  rot = orientation & 0b00011111; rot--; if(rot<  0 ) rot=5; orientation = rot | ( orientation & 0b11100000 );  break;
                case SDLK_LEFTBRACKET:  iShape ++; if( iShape>=nShapes ) iShape = 0;            break;
                case SDLK_RIGHTBRACKET: iShape --; if( iShape< 0       ) iShape = nShapes-1;    break;
                case SDLK_RETURN:         changeBlock( ix, iy, iz, orientation, shapes[iShape] ); break;
                case SDLK_BACKSPACE: changeBlock( ix, iy, iz, 0, 0 ); break;
                case SDLK_u:  qCamera.setOne();  break;
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
    AppSDL2OGL::eventHandling( event );
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
















