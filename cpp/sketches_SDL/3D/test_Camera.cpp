
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw.h"
#include "Draw2D.h"
#include "Draw3D.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"
//#include "Body.h"
#include "Camera.h"
#include "cameraOGL.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"
#include "SDL_utils.h"

char str[1024];


/*

TODO:
   - When switching perspective and orthographic, try to keep overall zoom

*/

class TestAppCamera : public AppSDL2OGL_3D {
	public:
    //MultiFight3DWorld world;
    double dvel = 10.0;

    //std::vector<KinematicBody*> objects;
    int nobject = 100;
    Vec3d* objects;

    int nCPs = 0;
    Vec3f* CPs=0;

    Camera cam;

    int fontTex=0;

    int defaultObjectShape, defaultObjectHitShape;

	virtual void draw   ();
	virtual void drawHUD();
	virtual void camera();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	virtual void keyStateHandling( const Uint8 *keys );

	TestAppCamera( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppCamera::TestAppCamera( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );

    defaultObjectShape = glGenLists(1);
    glNewList( defaultObjectShape , GL_COMPILE );
        glEnable( GL_LIGHTING );
        glColor3f( 0.8f, 0.8f, 0.8f ); Draw3D::drawPolygons( Solids::Icosahedron_nfaces, Solids::Icosahedron_ngons, Solids::Icosahedron_faces,  Solids::Icosahedron_verts );
    glEndList();

    defaultObjectHitShape = glGenLists(1);
    glNewList( defaultObjectHitShape , GL_COMPILE );
        glDisable ( GL_LIGHTING );
        Draw3D::drawAxis ( 3.0f );
        glColor3f( 0.8f, 0.0f, 0.8f ); Draw3D::drawSphereOctLines( 16, 2.0, (Vec3f){0.0,0.0,0.0} );
    glEndList();

    objects = new Vec3d[ nobject ];
    float Lspan = 50.0;
    for( int i=0; i<nobject; i++ ){
        objects[i].set( randf(-Lspan,Lspan), randf(-Lspan,Lspan), randf(-Lspan,Lspan) );
    }
    zoom = 16.0;

    cam.zoom   = zoom;
    cam.aspect = ASPECT_RATIO;

    nCPs = 4;
    CPs = new Vec3f[nCPs];

    CPs[0].set(0.0,0.0,0.0);
    CPs[1].set(1.0,0.0,0.0);
    CPs[2].set(0.0,2.0,0.0);
    CPs[3].set(0.0,0.0,3.0);

}

void TestAppCamera::camera(){
    // AppSDL2OGL_3D::camera();
    //Quat4f qCam = (Quat4f)qCamera;
    ((Quat4f)qCamera).toMatrix(cam.rot);
    //Cam::ortho( cam, true );
    //Cam::perspective( cam );
    if (perspective){ Cam::perspective( cam ); }
    else            { Cam::ortho( cam, true ); }
}

void TestAppCamera::draw(){
    printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	/*
    for( int i=0; i<nobject; i++ ) {
        glPushMatrix();
        Vec3d * o = objects + i;
        glTranslatef( o->x, o->y, o->z );
        double t = raySphere( (Vec3d)cam.pos, (Vec3d)cam.rot.c, 2.0, *o );
        if( ( t>0 ) && (t < 1000.0 ) ){
            glCallList( defaultObjectHitShape );
        }
        glCallList( defaultObjectShape );
        glPopMatrix();
    }
    */

    glColor3f(0.7,0.7,0.7);
    Draw3D::drawPanel( Vec3fZero, Mat3fIdentity, {20.0,10.0} );
    Draw3D::drawAxis( 100.0 );
    Draw3D::drawMatInPos( cam.rot, cam.pos );

    glColor3f(0.5,1.0,0.5);
    for( int i=0; i<nCPs; i++ ) {
        Draw3D::drawPointCross( CPs[i], 0.2 );
        //Draw::drawText( );
    }

};

void TestAppCamera::drawHUD(){
    glDisable ( GL_LIGHTING );
    glDisable ( GL_DEPTH_TEST );

    glColor3f( 0.0f, 1.0f, 0.0f );
    glBegin( GL_LINES );
    float whalf = WIDTH *0.5;
    float hhalf = HEIGHT*0.5;
    glVertex3f( whalf-10,hhalf, 0 ); glVertex3f( whalf+10,hhalf, 0 );
    glVertex3f( whalf,hhalf-10, 0 ); glVertex3f( whalf,hhalf+10, 0 );
    glEnd();

    Vec3f pScreen;
    for( int i=0; i<nCPs; i++ ) {
        glColor3f(1.0,0.5,1.0);

        Vec2f pix;

        if (perspective){ pix = cam.word2pixPersp( CPs[i], {(float)WIDTH,(float)HEIGHT} ); }
        else            { pix = cam.word2pixOrtho( CPs[i], {(float)WIDTH,(float)HEIGHT} ); }
        //cam.word2screenPersp( CPs[i], pScreen );
        //Draw2D::drawPointCross( (Vec2f){pScreen.x*WIDTH,pScreen.y*HEIGHT} + (Vec2f){WIDTH*0.5f,HEIGHT*0.5f}, 5 );
        Draw2D::drawPointCross( pix, 5 );

        //cam.word2screenOrtho( CPs[i], pScreen );
        //Draw2D::drawPointCross( (Vec2f){pScreen.x*WIDTH,pScreen.y*HEIGHT} + (Vec2f){WIDTH*0.5f,HEIGHT*0.5f}, 5 );
        //Draw2D::drawPointCross( cam.word2pixOrtho( CPs[i], {(float)WIDTH,(float)HEIGHT} ), 5 );

        glColor3f(0.0,0.0,0.0);
        sprintf( str, "%i", i );
        Draw2D::drawText( str, 0, (Vec2d)pix, 0.0, fontTex, fontSizeDef );
    }
    //glColor3f(0.0,0.0,0.0);
    //Draw2D::drawText( "abcdefghijklmnopqrstuvwxyz \n0123456789 \nABCDEFGHIJKLMNOPQRTSTUVWXYZ \nxvfgfgdfgdfgdfgdfgdfg", 0, {100.0, 300.0}, 0.0, fontTex, fontSizeDef );
    //Draw::drawText( "abcdefghijklmnopqrstuvwxyz \n0123456789 \nABCDEFGHIJKLMNOPQRTSTUVWXYZ \nxvfgfgdfgdfgdfgdfgdfg", fontTex, fontSizeDef, {10,5} );
}

void TestAppCamera::keyStateHandling( const Uint8 *keys ){
    if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.dyaw  (  keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.dyaw  ( -keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.dpitch(  keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.dpitch( -keyRotSpeed ); }

    float moveSpeed = 1.0f;
	if( keys[ SDL_SCANCODE_A ] ){ cam.pos.add_mul( cam.rot.a, -moveSpeed ); }
	if( keys[ SDL_SCANCODE_D ] ){ cam.pos.add_mul( cam.rot.a,  moveSpeed ); }
    if( keys[ SDL_SCANCODE_W ] ){ cam.pos.add_mul( cam.rot.b,  moveSpeed ); }
	if( keys[ SDL_SCANCODE_S ] ){ cam.pos.add_mul( cam.rot.b, -moveSpeed ); }
    if( keys[ SDL_SCANCODE_Q ] ){ cam.pos.add_mul( cam.rot.c, -moveSpeed ); }
	if( keys[ SDL_SCANCODE_E ] ){ cam.pos.add_mul( cam.rot.c,  moveSpeed ); }

};

void TestAppCamera::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;

                case SDLK_ESCAPE:   quit(); break;
                case SDLK_KP_MINUS: cam.zoom*=VIEW_ZOOM_STEP; break;
                case SDLK_KP_PLUS:  cam.zoom/=VIEW_ZOOM_STEP; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            break;
    };
    //AppSDL2OGL::eventHandling( event );
}



// ===================== MAIN

TestAppCamera * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppCamera( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















