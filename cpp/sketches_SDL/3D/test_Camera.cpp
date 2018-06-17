
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
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

class TestAppCamera : public AppSDL2OGL_3D {
	public:
    //MultiFight3DWorld world;
    double dvel = 10.0;

    //std::vector<KinematicBody*> objects;
    int nobject = 100;
    Vec3d* objects;
    Camera cam;

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

}

void TestAppCamera::camera(){
   // AppSDL2OGL_3D::camera();
   //Quat4f qCam = (Quat4f)qCamera;
   ((Quat4f)qCamera).toMatrix(cam.rot);
   //Cam::ortho( cam, true );
   Cam::perspective( cam );
}

void TestAppCamera::draw(){
    printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

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
    Draw3D::drawPanel( Vec3fZero, Mat3fIdentity, {20.0,10.0} );
    Draw3D::drawAxis( 100.0 );
    Draw3D::drawMatInPos( cam.rot, cam.pos );

};

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
                case SDLK_p:  first_person = !first_person; break;
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


void TestAppCamera::drawHUD(){
    glDisable ( GL_LIGHTING );
    glColor3f( 0.0f, 1.0f, 0.0f );
    glBegin( GL_LINES );
    float whalf = WIDTH *0.5;
    float hhalf = HEIGHT*0.5;
    glVertex3f( whalf-10,hhalf, 0 ); glVertex3f( whalf+10,hhalf, 0 );
    glVertex3f( whalf,hhalf-10, 0 ); glVertex3f( whalf,hhalf+10, 0 );
    glEnd();
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
















