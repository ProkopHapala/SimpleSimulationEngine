
/// @file @brief  This is a focused demo on object selection in a 3D scene. It uses the `raytrace.h` library to cast a ray from the camera through the mouse cursor to determine which object is being pointed at. The test specifically highlights `raySphere()` for picking spherical objects. When an object is "picked," it is typically highlighted. This is a fundamental mechanic for any interactive 3D application.
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

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

class TestAppMousePicking : public AppSDL2OGL_3D {
	public:
    //MultiFight3DWorld world;
    double dvel = 10.0;

    //std::vector<KinematicBody*> objects;
    int nobject = 100;
    Vec3d* objects;

    int defaultObjectShape, defaultObjectHitShape;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppMousePicking( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppMousePicking::TestAppMousePicking( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    defaultObjectShape = glGenLists(1);
    glNewList( defaultObjectShape , GL_COMPILE );
        //glPushMatrix();
        //glDisable ( GL_LIGHTING );
        //Draw3D::drawAxis ( 3.0f );
        //glColor3f( 1.0f, 0.0f, 1.0f ); Draw3D::drawLines   ( Solids::Icosahedron_nedges, Solids::Icosahedron_edges, Solids::Icosahedron_verts                             );
        glEnable( GL_LIGHTING );
        glColor3f( 0.8f, 0.8f, 0.8f ); Draw3D::drawPolygons( Solids::Icosahedron_nfaces, Solids::Icosahedron_ngons, Solids::Icosahedron_faces,  Solids::Icosahedron_verts );
        //glEnable( GL_LIGHTING ); glColor3f( 0.8f, 0.8f, 0.8f );   Draw3D::drawSphere_oct( 3, 1.0, {0.0,0.0,0.0} );
        //glPopMatrix();
    glEndList();
    //objects = new KinematicBody();

    defaultObjectHitShape = glGenLists(1);
    glNewList( defaultObjectHitShape , GL_COMPILE );
        //glPushMatrix();
        glDisable ( GL_LIGHTING );
        Draw3D::drawAxis ( 3.0f );
        //glColor3f( 1.0f, 0.0f, 1.0f ); Draw3D::drawLines   ( Solids::Icosahedron_nedges, Solids::Icosahedron_edges, Solids::Icosahedron_verts                             );
        glColor3f( 0.8f, 0.0f, 0.8f ); Draw3D::drawSphereOctLines( 16, 2.0, (Vec3f){0.0,0.0,0.0} );
        //glPopMatrix();
    glEndList();
    //objects = new KinematicBody();

    objects = new Vec3d[ nobject ];
    float Lspan = 50.0;
    for( int i=0; i<nobject; i++ ){
        //KinematicBody* obj = new KinematicBody();
        //obj->lpos.set( randf(-Lspan,Lspan), randf(-Lspan,Lspan), randf(-Lspan,Lspan) );
        //obj->lrot.setOne();
        //objects.push_back( obj );
        objects[i].set( randf(-Lspan,Lspan), randf(-Lspan,Lspan), randf(-Lspan,Lspan) );
    }

    //printf( " world.defaultObjectShape %i \n", world.defaultObjectShape );
    zoom = 16.0;

}

void TestAppMousePicking::draw(){
    printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    //glEnable( GL_LIGHTING );
    //glCallList( world.defaultObjectShape );
	//glDisable ( GL_LIGHTING );
	//Draw3D::drawAxis ( 3.0f );

    //printf( "camMat.a (%3.3f,%3.3f,%3.3f) \n", camMat.a.x, camMat.a.y, camMat.a.z );

    //for( auto o : world.objects ) {
    for( int i=0; i<nobject; i++ ) {
        //float glMat[16];
        glPushMatrix();
        //Draw3D::toGLMat( o->lpos, o->lrot, glMat );
        //glMultMatrixf( glMat );
        Vec3d * o = objects + i;
        glTranslatef( o->x, o->y, o->z );

        //double t = raySphere( {0.0d,0.0d,0.0d}, camMat.c, 2.0, o->lpos );
        double t = raySphere( (Vec3d)cam.pos, (Vec3d)cam.rot.c, 2.0, *o );
        if( ( t>0 ) && (t < 1000.0 ) ){
            //printf( " t %f  pos (%3.3f,%3.3f,%3.3f) \n", t, o->lpos.x, o->lpos.y, o->lpos.z );
            glCallList( defaultObjectHitShape );
        }

        glCallList( defaultObjectShape );
        glPopMatrix();
    }

};

/*
void TestAppMousePicking::mouseHandling( ){
    int mx,my;
    //Uint32 buttons = SDL_GetMouseState( &mx, &my );
    //SDL_WarpMouse(0, 0);
    SDL_GetRelativeMouseState( &mx, &my);
    //printf( " %i %i \n", mx,my );
    Quat4d q; q.fromTrackball( 0, 0, mx*0.001, my*0.001 );
    //qCamera.qmul( q );
    qCamera.qmul_T( q );
    camMat( );

    //defaultMouseHandling( mouseX, mouseY );
    //world.anchor.set( mouse_begin_x, mouse_begin_y );
}
*/

/*
void TestAppMousePicking::keyStateHandling( const Uint8 *keys ){
};
*/


void TestAppMousePicking::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
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


void TestAppMousePicking::drawHUD(){
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

TestAppMousePicking * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppMousePicking( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
