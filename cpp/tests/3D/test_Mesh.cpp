
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
#include "Mesh.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"


// ============= Application

class TestAppMesh : public AppSDL2OGL_3D {
	public:

    Mesh mesh;

    bool dragging;
    Vec2f mouse0;
    int ipicked;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppMesh( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppMesh::TestAppMesh( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    mesh.fromFileOBJ( "common_resources/turret.obj" );
    printf("initialization DONE !");

    mesh.rendered_shape = glGenLists(1);
    glNewList( mesh.rendered_shape , GL_COMPILE );
        glEnable( GL_LIGHTING );
        glColor3f( 0.8f, 0.8f, 0.8f );
        Draw3D::drawMesh( mesh );
    glEndList();

}

void TestAppMesh::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    glEnable( GL_LIGHTING );
    glEnable(GL_DEPTH_TEST);
    glCallList(mesh.rendered_shape);

    glDisable( GL_LIGHTING );
    Draw3D::drawAxis(1.0);

    glDisable( GL_DEPTH_TEST );
    /*
    // this should work for perspective
    Vec3d hRay; //hRay.set_lincomb(  camMat.a, camMat.b, camMat.b, );
    camMat.dot_to( {mouseX-WIDTH*0.5f, mouseY-HEIGHT*0.5f, zoom }, hRay);
    Draw3D::drawLine      ( camPos, camPos+hRay );
    Draw3D::drawPointCross( camPos+hRay, 0.2 );
    int ip = mesh.pickVertex( camPos, camMat.c );
    */


    glColor3f(0,0,0); Draw3D::drawPointCross( mesh.points[ipicked], 0.2 );

};


void TestAppMesh::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    Vec3d ray0;
                    ray0.set_lincomb( 1, mouse_begin_x, mouse_begin_y, camPos, camMat.a, camMat.b );
                    //glColor3f(1,1,1); Draw3D::drawPointCross( ray0, 0.2 );
                    ipicked= mesh.pickVertex( ray0, camMat.c );
                    mouse0.set(mouse_begin_x, mouse_begin_y);
                    dragging=true;
                    break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    if(dragging){
                        Vec3d d;
                        d.set_lincomb( 1, mouse_begin_x, mouse_begin_y, camPos, camMat.a, camMat.b );
                        d.sub(mesh.points[ipicked]);
                        double c = d.dot(camMat.c);  d.add_mul(camMat.c, -c);
                        mesh.points[ipicked].add(d);

                        glDeleteLists(mesh.rendered_shape, 1);
                        mesh.rendered_shape = glGenLists(1);
                        glNewList( mesh.rendered_shape , GL_COMPILE );
                            glEnable( GL_LIGHTING );
                            glColor3f( 0.8f, 0.8f, 0.8f );
                            Draw3D::drawMesh( mesh );
                        glEndList();
                    }
                    dragging=false;
                    break;
            }
    };
    AppSDL2OGL::eventHandling( event );
}

void TestAppMesh::drawHUD(){
    glDisable ( GL_LIGHTING );

}

// ===================== MAIN

TestAppMesh * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppMesh( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















