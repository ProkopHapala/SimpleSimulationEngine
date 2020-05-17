
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

#include "EditorGizmo.h"
#include "GUI.h"


#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

class TestAppMousePicking : public AppSDL2OGL_3D {
	public:
    //MultiFight3DWorld world;
    double dvel = 10.0;

    //std::vector<KinematicBody*> objects;
    int nobject = 100;
    EditorGizmo gizmo;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppMousePicking( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppMousePicking::TestAppMousePicking( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {
    gizmo.cam = &cam;

    int np = 30;
    gizmo.npoint = np;
    gizmo.points = new Vec3d[np];
    double sz = 5.0;
    for(int i=0; i<gizmo.npoint; i++){
        gizmo.points[i].fromRandomBox( {-sz,-sz,-sz}, {sz,sz,sz} );
    }
}

void TestAppMousePicking::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	gizmo.draw();

	Draw3D::drawPoints( gizmo.npoint, gizmo.points ,0.1 );


	for(auto& it: gizmo.selection){
        Vec3d p = gizmo.points[it.first];
        Draw::color_of_hash( it.second + 15454 );
        Draw3D::drawPointCross( p, 0.2 );
	};

};




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
    };
    AppSDL2OGL::eventHandling( event );
    Vec2f pix = { 2*mouseX/float(HEIGHT) - ASPECT_RATIO,
                  2*mouseY/float(HEIGHT) - 1              };
    cam.persp = perspective;
    cam.zoom  = zoom;
    gizmo.onEvent( pix, event );
}


void TestAppMousePicking::drawHUD(){
    /*
    glDisable ( GL_LIGHTING );
    glColor3f( 0.0f, 1.0f, 0.0f );
    glBegin( GL_LINES );
    float whalf = WIDTH *0.5;
    float hhalf = HEIGHT*0.5;
    glVertex3f( whalf-10,hhalf, 0 ); glVertex3f( whalf+10,hhalf, 0 );
    glVertex3f( whalf,hhalf-10, 0 ); glVertex3f( whalf,hhalf+10, 0 );
    glEnd();
    */
}

// ===================== MAIN

TestAppMousePicking * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppMousePicking( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















