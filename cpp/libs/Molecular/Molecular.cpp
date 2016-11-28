
#include <stdlib.h>
#include <stdio.h>

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"

#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"



class TestAppMesh : public AppSDL2OGL_3D {
	public:

    Mesh mesh;

    bool dragging;
    Vec2f mouse0;
    int ipicked;

	//virtual void draw   ();
	//virtual void drawHUD();
	//virtual void mouseHandling( );
	//virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppMesh( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppMesh::TestAppMesh( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

};


TestAppMesh * thisApp;


extern "C"{

    void printHello(){
        printf("Hello!\n");
    }

    void initWindow(){
        SDL_Init(SDL_INIT_VIDEO);
        SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
        int junk;
        thisApp = new TestAppMesh( junk , 800, 600 );
        thisApp->loop( 1000000 );
    }

}
