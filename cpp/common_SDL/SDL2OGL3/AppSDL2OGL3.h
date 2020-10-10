#ifndef  AppSDL2OGL3_h
#define  AppSDL2OGL3_h


#include <vector>
#include <GL/glew.h>
//#define GL_GLEXT_PROTOTYPES
//#include <GL/gl.h>
//#include <SDL2/SDL.h>

#include <fastmath.h>
#include <Vec2.h>
#include <Vec3.h>
#include <Mat3.h>
#include <quaternion.h>
#include <raytrace.h>
//#include <Body.h>

#include "Shader.h"
#include "GLObject.h"
#include "SceneOGL3.h"
#include "ScreenSDL2OGL3.h"
#include "AppSDL2OGL3.h"

class AppSDL2OGL3{ public:

    int frameCount = 0;
    bool STOP = false;

    float zoomStep      = 0.01;
    float mouseRotSpeed = 0.001;
    float keyRotSpeed   = 0.005;
    float keyMoveSpeed  = 1.0;
    //Quat4f qCamera;

    //SceneNode3D    * thisNode;
    //int nscreens=0;
    //ScreenSDL2OGL3 * thisScreen;

    int mouseX, mouseY;
    const Uint8 *keys;
    Uint32 mouseButtons;

    ScreenSDL2OGL3* screen = NULL;
    std::vector<ScreenSDL2OGL3*> screens;


    //static const int nBodies = 16;
    //PointBody bodies[nBodies];      // this would require to decouple "class Body" from SDL2OGL

    virtual void update();
    virtual void eventHandling(const SDL_Event& event);
    //virtual void mouseHandling();
    //void init();
    //void draw();
    void loop( int nframes );
    void quit();
    void die( char const *msg );

    void initSDL(int W, int H);
    AppSDL2OGL3 (int W, int H);
    //AppSDL2OGL3(){ printf("DEBUG AppSDL2OGL3()\n"); };

};

void AppSDL2OGL3::initSDL( int W, int H){
    //printf("APP DEBUG 1\n");
    if (SDL_Init(SDL_INIT_VIDEO) < 0) die( "Unable to initialize SDL" );
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3); // Opengl 3.2
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 2); // Opengl 3.2
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER,  1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE,   24);
    //printf("APP DEBUG 2\n");
    //https://forums.libsdl.org/viewtopic.php?p=41421
    SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 0);
    screens.push_back( new ScreenSDL2OGL3( W, H) );
    screen=screens[0];
    //printf("APP DEBUG 3\n");
    printf( "GL_VENDOR  : %s \n", glGetString(GL_VENDOR  ) );
	printf( "GL_VERSION : %s \n", glGetString(GL_VERSION ) );
	//printf("APP DEBUG 4\n");
    glewExperimental = true; // Needed for core profile
	if (glewInit() != GLEW_OK) {
		fprintf(stderr, "Failed to initialize GLEW\n");
		getchar();
		quit();
		//return -1;
	}
	//printf("APP DEBUG initSDL DONE\n");
    //window = SDL_CreateWindow("Tutorial2", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WIDTH, HEIGHT, SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN);
    //if ( !window ) die("Unable to create window");
}

AppSDL2OGL3::AppSDL2OGL3(int W, int H){
    initSDL(W, H);
    //qCamera.setOne();
}

// FUNCTION ======	inputHanding
void AppSDL2OGL3::update(){
	SDL_Event event;
	while(SDL_PollEvent(&event)){ eventHandling( event ); }
	keys = SDL_GetKeyboardState(NULL);
	mouseButtons = SDL_GetMouseState( &mouseX, &mouseY );
    if( screen  ){
        Quat4f& qCamera = screen->qCamera;
        Camera& cam = screen->cam;
        //if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.dyaw  (  keyRotSpeed ); }
        //if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.dyaw  ( -keyRotSpeed ); }
        if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.droll (  keyRotSpeed ); }
        if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.droll ( -keyRotSpeed ); }
        if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.dpitch(  keyRotSpeed ); }
        if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.dpitch( -keyRotSpeed ); }
        if( keys[ SDL_SCANCODE_W     ] ){ cam.pos.add_mul( cam.rot.b,  keyMoveSpeed ); printf(  "cam.pos %f %f %f \n", cam.pos.x, cam.pos.y, cam.pos.z ); }
        if( keys[ SDL_SCANCODE_S     ] ){ cam.pos.add_mul( cam.rot.b, -keyMoveSpeed ); }
        if( keys[ SDL_SCANCODE_E     ] ){ cam.pos.add_mul( cam.rot.c,  keyMoveSpeed ); }
        if( keys[ SDL_SCANCODE_Q     ] ){ cam.pos.add_mul( cam.rot.c, -keyMoveSpeed ); }
        if( keys[ SDL_SCANCODE_A     ] ){ cam.pos.add_mul( cam.rot.a,  keyMoveSpeed ); }
        if( keys[ SDL_SCANCODE_D     ] ){ cam.pos.add_mul( cam.rot.a, -keyMoveSpeed ); }

        if( keys[ SDL_SCANCODE_KP_PLUS  ] ){ cam.zoom*=(1.+zoomStep); }
        if( keys[ SDL_SCANCODE_KP_MINUS ] ){ cam.zoom/=(1.+zoomStep); }
        int dmx,dmy;
        Uint32 buttons = SDL_GetRelativeMouseState( &dmx, &dmy);
        bool polarCamera = true;
        if(polarCamera){
            if ( buttons & SDL_BUTTON(SDL_BUTTON_RIGHT)) {
                screen->cam.rot.rotate( dmx*mouseRotSpeed, {0.0,1.0,0.0} );
                screen->cam.rot.rotate( dmy*mouseRotSpeed,screen->cam.rot.a );
                qCamera.fromMatrix(screen->cam.rot);
            }
        }else{
            if ( buttons & SDL_BUTTON(SDL_BUTTON_RIGHT)) {
                Quat4f q; q.fromTrackball( 0, 0, -dmx*mouseRotSpeed, dmy*mouseRotSpeed );
                //qCamera.qmul_T( q );
                qCamera.qmul( q );
            }
            qCamera.toMatrix_T( screen->cam.rot );
        }
        //printf( "screen->camLookAt %i \n", screen->camLookAt );
        if(screen->camLookAt) cam.lookAt( *screen->camLookAt, screen->camDist );
    }
}

void AppSDL2OGL3::eventHandling( const SDL_Event& event ){
		if( event.type == SDL_KEYDOWN ){
			if(event.key.keysym.sym == SDLK_ESCAPE ) { quit(); }
			if(event.key.keysym.sym == SDLK_SPACE  ) { STOP=!STOP; }
		}
		if( event.type == SDL_QUIT){ quit();  };
};

void AppSDL2OGL3::quit(){
	//glDeleteVertexArrays(1, &vao);
    //if( context != NULL ) SDL_GL_DeleteContext( context );
    //if( window  != NULL ) SDL_DestroyWindow   ( window  );
    SDL_Quit();
	exit(0);
};

void AppSDL2OGL3::die( char const *msg ){
    printf("%s: %s\n", msg, SDL_GetError());
    SDL_Quit();
    exit(1);
}

void AppSDL2OGL3::loop( int nframes ){
    for ( int iframe=1; iframe<nframes; iframe++)    {
        update();
 		if(!STOP) for(auto screen:screens) screen->draw();
		frameCount++;
        SDL_Delay(10);
    }
}

#endif
