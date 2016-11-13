

//#define SPEED_TEST

#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "Draw3D.h"
//#include "testUtils.h"

#include "Body.h"

#include "AeroSurf.h"
#include "AeroCraft.h"
#include "FieldPatch.h"
#include "AeroCraftWorld.h"
#include "AeroCraftGUI.h"

// ===============================
// ===== GLOBAL CONSTAMNTS
// ===============================

/*
const int    WIDTH        = 800; const int  WIDTH_HALF = WIDTH/2;
const int    HEIGHT       = 600; const int  HEIGHT_HALF = HEIGHT/2;
const int    VIEW_DEPTH   = 10000;
const float  ASPECT_RATIO = WIDTH/(float)HEIGHT ;

const float VIEW_ZOOM_STEP    = 1.2f;
const float VIEW_ZOOM_DEFAULT = 0.5f;

const int   PHYS_STEPS_PER_FRAME = 1;
const float PHYS_TIME_PER_FRAME  = 0.01;
const float PHYS_dt              = PHYS_TIME_PER_FRAME/((float)PHYS_STEPS_PER_FRAME);

*/

bool loopEnd = false;
bool STOP    = false;
//#define COLDEPTH 16

// ===============================
// ===== GLOBAL VARIABLES
// ===============================

//SDL_Surface *screen;
SDL_Event event;

int frameCount=0;

AeroCraftWorld world;
AeroCraftGUI* thisScreen;

// ====================================
// ===== FUNCTION FORWARD DECLARATIONS
// ====================================

void quit(){SDL_Quit(); exit(1);};
void setup();
void inputHanding();

// ===============================
// ===== FUNCTION IMPLEMENTATION
// ===============================


// FUNCTION ======	camera


double tickSum=0;
int    stepSum=0;

// FUNCTION ======	inputHanding
void inputHanding(){
	while(SDL_PollEvent(&event)){
		if( event.type == SDL_KEYDOWN ){
			switch( event.key.keysym.sym ){
				case SDLK_ESCAPE   : quit(); break;
				case SDLK_KP_PLUS  : thisScreen->zoom/=VIEW_ZOOM_STEP; printf("zoom: %f \n", thisScreen->zoom); break;
				case SDLK_KP_MINUS : thisScreen->zoom*=VIEW_ZOOM_STEP; printf("zoom: %f \n", thisScreen->zoom); break;
				case SDLK_SPACE    : STOP = !STOP; printf( STOP ? " STOPED\n" : " UNSTOPED\n"); break;

				case SDLK_p    : thisScreen->first_person = !thisScreen->first_person; break;
			}
		}
		if( event.type == SDL_QUIT){ quit();  };
	}


	// see http://www.libsdl.org/release/SDL-1.2.15/docs/html/sdlkey.html
	// http://sdl.beuc.net/sdl.wiki/SDL_GetKeyState

	//Uint8 *keystate = SDL_GetKeyState(NULL);
	const Uint8 *keys = SDL_GetKeyboardState(NULL);

	if ( keys[ SDL_SCANCODE_DOWN  ] ) { thisScreen->qCamera.pitch( -0.005 );  }
	if ( keys[ SDL_SCANCODE_UP    ] ) { thisScreen->qCamera.pitch(  0.005  );  }
	if ( keys[ SDL_SCANCODE_RIGHT ] ) { thisScreen->qCamera.yaw  (  0.005  );  }
	if ( keys[ SDL_SCANCODE_LEFT  ] ) { thisScreen->qCamera.yaw  ( -0.005 ); }

	float dpitch = 0.01;
	float droll  = 0.01;
	float dyaw   = 0.01;

	if      ( keys[ SDL_SCANCODE_A ] ){ world.myCraft->panels[0].lrot.rotate( -droll, { 1,0,0 } );  world.myCraft->panels[1].lrot.rotate( +0.01, { 1,0,0 } );    }
	else if ( keys[ SDL_SCANCODE_D ] ){ world.myCraft->panels[0].lrot.rotate( +droll, { 1,0,0 } );  world.myCraft->panels[1].lrot.rotate( -0.01, { 1,0,0 } );    }

    if      ( keys[ SDL_SCANCODE_W ] ){ world.myCraft->panels[2].lrot.rotate( -dpitch, { 1,0,0 } );  }
	else if ( keys[ SDL_SCANCODE_S ] ){ world.myCraft->panels[2].lrot.rotate( +dpitch, { 1,0,0 } );  }

    if      ( keys[ SDL_SCANCODE_Q ] ){ world.myCraft->panels[3].lrot.rotate( -dyaw, { 0,1,0 } );  }
	else if ( keys[ SDL_SCANCODE_E ] ){ world.myCraft->panels[3].lrot.rotate( +dyaw, { 0,1,0 } );  }


	//if ( keystate[ SDL_SCANCODE_DOWN  ] ) { qmouse.pitch2( -0.005 ); }
	//if ( keystate[ SDL_SCANCODE_UP    ] ) { qmouse.pitch2( 0.005 ); }
	//if ( keystate[ SDL_SCANCODE_RIGHT ] ) { qmouse.yaw2  ( 0.005 ); }
	//if ( keystate[ SDL_SCANCODE_LEFT  ] ) { qmouse.yaw2  ( -0.005 ); }

	// mouse Camera
	int mx,my;
	SDL_GetMouseState(&mx,&my);
	int dmx = mx - thisScreen->WIDTH/2; 	int dmy = my - thisScreen->HEIGHT/2 ;
	//printf( " mx: %i  my: %i dmx: %i dmy: %i ",mx, my, dmx, dmy );
	//qmouse.pitch( 0.001* dmy );
	//qmouse.yaw  ( 0.001* dmx );

	thisScreen->qCamera.pitch( 0.001* dmy );
	thisScreen->qCamera.yaw  ( 0.001* dmx );
	//SDL_WarpMouse( thisScreen->WIDTH/2, thisScreen->HEIGHT/2 );
	SDL_WarpMouseInWindow( thisScreen->window, thisScreen->WIDTH/2, thisScreen->HEIGHT/2 );

	Mat3d matCam;
	thisScreen->qCamera.toMatrix( matCam );
	world.steerToDir( matCam.c );


}

// FUNCTION ======	setup
void setup(){
	srand(1234);
    world.init();
    thisScreen->world = &world;
    thisScreen->qCamera.setOne();

    thisScreen->VIEW_DEPTH = 10000.0f;
    thisScreen->first_person = false;
    printf( " world.init(); DONE! \n" );
}

/*
void loop( int n ){
	loopEnd = false;
	for( int iframe=0; iframe<n; iframe++ ){
		inputHanding();
		if(!STOP){
			thisScreen->update();
		}
		//printf(" %i \n", iframe );
		if( thisScreen->delay>0 ) SDL_Delay( thisScreen->delay );
		frameCount++;
		if(loopEnd) break;
	}
}
*/

void loop(int n ){
	loopEnd = false;
	for( int iframe=0; iframe<n; iframe++ ){
        //printf( " inputHanding(); \n" );
		inputHanding();
		if(!STOP){
			//update();
			//printf( " thisScreen->update(); \n" );
			thisScreen->update();
			//thisScreen->thisShip = thisShip; // FIXME
		}
		//printf(" %i \n", iframe );
		SDL_Delay( 10 );
		//SDL_Delay(  int(PHYS_TIME_PER_FRAME*1000) );
		frameCount++;
		if(loopEnd) break;
	}
}

// FUNCTION ======  main
int main(int argc, char *argv[]){

	// creating windows
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int sid;
	//thisScreen  = new Screen2D( sid, 800,600);
	thisScreen  = new AeroCraftGUI( sid, 800,600 );

	setup();

	//loop( 1 );
	loop( 1000000 );

	return 0;
}


