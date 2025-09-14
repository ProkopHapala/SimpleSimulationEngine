
/*
TO DO :
	*	each of a few runs ends up as segmentation fault:
		./test: line 6: 11233 Segmentation fault      (core dumped) ./program.x
		try to use valgrind ?
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

// ========= include from common
#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"

#include "GridMap2D.h"

#include "Draw2D.h"
//#include "drawUtils.h"

#include "Body.h"
#include "Body2D.h"
#include "AeroSurf2D.h"

#include "geom2D.h"

// ========= include from local app

#include "GameWorld.h"
#include "Projectile.h"
#include "Gun.h"
#include "Yacht2D.h"
#include "Frigate2D.h"

// ===============================
// ===== GLOBAL VARIABLES
// ===============================

bool  isLib             = false;
bool  loopEnd           = false;
int   frameCount		=	0;

SDL_Event		 event;

GameWorld world;

//std::vector<Projectile*> projectiles( 100 );

#include "GameScreen.h"
#include <AppSDL2OGL.h>

GameScreen* thisScreen;
Frigate2D*  thisShip;

// ===============================
// ===== FUNCTIONS
// ===============================

void quit(){
	SDL_Quit();
	exit(1);
};

void escape(){
	if(isLib){
		printf( " ending loop \n");
		loopEnd = true;
	}else{
		printf( " exiting \n");
		quit();
	}
};


void update(){
   // world->update();   // this is now inside GameScreen update
};

void setup(){
    printf( " main.setup \n" );
    world.init();
    thisShip = world.ships.front();
    thisScreen->world = &world;
    //thisScreen->zoom = 100;
    printf( " world.init DONE \n" );
}

void inputHanding(){
	while(SDL_PollEvent(&event)){ // be carefull to include all event handling instide the while loop !!!
		if( event.type == SDL_KEYDOWN ){
			switch( event.key.keysym.sym ){
				case SDLK_ESCAPE:   escape(); break;
				case SDLK_SPACE:    thisScreen->STOP = !thisScreen->STOP; printf( thisScreen->STOP ? " STOPED\n" : " UNSTOPED\n"); break;
				case SDLK_KP_MINUS: thisScreen->zoom*=VIEW_ZOOM_STEP; break;
				case SDLK_KP_PLUS:  thisScreen->zoom/=VIEW_ZOOM_STEP; break;
				case SDLK_KP_1:     thisShip->fire_left ( &world.projectiles ); break;
				case SDLK_KP_2:     thisShip->fire_right( &world.projectiles ); break;
			}

/*
			if( event.key.keysym.sym == SDLK_ESCAPE   ) { escape(); }
			if( event.key.keysym.sym == SDLK_SPACE    ) { STOP = !STOP; printf( STOP ? " STOPED\n" : " UNSTOPED\n"); }
			if( event.key.keysym.sym == SDLK_KP_MINUS ) { thisScreen->zoom*=VIEW_ZOOM_STEP; }
			if( event.key.keysym.sym == SDLK_KP_PLUS  ) { thisScreen->zoom/=VIEW_ZOOM_STEP; }
*/
		}
		if( event.type == SDL_QUIT){ quit();  };

	} // while(SDL_PollEvent(&event))

	const Uint8 *keys = SDL_GetKeyboardState(NULL);
	if( keys[ SDL_SCANCODE_LEFT  ] ){ thisShip->rudder.setAngle( thisShip->rudder.phi + 0.01 );  }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ thisShip->rudder.setAngle( thisShip->rudder.phi - 0.01 );  }

	if( keys[ SDL_SCANCODE_UP  ]  ){ thisShip->mast.setAngle( thisShip->mast.phi + 0.01 );  }
	if( keys[ SDL_SCANCODE_DOWN ] ){ thisShip->mast.setAngle( thisShip->mast.phi - 0.01 );  }

	SDL_GetMouseState( &thisScreen->mouseX, &thisScreen->mouseY );
	//printf( "frame %i mouseX moyseY  %i %i   \n", frameCount, mouseX, mouseY );
}


void loop(int n ){
	loopEnd = false;
	for( int iframe=0; iframe<n; iframe++ ){
		inputHanding();
		if(!thisScreen->STOP){
			update();
			thisScreen->update();
			thisScreen->thisShip = thisShip; // FIXME
		}
		//printf(" %i \n", iframe );
		SDL_Delay( 10 );
		frameCount++;
		if(loopEnd) break;
	}
}

// FUNCTION ======  main
int main(int argc, char *argv[]){

	// // creating windows
	// SDL_Init(SDL_INIT_VIDEO);
	// SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	// int sid;
	// //thisScreen  = new Screen2D( sid, 800,600);
	// thisScreen  = new GameScreen( sid, 800,600 );

    setbuf(stdout, NULL);                 // disable stdout buffering
    //setvbuf(stdout, NULL, _IONBF, 0);  //Or use the more flexible setvbuf:
    // example: use like : ./spaceCraftEditor -s data/ship_ICF_interceptor_1.lua
    printf( "argc %i \n", argc );
    SDL_DisplayMode dm = initSDLOGL( 8 );
    int junk;
    //thisApp = new SailWar_single( junk, dm.w-150, dm.h-100 );   
    thisScreen  = new GameScreen( junk, dm.w-150, dm.h-100 );

	setup();

	//loop( 1 );
	loop( 1000000 );

	return 0;
}


// ==========================================================
// ===== Export this functions to Dynamic library for Python
// ==========================================================

extern "C"{ }


