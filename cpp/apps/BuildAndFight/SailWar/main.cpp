
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
#include "math/fastmath.h"
#include "math/Vec2.h"
#include "math/Vec3.h"
#include "math/Mat3.h"
#include "math/quaternion.h"
#include "math/raytrace.h"
#include "SDL2OGL/drawMath.h"
#include "SDL2OGL/drawMath2D.h"
#include "SDL2OGL/drawUtils.h"
#include "dynamics/rigidBody.h"
#include "dynamics/Body2D.h"
#include "dynamics/AeroSurf2D.h"
#include "dynamics/Chain2D.h"

// ========= include from local app
#include "include/GameWorld.h"
#include "include/Projectile.h"
#include "include/Gun.h"

#include "include/Ship2D.h"
#include "include/Frigate2D.h"

// ===============================
// ===== GLOBAL CONSTAMNTS
// ===============================

const float	VIEW_ZOOM_STEP     = 1.2f;
const float	VIEW_ZOOM_DEFAULT  = 10.0f;
const float	VIEW_DEPTH_DEFAULT = 1000.0;

// ===============================
// ===== GLOBAL VARIABLES
// ===============================

bool  isLib             = false;
bool  STOP          	= false; 
bool  loopEnd           = false;
int   frameCount		=	0;

SDL_Event		 event; 
int perFrame = 10;

double dt = 0.0001;


GameWorld world;

Vec2d windSpeed, watterSpeed;

Frigate2D   ship1;
Frigate2D   ship2;

//std::vector<Projectile*> projectiles( 100 );  
std::vector<Projectile*> projectiles;  // see http://stackoverflow.com/questions/11457571/how-to-set-initial-size-of-stl-vector

const int npts = 4;
static double poss[npts*2] = { -1.0, 0.0,   0.0, -0.1,   0.0, +0.1,   +1.0, 0.0  };
static double mass[npts  ] = {  10.0, 50.0, 50.0, 10.0  };


#include "SDL2OGL/Screen2D.h"
#include "include/GameScreen.h"
GameScreen* thisScreen;

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
	
	//ship1.applySailForces(  { 0.0, 0.0 },  { 0.0, 1.0 }  );
	//ship1.move( dt );

/*
	for( int i=0; i<perFrame; i++ ){
		ship1.clean_temp( );
		ship1.applySailForces(  windSpeed,  watterSpeed );
		ship1.move( dt );
	}
*/

	if( ship1.phi != ship1.phi ) STOP = true;

	//STOP = true;
};

void setup(){
	int ifree,igl,nvert,ndiv;

	world.ground_level = 0.0d;
	world.watter_speed.set(   0.0, 0.0      );
	world.wind_speed  .set( -10.0, 0.0, 0.0 );

	windSpeed  .set( -10.0, 0.0 );
	watterSpeed.set(  0.0, 0.0 );

	int FigateShape = glGenLists(1);
	glNewList( FigateShape , GL_COMPILE );
	glBegin   (GL_TRIANGLE_FAN);	       
		glNormal3f( 0.0f, 0.0f, 1.0f ); 
		glVertex3f( +1.5,  0.0, 0 );
 		glVertex3f( +0.5,  0.2, 0 );
		glVertex3f( -1.0,  0.2, 0 );
 		glVertex3f( -1.0, -0.2, 0 );
		glVertex3f( +0.5, -0.2, 0 );
		glVertex3f( +1.5,  0.0, 0 );
	glEnd();
	glEndList();

	printf( " >>> Setup  ship1: \n" );
	ship1.loadFromFile( "data/FrigateType.txt" );
	ship1.from_mass_points( 2, mass, (Vec2d*)poss );  printf( " I invI  %f %f \n", ship1.I, ship1.invI );
	ship1.setDefaults();
	ship1.setAngle( M_PI*0.6   );
	ship1.pos.set ( {0.0, 0.0} );
	ship1.omega = 0.0;
	ship1.shape = FigateShape;
	ship1.initAllGuns( 1 );

	printf( " >>> Setup  ship2: \n" );
	ship2.loadFromFile( "data/FrigateType.txt" );
	ship2.from_mass_points( 2, mass, (Vec2d*)poss );  printf( " I invI  %f %f \n", ship1.I, ship1.invI );
	ship2.setDefaults();
	ship2.setAngle( M_PI*0.6   );
	ship2.pos.set ( {1., 1.0} );
	ship2.omega = 0.0;
	ship2.shape = FigateShape;

	printf( " >>> Setup  ship1 DONE \n" );

	projectiles.reserve(100); 

}

void inputHanding(){
	while(SDL_PollEvent(&event)){ // be carefull to include all event handling instide the while loop !!!
		if( event.type == SDL_KEYDOWN ){ 
			switch( event.key.keysym.sym ){
				case SDLK_ESCAPE:   escape(); break;
				case SDLK_SPACE:    STOP = !STOP; printf( STOP ? " STOPED\n" : " UNSTOPED\n"); break;
				case SDLK_KP_MINUS: thisScreen->zoom*=VIEW_ZOOM_STEP; break;
				case SDLK_KP_PLUS:  thisScreen->zoom/=VIEW_ZOOM_STEP; break;
				case SDLK_KP_1:     ship1.fire_left ( &projectiles ); break;
				case SDLK_KP_2:     ship1.fire_right( &projectiles ); break;
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
	if( keys[ SDL_SCANCODE_LEFT  ] ){  ship1.rudder.setAngle( ship1.rudder.phi + 0.01 );  }
	if( keys[ SDL_SCANCODE_RIGHT ] ){  ship1.rudder.setAngle( ship1.rudder.phi - 0.01 );  }

	if( keys[ SDL_SCANCODE_UP  ]  ){  ship1.mast.setAngle( ship1.mast.phi + 0.01 );  }
	if( keys[ SDL_SCANCODE_DOWN ] ){  ship1.mast.setAngle( ship1.mast.phi - 0.01 );  }

	SDL_GetMouseState( &thisScreen->mouseX, &thisScreen->mouseY );
	//printf( "frame %i mouseX moyseY  %i %i   \n", frameCount, mouseX, mouseY );
}


void loop(int n ){
	loopEnd = false;
	for( int iframe=0; iframe<n; iframe++ ){
		inputHanding(); 
		if(!STOP){
			update();
			thisScreen->update();
		}  
		//printf(" %i \n", iframe );
		SDL_Delay( 10 );
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
	thisScreen  = new GameScreen( sid, 800,600 ); 

	setup();

	printf( " setup done \n" );

	//loop( 1 );
	loop( 1000000 );

	return 0;
}


// ==========================================================
// ===== Export this functions to Dynamic library for Python
// ==========================================================

extern "C"{ }


