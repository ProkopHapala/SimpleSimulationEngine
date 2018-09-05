
//#define SPEED_TEST

//#include <ctime>
//#include <iostream>
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
//#include "IO_utils.h"
//#include "DynamicControl.h"

#include "SoftPolyLine2D.h"
//#include "SoftPolyLine3D.h"
#include "SwordArm.h"
#include "CombatMoves.h"

#include "SDL_utils.h"
#include "Draw.h"
#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"

/*
NOTES:
 - Colorize by distance - distance perception is missing in most of first-person figting games
 - camera is automatically rotated toward enemy

*/



// ====================================
//      SwordPlayGUI
// ====================================

class SwordPlayGUI : public AppSDL2OGL_3D { public:

    const Uint8 *scanKeys;
    Uint32 mouseButtons;

    int      fontTex;

    std::vector<SwordArm> warriors;

    SwordArm* me       = 0;
    SwordArm* opponent = 0;

    // ==== function declarations

    virtual void camera     ();
    //virtual void cameraHUD();
    virtual void draw   ();
    virtual void drawHUD();
    virtual void eventHandling   ( const SDL_Event& event  );
    virtual void keyStateHandling( const Uint8 *keys );
    virtual void mouseHandling   ( );

    SwordPlayGUI( int& id, int WIDTH_, int HEIGHT_ );

};

void SwordPlayGUI::camera (){
    cam.rot.c =  (Vec3f)( opponent->head.p - me->head.p );
    cam.rot.c.normalize();
    cam.rot.b = Vec3fY;
    //cam.rot.c = (opponent->head.pos.p-me->head.pos.p).normalize();
    //cam.rot.fromDirUp( , Vec3dY );
    cam.rot.a.set_cross( cam.rot.b, cam.rot.c );
    cam.rot.a.normalize();
    Cam::perspective(cam);
    //Cam::ortho(cam, true);
}

void SwordPlayGUI::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);

}

void SwordPlayGUI::drawHUD(){

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);

}

SwordPlayGUI:: SwordPlayGUI( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {
    upTime=0;
    SDL_SetRelativeMouseMode( SDL_TRUE );
    VIEW_DEPTH = 10000;
    zoom = 30.0;

    warriors.resize(2);
    me       = &warriors[0];
    opponent = &warriors[1];

};

void SwordPlayGUI:: eventHandling   ( const SDL_Event& event  ){
    switch( event.type ){
        case SDL_KEYDOWN :
        switch( event.key.keysym.sym ){
            case SDLK_ESCAPE   : SDL_Quit(); exit(1); break;
            case SDLK_KP_PLUS  : zoom/=VIEW_ZOOM_STEP; printf("zoom: %f \n", zoom); break;
            case SDLK_KP_MINUS : zoom*=VIEW_ZOOM_STEP; printf("zoom: %f \n", zoom); break;
            //case SDLK_SPACE    : STOP = !STOP; printf( STOP ? " STOPED\n" : " UNSTOPED\n"); break;
            //case SDLK_SPACE    : SimOn = !SimOn; printf( SimOn ? " STOPED\n" : " UNSTOPED\n"); break;
        }; break;
        case SDL_QUIT: SDL_Quit(); exit(1); break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_RIGHT:
                    break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_RIGHT:
                    break;
            }
            break;
    }
};

void SwordPlayGUI:: keyStateHandling( const Uint8 *keys ){
    scanKeys = keys;

};

void SwordPlayGUI:: mouseHandling   ( ){
	// mouse Camera
	int mx,my;
	mouseButtons = SDL_GetMouseState(&mx,&my);
	int dmx = mx - WIDTH/2; 	int dmy = my - HEIGHT/2 ;
	mouseX = dmx;
	mouseY = -dmy;
	SDL_GetRelativeMouseState(&dmx,&dmy);
	qCamera.pitch( 0.005* dmy );
	qCamera.yaw  ( 0.005* dmx );
};

// =========== MAIN

SwordPlayGUI * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	SDL_ShowCursor(SDL_DISABLE);
	SDL_DisplayMode dm;
    SDL_GetDesktopDisplayMode(0, &dm);
	int junk;
	thisApp = new SwordPlayGUI( junk , dm.w-150, dm.h-100 );
	SDL_SetWindowPosition(thisApp->window, 100, 0 );
	thisApp->loop( 1000000 );
	return 0;
}


