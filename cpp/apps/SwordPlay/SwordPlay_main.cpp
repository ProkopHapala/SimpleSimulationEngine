
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

#include "MusculeSkelet.h"

#include "SDL_utils.h"
#include "Draw.h"
#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"

/*
NOTES:
 - Colorize by distance - distance perception is missing in most of first-person figting games
 - camera is automatically rotated toward enemy

*/

void drawMusculeSkelet(const MusculeSkelet& body){
    glColor3f(1,1,1);
    Draw3D::drawPoints( body.npoints, body.points, 0.2);

    glBegin(GL_LINES);
    glColor3f(0,0,0);
    for(int i=0; i<body.nbonds; i++){
        const Bond&  b  = body.bonds [i];
        const Vec3d& p1 = body.points[b.i];
        const Vec3d& p2 = body.points[b.j];
        glVertex3f(p1.x,p1.y,p1.z);   glVertex3f(p2.x,p2.y,p2.z);
    }
    glColor3f(1,0,0);

    for(int i=0; i<body.nmuscules; i++){
        Vec3d p1,p2;
        body.getMusculeAnchors(i,p1,p2);
        glVertex3f(p1.x,p1.y,p1.z);   glVertex3f(p2.x,p2.y,p2.z);
    }

    glEnd();
}

void makeWarrior_type1(MusculeSkelet& body){

    const int np=13;
    Vec3d ps[np]{
        {  2,0, 1  }, // 0 foot R
        { -2,0, 1  }, // 1 foot L
        {  0,0,-1 }, // 2 tail tip (?)

        {  2,2, 1 }, // 3 hip R
        { -2,2, 1 }, // 4 hip L
        {  0,2,-1 }, // 5 loin

        {  2,4, 1 }, // 6 shoulder R
        { -2,4, 1 }, // 7 shoulder L
        {  0,4,-1 }, // 8 neck

        {  2,4, 2.5 }, // 9  elbow R
        { -2,4, 2.5 }, // 10 elbow L
        // sword
        {  0,5, 3.5 }, // 11 hilt
        {  0,7, 5.0 }, // 12 point
    };

    const int nb=15;
    Vec2i links[nb]{
        // Long
        { 0,3 }, // 0 leg R
        { 1,4 }, // 1 leg L
        { 2,5 }, // 2 tail
        { 5,8 }, // 3 spine
        { 6,9 }, // 4 arm R
        { 7,10}, // 5 arm L
        { 9,11 }, // 6 forearm R
        { 10,11}, // 7 forearm L
        { 11,12}, // 8 sword
        // Flat bones
        // pelvis
        {3,4}, // 9 pubic (pelvis-fw)
        {3,5}, // 10 ilium right
        {4,5}, // 11 ilium Left
        // Shoulders
        {6,7}, // 12 colar bone
        {6,8}, // 13 Scapula Right
        {7,8}, // 14 Scapula Left

        /*
        // push-muscules
        // between feets
        { 0,1 }, // 0 L-R foots
        { 0,2 }, // 1 Rfoot-tail
        { 1,2 }, // 2 Lfoot-tail
        // Arm muscules
        { 7,9 }, // 3 arm R Pectoralis
        { 6,10}, // 4 arm L Pectoralis

        { 8,9 }, // 5 arm R Deltoid
        { 8,10}, // 6 arm L Deltoid

        { 6,11 }, // 7 arm R Deltoid
        { 7,11 }, // 8 arm L Deltoid

        { 9,12 },  // 9  arm R Deltoid
        { 10,12 }, // 10 arm L Deltoid
        */
    };

    const int nmus=15;
    Muscule muscules[nmus]{
        // push-muscules
        // between legs
        { {0,1}, {0.5,0.5}, },  // 0  Between Right-Left leg
        { {0,2}, {0.5,0.5}, },  // 1  Right-tail
        { {1,2}, {0.5,0.5}, },  // 2  Left-tail
        // torso
        { {10,13}, {0.5,0.5}, },  // 3 Trapz Right
        { {11,14}, {0.5,0.5}, },  // 4 Trapz Left
        // Pectoralis
        { {4,12}, {0.5,0.5}, },  // 5 Pectoralis Right
        { {5,12}, {0.5,0.5}, },  // 6 Pectoralis Left
        // Deltoid
        { {4,13}, {0.5,0.5}, },  // 7  Deltoid Right
        { {5,14}, {0.5,0.5}, },  // 8  Deltoid Left
        // elbow
        { {4,6}, {0.5,0.5}, },  // 9   elbow Right
        { {5,7}, {0.5,0.5}, },  // 10  elbow Left
        // wrist
        { {6,8}, {0.5,0.5}, },  // 11  wrist Right
        { {7,8}, {0.5,0.5}, },  // 12  wrist Left
    };

    body.allocate( np , nb, 0, nmus );
    for(int i=0;i<np;i++  ){ body.points  [i]=ps[i];       };
    for(int i=0;i<nmus;i++){ body.muscules[i]=muscules[i]; };

    for(int i=0;i<nb;i++){
        const Vec2i& l = links[i];
        Bond&  b       = body.bonds[i];
        b.id=i;
        b.i=l.i;
        b.j=l.j;
    };

}


// ====================================
//      SwordPlayGUI
// ====================================

class SwordPlayGUI : public AppSDL2OGL_3D { public:

    const Uint8 *scanKeys;
    Uint32 mouseButtons;

    int      fontTex;

    std::vector<MusculeSkelet> warriors;

    //SwordArm* me       = 0;
    //SwordArm* opponent = 0;


    //MusculeSkelet body1;
    //MusculeSkelet body2;

    MusculeSkelet* me;
    MusculeSkelet* opponent;

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
    /*
    cam.rot.c =  (Vec3f)( opponent->head.p - me->head.p );
    cam.rot.c.normalize();
    cam.rot.b = Vec3fY;
    //cam.rot.c = (opponent->head.pos.p-me->head.pos.p).normalize();
    //cam.rot.fromDirUp( , Vec3dY );
    cam.rot.a.set_cross( cam.rot.b, cam.rot.c );
    cam.rot.a.normalize();
    Cam::perspective(cam);
    //Cam::ortho(cam, true);
    */
    AppSDL2OGL_3D::camera();
}

void SwordPlayGUI::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);

    for(const MusculeSkelet& body: warriors){
        drawMusculeSkelet( body );
    }

    Draw3D::drawAxis(5);

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

    makeWarrior_type1( warriors[0] );
    makeWarrior_type1( warriors[1] );

    for(int i=0;i<opponent->npoints;i++){
        opponent->points[i].mul({-1.0,1.0,-1.0});
        opponent->points[i].add(0,0,15.0);
    }

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


