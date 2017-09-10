

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

#include "Shooter.h"

#include "AeroSurf.h"
#include "AeroCraft.h"
#include "AeroCraftWarrior.h"

#include "FieldPatch.h"
//#include "AeroCraftWorld.h"
//#include "AeroCraftGUI.h"

#include "SDL_utils.h"
#include "Draw.h"
#include "Draw3D.h"
#include "AeroDraw.h"
#include "AeroTest.h"

#include "GUI.h"
#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

// ===============================
// ===== GLOBAL CONSTAMNTS
// ===============================

AeroTester         * tester      = NULL;
AeroCraftControler * autoPilot1  = NULL;
Shooter            * world       = NULL;

//AeroTester         tester     ;
//AeroCraftControler autoPilot1 ;
//Shooter            world      ;

//SDL_Surface *screen;
SDL_Event event;

int frameCount=0;
double tickSum=0;
int    stepSum=0;

bool loopEnd = false;
bool STOP    = false;

// ====================================
//      AeroCraftGUI
// ====================================

class AeroCraftGUI : public AppSDL2OGL_3D { public:
	//AeroCraftWorld * world;
	//Shooter * world;

    int      fontTex;
    GUIPanel   panel;
    MultiPanel mpanel;
    GUITextInput txt;

    GUIAbstractPanel*  focused = NULL;

    bool mouseSteer = false;
    bool autoPilot  = false;

	int perFrame = 10;
	//double dt = 0.001;

	AeroCraftWarrior * myCraft;
	//AeroCraft * myCraft;
	AeroCraft * myCraft_bak;

	FieldPatch fieldPatch;
	int buildings_shape = -1;
	int terrain_shape   = -1;

    bool staticTest = true;
    // - put to Spline manager ? ... make indepemented AeroCraft Test ?

	//int fontTex_DEBUG;

	// ==== function declarations

	//void reallocateTrj(int n);
	void resetSteer( );
	void steerToDir( const Vec3d& dir );

	//void renderSkyBox( float x0, float y0, float z0 );
	//void drawStaticTest2D();

	virtual void camera     ();
	//virtual void cameraHUD();
	virtual void draw   ();
	virtual void drawHUD();

    virtual void eventHandling   ( const SDL_Event& event  );
	virtual void keyStateHandling( const Uint8 *keys );
    virtual void mouseHandling   ( );

	AeroCraftGUI( int& id, int WIDTH_, int HEIGHT_ );

};

void AeroCraftGUI::camera (){
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    float fov = VIEW_ZOOM_DEFAULT/zoom;
    glFrustum( -ASPECT_RATIO, ASPECT_RATIO, -1, 1, 1*fov, VIEW_DEPTH*fov );
    Mat3d camMat;
    Vec3f camPos;
    convert(myCraft->pos, camPos );
    float camDist = 10.0;
    if(first_person){
        // third person camera attached to aero-craft
        camMat.setT( myCraft->rotMat );
        //glTranslatef ( -camPos.x, -camPos.y, -camPos.z );
    }else{
        // third person camera attached to aero-craft
        qCamera.toMatrix( camMat );
        camMat.T();
	}
	float glMat[16];
	Draw3D::toGLMatCam( { 0.0f, 0.0f, 0.0f}, camMat, glMat );
	glMultMatrixf( glMat );
    glTranslatef ( -camPos.x+camMat.cx*camDist, -camPos.y+camMat.cy*camDist, -camPos.z+camMat.cz*camDist );
    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity();

}

void AeroCraftGUI::resetSteer( ){
    myCraft->panels[0].lrot = myCraft_bak->panels[0].lrot;
    myCraft->panels[1].lrot = myCraft_bak->panels[1].lrot;
    myCraft->panels[2].lrot = myCraft_bak->panels[2].lrot;
    myCraft->panels[3].lrot = myCraft_bak->panels[3].lrot;
}

void AeroCraftGUI::steerToDir( const Vec3d& dir ){
    Mat3d rotMatT;
    rotMatT.setT(myCraft->rotMat);
    Draw3D::drawMatInPos( rotMatT, myCraft->pos );
    double dyaw   = rotMatT.a.dot( dir );
    double dpitch = rotMatT.b.dot( dir );
    const double acut = 0.1;
    //double droll = (a>acut)?(a-acut):((a<-acut)?(a+acut):0.0d);
    double droll = dyaw;
    resetSteer();
    myCraft->steerTo( 0.1*droll, 0.5*dpitch+0.5*fabs(dyaw), 0.5*dyaw);
};

void AeroCraftGUI::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	if(staticTest) return;

	world->update_world(); // ALL PHYSICS COMPUTATION DONE HERE
	camera ();

	renderSkyBox(myCraft->pos.x, myCraft->pos.y-1000, myCraft->pos.z, VIEW_DEPTH*0.25 );
	glEnable(GL_DEPTH_TEST);

	glEnable    (GL_LIGHTING);
	glShadeModel(GL_FLAT);

	//world->myCraft->render();
	renderAeroCraft( *myCraft );

	//glDisable (GL_LIGHTING);
	glShadeModel(GL_SMOOTH);

	if ( buildings_shape >0 ) glCallList( buildings_shape );
	if ( terrain_shape >0)  {
		glCallList( terrain_shape   );
	} else {
	 	// terrain
		float groundsz = VIEW_DEPTH;
		glBegin(GL_QUADS);
			glColor3f( 0.3, 0.6, 0.1 );
			glNormal3f(0,1,0);
			glVertex3f( -groundsz, 0, -groundsz );
			glVertex3f( +groundsz, 0, -groundsz );
			glVertex3f( +groundsz, 0, +groundsz );
			glVertex3f( -groundsz, 0, +groundsz );
		glEnd();
	};

	//Draw3D::drawAxis( 1000 );

	/*
	glColor4f(1.0f,1.0f,1.0f,0.9f);
	char str[256];
	sprintf(str, "speed %3.3f\0",world->myCraft->vel.norm());
	Draw3D::drawText(str, world->myCraft->pos, fontTex, 0.2, 0, 0 );
	//Draw3D::drawText( "AHOJ!\0", world->myCraft->pos, fontTex, 0.2, 0, 0 );
	//Draw3D::drawText( "AHOJ!\0", {0.0,0.0,0.0}, fontTex, 0.5, 0, 0 );
	*/

	if(autoPilot){
        //printf("autoPiloting frame %i\n", frameCount);
        autoPilot1->control(world->dt); return;
    }

    if(mouseSteer){
        if (first_person){
            double dpitch=mouseY*0.005;
            double dyaw  =mouseX*0.002;
            double droll =0.2*dyaw;
            resetSteer( );
            myCraft->steerTo(droll, dpitch , dyaw);
        }else{
            Mat3d matCam;
            qCamera.toMatrix_T( matCam );
            Draw3D::drawMatInPos(matCam, myCraft->pos);
            steerToDir( matCam.c );
        }
    }

};

void AeroCraftGUI::drawHUD(){
	//panel .tryRender();  panel.draw();
	//mpanel.tryRender(); mpanel.draw();
	//if(focused) Draw2D::drawRectangle(focused->xmin,focused->ymin,focused->xmax,focused->ymax,false);

	if(staticTest){ drawStaticTest2D( *tester, fontTex, WIDTH, HEIGHT ); return; }

	char str[256];
	//sprintf(str, "speed %3.3f attitude %4.3f glideRatio %3.3f \0",world->myCraft->vel.norm(), world->myCraft->pos.y,   -sqrt(sq(world->myCraft->vel.x)+sq(world->myCraft->vel.z))/world->myCraft->vel.y );
	double vtot   = myCraft->vel.norm();
	//double thrust = world->myCraft->propelers[0].getThrust( vtot );
	double thrust = myCraft->totalThrust.norm();
	sprintf(str, "attitude %4.3f speed %3.3f vVert %3.3f tgAlfa %3.3f thrust %3.3f \0", myCraft->pos.y, vtot, myCraft->vel.y, myCraft->vel.y/vtot, thrust );
	glColor4f(1.0f,1.0f,1.0f,0.9f); Draw::drawText( str, fontTex, 10, 0,0 );

	if(first_person){ glColor4f(1.0f,1.0f,1.0f,0.9f); Draw2D::drawPointCross({mouseX+WIDTH*0.5,mouseY+HEIGHT*0.5},5.0); }

	/*
	int npol = 101;
	double phi0 =0.0;
	double dphi =3.14/(npol-1);
    for( int i=0; i<npol; i++){
        double phi = phi0 + dphi*i;
        double sa = sin( phi );
        double ca = cos( phi );
        double CD,CL;
        world->myCraft->panels[0].polarModel( ca, sa, CD, CL );
        printf("%i: %3.3f (%3.3f,%3.3f) (%3.3f,%3.3f) \n", i, phi, ca, sa, CD, CL);
    };
    exit(0);
    */

}

//void AeroCraftGUI:: drawHUD(){};

//AeroCraftGUI:: AeroCraftGUI( int& id, int WIDTH_, int HEIGHT_ ) : ScreenSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {
AeroCraftGUI:: AeroCraftGUI( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {
    //fontTex = makeTexture( "common_resources/dejvu_sans_mono.bmp" );
    //fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA.bmp" );
    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    //fontTex = makeTexture( "common_resources/dejvu_sans_mono_Alpha.bmp" );

    panel.init( 5,5,105,35,  fontTex );
    panel.caption   = "rotation [Rad]"; panel.vmin = -3.14159265359; panel.vmax = 3.14159265359;

    mpanel.initMulti( 120,5,200,120, fontTex , 4 );
    mpanel.caption="MultiPanel_1";

    txt.inputText = "insert number using =+-*/";

    SDL_StartTextInput ();

    // ======= from AeroCraftWorld

	//makeEnvironment( 20000.0f );
	printf( " Environment DONE! \n" );
	//makeAeroCraft();

	world = new Shooter();

    world->perFrame = 1;
    world->dt       = 0.005d;

	printf( "DEBUG 1 " );


	//AeroCraftWarrior = AeroCraftWarrior();

	myCraft_bak = new AeroCraft();          myCraft_bak->fromFile("data/AeroCraft1.ini");
	myCraft     = new AeroCraftWarrior();   myCraft    ->fromFile("data/AeroCraft1.ini");
	//myCraft     = new AeroCraft();   myCraft->fromFile("data/AeroCraft1.ini");
    myCraft->pos.y=200.0;
    myCraft->vel.set_mul( myCraft->rotMat.c, 100.0 );

    world->registrWarrior(myCraft);

    printf( "DEBUG 2 " );

    autoPilot1  = new AeroCraftControler();
    autoPilot1->craft=myCraft;

    printf( "DEBUG 3 " );

    tester      = new AeroTester();
    tester->autoPilot1 = autoPilot1;
    tester->myCraft    = myCraft;
    tester->gravityG   = world->gravity;
    //tester->reallocateTrj(int n);

    staticTest=false;
    if( staticTest ) tester->doStaticTesting( 500, 0.01, 300.0, 5.0 );

    printf( " AeroCraft DONE! \n" );

    // ======= fmakeEnvironment( 20000.0f );

	//buildings_shape = makeBuildings( 10, 10, 100, 100, 0.5, 0.5, 50, 100 );
	float sz = 1000.0;
	//buildings_shape = makeBuildingsClusters( 30, 3, 10,   -sz,         sz,         -sz,          sz,            0, 500,   20, 100,   10, 100 );
	//buildings_shape = makeBuildingsClusters( 30, 3, 10, -VIEW_DEPTH/2, VIEW_DEPTH/2, -VIEW_DEPTH/2, VIEW_DEPTH/2,    0, 500,   20, 100,   10, 100 );
	//buildings_shape= makeBuildingsClusters( 100, 5, 5, -VIEW_DEPTH/2, VIEW_DEPTH/2, -VIEW_DEPTH/2, VIEW_DEPTH/2,    0, 500,   100, 100,   10, 100 );

	double h0    = 1;
	//float tersz = VIEW_DEPTH/2;
	//terrain   = FieldPatch::makeList( 15, { 0.5,   -tersz,-tersz,h0,   tersz,-tersz,h0,  -tersz,tersz,h0,   tersz,tersz,h0   }   );

	//terrain_shape     = fieldPatch.makeList( 5, { 0.5d,   {-sz,-sz,h0},  { sz,-sz,h0},  {-sz,sz,h0},   {sz,sz,h0}   }   );
	//terrain_shape   = FieldPatch::makeList( 15, { 0.5,   Vec3d(-tersz,-tersz,h0),  Vec3d( tersz,-tersz,h0),  Vec3d(-tersz,tersz,h0),   Vec3d(tersz,tersz,h0)   }   );

    printf( " Environment DONE! \n" );

    //panel.nChars = 6;
};

// ===============================
// ===== GLOBAL VARIABLES
// ===============================

//AeroCraftGUI* thisScreen;

// ====================================
// ===== FUNCTION FORWARD DECLARATIONS
// ====================================



// ===============================
// ===== FUNCTION IMPLEMENTATION
// ===============================

// FUNCTION ======	camera

void AeroCraftGUI:: eventHandling   ( const SDL_Event& event  ){
    switch( event.type ){
        case SDL_KEYDOWN :
        switch( event.key.keysym.sym ){
            case SDLK_ESCAPE   : SDL_Quit(); exit(1); break;
            case SDLK_KP_PLUS  : zoom/=VIEW_ZOOM_STEP; printf("zoom: %f \n", zoom); break;
            case SDLK_KP_MINUS : zoom*=VIEW_ZOOM_STEP; printf("zoom: %f \n", zoom); break;
            case SDLK_SPACE    : STOP = !STOP; printf( STOP ? " STOPED\n" : " UNSTOPED\n"); break;

            case SDLK_u : autoPilot    = !autoPilot;    break;
            case SDLK_p : first_person = !first_person; break;
            case SDLK_m : mouseSteer   = !mouseSteer;   break;
            case SDLK_c : resetSteer( );
        }; break;
        case SDL_QUIT: SDL_Quit(); exit(1); break;

        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_RIGHT:
                    mouseSteer = true;
                    break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_RIGHT:
                    mouseSteer = false;
                    resetSteer();
                    break;
            }
            break;
    }
};

void AeroCraftGUI:: keyStateHandling( const Uint8 *keys ){
	//Uint8 *keystate = SDL_GetKeyState(NULL);
	//const Uint8 *keys = SDL_GetKeyboardState(NULL);

	if ( keys[ SDL_SCANCODE_DOWN  ] ) { qCamera.pitch( -0.005 );  }
	if ( keys[ SDL_SCANCODE_UP    ] ) { qCamera.pitch(  0.005  );  }
	if ( keys[ SDL_SCANCODE_RIGHT ] ) { qCamera.yaw  (  0.005  );  }
	if ( keys[ SDL_SCANCODE_LEFT  ] ) { qCamera.yaw  ( -0.005 ); }

	float dpitch = 0.01;
	float droll  = 0.002;
	float dyaw   = 0.01;

	if      ( keys[ SDL_SCANCODE_A ] ){ myCraft->panels[0].lrot.rotate( +droll, { 1,0,0 } );  myCraft->panels[1].lrot.rotate( -droll, { 1,0,0 } );    }
	else if ( keys[ SDL_SCANCODE_D ] ){ myCraft->panels[0].lrot.rotate( -droll, { 1,0,0 } );  myCraft->panels[1].lrot.rotate( +droll, { 1,0,0 } );    }

    if      ( keys[ SDL_SCANCODE_W ] ){ myCraft->panels[2].lrot.rotate( +dpitch, { 1,0,0 } );  }
	else if ( keys[ SDL_SCANCODE_S ] ){ myCraft->panels[2].lrot.rotate( -dpitch, { 1,0,0 } );  }

    if      ( keys[ SDL_SCANCODE_Q ] ){ myCraft->panels[3].lrot.rotate( +dyaw, { 0,1,0 } );  }
	else if ( keys[ SDL_SCANCODE_E ] ){ myCraft->panels[3].lrot.rotate( -dyaw, { 0,1,0 } );  }
    //if ( keystate[ SDL_SCANCODE_DOWN  ] ) { qmouse.pitch2( -0.005 ); }
	//if ( keystate[ SDL_SCANCODE_UP    ] ) { qmouse.pitch2( 0.005 ); }
	//if ( keystate[ SDL_SCANCODE_RIGHT ] ) { qmouse.yaw2  ( 0.005 ); }
	//if ( keystate[ SDL_SCANCODE_LEFT  ] ) { qmouse.yaw2  ( -0.005 ); }
};

void AeroCraftGUI:: mouseHandling   ( ){
	// mouse Camera
	int mx,my;
	SDL_GetMouseState(&mx,&my);
	int dmx = mx - WIDTH/2; 	int dmy = my - HEIGHT/2 ;
	mouseX = dmx;
	mouseY = -dmy;
	//printf( " mx: %i  my: %i dmx: %i dmy: %i ",mx, my, dmx, dmy );
	//qmouse.pitch( 0.001* dmy );
	//qmouse.yaw  ( 0.001* dmx );

	SDL_GetRelativeMouseState(&dmx,&dmy);
	qCamera.pitch( 0.005* dmy );
	qCamera.yaw  ( 0.005* dmx );
};

/*

void quit(){SDL_Quit(); exit(1);};
void setup();
void inputHanding();

// FUNCTION ======	setup
void setup(){
	srand(1234);
    world.init();
    thisScreen->world = &world;
    thisScreen->qCamera.setOne();
    thisScreen->VIEW_DEPTH = 100000;

    thisScreen->VIEW_DEPTH = 10000.0f;
    thisScreen->first_person = false;
    printf( " world.init(); DONE! \n" );

    world.fontTex_DEBUG = thisScreen->fontTex;
    SDL_ShowCursor( SDL_FALSE );
}

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
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int sid;
	//thisScreen  = new Screen2D( sid, 800,600);
	thisScreen  = new AeroCraftGUI( sid, 800,600 );

	setup();

	//loop( 1 );
	loop( 1000000 );

	return 0;
}
*/

AeroCraftGUI * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new AeroCraftGUI( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}


