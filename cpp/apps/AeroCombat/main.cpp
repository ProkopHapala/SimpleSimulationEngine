

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

//#include "SimplexRuler.h"
//#include "Ruler2DFast.h"
//#include "TerrainHydraulics.h"

#include "Terrain25D.h"
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

Shooter            * world      = NULL;
AeroCraftControler * autoPilot  = NULL;
AeroTester         * tester     = NULL;

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

// ===============================
// ===== Free Functions
// ===============================

void rotateTo( int pivot, Mat3d& rot, const Mat3d& rot0, const Vec3d& xhat, const Vec3d& yhat, double sin_angle ){
    //rot.add_mul( rot0, coef ); rot.normalize();
    double ca,sa;
    int i3 = pivot*3;
    Vec3d& piv  = *(Vec3d*)(rot .array+i3);
    Vec3d& piv0 = *(Vec3d*)(rot0.array+i3);
    piv.getInPlaneRotation( piv0, xhat, yhat, ca, sa );
    if( sa < sin_angle ){
        Vec3d uaxis; uaxis.set_cross( xhat, yhat );
        rot.rotate_csa( ca, sa, uaxis );
    }else{
        rot.set(rot0);
    }
}

Terrain25D * prepareTerrain( int nsz, int nsub, double step, double hmax ){
    Terrain25D_bicubic * terrain = new Terrain25D_bicubic();
    terrain->ruler.setup( (Vec2d){nsz*0.5,nsz*0.5}*-step, (Vec2d){step,step} );
    terrain->allocate( {nsz,nsz} );
    terrain->makeRandom( 0.0, hmax );

    terrain->shape = glGenLists(1);
    glNewList( terrain->shape , GL_COMPILE );
    //int na=100,nb=100;

    int na = (terrain->ruler.n.a - 3)*nsub;
    int nb = (terrain->ruler.n.b - 3)*nsub;
    float da=terrain->ruler.step.a/float(nsub);
    float db=terrain->ruler.step.b/float(nsub);
    float x0=terrain->ruler.pos0.x;
    float y0=terrain->ruler.pos0.y;

    glEnable(GL_LIGHTING);
    glColor3f (0.5f,0.5f,0.5f);
    glNormal3f(0.0f,1.0f,0.0f);

    float * oldvals = new float[na*3];
    for(int ia=0; ia<na; ia++){
        glBegin(GL_TRIANGLE_STRIP);
        for(int ib=0; ib<nb; ib++){
            int i3 = 3*ib;
            Vec2d dv1,dv2;
            Vec2d p1; p1.set( (ia  )*da+x0, ib*db+y0 );
            Vec2d p2; p2.set( (ia+1)*da+x0, ib*db+y0 );
            float v1,v2;
            if( ia == 0 ){
                v1 = (float)terrain->eval( p1, dv1 );
            }else{
                v1 = oldvals[i3]; dv1.x=oldvals[i3+1]; dv1.y=oldvals[i3+2];
            }
            v2 = (float)terrain->eval( p2, dv2 );
            oldvals[i3] = v2; oldvals[i3+1] = dv2.x; oldvals[i3+2] = dv2.y;
            glNormal3f(-dv1.x,1.0,-dv1.y); glVertex3f( (float)p1.x,  v1, (float)p1.y );
            glNormal3f(-dv2.x,1.0,-dv2.y); glVertex3f( (float)p2.x,  v2, (float)p2.y );

            //glColor3f(v1,0.5,-v1); glVertex3f( (float)p1.x,  v1, (float)p1.y );
            //glColor3f(v2,0.5,-v2); glVertex3f( (float)p2.x,  v2, (float)p2.y );
            //printf( " %i (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f)\n", p1.x, p1.y, v1 ,  p2.x, p2.y, v2  );
        }
        glEnd();
    }

    // Normals
    /*
    glBegin(GL_LINES);
    for(int ia=0; ia<na; ia++){
        for(int ib=0; ib<nb; ib++){
            int i3 = 3*ib;
            Vec2d p,dv; p.set( ia*da+x0, ib*db+y0 );
            double v = (float)terrain->eval( p, dv );
            glVertex3f( (float)p.x,         v, (float)p.y );
            glVertex3f( (float)(p.x-dv.x),  v+1.0, (float)(p.y-dv.y) );
        }

    }
    glEnd();
    */
    glEndList();
    delete [] oldvals;
    return terrain;
}


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

    bool mouseSteer   = false;
    bool autoRetractAirelon = true;
    bool useAutoPilot = false;

	int perFrame = 10;
	//double dt = 0.001;

	AeroCraftWarrior * myCraft;
	//AeroCraft * myCraft;
	AeroCraft * myCraft_bak;

	//FieldPatch fieldPatch;
	//int buildings_shape = -1;
	//int terrain_shape   = -1;

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

	//renderSkyBox(myCraft->pos.x, myCraft->pos.y-1000, myCraft->pos.z, VIEW_DEPTH*0.25 );
	glEnable(GL_DEPTH_TEST);

	glEnable    (GL_LIGHTING);
	glShadeModel(GL_FLAT);

	//world->myCraft->render();
	renderAeroCraft( *myCraft );

	//glDisable (GL_LIGHTING);
	glShadeModel(GL_SMOOTH);

	/*

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

	*/

	if(world->terrain) glCallList(world->terrain->shape);

	//Draw3D::drawAxis( 1000 );

	if(useAutoPilot){
        //printf("autoPiloting frame %i\n", frameCount);
        autoPilot->control(world->dt); return;
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
}

AeroCraftGUI:: AeroCraftGUI( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    printf( " === GUI \n" );

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

    printf( " === world  \n" );

	world = new Shooter();
    world->perFrame = 1;
    world->dt       = 0.005d;

    printf( " === Environment \n" );
    world->terrain =  prepareTerrain( 128, 2, 100.0, 50 );

    printf( " === aerocraft \n" );

	myCraft_bak = new AeroCraft();          myCraft_bak->fromFile("data/AeroCraft1.ini");
	myCraft     = new AeroCraftWarrior();   myCraft    ->fromFile("data/AeroCraft1.ini");
	//myCraft     = new AeroCraft();   myCraft->fromFile("data/AeroCraft1.ini");
    myCraft->pos.y=200.0;
    myCraft->vel.set_mul( myCraft->rotMat.c, 100.0 );
    world->registrWarrior(myCraft);

    printf( " === autoPilot1 \n" );

    autoPilot  = new AeroCraftControler();
    autoPilot->craft=myCraft;

    printf( " === tester \n" );

    tester      = new AeroTester();
    tester->autoPilot = autoPilot;
    tester->craft     = myCraft;
    tester->gravityG  = world->gravity;
    //tester->reallocateTrj(int n);
    staticTest=false;
    if( staticTest ) tester->doStaticTesting( 500, 0.01, 300.0, 5.0 );

};

void AeroCraftGUI:: eventHandling   ( const SDL_Event& event  ){
    switch( event.type ){
        case SDL_KEYDOWN :
        switch( event.key.keysym.sym ){
            case SDLK_ESCAPE   : SDL_Quit(); exit(1); break;
            case SDLK_KP_PLUS  : zoom/=VIEW_ZOOM_STEP; printf("zoom: %f \n", zoom); break;
            case SDLK_KP_MINUS : zoom*=VIEW_ZOOM_STEP; printf("zoom: %f \n", zoom); break;
            case SDLK_SPACE    : STOP = !STOP; printf( STOP ? " STOPED\n" : " UNSTOPED\n"); break;

            case SDLK_u : useAutoPilot = !useAutoPilot;    break;
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

    if( autoRetractAirelon ){
        //myCraft->panels[0].lrot.rotateTo( myCraft_bak->panels[0].lrot, 0.01 );
        //myCraft->panels[1].lrot.rotateTo( myCraft_bak->panels[1].lrot, 0.01 );
        rotateTo( 1, myCraft->panels[0].lrot, myCraft_bak->panels[0].lrot, { 0,0,1 }, { 0,1,0 }, 0.01 );
        rotateTo( 1, myCraft->panels[0].lrot, myCraft_bak->panels[0].lrot, { 0,0,1 }, { 0,1,0 }, 0.01 );
    }

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

};

void AeroCraftGUI:: mouseHandling   ( ){
	// mouse Camera
	int mx,my;
	SDL_GetMouseState(&mx,&my);
	int dmx = mx - WIDTH/2; 	int dmy = my - HEIGHT/2 ;
	mouseX = dmx;
	mouseY = -dmy;

	SDL_GetRelativeMouseState(&dmx,&dmy);
	qCamera.pitch( 0.005* dmy );
	qCamera.yaw  ( 0.005* dmx );
};

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


