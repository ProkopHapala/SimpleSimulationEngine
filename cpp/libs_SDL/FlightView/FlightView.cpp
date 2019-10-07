
// see
// https://www.opengl.org/discussion_boards/showthread.php/171319-glFlush-or-glFinish-with-mulithreading

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <thread>
#include <vector>

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "SDL_utils.h"
#include "Draw.h"
#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"

#include "AeroSurf.h"
#include "AeroCraft.h"
#include "geom3D.h"

#include "ShooterCommon.h"
#include "Projectile3D.h"

#include "Draw.h"
#include "Draw3D.h"
#include "AeroDraw.h"

#include "libFlight.h"

char* work_dir = 0;

int camType = 2;

FlightWorld* world = 0;
std::vector<AeroSurfaceDebugRecord> dbgRects;
char str[2<<16];

void drawBurst( Burst3d& burst, int nseg, float shotSz, bool bBox, double dt ){
    glColor3f(1.0,1.0,0.0);
    int n = burst.shots.size();
    Vec3d op;
    for(int i=0; i<n; i++){
        Draw3D::drawPointCross( burst.shots[i].pos, shotSz );
        //if( tmpPos ) Draw3D::drawLine( tmpPos[i], burst.shots[i].pos );
        burst.shots[i].getOldPos(dt,op);
        Draw3D::drawLine( op, burst.shots[i].pos );
    }
    //if( tmpPos ){ glColor3f(1.0,1.0,0.0); Draw3D::drawPointCross( tmpPos[n-1], 0.5 ); }
    if(bBox){
        burst.shots[0].getOldPos(dt,op);
        glColor3f(1.0,1.0,0.0); Draw3D::drawPointCross( op, 0.5 );
        glColor3f(0.0,1.0,1.0); Draw3D::drawPointCross( burst.shots[0].pos, 0.5 );
        glColor3f(0.0,0.5,0.0);
        Vec3d b = burst.bbox.p+burst.bbox.hdir*burst.bbox.l;
        Draw3D::drawPointCross( burst.bbox.p, 3.0 );
        Draw3D::drawPointCross( b           , 3.0 );
        Draw3D::drawCylinderStrip_wire( nseg, burst.bbox.r, burst.bbox.r, (Vec3f)burst.bbox.p, (Vec3f)b );
    }
}


class AppFlightView : public AppSDL2OGL_3D { public:

    int fontTex   = 0;
    int txGround  = 0;
    int gloTarget = 0;
    AeroCraft *myCraft = 0;
    //GUI gui;




	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	//virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	void writeAeroCraftStateHUD();

	AppFlightView( int& id, int WIDTH_, int HEIGHT_ );

};

AppFlightView::AppFlightView( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    myCraft = &world->craft1;
    dbgRects.resize( myCraft->nPanels );
    for(int i=0; i<myCraft->nPanels; i++){
        myCraft->panels[i].dbgRec=&dbgRects[i];
    };

    gloTarget=glGenLists(1);
	glNewList( gloTarget, GL_COMPILE );
        //glEnable    ( GL_LIGHTING );
        //glShadeModel( GL_FLAT     );
        glDisable( GL_LIGHTING );
        //Draw3D::drawSphereOctLines( 2, 1.0, {0.0,0.0,0.0} );
        Draw3D::drawSphere_oct( 2, 1.0, {0.0,0.0,0.0}, true );
	glEndList();

	sprintf(str, "%s/common_resources/dejvu_sans_mono_RGBA_pix.bmp", work_dir );
	printf( "FlightView load from:  %s \n", str  );
    fontTex = makeTextureHard( str );
    //GUI_fontTex = fontTex;

    int imgW=1024,imgH=1024;
    uint8_t* pixels = new uint8_t[imgH*imgW*4];
    for(int iy=0; iy<imgH; iy++){ for(int ix=0; ix<imgW; ix++){
        int i = (iy*imgW+ix)*4;
        pixels[i+0] = ix ^ iy;
        pixels[i+1] = 128;
        pixels[i+2] = ix ^ iy;
        pixels[i+3] = 255;
    } }

    glEnable(GL_TEXTURE_2D);
    glActiveTexture(GL_TEXTURE0);
    glGenTextures  ( 1, (GLuint*)&txGround  );
    glBindTexture  ( GL_TEXTURE_2D, txGround  );
    //glTexImage2D(GL_TEXTURE_2D, 0, GL_R, imgW, imgH, 0,  GL_R, GL_UNSIGNED_BYTE, pixels );
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, imgW, imgH, 0,  GL_RGBA, GL_UNSIGNED_BYTE, pixels );
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    //glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    //glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glDisable(GL_TEXTURE_2D);
    delete [] pixels;

    zoom = 30.0;
    VIEW_DEPTH = 50000;
};

void AppFlightView::draw(  ) {
    //glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    glClearColor( 0.5f, 0.5f, 0.7f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);

    glMatrixMode(GL_PROJECTION); //glLoadIdentity();

    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    glDisable(GL_CULL_FACE);

    //AeroCraft* craft = &world->craft1;

    //camera_FPS ( craft->pos, craft->rotMat );

    /*
    switch( camType ){
        case 1:
            qCamera.toMatrix( cam.rot );
            cam.pos = (Vec3f)myCraft->pos;
            camera_FwUp( myCraft->pos, (Vec3d)cam.rot.c, {0.0,1.0,0.0}, false );
            break;
        case 2:
            camera_FPS( myCraft->pos, myCraft->rotMat );
            break;
    }
    */

    int mx,my; Uint32 buttons = SDL_GetRelativeMouseState( &mx, &my);
    //printf( " %i %i \n", mx,my );
    if ( buttons & SDL_BUTTON(SDL_BUTTON_RIGHT) ) {
        Quat4f q; q.fromTrackball( 0, 0, -mx*mouseRotSpeed, my*mouseRotSpeed );
        qCamera.qmul_T( q );
        qCamera.toMatrix( cam.rot );
        cam.pos = (Vec3f)myCraft->pos;
        camera_FwUp( myCraft->pos, (Vec3d)cam.rot.c, {0.0,1.0,0.0}, false );
    }else{
        camera_FPS( myCraft->pos, myCraft->rotMat );
        qCamera.fromMatrix( (Mat3f)myCraft->rotMat );
    }


    glMatrixMode(GL_MODELVIEW); glLoadIdentity();

    //Draw3D::drawMatInPos( Mat3dIdentity*30.0, myCraft->pos );

    renderAeroCraft       ( *myCraft, true, 1.0 );
    //renderAeroCraftVectors( *myCraft, 0.01, 1.0, true );

    glColor3f( 0.5, 0.7, 0.5 );
    glEnable     ( GL_TEXTURE_2D );
    glBindTexture( GL_TEXTURE_2D, txGround );
    Draw3D::drawPanel( Vec3dZero, Mat3dIdentity, {50000.0,50000.0} ); // ground
    glDisable    ( GL_TEXTURE_2D );

    //printf("nTargets %i \n", world->nTargets);
    glColor3d(0.0,0.0,0.0);
    for(int i=0; i<world->nTargets; i++ ){
        //printf("%i (%g,%g,%g) %g \n", i, world->targets[i].p.x, world->targets[i].p.y, world->targets[i].p.z, world->targets[i].r );
        Draw3D::drawShape( world->targets[i].p, Mat3dIdentity*world->targets[i].r, gloTarget, false );
    }

    //glEnable(GL_BLEND);
    for(Burst3d& burst : world->bursts){
        //drawBurst( burst, 8, 1.0, true, world->dt );
        drawBurst( burst, 8, 1.0, false, world->dt );
    }

};

void AppFlightView::drawHUD(  ) {
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);

    //gui.draw();
    //drawControlSurfaceStateHUD();
    //drawPolarPlotHUD();
    //drawAeroCraftTrjHUD();
    writeAeroCraftStateHUD();
    //drawCompassHUD( myCraft->pos.xz(), myCraft->rotMat.c.xz() );
}


void AppFlightView::writeAeroCraftStateHUD(){
	Vec3d& pos    = myCraft->pos;
	Vec3d& vel    = myCraft->vel;
    double vtot   = vel.norm();
    double thrust = myCraft->totalThrust.norm();
    double power  = myCraft->getTotalPower();
    Vec3d vdir  = vel*(1.0/vtot);
    char* s=str;
    s+=sprintf(s,"Time[s]: %3.1f\n",              world->time );
    s+=sprintf(s,"altitude[m  ] %g\n",            pos.y  );
    s+=sprintf(s,"speed   [m/s] %g [km/h] %g \n", vtot, vtot*3.6 );
    s+=sprintf(s,"climb   [m/s] %g           \n", vel.y  );
    s+=sprintf(s,"thrust  [N  ] %g [kg]   %g \n", thrust, thrust/9.81   );
    s+=sprintf(s,"thrust0 [N  ] %g P[kW]  %g\n", power/vtot, power*1e-3 );
    s+=sprintf(s,"pos  (%g,%g,%g) \n", pos .x, pos .y, pos .z  );
    s+=sprintf(s,"vel  (%g,%g,%g) \n", vel .x, vel .y, vel .z  );
    s+=sprintf(s,"vdir (%g,%g,%g) \n", vdir.x, vdir.y, vdir.z );

    s+=sprintf(s,"--- control --- \n", vdir.x, vdir.y, vdir.z );
    s+=sprintf(s,"pitch    %g \n", world->controlsState[iPitch].x );
    s+=sprintf(s,"yaw      %g \n", world->controlsState[iYaw  ].x );
    s+=sprintf(s,"roll     %g \n", world->controlsState[iRoll ].x );
    s+=sprintf(s,"throttle %g \n", world->controlsState[iThrottle ].x );

    //s+=sprintf(s,"--- Pro --- \n", vdir.x, vdir.y, vdir.z );
    int nprj=0;
    for(Burst3d& burst : world->bursts){ nprj += burst.shots.size(); };
    s+=sprintf(s," N burst       %i \n",  world->bursts.size() );
    s+=sprintf(s," N projectiles %i \n", nprj );

    //s+=sprintf(s,"vdir (%g,%g,%g) \n", vdir.x, vdir.y, vdir.z );
    glPushMatrix();
    glColor4f(1.0f,1.0f,1.0f,0.9f);
    glTranslatef( 10.0, HEIGHT-20, 0.0 );
    Draw::drawText( str, fontTex, fontSizeDef, {40,40} );
    glPopMatrix();
}


AppFlightView * thisApp = 0;

extern "C"{

void init( int w, int h, void* world_, char* work_dir_ ){
    world = (FlightWorld*)world_;
    work_dir = work_dir_;

    SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	//SDL_ShowCursor(SDL_DISABLE);
    int junk;
	thisApp = new AppFlightView( junk , w, h );
	SDL_SetWindowPosition(thisApp->window, 100, 0 );

}

bool draw(){
    //thisApp->update();
    if( thisApp->GL_LOCK ) return true;
    SDL_GL_MakeCurrent(thisApp->window, thisApp->glctx);
    thisApp->inputHanding();
    thisApp->GL_LOCK = true;
    thisApp->draw();
    thisApp->cameraHUD();
    thisApp->drawHUD();
    //SDL_RenderPresent(thisApp->renderer);
    SDL_GL_SwapWindow(thisApp->window);
    thisApp->GL_LOCK = false;
    return thisApp->GL_LOCK;

}

}
