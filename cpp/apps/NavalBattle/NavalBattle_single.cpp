
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw3D.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "geom3D.h"
#include "Body2D.h"
#include "raytrace.h"

#include "Draw2D.h"
#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

#include "Object3D.h"
//#include "Warrior3D.h"
#include "Ship2D.h"
#include "Battleship.h"
#include "Projectile3D.h"
#include "Shooter.h"


const int BB_Type = 1;

class NavalBattle_single : public AppSDL2OGL_3D {
	public:
    Shooter world;
    double dvel = 10.0;

    Ship2D *ship1 = NULL;

    Vec3d crosshair;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	virtual void keyStateHandling( const Uint8 *keys );
    virtual void mouseHandling( );
    void         camera();

	NavalBattle_single( int& id, int WIDTH_, int HEIGHT_ );

};

NavalBattle_single::NavalBattle_single( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {
    world.init_world();
    world.wind_speed.set(0.0);
    world.watter_speed.set(0.0);

    ship1 = new Battleship();
    // ShipHull.obj exported from blender with Forward=-Y  Up=-Z
    ((Battleship*)ship1)->loadFromFile( "data/USS_Arizona.ini" );
    ((Battleship*)ship1)->render( );
    ship1->setAngle ( 1.0 );
    ship1->throttle = 0.0;

    world.warriors25D.push_back(ship1);

	printf( "==== INITIALIZATION DONE ====\n" );


	zoom = 500.0;
}

void NavalBattle_single::draw(){
    //delay = 200;
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    //glEnable( GL_LIGHTING );
    //glCallList( world.defaultObjectShape );

    Draw3D::drawPointCross( crosshair , 10.0 );

	delay          = 100;
	world.dt       = 0.1;
	world.perFrame = 1;
    world.update_world( );

    glEnable( GL_LIGHTING );

    float glMat[16];
    if (ship1){
        printf( "drawing ship! rudder %f pos (%3.3f,%3.3f) rot (%3.3f,%3.3f) vel (%3.3f,%3.3f) \n", ship1->rudder.phi, ship1->pos.x, ship1->pos.y,  ship1->rot.x, ship1->rot.y, ship1->vel.x, ship1->vel.y );
        //glCallList(ship1->shape);
        //Draw2D::drawShape( ship1->pos, ship1->rot, ship1->shape );
        //Draw3D::drawShape( const Vec3d& pos, const Mat3d& rot, int shape );

        //for( Turret* tur : ((Battleship*)ship1)->turrets ){
        //    Draw2D::drawVecInPos( {tur->gpos.x,tur->gpos.z}, {tur->grot.c.x,tur->grot.c.z} );
        //}

        ship1->draw();

        /*
        glPushMatrix();
        gpos
        grot
        Draw3D::toGLMat( ship1->gpos, ship1->grot, o->span, glMat );
        glMultMatrixf( glMat );
        glCallList( ship1->mesh->rendered_shape );
        glPopMatrix();
        */
    }

    /*
    for( auto p : world.projectiles ) {
        Draw3D::drawLine(p->pos, p->old_pos);
        Draw3D::drawPointCross(p->pos,0.1);
    }
    */

    glDisable ( GL_LIGHTING );
	Draw3D::drawAxis ( 3.0f );

};

void NavalBattle_single::keyStateHandling( const Uint8 *keys ){

    if( ship1 != NULL ){
       // if( keys[ SDL_SCANCODE_W ] ){ ship1->vel.add_mul( camMat.c, +0.1 ); }
       // if( keys[ SDL_SCANCODE_S ] ){ ship1->vel.add_mul( camMat.c, -0.1 ); }
       //if( keys[ SDL_SCANCODE_A ] ){ ship1->vel.add_mul( camMat.a, -0.1 ); }
       //if( keys[ SDL_SCANCODE_D ] ){ ship1->vel.add_mul( camMat.a, +0.1 ); }

        //if( keys[ SDL_SCANCODE_LEFT  ] ){ ship1->rudder.setAngle( ship1->rudder.phi + 0.01 );  }
        //if( keys[ SDL_SCANCODE_RIGHT ] ){ ship1->rudder.setAngle( ship1->rudder.phi - 0.01 );  }

        if( keys[ SDL_SCANCODE_LEFT  ] ){ ship1->setAngle( ship1->phi + 0.01 );  }
        if( keys[ SDL_SCANCODE_RIGHT ] ){ ship1->setAngle( ship1->phi - 0.01 );  }

        //if( keys[ SDL_SCANCODE_SPACE ] ){ warrior1->vel.mul( 0.9 ); }
    }

    //if( keys[ SDL_SCANCODE_W ] ){ camPos.add_mul( camMat.c, +0.1 ); }
	//if( keys[ SDL_SCANCODE_S ] ){ camPos.add_mul( camMat.c, -0.1 ); }
	//if( keys[ SDL_SCANCODE_A ] ){ camPos.add_mul( camMat.a, -0.1 ); }
	//if( keys[ SDL_SCANCODE_D ] ){ camPos.add_mul( camMat.a, +0.1 ); }
	//if( keys[ SDL_SCANCODE_W ] ){ camPos.z += +0.1; }
	//if( keys[ SDL_SCANCODE_S ] ){ camPos.z += -0.1; }
	//if( keys[ SDL_SCANCODE_A ] ){ camPos.x += +0.1; }
	//if( keys[ SDL_SCANCODE_D ] ){ camPos.x += -0.1; }
    if( keys[ SDL_SCANCODE_Q ] ){ qCamera.droll2(  +0.01 ); }
	if( keys[ SDL_SCANCODE_E ] ){ qCamera.droll2(  -0.01 ); }

};

void NavalBattle_single::mouseHandling( ){
    int mx,my;
    //Uint32 buttons = SDL_GetRelativeMouseState( &mx, &my);
    Uint32 buttons = SDL_GetMouseState(&mx, &my);
    crosshair.set((mx-WIDTH*0.5)*(2*zoom/HEIGHT),0.0,(my-HEIGHT*0.5)*(2*zoom/HEIGHT));
    ((Battleship*)ship1)->aim( crosshair, world.gravity );


    /*
    camPhi   += mx*0.001;
    camTheta += my*0.001;
    camMat.a.set( sin(camPhi),0.0,cos(camPhi) );
    double ct=cos(camTheta),st=sin(camTheta);
    camMat.b.set(-camMat.a.z*st, ct, camMat.a.x*st);
    camMat.c.set(-camMat.a.z*ct,-st, camMat.a.x*ct); // up vector
    //printf("camMat:\n",);
    printf("camMat: (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f)\n",camMat.ax,camMat.ay,camMat.az,  camMat.bx,camMat.by,camMat.bz, camMat.cx,camMat.cy,camMat.cz);
    //qCamera.fromAngleAxis( camTheta, {sin(camPhi),0, cos(camPhi)} );
    //qCamera.fromAngleAxis( camTheta, {0,sin(camPhi), cos(camPhi)} );
    */

}


void NavalBattle_single::camera(){

    //float camMatrix[16];
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity();
    float fov = VIEW_ZOOM_DEFAULT/zoom;
    //glFrustum( -ASPECT_RATIO, ASPECT_RATIO, -1, 1, 1*fov, VIEW_DEPTH*fov );
    glOrtho  ( -zoom*ASPECT_RATIO, zoom*ASPECT_RATIO, -zoom, zoom, -VIEW_DEPTH, +VIEW_DEPTH );
    //Draw3D::toGLMatCam( {0.0d,0.0d,0.0d}, camMat, camMatrix );

    glRotatef(90,1.0,0.0,0.0);
    //glMultMatrixf( camMatrix );
    //glTranslatef ( -camPos.x, -camPos.y, -camPos.z );

    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity();


}

void NavalBattle_single::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_f:  warrior1->tryJump(); break;
                //case SDLK_h:  warrior1->tryJump(); break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
                case SDLK_ESCAPE:   quit(); break;
                //case SDLK_SPACE:    STOP = !STOP; printf( STOP ? " STOPED\n" : " UNSTOPED\n"); break;
                case SDLK_KP_MINUS: zoom*=VIEW_ZOOM_STEP; break;
                case SDLK_KP_PLUS:  zoom/=VIEW_ZOOM_STEP; break;
            }
            break;

        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    ship1->trigger = true;
                    break;
            }
            break;

        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    ship1->trigger = false;
                    break;
            }
            break;
    };
    //AppSDL2OGL::eventHandling( event );
}




void NavalBattle_single::drawHUD(){
    glDisable ( GL_LIGHTING );
    glColor3f( 0.0f, 1.0f, 0.0f );
    drawCrosshair( 10 );
}

// ===================== MAIN

NavalBattle_single * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new NavalBattle_single( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















