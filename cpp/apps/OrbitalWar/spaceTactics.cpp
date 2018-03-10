
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "Draw3D.h"
#include "SDL_utils.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"

#include "spline_hermite.h"
#include "SplineManager.h"


#include "SpaceBodies.h"
#include "SpaceWorld.h"

#include "AppSDL2OGL_3D.h"
#include "GUI.h"
#include "testUtils.h"

/*


#### Ideas:
 -

 - Optimal placement of units is given by
    - effective range of weapons and
    - narrow angle of shields

    *       *       *        *   //army 1
    ..
    . .
    .  .
    .   .
    .    .
    .     .
    .      .
    *       *       *        *   //army 1

*/

void drawPolyLineRef( int n, Vec3d * ps, Vec3d * ps_ref ){   // closed=false
    //printf("%i %i\n", n, closed );
    glBegin(GL_LINE_STRIP);
    for(int i=0; i<n; i++){
        glVertex3d( ps[i].x-ps_ref[i].x, ps[i].y-ps_ref[i].y, ps[i].z-ps_ref[i].z );
    };
    glEnd();
};


class SpaceTactics : public AppSDL2OGL_3D { public:
    SplineManager splines;

    SpaceWorld world;

    bool dragging=false;
    int     nsamp=0;
    double tstart, tend, dt;

    /*
    bool dragging = false;
    int iedit=0, ipoint=0;
    double mouse_t,mouse_val;
    double timeScale=100.0,valScale=100.0f;
    */


    SpaceBody* referenceBody = 0;

    int          fontTex;
    GUITextInput txtStatic;
    //GUITextInput txtPop;
    char curCaption[100];

    //DEGUB
    //Vec3d *  dpos;
    //Vec3d * ddpos;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	//virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );
    virtual void mouseHandling( );

	SpaceTactics( int& id, int WIDTH_, int HEIGHT_ );

};

SpaceTactics::SpaceTactics( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    world.addPlanet( "Sun"  , 100000000000.0, 1.0, { 0.0, 0.0, 0.0}, {0.0,0.0,0.0} );
    world.addPlanet( "Earth",   1000000000.0, 0.1, {10.0, 0.0, 0.0}, {0.0,1.0,0.0} );

    world.allocateODE ();
    world.allocateTrjs( 5000 );
    //world.allocateTrjs( 20 );

    referenceBody = &world.planets[0];

    world.predictTrjs( 100000000, 0.1 );
    //exit(0);

    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    //txtPop   .inputText = "txtPop";
    txtStatic.inputText = "txtStatic";

    zoom = 20;

}

void SpaceTactics::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glColor3f(0.0f,0.0f,0.0f);
	int i=0;
	for(SpaceBody& b : world.planets ){
        //Draw3D::drawPolyLine( world.trj_n, b.trjPos, false );
        drawPolyLineRef( world.trj_n, b.trjPos, referenceBody->trjPos );
        //if(i==1) for(int j=0; j<world.trj_n; j++ ) printf( " %i %f %f %f \n", j, b.trjPos[j].x, b.trjPos[j].y, b.trjPos[j].z );
        i++;
	}
    for(SpaceBody& b : world.ships ){
        drawPolyLineRef( world.trj_n, b.trjPos, referenceBody->trjPos );
	}
    //exit(0);

};

void SpaceTactics::drawHUD(){
    glDisable ( GL_LIGHTING );

    glColor3f(1.0f,1.0f,1.0f);   txtStatic.view3D( {5,5}, fontTex, 8 );

}

//void SpaceTactics::keyStateHandling( const Uint8 *keys ){ };

void SpaceTactics::mouseHandling( ){
    /*
    SDL_GetMouseState( &mouseX, &mouseY ); mouseY=HEIGHT-mouseY;
    mouse_t   = (mouseX-20 )/timeScale + tstart;
    mouse_val = (mouseY-200)/valScale;
    sprintf( curCaption, "%f %f\0", mouse_t, mouse_val );
    int ipoint_ = binSearchFrom<double>(mouse_t,splines.n,splines.ts);
    if( (splines.ts[ipoint_+1]-mouse_t)<(mouse_t-splines.ts[ipoint_]) ) ipoint_++;
    if(ipoint_!=ipoint){
        ipoint=ipoint_;
        char buff[100];
        Vec3d r,v;
        r.set( splines.CPs[0][ipoint],splines.CPs[1][ipoint],splines.CPs[2][ipoint] );
        v.set( splines.getPointDeriv(ipoint,0), splines.getPointDeriv(ipoint,1), splines.getPointDeriv(ipoint,2) );
        sprintf(buff, "%i %f r(%3.3f,%3.3f,%3.3f) v(%3.3f,%3.3f,%3.3f)", ipoint, splines.ts[ipoint], r.x,r.y,r.z,   v.x,v.y,v.z );
        txtStatic.inputText = buff;
    }
    //printf("%i %i %3.3f %3.3f \n", mouseX,mouseY, mouse_t, mouse_val );
    */
}

/*
void SpaceTactics::eventHandling ( const SDL_Event& event  ){
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

                //case SDLK_x:  iedit=0; break;
                //case SDLK_y:  iedit=1; break;
                //case SDLK_z:  iedit=2; break;
            }
            break;

        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                case SDL_BUTTON_RIGHT:
                    //float t   = (mouseX+20 )/timeScale + tstart;
                    //float val = (mouseY-200)/valScale;
                    //ipoint = binSearchFrom<double>(mouse_t,splines.n,splines.ts);
                    dragging=true;
                    break;
            }
            break;

        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    if(dragging){
                    }
                    dragging=false;
                    break;
                case SDL_BUTTON_RIGHT:
                    if(dragging){

                    }
                    dragging=false;
                    break;
            }
            break;
    };
    //AppSDL2OGL::eventHandling( event );
}
*/



// ===================== MAIN

SpaceTactics * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new SpaceTactics( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















