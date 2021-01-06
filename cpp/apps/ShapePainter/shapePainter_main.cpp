
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "testUtils.h"

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "Draw3D.h"
#include "SDL_utils.h"

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
//#include "raytrace.h"

#include "lineSearch.h"

#include "Spline2d.h"
#include "geom2D.h"
#include "spline_Circle.h"


//#include "SphereSampling.h"
//#include "DrawSphereMap.h"
//#include "Draw3D_Surf.h"

#include "AppSDL2OGL_3D.h" // should we use 2D instead ? - no, we wan to draw projections later
#include "GUI.h"
#include "IO_utils.h"

#include "ShapePainter.h"

//#include "Table.h"
//#include "Tree.h"
//#include "IntersectionCurve.h"

//using namespace Painter;


int fontTex = 0 ;

class PainterGUI : public AppSDL2OGL_3D { public:

	//DropDownList lstLuaFiles;
    GUI gui;

    //Spline2d     spline1;
    //CircleSpline cspline;
    //IntersectionCurve impCurve;
    //DropDownList* lstLuaFiles=0;
    //OnSelectLuaShipScript onSelectLuaShipScript;

    //double gridStep = 1.0;

    //int ipick = -1;
    Vec2f  opmouse;


    CornerBrush cornerBrush;

    Brush* brush;

    //int ogl;

    // ======= Functions

    //void selectCompGui();

    // https://stackoverflow.com/questions/29145476/requiring-virtual-function-overrides-to-use-override-keyword

	virtual void draw   () override;
	virtual void drawHUD() override;
	//virtual void mouseHandling( );
	//virtual void camera();
	virtual void eventHandling   ( const SDL_Event& event  ) override;
	virtual void keyStateHandling( const Uint8 *keys )       override;
    //virtual void mouseHandling( );

	PainterGUI( int& id, int WIDTH_, int HEIGHT_ );

};

PainterGUI::PainterGUI( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {
    brush=&cornerBrush;
    //brush;
}



void PainterGUI::draw(){
    glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glDisable(GL_DEPTH_TEST);
	glEnable (GL_BLEND     );
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	brush->view();

    //glCallList(ogl);

};

void PainterGUI::drawHUD(){
    glDisable( GL_LIGHTING );
    glDisable(GL_DEPTH_TEST);

    //gui.draw();
    //glColor3f(1.0f,1.0f,1.0f);   txtStatic.view3D( {5,5}, fontTex, 8 );
    //glPopMatrix();
}

void PainterGUI::keyStateHandling( const Uint8 *keys ){

    //Mat3d camMat;
    //qCamera.toMatrix_T(camMat);
    //Draw3D::drawMatInPos(  (Mat3f)camMat, (Vec3f){0.0,0.0,0.0} );
    //qCamera.toMatrix(camMat);

    if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.dyaw  (  keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.dyaw  ( -keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.dpitch(  keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.dpitch( -keyRotSpeed ); }

    if( keys[ SDL_SCANCODE_W ] ){ cam.pos.add_mul( cam.rot.b, +0.05*zoom ); }
	if( keys[ SDL_SCANCODE_S ] ){ cam.pos.add_mul( cam.rot.b, -0.05*zoom );  }
	if( keys[ SDL_SCANCODE_A ] ){ cam.pos.add_mul( cam.rot.a, -0.05*zoom );  }
	if( keys[ SDL_SCANCODE_D ] ){ cam.pos.add_mul( cam.rot.a, +0.05*zoom );  }

};

void PainterGUI::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    //if(event.type == SDL_MOUSEWHEEL){
    //    if     (event.wheel.y > 0){ zoom*=VIEW_ZOOM_STEP; }
    //    else if(event.wheel.y < 0){ zoom/=VIEW_ZOOM_STEP; }
    //}

    mouseState = SDL_GetMouseState(0,0);


    if(event.type == SDL_MOUSEWHEEL){
        if     (event.wheel.y > 0){ zoom/=VIEW_ZOOM_STEP; }
        else if(event.wheel.y < 0){ zoom*=VIEW_ZOOM_STEP; }
    }
    gui.onEvent(mouseX,mouseY,event);
    Vec2f pmouse=(Vec2f){mouse_begin_x,mouse_begin_y};


    brush->eventHandling(event, pmouse, mouseState );

    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_m:  edit_mode = (EDIT_MODE)((((int)edit_mode)+1)%((int)EDIT_MODE::size)); printf("edit_mode %i\n", (int)edit_mode); break;
                //case SDLK_h:  warrior1->tryJump(); break;
                //case SDLK_l:    onSelectLuaShipScript.GUIcallback(lstLuaFiles); break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:break;
                case SDL_BUTTON_RIGHT: break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:break;
                case SDL_BUTTON_RIGHT:break;
            }
            break;
        case SDL_MOUSEMOTION:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:break;
                case SDL_BUTTON_RIGHT:break;
            }
        break;
    };
    AppSDL2OGL::eventHandling( event );

    //printf( "compGui %li \n", compGui );
    //if(compGui) if( compGui->check() ){ renderShip(); }

}

// ===================== MAIN

PainterGUI * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);


	// https://www.opengl.org/discussion_boards/showthread.php/163904-MultiSampling-in-SDL
	//https://wiki.libsdl.org/SDL_GLattr
	//SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 8);
    //SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 8);
    //SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 8);
    //SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE, 8);
    //SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    //SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    //SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 1);
    //SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 2);
    //SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 4);
    //SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 8);
    //SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 16);
    //glEnable(GL_MULTISAMPLE);


	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
    SDL_DisplayMode dm;
    SDL_GetDesktopDisplayMode(0, &dm);
	thisApp = new PainterGUI( junk , dm.w-150, dm.h-100 );
	//thisApp = new PainterGUI( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















