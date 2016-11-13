
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "SDL_utils.h"

#include "fastmath.h"
#include "Vec2.h"
#include "Solids.h"

#include "Draw.h"
#include "Draw2D.h"
#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"

#include "GUI.h"

// ======================  TestApp


void command_example( double value ){
    printf( "I'm writting value : %6.6f \n", value );
}

class TestAppGUI : public AppSDL2OGL_3D {
	public:
    int      fontTex;
    GUIPanel   panel;
    MultiPanel mpanel;
    GUITextInput txt;

    GUIAbstractPanel*  focused = NULL;

	// ---- function declarations
	virtual void draw   ();
	virtual void drawHUD();
	virtual void eventHandling ( const SDL_Event& event  );

	TestAppGUI( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppGUI::TestAppGUI( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    //fontTex = makeTexture( "common_resources/dejvu_sans_mono.bmp" );
    //fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA.bmp" );
    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    //fontTex = makeTexture( "common_resources/dejvu_sans_mono_Alpha.bmp" );

    panel.init( 5,5,105,35,  fontTex );
    panel.caption   = "rotation [Rad]"; panel.vmin = -3.14159265359; panel.vmax = 3.14159265359;
    //sprintf(panel.val_text, "NA" );
    panel.command = &command_example; panel.isButton = true;

    mpanel.initMulti( 120,5,200,120, fontTex , 4 );
    mpanel.caption="MultiPanel_1";

    txt.inputText = "insert number using =+-*/";

    SDL_StartTextInput ();
    //panel.nChars = 6;
}

void TestAppGUI::draw   (){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glDisable ( GL_LIGHTING );

	//Draw3D::drawText( txt.inputText.c_str(), {1.0,1.0,1.0}, fontTex, 0.3, 0, 0 );

	Draw3D::drawAxis ( 3.0f );

	glColor4f(1.0f,1.0f,1.0f,0.9f);
	txt.view3D( {1.0,1.0,1.0}, fontTex, 0.3 );

	//Draw3D::drawText( "AHOJ!\0", {1.0,1.0,1.0}, fontTex, 5.0, 0, 0 );

};

void TestAppGUI::drawHUD(){

	panel .tryRender();  panel.draw();
	mpanel.tryRender(); mpanel.draw();
	if(focused) Draw2D::drawRectangle(focused->xmin,focused->ymin,focused->xmax,focused->ymax,false);

	//glColor3f(1.0f,1.0f,1.0f);
	//panel.render();
	//Draw2D::drawRectangle ( 100, 100, 300, 200, false    );
	//Draw2D::drawLine      ( {0.0f,0.0f},{100.0f, 200.0f} );
}

void TestAppGUI::eventHandling ( const SDL_Event& event  ){
    //printf( "NBodyWorldApp::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN:
            if(focused){ focused->onKeyDown( event ); }else{ txt.onKeyDown(  event ); }; break;
        case SDL_TEXTINPUT:
            if(focused){ focused->onText   ( event ); }else{ txt.onText   ( event );  }; break;
        case SDL_MOUSEBUTTONDOWN:
            //printf( "%i %i\n", mouseX, mouseY );
            GUIAbstractPanel* active = NULL; focused=NULL;
            active = mpanel.onMouse( mouseX, mouseY, event ); if(active) focused=active;
            active = panel .onMouse( mouseX, mouseY, event ); if(active) focused=active;
            break;
    };
    AppSDL2OGL::eventHandling( event );
    //camStep = zoom*0.05;
}

// ===================== MAIN

TestAppGUI * testApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestAppGUI( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}
















