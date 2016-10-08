
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
    GUIPanel panel;

	// ---- function declarations
	virtual void draw   ();
	virtual void drawHUD();
	virtual void eventHandling ( const SDL_Event& event  );

	TestAppGUI( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppGUI::TestAppGUI( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    fontTex = makeTexture( "/home/prokop/git/SimpleSimulationEngine/cpp/common_resources/dejvu_sans_mono.bmp" );

    panel.init( 5,5,100,40,  fontTex );
    panel.caption   = "rotation [Rad]"; panel.vmin = -3.14159265359; panel.vmax = 3.14159265359;
    //sprintf(panel.val_text, "NA" );
    panel.command = &command_example; panel.isButton = true;

    SDL_StartTextInput ();
    //panel.nChars = 6;
}

void TestAppGUI::draw   (){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glDisable ( GL_LIGHTING );

	Draw3D::drawAxis ( 3.0f );

};

void TestAppGUI::drawHUD(){
	if( panel.redraw ) panel.render();
	glCallList( panel.gllist );

	//glColor3f(1.0f,1.0f,1.0f);
	//panel.render();
	//Draw2D::drawRectangle ( 100, 100, 300, 200, false    );
	//Draw2D::drawLine      ( {0.0f,0.0f},{100.0f, 200.0f} );
}

void TestAppGUI::eventHandling ( const SDL_Event& event  ){
    //printf( "NBodyWorldApp::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN:
            panel.onKeyDown( event ); break;
        case SDL_TEXTINPUT:
            panel.onText( event ); break;
        case SDL_MOUSEBUTTONDOWN:
            //printf( "%i %i\n", mouseX, mouseY );

            bool button_clicked = panel.onMouse( mouseX, mouseY, event );
            if  (button_clicked) printf("button_clicked!!! panel.value= %3.6f\n", panel.value );
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
















