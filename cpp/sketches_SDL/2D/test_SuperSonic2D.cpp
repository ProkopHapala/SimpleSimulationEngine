
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "AppSDL2OGL.h"

#include "fastmath.h"
#include "Vec2.h"
//#include "geom2D.h"
#include "SuperSonic2D.h"
//#include "Fluid2D.h"
//#include "Fluid2D.cpp"
#include "testUtils.h"


#include "SDL_utils.h"
#include "GUI.h"



class TestAppSuperSonic2D : public AppSDL2OGL{
	public:
    int job_type = 1;
    int perframe = 3;
    //TerrainHydraulics terrain;
    int shape;
    bool running = true;

    SuperSonic2D solver;

	// ---- function declarations
	virtual void draw   ();
	virtual void drawHUD();
    virtual void eventHandling( const SDL_Event& event );
	//virtual int tileToList( float x0, float y0, float x1, float y1 );

    //void renderMapContent ( float x0, float y0, float scale, float csc, float hsc );
    //double terrain_color( int i );
	//void drawSimplexGrid( int n, float step );

	TestAppSuperSonic2D ( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppSuperSonic2D ::TestAppSuperSonic2D ( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

    GUI_fontTex = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );

    //((GUIPanel*)gui.addPanel( new GUIPanel( "dt [Rad]", 5,5,105,35, true, true ) )) -> setRange(0.0, 0.2);
    //    ->command = &command_example;
    //panel_dt.initPanel( "dt", 5,5,90,35 );
    //panel_dt.setRange(0.0, 0.2);
    //panel_dt.value = 0.05;
    //gui.addPanel( &panel_dt );



}

void TestAppSuperSonic2D ::draw(){

    //delay = 100;
    glClearColor( 1.0f, 1.0f, 1.0f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable(GL_LIGHTING);

};

void TestAppSuperSonic2D ::drawHUD(){
    glDisable( GL_LIGHTING );
    glDisable(GL_DEPTH_TEST);

    //Draw2D::drawText( "AHOJ !!!!", 0, {10, 100}, 0.0, GUI_fontTex, fontSizeDef );
    //gui.draw();
}


void TestAppSuperSonic2D ::eventHandling( const SDL_Event& event ){
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
               // case SDLK_r:  terrain.initErrosion( 0.8 ); running=true;  break;
            }
            break;
         /*
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //paintSimplex( mouse_begin_x, mouse_begin_y );
                    mouse_left  = true;
                    break;
                case SDL_BUTTON_RIGHT:
                    mouse_right = true;
                    //eraseSimplex( mouse_begin_x, mouse_begin_y );
                    break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    mouse_left = false;
                    break;
                case SDL_BUTTON_RIGHT:
                    mouse_right = false;
                    break;
            }
            break;
        */
    };
    //gui.onEvent( mouseX, HEIGHT-mouseY, event );
    AppSDL2OGL::eventHandling( event );
};

// ===================== MAIN

TestAppSuperSonic2D * testApp;

int main(int argc, char *argv[]){

	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestAppSuperSonic2D ( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}
