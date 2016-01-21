
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "drawMath2D.h"

#include "AppSDL2OGL.h"
#include "testUtils.h"

class NBodyWorldApp : public AppSDL2OGL {
	public:

	virtual void draw   ();
	virtual void drawHUD();
	NBodyWorldApp( int& id, int WIDTH_, int HEIGHT_ );

};

NBodyWorldApp::NBodyWorldApp( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

}

void NBodyWorldApp::draw(){

};

void NBodyWorldApp::drawHUD(){

}

// ===================== MAIN

NBodyWorldApp * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	thisApp = new NBodyWorldApp( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















