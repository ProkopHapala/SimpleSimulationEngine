
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "AppSDL2OGL.h" // THE HEADER

void AppSDL2OGL::quit(){ SDL_Quit(); exit(1); };

void AppSDL2OGL::loop( int n ){
	loopEnd = false;
	for( int iframe=0; iframe<n; iframe++ ){
		inputHanding();
		if(!STOP){
			update();
		}
		//printf(" %i \n", iframe );
		SDL_Delay( 10 );
		frameCount++;
		if(loopEnd) break;
	}
}

void AppSDL2OGL::inputHanding(){
	SDL_Event		 event;
	while(SDL_PollEvent(&event)){ // be carefull to include all event handling instide the while loop !!!
		if( event.type == SDL_KEYDOWN ){
			switch( event.key.keysym.sym ){
				case SDLK_ESCAPE:   quit(); break;
				case SDLK_SPACE:    STOP = !STOP; printf( STOP ? " STOPED\n" : " UNSTOPED\n"); break;
				case SDLK_KP_MINUS: zoom*=VIEW_ZOOM_STEP; break;
				case SDLK_KP_PLUS:  zoom/=VIEW_ZOOM_STEP; break;
			}
		}
		if( event.type == SDL_QUIT){ quit();  };

	} // while(SDL_PollEvent(&event))
	const Uint8 *keys = SDL_GetKeyboardState(NULL);
	if( keys[ SDL_SCANCODE_LEFT  ] ){ camX0 -= camStep; }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ camX0 += camStep; }
	if( keys[ SDL_SCANCODE_UP    ] ){ camY0 += camStep; }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ camY0 -= camStep; }
	SDL_GetMouseState( &mouseX, &mouseY );
	mouse_begin_x = mouseRight( mouseX ) + camX0;
	mouse_begin_y = mouseUp   ( mouseY ) + camY0;
    fWIDTH  = zoom*ASPECT_RATIO;
	fHEIGHT = zoom;
	camXmin = camX0 - fWIDTH; camYmin = camY0 - fHEIGHT;
	camXmax = camX0 + fWIDTH; camYmax = camY0 + fHEIGHT;
}

AppSDL2OGL::AppSDL2OGL( int& id, int WIDTH_, int HEIGHT_ ) : ScreenSDL2OGL( id, WIDTH_, HEIGHT_ ) {

}

