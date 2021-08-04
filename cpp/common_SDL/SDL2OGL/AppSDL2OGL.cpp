
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "AppSDL2OGL.h" // THE HEADER

void AppSDL2OGL::quit(){ SDL_Quit(); loopEnd=true; /*exit(1);*/ };

void AppSDL2OGL::wait(int ms){
    for( int i=0; i<ms; i+=timeSlice ){
        uint32_t tnow=SDL_GetTicks();
        if(tnow>=(upTime+ms)){upTime=tnow; break; }
        SDL_Delay(timeSlice);
    }
};


void AppSDL2OGL::loop( int n ){
	loopEnd = false;
	for( int iframe=0; iframe<n; iframe++ ){
		inputHanding();
		//if(!STOP){update();} // DEPRECATED: usually we want to stop physics, not drawing
		update();
		for(ScreenSDL2OGL* w : child_windows){
            if(w) w->update();
        }
		//printf(" %i \n", iframe );
        wait(delay);
		//if( delay>0 ) SDL_Delay( delay );
		//frameCount++;
		if(loopEnd) break;
	}
}

void AppSDL2OGL::inputHanding(){
    // https://stackoverflow.com/questions/24918626/sdl-window-input-focus-and-sdl-window-mouse-focus
    // SDL_WINDOW_INPUT_GRABBED set with SDL_SetWindowGrab and forces the mouse to stay inside the window (so the window has both mouse and keyboard focus).
    // SDL_WINDOW_INPUT_FOCUS flag indicates whether the window is active or not (has keyboard input, and probably other controller inputs too).
    // SDL_WINDOW_MOUSE_FOCUS indicates whether the mouse is hovering over the window, even if the window is not active.
    //printf("DEBUG AppSDL2OGL::inputHanding 0 \n");
    wflags = SDL_GetWindowFlags(window);
    bool bFocus =  wflags & SDL_WINDOW_INPUT_FOCUS;
    //printf( "window[%i] inputHanding   %i %i \n", id, wflags&SDL_WINDOW_INPUT_FOCUS, SDL_WINDOW_INPUT_FOCUS );
    //if( wflags & SDL_WINDOW_INPUT_FOCUS ){
    //printf( "window[%i] inputHanding \n", id );
    const Uint8 *keys = SDL_GetKeyboardState(NULL);
    if( bFocus ){
        //printf( "window[%i] inputHanding has focus \n", id );
        keyStateHandling( keys );
        mouseHandling( );
    }
    //printf("DEBUG AppSDL2OGL::inputHanding 1 \n");
    bool focuses[child_windows.size()];
    for(int i=0; i<child_windows.size(); i++){
        ScreenSDL2OGL* w = child_windows[i];
        if(w==0){ focuses[i] = 0; continue; }
        focuses[i] = SDL_WINDOW_INPUT_FOCUS & SDL_GetWindowFlags(w->window);
        if( focuses[i] ){
            w->keyStateHandling ( keys );
            w->mouseHandling( );
        }
    }
    SDL_Event		 event;
    //printf("DEBUG AppSDL2OGL::inputHanding 2 \n");
    while(SDL_PollEvent(&event)){
        if( bFocus ){ eventHandling( event ); }
        for(int i=0; i<child_windows.size(); i++){
            if( focuses[i] ){
                child_windows[i]->eventHandling( event );
                if( !child_windows[i] ) focuses[i] = 0;
            }
        }
    }
    //printf("DEBUG AppSDL2OGL::inputHanding 3 \n");
}


void AppSDL2OGL::eventHandling( const SDL_Event& event ){
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_ESCAPE:   quit(); break;
            } break;
        case SDL_WINDOWEVENT:
            switch (event.window.event) {
                case SDL_WINDOWEVENT_CLOSE: quit(); break;
            } break;
        case SDL_QUIT: quit(); break;
    };
    ScreenSDL2OGL::eventHandling(event);
};

/*
void AppSDL2OGL::eventHandling( const SDL_Event& event ){
    switch( event.type ){
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:  LMB = true;  break;
                case SDL_BUTTON_RIGHT: RMB = true;  break;
            };  break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:  LMB = false; break;
                case SDL_BUTTON_RIGHT: RMB = false; break;
            }; break;
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_ESCAPE:   quit(); break;
                //case SDLK_SPACE:    STOP = !STOP; printf( STOP ? " STOPED\n" : " UNSTOPED\n"); break;
                case SDLK_KP_MINUS: zoom*=VIEW_ZOOM_STEP; break;
                case SDLK_KP_PLUS:  zoom/=VIEW_ZOOM_STEP; break;
            } break;
        case SDL_QUIT: quit(); break;
    };
};

void AppSDL2OGL::keyStateHandling( const Uint8 *keys ){
    if( keys[ SDL_SCANCODE_LEFT  ] ){ camX0 -= camStep; }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ camX0 += camStep; }
	if( keys[ SDL_SCANCODE_UP    ] ){ camY0 += camStep; }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ camY0 -= camStep; }
};

void AppSDL2OGL::mouseHandling( ){
    SDL_GetMouseState( &mouseX, &mouseY ); //mouseY=HEIGHT-mouseY; // this is done in mouseUp()
    defaultMouseHandling( mouseX, mouseY );
}
*/

void AppSDL2OGL::removeChild(ScreenSDL2OGL* child){
    for(int i=0; i<child_windows.size(); i++){
        if(child_windows[i]==child){
            child_windows[i] = 0;
        }
    }
}

AppSDL2OGL::AppSDL2OGL( int& id, int WIDTH_, int HEIGHT_, const char* name ) : ScreenSDL2OGL( id, WIDTH_, HEIGHT_, name ) {

}

