
#ifndef  AppSDL2OGL_h
#define  AppSDL2OGL_h

#include <vector>
#include "ScreenSDL2OGL.h"



inline SDL_DisplayMode initSDLOGL( int multi_sample=8 ){
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
    if( multi_sample>0 ){ 
        SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, multi_sample );
        glEnable(GL_MULTISAMPLE);
    }
	//SDL_SetRelativeMouseMode( SDL_TRUE );
    SDL_DisplayMode dm;
    SDL_GetDesktopDisplayMode(0, &dm);
    return dm;
}



class AppSDL2OGL : public ScreenSDL2OGL{ public:
	bool LMB=false,RMB=false;
	int  upTime=0,delay=20,timeSlice=5; //,frameCount=0;
	bool loopEnd    = false, STOP = false;
	//float camStep   = VIEW_MOVE_STEP;

    std::vector<ScreenSDL2OGL*> child_windows;

// ============ function declarations

    void wait(float ms);
	virtual void quit(       );
	void         wait(int ms);
	virtual void loop( int n );
	virtual void inputHanding();
	virtual void eventHandling   ( const SDL_Event& event               );
	//virtual void keyStateHandling( const Uint8 *keys                    );
	//virtual void mouseHandling   ( );
	//void defaultMouseHandling    ( const int& mouseX, const int& mouseY );
    virtual void removeChild(ScreenSDL2OGL* child);

	AppSDL2OGL( int& id, int WIDTH_, int HEIGHT_, const char* name=0 );

};

#endif
