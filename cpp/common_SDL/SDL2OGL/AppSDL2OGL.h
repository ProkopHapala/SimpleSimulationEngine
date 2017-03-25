
#ifndef  AppSDL2OGL_h
#define  AppSDL2OGL_h

#include "ScreenSDL2OGL.h"

const float	VIEW_MOVE_STEP     = 0.2f;

class AppSDL2OGL : public ScreenSDL2OGL{
	public:
	bool LMB=false,RMB=false;
	int frameCount = 0;
	bool loopEnd   = false, STOP = false;
	float camStep  = VIEW_MOVE_STEP;
	int delay = 10;

// ============ function declarations

	virtual void quit(       );
	virtual void loop( int n );
	virtual void inputHanding();
	virtual void eventHandling   ( const SDL_Event& event               );
	virtual void keyStateHandling( const Uint8 *keys                    );
	virtual void mouseHandling   ( );
	void defaultMouseHandling    ( const int& mouseX, const int& mouseY );

	AppSDL2OGL( int& id, int WIDTH_, int HEIGHT_ );

};

#endif
