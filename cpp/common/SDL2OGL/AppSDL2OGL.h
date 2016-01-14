
#ifndef  AppSDL2OGL_h
#define  AppSDL2OGL_h

#include "ScreenSDL2OGL.h"

const float	VIEW_MOVE_STEP     = 0.2f;

class AppSDL2OGL : public ScreenSDL2OGL{
	public:
	int frameCount = 0;
	bool loopEnd   = false, STOP = false;
	float camStep  = VIEW_MOVE_STEP; 

// ============ function declarations

	virtual void quit(       );
	virtual void loop( int n );
	virtual void inputHanding();

	AppSDL2OGL( int& id, int WIDTH_, int HEIGHT_ );

};

#endif
