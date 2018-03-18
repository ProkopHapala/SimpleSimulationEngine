
#ifndef  AppSDL2OGL_h
#define  AppSDL2OGL_h

#include "ScreenSDL2OGL.h"

const float	VIEW_MOVE_STEP     = 0.2f;

class AppSDL2OGL : public ScreenSDL2OGL{
	public:
	bool LMB=false,RMB=false;
	int  upTime=0,delay=20,timeSlice=5,frameCount=0;
	bool loopEnd    = false, STOP = false;
	float camStep   = VIEW_MOVE_STEP;


// ============ function declarations

    void wait(float ms);
	virtual void quit(       );
	void         wait(int ms);
	virtual void loop( int n );
	virtual void inputHanding();
	virtual void eventHandling   ( const SDL_Event& event               );
	virtual void keyStateHandling( const Uint8 *keys                    );
	virtual void mouseHandling   ( );
	void defaultMouseHandling    ( const int& mouseX, const int& mouseY );

	AppSDL2OGL( int& id, int WIDTH_, int HEIGHT_ );

};

#endif
