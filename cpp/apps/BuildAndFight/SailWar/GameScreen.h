
#ifndef GameScreen_h
#define GameScreen_h

#include "ScreenSDL2OGL.h"

#include "GameWorld.h"
#include "Frigate2D.h"

class GameScreen : public ScreenSDL2OGL {
	public:
	GameWorld * world    = NULL;
	Frigate2D * thisShip = NULL; // without NULL there were segfaults

	// ==== function declarations

	virtual void draw   ();
	virtual void drawHUD();
	GameScreen( int& id, int WIDTH_, int HEIGHT_ );

};

#endif  // #ifndef GameScreen_h

