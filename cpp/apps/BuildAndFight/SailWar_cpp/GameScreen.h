
#ifndef GameScreen_h
#define GameScreen_h

#include "Screen2D.h"

#include "GameWorld.h"
#include "Frigate2D.h"

class GameScreen : public Screen2D {
	public:
	GameWorld * world;
	Frigate2D * thisShip;

	// ==== function declarations

	virtual void draw   ();
	virtual void drawHUD();
	GameScreen( int& id, int WIDTH_, int HEIGHT_ );

};

#endif  // #ifndef GameScreen_h

