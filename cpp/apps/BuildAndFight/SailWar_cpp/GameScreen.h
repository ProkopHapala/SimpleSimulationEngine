
#ifndef GameScreen_h
#define GameScreen_h

#include "Screen2D.h"

#include "GameWorld.h"

class GameScreen : public Screen2D {
	public:
	GameWorld * world;

	// ==== function declarations

	virtual void draw();
	GameScreen( int& id, int WIDTH_, int HEIGHT_ );

};

#endif  // #ifndef GameScreen_h

