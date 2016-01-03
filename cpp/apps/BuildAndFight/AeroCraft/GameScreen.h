
#ifndef GameScreen_h
#define GameScreen_h

#include "ScreenSDL2OGL_3D.h"

#include "GameWorld.h"
#include "AeroCraft.h"

class GameScreen : public ScreenSDL2OGL_3D {
	public:
	GameWorld * world;
	AeroCraft * thisShip;

	// ==== function declarations

	void renderSkyBox();

	virtual void camera     ();
	//virtual void cameraHUD();
	virtual void draw       ();
	//virtual void drawHUD  ();
	GameScreen( int& id, int WIDTH_, int HEIGHT_ );




};

#endif  // #ifndef GameScreen_h

