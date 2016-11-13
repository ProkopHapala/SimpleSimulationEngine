
#ifndef AeroCraftGUI_h
#define AeroCraftGUI_h

#include "ScreenSDL2OGL_3D.h"

#include "AeroCraftWorld.h"
#include "AeroCraft.h"

class AeroCraftGUI : public ScreenSDL2OGL_3D {
	public:
	int default_font_texture;
	AeroCraftWorld * world;

	// ==== function declarations

	void renderSkyBox();

	virtual void camera     ();
	//virtual void cameraHUD();
	virtual void draw       ();
	//virtual void drawHUD  ();
	AeroCraftGUI( int& id, int WIDTH_, int HEIGHT_ );




};

#endif  // #ifndef GameScreen_h

