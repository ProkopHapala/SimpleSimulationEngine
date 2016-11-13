
#ifndef AeroCraftGUI_h
#define AeroCraftGUI_h

#include "ScreenSDL2OGL_3D.h"
#include "AppSDL2OGL_3D.h"

#include "GUI.h"

#include "AeroCraftWorld.h"
#include "AeroCraft.h"

//class AeroCraftGUI : public ScreenSDL2OGL_3D {
class AeroCraftGUI : public AppSDL2OGL_3D {
	public:
	AeroCraftWorld * world;

    int      fontTex;
    GUIPanel   panel;
    MultiPanel mpanel;
    GUITextInput txt;

    GUIAbstractPanel*  focused = NULL;

    bool mouseSteer = false;

	// ==== function declarations

	void renderSkyBox();

	virtual void camera     ();
	//virtual void cameraHUD();
	virtual void draw   ();
	virtual void drawHUD();
	AeroCraftGUI( int& id, int WIDTH_, int HEIGHT_ );

};

#endif  // #ifndef GameScreen_h

