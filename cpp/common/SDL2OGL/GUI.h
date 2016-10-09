
#ifndef  GUI_h
#define  GUI_h

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "Draw.h"
#include "Draw2D.h"
#include "Draw3D.h"

#include <string>

// ==============================
//    class GUITextInput
// ==============================

class GUITextInput{
    public:

    int         curPos=0;
	std::string inputText;

    bool     isNumber=true,checkRange=false;
	float    vmin=0.0f, vmax=1.0f;
	double   value=0.0d;
	SDL_Keycode num_op = 0;

	bool     modified=true,entered=false;

	// ==== functions

	void         applyVal( float f );

    virtual void view3D( const Vec3d& pos, int fontTex, float textSize );
    virtual void onKeyDown( SDL_Event e );
	virtual void onText( SDL_Event e );

};

// ==============================
//    class GUIAbstractPanel
// ==============================

class GUIAbstractPanel{
    public:
	int  xmin=256,xmax=128,ymin=0,ymax=0;
	bool visible=true, disabled=false;

	uint32_t bgColor=0xA0A0A0, textColor=0x000000;

	bool     redraw=true;
	int      gllist=0;

	int      fontTex=0;
    char*    caption=NULL;

	// ==== functions

    inline bool check      ( int  x, int  y ){  return (x>xmin)&&(x<xmax)&&(y>ymin)&&(y<ymax); }
	inline void toRelative ( int& x, int& y ){ x-=xmin; y-=ymin; }

	virtual void              onKeyDown( SDL_Event e )                 = 0;
	virtual GUIAbstractPanel* onMouse( int x, int y, SDL_Event event ) = 0;
    virtual void              onText( SDL_Event e )                    = 0;

    virtual void draw  ( );
    virtual void tryRender();
    virtual void init( int xmin_, int ymin_, int xmax_, int ymax_, int fontTex_ );

};

// ==============================
//       class GUIPanel
// ==============================

class GUIPanel : public GUIAbstractPanel {
    public:
	bool isSlider=true, isButton=false;

	uint32_t barColor=0x00FF00;

	bool     executed=false;
	int      curPos=0;
	std::string inputText;

	float    vmin=0.0f, vmax=1.0f;
	double   value=0.0d;

	void (*command)(double) = NULL;

    // ==== functions

    virtual void draw  ( );
	virtual void tryRender();
    virtual void              onKeyDown( SDL_Event e );
    virtual void              onText( SDL_Event e );
    virtual GUIAbstractPanel* onMouse( int x, int y, SDL_Event event );

	// ===== inline functions
	inline double x2val( float  x   ){ return ( x*(vmax-vmin)/(xmax-xmin) )+ vmin; };
	inline float  val2x( double val ){ return (val-vmin)*(xmax-xmin)/(vmax-vmin);  };

};

// ==============================
//     class  MultiPanel
// ==============================

class MultiPanel : public GUIAbstractPanel {
    public:
    int nsubs;
    GUIPanel ** subs;

    // ==== functions

    void initMulti( int xmin_, int ymin_, int xmax_, int ymax_, int fontTex_, int nsubs_ );

    virtual void draw  ( );
    virtual void tryRender( );
    virtual GUIAbstractPanel* onMouse  ( int x, int y, SDL_Event event );

    virtual void onKeyDown( SDL_Event e ){};
    virtual void onText   ( SDL_Event e ){};

};

#endif
