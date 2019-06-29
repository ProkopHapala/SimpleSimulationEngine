
#ifndef  ScreenSDL2OGL_h
#define  ScreenSDL2OGL_h

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

const float	VIEW_ZOOM_STEP     = 1.2f;
const float	VIEW_ZOOM_DEFAULT  = 10.0f;
const float	VIEW_DEPTH_DEFAULT = 1000.0;

class ScreenSDL2OGL{
	public:
	// World2D* scene;   // TODO
	int   WIDTH;
	int   HEIGHT;
	float VIEW_DEPTH;
	float ASPECT_RATIO;
	float zoom;

	float camX0=0.0f,camY0=0.0f;
	float fWIDTH, fHEIGHT, camXmin, camYmin, camXmax, camYmax;

	int   mouseX, mouseY;
	float mouse_begin_x;
	float mouse_begin_y;
	//float mouse_end_x;
	//float mouse_end_y;

	bool GL_LOCK = false;

	bool hasFocus;
	SDL_Window*      window;
	SDL_Renderer*    renderer;

// ============ function declarations

	virtual void camera      ();
	virtual void cameraHUD   ();
	virtual void draw        ();
    virtual void drawHUD     ();
	virtual void inputHanding();
    virtual void updateMousePos ( int x, int y );

	virtual void setupRenderer ();
	virtual void setDefaults   ();

	void update( );

	void init( int& id, int WIDTH_, int HEIGHT_ );
	ScreenSDL2OGL ( int& id, int WIDTH_, int HEIGHT_ );

// === inline functions

	inline float mouseUp     ( float mY ){ return 2*zoom*( 0.5 -mY/float(HEIGHT)                    ); };
	inline float mouseUp_    ( float mY ){ return 2*zoom*(      mY/float(HEIGHT) - 0.5              ); };
	inline float mouseRight  ( float mX ){ return 2*zoom*(      mX/float(HEIGHT) - 0.5*ASPECT_RATIO ); };

	bool wait_LOCK( int n, int ms ){ if(!GL_LOCK) return true; for(int i=0; i<n; i++){ SDL_Delay(ms); if(!GL_LOCK) return true; } return false; }

};

#endif
