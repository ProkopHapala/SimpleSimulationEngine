
#ifndef  ScreenSDL2OGL_h
#define  ScreenSDL2OGL_h

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

const float	VIEW_ZOOM_STEP     = 1.2f;
const float	VIEW_ZOOM_DEFAULT  = 10.0f;
const float	VIEW_DEPTH_DEFAULT = 1000.0;

const float	VIEW_MOVE_STEP     = 0.2f;


//void setupOpenGLglobals();
void setLightingRGB();
void setLightingNormal();

class ScreenSDL2OGL{
	public:
	// World2D* scene;   // TODO
	int   id;
	ScreenSDL2OGL* parent = 0;
	int   iinparent=0;
	int frameCount=0;

	int   WIDTH;
	int   HEIGHT;
	float VIEW_DEPTH;
	float ASPECT_RATIO;
	float zoom;

    bool LMB=false,RMB=false;
	float camStep   = VIEW_MOVE_STEP;

	float camX0=0.0f,camY0=0.0f;
	float fWIDTH, fHEIGHT, camXmin, camYmin, camXmax, camYmax;

	int   mouseX, mouseY;
	float mouse_begin_x;
	float mouse_begin_y;
	//float mouse_end_x;
	//float mouse_end_y;

	bool GL_LOCK = false;

	//bool hasFocus;
	Uint32         wflags;
	SDL_Window*    window;
	SDL_GLContext  glctx;
	//SDL_Renderer*    renderer;

// ============ function declarations

	virtual void camera      ();
	virtual void cameraHUD   ();
	virtual void draw        ();
    virtual void drawHUD     ();
	//virtual void inputHanding();
	virtual void eventHandling   ( const SDL_Event& event               );
	virtual void keyStateHandling( const Uint8 *keys                    );
	virtual void mouseHandling   ( );

    //virtual void updateMousePos ( int x, int y );

	//virtual void setupRenderer ();
	virtual void setDefaults   ();

	void update( );

	void init( int& id, int WIDTH_, int HEIGHT_ );
	ScreenSDL2OGL ( int& id, int WIDTH_, int HEIGHT_ );
	virtual ~ScreenSDL2OGL();
	virtual void removeChild(ScreenSDL2OGL* child);

// === inline functions

	inline float mouseUp     ( float mY ){ return 2*zoom*( 0.5 -mY/float(HEIGHT)                    ); };
	inline float mouseUp_    ( float mY ){ return 2*zoom*(      mY/float(HEIGHT) - 0.5              ); };
	inline float mouseRight  ( float mX ){ return 2*zoom*(      mX/float(HEIGHT) - 0.5*ASPECT_RATIO ); };

	void defaultMouseHandling    ( const int& mouseX, const int& mouseY );

	bool wait_LOCK( int n, int ms ){ if(!GL_LOCK) return true; for(int i=0; i<n; i++){ SDL_Delay(ms); if(!GL_LOCK) return true; } return false; }

};

#endif
