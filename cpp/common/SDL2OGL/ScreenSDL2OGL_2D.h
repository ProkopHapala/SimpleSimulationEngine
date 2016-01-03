
#ifndef  ScreenSDL2OGL_2D_h
#define  ScreenSDL2OGL_2D_h

const float	VIEW_ZOOM_STEP     = 1.2f;
const float	VIEW_ZOOM_DEFAULT  = 10.0f;
const float	VIEW_DEPTH_DEFAULT = 1000.0;

class ScreenSDL2OGL_2D : public ScreenSDL2OGL {
	public:
	// World2D* scene;   // TODO
	int   WIDTH;
	int   HEIGHT;
	float VIEW_DEPTH;
	float ASPECT_RATIO;
	float zoom;

	int   mouseX, mouseY;
	float mouse_begin_x ;
	float mouse_begin_y ;

	bool hasFocus;
	SDL_Window*      window;
	SDL_Renderer*    renderer;

// ============ function declarations

	virtual void camera       ();
	virtual void cameraHUD    ();
	virtual void draw         ();
    virtual void drawHUD      ();
	virtual void inputHanding ();

	virtual void setupRenderer();
	virtual void setDefaults  ();

	void update( );

	void init( int& id, int WIDTH_, int HEIGHT_ );
	ScreenSDL2OGL_2D ( int& id, int WIDTH_, int HEIGHT_ );

// === inline functions

	inline float mouseUp     ( float mY ){ return 2*zoom*( 0.5 -mY/float(HEIGHT)                    ); };
	inline float mouseRight  ( float mX ){ return 2*zoom*(      mX/float(HEIGHT) - 0.5*ASPECT_RATIO ); };

};

#endif
