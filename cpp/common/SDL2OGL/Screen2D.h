
#ifndef  Screen2D_h
#define  Screen2D_h

const float	VIEW_ZOOM_STEP     = 1.2f;
const float	VIEW_ZOOM_DEFAULT  = 10.0f;
const float	VIEW_DEPTH_DEFAULT = 1000.0;

class Screen2D{
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

	void camera   ();
	void cameraHUD();

	void update( );
	virtual void draw   ();
    virtual void drawHUD();
	void inputHanding   ();

	void setupRenderer();
	void setDefaults();
	void init( int& id, int WIDTH_, int HEIGHT_ );
	Screen2D ( int& id, int WIDTH_, int HEIGHT_ );

// === inline functions

	inline float mouseUp     ( float mY ){ return 2*zoom*( 0.5 -mY/float(HEIGHT)                    ); };
	inline float mouseRight  ( float mX ){ return 2*zoom*(      mX/float(HEIGHT) - 0.5*ASPECT_RATIO ); };

};

#endif
