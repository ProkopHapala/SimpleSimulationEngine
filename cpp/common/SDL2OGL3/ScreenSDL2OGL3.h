
#ifndef  ScreenSDL2OGL3_h
#define  ScreenSDL2OGL3_h

#include <fastmath.h>
#include <Vec2.h>
#include <Vec3.h>
#include <Mat3.h>
#include <quaternion.h>

#define GL_GLEXT_PROTOTYPES
#include <GL/gl.h>
#include <SDL2/SDL.h>

#include <Shader.h>
#include <GLObject.h>
#include <SceneNode.h>

const float	VIEW_ZOOM_STEP     = 1.2f;
const float	VIEW_ZOOM_DEFAULT  = 10.0f;
const float	VIEW_DEPTH_DEFAULT = 1000.0;

class ScreenSDL2OGL3{
	public:

    Quat4d qCam;
	Vec3d  camPos;
	Mat3d  camRot;
	Vec2d  tgFrustrum; // tangent of scene frustrum

    bool deffered;
	Shader   * canvasShader;
	GLObject * canvasQuad;
	GLuint canvas_depth_tex,canvas_color_tex;
    GLuint FramebufferName = 0;
    GLuint depthrenderbuffer;

    SceneNode3D * scene;

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

	bool hasFocus;
	SDL_Window*      window;
	SDL_Renderer*    renderer;
	SDL_GLContext    context = NULL;

// ============ function declarations

	virtual void draw        ();
	virtual void drawDeffered();

	virtual void setupDefferedRender();
	virtual void setDefaults  ();

	void update( );
	bool checkFramebufferStatus();

	int init( int& id, int WIDTH_, int HEIGHT_ );
	ScreenSDL2OGL3( int& id, int WIDTH_, int HEIGHT_ );

};

#endif
