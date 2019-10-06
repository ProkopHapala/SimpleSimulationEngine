
#ifndef  ScreenSDL2OGL_3D_h
#define  ScreenSDL2OGL_3D_h

#include "Vec3.h"
#include "Vec2.h"
#include "Mat3.h"
#include "quaternion.h"

#include "Camera.h"
#include "cameraOGL.h"
#include "ScreenSDL2OGL.h"

class ScreenSDL2OGL_3D : public ScreenSDL2OGL{ public:

	bool mouseSpinning = false;
	//Quat4f qCamera;
	//Mat3f  camMat;
	float  mouseRotSpeed   = 0.001;
	float  keyRotSpeed     = 0.01;
    float  cameraMoveSpeed = 0.2f;
    Quat4f qCamera;
	//Mat3d  camMat;
	//Vec3d  camPos;
	Camera cam;

	float camDist = 50.0;
	Vec2i spinning_start;

	bool perspective  = false;
	bool first_person = false;

	bool  mouse_spinning;

	// ==== function declarations
	virtual void draw     ();
	virtual void camera      ();
    virtual void eventHandling   ( const SDL_Event& event               );
	virtual void keyStateHandling( const Uint8 *keys                    );
	virtual void mouseHandling   ( );

	// TODO
	//void getCameraDirections();
	//void mouse_camera ( float x, float y );

	ScreenSDL2OGL_3D ( int& id, int WIDTH_, int HEIGHT_ );


	// ==== inline functions

	void  startSpining ( float x, float y              ){ mouse_spinning = true; mouse_begin_x  = x; mouse_begin_y  = y;	}
	void  endSpining   (                               ){ mouse_spinning = false;	                                    }
	//void  projectMouse ( float mX, float mY, Vec3d& mp ){ mp.set_lincomb( mouseRight(mX), camRight,  mouseUp(mY), camUp ); };
	//void  projectMouse ( float mX, float mY, Vec3f& mp ){ mp.set_lincomb( mouseRight(mX), camMat.a,  mouseUp(mY), camMat.b ); };

};

#endif







