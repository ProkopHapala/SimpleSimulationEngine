
#ifndef  ScreenSDL2OGL_3D_h
#define  ScreenSDL2OGL_3D_h

#include "Vec3.h"
#include "Vec2.h"
#include "Mat3.h"
#include "quaternion.h"

//#include "Camera.h"
//#include "cameraOGL.h"

#include "ScreenSDL2OGL.h"

class ScreenSDL2OGL_3D : public ScreenSDL2OGL{
	public:

	//float qCamera   [4];
	//float qCameraOld[4];
	//Vec3d camDir, camUp, camRight;

	Quat4f qCamera   ;
	Quat4f qCameraOld;
	Mat3f  camMat;

	bool  mouse_spinning;

	// ==== function declarations

	virtual void camera      ();

	// TODO
	//void getCameraDirections();
	//void mouse_camera ( float x, float y );

	ScreenSDL2OGL_3D ( int& id, int WIDTH_, int HEIGHT_ );


	// ==== inline functions

	void  startSpining ( float x, float y              ){ mouse_spinning = true; mouse_begin_x  = x; mouse_begin_y  = y;	}
	void  endSpining   (                               ){ mouse_spinning = false;	                                    }
	//void  projectMouse ( float mX, float mY, Vec3d& mp ){ mp.set_lincomb( mouseRight(mX), camRight,  mouseUp(mY), camUp ); };
	void  projectMouse ( float mX, float mY, Vec3f& mp ){ mp.set_lincomb( mouseRight(mX), camMat.a,  mouseUp(mY), camMat.b ); };

};

#endif







